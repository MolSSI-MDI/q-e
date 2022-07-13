MODULE MDI_IMPLEMENTATION
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_bcast, mp_sum
  USE mp_pools,         ONLY : intra_pool_comm, me_pool
  USE mdi,              ONLY : MDI_Init, MDI_Send, MDI_INT, MDI_CHAR, MDI_NAME_LENGTH, &
       MDI_Accept_communicator, MDI_Recv_command, MDI_Recv, &
       MDI_Set_execute_command_func, MDI_MPI_get_world_comm, MDI_DOUBLE, MDI_BYTE, &
       MDI_ENGINE, MDI_Get_role, MDI_Register_command, MDI_Register_node, &
       MDI_Register_callback, MDI_COMMAND_LENGTH, MDI_MPI_get_world_comm, &
       MDI_Plugin_get_argc, MDI_Plugin_get_arg
  USE mdi_engine,       ONLY : is_mdi, mdi_forces, &
       rid, rid_old, get_mdi_options, socket, mdi_exit_flag

  IMPLICIT NONE

  LOGICAL :: cell_changed = .true.
  LOGICAL :: mdi_first_scf = .true.
  REAL*8 :: omega_reset

  CONTAINS

  FUNCTION MDI_Plugin_init_qemdi() bind ( C, name="MDI_Plugin_init_qemdi" )
  !USE mpi,               ONLY : MPI_Comm_rank
  !USE parallel_include,  ONLY : MPI_Comm_rank, MPI_COMM_WORLD, MPI_Allreduce, MPI_MAX, &
  !        MPI_IN_PLACE
  USE ISO_C_binding
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: set_command_line
  USE input_parameters,  ONLY : prefix
  USE parallel_include

  ! MDI Communicator to the driver
  INTEGER :: comm

  ! MPI intra-communicator for this code
  INTEGER :: world_comm

  ! Flag to terminate MDI response function
  LOGICAL :: terminate_flag = .false.

    INTEGER :: MDI_Plugin_init_qemdi
    INTEGER :: ierr
    INTEGER :: argc
    INTEGER :: iarg
    CHARACTER(LEN=1024) :: option
    CHARACTER(LEN=1024) :: mdi_options
    CHARACTER(LEN=1024) :: input_file
    LOGICAL :: mdi_options_found = .false.
    LOGICAL :: input_file_found = .false.

    INTEGER                :: nim, npt, npl, nta, nbn, ndg
    !CHARACTER(LEN=80)      :: infile
    INTEGER                :: retval
    INTEGER                :: world_rank, world_rank_reduce
    CHARACTER(LEN=80)      :: suffix


    WRITE(6,*)"IN PLUGIN_INIT"
    FLUSH(6)

    ! Get the command-line options from the driver
    mdi_options_found = .false.
    CALL MDI_Plugin_get_argc(argc, ierr)
    DO iarg=0, argc-1
       CALL MDI_Plugin_get_arg(iarg, option, ierr)
       IF ( (TRIM(option) .eq. "-mdi") .or. (TRIM(option) .eq. "--mdi") ) THEN
          IF ( argc .gt. (iarg+1) ) THEN
             CALL MDI_Plugin_get_arg(iarg+1, mdi_options, ierr)
             mdi_options_found = .true.
          ELSE
             WRITE(6,*)'ERROR: argument to -mdi option not provided'
             MDI_Plugin_init_qemdi = 1
             RETURN
          END IF
       END IF
       IF ( (TRIM(option) .eq. "-in") .or. (TRIM(option) .eq. "--in") ) THEN
          IF ( argc .gt. (iarg+1) ) THEN
             CALL MDI_Plugin_get_arg(iarg+1, input_file, ierr)
             input_file_found = .true.
          ELSE
             WRITE(6,*)'ERROR: argument to -in option not provided'
             MDI_Plugin_init_qemdi = 1
             RETURN
          END IF
       END IF
    END DO
    IF ( .not. mdi_options_found ) THEN
       WRITE(6,*)'ERROR: --mdi option not provided'
       MDI_Plugin_init_qemdi = 1
       RETURN
    END IF
    IF ( .not. input_file_found ) THEN
       WRITE(6,*)'ERROR: --in option not provided'
       MDI_Plugin_init_qemdi = 1
       RETURN
    END IF

    ! Call MDI_Init
    CALL MDI_Init(mdi_options, ierr)

    ! Get the MPI intra-communicator over which this plugin will run
    CALL MDI_MPI_get_world_comm(world_comm, ierr);


    ! Register all supported MDI Commands
    CALL MDI_Register_node("@DEFAULT", ierr)
    CALL MDI_Register_command("@DEFAULT", "<CELL", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CELL", ierr)
    CALL MDI_Register_command("@DEFAULT", ">COORDS", ierr)


    ! SHOULD GET THESE THROUGH GET ARG
    nim = 1
    npl = 1
    nta = 1
    nbn = 1
    ndg = 1
    !
    CALL set_command_line( nimage=nim, npool=npl, ntg=nta, &
         nband=nbn, ndiag=ndg )
    CALL mp_startup ( my_world_comm=world_comm )
    !CALL mp_startup ( mdi_initialization = mdi_options )
    CALL environment_start ( 'PWSCF' )
    !
    !infile = '/repo/tests/water/qe.in'
    CALL read_input_file ('PW+iPi', input_file )
    !
    ! Append a label to the file prefix, so that each plugin instance uses different files
    !
    CALL MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, world_rank, 1, MPI_INT, MPI_MIN, world_comm, ierr)
    WRITE(suffix, '(I0)')world_rank
    prefix = TRIM(prefix) // '_'
    prefix = TRIM(prefix) // TRIM(suffix)
    !
    ! Start a PW calculation, which will listen for MDI commands
    !
    retval = 0
    CALL mdi_listen( retval )
    WRITE(6,*)'AFTER MDI LISTEN'
    !
    MDI_Plugin_init_qemdi = retval
    !
  END FUNCTION MDI_Plugin_init_qemdi



  SUBROUTINE mdi_listen ( exit_status )
    !!
    !! Driver for IPI
    !!
    USE ISO_C_BINDING
    !
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE parameters,       ONLY : ntypx, npk
    USE check_stop,       ONLY : check_stop_init
    USE mp_global,        ONLY : mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : gamma_only, conv_elec, istep, ethr, lscf, lmd, &
                                 lforce => tprnfor, tstress
    USE ions_base,        ONLY : tau, nat
    USE cell_base,        ONLY : alat, at, omega, bg
    USE cellmd,           ONLY : omega_old, at_old, calc
    USE force_mod,        ONLY : force
    USE f90sockets,       ONLY : readbuffer, writebuffer
    USE extrapolation,    ONLY : update_file, update_pot
    USE io_files,         ONLY : iunupdate, nd_nmbr, prefix, tmp_dir, postfix, &
         wfc_dir, delete_if_present, seqopn
    USE scf,              ONLY : rho
    USE lsda_mod,         ONLY : nspin
    USE fft_base,         ONLY : dfftp
    !<<<
    USE mdi,              ONLY : MDI_Accept_communicator, MDI_Set_execute_command_func, &
                                 MDI_Recv_command
    !USE command_line_options, ONLY : command_line
    !>>>
    !
    IMPLICIT NONE
    !
    !! Gives the exit status at the end
    INTEGER, INTENT(OUT) :: exit_status
    !
    ! Local variables
    INTEGER, PARAMETER :: MSGLEN=12
    REAL*8, PARAMETER :: gvec_omega_tol=1.0D-1
    LOGICAL :: isinit=.false., hasdata=.false., exst
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER :: ccmd, i, info
    REAL*8 :: sigma(3,3), at_reset(3,3), dist_reset, ang_reset
    REAL *8 :: cellih(3,3), vir(3,3), pot
    REAL*8 :: dist_ang(6), dist_ang_reset(6)
    INTEGER :: ierr
    
    ! MDI Plugin callback function
    PROCEDURE(mdi_execute_command), POINTER :: mdi_execute_command_func => null()
    TYPE(C_PTR)                         :: class_obj
    mdi_execute_command_func => mdi_execute_command

    !----------------------------------------------------------------------------
    !
    lscf      = .true.
    lforce    = .true.
    tstress   = .true.
    lmd       = .true.
    omega_reset = 0.d0
    !
    exit_status = 0
    !
    IF (ionode) CALL plugin_arguments()
    CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
    !
    ! ... needs to come before iosys() so some input flags can be
    !     overridden without needing to write PWscf specific code.
    !
    ! ... convert to internal variables
    !
    CALL iosys()
    !
    IF ( gamma_only ) WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
    !
    ! call to void routine for user defined / plugin patches initializations
    !
    CALL plugin_initialization()
    !
    CALL check_stop_init()
    CALL setup()
    !
    ! ... Initializations
    !
    CALL init_run()
    !
    ! Set MDI variables
    !
    is_mdi = .true.
    !
    IF (is_mdi) THEN
       CALL allocate_nat_arrays()
    END IF
    !
    CALL MDI_Accept_Communicator( socket, ierr )
    CALL MDI_Set_execute_command_func(mdi_execute_command_func, class_obj, ierr)
    !
    DO WHILE ( .not. mdi_exit_flag )
       !
       CALL MDI_Recv_Command( header, socket, ierr )
       CALL mp_bcast( header, ionode_id, intra_image_comm )
       !
       IF ( ionode ) write(*,*) " @ DRIVER MODE: Message from server: ", trim( header )
       !
       call mdi_execute_command(header, socket, exit_status)
       IF ( exit_status .ne. 0 ) THEN
          CALL errore('mdi_listen','execute command failed',1)
       END IF
       !
    END DO
    !
  END SUBROUTINE mdi_listen




  !<<<
  !---------------------------------------------------------------------!
  !
  !
  SUBROUTINE read_mm_charge(socketfd)
    USE qmmm, ONLY : charge_mm, tau_mask, nat_mm
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: i, ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM charges"
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( charge_mm, nat_mm, MDI_DOUBLE, socketfd, ierr )
    !
#if defined(__MPI)
    CALL mp_bcast(charge_mm, ionode_id, world_comm)
#endif
    !
    ! clear charge for QM atoms
    DO i = 1, nat_mm
       IF(tau_mask(i) .eq. -1)CYCLE
       charge_mm(i) = 0.0d0
    ENDDO
    !
    !IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: nat_all",nat_all
    !
    !IF (ionode) THEN
    !   WRITE(stdout,*)
    !   DO i = 1, nat_all
    !       WRITE(stdout,'(5X,A,3F10.6,2X,A,F10.6,2X,A,I2)') &
    !            'QMMM: tau_mm ',tau_mm(:,i),' charge_mm ',charge_mm(i),' QA ',tau_mask(i)
    !   END DO
    !END IF
    !
  END SUBROUTINE read_mm_charge
  !
  !
  SUBROUTINE read_mm_mask(socketfd)
    USE qmmm, ONLY : tau_mask, nat_mm
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM mask"
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( tau_mask, nat_mm, MDI_INT, socketfd, ierr )
    !
#if defined(__MPI)
    CALL mp_bcast(tau_mask, ionode_id, world_comm)
#endif
    !
  END SUBROUTINE read_mm_mask
  !
  !
  SUBROUTINE read_mm_coord(socketfd)
    USE qmmm, ONLY : nat_mm, tau_mm
    USE cell_base, ONLY : alat
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: ierr
    REAL(DP) :: buf(3*nat_mm)
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM coordinates"
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( buf, 3*nat_mm, MDI_DOUBLE, socketfd, ierr )
    tau_mm = RESHAPE(buf, (/3,nat_mm/))
    !
    tau_mm = tau_mm / alat ! internally positions are in alat
    !
#if defined(__MPI)
    CALL mp_bcast(tau_mm, ionode_id, world_comm)
#endif
    !
  END SUBROUTINE read_mm_coord
  !
  !
  SUBROUTINE read_types(socketfd)
    USE qmmm, ONLY : nat_mm, types
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM types"
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( types, nat_mm, MDI_INT, socketfd, ierr )
    !
#if defined(__MPI)
    CALL mp_bcast(types, ionode_id, world_comm)
#endif
    !
  END SUBROUTINE read_types
  !
  !
  SUBROUTINE read_mass(socketfd)
    USE qmmm,      ONLY : aradii, mass, rc_mm, nat_mm, types, ntypes
    USE cell_base, ONLY : alat
    USE constants, ONLY : bohr_radius_angs
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: ierr
    INTERFACE
       SUBROUTINE ec_fill_radii ( aradii, nat_mm, mass, types, ntypes, flag ) &
            BIND(C,name="ec_fill_radii")
         USE ISO_C_BINDING
         REAL(kind=c_double), INTENT(OUT) :: aradii(*)
         REAL(kind=c_double), INTENT(IN) :: mass(*)
         INTEGER(kind=c_int), INTENT(IN) :: types(*)
         INTEGER(kind=c_int), INTENT(IN) :: nat_mm, ntypes, flag
       END SUBROUTINE ec_fill_radii
    END INTERFACE
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM mass"
    !
    ! ... Read the dimensions of the MM cell
    !
    IF ( ionode ) CALL MDI_Recv( mass, ntypes+1, MDI_DOUBLE, socketfd, ierr )
    !
#if defined(__MPI)
    CALL mp_bcast(mass, ionode_id, world_comm)
#endif
    !
    ! do pre-forces work
    IF (ionode) THEN

        !CALL qmmm_center_molecule
        !CALL qmmm_minimum_image

        ! set atomic radii
        CALL ec_fill_radii( aradii, nat_mm, mass, types, ntypes, 1 )

    END IF
    CALL mp_bcast(aradii, ionode_id, world_comm)    
    rc_mm = aradii
    ! Convert radii to Bohr units
    rc_mm = rc_mm / (alat * bohr_radius_angs)
    !
  END SUBROUTINE read_mass
  !
  !
  SUBROUTINE read_aradii(socketfd)
    USE qmmm,      ONLY : aradii, rc_mm, nat_mm
    USE cell_base, ONLY : alat
    USE constants, ONLY : bohr_radius_angs
    INTEGER, INTENT(IN) :: socketfd
    INTEGER :: ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM aradii"
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( aradii, nat_mm, MDI_DOUBLE, socketfd, ierr )
    !
#if defined(__MPI)
    CALL mp_bcast(aradii, ionode_id, world_comm)
#endif
    !
    rc_mm = aradii
    ! Convert radii to Bohr units
    rc_mm = rc_mm / (alat * bohr_radius_angs)
    !
  END SUBROUTINE read_aradii

  !---------------------------------------------------------------------!
  ! communicate forces of the QM system to MM-master
  !
  SUBROUTINE write_ec_force( sockfd )
    USE qmmm, ONLY : qmmm_mode, nat_qm, force_qm
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sockfd
    INTEGER :: ierr
    REAL(DP) :: buf(3*nat_qm)
    

    IF (qmmm_mode .ne. 2) RETURN

    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Writing EC forces on QM atoms"
    !
    buf=RESHAPE(force_qm, (/ 3 * nat_qm /) ) * 0.5   ! force in a.u.
    !
    CALL MDI_Send( buf, 3*nat_qm, MDI_DOUBLE, sockfd, ierr )
    !

  END SUBROUTINE write_ec_force

  !---------------------------------------------------------------------!
  ! communicate forces of the QM system to MM-master
  !
  SUBROUTINE write_mm_force( sockfd, rho, nspin, dfftp )
    USE qmmm, ONLY : qmmm_mode, nat_mm, force_mm, qmmm_force_esf
    !
    USE fft_types,          ONLY : fft_type_descriptor
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sockfd
    INTEGER :: ierr
    REAL(DP) :: rho(:,:)
    INTEGER  :: nspin
    TYPE(fft_type_descriptor) :: dfftp
    REAL(DP) :: buf(3*nat_mm)

    IF (qmmm_mode .ne. 2) RETURN

    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Writing EC forces on MM atoms"

    !IF (ionode .and. (qmmm_verb > 0)) &
    !     WRITE(stdout,'(/,5X,A)') 'QMMM: compute EC forces'
    CALL qmmm_force_esf( rho, nspin, dfftp )

    !IF (qmmm_verb > 0) WRITE(stdout,'(5X,A)') 'QMMM: update forces'
    !
    !!!! Note, not used if ec_alg is false. Optimize excluding this send as well
    buf=RESHAPE(force_mm, (/ 3 * nat_mm /) ) * 0.5   ! force in a.u.
    CALL MDI_Send( buf, 3*nat_mm, MDI_DOUBLE, sockfd, ierr )

  END SUBROUTINE write_mm_force

  !---------------------------------------------------------------------!
  ! send the number of grid points used to represent the density
  !
  SUBROUTINE send_ndensity( sockfd, dfftp )
    USE fft_types,          ONLY : fft_type_descriptor
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sockfd
    TYPE(fft_type_descriptor) :: dfftp
    INTEGER :: ierr
    !
    INTEGER :: idx, i, j, k, j0, k0
    INTEGER :: ir
    INTEGER :: ngrid
    !
    ngrid = 0
    !
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    !
    DO ir = 1, dfftp%nnr
       idx = ir -1
       k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
       idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
       k   = k + k0
       IF ( k .GE. dfftp%nr3 ) CYCLE
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*j
       j   = j + j0
       IF ( j .GE. dfftp%nr2 ) CYCLE
       i   = idx
       IF ( i .GE. dfftp%nr1 ) CYCLE
       ngrid = ngrid + 1
    END DO
    !
    CALL mp_sum(ngrid, intra_pool_comm)
    !
    CALL MDI_Send( ngrid, 1, MDI_INT, sockfd, ierr )
    !
  END SUBROUTINE send_ndensity

  !---------------------------------------------------------------------!
  ! send the coordinates of the grid points used to represent the density
  !
  SUBROUTINE send_cdensity( sockfd, dfftp )
    USE fft_types,          ONLY : fft_type_descriptor
    USE cell_base,          ONLY : at, alat
    USE mp,                 ONLY : mp_gather, mp_gather
    USE mp_pools,           ONLY : intra_pool_comm, npool, nproc_pool, me_pool
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sockfd
    TYPE(fft_type_descriptor) :: dfftp
    INTEGER :: ierr
    !
    INTEGER :: idx, i, j, k, j0, k0
    INTEGER :: ir
    INTEGER :: ngrid, igrid, mygrid
    REAL(DP) :: s(3),r(3)
    REAL(DP), ALLOCATABLE :: cdensity(:), cdensity_all(:)
    !
    INTEGER :: mygrid_all(nproc_pool)
    INTEGER :: mygrid_displ(nproc_pool)
    !
    mygrid = 0
    !
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    !
    ! get the number of grid points
    DO ir = 1, dfftp%nnr
       idx = ir -1
       k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
       idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
       k   = k + k0
       IF ( k .GE. dfftp%nr3 ) CYCLE
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*j
       j   = j + j0
       IF ( j .GE. dfftp%nr2 ) CYCLE
       i   = idx
       IF ( i .GE. dfftp%nr1 ) CYCLE
       mygrid = mygrid + 1
    END DO
    !
    ! determine the number of grid points on each process
    CALL mp_gather(mygrid, mygrid_all, 0, intra_pool_comm)
    CALL mp_bcast(mygrid_all, 0, intra_pool_comm)
    DO ir = 1, nproc_pool
       if(ionode)WRITE(6,*)'   NGRID: ',ir,mygrid_all(ir)
    END DO
    !
    ngrid = mygrid
    CALL mp_sum(ngrid, intra_pool_comm)
    !
    ! construct the array of grid coordinates
    ALLOCATE( cdensity( 3*mygrid ) )
    IF ( me_pool .EQ. 0 ) ALLOCATE( cdensity_all( 3*ngrid ) )
    igrid = 1
    DO ir = 1, dfftp%nnr
       idx = ir -1
       k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
       idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
       k   = k + k0
       IF ( k .GE. dfftp%nr3 ) CYCLE
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*j
       j   = j + j0
       IF ( j .GE. dfftp%nr2 ) CYCLE
       i   = idx
       IF ( i .GE. dfftp%nr1 ) CYCLE
       !
       s(1) = DBLE(i)/DBLE(dfftp%nr1)
       s(2) = DBLE(j)/DBLE(dfftp%nr2)
       s(3) = DBLE(k)/DBLE(dfftp%nr3)
       !
       r=matmul(at,s)
       !
       cdensity( 3*(igrid-1) + 1 ) = r(1)*alat
       cdensity( 3*(igrid-1) + 2 ) = r(2)*alat
       cdensity( 3*(igrid-1) + 3 ) = r(3)*alat
       !
       igrid = igrid + 1
       !
    END DO
    !
    ! gather the grid coordinates across all processes
    !
    mygrid_displ(1) = 0
    DO ir = 2, nproc_pool
       mygrid_displ(ir) = mygrid_displ(ir-1) + 3*mygrid_all(ir-1)
    END DO
    !DO ir = 1, nproc_pool
    !   WRITE(6,*)'   DISPL: ',ir,mygrid_displ(ir)
    !END DO
    CALL mp_gather(cdensity, cdensity_all, 3*mygrid_all, mygrid_displ, 0, intra_pool_comm)
    !IF ( ionode ) THEN
    !   DO ir = 1, ngrid
    !      WRITE(6,*)'CDEN: ',ir,cdensity_all(ir)
    !   END DO
    !END IF
    !
    CALL MDI_Send( cdensity_all, 3*ngrid, MDI_DOUBLE, sockfd, ierr )
    DEALLOCATE( cdensity )
    IF ( me_pool .EQ. 0 ) DEALLOCATE( cdensity_all )
    !
  END SUBROUTINE send_cdensity


  SUBROUTINE send_density( sockfd, rho, nspin, dfftp )
    USE fft_types,          ONLY : fft_type_descriptor
    USE cell_base,          ONLY : at, alat
    USE mp,                 ONLY : mp_gather, mp_gather
    USE mp_pools,           ONLY : intra_pool_comm, npool, nproc_pool, me_pool
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sockfd
    REAL(DP) :: rho(:,:)
    INTEGER  :: nspin
    TYPE(fft_type_descriptor) :: dfftp
    INTEGER :: ierr
    !
    INTEGER :: idx, i, j, k, j0, k0
    INTEGER :: ir, ispin
    INTEGER :: ngrid, igrid, mygrid
    REAL(DP) :: s(3),r(3)
    REAL(DP), ALLOCATABLE :: density(:), density_all(:)
    !
    INTEGER :: mygrid_all(nproc_pool)
    INTEGER :: mygrid_displ(nproc_pool)
    !
    mygrid = 0
    !
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    !
    ! get the number of grid points
    DO ir = 1, dfftp%nnr
       idx = ir -1
       k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
       idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
       k   = k + k0
       IF ( k .GE. dfftp%nr3 ) CYCLE
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*j
       j   = j + j0
       IF ( j .GE. dfftp%nr2 ) CYCLE
       i   = idx
       IF ( i .GE. dfftp%nr1 ) CYCLE
       mygrid = mygrid + 1
    END DO
    !
    ! determine the number of grid points on each process
    CALL mp_gather(mygrid, mygrid_all, 0, intra_pool_comm)
    CALL mp_bcast(mygrid_all, 0, intra_pool_comm)
    DO ir = 1, nproc_pool
       if(ionode)WRITE(6,*)'   NGRID: ',ir,mygrid_all(ir)
    END DO
    !
    ngrid = mygrid
    CALL mp_sum(ngrid, intra_pool_comm)
    !
    ! construct the array of grid coordinates
    ALLOCATE( density( mygrid ) )
    IF ( me_pool .EQ. 0 ) ALLOCATE( density_all( ngrid ) )
    igrid = 1
    DO ir = 1, dfftp%nnr
       idx = ir -1
       k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
       idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
       k   = k + k0
       IF ( k .GE. dfftp%nr3 ) CYCLE
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*j
       j   = j + j0
       IF ( j .GE. dfftp%nr2 ) CYCLE
       i   = idx
       IF ( i .GE. dfftp%nr1 ) CYCLE
       !
       !s(1) = DBLE(i)/DBLE(dfftp%nr1)
       !s(2) = DBLE(j)/DBLE(dfftp%nr2)
       !s(3) = DBLE(k)/DBLE(dfftp%nr3)
       !
       !r=matmul(at,s)
       !
       !density( igrid ) = r(1)*alat
       density( igrid ) = 0.0_DP
       DO ispin = 1, nspin
          density( igrid ) = density( igrid ) + rho( ir, ispin )
       END DO
       !
       igrid = igrid + 1
       !
    END DO
    !
    ! gather the grid coordinates across all processes
    !
    mygrid_displ(1) = 0
    DO ir = 2, nproc_pool
       mygrid_displ(ir) = mygrid_displ(ir-1) + mygrid_all(ir-1)
    END DO
    CALL mp_gather(density, density_all, mygrid_all, mygrid_displ, 0, intra_pool_comm)
    !
    CALL MDI_Send( density_all, ngrid, MDI_DOUBLE, sockfd, ierr )
    DEALLOCATE( density )
    IF ( me_pool .EQ. 0 ) DEALLOCATE( density_all )
    !
  END SUBROUTINE send_density
  !>>>






  !
  SUBROUTINE recv_npotential( mdi_comm )
    USE mdi_engine, ONLY : npotential, potential
    !
    INTEGER :: mdi_comm, ierr
    !
    IF ( ionode ) WRITE(stdout,*)'In recv_npotential'
    !
    CALL MDI_Recv( npotential, 1, MDI_INT, mdi_comm, ierr )
    CALL mp_bcast( npotential, ionode_id, world_comm )
    !
    ! Allocate the potential array
    IF ( ALLOCATED( potential ) ) DEALLOCATE( potential )
    ALLOCATE( potential(npotential) )
    potential(:) = 0.0_DP
    !
  END SUBROUTINE recv_npotential
  !
  SUBROUTINE recv_potential( mdi_comm )
    USE mdi_engine, ONLY : npotential, potential
    !
    INTEGER :: mdi_comm, ierr
    !
    IF ( .not. ALLOCATED( potential ) ) THEN
       IF ( ionode ) WRITE(*,*) "MDI received the >POTENTIAL command before receiving the >NPOTENTIAL command"
       STOP 270
    END IF
    !
    CALL MDI_Recv( potential, npotential, MDI_DOUBLE, mdi_comm, ierr )
    potential(:) = potential(:) * 2.0_DP
    CALL mp_bcast( potential, ionode_id, world_comm )
    !
  END SUBROUTINE recv_potential
  !
  !
  SUBROUTINE read_qmmm_mode()
    USE qmmm,             ONLY : qmmm_mode, qmmm_initialization
    USE mp_global,        ONLY : intra_image_comm
    !
    INTEGER :: ierr
    !
    ! ... Reads the number of atoms
    !
    CALL MDI_Recv( qmmm_mode, 1, MDI_INT, socket, ierr )
    CALL mp_bcast( qmmm_mode, ionode_id, intra_image_comm )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read qmmm mode: ",qmmm_mode
    !
    CALL qmmm_initialization()
    !
  END SUBROUTINE read_qmmm_mode
  !
  !
  SUBROUTINE recv_nat_mm()
    USE qmmm, ONLY : set_mm_natoms
    !
    INTEGER :: natoms_in
    INTEGER :: ierr
    !
    ! ... Reads the number of mm atoms
    !
    CALL MDI_Recv( natoms_in, 1, MDI_INT, socket, ierr )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read mm natoms: ",natoms_in
    !
    CALL set_mm_natoms(natoms_in)
    !
  END SUBROUTINE recv_nat_mm
  !
  !






  !
  !
  SUBROUTINE set_replica_id()
    USE io_global,        ONLY : ionode, ionode_id
    USE mp_global,        ONLY : intra_image_comm
    USE mdi_engine,       ONLY : rid, rid_old, socket
    USE mdi,              ONLY : MDI_INT, MDI_Recv
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    ! ... Check if the replica id (rid) is the same as in the last run
    !
    CALL MDI_Recv( rid, 1, MDI_INT, socket, ierr )
    CALL mp_bcast( rid, ionode_id, intra_image_comm )
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Receiving replica", rid, rid_old
    IF ( rid .NE. rid_old .AND. .NOT. mdi_first_scf ) THEN
       !
       ! ... If a different replica reset the history
       !
       IF ( .NOT. mdi_first_scf) CALL reset_history_for_extrapolation()
    END IF
    !
    rid_old = rid
    !
  END SUBROUTINE set_replica_id
  !
  !
  SUBROUTINE set_nat()
    USE io_global,        ONLY : ionode, ionode_id
    USE mp_global,        ONLY : intra_image_comm
    USE ions_base,        ONLY : nat
    USE mdi_engine,       ONLY : socket
    USE mdi,              ONLY : MDI_INT, MDI_Recv
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    ! ... Reads the number of atoms
    !
    CALL MDI_Recv( nat, 1, MDI_INT, socket, ierr )
    CALL mp_bcast( nat, ionode_id, intra_image_comm )
    !
    !
    CALL allocate_nat_arrays()
    !
  END SUBROUTINE set_nat


















  SUBROUTINE mdi_execute_command(command, mdi_comm, ierr)
    USE io_global,        ONLY : ionode
    USE scf,              ONLY : rho
    USE lsda_mod,         ONLY : nspin
    USE fft_base,         ONLY : dfftp
    USE mdi_engine,       ONLY : scf_current, socket, mdi_exit_flag
    USE qmmm,             ONLY : qmmm_minimum_image, qmmm_center_molecule

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: command
    INTEGER, INTENT(IN)           :: mdi_comm
    INTEGER, INTENT(OUT)          :: ierr

    ! Local variables
    INTEGER, PARAMETER :: MSGLEN=12
    REAL*8, PARAMETER :: gvec_omega_tol=1.0D-1
    LOGICAL :: isinit=.false., hasdata=.false., exst
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER :: ccmd, i, info
    REAL*8 :: sigma(3,3), at_reset(3,3), dist_reset, ang_reset
    REAL *8 :: cellih(3,3), vir(3,3), pot
    REAL*8 :: dist_ang(6), dist_ang_reset(6)

    ierr = 0

    SELECT CASE ( trim( command ) )
    CASE( ">RID" )
       CALL set_replica_id()
       isinit=.true.
       !
    CASE( ">NATOMS" )
       scf_current = .false.
       CALL set_nat()
       !
    CASE( ">MM_NATOMS" )
       scf_current = .false.
       CALL recv_nat_mm()
       !
    CASE( ">NTYPES" )
       scf_current = .false.
       CALL read_ntypes()
       !
    CASE( ">CELL" )
       scf_current = .false.
       CALL read_cell()
       CALL update_cell()
       !
    CASE( "<CELL" )
       CALL send_cell()
       !
    CASE( ">COORDS" )
       scf_current = .false.
       CALL read_coordinates()
       !
    CASE( "<COORDS" )
       CALL send_coordinates()
       !
    CASE( ">QMMM_MODE" )
       scf_current = .false.
       CALL read_qmmm_mode()
       !
    CASE( ">MM_CELL" )
       scf_current = .false.
       CALL read_cell_mm()
       !
    CASE( ">MM_CHARGES" )
       scf_current = .false.
       CALL read_mm_charge(socket)
       !
    CASE( ">MM_MASK" )
       scf_current = .false.
       CALL read_mm_mask(socket)
       !
    CASE( ">MM_COORDS" )
       scf_current = .false.
       CALL read_mm_coord(socket)
       !
    CASE( ">MM_TYPES" )
       scf_current = .false.
       CALL read_types(socket)
       !
    CASE( ">MM_MASSES" )
       scf_current = .false.
       CALL read_mass(socket)
       !
    CASE( ">ARADII" )
       scf_current = .false.
       CALL read_aradii(socket)
       !
    CASE( "RECENTER" )
       scf_current = .false.
       CALL qmmm_center_molecule()
       CALL qmmm_minimum_image()
       !
    CASE( "SCF" )
       CALL run_scf()
       !
    CASE( "<ENERGY" )
       CALL write_energy()
       !
    CASE( "<PE" )
       CALL write_energy()
       !
    CASE( "<FORCES" )
       CALL write_forces()
       !
    CASE( "<EC_FORCES" )
       CALL write_ec_force(socket)
       !
    CASE( "<MM_FORCES" )
       CALL write_mm_force(socket, rho%of_r, nspin, dfftp)
       !
    CASE( "<NATOMS" )
       CALL send_natoms()
       !
    CASE( "<NDENSITY" )
       CALL send_ndensity(socket, dfftp)
       !
    CASE( "<CDENSITY" )
       CALL send_cdensity(socket, dfftp)
       !
    CASE( "<DENSITY" )
       IF ( .not. scf_current ) THEN
          CALL run_scf()
       END IF
       CALL send_density(socket, rho%of_r, nspin, dfftp)
       !
    CASE( "<CHARGES" )
       CALL send_charges()
       !
    CASE( ">NPOTENTIAL" )
       scf_current = .false.
       CALL recv_npotential(socket)
       !
    CASE( ">POTENTIAL" )
       scf_current = .false.
       CALL recv_potential(socket)
       !
    CASE( "<STRESS" )
       CALL send_stress()
       !
    CASE( "EXIT" )
       mdi_exit_flag = .true.
       !
    CASE DEFAULT
       IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Unrecognized command: ",trim(command)
       ierr = 130
    END SELECT
    !
  END SUBROUTINE MDI_EXECUTE_COMMAND



  
  !
  !
  SUBROUTINE allocate_nat_arrays()
    USE io_global,        ONLY : ionode
    USE mdi_engine,       ONLY : combuf
    USE ions_base,        ONLY : nat
    USE qmmm,             ONLY : set_qm_natoms
    !
    IMPLICIT NONE
    !
    ! ... Allocate the dummy array for the atoms coordinates
    !
    IF ( ALLOCATED( combuf ) ) DEALLOCATE( combuf )
    ALLOCATE( combuf( 3 * nat ) )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read number of atoms: ",nat
    !
    CALL set_qm_natoms(nat)
    !
  END SUBROUTINE allocate_nat_arrays
  !
  !
  SUBROUTINE read_ntypes()
    USE io_global,        ONLY : ionode
    USE mdi_engine,       ONLY : socket
    USE mdi,              ONLY : MDI_INT, MDI_Recv
    USE qmmm,             ONLY : set_ntypes
    !
    IMPLICIT NONE
    !
    INTEGER :: ntypes_in
    INTEGER :: ierr
    !
    ! ... Reads the number of atom types
    !
    CALL MDI_Recv( ntypes_in, 1, MDI_INT, socket, ierr )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read ntypes: ",ntypes_in
    !
    CALL set_ntypes(ntypes_in)
    !
  END SUBROUTINE read_ntypes
  !
  !
  SUBROUTINE read_cell()
    USE io_global,        ONLY : ionode, ionode_id
    USE mp_global,        ONLY : intra_image_comm
    USE cellmd,           ONLY : omega_old, at_old
    USE cell_base,        ONLY : alat, at, omega
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Recv
    USE mdi_engine,       ONLY : socket
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    REAL *8 :: cellh(3,3)
    REAL *8 :: mtxbuffer(9)
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading cell "
    !
    IF ( .NOT. mdi_first_scf) THEN
       at_old = at
       omega_old = omega
    END IF
    !
    ! ... First reads cell and the number of atoms
    !
    CALL MDI_Recv( mtxbuffer, 9, MDI_DOUBLE, socket, ierr )
    cellh = RESHAPE(mtxbuffer, (/3,3/))
    !
    ! ... Share the received data 
    !
    CALL mp_bcast( cellh,  ionode_id, intra_image_comm )
    !
    ! ... Convert the incoming configuration to the internal pwscf format
    !
    cellh  = TRANSPOSE(  cellh )                 ! row-major to column-major
    at = cellh / alat                            ! internally cell is in alat
    !
    cell_changed = .true.
    !
  END SUBROUTINE read_cell
  !
  !
  SUBROUTINE send_cell()
    USE cell_base,        ONLY : alat, at
    USE io_global,        ONLY : ionode
    !
    USE kinds,            ONLY : DP
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE mdi_engine,       ONLY : socket
    !
    IMPLICIT NONE
    !
    REAL(DP) :: cell_mdi(9)
    INTEGER :: ierr
    REAL *8 :: cellh(3,3)
    !
    cellh = at * alat
    cellh = TRANSPOSE( cellh )
    !
    cell_mdi(1:9) = RESHAPE( cellh, (/9/))
    !cell_mdi(10:12) = 0.0_DP
    !
    CALL MDI_Send( cell_mdi, 9, MDI_DOUBLE, socket, ierr)
    !
  END SUBROUTINE send_cell
  !
  !
  SUBROUTINE update_cell()
    USE cell_base,        ONLY : alat, at, omega, bg
    USE io_global,        ONLY : ionode
    !
    IMPLICIT NONE
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Updating cell "
    !
    ! ... Recompute cell data
    !
    CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
    CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
    !
    !! ... Check if the cell is changed too much and in that case reset the
    !! ... g-vectors
    !!
    !lgreset = ( ABS ( omega_reset - omega ) / omega .GT. gvec_omega_tol )
    !
  END SUBROUTINE update_cell
  !
  !
  SUBROUTINE read_cell_mm()
    USE io_global,        ONLY : ionode
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Recv
    USE mdi_engine,       ONLY : socket
    USE qmmm,             ONLY : set_cell_mm
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    REAL *8 :: mtxbuffer(9)
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM cell "
    !
    ! ... Read the dimensions of the MM cell
    !
    CALL MDI_Recv( mtxbuffer, 9, MDI_DOUBLE, socket, ierr )
    !
    CALL set_cell_mm(mtxbuffer)
    !
  END SUBROUTINE read_cell_mm
  !
  !
  SUBROUTINE read_coordinates()
    USE cell_base,        ONLY : alat
    USE io_global,        ONLY : ionode, ionode_id
    USE mp_global,        ONLY : intra_image_comm
    USE ions_base,        ONLY : tau
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Recv
    USE mdi_engine,       ONLY : socket, combuf
    USE ions_base,        ONLY : nat
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER :: i, ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading coordinates "
    !
    ! ... Read the atoms coordinates and share them
    !
    CALL MDI_Recv( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    CALL mp_bcast( combuf, ionode_id, intra_image_comm)
    !
    ! ... Convert the incoming configuration to the internal pwscf format
    !
    tau = RESHAPE( combuf, (/ 3 , nat /) )/alat  ! internally positions are in alat 
    !
  END SUBROUTINE read_coordinates
  !
  !
  SUBROUTINE send_coordinates()
    USE cell_base,        ONLY : alat
    USE io_global,        ONLY : ionode
    USE ions_base,        ONLY : tau
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE mdi_engine,       ONLY : socket, combuf
    USE ions_base,        ONLY : nat
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Sending coordinates "
    !
    ! ... Convert the internal pwscf format to the MDI standard
    !
    combuf = RESHAPE( tau, (/ 3*nat /) )*alat
    !
    ! ... Read the atoms coordinates and share them
    !
    CALL MDI_Send( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE send_coordinates
  !
  !
  SUBROUTINE send_charges()
    !
    USE ions_base,        ONLY : nat
    USE io_global,        ONLY : ionode
    USE kinds,            ONLY : DP
    USE ions_base,        ONLY : zv, ityp
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE mdi_engine,       ONLY : socket
    !
    IMPLICIT NONE
    !
    REAL(DP) :: charges(nat)
    INTEGER  :: i, ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Sending charges "
    !
    ! ... Send the charges of each ion
    !
    DO i=1, nat
       charges(i) = zv(ityp(i))
    END DO
    CALL MDI_Send( charges, nat, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE send_charges
  !NOTE:
  !ALSO NEED qm_charge, mm_charge_all, mm_coord_all, mm_mask_all, type, mass
  !
  !
  SUBROUTINE run_scf()
    USE cell_base,        ONLY : alat, at, omega, bg
    USE extrapolation,    ONLY : update_pot
    USE control_flags,    ONLY : conv_elec, treinit_gvecs
    USE mdi_engine,       ONLY : scf_current
    !
    IMPLICIT NONE
    REAL*8, PARAMETER :: gvec_omega_tol=1.0D-1
    !
    ! ... Recompute cell data
    !
    IF ( cell_changed ) THEN
       CALL update_cell()
       !
       ! ... If the cell is changes too much, reinitialize G-Vectors
       ! ... also extrapolation history must be reset
       ! ... If firststep, it will also be executed (omega_reset equals 0),
       ! ... to make sure we initialize G-vectors using positions from I-PI
       IF ( ((ABS( omega_reset - omega ) / omega) .GT. gvec_omega_tol) &
                                   .AND. (gvec_omega_tol .GE. 0.d0) ) THEN
          IF ( ionode ) THEN
             IF ( mdi_first_scf ) THEN
                WRITE(*,*) " @ DRIVER MODE: initialize G-vectors "
             ELSE
                WRITE(*,*) " @ DRIVER MODE: reinitialize G-vectors "
             END IF
          END IF
          CALL initialize_g_vectors()
          CALL reset_history_for_extrapolation()
          !
       ELSE
          !
          ! ... Update only atomic position and potential from the history
          ! ... if the cell did not change too much
          !
          IF (.NOT. mdi_first_scf) THEN
             IF ( treinit_gvecs ) THEN
                IF ( cell_changed ) CALL scale_h()
                CALL reset_gvectors ( )
             ELSE
                CALL update_pot()
                CALL hinit1()
             END IF
          END IF
       END IF
    ELSE
       IF (.NOT. mdi_first_scf ) THEN
           CALL update_pot()
           CALL hinit1()
       ENDIF
    END IF
    !
    ! ... Run an scf calculation
    !
    CALL electrons()
    IF ( .NOT. conv_elec ) THEN
       CALL punch( 'all' )
       CALL stop_run( conv_elec )
    ENDIF
    !
    ! Reset SCF-related flags
    !
    mdi_first_scf = .false.
    scf_current = .true.
    cell_changed = .false.
    !
  END SUBROUTINE run_scf
  !
  !
  SUBROUTINE write_energy()
    USE io_global,        ONLY : ionode
    USE ener,             ONLY : etot
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE mdi_engine,       ONLY : scf_current, socket
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    ! ... Run an SCF calculation
    !
    IF ( .not. scf_current ) THEN
       CALL run_scf()
    END IF
    !
    ! ... Write the total energy in a.u.
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Sending energy: ",0.5*etot
    CALL MDI_Send( 0.5*etot, 1, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE write_energy
  !
  !
  SUBROUTINE write_forces()
    USE io_global,        ONLY : ionode
    USE mdi_engine,       ONLY : mdi_forces, combuf, scf_current, socket
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE ions_base,        ONLY : nat
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    ! ... Run an SCF calculation
    !
    IF ( .not. scf_current ) THEN
       CALL run_scf()
    END IF
    !
    ! ... Compute forces
    !
    CALL forces()
    !
    ! ... Converts forces to a.u.
    ! ... (so go from Ry to Ha)
    !
    combuf=RESHAPE(mdi_forces, (/ 3 * nat /) ) * 0.5   ! force in a.u.
    !
    ! ... Write the forces
    !
    CALL MDI_Send( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE write_forces
  !
  !
  SUBROUTINE send_natoms()
    USE ions_base,        ONLY : nat
    USE io_global,        ONLY : ionode
    USE mdi,              ONLY : MDI_INT, MDI_Send
    USE mdi_engine,       ONLY : socket
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    ! ... Send the number of atoms in the system
    !
    CALL MDI_Send( nat, 1, MDI_INT, socket, ierr )
    !
  END SUBROUTINE send_natoms
  !
  !
  SUBROUTINE send_stress()
    USE cell_base,        ONLY : alat, at, omega
    USE io_global,        ONLY : ionode
    !
    USE kinds,            ONLY : DP
    USE mdi,              ONLY : MDI_DOUBLE, MDI_Send
    USE mdi_engine,       ONLY : scf_current, socket
    !
    IMPLICIT NONE
    !
    REAL(DP) :: stress_mdi(9)
    INTEGER :: ierr
    REAL*8 :: sigma(3,3), vir(3,3)
    !
    ! ... Run an SCF calculation
    !
    IF ( .not. scf_current ) THEN
       CALL run_scf()
    END IF
    !
    ! ... Compute stress
    !
    CALL stress()
    !
    vir=TRANSPOSE( sigma ) * omega * 0.5          ! virial in a.u & no omega scal.
    !
    stress_mdi(1:9) = RESHAPE( vir, (/9/))
    !
    CALL MDI_Send( stress_mdi, 9, MDI_DOUBLE, socket, ierr)
    !
  END SUBROUTINE send_stress
  !
  !
  SUBROUTINE initialize_g_vectors()
    USE io_global,        ONLY : ionode, ionode_id
    USE cellmd,           ONLY : omega_old, at_old
    USE mp_global,        ONLY : intra_image_comm
    USE cell_base,        ONLY : at, bg, omega
    USE mp,               ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    CALL clean_pw( .FALSE. )
    CALL init_run()
    !
    CALL mp_bcast( at,        ionode_id, intra_image_comm )
    CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
    CALL mp_bcast( omega,     ionode_id, intra_image_comm )
    CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
    CALL mp_bcast( bg,        ionode_id, intra_image_comm )
    !
    omega_reset = omega
    !<<<
    !
    !lgreset = .false.
    !>>>
    !
  END SUBROUTINE initialize_g_vectors
  !
  SUBROUTINE reset_history_for_extrapolation()
    USE io_global,        ONLY : ionode
    USE io_files,         ONLY : iunupdate, nd_nmbr, prefix, tmp_dir, postfix, &
                                 wfc_dir, delete_if_present, seqopn
    USE extrapolation,    ONLY : update_file
    !
    ! ... Resets history of wavefunction and rho as if the
    ! ... previous step was the first one in the calculation.
    ! ... To this end, files with rho and wfc from previous steps
    ! ... must be deleted, and iunupdate unit wiped. The latter
    ! ... is achieved by deleting the file and recreating it using
    ! ... update_file() routine.
    !
    IMPLICIT NONE
    LOGICAL :: exst
    !
    ! ... Delete history files, names correspond to the ones
    ! ... in the update_pot() routine.
    !
    CALL delete_if_present(TRIM( wfc_dir ) // TRIM( prefix ) // '.oldwfc' // nd_nmbr)
    CALL delete_if_present(TRIM( wfc_dir ) // TRIM( prefix ) // '.old2wfc' // nd_nmbr)
    IF ( ionode ) THEN
       CALL delete_if_present(TRIM( tmp_dir ) // TRIM( prefix ) // postfix // 'charge-density.old.dat')
       CALL delete_if_present(TRIM( tmp_dir ) // TRIM( prefix ) // postfix // 'charge-density.old2.dat')
       !
       ! ... The easiest way to wipe the iunupdate unit, is to delete it
       ! ... and run update_file(), which will recreate the file
       !
       CALL seqopn( iunupdate, 'update', 'FORMATTED', exst )
       CLOSE(UNIT=iunupdate, STATUS='DELETE')
    END IF
    !
    CALL update_file()
    !
  END SUBROUTINE




END MODULE MDI_IMPLEMENTATION
