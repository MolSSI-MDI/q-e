!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE run_mdi
  !
  USE io_global,        ONLY : ionode, ionode_id
  USE cell_base,        ONLY : alat, at, omega, bg
  USE mdi,              ONLY : MDI_Send, MDI_Recv, MDI_Recv_Command, &
                               MDI_Accept_Communicator, &
                               MDI_CHAR, MDI_DOUBLE, MDI_INT, &
                               MDI_Set_execute_command_func
  use mdi_engine,       ONLY : socket, scf_current, get_mdi_options, &
                               firststep, combuf, &
                               recv_nat_mm
  !
  USE ISO_C_BINDING
  !
  IMPLICIT NONE
  SAVE
  !
  PRIVATE
  !
  INTEGER             :: nat
  INTEGER             :: rid, rid_old=-1

  PUBLIC :: mdi_listen

CONTAINS


  SUBROUTINE mdi_execute_command(command, mdi_comm, ierr)
    USE io_global,        ONLY : ionode
    USE scf,              ONLY : rho
    USE lsda_mod,         ONLY : nspin
    USE fft_base,         ONLY : dfftp

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
    INTEGER :: nat, rid, ccmd, i, info, rid_old=-1
    REAL*8 :: sigma(3,3), at_reset(3,3), dist_reset, ang_reset
    REAL *8 :: cellih(3,3), vir(3,3), pot
    REAL*8 :: dist_ang(6), dist_ang_reset(6)
    WRITE(6,*)'--------------------------------------'
    WRITE(6,*)'--------------------------------------'
    WRITE(6,*)'IN MDI_EXECUTE_COMMAND'
    WRITE(6,*)'COMMAND: ',trim(command)
    WRITE(6,*)'--------------------------------------'
    WRITE(6,*)'--------------------------------------'
    FLUSH(6)

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
     CASE( "EXIT" )
        ierr = 0
        RETURN
        !
     !>>>
     CASE DEFAULT
        IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Unrecognized command: ",trim(header)
        ierr = 130
        RETURN
     END SELECT

  END SUBROUTINE MDI_EXECUTE_COMMAND

  
  SUBROUTINE mdi_listen ( srvaddress, exit_status, mdi_options ) 
    !!
    !! Driver for IPI
    !!
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE parameters,       ONLY : ntypx, npk, lmaxx
    USE check_stop,       ONLY : check_stop_init
    USE mp_global,        ONLY : mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : gamma_only, conv_elec, istep, ethr, lscf, lmd
    USE cellmd,           ONLY : lmovecell
    USE force_mod,        ONLY : lforce, lstres
    USE ions_base,        ONLY : tau, nat_input => nat
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
    USE mdi_engine,       ONLY : is_mdi, recv_npotential, recv_potential, mdi_forces
    !USE command_line_options, ONLY : command_line
    !>>>
    !
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: exit_status
    !! Gives the exit status at the end
    CHARACTER(*), INTENT(IN) :: srvaddress
    CHARACTER(len=1024), OPTIONAL :: mdi_options
    !! Gives the socket address 
    !
    ! Local variables
    INTEGER, PARAMETER :: MSGLEN=12
    REAL*8, PARAMETER :: gvec_omega_tol=1.0D-1
    LOGICAL :: isinit=.false., hasdata=.false., exst
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER :: nat, rid, ccmd, i, info, rid_old=-1
    REAL*8 :: sigma(3,3), at_reset(3,3), dist_reset, ang_reset
    REAL *8 :: cellih(3,3), vir(3,3), pot
    REAL*8 :: dist_ang(6), dist_ang_reset(6)
    INTEGER :: ierr
    
    ! MDI Plugin callback function
    PROCEDURE(mdi_execute_command), POINTER :: mdi_execute_command_func => null()
    TYPE(C_PTR)                         :: class_obj
    mdi_execute_command_func => mdi_execute_command
    
    !>>>
    !----------------------------------------------------------------------------
    !
    lscf      = .true.
    lforce    = .true.
    lstres    = .true.
    lmd       = .true.
    lmovecell = .true.
    firststep = .true.
    !omega_reset = 0.d0
    !
    exit_status = 0
    IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
    !
    !<<<
    WRITE(6,*)'At start of run_driver'
    !>>>
    IF (ionode) CALL plugin_arguments()
    !<<<
    WRITE(6,*)'Calling plugin_arguments_bcast'
    !>>>
    CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
    !
    ! ... needs to come before iosys() so some input flags can be
    !     overridden without needing to write PWscf specific code.
    !
    ! ... convert to internal variables
    !
    !<<<
    WRITE(6,*)'Calling iosys'
    !>>>
    CALL iosys()
    !
    IF ( gamma_only ) WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
    !
    ! call to void routine for user defined / plugin patches initializations
    !
    !<<<
    WRITE(6,*)'Calling plugin_initialization'
    !>>>
    CALL plugin_initialization()
    !
    CALL check_stop_init()
    CALL setup()
    ! ... Initializations
    CALL init_run()
    !<<<
    !
    IF ( .not. PRESENT( mdi_options ) ) THEN
       mdi_options = get_mdi_options( )
    END IF
    IF ( .not. trim(mdi_options) == ' ' ) is_mdi = .true.
    WRITE(6,*)'MDI options: ',mdi_options
    WRITE(6,*)'is_mdi: ',is_mdi
    WRITE(6,*)'Calling create socket'
    IF (is_mdi) THEN
       nat = nat_input
       CALL allocate_nat_arrays()
    END IF
    !
    IF (ionode) THEN
       CALL MDI_Accept_Communicator( socket, ierr )
       CALL MDI_Set_execute_command_func(mdi_execute_command_func, class_obj, ierr)
    END IF
    WRITE(6,*)'Finished calling create socket'
    !>>>
    !
    driver_loop: DO
       !
       IF ( ionode ) CALL MDI_Recv_Command( header, socket, ierr )
       WRITE(6,*)'==============================='
       WRITE(6,*)'New command: ',trim(header)
       WRITE(6,*)'==============================='
       FLUSH(6)
       CALL mp_bcast( header, ionode_id, intra_image_comm )
       !
       IF ( ionode ) write(*,*) " @ DRIVER MODE: Message from server: ", trim( header )
       !
       
       
       call mdi_execute_command(header, socket, ierr)
       exit_status = ierr
       
     
       !
    END DO driver_loop
    !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
          & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
          & /,5X,'Max number of k-points (npk) = ',I6,&
          & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
    !
  END SUBROUTINE mdi_listen
  !
  !
  SUBROUTINE set_replica_id()
    USE mp_global,        ONLY : intra_image_comm
    !
    INTEGER :: ierr
    !
    ! ... Check if the replica id (rid) is the same as in the last run
    !
    IF ( ionode ) CALL MDI_Recv( rid, 1, MDI_INT, socket, ierr )
    CALL mp_bcast( rid, ionode_id, intra_image_comm )
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Receiving replica", rid, rid_old
    IF ( rid .NE. rid_old .AND. .NOT. firststep ) THEN
       !
       ! ... If a different replica reset the history
       !
       IF ( .NOT. firststep) CALL reset_history_for_extrapolation()
    END IF
    !
    rid_old = rid
    !
  END SUBROUTINE set_replica_id
  !
  !
  SUBROUTINE set_nat()
    USE mp_global,        ONLY : intra_image_comm
    USE ions_base,        ONLY : nat_input => nat
    !
    INTEGER :: ierr
    !
    ! ... Reads the number of atoms
    !
    IF ( ionode ) CALL MDI_Recv( nat, 1, MDI_INT, socket, ierr )
    CALL mp_bcast(    nat, ionode_id, intra_image_comm )
    nat_input = nat
    !
    !
    CALL allocate_nat_arrays()
    !
  END SUBROUTINE set_nat
  !
  !
  SUBROUTINE allocate_nat_arrays()
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
    !
    INTEGER :: ntypes_in
    INTEGER :: ierr
    !
    ! ... Reads the number of atom types
    !
    IF ( ionode ) CALL MDI_Recv( ntypes_in, 1, MDI_INT, socket, ierr )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read ntypes: ",ntypes_in
    !
    CALL set_ntypes(ntypes_in)
    !
  END SUBROUTINE read_ntypes
  !
  !
  SUBROUTINE read_cell()
    USE mp_global,        ONLY : intra_image_comm
    USE cellmd,           ONLY : omega_old, at_old
    !
    INTEGER :: ierr
    REAL *8 :: cellh(3,3)
    REAL *8 :: mtxbuffer(9)
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading cell "
    !
    IF ( .NOT. firststep) THEN
       at_old = at
       omega_old = omega
    END IF
    !
    ! ... First reads cell and the number of atoms
    !
    IF ( ionode ) CALL MDI_Recv( mtxbuffer, 9, MDI_DOUBLE, socket, ierr )
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Received cell "
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: mtxbuffer: ",mtxbuffer
    cellh = RESHAPE(mtxbuffer, (/3,3/))
    !
    ! ... Share the received data 
    !
    CALL mp_bcast( cellh,  ionode_id, intra_image_comm )
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: cellh: ",cellh
    !
    ! ... Convert the incoming configuration to the internal pwscf format
    !
    cellh  = TRANSPOSE(  cellh )                 ! row-major to column-major
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: alat ",alat
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: at before ",at
    at = cellh / alat                            ! internally cell is in alat
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: at after ",at
    !
  END SUBROUTINE read_cell
  !
  !
  SUBROUTINE send_cell()
    !
    USE kinds,            ONLY : DP
    !
    REAL(DP) :: cell_mdi(9)
    INTEGER :: ierr
    REAL *8 :: cellh(3,3)
    !
    IF ( ionode ) THEN
       !
       cellh = at * alat
       cellh = TRANSPOSE( cellh )
       !
       cell_mdi(1:9) = RESHAPE( cellh, (/9/))
       !cell_mdi(10:12) = 0.0_DP
       !
       CALL MDI_Send( cell_mdi, 9, MDI_DOUBLE, socket, ierr)
       !
    END IF
    !
  END SUBROUTINE send_cell
  !
  !
  SUBROUTINE update_cell()
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
    !
    INTEGER :: ierr
    REAL *8 :: mtxbuffer(9)
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading MM cell "
    !
    ! ... Read the dimensions of the MM cell
    !
    IF ( ionode ) CALL MDI_Recv( mtxbuffer, 9, MDI_DOUBLE, socket, ierr )
    !
    CALL set_cell_mm(mtxbuffer)
    !
  END SUBROUTINE read_cell_mm
  !
  !
  SUBROUTINE read_coordinates()
    USE mp_global,        ONLY : intra_image_comm
    USE ions_base,        ONLY : tau
    !
    INTEGER :: i, ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Reading coordinates "
    !
    ! ... Read the atoms coordinates and share them
    !
    IF ( ionode ) CALL MDI_Recv( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    CALL mp_bcast( combuf, ionode_id, intra_image_comm)
    !
    WRITE(6,*)" @ DRIVER MODE: input coordinates"
    DO i=1, nat
       WRITE(6,*)i,combuf((i-1)*3+1),combuf((i-1)*3+2),combuf((i-1)*3+3)
    END DO
    !
    ! ... Convert the incoming configuration to the internal pwscf format
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: old coords ",tau
    tau = RESHAPE( combuf, (/ 3 , nat /) )/alat  ! internally positions are in alat 
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: new coords ",tau
    !
  END SUBROUTINE read_coordinates
  !
  !
  SUBROUTINE send_coordinates()
    USE ions_base,        ONLY : tau
    !
    INTEGER :: ierr
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: Sending coordinates "
    !
    ! ... Convert the internal pwscf format to the MDI standard
    !
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: old coords ",tau
    combuf = RESHAPE( tau, (/ 3*nat /) )*alat
    IF ( ionode ) WRITE(*,*) " @ DRIVER MODE: new coords ",combuf
    !
    ! ... Read the atoms coordinates and share them
    !
    IF ( ionode ) CALL MDI_Send( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE send_coordinates
  !
  !
  SUBROUTINE send_charges()
    !
    USE kinds,            ONLY : DP
    USE ions_base,        ONLY : zv, ityp
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
    IF ( ionode ) CALL MDI_Send( charges, nat, MDI_DOUBLE, socket, ierr )
    !
    WRITE(6,*)" @ DRIVER MODE: sent charges: "
    DO i=1, nat
       WRITE(6,*)i,charges(i)
    END DO
    !
  END SUBROUTINE send_charges
  !NOTE:
  !ALSO NEED qm_charge, mm_charge_all, mm_coord_all, mm_mask_all, type, mass
  !
  !
  SUBROUTINE run_scf()
    USE control_flags,    ONLY : conv_elec
    !
    ! ... Initialize the G-Vectors when needed
    !
    !IF ( lgreset ) THEN
    !   !
    !   ! ... Reinitialize the G-Vectors if the cell is changed
    !   !
    !   CALL initialize_g_vectors()
    !   !<<<
    !   !lgreset = .false.
    !   !>>>
    !   !
    !ELSE
       !
       ! ... Update only atomic position and potential from the history
       ! ... if the cell did not change too much
       !
       CALL update_pot()
       CALL hinit1()
    !END IF
    !
    ! ... Run an scf calculation
    !
    CALL electrons()
    IF ( .NOT. conv_elec ) THEN
       CALL punch( 'all' )
       CALL stop_run( conv_elec )
    ENDIF
    scf_current = .true.
    !
  END SUBROUTINE run_scf
  !
  !
  SUBROUTINE write_energy()
    USE ener,             ONLY : etot
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
    IF ( ionode ) CALL MDI_Send( 0.5*etot, 1, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE write_energy
  !
  !
  SUBROUTINE write_forces()
    USE mdi_engine,       ONLY : mdi_forces
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
    IF ( ionode ) CALL MDI_Send( combuf, 3*nat, MDI_DOUBLE, socket, ierr )
    !
  END SUBROUTINE write_forces
  !
  !
  SUBROUTINE send_natoms()
    INTEGER :: ierr
    !
    ! ... Send the number of atoms in the system
    !
    IF ( ionode ) CALL MDI_Send( nat, 1, MDI_INT, socket, ierr )
    !
  END SUBROUTINE send_natoms
  !>>>
  !
  !
  SUBROUTINE initialize_g_vectors()
    USE cellmd,           ONLY : omega_old, at_old
    USE mp_global,        ONLY : intra_image_comm
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
    !omega_reset = omega
    !<<<
    !
    !lgreset = .false.
    !>>>
    !
  END SUBROUTINE initialize_g_vectors
  !
  SUBROUTINE reset_history_for_extrapolation()
    USE io_files,         ONLY : iunupdate, nd_nmbr, prefix, tmp_dir, postfix, &
                                 wfc_dir, delete_if_present, seqopn
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
  !
  !
END MODULE run_mdi
