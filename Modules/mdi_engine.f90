!
! Copyright (C) 2013 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==-----------------------------------------------------------------------==!
MODULE mdi_engine
  !==---------------------------------------------------------------------==!
  USE cell_base,        ONLY : alat
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE ions_base,        ONLY : nat
  USE kinds,            ONLY : DP
  USE mdi,              ONLY : MDI_Send, MDI_Recv, MDI_Recv_Command, &
                               MDI_Accept_Communicator, &
                               MDI_CHAR, MDI_DOUBLE, MDI_INT
  USE mp_global,        ONLY : mp_bcast, mp_sum, mp_barrier, intra_pool_comm
  USE mp_pools,         ONLY : me_pool, nproc_pool
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_gather
  USE qmmm,             ONLY : qmmm_mode, qmmm_initialization, set_mm_natoms, &
       set_qm_natoms, set_ntypes, set_cell_mm, &
       read_mm_charge, read_mm_mask, read_mm_coord, &
       read_types, read_mass, write_ec_force, &
       write_mm_force, qmmm_center_molecule, &
       qmmm_minimum_image, read_aradii, send_ndensity, &
       send_cdensity, send_density
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! Indicates whether or not the code is running with MDI
  LOGICAL :: is_mdi = .false.
  !
  ! Return value of MDI functions
  INTEGER :: ierr
  !
  ! Number of points used to represent an external potential on a grid
  INTEGER :: npotential = 0
  !
  ! External potential at each grid point
  REAL(DP), ALLOCATABLE :: potential(:)
  !
  ! Forces used in response to the <FORCES command
  REAL(DP), ALLOCATABLE :: mdi_forces(:,:)
  !
  ! The mdi communicator
  INTEGER             :: socket
  !
  ! Flag whether the most recent SCF results are still valid
  LOGICAL :: scf_current=.false.
  !
  REAL*8, ALLOCATABLE :: combuf(:)
  LOGICAL :: firststep
  INTEGER             :: rid, rid_old=-1
  LOGICAL :: mdi_exit_flag = .false.
  !
  PUBLIC :: is_mdi, mdi_forces, socket, scf_current
  PUBLIC :: firststep, combuf, rid, rid_old
  PUBLIC :: mdi_exit_flag
  !
  PUBLIC :: recv_npotential, recv_potential
  PUBLIC :: mdi_add_potential, set_mdi_forces
  PUBLIC :: get_mdi_options
  PUBLIC :: read_qmmm_mode
  PUBLIC :: recv_nat_mm
  !
CONTAINS
  !
  SUBROUTINE recv_npotential( mdi_comm )
    !
    INTEGER :: mdi_comm
    !
    IF ( ionode ) WRITE(stdout,*)'In recv_npotential'
    !
    IF ( ionode ) CALL MDI_Recv( npotential, 1, MDI_INT, mdi_comm, ierr )
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
    !
    INTEGER :: mdi_comm
    !
    IF ( .not. ALLOCATED( potential ) ) THEN
       IF ( ionode ) WRITE(*,*) "MDI received the >POTENTIAL command before receiving the >NPOTENTIAL command"
       STOP 270
    END IF
    !
    IF ( ionode ) CALL MDI_Recv( potential, npotential, MDI_DOUBLE, mdi_comm, ierr )
    potential(:) = potential(:) * 2.0_DP
    CALL mp_bcast( potential, ionode_id, world_comm )
    !
  END SUBROUTINE recv_potential
  !
  SUBROUTINE mdi_add_potential( vltot, dfftp )

    !
    !   This routine adds an electrostatic field due to MM atoms to the 
    !   local potential.
    !
    USE cell_base,          ONLY : alat, at, omega
    !USE ions_base,          ONLY : zv, tau
    !USE constants,          ONLY : e2, eps8, bohr_radius_angs
    !USE io_global,          ONLY : stdout,ionode
    USE fft_types,          ONLY : fft_type_descriptor
    !
    !USE constraints_module, ONLY : pbc
    !
    IMPLICIT NONE
    !
    REAL(DP) :: vltot(:)
    TYPE(fft_type_descriptor) :: dfftp
    !
    ! local variables
    !
    INTEGER :: idx, i, j, k, j0, k0
    INTEGER :: ir, ipot
    INTEGER :: ngrid, mygrid, mygrid_displ
    INTEGER :: mygrid_all(nproc_pool)
    !
    !INTEGER :: i_mm, i_qm, ipol, ii_qm
    ! r_nn is the cutoff for the nearest neighbour
    REAL(DP) :: s(3),r(3), dist, r_nn, fder
    !
    !REAL(DP) :: esfcontrib
    !
    !!!!
    !
    IF ( ionode ) WRITE(stdout,*)'In mdi_add_potential'
    IF ( ionode ) WRITE(stdout,*)'VLTOT size: ',size(vltot)
    !
    IF ( .not.is_mdi .or. npotential.eq.0 ) RETURN
    !
    !ngrid = dfftp%nnr
    mygrid = 0
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
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
    CALL mp_gather(mygrid, mygrid_all, 0, intra_pool_comm)
    CALL mp_bcast(mygrid_all, 0, intra_pool_comm)
    mygrid_displ = 0
    DO ir = 0, me_pool - 1
       mygrid_displ = mygrid_displ + mygrid_all(ir+1)
    END DO
    ngrid = mygrid
    CALL mp_sum( ngrid, intra_pool_comm )
    IF ( ionode ) WRITE(stdout,*)'VLTOT NGRID: ',ngrid
    IF ( ngrid .ne. npotential ) THEN
       WRITE(*,*)'QE currently only supports the >NPOTENTIAL command when NPOTENTIAL equals DFFTP%NNR'
       WRITE(*,*)'   NPOTENTIAL, DFFTP%NNR: ',npotential,ngrid
       STOP 270
    END IF
    !
    !!!DO ir = 1, dfftp%nnr
    !!!   !IF ( ionode ) WRITE(stdout,*)'VLTOT: ',ir,vltot(ir),potential(ir)
    !!!   vltot(ir) = vltot(ir) + potential(ir + me_pool*dfftp%nnr)
    !!!END DO
    ipot = 0
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
       ipot = ipot + 1
       vltot(ir) = vltot(ir) + potential(ipot + mygrid_displ)
       !
    END DO
    !
    !!!!
    !
    ! if either the MS2 or EC aren't enabled, exit immediately
    !IF( qmmm_mode /= 2 ) RETURN
    !
    ! Index for parallel summation
    !
    !j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    !
    !r_nn = 50000.d0 ! cut-off for the nearest neighbour
    !
    !r(:) = 0.d0
    ! 
    !DO ir = 1, dfftp%nnr
    !   idx = ir -1
    !   k   = idx / (dfftp%nr1x * dfftp%my_nr2p )
    !   idx = idx - (dfftp%nr1x * dfftp%my_nr2p ) * k
    !   k   = k + k0
    !   IF ( k .GE. dfftp%nr3 ) CYCLE
    !   j   = idx / dfftp%nr1x
    !   idx = idx - dfftp%nr1x*j
    !   j   = j + j0
    !   IF ( j .GE. dfftp%nr2 ) CYCLE
    !   i   = idx
    !   IF ( i .GE. dfftp%nr1 ) CYCLE
    !   !
    !   s(1) = DBLE(i)/DBLE(dfftp%nr1)
    !   s(2) = DBLE(j)/DBLE(dfftp%nr2)
    !   s(3) = DBLE(k)/DBLE(dfftp%nr3)
    !   !
    !   r=matmul(at,s)
    !   ! 
    !   ! Clear the contribute (it's an accumulator)
    !   esfcontrib = 0.0D0
    !   !
    !   DO i_mm = 1, nat_mm 

    !      if(tau_mask(i_mm) .ne. -1)cycle ! only MM atoms contribute to ESF

    !      dist=sqrt((tau_mm(1, i_mm)-r(1))**2 + (tau_mm(2, i_mm)-r(2))**2 + (tau_mm(3, i_mm)-r(3))**2)
          !
    !      if(dist .LE. r_nn) then
    !          esfcontrib = esfcontrib - e2*charge_mm(i_mm)*(rc_mm(i_mm)**4 -dist**4)/(rc_mm(i_mm)**5 -dist**5) / alat
    !      end if
    !   ENDDO
       !
       ! Add the contribute
    !   vltot(ir) = vltot(ir) + esfcontrib
       !
    !END DO
    !
    !r(:) = 0.D0
    !force_qm = 0.D0

    !ii_qm = 1
    !DO i_qm = 1, nat_mm
    !   if(tau_mask(i_qm) .eq. -1)cycle
    !   DO i_mm = 1, nat_mm
    !      IF(tau_mask(i_mm) .ne. -1)CYCLE
    !      dist = sqrt((tau_mm(1, i_mm) - tau_mm(1, i_qm))**2 +   &
    !                  (tau_mm(2, i_mm) - tau_mm(2, i_qm))**2 +   &
    !                  (tau_mm(3, i_mm) - tau_mm(3, i_qm))**2)

    !            fder = ( 5.d0*(dist**4)*( rc_mm(i_mm)**4 - dist**4 ) -   &
    !                     4.d0*(dist**3)*( rc_mm(i_mm)**5 - dist**5 ) ) / & 
    !                        ( ( rc_mm(i_mm)**5 - dist**5 )**2 )
    !            DO ipol = 1,3
    !                  force_qm(ipol,ii_qm) = force_qm(ipol,ii_qm) -  &
    !                  e2*charge_mm(i_mm)*zv(tau_mask(i_qm)) *       &
    !                  fder*(tau_mm(ipol, i_qm)-tau_mm(ipol, i_mm))/dist
    !            ENDDO
    !   ENDDO
    !   ii_qm = ii_qm + 1
    !ENDDO

    !force_qm=force_qm/(alat**2)

    RETURN

  END SUBROUTINE mdi_add_potential
  !
  SUBROUTINE set_mdi_forces( forces, ipol )
    !
    REAL(DP), INTENT(IN) :: forces(:,:)
    INTEGER, INTENT(IN) :: ipol
    INTEGER :: iatom
    !
    IF (.not.ALLOCATED(mdi_forces)) ALLOCATE(mdi_forces(3,nat))
    !
    DO iatom = 1, nat
       mdi_forces(ipol,iatom) = forces(ipol,iatom)
    END DO
    !
  END SUBROUTINE set_mdi_forces
  !
  FUNCTION get_mdi_options ( ) RESULT ( options )
    ! 
    ! checks for the presence of a command-line option of the form
    ! -mdi "options" or --mdi "options";
    ! returns "options", used to run pw.x in driver mode.
    ! On input, "commmand_line" must contain the unprocessed part of the command
    ! line, on all processors, as returned after a call to "get_cammand_line"
    !
    USE command_line_options, ONLY : my_iargc, my_getarg
    IMPLICIT NONE
    !CHARACTER(LEN=*), INTENT(IN) :: command_line
    CHARACTER(LEN=1024) :: options
    !
    INTEGER  :: nargs, narg
    CHARACTER (len=1024) :: arg
    !
    nargs = command_argument_count()
    options = ' '
    !IF ( command_line == ' ' ) RETURN
    !
    !nargs = my_iargc ( command_line )
    !
    narg = 0
10  CONTINUE
    CALL get_command_argument(narg, arg)
    !CALL my_getarg ( command_line, narg, arg )
    IF ( TRIM (arg) == '-mdi' .OR. TRIM (arg) == '--mdi' ) THEN
       IF ( options == ' ' ) THEN
          narg = narg + 1
          IF ( narg > nargs ) THEN
             CALL infomsg('get_server_address','missing server IP in command line')
             RETURN
          ELSE
             CALL get_command_argument(narg, options)
             !CALL my_getarg ( command_line, narg, options )
          END IF
       ELSE
          CALL infomsg('get_server_address','duplicated server IP in command line')
       END IF
    END IF
    narg = narg + 1
    IF ( narg > nargs ) RETURN
    GO TO 10
    !
  END FUNCTION get_mdi_options
  !
  !
  SUBROUTINE read_qmmm_mode()
    USE qmmm,             ONLY : qmmm_mode
    USE mp_global,        ONLY : intra_image_comm
    !
    INTEGER :: ierr
    !
    ! ... Reads the number of atoms
    !
    IF ( ionode ) CALL MDI_Recv( qmmm_mode, 1, MDI_INT, socket, ierr )
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
    !
    INTEGER :: natoms_in
    INTEGER :: ierr
    !
    ! ... Reads the number of mm atoms
    !
    IF ( ionode ) CALL MDI_Recv( natoms_in, 1, MDI_INT, socket, ierr )
    !
    IF ( ionode ) write(*,*) " @ DRIVER MODE: Read mm natoms: ",natoms_in
    !
    CALL set_mm_natoms(natoms_in)
    !
  END SUBROUTINE recv_nat_mm
  !
  !







!
END MODULE mdi_engine
