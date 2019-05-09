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
  USE kinds,            ONLY : DP
  USE mdi,              ONLY : MDI_Send, MDI_Recv, MDI_Recv_Command, &
                               MDI_Accept_Communicator, &
                               MDI_CHAR, MDI_DOUBLE, MDI_INT
  USE mp_global,        ONLY : mp_bcast, mp_sum, mp_barrier, intra_pool_comm
  USE mp_pools,         ONLY : me_pool
  USE mp_world,         ONLY : world_comm
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
  PUBLIC :: is_mdi
  PUBLIC :: recv_npotential, recv_potential
  PUBLIC :: mdi_add_potential
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
    !USE cell_base,          ONLY : alat, at, omega
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
    !INTEGER :: idx, i, j, k, j0, k0
    INTEGER :: ir
    INTEGER :: ngrid
    !
    !INTEGER :: i_mm, i_qm, ipol, ii_qm
    ! r_nn is the cutoff for the nearest neighbour
    !REAL(DP) :: s(3),r(3), dist, r_nn, fder
    !
    !REAL(DP) :: esfcontrib
    !
    !!!!
    !
    IF ( ionode ) WRITE(stdout,*)'In mdi_add_potential'
    !
    IF ( .not.is_mdi .or. npotential.eq.0 ) RETURN
    !
    ngrid = dfftp%nnr
    CALL mp_sum( ngrid, intra_pool_comm )
    IF ( ngrid .ne. npotential ) THEN
       WRITE(*,*)'QE currently only supports the >NPOTENTIAL command when NPOTENTIAL equals DFFTP%NNR'
       STOP 270
    END IF
    !
    DO ir = 1, dfftp%nnr
       IF ( ionode ) WRITE(stdout,*)'VLTOT: ',ir,vltot(ir),potential(ir)
       vltot(ir) = vltot(ir) + potential(ir + me_pool*dfftp%nnr)
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
END MODULE mdi_engine
