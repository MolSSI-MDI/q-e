PROGRAM QEMDI

  USE mdi,              ONLY : MDI_Init, MDI_MPI_get_world_comm, &
                               MDI_Register_command, MDI_Register_Node
  USE mdi_implementation, ONLY : MDI_Plugin_init_qemdi, mdi_listen

  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: set_command_line
  USE io_global, ONLY: stdout
  USE parallel_include, ONLY : MPI_THREAD_MULTIPLE

  IMPLICIT NONE

  ! MPI intra-communicator for this code
  INTEGER :: world_comm

  INTEGER  :: narg
  INTEGER :: iarg, ierr
  CHARACTER(len=1024) :: arg, mdi_options, input_file
  INTEGER             :: nim, npt, npl, nta, nbn, ndg
  INTEGER                :: retval
#if defined(_OPENMP)
  INTEGER :: PROVIDED
#endif


  WRITE(stdout,*)'HELLO WORLD!'

  ! Get the MDI command-line options
  narg = command_argument_count()
  mdi_options = ' '
  input_file = ' '

  DO iarg=1, narg
    CALL get_command_argument(iarg, arg)
    IF ( TRIM (arg) == '-mdi' .OR. TRIM (arg) == '--mdi' ) THEN
       IF ( mdi_options == ' ' ) THEN
          IF ( iarg+1 > narg ) THEN
             CALL infomsg('main','missing argument for -mdi command-line argument')
!             EXIT 1
          ELSE
             CALL get_command_argument(iarg+1, mdi_options)
          END IF
       ELSE
          CALL infomsg('main','duplicated -mdi command-line argument')
!          EXIT 1
       END IF
    END IF
    IF ( TRIM (arg) == '-in' ) THEN
       IF ( input_file == ' ' ) THEN
          IF ( iarg+1 > narg ) THEN
             CALL infomsg('main','missing argument for -in command-line argument')
!             EXIT 1
          ELSE
             CALL get_command_argument(iarg+1, input_file)
          END IF
       ELSE
          CALL infomsg('main','duplicated -in command-line argument')
!          EXIT 1
       END IF
    END IF
  END DO
  IF ( mdi_options == ' ' ) THEN
    CALL infomsg('main','missing -mdi command-line argument')
!    EXIT 1
  END IF
  IF ( input_file == ' ' ) THEN
    CALL infomsg('main','missing -in command-line argument')
!    EXIT 1
  END IF


#if defined(_OPENMP)
   CALL MPI_Init_thread(MPI_THREAD_MULTIPLE, PROVIDED, ierr)
#else
   CALL MPI_Init(ierr)
#endif


  ! Call MDI_Init

  !mdi_options = "-name QM -method TCP -role ENGINE -hostname localhost -port 8021"
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
  WRITE(stdout,*)'HERE 1'
  CALL set_command_line( nimage=nim, npool=npl, ntg=nta, &
       nband=nbn, ndiag=ndg )
  WRITE(stdout,*)'HERE 2'
  CALL mp_startup ( my_world_comm=world_comm )
  WRITE(stdout,*)'HERE 3'
  CALL environment_start ( 'PWSCF' )
  WRITE(stdout,*)'HERE 4'
  !
  !input_file = '/repo/tests/water/qe.in'
  CALL read_input_file ('PW+iPi', input_file )
  WRITE(stdout,*)'HERE 5'

  !
  ! Start a PW calculation, which will listen for MDI commands
  !
  retval = 0
  CALL mdi_listen( retval )
  WRITE(stdout,*)'AFTER MDI LISTEN'

END PROGRAM QEMDI
