MODULE MDI_IMPLEMENTATION

  IMPLICIT NONE

  CONTAINS

  FUNCTION MDI_Plugin_init_qecouple() bind ( C, name="MDI_Plugin_init_qecouple" )
  USE mpi
  USE ISO_C_binding
  USE mdi,              ONLY : MDI_Init, MDI_Send, MDI_INT, MDI_CHAR, MDI_NAME_LENGTH, &
       MDI_Accept_communicator, MDI_Recv_command, MDI_Recv, &
       MDI_Set_execute_command_func, MDI_MPI_get_world_comm, MDI_DOUBLE, MDI_BYTE, &
       MDI_ENGINE, MDI_Get_role, MDI_Register_command, MDI_Register_node, &
       MDI_Register_callback, MDI_COMMAND_LENGTH, MDI_MPI_get_world_comm, &
       MDI_Plugin_get_argc, MDI_Plugin_get_arg

  ! MDI Communicator to the driver
  INTEGER :: comm

  ! MPI intra-communicator for this code
  INTEGER :: world_comm

  ! Flag to terminate MDI response function
  LOGICAL :: terminate_flag = .false.

    INTEGER :: MDI_Plugin_init_qecouple
    INTEGER :: ierr
    INTEGER :: argc
    INTEGER :: iarg
    CHARACTER(LEN=1024) :: option
    CHARACTER(LEN=1024) :: mdi_options
    LOGICAL :: mdi_options_found

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
             MDI_Plugin_init_qecouple = 1
             RETURN
          END IF
       END IF
    END DO
    IF ( .not. mdi_options_found ) THEN
       WRITE(6,*)'ERROR: -mdi option not provided'
       MDI_Plugin_init_qecouple = 1
       RETURN
    END IF

    ! Call MDI_Init
    !world_comm = MPI_COMM_WORLD
    CALL MDI_Init(mdi_options, ierr)

    ! Get the MPI intra-communicator over which this plugin will run
    CALL MDI_MPI_get_world_comm(world_comm, ierr);

    MDI_Plugin_init_qecouple = 0

  END FUNCTION MDI_Plugin_init_qecouple
END MODULE MDI_IMPLEMENTATION
