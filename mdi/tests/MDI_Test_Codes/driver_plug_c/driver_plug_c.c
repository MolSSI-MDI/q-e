#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdi.h"


void mpi_error(const char* errormsg) {
  MPI_Abort(MPI_COMM_WORLD, 1);
}


int code_for_plugin_instance(void* mpi_comm_ptr, MDI_Comm mdi_comm, void* class_object) {
  MPI_Comm mpi_comm = *(MPI_Comm*) mpi_comm_ptr;
  int my_rank;
  MPI_Comm_rank(mpi_comm, &my_rank);

  // Determine the name of the engine
  char* engine_name = malloc( MDI_NAME_LENGTH * sizeof(char) );
  if ( MDI_Send_command("<NAME", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }
  if ( MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, mdi_comm) != 0 ) {
    mpi_error("MDI_Recv returned non-zero exit code.");
  }

  if ( my_rank == 0 ) {
    printf("Engine name: %s\n", engine_name);
  }
  free( engine_name );


  // Get the number of atoms
  int natoms;
  if ( MDI_Send_command("<NATOMS", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }
  if ( MDI_Recv(&natoms, 1, MDI_INT, mdi_comm) != 0 ) {
    mpi_error("MDI_Recv returned non-zero exit code.");
  }

  // Broadcast the number of atoms to all ranks within this plugin instance
  //MPI_Bcast(&natoms, 1, MPI_INT, 0, mpi_comm);
  if ( my_rank == 0 ) {
    printf("Engine natoms: %d\n", natoms);
  }

  /*
  // Get the nuclear coordinates
  double* coords = new double [3*natoms];
  if ( MDI_Send_command("<COORDS", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }
  if ( MDI_Recv(coords, 3*natoms, MDI_DOUBLE, mdi_comm) != 0 ) {
    mpi_error("MDI_Recv returned non-zero exit code.");
  }
  delete[] coords;
  */

  // Send the "EXIT" command to the engine
  if ( MDI_Send_command("EXIT", mdi_comm) != 0 ) {
    mpi_error("MDI_Send_command returned non-zero exit code.");
  }

  return 0;
}


int main(int argc, char **argv) {
  int ret;

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  // Initialize MDI
  ret = MDI_Init(&argc, &argv);

  // Confirm that MDI was initialized successfully
  int initialized_mdi;
  ret = MDI_Initialized(&initialized_mdi);

  // Get the correct MPI intra-communicator for this code
  ret = MDI_MPI_get_world_comm(&world_comm);

  // Number of ranks that will run the driver
  // This is the number of ranks that will NOT run plugin instances
  // The value of this variable is read from the command-line options
  int driver_nranks = 0;

  // Number of ranks running EACH plugin instance
  // The value of this variable is read from the command-line options
  int plugin_nranks = 1;

  // Name of the plugin to use
  // The value of this variable is read from the command-line options
  char* plugin_name = NULL;

  // Read through all the command line options
  int iarg = 1;
  while ( iarg < argc ) {

    if ( strcmp(argv[iarg],"-plugin_name") == 0 ) {

      // Ensure that the argument to the -plugin_name option was provided
      if ( argc-iarg < 2 ) {
	mpi_error("The -plugin_name argument was not provided.");
      }

      // Set driver_nranks
      plugin_name = argv[iarg+1];
      iarg += 2;

    }
    else {
      mpi_error("Unrecognized option.");
    }

  }

  // Verify the value of driver_nranks
  if ( driver_nranks < 0 ) {
    mpi_error("Invalid value for driver_nranks [0, inf).");
  }

  // Verify the value of plugin_nranks
  if ( plugin_nranks <= 0 ) {
    mpi_error("Invalid value for plugin_nranks (0, inf).");
  }

  // Verify that the value of driver_nranks and plugin_nranks is consistent with world_size
  int world_size;
  MPI_Comm_size(world_comm, &world_size);
  if ( (world_size - driver_nranks) % plugin_nranks != 0 ) {
    mpi_error("Invalid values for driver_nranks and plugin_nranks: world_size - driver_nranks must be divisible by plugin_nranks.");
  }

  // Verify the value of plugin_name
  if ( plugin_name == NULL ) {
    mpi_error("Plugin name was not provided.");
  }
  
  // Split world_comm into MPI intra-comms for the driver and each plugin
  MPI_Comm intra_comm;
  int my_rank, color, intra_rank;
  MPI_Comm_rank(world_comm, &my_rank);
  if ( my_rank < driver_nranks ) {
    color = 0;
  }
  else {
    color = ( ( my_rank - driver_nranks ) / plugin_nranks ) + 1;
  }
  MPI_Comm_split(world_comm, color, my_rank, &intra_comm);
  MPI_Comm_rank(intra_comm, &intra_rank);

  if ( color == 0 ) { // Driver intra-comm

    if (intra_rank == 0 ) {
      printf("I am the driver\n");
    }

  }
  else { // Engine instance intra-comm

    if ( intra_rank == 0 ) {
      printf("I am engine instance: \n");
    }

    // Initialize and run an instance of the engine library
    if ( MDI_Launch_plugin(plugin_name,
			   "-mdi \"-name MM -role ENGINE -method LINK\"",
			   &intra_comm,
			   code_for_plugin_instance,
			   NULL) != 0 ) {
      mpi_error("MDI_Launch_plugin returned non-zero exit code.");
    }
  }

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
