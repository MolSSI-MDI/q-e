#ifndef QMMM_MESSAGES
#define QMMM_MESSAGES

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include <errno.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <netdb.h>
#include "sockets.h"

//#define SOCKET_NAME "./9Lq7BNBnBycd6nxy.socket"
//#define BUFFER_SIZE 12

int initialize_client();

int driver_socket;
int qm_socket, qm_socket_in;
int mm_socket, mm_socket_in;
int mm_subset_socket, mm_subset_socket_in;
struct sockaddr_un qm_server, qm_client;
char buffer[BUFFER_SIZE];

int natoms = 3;
int num_qm = 2;
int num_mm = 2;
int ntypes = 2;

float boxlo0 = 0.0;
float boxlo1 = 0.0;
float boxlo2 = 0.0;
float boxhi0 = 0.0;
float boxhi1 = 0.0;
float boxhi2 = 0.0;
float cellxy = 0.0;
float cellxz = 0.0;
float cellyz = 0.0;

double *qm_coord;
double *qm_charge;
double *mm_charge_all;
double *mm_coord_all;
int *mm_mask_all;
int *type;
double *mass;

double *qm_force;
double *mm_force_all;
double *mm_force_on_qm_atoms;

#endif
