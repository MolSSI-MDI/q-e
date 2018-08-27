#ifndef QMMM_SOCKETS
#define QMMM_SOCKETS

#include <stdio.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <errno.h>

#define BUFFER_SIZE 12

void error(char*);
int send_label(int, char*);
int read_label(int, char*);
int send_array(int, void*, int);
int receive_array(int, void*, int);


#endif
