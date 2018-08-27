/* Library for sending and receiving messagess through sockets */
//#include <cstdlib>
#include <stdio.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <errno.h>
#include <unistd.h>
#include "sockets.h"


/* Print an error message */
void error(char *msg)
{
  perror(msg);
  exit(1);
}



/* Send text through the socket */
int send_label(int socket, char *msg)
{
  int ret;
  int remaining;
  char buffer[BUFFER_SIZE];

  strcpy(buffer, msg);

  char *buf = (char*)&buffer;
  remaining = sizeof(buffer);
  do {
    ret = write(socket, buf, remaining);
    if (ret < 0) {
      error("Could not write to socket");
    }
    else if (ret == 0) {
      error("Wrote label message of size zero");
    }
    else {
      buf += ret;
      remaining -= ret;
    }
    
  }
  while (remaining > 0);

}



/* Read a label from the socket */
int read_label(int socket, char *buf)
{
  int ret;
  int remaining;

  remaining = BUFFER_SIZE;
  do {
    ret = read(socket, buf, remaining);
    if (ret < 0) {
      error("Could not read message");
    }
    else if (ret == 0) {
      error("Read label message of size zero");
    }
    else {
      buf += ret;
      remaining -= ret;
    }
  }
  while (remaining > 0);

}



/* Send an array through the socket */
int send_array(int socket, void *data, int size)
{
  int ret;
  int remaining;

  remaining = size;
  do {
    ret = write(socket, data, remaining);
    if (ret < 0) {
      error("Could not write to socket");
    }
    else if (ret == 0) {
      error("Wrote array message of size zero");
    }
    else {
      data += ret;
      remaining -= ret;
    }
    
  }
  while (remaining > 0);

}



/* Read an array from the socket */
int receive_array(int socket, void *data, int size)
{
  int ret;
  int remaining;

  remaining = size;
  do {
    ret = read(socket, data, remaining);
    if (ret < 0) {
      error("Could not read message");
    }
    else if (ret == 0) {
      error("Read array message of size zero");
    }
    else {
      data += ret;
      remaining -= ret;
    }
  }
  while (remaining > 0);

}
