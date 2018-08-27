!F90 ISO_C_BINGING wrapper for socket communication.

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and
!      port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

   MODULE F90MESSAGES
   USE ISO_C_BINDING
   
   IMPLICIT NONE

  INTERFACE
     SUBROUTINE initialize_client_c() BIND(C, name="initialize_client")
       USE ISO_C_BINDING
     END SUBROUTINE initialize_client_c
  END INTERFACE

    
!    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
!      USE ISO_C_BINDING
!    INTEGER(KIND=C_INT)                      :: psockfd
!    TYPE(C_PTR), VALUE                       :: pdata
!    INTEGER(KIND=C_INT)                      :: plen
!
!    END SUBROUTINE writebuffer_csocket       

!    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
!      USE ISO_C_BINDING
!    INTEGER(KIND=C_INT)                      :: psockfd
!    TYPE(C_PTR), VALUE                       :: pdata
!    INTEGER(KIND=C_INT)                      :: plen
!
!    END SUBROUTINE readbuffer_csocket   
!  END INTERFACE

   CONTAINS
   
   SUBROUTINE initialize_client()
     CALL initialize_client_c()
   END SUBROUTINE

  END MODULE
