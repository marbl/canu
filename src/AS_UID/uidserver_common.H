#include "AS_global.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <syslog.h>

#include <netdb.h>

#include <netinet/in.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include <signal.h>
#include <errno.h>



char   *databaseName = NULL;
char   *serverName   = NULL;
char   *serverPort   = NULL;

int32   serverSocket = 0;

uint64  curUID       = 0;       //  Next UID to serve
uint64  maxUID       = 0;       //  UIDs we have checked out from the database
uint64  incUID       = 100000;  //  Number of UIDs to grab per server block


#define UID_CLIENT_SEND_TIMEOUT   20
#define UID_CLIENT_RECV_TIMEOUT   20
#define UID_SERVER_SEND_TIMEOUT   20
#define UID_SERVER_RECV_TIMEOUT   20


#define UIDserverMessage_UNDEFINED  0x00
#define UIDserverMessage_SHUTDOWN   0x01
#define UIDserverMessage_REQUEST    0x02
#define UIDserverMessage_GRANTED    0x03
#define UIDserverMessage_ERROR      0xff

struct UIDserverMessage {
  uint64    message;
  uint64    bgnUID;
  uint64    endUID;
  uint64    numUID;
};

////////////////////////////////////////
//
//  Logging and misc.
//
static
void
logError(int level, const char *fmt, ...) {
  va_list   ap;

  va_start(ap, fmt);

  //  In decreasing order of severity
  //
  //  LOG_EMERG   - panic, kaboom
  //  LOG_ALERT   - correct NOW
  //  LOG_CRIT    - immenient death
  //  LOG_ERR     - errors
  //  LOG_WARNING - warnings
  //  LOG_NOTICE  - information, correctable
  //  LOG_INFO    - information, fyi
  //  LOG_DEBUG   - debugging
  //
  //openlog("uid_server", LOG_CONS | LOG_PID, LOG_DAEMON);
  //vsyslog(level, fmt, ap);
  //closelog();

  vfprintf(stderr, fmt, ap);

  if (level == LOG_ERR)
    assert(0);
}


static
void
sigHandler(int signal) {
  logError(LOG_WARNING, "ignored signal %s (%d)\n", strsignal(signal), signal);
}





////////////////////////////////////////
//
//  Database functions
//
static
int32
writeDatabase(void) {

  errno=0;

  FILE *F = fopen(databaseName, "w");
  if (errno) {
    logError(LOG_ERR, "could not open database '%s' for writing: %s\n", databaseName, strerror(errno));
    return(1);
  }

  logError(LOG_INFO, "checked out up to UID %d\n", maxUID);

  fwrite(&maxUID, sizeof(uint64), 1, F);
  if (errno) {
    logError(LOG_ERR, "failed to write maxUID %d to database '%s': %s\n", maxUID, databaseName, strerror(errno));
    return(1);
  }

  fclose(F);
  if (errno) {
    logError(LOG_ERR, "failed to close database '%s' after writing maxUID %d: %s\n", databaseName, maxUID, strerror(errno));
    return(1);
  }

  return(0);
}


static
int32
readDatabase(void) {

  errno=0;

  FILE *F = fopen(databaseName, "r");
  if (errno) {
    logError(LOG_ERR, "could not open database '%s' for reading: %s\n", databaseName, strerror(errno));
    return(1);
  }

  fread(&curUID, sizeof(uint64), 1, F);
  if (errno) {
    logError(LOG_ERR, "failed to read curUID from database '%s': %s\n", databaseName, strerror(errno));
    return(1);
  }

  fclose(F);

  maxUID = curUID + incUID;

  logError(LOG_INFO, "restarting from UID %d\n", curUID);

  return(writeDatabase());
}


////////////////////////////////////////
//
//
//
static
int32
recvMessage(int32 fd, UIDserverMessage *mesg) {
  XDR     xdr;
  char    buffer[256];
  int32   nleft = sizeof(UIDserverMessage);
  int32   nxfer = 0;

  //  Read the message from the network.

  while (nleft > 0) {
    errno = 0;

    int32 n = recv(fd, buffer + nxfer, nleft, 0);

    if (n < 0) {
      logError(LOG_ERR, "error reading from network: %s\n", strerror(errno));
      return(1);
    }

    if (n == 0)
      break;

    nleft -= n;
    nxfer += n;
  }

  if (nleft > 0) {
    logError(LOG_ERR, "short read from network, read %d out of %d bytes\n", nxfer, sizeof(UIDserverMessage));
    return(1);
  }

  //  Decode the network message into a UIDserverMessage.

  xdrmem_create(&xdr, buffer, 256, XDR_DECODE);

  int32 success = 0;

  success += xdr_u_hyper(&xdr, &mesg->message);
  success += xdr_u_hyper(&xdr, &mesg->bgnUID);
  success += xdr_u_hyper(&xdr, &mesg->endUID);
  success += xdr_u_hyper(&xdr, &mesg->numUID);

  if (success != 4) {
    logError(LOG_ERR, "failed to unpack a message: %d %s\n", success, strerror(errno));
    exit(1);
  }

  xdr_destroy(&xdr);

  return(0);
}


static
int32
sendMessage(int32 fd, UIDserverMessage *mesg) {
  XDR     xdr;
  char    buffer[256];
  int32   nleft = sizeof(UIDserverMessage);
  int32   nxfer = 0;

  errno = 0;

  //  Encode the UIDserverMessage into a network message

  xdrmem_create(&xdr, buffer, 256, XDR_ENCODE);

  int32 success = 0;

  success += xdr_u_hyper(&xdr, &mesg->message);
  success += xdr_u_hyper(&xdr, &mesg->bgnUID);
  success += xdr_u_hyper(&xdr, &mesg->endUID);
  success += xdr_u_hyper(&xdr, &mesg->numUID);

  if (success != 4) {
    logError(LOG_ERR, "failed to pack response: %d %s\n", success, strerror(errno));
  }

  xdr_destroy(&xdr);

  //  Send the message to the netowrk.

  while (nleft > 0) {
    int32 n = send(fd, buffer + nxfer, nleft, 0);

    if (n <= 0) {
      logError(LOG_ERR, "error sending to network: %s\n", strerror(errno));
      return(1);
    }

    nleft -= n;
    nxfer += n;
  }

  if (nleft > 0) {
    logError(LOG_ERR, "short write to network, wrote %d out of %d bytes\n", nxfer, sizeof(UIDserverMessage));
    return(1);
  }

  return(0);
}




////////////////////////////////////////
//
//
//
int32
connectToServer(char *serverName, char *serverPort) {
  int32                   scid     = 0;
  struct addrinfo         hints;
  struct addrinfo        *servinfo = NULL;
  struct addrinfo        *servloop = NULL;

  memset(&hints, 0, sizeof hints);

  hints.ai_family   = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;

  //fprintf(stderr, "connectToServer %s port %s\n", serverName, serverPort);

  errno = getaddrinfo(serverName, serverPort, &hints, &servinfo);
  if (errno != 0) {
    fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(scid));
    exit(1);
  }

  // loop through all the results and connect to the first we can
  for (servloop = servinfo; servloop != NULL; servloop = servloop->ai_next) {
    scid = socket(servloop->ai_family, servloop->ai_socktype, servloop->ai_protocol);
    if (errno) {
      fprintf(stderr, "socket() failed: %s\n", strerror(errno));
      continue;
    }

    connect(scid, servloop->ai_addr, servloop->ai_addrlen);
    if (errno) {
      close(scid);
      fprintf(stderr, "connect() failed: %s\n", strerror(errno));
      continue;
    }

    // if we get here, we must have connected successfully
    break;
  }

  if (servloop == NULL) {
    // looped off the end of the list with no connection
    fprintf(stderr, "failed to connect\n");
    exit(1);
  }

  freeaddrinfo(servinfo); // all done with this structure

  return(scid);
}


