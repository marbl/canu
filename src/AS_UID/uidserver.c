
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: uidserver.c,v 1.1 2009-11-23 00:31:38 brianwalenz Exp $";

#include "uidserver_common.h"

#define ACTION_NONE  0
#define ACTION_INIT  1
#define ACTION_KILL  2
#define ACTION_RESU  3




static
void
initializeServer(void) {

  errno = 0;

  struct addrinfo       hints    = {0};
  struct addrinfo      *servinfo = NULL;
  struct addrinfo      *servloop = NULL;

  //  Create a connection for communication

  hints.ai_family    = AF_UNSPEC;
  hints.ai_socktype  = SOCK_STREAM;
  hints.ai_flags     = AI_PASSIVE; // use my IP address

  errno = getaddrinfo(NULL, serverPort, &hints, &servinfo);
  if (errno != 0) {
    fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(errno));
    exit(1);
  }

  //  Loop through all the results and bind to the first we can

  for (servloop = servinfo; servloop != NULL; servloop = servloop->ai_next) {
    serverSocket = socket(servloop->ai_family, servloop->ai_socktype, servloop->ai_protocol);
    if (errno) {
      perror("socket");
      continue;
    }

    bind(serverSocket, servloop->ai_addr, servloop->ai_addrlen);
    if (errno) {
      close(serverSocket);
      perror("bind");
      continue;
    }

    //  If we get here, we must have connected successfully
    break;
  }

  if (servloop == NULL) {
    // looped off the end of the list with no successful bind
    fprintf(stderr, "failed to bind socket\n");
    exit(1);
  }

  freeaddrinfo(servinfo); // all done with this structure

  //  Set some timeouts

  {
    struct timeval to;

    to.tv_sec = UID_SERVER_SEND_TIMEOUT;
    to.tv_usec = 0;
    if (setsockopt(serverSocket,
                   SOL_SOCKET,
                   SO_SNDTIMEO,
                   (const void *)&to,
                   sizeof(struct timeval)) < 0) {
      logError(LOG_ERR, "failed to setsockopt() for sends: %s\n", strerror(errno));
      exit(1);
    }
  }

  {
    struct timeval to;
    to.tv_sec = UID_SERVER_RECV_TIMEOUT;
    to.tv_usec = 0;
    if (setsockopt(serverSocket,
                   SOL_SOCKET,
                   SO_RCVTIMEO,
                   (const void *)&to,
                   sizeof(struct timeval)) < 0) {
      logError(LOG_ERR, "failed to setsockopt() for receives: %s\n", strerror(errno));
      exit(1);
    }
  }

  //  Active the connection

  {
    errno = 0;

    listen(serverSocket, SOMAXCONN/2);
    if (errno) {
      close(serverSocket);
      logError(LOG_ERR, "failed to listen(): %s\n", strerror(errno));
      exit(1);
    }
  }

  //  Last step, read the current highest UID number from the database file (the error is logged in
  //  readDatabase()).

  if (readDatabase() == 1)
    exit(1);
}




static
int
runServer(void) {
  int32   clientSocket = 0;

  errno = 0;

  //  Block until we get a connection

  {
    struct sockaddr_in   sockdata = {0};
    socklen_t            socksize = sizeof(struct sockaddr_in);

    clientSocket = accept(serverSocket, (struct sockaddr *)&sockdata, &socksize);
  }

  if (errno) {
    logError(LOG_ERR, "runServer failed to accept(): %s\n", strerror(errno));
    exit(1);
  }

  //  Read the client request

  UIDserverMessage  request  = {0};
  UIDserverMessage  response = {0};

  if (recvMessage(clientSocket, &request)) {
    logError(LOG_ERR, "runServer bogus message %d %d %d %d\n",
             request.message, request.bgnUID, request.endUID, request.numUID);
    return(1);
  }

  //  Process the request


  //  Shutdown.  Close connections.  Prune out UIDs the server has allocated but has not yet
  //  assigned.  Failing to prune isn't an error, we just lose a little bit of our UID space.
  //
  if (request.message == UIDserverMessage_SHUTDOWN) {
    logError(LOG_INFO, "port %s shutting down\n", serverPort);

    if (clientSocket)  close(clientSocket);
    if (serverSocket)  close(serverSocket);

#warning PRUNING DOES NOT SEEM TO WORK

    maxUID = curUID;

    if (writeDatabase()) {
      logError(LOG_WARNING, "port %d failed to prune database\n", serverPort);
      response.message = UIDserverMessage_ERROR;
    }

    fprintf(stderr, "SHUTDOWN!\n");
    return(0);
  }


  //  Client is asking for more UIDs.
  //
  if (request.message == UIDserverMessage_REQUEST) {

    if (curUID + request.numUID >= maxUID) {
      maxUID = curUID + request.numUID + incUID;

      if (writeDatabase()) {
        response.message = UIDserverMessage_ERROR;
      }
    }

    if (response.message == 0) {
      response.message = UIDserverMessage_GRANTED;
      response.bgnUID  = curUID;
      response.endUID  = curUID + request.numUID;
      response.numUID  = request.numUID;

      curUID += request.numUID;
    }

    goto sendMessageToClient;
  }


  //  Huh??  This is OUR response!
  if (request.message == UIDserverMessage_GRANTED) {
    response.message = UIDserverMessage_ERROR;
  }


  //  Unknown client request.
  response.message = UIDserverMessage_ERROR;

  //  Send the response

 sendMessageToClient:

  sendMessage(clientSocket, &response);
  close(clientSocket);

  return(1);
}



int32
main(int32 argc, char** argv) {
  int action = ACTION_NONE;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-i") == 0) {
      action        = ACTION_INIT;
      databaseName  = argv[++arg];
      maxUID        = strtoul(argv[++arg], (char**)NULL, 10);
    
    } else if (strcmp(argv[arg], "-k") == 0) {
      action        = ACTION_KILL;
      serverName = argv[++arg];
      serverPort = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      action       = ACTION_RESU;
      databaseName = argv[++arg];
      serverPort   = argv[++arg];

    } else {
      fprintf(stderr, "Unknown option %s\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) || (action == ACTION_NONE)) {
    fprintf(stderr, "usage: %s <command>\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "-i <filename> <initial_uid>\n", argv[0]);
    fprintf(stderr, "    Initialize a UID database file. Should only\n");
    fprintf(stderr, "    be done ONCE when a machine is configured.  Should NOT be re-run if\n");;
    fprintf(stderr, "    the machine is turned off or server re-booted or may risk\n");
    fprintf(stderr, "    corrupting UID uniqueness.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-r <filename> <start#> <size> <max_block> <update_freq> <port>\n", argv[0]);
    fprintf(stderr, "    Start the UID server from an initialized database file.\n");
    fprintf(stderr, "                                                                   \n");
    fprintf(stderr, "-k <hostname> <port>\n");
    fprintf(stderr, "    Kill a running UID server.\n");
    fprintf(stderr, "\n");

    exit(0);
  }



  //   Send a kill message to the server specified on the command line.
  if (action == ACTION_KILL) {
    UIDserverMessage        mesg     = {0};

    mesg.message = UIDserverMessage_SHUTDOWN;
    mesg.bgnUID  = 0;
    mesg.endUID  = 0;
    mesg.numUID  = 0;

    int32 scid = connectToServer(serverName, serverPort);
    sendMessage(scid, &mesg);
    close(scid);

    exit(0);
  }



  if (action == ACTION_INIT) {
    if (writeDatabase()) {
      fprintf(stderr, "Could not initialize database file '%s': %s\n", databaseName, strerror(errno));
      exit(1);
    }

    exit(0);
  }



  //  Fork off a server.  Most comes from Steven's "UNIX Network Programming".
#if 0
  if (fork() > 0)
    exit(0);

  if (setsid() == -1)
    fprintf(stderr, "setsid() failed: %s\n", strerror(errno));


  //  Ignore terminal disconnects.
  signal(SIGHUP, SIG_IGN);

  //  Forg again to complete resetting of context.
  if (fork() > 0)
    exit(0);
#endif

  //  Handle signals.

  //chdir("/");
  //umask(0);

  {
    struct sigaction sigact = {0};

    sigact.sa_handler = sigHandler;
    sigact.sa_flags   = 0;

    for (int sn = SIGHUP; sn <= SIGUSR2; sn++)
      sigaction(sn, &sigact, NULL);
  }

  initializeServer();

  while (runServer()) {
    //logError(LOG_INFO, "port %s starting again\n", serverPort);
  }

  logError(LOG_INFO, "port %s exiting\n", serverPort);

  exit(0);
}

