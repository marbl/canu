
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

#include "SYS_UIDcommon.h"
#include "SYS_UIDerror.h"
#include "SYS_UIDserver.h"

static void SYS_UIDsignalHandlerFunction(int signal);

/*******************************************************************************

Description: main() for the UID server

*******************************************************************************/
cds_int32 main(cds_int32 argc, char** argv)
{
  int i;
  struct sigaction sigact;
  int startsignal = SIGHUP;
  int endsignal = SIGUSR2;
  int signalnum;

   SYS_UIDtype = UID_SERVER_TYPE;

   SYS_UIDparseOptions(argc, argv);

  /* first, must test for run mode and handle daemon configuration *
   * if necessary                                                  */

   if (argc >= 8 && argv[1][1] == 'r')
   {
      if (fork() > 0)           /* fork to establish daemon as session leader */
         exit(0);

      if (setsid() == -1)
	SYS_UIDerrorMsg("UID server: setsid ERROR");
   }

   /* this section taken straight from Steven's "UNIX Network Programming" */
 
   signal(SIGHUP, SIG_IGN); /* prepare to ignore terminal disconnects */

   if (fork() > 0) /* fork again for complete resetting of context */
       exit(0);

   chdir("/");

   umask(0);

   /* NOTE: shouldn't need this code b/c no file activity to this point in parent *
    *   for (i=0;i<64;i++)                                                        *
    *      close(i);                                                              *
    */

   /* past final fork - now set universal signal handling */

   /* initialze data structure for using the            *
    * sig_handler_function                              */
   memset(&sigact, 0, sizeof(struct sigaction));
   sigact.sa_handler = SYS_UIDsignalHandlerFunction;
   sigact.sa_flags = 0;

   /* next step is to register the signal handler with  *
    * each of the signals of interest                   */
   for (signalnum = startsignal; signalnum <= endsignal; signalnum++)
   {
     sigaction(signalnum, &sigact, NULL);
   }

   if (SYS_UIDserverInitialize(argc,argv) == UID_FAILS)  /* perform initialization */
   {
      SYS_UIDerrorMsg("Could not initialize UID Server\n");
      exit(1);
   }

   if (SYS_UIDserverStart() == UID_FAILS) // should block indefinitely unless error
   {
      SYS_UIDerrorMsg("UID Server failure\n");
      exit(1);
   }

   return(0);
}

/*******************************************************************************

Description: catches signals and handles debug messages

*******************************************************************************/
static void
SYS_UIDsignalHandlerFunction(int signal)
{
   if (SYS_UIDtype == UID_SERVER_TYPE)
     {
       switch(signal) {
       case SIGHUP:
         SYS_UIDerrorMsg("UID server: received SIGHUP signal");
         break;
       case SIGINT:
         SYS_UIDerrorMsg("UID server: received SIGINT signal");
         break;
       case SIGQUIT:
         SYS_UIDerrorMsg("UID server: received SIGQUIT signal");
         break;
       case SIGILL:
         SYS_UIDerrorMsg("UID server: received SIGILL signal");
         break;
       case SIGTRAP:
         SYS_UIDerrorMsg("UID server: received SIGTRAP signal");
         break;
       case SIGABRT:
         SYS_UIDerrorMsg("UID server: received SIGABRT signal");
         break;
       /*
       case SIGEMT:
         SYS_UIDerrorMsg("UID server: received SIGEMT signal");
         break;
       */
       case SIGFPE:
         SYS_UIDerrorMsg("UID server: received SIGFPE signal");
         break;
       case SIGKILL:
         SYS_UIDerrorMsg("UID server: received SIGKILL signal");
         break;
       case SIGBUS:
         SYS_UIDerrorMsg("UID server: received SIGBUS signal");
         break;
       case SIGSEGV:
         SYS_UIDerrorMsg("UID server: received SIGSEGV signal");
         break;
       case SIGSYS:
         SYS_UIDerrorMsg("UID server: received SIGSYS signal");
         break;
       case SIGPIPE:
         SYS_UIDerrorMsg("UID server: received SIGPIPE signal");
         break;
       case SIGALRM:
         SYS_UIDerrorMsg("UID server: received SIGALRM signal");
         break;
       case SIGTERM:
         SYS_UIDerrorMsg("UID server: received SIGTERM signal");
         break;
       case SIGURG:
         SYS_UIDerrorMsg("UID server: received SIGURG signal");
         break;
       case SIGSTOP:
         SYS_UIDerrorMsg("UID server: received SIGSTOP signal");
         break;
       case SIGTSTP:
         SYS_UIDerrorMsg("UID server: received SIGTSTP signal");
         break;
       case SIGCONT:
         SYS_UIDerrorMsg("UID server: received SIGCONT signal");
         break;
       case SIGCHLD:
         SYS_UIDerrorMsg("UID server: received SIGCHLD signal");
         break;
       case SIGTTIN:
         SYS_UIDerrorMsg("UID server: received SIGTTIN signal");
         break;
       case SIGTTOU:
         SYS_UIDerrorMsg("UID server: received SIGTTOU signal");
         break;
#ifdef SIGPOLL
       case SIGPOLL:
         SYS_UIDerrorMsg("UID server: received SIGPOLL signal");
         break;
#endif
       case SIGXCPU:
         SYS_UIDerrorMsg("UID server: received SIGXCPU signal");
         break;
       case SIGXFSZ:
         SYS_UIDerrorMsg("UID server: received SIGXFSZ signal");
         break;
       case SIGVTALRM:
         SYS_UIDerrorMsg("UID server: received SIGVTALRM signal");
         break;
       case SIGPROF:
         SYS_UIDerrorMsg("UID server: received SIGPROF signal");
         break;
       case SIGWINCH:
         SYS_UIDerrorMsg("UID server: received SIGWINCH signal");
         break;
       /*
       case SIGINFO:
         SYS_UIDerrorMsg("UID server: received SIGINFO signal");
         break;
       */
       case SIGUSR1:
         SYS_UIDerrorMsg("UID server: received SIGUSR1 signal");
         break;
       case SIGUSR2:
         SYS_UIDerrorMsg("UID server: received SIGUSR2 signal");
         break;
       default:
         sprintf(SYS_UIDerr_str,"UID server: received signal %d", signal);
         SYS_UIDerrorMsg(SYS_UIDerr_str);
       }
     } else {
       ; /* do nothing */
     }
}
