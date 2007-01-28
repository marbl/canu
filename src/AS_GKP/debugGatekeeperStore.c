
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
static char CM_ID[] = "$Id: debugGatekeeperStore.c,v 1.7 2007-01-28 21:52:24 brianwalenz Exp $";

/* Dump the gatekeeper stores for debug */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

int  main(int argc, char * argv [])
{
  int i;
  char * gkpStorePath = NULL;
  GateKeeperStore gkpStore;

  /**************** Process Command Line Arguments *********************/
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "g:")) != EOF))
    {
      switch(ch)
      {
        case 'g':
          gkpStorePath = optarg;
          break;
        case '?':
          fprintf(stderr,"Unrecognized option -%c",optopt);
        default :
          errflg++;
      }
    }

    if(gkpStorePath == NULL)
    {
      fprintf (stderr, "USAGE: %s [-g gatekeeperStorePath]\n", argv[0]);
      exit (EXIT_FAILURE);
    }
  }
   

  InitGateKeeperStore(&gkpStore, gkpStorePath);
  OpenReadOnlyGateKeeperStore(&gkpStore);
  {
    volatile PHashValue_AS value;
    volatile GateKeeperBatchRecord batch;
    volatile GateKeeperFragmentRecord frag1, frag2;
    volatile GateKeeperLocaleRecord locale;
    volatile GateKeeperSequenceRecord seq;
    volatile GateKeeperBactigRecord bactig;
    volatile GateKeeperDistanceRecord dist;
    volatile GateKeeperLinkRecord link;
    volatile StoreStat stat;
    
    fprintf(stderr, "Press Ctrl-C to interrupt in debugger!\n");
    while(TRUE)
    {
    }
  }
}
