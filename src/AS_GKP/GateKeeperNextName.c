
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
static char CM_ID[] = "$Id: GateKeeperNextName.c,v 1.1.1.1 2004-04-14 13:51:40 catmandew Exp $";

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

/*********************************************************************
 * Module:   GateKeeperNextName
 *
 *
 * Command Line Interface
 *      gateKeeperNextName <projName> <gkpStorePath>
 *
 * Writes the name of the next batch output of the gatekeeper to stdout
 *
 *    Programmer:  S. Kravitz
 *       Written:  Mar 2000
 *
 *********************************************************************/

#define EXE_NAME "GateKeeperNextName"
#define INCREMENT 1
#define ROLE "next"

void usage(void);

int  main(int argc, char * argv []){

  GateKeeperStore GkpStore;
  char *projName;
  char *gkpStorePath;
  int batchNum;
  int store_exists;

  if(argc != 3){
    fprintf(stderr,"*usage: " EXE_NAME " <projName> <gkpStorePath>\n");
    fprintf(stderr,"        Writes the name of the " ROLE " gatekeeper output to stdout\n");
    exit(1);
  }

  projName = argv[1];
  gkpStorePath = argv[2];
  
  fprintf(stderr,"* " EXE_NAME " proj:%s store:%s\n",
	  projName, gkpStorePath);

  
  InitGateKeeperStore(&GkpStore, gkpStorePath);
  store_exists = TestOpenGateKeeperStore(&GkpStore);

  if(store_exists != 1){
    fprintf(stderr,"* Store %s does not exist or is not complete (%d)...exiting\n",
	    gkpStorePath, store_exists);
    exit(1);
  }

  OpenGateKeeperStore(&GkpStore); // this is overkill, since we need one file...what the hey...
  batchNum = getNumGateKeeperBatchs(GkpStore.batStore);

  fprintf(stdout,"%s_%05d.inp\n", projName, batchNum + INCREMENT);
  CloseGateKeeperStore(&GkpStore);
  exit(0);
  return 0;
}

