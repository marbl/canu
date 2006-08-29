
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

/*********************************************************************
 * Writes the name of the current / next batch output of the
 * gatekeeper to stdout
 *
 *    Programmer:  S. Kravitz, B. Walenz
 *       Written:  Mar 2000,   Aug 2006
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"

int
main(int argc, char **argv){
  GateKeeperStore GkpStore;
  int             increment = 0;
  char           *projName = NULL;
  char           *gkpStorePath = NULL;
  int             batchNum  = 0;

  if ((argc != 4) && (argc != 3)) {
    fprintf(stderr,"*usage: %s [-c | -n] <projName> <gkpStorePath>\n");
    fprintf(stderr,"        -c   Writes the name of the current gatekeeper output to stdout\n");
    fprintf(stderr,"        -n   Writes the name of the next gatekeeper output to stdout\n");
    exit(1);
  }
  if        (strcmp(argv[1], "-c") == 0) {
    increment = 0;
  } else if (strcmp(argv[1], "-n") == 0) {
    increment = 1;
  }
  projName     = argv[argc-2];
  gkpStorePath = argv[argc-1];

  InitGateKeeperStore(&GkpStore, gkpStorePath);

  if (1 != TestOpenGateKeeperStore(&GkpStore))
    fprintf(stderr, "* ERROR: Store %s does not exist or is not complete.\n",
            gkpStorePath), exit(1);

  OpenGateKeeperStore(&GkpStore);
  batchNum = getNumGateKeeperBatchs(GkpStore.batStore);
  CloseGateKeeperStore(&GkpStore);

  fprintf(stdout, "%s_%05d.inp\n", projName, batchNum + increment);

  return(0);
}

