
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

static char CM_ID[] = "$Id: LastFragInStoreOVL.c,v 1.5 2007-02-12 22:16:58 brianwalenz Exp $";

#include  <stdio.h>
#include "AS_global.h"
#include "AS_PER_gkpStore.h"

int
main(int argc, char **argv) {

  if (testOpenGateKeeperStore(argv[1], FALSE)) {
    GateKeeperStore  *g = openGateKeeperStore(argv[1], FALSE);
    fprintf(stdout, "Last frag in store is iid = "F_IID"\n", getLastElemFragStore(g));
    closeGateKeeperStore(g);
    exit(0);
  }

  if (strcmp(argv[1], "-h") == 0) {
    fprintf(stderr, "usage: %s [-h] some.gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Prints the internal ID of the last fragment in the fragment store\n");
    fprintf(stderr, "specified on the command line.\n");
  } else {
    fprintf(stderr, "Failed to open GateKeeperStore '%s'.\n", argv[1]);
  }

  exit(1);
}
