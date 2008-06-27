
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"

void outputBAC(CDS_UID_t distacnum, CDS_UID_t locacnum)
{

  GenericMesg outMesg;
  BacMesg bacMesg;
  char buffer[200];


  bacMesg.action = AS_ADD;
  bacMesg.elength = distacnum;
  bacMesg.ebac_id = locacnum;
  bacMesg.eseq_id = 0;
  bacMesg.iseq_id = 0;
  bacMesg.type = AS_ENDS;
  // reproducible timestamp for regression tests
  bacMesg.entry_time = CREATION_TIME;
  bacMesg.num_bactigs = 0;
  bacMesg.bactig_list = NULL;
  bacMesg.ilength = 0;

  sprintf(buffer,"BAC for Bac Ends " F_UID "\n", locacnum );
  bacMesg.source = buffer;

  outMesg.m = &bacMesg;
  outMesg.t = MESG_BAC;

  WriteProtoMesg_AS(stdout,&outMesg);
  fflush(NULL);
}


/* >>>> MAIN / TOP <<<< */

int main(int argc, char *argv[])
{


  if (argc != 3)
    {
      exit (1);
    }

  outputBAC(STR_TO_UID(argv[1], NULL, 10), STR_TO_UID(argv[2], NULL, 10));

  exit(0);
}
