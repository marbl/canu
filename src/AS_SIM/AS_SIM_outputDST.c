
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


void outputDST(Distance_ID acnum, float32 mean, float32 delta)
{
  
  GenericMesg outMesg;
  DistanceMesg distanceMesg;
  float32 stddev = delta/3.0;  // Use delta = 3 sigma relationship

  // fprintf(stderr,"acnum = " F_U64 " mean = %g std =%g\n", acnum, mean, stddev);

  distanceMesg.action = AS_ADD;
  distanceMesg.eaccession = acnum;
  distanceMesg.mean = mean;
  distanceMesg.stddev = stddev;

  outMesg.m = &distanceMesg;
  outMesg.t = MESG_DST;

  WriteProtoMesg_AS(stdout,&outMesg);

}


/* >>>> MAIN / TOP <<<< */

int main(int argc, char *argv[])
{ 

  //      fprintf(stderr,"argc = %d  argv[1] = %s argv[2] = %s argv[3] = %s\n",
  //	      argc, argv[1], argv[2], argv[3]);

  if (argc != 4)
    { 
      exit (1);
    }

  outputDST(STR_TO_UID(argv[1], NULL, 10), atof(argv[2]), atof(argv[3]));

  exit(0);
}
