
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
/* $Id: ofg2asm.c,v 1.1 2005-11-17 21:54:49 catmandew Exp $   */

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"

int main(int argc, char ** argv)
{
  FILE * fpIn;
  char fnIn[2048];
  FILE * fpOut;
  char fnOut[2048];
  MesgReader reader;
  GenericMesg * gen;
  OFGMesg * ofg;
  AugFragMesg afg;
  InternalDistMesg * idt;
  SnapMateDistMesg mdi;

  if( argc != 2 )
  {
    fprintf( stderr, "Usage: %s assemblyPrefix\n", argv[0] );
    return 1;
  }
  sprintf(fnIn, "%s.ofg", argv[1]);
  sprintf(fnOut, "%s.asm", argv[1]);

  if((fpIn = fopen(fnIn, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open file %s for reading.\n", fnIn);
    return 1;
  }
  reader = InputFileType_AS(fpIn);

  if((fpOut = fopen(fnOut, "w")) == NULL)
  {
    fprintf(stderr, "Failed to open file %s for writing.\n", fnOut);
    return 2;
  }

  // initialize unchanging attributes of afg messages
  afg.screened = NULL;
  afg.mate_status = UNRESOLVED_MATE;
  afg.chimeric = 0;
  afg.chaff = 1;

  // initialize unchanging attributes of mdi messages
  mdi.min = CDS_COORD_MIN;
  mdi.max = CDS_COORD_MAX;
  mdi.num_buckets = 0;
  mdi.histogram = NULL;

  // do messages
  while(reader(fpIn, &gen) != EOF)
  {
    int skipIt = 0;
    switch(gen->t)
    {
      case MESG_OFG:
        ofg = (OFGMesg *) gen->m;
        // copy ofg to afg
        afg.eaccession = ofg->eaccession;
        afg.iaccession = ofg->iaccession;
        afg.clear_rng.bgn = ofg->clear_rng.bgn;
        afg.clear_rng.end = ofg->clear_rng.end;
        gen->t = MESG_AFG;
        gen->m = &afg;
        break;
      case MESG_IDT:
        idt = (InternalDistMesg *) gen->m;
        // copy idt to mdi
        mdi.erefines = idt->eaccession;
        mdi.irefines = idt->iaccession;
        mdi.mean = idt->mean;
        mdi.stddev = idt->stddev;
        gen->t = MESG_MDI;
        gen->m = &mdi;
        break;
      case MESG_ADT:
        // do nothing
        break;
      default:
        skipIt = 1;
        break;
    }
    if(!skipIt)
      WriteProtoMesg_AS(fpOut, gen);
  }
  fclose(fpIn);
  fclose(fpOut);

  return 0;
}
