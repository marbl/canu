
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
/* $Id: null_scaffolder.c,v 1.3 2005-03-22 19:04:47 jason_miller Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "AS_global.h"

static void Rewind( FILE * fp, MesgReader reader )
{
  if( fp != NULL )
  {
    if( reader == ReadProtoMesg_AS )
      rewind( fp );
    else
      CDS_FSEEK( fp, (off_t) sizeof( int ), SEEK_SET );
  }

}

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  MesgReader   reader;
  IntAugFragMesg af_mesg;
  IntConConMesg contig;
  IntUnitigPos iup;
  IntContigPairs icp;
  IntScaffoldMesg scaffold;
  IntDegenerateScaffoldMesg dscaffold;
  IntUnitigMesg *unitig;
  IntMultiPos *f_list;
  int i,num_frags;
  int32 scaffid=1;
  float coverage_threshold=10.0;
  int pass=0;
  int in_iums=0;
  int nread;
  if (argc > 1) coverage_threshold = atof(argv[1]);
  reader = InputFileType_AS(stdin);

  while ( 1 ){
    nread = reader(stdin,&pmesg);
    if (nread == EOF && pass < 3 ) {
       in_iums = 0;
       pass++;
       Rewind(stdin,reader);
    } else if (nread == EOF) {
       break;
    } else if (pmesg->t == MESG_IUM) {
      in_iums=1;
      if (pass == 0) {
         f_list = ((IntUnitigMesg *)pmesg->m)->f_list;
         num_frags = ((IntUnitigMesg *)pmesg->m)->num_frags;
         pmesg->m = &af_mesg;
         pmesg->t = MESG_IAF;
         af_mesg.mate_status =  UNRESOLVED_MATE;
         af_mesg.chimeric = 0;
         af_mesg.clear_rng.bgn = -1;
         af_mesg.clear_rng.end = -1;
         for (i=0;i<num_frags;i++) {
           af_mesg.iaccession = f_list[i].ident;
           af_mesg.type = f_list[i].type;
           WriteProtoMesg_AS(stdout,pmesg);
         }
      // go through unitig f_list and generate a frag message for each
      } else if (pass==1) {
      // just write the unitig message on the second pass
         WriteProtoMesg_AS(stdout,pmesg);
      } else if (pass==2) {
         unitig = (IntUnitigMesg *) pmesg->m;
      // produce a dummy contig, promoting each unitig to contig.
         contig.iaccession = unitig->iaccession;
         contig.placed = ( unitig->coverage_stat >= coverage_threshold ) ? AS_PLACED: AS_UNPLACED;
         contig.length = unitig->length;
         contig.consensus = unitig->consensus;
         contig.quality = unitig->quality;
         contig.forced = unitig->forced;
         contig.num_pieces = unitig->num_frags;
         contig.num_unitigs =  1;
         contig.pieces = unitig->f_list;
         iup.type = (contig.placed)?AS_UNIQUE_UNITIG:AS_SINGLE_UNITIG;
         iup.ident = unitig->iaccession;
         iup.position.bgn = 0;
         iup.position.end = contig.length;
         iup.delta_length = 0;
         iup.delta = NULL;
         contig.unitigs = &iup;
         pmesg->t = MESG_ICM;
         pmesg->m = &contig;
         WriteProtoMesg_AS(stdout,pmesg);
      } else {
         unitig = (IntUnitigMesg *) pmesg->m;
      // produce a scaffold message, placing each D-unique in a singleton scaffold
        if ( unitig->coverage_stat >= coverage_threshold ) {
         icp.contig1 = unitig->iaccession;
         icp.contig2 = unitig->iaccession;
         icp.mean = 0;
         icp.stddev = 0;
         icp.orient = AB_AB;
         scaffold.iaccession = scaffid++;
         scaffold.num_contig_pairs = 0;
         scaffold.contig_pairs = &icp;
         pmesg->t = MESG_ISF;
         pmesg->m = &scaffold;
         WriteProtoMesg_AS(stdout,pmesg);
        } else {
         dscaffold.icontig = unitig->iaccession;
         pmesg->t = MESG_IDS;
         pmesg->m = &dscaffold;
         WriteProtoMesg_AS(stdout,pmesg);
        }
      }   
    } else {
      if ( pass == 3) {
        if ( in_iums ) {
          break;
        }
      } else { // pass < 3
        if ( !in_iums ) {  // pass through the initial bunch of messages
          if ( pass == 0 ) WriteProtoMesg_AS(stdout,pmesg);
        } else { // beyond the stream of unitigs, rewind and now output the unitigs
          pass++;
          in_iums=0;
          Rewind(stdin,reader);
        }
      }
    }
 }
 exit (0);
}
