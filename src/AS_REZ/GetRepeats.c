
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
   Module: get_repeats.c
   Description: displays multialignment and the repeat instances
                a fragment in the alignment comes from.
 *********************************************************************/

/**********************************************************************
 $Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_REZ/Attic/GetRepeats.c,v $
 $Revision: 1.4 $
***********************************************************************/

/*************************************************************************
   An interim translator from chunk messages to pseudo-contig messages
   for nuggetizer.
 *************************************************************************/
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"



int main (int argc, char *argv[]) {
  GenericMesg *pmesg;  
  ChunkMesg *unitig;
  IntContigMesg *icontig;
  MessageType imesgtype;
  MesgReader   reader;
  int i;
  int chunkId = -1;
  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "c")) != EOF))
      switch(ch) 
	{
	case 'c':
	  chunkId = atoi(argv[2]);
	  break;
	case '?':
	  fprintf(stderr,"Unrecognized option -%c\n",optopt);
	}
  }

  if( chunkId == -1 ){
    fprintf(stderr,"Missing chunk number. Provide it with -c flag.\n");
    exit(1);
  }

  //  icontig = (IntContigMesg *) malloc(sizeof(IntContigMesg));
  i=0;
  reader = InputFileType_AS(stdin);
  while ( reader(stdin,&pmesg) != EOF) {
    imesgtype = pmesg->t;
    switch(imesgtype){
    case MESG_CHK:
      {
	unitig = pmesg->m;
	if( /* unitig->source[5] == 'r' && */unitig-> iaccession == chunkId){
	 WriteProtoMesg_AS(stdout, pmesg);
	 //	 free((LayoutPos *)icontig->reads);
	 fflush(stdout);

       }
      } 
      i++;
      break;
    default:
      {
//       Swallowing all messages normally intended for CGW
//       WriteBinaryMesg_AS(stdout, pmesg);
//       fprintf(stderr,"Ignoring message... Type: %d\n",imesgtype);
      }
    }
  }
  fprintf(stderr,"getrep searched %d chunks...\n",i);
  exit(0);
}

