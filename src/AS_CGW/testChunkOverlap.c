
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
/* 	$Id: testChunkOverlap.c,v 1.1.1.1 2004-04-14 13:51:06 catmandew Exp $	 */
#include "AS_global.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "AS_UTL_Var.h"
#include "AS_PER_SafeIO.h"
#include "ChunkStore_AS.h"
#include "AS_ALN_aligners.h"
#include "AS_global.h"


ChunkFileT *ChunkFile;


void ComputeOverlap(CDS_CID_t cidA, char *consensusA,
                    CDS_CID_t cidB, char *consensusB){

          CDS_COORD_t A_hang_estimate,A_hang_lower,A_hang_upper;
          InternalFragMesg AFR, BFR;
          OverlapMesg     *O;
          int where;
	  CDS_COORD_t beg, end;

          AFR.sequence   = consensusA;
          AFR.quality    = NULL;
          AFR.eaccession = 0;
          AFR.iaccession = cidA;
          BFR.sequence   = consensusB;
          BFR.quality    = NULL;
          BFR.eaccession = 0;
          BFR.iaccession = cidB;

	  end = strlen(consensusA);
	  beg = -end;

	  fprintf(stderr,"* Calling DP_Compare with %s sequences [" F_COORD "," F_COORD "] %c %c\n",
		  (!strcmp(AFR.sequence, BFR.sequence)?"EQUAL":"DIFFERENT"),
		  beg, end, *AFR.sequence, *BFR.sequence);

#define CGW_DP_ERATE .10
#define CGW_DP_THRESH 1e-6
#define CGW_DP_MINLEN 20

          O = DP_Compare_AS(&AFR,&BFR, 
                             beg, end, FALSE,
			    CGW_DP_ERATE,CGW_DP_THRESH, CGW_DP_MINLEN,
#if 0
			    AS_FIND_OVERLAP,  // Slower but more accurate?
#else
			    AS_FIND_ALIGN,    // Faster but less accurate?
#endif
			    &where);


	  if(O == NULL){
	    fprintf(stderr,"* NO Overlap found between %d and %d\n",
		    cidA, cidB);
	  }else{
	    fprintf(stderr,"* Overlap FOUND between %d and %d\n",
		    cidA, cidB);
	  }

}

int Debug = FALSE; 

int main(int argc, char *argv[]){
  int verbose = 0;
  char *tempPath;
  char chunkStorePath[FILENAME_MAX];

   { /* Parse the argument list using "man 3 getopt". */ 
     int ch,errflg=0;
     optarg = NULL;
     while (!errflg && ((ch = getopt(argc, argv, "vd")) != EOF))
       switch(ch) {
       case 'd':
	 Debug = TRUE;
	 break;
       case 'v':
	 verbose = TRUE;
	 break;
       case '?':
	 fprintf(stderr,"Unrecognized option -%c",optopt);
       default :
	 errflg++;
       }

     if((argc - optind != 1 ))
       {
	 fprintf (stderr, "USAGE:  testChunkOverlap [-v] "
		  "<ChunkStorePath> \n"
		  "Uses consensus data in <ChunkStorePath> to test aligner\n"
		  "Use -v to specify verbose output\n");
	 exit (EXIT_FAILURE);
       }

     tempPath = argv[optind++];
     ChunkFile = OpenChunkFile(tempPath);

     /* End of command line parsing */
   }

   {
     int i,j;
     char *consensusI, *qualityI;
     ChunkFragPosT *fragments;
     fprintf(stderr,"* Forward *\n");

     if(Debug){
       i = 3;
       GetChunk(ChunkFile, i, &consensusI, &qualityI, &fragments);
       ComputeOverlap(i,consensusI,i, consensusI);
     }

     for(i = 1; i < ChunkFile->indexLength; i++){
	 GetChunk(ChunkFile, i, &consensusI, &qualityI, &fragments);
	 ComputeOverlap(i,consensusI,i, consensusI);
     }

    }
}

