
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
#include <stdlib.h> 
#include <stdio.h> 
#include <assert.h>

#include "AS_global.h" 
#include "AS_ALN_aligners.h"
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_SDB_SequenceDB.h"
#include "Array_CNS.h"
#include "MultiAlignStore_CNS.h"
#include "ChunkOverlap_CGW.h"
#include "GreedyOverlapREZ.h"
#include "UtilsREZ.h"
#include "dpc_CNS.h"
#include "CommonREZ.h"
#include "GraphCGW_T.h"

Overlap* OverlapSequences1( char *seq1, char *seq2,
                            ChunkOrientationType orientation, 
                            int min_ahang, int max_ahang,
                            double erate, double thresh, int minlen,
                            CompareOptions what);

FILE *myerr=NULL;

#if 0
#define CGW_FUDGE_FACTOR (0.026)
#define CGW_DP_ERATE .10
#define CGW_DP_THRESH 1e-6
#define CGW_DP_MINLEN 20
#define CGW_DP_TRY_HARDER_MINLEN 20
#define CGW_DP_DESPERATION_MINLEN 10
#endif

int main(int argc, char *argv[]){
  
// args: seqStore revision_no  contig_id1 contig_id2 orientation
#if 0
  VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(500);
  FragStoreHandle global_fragStore = openFragStore(argv[3], "rb");
#endif
  tSequenceDB *sequenceDB = OpenSequenceDB(argv[1], FALSE, atoi(argv[2]));	
  int id1 = atoi(argv[3]);
  int id2 = atoi(argv[4]);
  char orientC = argv[5][0];
  Overlap * olap;
  int length1, length2;
  
  MultiAlignT *ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, id1, FALSE);
  MultiAlignT *ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, id2, FALSE);
  VA_TYPE(char) *seq1 = CreateVA_char(10000);
  VA_TYPE(char) *seq2 = CreateVA_char(10000);
  VA_TYPE(char) *qua = CreateVA_char(10000);
  ChunkOrientationType orient ;
  myerr = stderr;
  
  switch(orientC){
    case 'A':
    case 'N':
    case 'I':
    case 'O':
      orient = (ChunkOrientationType) orientC;
      break;
      
    default:
      fprintf(stderr,"* Illegal orientation %c...must be [AINO]\n", orientC);
      exit(1);
  }
  
  fprintf( myerr, "Aligning (%d,%d,%c)\n", id1, id2, orient);
  
  
  GetMultiAlignUngappedConsensus(ma1, seq1, qua); 
  GetMultiAlignUngappedConsensus(ma2, seq2, qua); 
  
  length1 = strlen(Getchar(ma1->consensus,0));
  length2 = strlen(Getchar(ma2->consensus,0));
  
  olap = OverlapSequences1( Getchar(seq1,0), Getchar(seq2,0), orient, 
                            -length2, length1,
                            0.06, CGW_DP_THRESH, CGW_DP_MINLEN, AS_FIND_ALIGN);
  
  fprintf( stderr, "olap->begpos = %d\n", olap->begpos);
  fprintf( stderr, "olap->endpos = %d\n", olap->endpos);
  fprintf( stderr, "olap->length = %d\n", olap->length);
  
  fprintf( stderr, "You can stop the debugger here if helpful.\n");
  return 0;
}

Overlap* OverlapSequences1( char *seq1, char *seq2,
                            ChunkOrientationType orientation, 
                            int min_ahang, int max_ahang,
                            double erate, double thresh, int minlen,
                            CompareOptions what)
{
  Overlap *omesg;
  int flip = 0;
  
  // if the orientation is BA_AB or BA_BA, we need to reverse complement the first contig
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip variable for the second contig
  if (orientation == AB_BA || orientation == BA_BA)
    flip = 1;
  
  // min_ahang and end are essentially bounds on the a-hang
  omesg = DP_Compare(seq1, seq2,
                     min_ahang, max_ahang, (int) flip,
                     erate, thresh, minlen,
                     what);
  
  Print_Overlap(myerr, seq1, seq2, omesg);
  
  // return seq1 to its original state
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq( seq1 );
  
  // omesg->begpos is the a-hang, omesg->endpos is the b-hang
  return omesg;
}
