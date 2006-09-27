
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

//  A collection of SDB related utility programs.

#include <stdlib.h> 
#include <stdio.h> 
#include <errno.h>
#include <assert.h>

#include "AS_global.h" 
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_CGW_dataTypes.h"
#include "AS_ALN_aligners.h"  //  Defines Overlap

//  Fetch two contigs from an SDB an output them to stdout in a format
//  suitable for readying by test-laligner_ca
//
void
dumpTwoContigs(char   *seqstore,
               int     version,
               int     contigiid1,
               int     contigiid2) {

  tSequenceDB *sequenceDB = OpenSequenceDB(seqstore, FALSE, version);	

  VA_TYPE(char) *seq1 = CreateVA_char(16384);
  VA_TYPE(char) *seq2 = CreateVA_char(16384);
  VA_TYPE(char) *qua  = CreateVA_char(16384);

  MultiAlignT   *ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, contigiid1, FALSE);
  MultiAlignT   *ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, contigiid2, FALSE);

  GetMultiAlignUngappedConsensus(ma1, seq1, qua); 
  GetMultiAlignUngappedConsensus(ma2, seq2, qua); 

  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n", contigiid1, contigiid1, Getchar(seq1,0));
  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n", contigiid2, contigiid2, Getchar(seq2,0));

  exit(0);
}



void
overlapTwoContigs(char   *seqstore,
                  int     version,
                  int     contigiid1,
                  int     contigiid2,
                  char    orientation) {

  switch (orientation) {
    case 'A':
    case 'N':
    case 'I':
    case 'O':
      break;
    default:
      fprintf(stderr,"invalid orientation.  Must be [AINO].\n");
      exit(1);
  }

  tSequenceDB *sequenceDB = OpenSequenceDB(seqstore, FALSE, version);	

  VA_TYPE(char) *seq1 = CreateVA_char(16384);
  VA_TYPE(char) *seq2 = CreateVA_char(16384);
  VA_TYPE(char) *qua  = CreateVA_char(16384);

  MultiAlignT   *ma1 =  LoadMultiAlignTFromSequenceDB(sequenceDB, contigiid1, FALSE);
  MultiAlignT   *ma2 =  LoadMultiAlignTFromSequenceDB(sequenceDB, contigiid2, FALSE);

  GetMultiAlignUngappedConsensus(ma1, seq1, qua); 
  GetMultiAlignUngappedConsensus(ma2, seq2, qua); 

  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n", contigiid1, contigiid1, Getchar(seq1,0));
  fprintf(stdout,"%d > Sequence %d\n%s\njunk\n", contigiid2, contigiid2, Getchar(seq2,0));

  int length1 = strlen(Getchar(ma1->consensus,0));
  int length2 = strlen(Getchar(ma2->consensus,0));
  
  // if the orientation is BA_AB or BA_BA, we need to reverse
  // complement the first contig
  //
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq(Getchar(seq1,0));
  
  // if the orientation is AB_BA or BA_BA, we need to set the flip
  // variable for the second contig
  //
  int flip = 0;
  if (orientation == AB_BA || orientation == BA_BA)
    flip = 1;
  
  // min_ahang and end are essentially bounds on the a-hang
  //
  Overlap *olap = DP_Compare(Getchar(seq1, 0),
                             Getchar(seq2, 0),
                             -length2, //  min_ahang
                             length1,  //  max_ahang
                             flip,
                             0.06,
                             CGW_DP_THRESH,
                             CGW_DP_MINLEN,
                             AS_FIND_ALIGN);

  Print_Overlap(stderr, Getchar(seq1, 0), Getchar(seq2, 0), olap);

  // return seq1 to its original state
  //
  if (orientation == BA_AB || orientation == BA_BA)
    Complement_Seq(Getchar(seq1,0));

  fprintf( stderr, "olap->begpos = %d\n", olap->begpos);
  fprintf( stderr, "olap->endpos = %d\n", olap->endpos);
  fprintf( stderr, "olap->length = %d\n", olap->length);
  
  exit(0);
}
  

void
dumpLength(char   *seqstore,
           int     version,
           int     partition) {

  tSequenceDBPartition *SDBPartition = openSequenceDBPartition(seqstore, version, partition);
  VA_TYPE(int32) *members = GetContentSequenceDBPartition(SDBPartition);

  MultiAlignT *ma;
  int i;

  for (i=0; i<GetNumint32s(members); i++) {
    CDS_CID_t    id = *Getint32(members,i);
    MultiAlignT *ma =  loadFromSequenceDBPartition(SDBPartition,id);

    fprintf(stderr, "id "F_CID" has "F_SIZE_T" fragments.\n", id, GetNumIntMultiPoss(ma->f_list));
  }

  exit(0);
}


void
viewMultiAlign(char   *seqstore,
               int     version,
               int     partition,
               int  showDots,
               int  showQual,
               int   isUnitig,  //  if not, isContig
               int   multialignID,
               char  *frgStore) {

  tSequenceDB *sequenceDB = OpenSequenceDB(seqstore, FALSE, version);

  //  XXXXXXX!  Finish me!

#if 0
  FragStoreHandle global_fragStore = openFragStore(frgFilename, "rb");

  VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(500);

  MultiAlignT *ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, multialign_id, is_unitig);

  MultiAlignT *newma;

  num_tigs = GetNumIntUnitigPoss(ma->u_list);

  for (i=0;i<num_tigs;i++) {
    IntUnitigPos *upos = GetIntUnitigPos(ma->u_list,i);
    IntElementPos pos;
    pos.type = AS_UNITIG;
    pos.ident = upos->ident;
    pos.position.bgn = upos->position.bgn;
    pos.position.end = upos->position.end;
    SetVA_IntElementPos(positions,i,&pos);
  }
    
  newma = MergeMultiAlignsFast_new(sequenceDB, NULLFRAGSTOREHANDLE, positions, 0, 1, NULL, NULL);
  //MultiAlignT2Array(ma, global_fragStore, NULLFRAGSTOREHANDLE, &depth, &multia, &id_array);
  //fprintf(stderr,"Converted into a character array\n");
  PrintMultiAlignT(stderr,newma,global_fragStore,NULL, NULLFRAGSTOREHANDLE,show_qv,dots,READSTRUCT_LATEST);
#endif



  //  Or we can dump
#if 0
  VA_TYPE(char) *ungappedCNS = CreateVA_char(256*1024);
  VA_TYPE(char) *ungappedQLT = CreateVA_char(256*1024);
  MultiAlignT *ma =  LoadMultiAlignTFromSequenceDB(sequenceDB, multialign_id, is_unitig);

  GetMultiAlignUngappedConsensus(ma, ungappedCNS, ungappedQLT);

  char *sequence=Getchar(ungappedCNS,0);
  if(reverse)
    SequenceComplement( sequence, NULL);

  len = strlen(sequence);

  fprintf(stdout,">%s%d len=%d\n%s\n",
          (is_unitig ? "unitig" : "contig"),
          multialign_id,
          len,
          sequence);
#endif

 
  exit(0);
}




void
testSDB_ReadInput(tSequenceDB *sequenceDB,
                  char        *iumInput,
                  int          check) {

  FILE *inputFile = fopen(iumInput, "r");
  if (errno)
    fprintf(stderr, "Can't open '%s': %s\n", iumInput, strerror(errno)), exit(1);
  
  MesgReader     reader = (MesgReader)InputFileType_AS(inputFile);
  int32          totalFrags = 1;

  GenericMesg   *pmesg;

  while ((EOF != reader(inputFile, &pmesg))) {
    if (pmesg->t == MESG_IUM) {
      IntUnitigMesg *ium_mesg = (IntUnitigMesg *)pmesg->m;
      MultiAlignT *ma = CreateMultiAlignTFromIUM(ium_mesg, totalFrags, FALSE);
      totalFrags += ium_mesg->num_frags;

      if(check) {
        MultiAlignT *loadma = LoadMultiAlignTFromSequenceDB(sequenceDB, ium_mesg->iaccession, TRUE);
        if (CompareMultiAlignT(ma, loadma)){
          fprintf(stderr, "MultiAlignment %d differ\n", ium_mesg->iaccession);
        }
        DeleteMultiAlignT(loadma);
        DeleteMultiAlignT(ma);
      }else{
        InsertMultiAlignTInSequenceDB(sequenceDB, ium_mesg->iaccession, TRUE, ma, FALSE);
      }
    }
  }
}


void
testSDB(char   *seqstore,
        int     version,
        int     partition,
        char   *iumInput) {

  tSequenceDB *sequenceDB = CreateSequenceDB(seqstore, 0, FALSE);

  testSDB_ReadInput(sequenceDB, iumInput, FALSE);
  SaveSequenceDB(sequenceDB);
  SaveSequenceDB(sequenceDB);
  DeleteSequenceDB(sequenceDB);
  sequenceDB = OpenSequenceDB(seqstore, FALSE, 1);
  testSDB_ReadInput(sequenceDB, iumInput, TRUE);
  testSDB_ReadInput(sequenceDB, iumInput, TRUE);
  ClearCacheSequenceDB(sequenceDB, TRUE);
  testSDB_ReadInput(sequenceDB, iumInput, TRUE);
  DeleteSequenceDB(sequenceDB);

  exit(0);
}




int
main(int argc, char **argv) {
  fprintf(stderr, "usage: %s [a ?useful? SDB utility]\n");
  exit(0);
}
