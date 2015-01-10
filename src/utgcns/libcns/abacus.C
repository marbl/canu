
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

static char *rcsid = "$Id$";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
//#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"









abacus::abacus() {

  //  This is called in ResetStores -- which is called before any
  //  consensus work is done.
  //
  //InitializeAlphTable(void)

  ALPHABET[0] = '-';  // These were lowercase, why?
  ALPHABET[1] = 'A';
  ALPHABET[2] = 'C';
  ALPHABET[3] = 'G';
  ALPHABET[4] = 'T';
  ALPHABET[5] = 'N';

  RALPHABET[ 0] = '-';
  RALPHABET[ 1] = 'A';
  RALPHABET[ 2] = 'C';
  RALPHABET[ 3] = 'G';
   RALPHABET[ 4] = 'T';
   RALPHABET[ 5] = 'N';
   RALPHABET[ 6] = 'a';  //  -A
   RALPHABET[ 7] = 'c';  //  -C
   RALPHABET[ 8] = 'g';  //  -G
   RALPHABET[ 9] = 't';  //  -T
   RALPHABET[10] = 'M';  //  AC
   RALPHABET[11] = 'R';  //  AG
   RALPHABET[12] = 'W';  //  AT
   RALPHABET[13] = 'S';  //  CG
   RALPHABET[14] = 'Y';  //  CT
   RALPHABET[15] = 'K';  //  GT
   RALPHABET[16] = 'm';  //  -AC
   RALPHABET[17] = 'r';  //  -AG
   RALPHABET[18] = 'w';  //  -AT
   RALPHABET[19] = 's';  //  -CG
   RALPHABET[20] = 'y';  //  -CT
   RALPHABET[21] = 'k';  //  -GT
   RALPHABET[22] = 'V';  //  ACG
   RALPHABET[23] = 'H';  //  ACT
   RALPHABET[24] = 'D';  //  AGT
   RALPHABET[25] = 'B';  //  CGT
   RALPHABET[26] = 'v';  //  -ACG
   RALPHABET[27] = 'h';  //  -ACT
   RALPHABET[28] = 'd';  //  -AGT
   RALPHABET[29] = 'b';  //  -CGT
   RALPHABET[30] = 'X';  //  ACGT
   RALPHABET[31] = 'x';  //  -ACGT

   RALPHABETC[ 0] = '-';
   RALPHABETC[ 1] = 'T';
   RALPHABETC[ 2] = 'G';
   RALPHABETC[ 3] = 'C';
   RALPHABETC[ 4] = 'A';
   RALPHABETC[ 5] = 'N';
   RALPHABETC[ 6] = 't';  //  -A
   RALPHABETC[ 7] = 'g';  //  -C
   RALPHABETC[ 8] = 'c';  //  -G
   RALPHABETC[ 9] = 'a';  //  -T
   RALPHABETC[10] = 'K';  //  AC
   RALPHABETC[11] = 'Y';  //  AG
   RALPHABETC[12] = 'W';  //  AT
   RALPHABETC[13] = 'S';  //  CG
   RALPHABETC[14] = 'R';  //  CT
   RALPHABETC[15] = 'M';  //  GT
   RALPHABETC[16] = 'k';  //  -AC
   RALPHABETC[17] = 'y';  //  -AG
   RALPHABETC[18] = 'w';  //  -AT
   RALPHABETC[19] = 's';  //  -CG
   RALPHABETC[20] = 'r';  //  -CT
   RALPHABETC[21] = 'm';  //  -GT
   RALPHABETC[22] = 'B';  //  ACG
   RALPHABETC[23] = 'D';  //  ACT
   RALPHABETC[24] = 'H';  //  AGT
   RALPHABETC[25] = 'V';  //  CGT
   RALPHABETC[26] = 'b';  //  -ACG
   RALPHABETC[27] = 'd';  //  -ACT
   RALPHABETC[28] = 'h';  //  -AGT
   RALPHABETC[29] = 'v';  //  -CGT
   RALPHABETC[30] = 'X';  //  ACGT
   RALPHABETC[31] = 'x';  //  -ACGT

   TAU_MISMATCH = 1.0 / (5.0 - 1.0);

   AMASK[0] = 013607700741;  //  -
   AMASK[1] = 015670707042;  //  a
   AMASK[2] = 016733131104;  //  c
   AMASK[3] = 017355252210;  //  g
   AMASK[4] = 017566464420;  //  t

  for (int32 i=0; i<256; i++)
    RINDEX[i] = 31;

  for (int32 i=0; i<CNS_NP; i++)
    RINDEX[(int)RALPHABET[i]] = i;

  RINDEX[(int)'n'] = RINDEX[(int)'N'];  //  Used in baseCount

  for (int32 i=0, qv=CNS_MIN_QV; i<CNS_MAX_QV - CNS_MIN_QV + 1; i++, qv++) {
    EPROB[i]= log(TAU_MISMATCH * pow(10, -qv/10.0));
    PROB[i] = log(1.0 - pow(10, -qv/10.0));
  }



};


abacus::~abacus() {
};















////////////////////////////////////////
//
//  Basic MANode operations
//

//external
int
abacus::GetMANodeConsensus(int32 mid, VA_TYPE(char) *sequence, VA_TYPE(char) *quality) {
}






////////////////////////////////////////
//
//  Bead Manipulation
//






////////////////////////////////////////
//
//  Columns
//




//external -- unused, but looks handy
#if 0
void
abacus::ShowColumn(abBeadID cid) {
  Column *column = GetColumn(columnStore,cid);
  Bead *call;
  Bead *bead;
  FragType type;
  UnitigType utype;
  ColumnBeadIterator ci;
  abBeadID bid;

  CreateColumnBeadIterator(cid,&ci);

  call = GetBead(beadStore,column->call);
  fprintf(stderr,"\nstore_index: %-20d ( prev: "F_BEADIDX" next: "F_BEADIDX")\n",column->lid,column->prev,column->next);
  fprintf(stderr,"ma_index:    %-20d\n",column->ma_index);
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"composition:\n");
  while ( (bid = NextColumnBead(&ci)).isValid() ) {
    bead = GetBead(beadStore,bid);
    type = GetFragment(fragmentStore,bead->frag_index)->type;
    utype = GetFragment(fragmentStore,bead->frag_index)->utype;
    fprintf(stderr,"             %c /%c (%10d) <-- %d iid:%d cid:%d UDLR:%d %d %d %d type:%c utype:%c\n",
            *Getchar(sequenceStore,bead->soffset),
            *Getchar(qualityStore,bead->soffset),
            bid,
            bead->frag_index,
            GetFragment(fragmentStore,bead->frag_index)->iid,
            bead->column_index,
            bead->up,
            bead->down,
            bead->prev,
            bead->next,
            type,
            utype);

    assert(bead->column_index == cid);
  }
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"call:        %c /%c\n",toupper(*Getchar(sequenceStore,call->soffset)),*Getchar(qualityStore,call->soffset));
}
#endif


////////////////////////////////////////
//
//  Data Management
//
void
abacus::ResetStores(int32 num_bases, int32 num_frags, int32 num_columns) {

  if (fragmentStore == NULL) {
    InitializeAlphTable();

    fragmentStore      = CreateVA_Fragment(num_frags);
    fragment_indices   = CreateVA_int32(num_frags);
    abacus_indices     = CreateVA_int32(50000);
    fragment_positions = CreateVA_CNS_AlignedContigElement(2 * num_frags);
    sequenceStore      = CreateVA_char(num_bases);
    qualityStore       = CreateVA_char(num_bases);
    columnStore        = CreateVA_Column(num_columns);
    beadStore          = CreateVA_Bead(num_bases + num_columns);
    manodeStore        = CreateVA_MANode(1);
  }

  ResetVA_Fragment(fragmentStore);
  MakeRoom_VA(fragmentStore, num_frags);

  ResetVA_int32(fragment_indices);
  MakeRoom_VA(fragment_indices, num_frags);

  ResetVA_int32(abacus_indices);

  ResetVA_CNS_AlignedContigElement(fragment_positions);
  MakeRoom_VA(fragment_positions, 2 * num_frags);

  ResetVA_char(sequenceStore);
  MakeRoom_VA(sequenceStore, num_bases);

  ResetVA_char(qualityStore);
  MakeRoom_VA(qualityStore, num_bases);

  ResetVA_Column(columnStore);
  MakeRoom_VA(columnStore, num_columns);

  ResetVA_Bead(beadStore);
  MakeRoom_VA(beadStore, num_bases + num_columns);

  ResetVA_MANode(manodeStore);

  gaps_in_alignment = 0;
}







////////////////////////////////////////
//
//  Variation Functions
//



//external
void
abacus::OutputReads(FILE *fout, Read *reads, int32 nr, int32 width) {
  int32 i, j;
  fprintf(fout, "\nReads =\n");

  for (i=0; i<nr; i++) {
    fprintf(fout, "%d   ", reads[i].allele_id);
    for (j=0; j<width; j++)
      fprintf(fout, "%c", reads[i].bases[j]);
    fprintf(fout, "\n");
  }
  fprintf(fout, "\n\n");
}

//external
void
abacus::OutputAlleles(FILE *fout, VarRegion *vreg) {
  int32 i, j;
  fprintf(fout,   "Outputting alleles:\n");
  fprintf(fout,   "nr= %d na= %d nca= %d\n", vreg->nr, vreg->na, vreg->nca);
  fprintf(fout,   "Num_reads= ");
  for (i=0; i<vreg->na; i++)
      fprintf(fout,   "%d ", vreg->alleles[i].num_reads);
  fprintf(fout,   "\n");
  fprintf(fout,   "Weights= ");
  for (i=0; i<vreg->na; i++)
      fprintf(fout,   "%d ", vreg->alleles[i].weight);
  fprintf(fout,   "\n");
  fprintf(fout,   "Reads= \n");
  for (i=0; i<vreg->na; i++) {
      fprintf(fout,   "   Allele order= %d, id= %d:\n", i, vreg->alleles[i].id);
      for (j=0; j<vreg->alleles[i].num_reads; j++) {
          int32 k, read_id = vreg->alleles[i].read_ids[j];
          int32 len = vreg->end-vreg->beg;
          fprintf(fout,   "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(fout,   "%c", vreg->reads[read_id].bases[k]);
          fprintf(fout,   "   %d\n", vreg->alleles[i].read_iids[j]);
        }
    }
}
