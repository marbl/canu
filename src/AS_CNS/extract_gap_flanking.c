
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
/* $Id: extract_gap_flanking.c,v 1.2 2004-09-23 20:25:20 mcschatz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  IntScaffoldMesg *isf;
  MesgReader   reader;
  IntScreenMatch *mat;
  IntContigPairs *pairs;
  IntConConMesg *contig;
  IntUnitigMesg *unitig;
  MultiAlignT *ma;
  MultiAlignT *ma1;
  MultiAlignT *ma2;
  MultiAlignStoreT *cstore = CreateMultiAlignStoreT();
  MultiAlignStoreT *ustore = CreateMultiAlignStoreT();
  IntMultiPos *frag1;
  IntMultiPos *frag2;
  int num_frag1;
  int num_frag2;
  int ma1_len;
  int ma2_len;
  VA_TYPE(int) *placed;
  int scaffold=-1;
  int num_pairs;
  int i;
  int isplaced = 1;
  int threshold=atoi(argv[1]);
  placed = CreateVA_int(800000);
  reader = InputFileType_AS( stdin );

 while (reader(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_IUM)  {
      unitig = pmesg->m;
      ma = CreateMultiAlignTFromIUM(unitig, unitig->iaccession,  0);
      SetMultiAlignInStore(ustore, unitig->iaccession, ma);
    }
    if (pmesg->t ==MESG_ICM)  {
      contig = pmesg->m;
      ma = CreateMultiAlignTFromICM(contig, contig->iaccession,  0);
      SetMultiAlignInStore(cstore, contig->iaccession, ma);
      if (contig->placed) {
        Setint(placed,contig->iaccession,&isplaced);
      }
    }
    if (pmesg->t ==MESG_ISF)  {
      scaffold++;
      isf = (IntScaffoldMesg *)pmesg->m; 
      pairs = isf->contig_pairs;
      num_pairs = isf->num_contig_pairs;
      for (i=0;i<num_pairs;i++){
        if (pairs[i].mean > threshold  )  {
           if ( Getint(placed,pairs[i].contig1)) {
             fprintf(stdout,"%d %d %d %f %c ",scaffold,pairs[i].contig1,pairs[i].contig2,pairs[i].mean,pairs[i].orient);
             ma1 = GetMultiAlignInStore(cstore,pairs[i].contig1);
             ma1_len = GetNumchars(ma1->consensus);
             ma2 = GetMultiAlignInStore(cstore,pairs[i].contig2);
             ma2_len = GetNumchars(ma2->consensus);
             // Check pieces list for end frags that are guides
             frag1 = GetIntMultiPos(ma1->f_list,0);
             frag2 = GetIntMultiPos(ma2->f_list,0);
             num_frag1 = GetNumIntMultiPoss(ma1->f_list);
             num_frag2 = GetNumIntMultiPoss(ma2->f_list);
             if (num_frag1 > 0 ) {
               int fi,leftend,rightend;
               if (pairs[i].orient == AB_AB || pairs[i].orient == AB_BA ) {
                 fi = num_frag1-1;
                 leftend = (frag1[fi].position.end > frag1[fi].position.bgn)?frag1[fi].position.end:
                                    frag1[fi].position.bgn;
                 while ( fi>0 && leftend > ma1_len-250) {
                  if (frag1[fi].type != AS_READ &&
                      frag1[fi].type != AS_B_READ &&
                      frag1[fi].type != AS_EXTR &&
                      frag1[fi].type != AS_TRNR) {
                     fprintf(stdout,"%d 1b %c ",
                          frag1[fi].ident,frag1[fi].type);
                  }
                  fi--;
                  leftend = (frag1[fi].position.end > frag1[fi].position.bgn)?frag1[fi].position.end:
                                    frag1[fi].position.bgn;
                 }
               } else {
                 fi = 0;
                 rightend = (frag1[fi].position.end < frag1[fi].position.bgn)?frag1[fi].position.end:
                                    frag1[fi].position.bgn;
                 while ( fi<num_frag1 && rightend < 250) {
                  if (frag1[fi].type != AS_READ &&
                      frag1[fi].type != AS_B_READ &&
                      frag1[fi].type != AS_EXTR &&
                      frag1[fi].type != AS_TRNR) {
                     fprintf(stdout,"%d 1e %c ",
                          frag1[fi].ident,frag1[fi].type);
                  }
                  fi++;
                  rightend = (frag1[fi].position.end < frag1[fi].position.bgn)?frag1[fi].position.end:
                                    frag1[fi].position.bgn;
                 }
               }
             }
             if (num_frag2 > 0 ) {
               int fi,leftend,rightend;
               if (pairs[i].orient == BA_BA || pairs[i].orient == AB_BA ) {
                 fi = num_frag2-1;
                 leftend = (frag2[fi].position.end > frag2[fi].position.bgn)?frag2[fi].position.end:
                                    frag2[fi].position.bgn;
                 while ( fi>0 && leftend > ma2_len-250) {
                  if (frag2[fi].type != AS_READ &&
                      frag2[fi].type != AS_B_READ &&
                      frag2[fi].type != AS_EXTR &&
                      frag2[fi].type != AS_TRNR) {
                     fprintf(stdout,"%d 2e %c ",
                          frag2[fi].ident,frag2[fi].type);
                  }
                  fi--;
                  leftend = (frag2[fi].position.end > frag2[fi].position.bgn)?frag2[fi].position.end:
                                    frag2[fi].position.bgn;
                 }
               } else {
                 fi = 0;
                 rightend = (frag2[fi].position.end < frag2[fi].position.bgn)?frag2[fi].position.end:
                                    frag2[fi].position.bgn;
                 while ( fi<num_frag2 && rightend < 250) {
                  if (frag2[fi].type != AS_READ &&
                      frag2[fi].type != AS_B_READ &&
                      frag2[fi].type != AS_EXTR &&
                      frag2[fi].type != AS_TRNR) {
                     fprintf(stdout,"%d 2b %c ",
                          frag2[fi].ident,frag2[fi].type);
                  }
                  fi++;
                  rightend = (frag2[fi].position.end < frag2[fi].position.bgn)?frag2[fi].position.end:
                                    frag2[fi].position.bgn;
                 }
               }
             }
           }
           fprintf(stdout,"\n");
           fflush(stdout);
        }
      }
   }
 }
 exit (0);
}
