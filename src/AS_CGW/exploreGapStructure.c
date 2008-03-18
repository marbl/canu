
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

#include "AS_global.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "MultiAlignment_CNS.h"

int
main (int argc , char * argv[] ) {
  int            ckptNum = NULLINDEX;

  int            DistFromGap = 20;
  MultiAlignT   *ma = NULL;

  int            sid = 0;

  GlobalData = CreateGlobal_CGW();

  argc = AS_configure(argc, argv);
 
  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      strcpy(GlobalData->File_Name_Prefix, argv[++arg]);
    } else if (strcmp(argv[arg], "-d") == 0) {
      DistFromGap = atoi(argv[++arg]);
      if (DistFromGap > 9999) {
        fprintf(stderr,"-d option too large.\n");
        exit(1);
      }
    } else if (strcmp(argv[arg], "-g") == 0) {
      strcpy(GlobalData->Gatekeeper_Store_Name, argv[++arg]);
    } else if (strcmp(argv[arg], "-n") == 0) {
      ckptNum = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err > 0) ||
      (GlobalData->File_Name_Prefix[0] == 0) ||
      (GlobalData->Gatekeeper_Store_Name[0] == 0)) {
    fprintf(stderr, "usage: %s -g <gkpStore> -c <ckptName> -n <ckptNum> -d distance_from_end\n",
            argv[0]);
    exit(1);
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint(GlobalData->File_Name_Prefix, ckptNum, FALSE);

  ma            = CreateEmptyMultiAlignT();

  for (sid = 0; sid < GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldTIterator CIs;
    CIScaffoldT * scaff;
    ContigT *contig;
    float endPrev;
    float varPrev;
	
    scaff = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    if ((scaff == NULL) ||
        (isDeadCIScaffoldT(scaff)) ||
        (scaff->type != REAL_SCAFFOLD))
      continue;

    fprintf(stdout,"Working on scaffold %d\n",sid);

    InitCIScaffoldTIterator( ScaffoldGraph, scaff, TRUE, FALSE, &CIs);
    endPrev=-1;
    while ( (contig = NextCIScaffoldTIterator( &CIs )) != NULL) {
      float begThis,varThis;

      if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
        begThis=contig->offsetAEnd.mean;
        varThis=contig->offsetAEnd.variance;
      } else {
        begThis=contig->offsetBEnd.mean;
        varThis=contig->offsetBEnd.variance;
      }

      if (endPrev != -1)
        fprintf(stdout,"  gap size %f , %f\n", begThis-endPrev, sqrt(varThis-varPrev));
		
      MultiAlignT  *newma = NULL;

      VA_TYPE(IntElementPos) *positions = CreateVA_IntElementPos(5000);
      MultiAlignT *ma =  loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contig->id, FALSE);
      int num_tigs = GetNumIntUnitigPoss(ma->u_list);
      int i;
      for (i=0;i<num_tigs;i++) {
        IntUnitigPos *upos = GetIntUnitigPos(ma->u_list,i);
        IntElementPos pos;
        pos.type = AS_UNITIG;
        pos.ident = upos->ident;
        pos.position.bgn = upos->position.bgn;
        pos.position.end = upos->position.end;
        SetVA_IntElementPos(positions,i,&pos);
      }
   
      newma = MergeMultiAlignsFast_new(ScaffoldGraph->sequenceDB, NULL, positions, 0, 1, NULL, NULL);

      int nfr = GetNumIntMultiPoss(newma->f_list);
      int len = GetMultiAlignLength(newma);
      int frontCnt = 0;
      int tailCnt  = 0;

      //fprintf(stdout,"Numfrgs in contig %d , len = %d\n",nfr,len);

      for(i=0;i<nfr;i++) {
        IntMultiPos *mpos = GetIntMultiPos(newma->f_list,i);
        int beg = mpos->position.bgn;
        int end = mpos->position.end;

        if(end< beg){
          int tmp=end;
          end=beg;
          beg=tmp;
        }

        if(beg<DistFromGap)
          frontCnt++;
        if(len-end<DistFromGap)
          tailCnt++;
      }

      DeleteVA_IntElementPos(positions);

      fprintf(stdout,"    %4d  <-- frgs within %9dbp of end --> %4d\n",
              frontCnt, DistFromGap, tailCnt);

      if(contig->offsetAEnd.mean<contig->offsetBEnd.mean){
        endPrev=contig->offsetBEnd.mean;
        varPrev=contig->offsetBEnd.variance;
      } else {
        endPrev=contig->offsetAEnd.mean;
        varPrev=contig->offsetAEnd.variance;
      }
    }
  }
}
