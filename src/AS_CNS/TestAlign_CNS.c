
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
   Module:       TestBaseCall_CNS.c
   Description:  
   Assumptions:  
                 
 *********************************************************************/

static char CM_ID[] = "$Id: TestAlign_CNS.c,v 1.14 2007-09-01 05:09:49 brianwalenz Exp $";

// Operating System includes:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Celera Assembler includes:
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_version.h"
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_ALN_forcns.h"

// Consensus includes:
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Globals_CNS.h"

float CNS_SEQUENCING_ERROR_EST = .02; // Used to calculate '-' probability
float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_QVS = 1;   // Used to direct basecalling to use quality value (versus strict majority rule)

void print_keys(void);

int main (int argc, char *argv[]) {
   char FORWARD_PIC[]="--->";
   char REVERSE_PIC[]="<---";
   int use_store=0;
   int align_unitigs=FALSE;
   int afragID=0, bfragID=0;
   int alid = 0,blid = 0;
   int aComplement=0, bComplement=0;
   int ahang_est=0;
   int ch,errflg=0,illegal_use=0,iflags=0;
   int olap_success=0;
   char frgStoreFileName[FILENAME_MAX];
   char SeqStoreFileName[FILENAME_MAX];
   int merge=0;
   int sdb_version=0;
   int local_aligner=0;
   tSequenceDB *sequenceDB=NULL;

   InitializeAlphTable();
   optarg = NULL;
   while (!errflg && ((ch = getopt(argc, argv, "ulf:s:V:h")) != EOF)) {
        switch(ch) {
        case 'f':
          use_store = 1;
          strcpy(frgStoreFileName, optarg);
          iflags++;
          iflags++;
          break;
        case 's':
          merge = 1;
          strcpy(SeqStoreFileName, optarg);
          iflags++;
          iflags++;
          break;
        case 'V':
          sdb_version = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'u':
         align_unitigs=TRUE;
         iflags++;
         break;
        case 'l':
         local_aligner=1;
         iflags++;
         break;
        case 'h':
        default:
          illegal_use = 1;
          break;
       }
   }
   if ( merge && ! use_store ) illegal_use=1;
   if ( merge && sdb_version == 0 ) illegal_use=1;
   if ( ! illegal_use && use_store ) {
    if ( (argc - iflags) != 6 ) {
      illegal_use = 1;   
    } else {
      afragID = atoi(argv[optind++]);
      aComplement = atoi(argv[optind++]);
      bfragID = atoi(argv[optind++]);
      bComplement = atoi(argv[optind++]);
      ahang_est = atoi(argv[optind++]);
      fprintf(stderr,"Testing alignment of fragments %d %s and %d %s with estimated ahang %d\n",
             afragID,(!aComplement)?FORWARD_PIC:REVERSE_PIC,
             bfragID,(!bComplement)?FORWARD_PIC:REVERSE_PIC,
             ahang_est);
    }
   }
   if ( illegal_use ) {
        fprintf(stderr,"\n%s runs two specified fragments through the alignment routine within consensus\n\n",argv[0]);
        fprintf(stderr,"Usage:\n\n");
        fprintf(stderr,"  %s -f frgStore Aiid Aflip Biid Bflip ahang_est\n",argv[0]);
        fprintf(stderr,"    -f frgStore       Path to valid fragStore\n");
        fprintf(stderr,"    Aiid Biid         assembler iids of two fragments to align\n");
        fprintf(stderr,"    Aflip Bflip       orientation of the A and B fragments w.r.t. unitig\n");
        fprintf(stderr,"    ahang_est         estimated ahang of the fragment overlap\n\n");
        fprintf(stderr,"Alternate usage (MergeMultiAligns on contig pair):\n\n"); 
        fprintf(stderr,"  %s -f frgStore -s SeqStore -V SDB_version [-l] Aiid Aflip Biid Bflip ahang_est\n",argv[0]);
        fprintf(stderr,"    -f frgStore       Path to valid fragStore\n");
        fprintf(stderr,"    -s SeqStore -V #  Path to valid SeqStore and SDB version (checkpoint #)\n");
        fprintf(stderr,"    Aiid Biid         assembler iids of two contigs to align\n");
        fprintf(stderr,"    Aflip Bflip       orientation of the A and B contig w.r.t. scaffold\n");
        fprintf(stderr,"    ahang_est         estimated ahang of the contig overlap\n\n");
        fprintf(stderr,"    -l                (optional) try alignment with Local_Overlap_forCNS\n\n");

       // fprintf(stderr,"  %s < fasta_frag_input\n\n", argv[0]);
       // fprintf(stderr,"  %s  -h            prints usage statement\n",argv[0]);
        exit(1);
   }
   ResetStores(LINE_MAX,20);
   if ( use_store ) {
     gkpStore = openGateKeeperStore(frgStoreFileName, FALSE);
     if (gkpStore == NULL)
       return 0;
 
    if ( ! merge ) {
      alid= AppendFragToLocalStore(AS_READ, afragID, 0,0, 0x0, 0, 0x0);
      blid= AppendFragToLocalStore(AS_READ, bfragID, 0,0, 0x0, 0, 0x0);
    } else {
      sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);
    }
   } else {
    // read the fragments from stdin
   }
   if (!merge)   {  // test "consensus" program alignments (fragment x fragment)
      VA_TYPE(int32) *trace=CreateVA_int32(AS_READ_MAX_LEN);
      OverlapType otype;
      int ovl=GetFragment(fragmentStore,alid)->length-ahang_est;
      olap_success = GetAlignmentTrace(alid, 0, blid, &ahang_est, ovl, trace, &otype,DP_Compare,1,1);
      if ( ! olap_success ) {
       // try again, perhaps with alternate overlapper
       olap_success = GetAlignmentTrace(alid, 0, blid, &ahang_est, ovl, trace, &otype,Local_Overlap_AS_forCNS,1,1);
      }
   } else { // test MergeMultiAligns alignments (contig x contig)
     VA_TYPE(IntMultiPos) *pos_va=CreateVA_IntMultiPos(2);
     IntMultiPos pos;
     MultiAlignT *afragMA= LoadMultiAlignTFromSequenceDB(sequenceDB, afragID, align_unitigs);
     MultiAlignT *bfragMA= LoadMultiAlignTFromSequenceDB(sequenceDB, bfragID, align_unitigs);
     pos.type = AS_CONTIG;
     pos.ident = afragID;
     if ( aComplement ) {
       pos.position.bgn=GetMultiAlignUngappedLength(afragMA);
       pos.position.end=0;
     } else {
       pos.position.bgn=0;
       pos.position.end=GetMultiAlignUngappedLength(afragMA);
     }
     pos.contained=0;
     AppendVA_IntMultiPos(pos_va,&pos);
     pos.ident = bfragID;
     if ( bComplement ) {
       pos.position.bgn=ahang_est+GetMultiAlignUngappedLength(bfragMA);
       pos.position.end=ahang_est;
     } else {
       pos.position.bgn=ahang_est;
       pos.position.end=ahang_est+GetMultiAlignUngappedLength(bfragMA);
     }
     AppendVA_IntMultiPos(pos_va,&pos);
     fprintf(stderr,"Aligning contig %d (flipped %d)\n",afragID,aComplement);
     {int i;
      int numu=GetNumIntUnitigPoss(afragMA->u_list);
      IntUnitigPos *tig=GetIntUnitigPos(afragMA->u_list,0);
      for (i=0;i<numu;i++,tig++) {
         fprintf(stderr,"\tCI %d, [%d, %d]\n",tig->ident,tig->position.bgn,tig->position.end);
      }
     }
     fprintf(stderr,"and contig %d (flipped %d)\n",bfragID,bComplement);
     {int i; 
      int numu=GetNumIntUnitigPoss(bfragMA->u_list);
      IntUnitigPos *tig=GetIntUnitigPos(bfragMA->u_list,0);
      for (i=0;i<numu;i++,tig++) {
         fprintf(stderr,"\tCI %d, [%d, %d]\n",tig->ident,tig->position.bgn,tig->position.end);
      }
     }
      
     MergeMultiAligns( sequenceDB, gkpStore, pos_va, 1, 0, DP_Compare, NULL);
     if ( local_aligner ) 
        MergeMultiAligns( sequenceDB, gkpStore, pos_va, 1, 0, Local_Overlap_AS_forCNS, NULL);
   }
   return 0;
}

