
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

static char CM_ID[] = "$Id: ScaffoldUnitigProfile_CNS.c,v 1.3 2005-03-22 19:04:34 jason_miller Exp $";

// Operating System includes:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

// Celera Assembler includes:
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "AS_UTL_version.h"
#include "AS_SDB_SequenceDBPartition.h"
#include "AS_ALN_forcns.h"

// Consensus includes:
#include "MultiAlignStore_CNS.h"
#include "MultiAlignment_CNS.h"
#include "Globals_CNS.h"

int ScaffoldUnitigProfile(FILE *profileFile,
                          int scaffoldID,int scaffoldOffset,
                          MultiAlignT *contig,
                          tSequenceDB *sequenceDBp,
                          VA_TYPE(UnitigData) *unitigData);

float CNS_SEQUENCING_ERROR_EST = .02; // Used to calculate '-' probability
float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_PUBLIC = 0;   // Used to direct basecalling to include public data
int   CNS_CALL_PUBLIC = 0;   // Used to direct basecalling to favor public data
int   CNS_USE_QVS = 1;   // Used to direct basecalling to use quality value (versus strict majority rule)

VA_DEF(IntContigPairs)
void print_keys(void);

int OutputScaffoldUnitigProfile(FILE *profileFile,ScaffoldData *scaff,IntContigPairs *cp,tSequenceDB *contigStore,
                          tSequenceDB *sequenceDB, VA_TYPE(UnitigData) *unitigData);
int UnitigDataCmp( const void *l, const void *m);

int main (int argc, char *argv[]) {
   int use_store=0;
   int cns_store=0;
   int ch,errflg=0,illegal_use=0,iflags=0;
   char frgStoreFileName[FILENAME_MAX];
   char SeqStoreFileName[FILENAME_MAX];
   char cnsStoreFileName[FILENAME_MAX];
   char tigStoreFileName[FILENAME_MAX];
   int partition_no=0;
   int merge=0;
   int sdb_version=0;
   int doALL=0;
   tSequenceDB *sequenceDB=NULL;
   tSequenceDB *contigStore = NULL;
   InitializeAlphTable();
   optarg = NULL;
   while (!errflg && ((ch = getopt(argc, argv, "f:s:V:hq:c:p:A")) != EOF)) {
        switch(ch) {
        case 'q':
          sscanf(optarg,"%f:%d:%f",&CNS_SEQUENCING_ERROR_EST,&CNS_HAPLOTYPES,&CNS_SNP_RATE);
          if (!(CNS_SEQUENCING_ERROR_EST > 0) || CNS_SEQUENCING_ERROR_EST > .10 ) {
             fprintf(stderr,"ERROR: Sequencing error estimate (-q flag) should be "
                           "within (0,.10) (%4f was specified\n",CNS_SEQUENCING_ERROR_EST);
             illegal_use = 1;
          }
          if (CNS_HAPLOTYPES < 1) {
            fprintf(stderr,"ERROR: Haplotypes sampled (-h flag) must be > 0 "
                           "(%d was specified\n",CNS_HAPLOTYPES);
            illegal_use = 1;
          }
          if ((CNS_SNP_RATE < 0) || CNS_SNP_RATE > .10 ) {
             fprintf(stderr,"ERROR: SNP rate estimate (-s flag) should be within [0,.10) "
                            "(%4f was specified\n",CNS_SNP_RATE);
            illegal_use = 1;
          }
          iflags++;
          iflags++;
          break;
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
        case 'c':
          cns_store=1;
          strcpy(cnsStoreFileName, optarg);
          iflags++;
          iflags++;
          break;
        case 'p':
          partition_no = atoi(optarg);
          iflags++;
          iflags++;
          break;
        case 'A':
          doALL=1;
          iflags++;
          break;
        case 'h':
        default:
          illegal_use = 1;
          break;
       }
   }
   if ( !merge ) illegal_use=1;
   if ( !use_store ) illegal_use=1;
   if ( !cns_store ) illegal_use=1;
   if ( sdb_version == 0 ) illegal_use=1;
   if ( ! illegal_use && use_store ) {
    if ( (argc - iflags) != 1 ) {
      illegal_use = 1;   
    } 
   }
   if ( illegal_use ) {
        fprintf(stderr,"\n%s builds or examines a post-consensus multialign store\n\n",argv[0]);
        fprintf(stderr,"Usage:\n\n");
        fprintf(stderr,"  %s -f frgStore -s SeqStore -V SDB_version -c cnsStore -p part_number\n",argv[0]);
        fprintf(stderr,"    -f frgStore       Path to valid fragStore\n");
        fprintf(stderr,"    -s SeqStore -V #  Path to valid SeqStore and SDB version (checkpoint #)\n");
        fprintf(stderr,"    -c cnsStore -p #  Path to cnsStore to create, and SDB version of that store\n");
        fprintf(stderr,"    -a                Analyze (versus build) the post-consensus multialign store\n");
        exit(1);
   }
   ResetStores(LINE_MAX,20);

   global_fragStore = openFragStore(frgStoreFileName, "rb");
   if (global_fragStore == NULLSTOREHANDLE) return 0;
 
   sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);

   contigStore = OpenSequenceDB(cnsStoreFileName, FALSE, partition_no);

   // to start with, want to build an index of the scaffolds
   {
     char scaffStoreFileName[FILENAME_MAX];
     FILE *scaffstream;
     VA_TYPE(IntContigPairs) *ctpStore;
     VA_TYPE(ScaffoldData) *scaffStore;
# if 0
     {
     ctpStore = CreateVA_IntContigPairs(1000000);
     scaffStore = CreateVA_ScaffoldData(100000);
     reader = InputFileType_AS( stdin );
     while (reader(stdin,&pmesg) != EOF){
      if ( pmesg->t == MESG_ISF ) {
        ScaffoldData scaffData;
        IntScaffoldMesg *scaff = pmesg->m;
        MultiAlignT *contig;
        int i,num_pairs=scaff->num_contig_pairs;
        IntContigPairs *cp=scaff->contig_pairs;
        scaffData.ident=scaff->iaccession;  
        scaffData.num_contig_pairs=scaff->num_contig_pairs;  
        scaffData.contig_pairs= GetNumIntContigPairss(ctpStore);
        if (num_pairs==0) num_pairs++;
        contig= LoadMultiAlignTFromSequenceDB(contigStore, cp[0].contig1, FALSE);
        scaffData.length=GetNumchars(contig->consensus);
        for (i=0;i<num_pairs;i++ ) {
           AppendVA_IntContigPairs(ctpStore,&cp[i]);
           scaffData.length+=(int) cp[i].mean; 
           contig= LoadMultiAlignTFromSequenceDB(contigStore, cp[i].contig2, FALSE);
           scaffData.length+=GetNumchars(contig->consensus); 
        }
        // set length properly in the special case of singleton scaffolds
        if ( scaff->num_contig_pairs == 0 ) {
           scaffData.length = scaffData.length/2;
        }
        SetVA_ScaffoldData(scaffStore,scaffData.ident,&scaffData);
      }
     }
     sprintf(scaffStoreFileName,"%s/scaff.dat",cnsStoreFileName);
     scaffstream=fopen(scaffStoreFileName,"w");
     CopyToFileVA_ScaffoldData(scaffStore,scaffstream);
     CopyToFileVA_ScaffoldData(ctpStore,scaffstream);
     fclose(scaffstream);
     }
#endif
     {
     int scaffoldID=0;
     FILE *tigstream;
     VA_TYPE(UnitigData) *unitigData;
     FILE *scaffOutput = NULL;
     char scaffFileName[FILENAME_MAX];
     int scaffold_basepairs=0;
     int scaffold_batch_basepairs=0;
     int scaffold_total_basepairs=0;
     int batch_no=0;

     sprintf(tigStoreFileName,"%s/unitig.dat",cnsStoreFileName);
     tigstream=fopen(tigStoreFileName,"r");
     unitigData = CreateFromFileVA_UnitigData(tigstream,0);

     sprintf(scaffStoreFileName,"%s/scaff.dat",cnsStoreFileName);
     scaffstream=fopen(scaffStoreFileName,"r");
     scaffStore = CreateFromFileVA_ScaffoldData(scaffstream,0);
     ctpStore = CreateFromFileVA_IntContigPairs(scaffstream,0);

     
     if ( doALL ) {
       sprintf(scaffFileName,"scaffold_unitigs.%d",batch_no);
       scaffOutput=fopen(scaffFileName,"w");
     }

     while ( 1 ) {
        // lookup scaffold, and produce output for it
        IntContigPairs *cp;
        ScaffoldData *scaff=NULL;
        if ( !doALL ) {
          if ( EOF == fscanf(stdin,"%d",&scaffoldID) ) {
            fprintf(stderr,"Hey, this is the EOF already\n");
            break;
          }
          sprintf(scaffFileName,"scaffold_unitigs.scaff_%d",scaffoldID);
          scaffOutput=fopen(scaffFileName,"w");
        }
        if ( doALL ) {
           while ( 1 ) {
             if ( scaffoldID == GetNumScaffoldDatas(scaffStore)) break;
             scaff = GetScaffoldData(scaffStore,scaffoldID++);
             if ( scaff != NULL ) break;
           }
        } else {
           scaff = GetScaffoldData(scaffStore,scaffoldID);
        }
        if ( doALL && scaff == NULL ) break;
        if ( scaff == NULL ) {
            fprintf(stderr,"Could not locate scaffold %d\n",scaffoldID);
            continue;
        } 
        if (doALL && scaffold_batch_basepairs > 10000000 ) {
             batch_no++;
             scaffold_batch_basepairs=0;
             fprintf(stderr,"Starting scaffold batch %d with scaffold %d...\n",batch_no,scaffoldID);
             sprintf(scaffFileName,"scaffold_unitigs.%d",batch_no);
             scaffOutput=fopen(scaffFileName,"w");
        }
        cp = GetIntContigPairs(ctpStore,scaff->contig_pairs); 
        scaffold_basepairs = OutputScaffoldUnitigProfile(scaffOutput,scaff,cp,contigStore,sequenceDB,unitigData);
        scaffold_batch_basepairs+=scaffold_basepairs;
        scaffold_total_basepairs+=scaffold_basepairs;
     }

     }
   }
   return 0;
}

int OutputScaffoldUnitigProfile(FILE *profileFile,ScaffoldData *scaff,IntContigPairs *cp,tSequenceDB *contigStore,
                          tSequenceDB *sequenceDB, VA_TYPE(UnitigData) *unitigData) {
        // produce output for the given scaffold
        int scaffoldID=scaff->ident;
        int gapped_length=0;
        MultiAlignT *contig;
        int num_pairs=scaff->num_contig_pairs;
        contig= LoadMultiAlignTFromSequenceDB(contigStore, cp[0].contig1, FALSE);
        assert(contig->id==cp[0].contig1);
        if ( num_pairs == 0 ) {
           ScaffoldUnitigProfile(profileFile,scaffoldID,gapped_length,contig,sequenceDB,unitigData);
           gapped_length+=GetNumchars(contig->consensus);
        } else {
           int i;
           if ( cp[0].orient==AS_ANTI  || cp[0].orient==AS_OUTTIE ) {
             contig = RevcomplMultiAlignT(contig);
           }
           ScaffoldUnitigProfile(profileFile,scaffoldID,gapped_length,contig,sequenceDB,unitigData);
           gapped_length+=GetNumchars(contig->consensus);
           for (i=0;i<num_pairs;i++) {
               // output N's
               int gapsize;
               gapsize = (int) cp[i].mean;
               if (gapsize < 20 ) gapsize = 20; 
               gapped_length+=gapsize;
               contig=LoadMultiAlignTFromSequenceDB(contigStore, cp[i].contig2, FALSE);
               assert(contig->id==cp[i].contig2);
               if ( cp[i].orient==AS_ANTI  || cp[i].orient==AS_INNIE ) {
                contig = RevcomplMultiAlignT(contig);
               }
               ScaffoldUnitigProfile(profileFile,scaffoldID,gapped_length,contig,sequenceDB,unitigData);
               gapped_length+=GetNumchars(contig->consensus);
           } 
        }
        return gapped_length;
}

int ScaffoldUnitigProfile(FILE *profileFile,int scaffoldID,int scaffoldOffset, MultiAlignT *contig,
                          tSequenceDB *sequenceDBp, VA_TYPE(UnitigData) *unitigData) {
       int i;       
       int num_unitigs=GetNumIntUnitigPoss(contig->u_list);
       UnitigData *gatheredUnitigData=(UnitigData *) safe_malloc(num_unitigs*sizeof(UnitigData));
       MultiAlignT *uma;
       for (i=0;i<num_unitigs;i++) {
        int left,right;
        IntUnitigPos *tig=GetIntUnitigPos(contig->u_list,i);
        uma =  LoadMultiAlignTFromSequenceDB(sequenceDBp,tig->ident, TRUE);
        gatheredUnitigData[i]=*GetUnitigData(unitigData,tig->ident);
        if ( tig->position.bgn < tig->position.end ) {
           left = tig->position.bgn;
           right = tig->position.end;
        } else {
           left = tig->position.end;
           right = tig->position.bgn;
        }
        gatheredUnitigData[i].left=left; 
        gatheredUnitigData[i].right=right; 
        gatheredUnitigData[i].type=tig->type; 
       }
       qsort((void *)gatheredUnitigData, num_unitigs, sizeof(UnitigData), UnitigDataCmp);

       for (i=0;i<num_unitigs;i++) {
         UnitigData ud=gatheredUnitigData[i];
         fprintf(profileFile, "%d\t%d\t%d\t%c\t%f\t%d\t%d\t%d\n", 
                 scaffoldID, contig->id, ud.ident, ud.type, ud.coverage_stat,ud.length,
                 scaffoldOffset+ud.left,scaffoldOffset+ud.right);
       }
       return GetNumchars(contig->consensus);
}
