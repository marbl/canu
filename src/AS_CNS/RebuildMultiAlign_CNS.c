
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

static char CM_ID[] = "$Id: RebuildMultiAlign_CNS.c,v 1.17 2007-09-01 05:09:49 brianwalenz Exp $";

// Operating System includes:
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // Linux has optarg here.

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
   int use_store=0;
   int cns_store=0;
   int createUnitigDat=FALSE;
   int ch,errflg=0,illegal_use=0,iflags=0;
   char frgStoreFileName[FILENAME_MAX];
   char SeqStoreFileName[FILENAME_MAX];
   char cnsStoreFileName[FILENAME_MAX];
   char tigStoreFileName[FILENAME_MAX];
   char unitigDefsFileName[FILENAME_MAX];
   int partition_no=0;
   int merge=0;
   int sdb_version=0;
   int SHOW_ALIGNMENT=0;
   tSequenceDB *sequenceDB=NULL;
   tSequenceDB *contigStore = NULL;
   GenericMesg *pmesg;
   int analyze=0;

   FILE *tigstream;
   MultiAlignT *ma;
   VA_TYPE(UnitigData) * unitigData;
   UnitigData dat;
   FILE *tigs;
   int32 acc;
   float   stat;
   int32 len;
   int counter=0;
   int n_frags;
   int contigID = 0;
   int contigFlip = 0;
   IntMultiPos *imps;
   
   InitializeAlphTable();
   optarg = NULL;
   while (!errflg && ((ch = getopt(argc, argv, "f:s:V:hq:c:p:au:v")) != EOF)) {
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
        case 'a':
          analyze = 1;
          iflags++;
          break;
	case 'u':
	  createUnitigDat=1;
	  strcpy(unitigDefsFileName, optarg);
          iflags++;
          iflags++;
          break;
        case 'v':
          SHOW_ALIGNMENT = 1;
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
        fprintf(stderr,"    -a [-u utgDat]    Analyze (versus build) the post-consensus multialign store\n");
        fprintf(stderr,"                      (or create unitig.dat file from the file of acc stat len\n");
        fprintf(stderr,"                      triplets in utgDat)\n");

        exit(1);
   }
   ResetStores(LINE_MAX,20);

   gkpStore = openGateKeeperStore(frgStoreFileName, FALSE);
   if (gkpStore == NULL)
     return 0;
 
   sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);

   if ( partition_no == 1 && ! analyze) {
     contigStore = CreateSequenceDB(cnsStoreFileName, 1000000, TRUE);
   } else {
     contigStore = OpenSequenceDB(cnsStoreFileName, !analyze , partition_no-2);
   }
   if ( !analyze) { // this block creates the cnsStore
     while (ReadProtoMesg_AS(stdin,&pmesg) != EOF){
      if ( pmesg->t == MESG_ICM ) {
	IntConConMesg *contig=pmesg->m;

	if(contig->length!=strlen(contig->consensus)||
	   contig->length!=strlen(contig->quality)){
	  fprintf(stderr,"%s: WARNING: skipping contig %d: length mismatch\n",
		  argv[0],contig->iaccession);
	  continue;
	}

        ma = CreateMultiAlignTFromICM(contig,contig->iaccession,0);
        InsertMultiAlignTInSequenceDB( contigStore, contig->iaccession, FALSE, ma, FALSE);
        if ( partition_no==1) {
          InsertMultiAlignTInSequenceDB( contigStore, contig->iaccession, TRUE, ma, FALSE);
        }
        //MultiAlignContig_NoCompute(ma,sequenceDB);
      }
     }
     SaveSequenceDB(contigStore);
   } else {
   // if we're in analyze mode (-a option) examine things in the store, don't write to it (except if we still have to construct the unitig data!)

     sprintf(tigStoreFileName,"%s/unitig.dat",cnsStoreFileName);

     if(createUnitigDat){
     
       tigs=fopen(unitigDefsFileName,"r");
       counter=0;


       unitigData=CreateVA_UnitigData(1000000);
       while (  fscanf(tigs,"%d %f %d",&acc,&stat,&len) !=EOF ) {
	 ma = LoadMultiAlignTFromSequenceDB(sequenceDB, acc, TRUE);
	 imps=GetIntMultiPos(ma->f_list,0);
	 n_frags=GetNumIntMultiPoss(ma->f_list);
	 dat.ident = acc;
	 dat.length = len;
	 dat.coverage_stat = stat;
	 // dat.status=is_tig_simple(imps,n_frags,len,gkpStore,ali,thresh,variant,&dat.prob);
	 SetVA_UnitigData(unitigData,acc,&dat);
	 counter++;
	 //	 if(counter == 1000 * (int) (counter/1000))
	 //	   fprintf(stderr,"working on unitig %d\n",acc);
       }
       tigstream=fopen(tigStoreFileName,"w");
       CopyToFileVA_UnitigData(unitigData,tigstream);
       fclose(tigstream);
       return(0);

     }
 
     tigstream=fopen(tigStoreFileName,"r");
     unitigData = CreateFromFileVA_UnitigData(tigstream,0);
 
    while ( fscanf(stdin,"%d %d",&contigID,&contigFlip) != EOF) {
     tSequenceDB *seqStore;
     seqStore = contigStore;
     ma = LoadMultiAlignTFromSequenceDB(seqStore, contigID, FALSE);
     // this should be replaced with something that's computed once and referenced later
     if ( contigFlip ) {
        ma = RevcomplMultiAlignT(ma);
     }
     //MultiAlignContig_NoCompute(stdout,-1,ma,sequenceDB,unitigData,1);
     if ( SHOW_ALIGNMENT ) {
       PrintMultiAlignT(stdout,
                        ma,
                        gkpStore,
			1,0,AS_READ_CLEAR_LATEST);
     }
    }
   }
   return 0;
}

