
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




/* ******  README : NOTES ON PREREQUISITES TO REBUILDING SCAFFOLDS  ***


This program can perform a variety of analyses on the multiple alignments of
one or more scaffolds.

It relies on having a post-final-consensus "cnsStore", which includes information from the ICMs, ISFs and IUMs of an assembly.

To build such a store, each of these three message times get analysed separately.

PROCESSING ICMs:

First, construct a file containing the post-final-consensus ICMs (or use the
entire cns file, but this means a lot of protoIO IO).  
Then, use rebuildalignment, e.g., assuming the ICMs are in oct01.cns_contigs:

rebuildalignment -f oct01.fStore -s oct01.sStore -V 23 -c oct01.cStore -p 1 < oct01.cns_contigs


PROCESSING IUMs:

Construct a file containing the post-consensus IUMs (again, can use entire
cns file, but wasteful).  
Then construct a "unitig data" file, each line giving acc, A-state, len for a unitig; e.g. 

awk 'substr($1,1,4)=="acc:"{acc=substr($1,5)}substr($1,1,4)=="len:"{len=substr($1,5);print acc,cov,len}substr($1,1,4)=="cov:"{cov=substr($1,5)}' oct01.iums > oct01.utgData

Then process the unitig data to get it into the store:

rebuildalignment -f oct01.fStore -s oct01.sStore -V 23 -c oct01.cStore -p 2 -a -u oct01.utgData


PROCESSING ISFs:

First, get all isfs into a file, e.g. "oct01.isfs"

Then put the scaffold data into the cnsStore:

rebuildscaffolds -f oct01.fStore -s oct01.sStore -V 23 -c oct01.cStore -p 0 -b < oct01.isfs


RUNNING ANALYSES:

Run the relevant flavor of rebuildscaffolds; e.g. to get the multiple sequence
alignment of scaffold 53:

echo 53 | rebuildscaffolds -f oct01.fStore -s oct01.sStore -V 23 -c oct01.cStore -p 0 -m

 *********************************************************************/

static char CM_ID[] = "$Id: RebuildScaffolds_CNS.c,v 1.17 2007-05-29 10:54:28 brianwalenz Exp $";

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
#include "colCorr_CNS.h"

float CNS_SEQUENCING_ERROR_EST = .02; // Used to calculate '-' probability
float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_PUBLIC = 0;   // Used to direct basecalling to include public data
int   CNS_CALL_PUBLIC = 0;   // Used to direct basecalling to favor public data
int   CNS_USE_QVS = 1;   // Used to direct basecalling to use quality value (versus strict majority rule)

VA_DEF(IntContigPairs)
void print_keys(void);

int OutputScaffoldProfile(FILE *,ScaffoldData *,IntContigPairs *, 
    tSequenceDB *, tSequenceDB *, VA_TYPE(UnitigData) *, int, int , int);

int PrintColumnCorrelation(FILE *, MultiAlignT *, GateKeeperStore *,
    GateKeeperStore *, int , int , uint32 , int );

void OutputCoordinateMap(FILE *,int ,MultiAlignT *,tSequenceDB *,int *,int);

int main (int argc, char *argv[]) {
   int use_store=0;
   int cns_store=0;
   int buildScaffDat=FALSE;
   int ch,errflg=0,illegal_use=0,iflags=0;
   char frgStoreFileName[FILENAME_MAX];
   char SeqStoreFileName[FILENAME_MAX];
   char cnsStoreFileName[FILENAME_MAX];
   char tigStoreFileName[FILENAME_MAX];
   int partition_no=0;
   int merge=0;
   int sdb_version=0;
   int doALL=0;
   int SHOW_MULTIALIGN=0;
   int SHOW_COLUMNCORRELATION=0;
   int SHOW_COORDINATE_MAP=0;
   tSequenceDB *sequenceDB=NULL;
   tSequenceDB *contigStore = NULL;
   GenericMesg *pmesg;

   InitializeAlphTable();
   optarg = NULL;

   
   while (!errflg && ((ch = getopt(argc, argv, "f:s:V:hq:c:p:bAmCN")) != EOF)) {
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
	case 'b':
	  buildScaffDat=1;
          iflags++;
          break;
        case 'A':
          doALL=1;
          iflags++;
          break;
        case 'm':
          SHOW_MULTIALIGN=1;
          iflags++;
          break;
        case 'C':
          SHOW_COLUMNCORRELATION=1;
          iflags++;
          break;
	case 'N':
	  SHOW_COORDINATE_MAP=1;
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
        fprintf(stderr,
        "  %s -f frgStore -s SeqStore -V SDB_version -c cnsStore -p part_number [-q haplo_criteria] [-A] [-b|-m|-C|-N]\n",argv[0]);
        fprintf(stderr,"    -f frgStore       Path to valid fragStore\n");
        fprintf(stderr,"    -s SeqStore -V #  Path to valid SeqStore and SDB version (checkpoint #)\n");
        fprintf(stderr,"    -c cnsStore -p #  Path to cnsStore and SDB version of that store\n");
	fprintf(stderr,"    -b                build (vs analyze) scaffold data from ISFs on stdin\n");
        fprintf(stderr,"    -A                Process all scaffolds\n");
        fprintf(stderr,"    -C                Analyze confirmed-mismatch correlations\n");
        fprintf(stderr,"    -N                Generate gapped->ungapped coordinate offsets\n");
        fprintf(stderr,"    -m                Show multiple sequence alignment\n");
        fprintf(stderr,"    -q %%f:%%d:%%f       Set haplotype test information (SEQ_ERROR_EST:HAPLOTYPES:SNP_RATE)\n");
        exit(1);
   }
   ResetStores(LINE_MAX,20);

   gkpStore = openGateKeeperStore(frgStoreFileName, FALSE);
   if (gkpStore == NULL) return 0;
 
   sequenceDB = OpenSequenceDB(SeqStoreFileName, FALSE, sdb_version);

   contigStore = OpenSequenceDB(cnsStoreFileName, FALSE, partition_no);

   {
     char scaffStoreFileName[FILENAME_MAX];
     FILE *scaffstream;
     VA_TYPE(IntContigPairs) *ctpStore;
     VA_TYPE(ScaffoldData) *scaffStore;
     
     // if called for, build an index of the scaffolds
     if(buildScaffDat){
       ctpStore = CreateVA_IntContigPairs(1000000);
       scaffStore = CreateVA_ScaffoldData(100000);
       while (ReadProtoMesg_AS(stdin,&pmesg) != EOF){
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
	   if(contig!=NULL) {
	     scaffData.length=GetNumchars(contig->consensus)-1;
	   } else {
	     fprintf(stderr,
                 "%s: WARNING: couldn't find contig %d in cnsStore -- skipping\n",
                 argv[0],cp[0].contig1);
	   }
	   for (i=0;i<num_pairs;i++ ) {
	     AppendVA_IntContigPairs(ctpStore,&cp[i]);
	     scaffData.length+=(int) cp[i].mean; 
	     contig= LoadMultiAlignTFromSequenceDB(contigStore, cp[i].contig2, FALSE);
	     if(contig!=NULL) {
	       scaffData.length+=GetNumchars(contig->consensus)-1; 
	     } else {
	       fprintf(stderr,
                   "%s: WARNING: couldn't find contig %d in cnsStore -- skipping\n",
                   argv[0],cp[i].contig2);
	     }

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
       return(0);
     }

     {
     int scaffoldID=-1;
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
       if ( SHOW_MULTIALIGN ) {
           sprintf(scaffFileName,"scaffold_alignment.%d",batch_no);
       } else {
	 if ( SHOW_COLUMNCORRELATION) {
           sprintf(scaffFileName,"scaffold_correlations.%d",batch_no);
	 } else {
	   if ( SHOW_COORDINATE_MAP) {
	     sprintf(scaffFileName,"scaffold_coordmap.%d",batch_no);
	   } else {
	     sprintf(scaffFileName,"scaffold_profile.%d",batch_no);
	   }
	 }
       }
       scaffOutput=fopen(scaffFileName,"w");
     }

     while ( 1 ) {
        // lookup scaffold, and produce output for it
        IntContigPairs *cp;
        ScaffoldData *scaff=NULL;
        if ( !doALL ) {
          if ( fscanf(stdin,"%d",&scaffoldID) == EOF ) {
            break;
          }
          if ( SHOW_MULTIALIGN ) {
            sprintf(scaffFileName,"scaffold_alignment.scaff_%d",scaffoldID);
          } else {
	    if ( SHOW_COLUMNCORRELATION) {
	      sprintf(scaffFileName,"scaffold_correlations.scaff_%d",scaffoldID);
	    } else {
	      if ( SHOW_COORDINATE_MAP) {
		sprintf(scaffFileName,"scaffold_coordmap.scaff_%d",scaffoldID);
	      } else {
		sprintf(scaffFileName,"scaffold_profile.scaff_%d",scaffoldID);
	      }
	    }
          }
          scaffOutput=fopen(scaffFileName,"w");
        }
        if ( doALL ) {
           scaffoldID++;
           scaff = GetScaffoldData(scaffStore,scaffoldID);
           while ( scaff != NULL && scaff->ident != scaffoldID ) {
             scaffoldID++;
             scaff = GetScaffoldData(scaffStore,scaffoldID);
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
             fprintf(stderr,
                 "Starting scaffold batch %d with scaffold %d...\n",
                 batch_no,scaffoldID);
             if ( SHOW_MULTIALIGN ) {
                sprintf(scaffFileName,"scaffold_profile.%d",batch_no);
             } else {

	       if ( SHOW_COLUMNCORRELATION) {
		 sprintf(scaffFileName,"scaffold_correlations.%d",batch_no);
	       } else {
		 if ( SHOW_COORDINATE_MAP) {
		   sprintf(scaffFileName,"scaffold_coordmap.%d",batch_no);
		 } else {
		   sprintf(scaffFileName,"scaffold_alignment.%d",batch_no);
		 }
	       }
             }
             scaffOutput=fopen(scaffFileName,"w");
        }
        cp = GetIntContigPairs(ctpStore,scaff->contig_pairs); 
        scaffold_basepairs = OutputScaffoldProfile(scaffOutput,scaff,cp,
            contigStore,sequenceDB,unitigData,SHOW_MULTIALIGN,
            SHOW_COLUMNCORRELATION,SHOW_COORDINATE_MAP);
	fclose(scaffOutput);
        scaffold_batch_basepairs+=scaffold_basepairs;
        scaffold_total_basepairs+=scaffold_basepairs;
     }

     }
   }
   return 0;
}

int PrintColumnCorrelation(FILE *out,
			   MultiAlignT *ma,  
			   GateKeeperStore *frag_store,
			   GateKeeperStore *bactig_store,
			   int show_qv, 
			   int dots,
			   uint32 clrrng_flag,
			   int offset){
  int i=0;
  ColumnCorrelationT *cc;
  cc = test_correlated_columns(ma,frag_store);
  fprintf(out,"<<< begin Contig %d >>>\n",ma->id);
  if(cc!=NULL){
    while( cc[i].col!=-1){
      fprintf(out,"%d %d %d\n",cc[i].col+offset,cc[i].col,cc[i].corr);
      i++;
    }
  }
  return (i);
}



int OutputScaffoldProfile(FILE *profileFile,
    ScaffoldData *scaff,
    IntContigPairs *cp,
    tSequenceDB *contigStore,
    tSequenceDB *sequenceDB, 
    VA_TYPE(UnitigData) *unitigData, 
    int show_ma, 
    int show_colcorr, 
    int show_coordmap) {

        // produce output for the given scaffold
        int scaffoldID=scaff->ident;
	int mapoffset=0;
        int gapped_length=0;
        MultiAlignT *contig;
	int contigOk = 1;
        int num_pairs=scaff->num_contig_pairs;
        contig= LoadMultiAlignTFromSequenceDB(contigStore, cp[0].contig1, FALSE);
        if(contig->id!=cp[0].contig1){
	  fprintf(stderr,
              "OutputScaffoldProfile: WARNING: couldn't load contig %d\n",
              cp[0].contig1);
	  contigOk=0;
	}

        if ( num_pairs == 0 ) {

	  if( contigOk ) {
	    if ( show_ma ) {
	      PrintMultiAlignT(profileFile,contig,gkpStore,
                               1,
                               1,AS_READ_CLEAR_LATEST);
	    } else {
	      if ( show_colcorr ) {
		PrintColumnCorrelation(profileFile,contig,gkpStore,
                    NULL, 1,
                    1,AS_READ_CLEAR_LATEST,gapped_length);
	      } else 
		if ( show_coordmap) {
		  OutputCoordinateMap(profileFile,scaffoldID,contig,sequenceDB,
                      &mapoffset,gapped_length);
		} else {
		  MultiAlignContig_NoCompute(profileFile,scaffoldID,contig,
                      sequenceDB,unitigData, NULL);
		}
	    }
	    gapped_length+=GetNumchars(contig->consensus)-1;
	  }

        } else {

           int i;

	   if( contigOk ) {
	     if ( cp[0].orient==AS_ANTI  || cp[0].orient==AS_OUTTIE ) {
	       contig = RevcomplMultiAlignT(contig);
	     }
	     if ( show_ma ) {
	       fprintf(profileFile,"Scaffold offset %d\n",gapped_length);
	       PrintMultiAlignT(profileFile,contig,gkpStore,
                   1,1,AS_READ_CLEAR_LATEST);
	     } else {
	       if ( show_colcorr ) {
		 PrintColumnCorrelation(profileFile,contig,gkpStore,
                      NULL, 
                      1,1,AS_READ_CLEAR_LATEST,gapped_length);
	       } else {
		 if ( show_coordmap) {
		   OutputCoordinateMap(profileFile,scaffoldID,contig,sequenceDB,
                       &mapoffset,gapped_length);
		 } else {
		   MultiAlignContig_NoCompute(profileFile,scaffoldID,contig,
                       sequenceDB,unitigData, NULL);
		 }
	       }
	     }

	     gapped_length+=GetNumchars(contig->consensus)-1;

	   }

           for (i=0;i<num_pairs;i++) {
               // output N's
#ifdef LINES_FOR_GAP
               int j;
#endif
               int gapsize;
               gapsize = (int) cp[i].mean;

	       if(gapsize < 20 && show_coordmap) {
		 int i;
		 for(i=gapsize;i<20;i++){
		   mapoffset --;
		   //		   fprintf(profileFile,"%d %d\n", gapped_length+i,mapoffset);
		 }
	       }

#ifdef LINES_FOR_GAP
               if (gapsize < 20 ) gapsize = 20; 
               for (j=0;j<gapsize;j++) {
                  if ( ! show_ma && ! show_colcorr ) 
                      fprintf(profileFile,"%d\t0\t0\t0\tN\n",scaffoldID);
               }
#else
	       if ( ! show_ma ) fprintf(profileFile,"gap: %d\n",gapsize);
#endif

	       if(show_coordmap){
		 fprintf(profileFile,"%d %d\n", gapped_length+1,mapoffset);
	       }

               gapped_length+=gapsize;

               contig=LoadMultiAlignTFromSequenceDB(contigStore, cp[i].contig2, 
                   FALSE);
	       if(contig->id==cp[i].contig2){
		 contigOk=1;
	       } else {
		 fprintf(stderr,
                     "OutputScaffoldProfile: WARNING: couldn't load contig %d\n",
                      cp[i].contig2);

		 continue;
	       }

               if ( cp[i].orient==AS_ANTI  || cp[i].orient==AS_INNIE ) {
                contig = RevcomplMultiAlignT(contig);
               }


               if ( show_ma ) {
                  fprintf(profileFile,"Scaffold offset %d\n",gapped_length);
                  PrintMultiAlignT(profileFile,contig,gkpStore,
                       1,1,AS_READ_CLEAR_LATEST);
               } else {
		 if ( show_colcorr ) {
		   PrintColumnCorrelation(profileFile,contig,gkpStore,
                       NULL, 
                       1,1,AS_READ_CLEAR_LATEST,gapped_length);
		 } else {
		   if ( show_coordmap) {
		     OutputCoordinateMap(profileFile,scaffoldID,contig,sequenceDB,
                       &mapoffset,gapped_length);
		   } else {
		     MultiAlignContig_NoCompute(profileFile,scaffoldID,contig,
                       sequenceDB,unitigData, NULL);
		   }
		 }
	       }


               gapped_length+=GetNumchars(contig->consensus)-1;

           } 
	}
        return gapped_length;
}


void OutputCoordinateMap(FILE *profileFile,int scaffoldID,MultiAlignT *contig,
    tSequenceDB *sequenceDB,int *mapoffset,int gapped_length)
{
  int i,n;
  char *seq;
  n=GetNumchars(contig->consensus)-1;
  seq=Getchar(contig->consensus,0);
  for(i=0;i<n;i++){
    if(seq[i]=='-'){
      (*mapoffset)++;
      fprintf(profileFile,"%d %d\n",gapped_length+i,*mapoffset);
    }
  }
}
