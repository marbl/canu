
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

static char CM_ID[] = "$Id: TestBaseCall_CNS.c,v 1.12 2006-09-18 18:48:17 gdenisov Exp $";

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

float CNS_SEQUENCING_ERROR_EST = .02; // Used to calculate '-' probability
float CNS_SNP_RATE   = 0.0003; // Used to calculate BIAS
int   CNS_HAPLOTYPES = 1;   // Used to calculate BIAS
int   CNS_USE_PUBLIC = 0;   // Used to direct basecalling to include public data
int   CNS_CALL_PUBLIC = 0;   // Used to direct basecalling to favor public data
int   CNS_USE_QVS = 1;   // Used to direct basecalling to use quality value (versus strict majority rule)

void print_keys(void);

int main (int argc, char *argv[]) 
{
   char seq[LINE_MAX];
   char qlt[LINE_MAX];
   char frag_type[LINE_MAX];
   char unitig_type[LINE_MAX];
   VarRegion  vr;  
   char base;
   int cid=-1;
   int display_keys=0;
   int ch,errflg=0,illegal_use=0,iflags=0;
   double var = 0.;

   vr.nr = 0;       
   InitializeAlphTable();
   optarg = NULL;
   while (!errflg && ((ch = getopt(argc, argv, "q:d:hPim")) != EOF)) {
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
        case 'd': // depth of coverage at which to start using public data, 0 == alway include public
          {
            CNS_USE_PUBLIC = atoi(optarg);
          }
          iflags++;
          iflags++;
          break;
        case 'P': // favor the public data, call public if discrepant with Celera reads
          CNS_CALL_PUBLIC = 1;
          iflags++;
          break;
        case 'm': // use majority rule (rather than Bayesian evaluation of quality values)
          CNS_USE_QVS = 0;
          iflags++;
          break;
        case 'i':
          display_keys = 1; 
          break;
        case 'h':
        default:
          illegal_use = 1;
          break;
       }
   }
   if ( display_keys ) {
        print_keys(); 
        exit(1);
   }
   if ( illegal_use ) {
        fprintf(stderr,"  Usage:\n\n");
        fprintf(stderr,"  %s [-d int] [-q string] [-P] [-h] [-i] < formatted_column_input_file\n",argv[0]);
        fprintf(stderr,"    -d int       Depth of Celera coverage below which to include external data in basecalling\n");
        fprintf(stderr,"                    0 (default) indicates that external data should always be used\n");
        fprintf(stderr,"                    1 yields the traditional behavior, which uses external only in absence of Celera\n");
        fprintf(stderr,"                  > 1 will include publice data is the Celera depth falls below the given value\n");
        fprintf(stderr,"    -q string    Override default quality call parameters\n");
        fprintf(stderr,"                    string is colon separated list of the form '%%f:%%d:%%f'\n");
        fprintf(stderr,"                    where first field is estimated sequencing error rate (default: .015)\n");
        fprintf(stderr,"                         second field is number of sequenced haplotypes (default: 1)\n");
        fprintf(stderr,"                          third field is estimated SNP rate (default: 1/1000)\n");
        fprintf(stderr,"    -P           Call public data if is conflicts with Celera reads\n");
        fprintf(stderr,"    -h           print usage statement\n");
        fprintf(stderr,"    -i           print QV conversion + polymorphism key\n\n");
        exit(1);
   }
   fscanf(stdin,"%s",seq);
   fscanf(stdin,"%s",qlt);
   fscanf(stdin,"%s",frag_type);
   fscanf(stdin,"%s",unitig_type);
   ResetStores(LINE_MAX,20);
   cid = SetupSingleColumn(seq,qlt,frag_type,unitig_type, NULL);
   
   BaseCall(cid,CNS_USE_QVS, &var, &vr, vr.best_allele, &base, 1, 0, NULL);
   ShowColumn(cid);
   fprintf(stdout,"\nparameters:\n");
   fprintf(stdout,"             CNS_USE_QVS = %d\n",CNS_USE_QVS); 
   fprintf(stdout,"CNS_SEQUENCING_ERROR_EST = %4.2f\n",CNS_SEQUENCING_ERROR_EST);
   fprintf(stdout,"            CNS_SNP_RATE = %6.4f\n",CNS_SNP_RATE);
   fprintf(stdout,"          CNS_HAPLOTYPES = %d\n",CNS_HAPLOTYPES);  
   fprintf(stdout,"          CNS_USE_PUBLIC = %d\n",CNS_USE_PUBLIC); 
   fprintf(stdout,"         CNS_CALL_PUBLIC = %d\n",CNS_CALL_PUBLIC); 
   return 0;
}

void print_keys(void) {
fprintf(stdout,"Polymorphism encoding table:\n");
fprintf(stdout,"-\t-\n");
fprintf(stdout,"A\tA\n");
fprintf(stdout,"C\tC\n");
fprintf(stdout,"G\tG\n");
fprintf(stdout,"T\tT\n");
fprintf(stdout,"a\t-/A\n");
fprintf(stdout,"c\t-/C\n");
fprintf(stdout,"g\t-/G\n");
fprintf(stdout,"t\t-/T\n");
fprintf(stdout,"M\tA/C\n");
fprintf(stdout,"R\tA/G\n");
fprintf(stdout,"W\tA/T\n");
fprintf(stdout,"S\tC/G\n");
fprintf(stdout,"Y\tC/T\n");
fprintf(stdout,"K\tG/T\n");
fprintf(stdout,"m\t-/A/C\n");
fprintf(stdout,"r\t-/A/G\n");
fprintf(stdout,"w\t-/A/T\n");
fprintf(stdout,"s\t-/C/G\n");
fprintf(stdout,"y\t-/C/T\n");
fprintf(stdout,"k\t-/G/T\n");
fprintf(stdout,"V\tA/C/G\n");
fprintf(stdout,"H\tA/C/T\n");
fprintf(stdout,"D\tA/G/T\n");
fprintf(stdout,"B\tC/G/T\n");
fprintf(stdout,"v\t-/A/C/G\n");
fprintf(stdout,"h\t-/A/C/T\n");
fprintf(stdout,"d\t-/A/G/T\n");
fprintf(stdout,"b\t-/C/G/T\n");
fprintf(stdout,"X\t/NA/C/G/T\n");
fprintf(stdout,"x\t/n-/A/C/G/T\n\n");

fprintf(stdout,"Quality value encoding table:\n");
{ int i;
  for (i=0;i<61;i++) {
    fprintf(stdout,"%c = %d\n",(char)i+'0',i);
  }
  fprintf(stdout,"\n");
}

}
