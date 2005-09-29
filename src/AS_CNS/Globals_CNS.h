
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
#ifndef GLOBALS_CNS_INCLUDE
#define GLOBALS_CNS_INCLUDE

#include "MultiAlignment_CNS.h"
#include "AS_UTL_Var.h"
#include "AS_PER_fragStorePartition.h"
#include "AS_SDB_SequenceDB.h"
#include "AS_SDB_SequenceDBPartition.h"
#include <math.h>
#include "AS_UTL_PHash.h"

//====================================================================
// Input parameters for tuning BaseCall
// Defaults can be overridden with command line flag 

   float CNS_SEQUENCING_ERROR_EST;// Used to calculate '-' probability
   float CNS_SNP_RATE;            // Used to calculate BIAS
   int   CNS_HAPLOTYPES;          // Used to calculate BIAS
   int   CNS_USE_PUBLIC;          // Used to include public data in basecalling
   int   CNS_CALL_PUBLIC;         // Used to favor public data in basecalling

//====================================================================
// Output file:

   int std_output;
   FILE *cnsout;
   char OutputFileName[FILENAME_MAX];
   char OutputFileNameTmp[FILENAME_MAX];

//====================================================================
// Log file:

   int std_error_log;
   FILE *cnslog;
   char LogFileName[FILENAME_MAX];

//====================================================================
// Persistent store of the fragment/bactig data (produced upstream)

   FragStoreHandle global_fragStore;
   tFragStorePartition *global_fragStorePartition;
   tSequenceDB *sequenceDB;
   tSequenceDBPartition *sequenceDB_part;
   int partitioned;
   FragStoreHandle global_bactigStore;
   PHashTable_AS *fragmentMap;
   PHashTable_AS *bactigMap;

//====================================================================
// Store for the multialignments of Unitigs (referenced by contigging)

   MultiAlignStoreT *unitigStore;

//====================================================================
// Stores for the sequence/quality/alignment information
// (reset after each multialignment)

   VA_TYPE(char) *sequenceStore;
   VA_TYPE(char) *qualityStore;
   VA_TYPE(Bead) *beadStore;

//====================================================================
// Local stores for 
//      fragment information: 
//                indices into sequence/quality stores
//                index into "bead" store for alignment information
//      column information: 
//                basecall, profile, index in multialignment
//                indexed pointers to next and previous columns
//      multialignment information: 
//                first and last column, profile, index in multialignment
//                VA of all component columns
//
// (All are reset after each multialignment)

   VA_TYPE(Fragment) *fragmentStore;
   VA_TYPE(Column) *columnStore;
   VA_TYPE(MANode) *manodeStore;

//====================================================================
// Convenience arrays for misc. fragment information
// (All are reset after each multialignment)
   typedef enum ALIGN_TYPE { AS_CONSENSUS = (int) 'C', AS_OVERLAY = (int) 'O', AS_MERGE = (int) 'M' } ALIGN_TYPE;
   ALIGN_TYPE ALIGNMENT_CONTEXT; 
   int USE_SDB;
   int USE_SDB_PART;


   VA_TYPE(int32) *fragment_indices;
   VA_TYPE(int32) *abacus_indices;
   VA_TYPE(PtrT) *fragment_positions;
   VA_TYPE(PtrT) *fragment_source;
   int64 gaps_in_alignment;
   int debug_out;
   int terminate_cond;
   int allow_forced_frags;
   int allow_neg_hang;
   VA_TYPE(int32) *bactig_delta_length;
   VA_TYPE(PtrT) *bactig_deltas;

static void CleanExit(char *mesg, int lineno, int rc) {
  char command[100+FILENAME_MAX];

  fprintf(stderr,"%s at line: %d, rc: %d\n",mesg,lineno,rc);
  if( cnsout != NULL && ! std_output ){
    fclose(cnsout);
    if ( debug_out ) {
      sprintf(command,"mv -f %s %s.dbg",OutputFileName,OutputFileName);
    } else {
      sprintf(command,"rm -f %s",OutputFileName);
    }
    fprintf(stderr,"%s\n",command);
    system(command);
    sprintf(command,"touch %s",OutputFileName);
    system(command);
  }
  if( cnslog != NULL && ! std_error_log) {
    fclose(cnslog);
    if ( debug_out ) {
      sprintf(command,"mv -f %s %s.dbg",LogFileName,LogFileName);
    } else {
      sprintf(command,"rm -f %s",LogFileName);
    }
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  exit(terminate_cond);
  //exit(rc);
}

static float CNS_eff_cov[37] = {
1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75,
5.00, 5.25, 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75,
9.00, 9.25, 9.50, 9.75, 10.00};

static int CNS_ctg_len[37] = {
942, 1103, 1280, 1475, 1694, 1944, 2232, 2568, 2965, 3436, 4000, 4675, 5487, 6468, 7653, 9089,
10831, 12948, 15523, 18661, 22488, 27160, 32869, 39853, 48404, 58880, 71727, 87491, 106850, 130640, 
159893, 195885, 240195, 294778, 362050, 445004, 547346};

static float EffectiveCoverage(int contig_length) {
   int index;
   if ( contig_length == 0 ) return 0;
   if ( contig_length < 11000 ) {
      index = (int) (6.53 * log(contig_length)) - 44;
   } else {
      index = (int) (5.08 * log(contig_length)) - 31;
   }
   if ( index < 0 ) index = 0;
   if ( index > 36 ) index = 36;
   return CNS_eff_cov[index];
}
#endif
