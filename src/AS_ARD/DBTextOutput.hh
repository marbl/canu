
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
#ifndef DBTextOutput_HH
#define DBTextOutput_HH

#include "IDBOutput.hh"

namespace AS_ARD {
   class DBTextOutput : public IDBOutput {
      private:
         // disable copy constructor
         DBTextOutput(DBTextOutput &);
         
      public:
         DBTextOutput();
         ~DBTextOutput();

         uint64 storeGenome(
                  const char * study,
                  const char * project,
                  const char * taxon);
         uint64 storeAssembly(
                  AS_UID assemblyEUID,
                  const char * date,
                  AS_UID genomeIID,
                  const char * op,
                  const char * genProg,
                  const char * ver,
                  const char status,
                  const char * notes);
         bool storeMDI2DB (
                  AS_UID erefines,  
                  IntDist_ID irefines,
                  float mean,
                  float stddev,
                  int32 min,
                  int32 max);
         bool storeAFG2DB (
                  AS_UID erefines,  
                  IntFragment_ID irefines,
                  MateStatType mate_status,
                  int32 chaff,
                  CDS_COORD_t begin,
                  CDS_COORD_t end);
         bool storeUTG2DB (
                  AS_UID eaccession,  
                  IntFragment_ID iaccession,
                  const char * source,
                  float mhp,
                  float coverage_stat,
                  UnitigStatus status,
                  CDS_COORD_t length,
                  const char * consensus,
                  const char * quality,
                  int32 forced,
                  int32 num_frags);
         bool storeMPS2DB (
                  AS_UID unitigID,                  
                  AS_UID afgID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         bool storeULK2DB (
                  AS_UID euid,
                  CDS_CID_t ciid,
                  ChunkOrientationType orientation,
                  UnitigOverlapType overlap_type,
                  int32 is_possible_chimera,
                  float mean_distance,
                  float std_deviation,
                  int32 num_contributing,
                  PlacementStatusType status);
         bool storeLKList2DB(int type, AS_UID ulkID, AS_UID utgID);
         bool storeJMP2DB(int jmpType, AS_UID jmpID, AS_UID ulkID, LinkType type);
         bool storeJMPList2DB(int jmpType, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID);
         
         bool storeCCO2DB (
                  AS_UID eaccession,  
                  IntFragment_ID iaccession,
                  ContigPlacementStatusType placed,
                  CDS_COORD_t length,
                  const char * consensus,
                  const char * quality,
                  int32 forced,
                  int32 num_pieces,
                  int32 num_unitigs,
                  int32 num_vars);
         bool storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,            
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         bool storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,            
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         bool storeVAR2DB(
                  AS_UID varID,
                  AS_UID ccoID,            
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  uint32 num_reads,
                  uint32 num_conf_alleles,
                  uint32 anchor_size,
                  CDS_COORD_t var_length,
                  int32 curr_var_id,
                  int32 phased_var_id);
         bool storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq);
         bool storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID);
         bool storeCLK2DB(
                  AS_UID euid,
                  CDS_CID_t ciid,
                  ChunkOrientationType orientation,
                  UnitigOverlapType overlap_type,
                  uint32 is_possible_chimera,
                  uint32 includes_guide,
                  float mean_distance,
                  float std_deviation,
                  uint32 num_contributing,
                  PlacementStatusType status);
         bool storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs);
         bool storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient);         
         bool storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID);
         bool storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd);

         bool commitMDI2DB() { return true; };
         bool commitAFG2DB() { return true; };
         bool commitUTG2DB() { return true; };
         bool commitMPS2DB() { return true; };
         bool commitULK2DB() { return true; };
         bool commitULKList2DB() { return true; };
         bool commitJMP2DB() { return true; };
         bool commitJMPList2DB() { return true; };
         bool commitCCO2DB() { return true; };
         bool commitCCOMPS2DB() { return true; };
         bool commitUPS2DB() { return true; };
         bool commitVAR2DB() { return true; };
         bool commitVARAllele2DB() { return true; };
         bool commitVARAFG2DB() { return true; };
         bool commitCLK2DB() { return true; };
         bool commitCLKList2DB() { return true; };
         bool commitCLKJMP2DB() { return true; };
         bool commitCLKJMPList2DB() { return true; };
         bool commitSCF2DB() { return true; };
         bool commitCTP2DB() { return true; };
         bool commitCTPList2DB() { return true; };
         bool commitCPS2DB() { return true; };
   };
}; 

#endif // DBTextOutput_HH
