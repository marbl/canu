
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

static const char *rcsid = "$Id: DBTextOutput.cc,v 1.7 2008-10-08 22:02:54 brianwalenz Exp $";

#include <iostream>
#include <assert.h>

#include "DBTextOutput.hh"

using AS_ARD::DBTextOutput;

DBTextOutput::DBTextOutput() :
      IDBOutput() {
}

DBTextOutput::~DBTextOutput() {
}

uint64 DBTextOutput::storeGenome(
         const char * study,
         const char * project,
         const char * taxon) {
   assert(study != NULL);
   assert(project != NULL);
   assert(taxon != NULL);

   std::cerr << "GEN: "
             << study << "\t"
             << project << "\t"
             << taxon << "\n";
   return 1;
}

uint64 DBTextOutput::storeAssembly(
         AS_UID assemblyEUID,
         const char * date,
         AS_UID genomeIID,
         const char * op,
         const char * genProg,
         const char * ver,
         const char status,
         const char * notes) {
   assert(date != NULL);
   assert(op != NULL);
   assert(genProg != NULL);
   assert(ver != NULL);
   assert(notes != NULL);

   std::cerr << "ASM: "
             << AS_UID_toInteger(assemblyEUID) << "\t"
             << date << "\t"
             << AS_UID_toInteger(genomeIID) << "\t"
             << op << "\t"
             << genProg << "\t"
             << ver << "\t"
             << status << "\t"
             << notes << "\n";

   return 1;
}

bool DBTextOutput::storeMDI2DB (
         AS_UID erefines,
         IntDist_ID irefines,
         float mean,
         float stddev,
         int32 min,
         int32 max) {
   std::cerr << "MDI: "
          << AS_UID_toInteger(erefines) << "\t"
          << irefines << "\t"
          << mean << "\t"
          << stddev << "\t"
          << min << "\t"
          << max  << "\n";

   return true;
}

bool DBTextOutput::storeAFG2DB (
         AS_UID eaccession,
         IntFragment_ID iaccession,
         MateStatType mate_status,
         int32 chaff,
         CDS_COORD_t bgn,
         CDS_COORD_t end) {
   std::cerr << "AFG: "
          << AS_UID_toInteger(eaccession) << "\t"
          << iaccession << "\t"
          << static_cast<char>(mate_status) << "\t"
          << chaff << "\t"
          << bgn  << "\t"
          << end  << "\n";

   return true;
}

bool DBTextOutput::storeUTG2DB (
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
         int32 num_frags) {
      std::cerr << "UTG: "
            << AS_UID_toInteger(eaccession) << "\t"
            << iaccession << "\t"
            << source << "\t"
            << mhp << "\t"
            << coverage_stat << "\t"
            << static_cast<char>(status) << "\t"
            // what are abp and bbp in CARD spec?
            << length << "\t"
            << consensus << "\t"
            << quality << "\t"
            << forced << "\t"
            << num_frags << "\n";

   return true;
}

bool DBTextOutput::storeMPS2DB (
         AS_UID unitigID,
         AS_UID afgID,
         FragType type,
         const char * source,
         CDS_COORD_t bgn,
         CDS_COORD_t end,
         int32 delta_length,
         std::string delta) {
   std::cerr << "MPS: "
         << AS_UID_toInteger(unitigID) << "\t"
         << AS_UID_toInteger(afgID) << "\t"
         << static_cast<char>(type) << "\t"
         << source << "\t"
         << bgn << "\t"
         << end << "\t"
         << delta_length << "\t"
         << delta << "\n";

   return true;
}

bool DBTextOutput::storeULK2DB (
         AS_UID euid,
         CDS_CID_t ciid,
         ChunkOrientationType orientation,
         UnitigOverlapType overlap_type,
         int32 is_possible_chimera,
         float mean_distance,
         float std_deviation,
         int32 num_contributing,
         PlacementStatusType status) {
   std::cerr << "ULK: "
         << AS_UID_toInteger(euid) << "\t"
         << ciid << "\t"
         << orientation << "\t"
         << overlap_type << "\t"
         << is_possible_chimera << "\t"
         //smp->includes_guide, - why no guide included?
         << mean_distance << "\t"
         << std_deviation << "\t"
         << num_contributing << "\t"
         << static_cast<char>(status) << "\n";

   return true;
}

bool DBTextOutput::storeLKList2DB(int type, AS_UID ulkID, AS_UID utgID) {
   if (type == ULK_TYPE) {
      std::cerr << "ULK-LIST" << "\t";
   }
   else if (type == CLK_TYPE) {
      std::cerr << "CLK-LIST" << "\t";
   }

   std::cerr
         << AS_UID_toInteger(ulkID)  << "\t"
         << AS_UID_toInteger(utgID)  << "\n";

   return true;
}

bool DBTextOutput::storeJMP2DB(int jmpType, AS_UID jmpID, AS_UID ulkID, LinkType type) {
   if (jmpType == ULK_TYPE) {
      std::cerr << "JMP: " << "\t";
   }
   else if (jmpType == CLK_TYPE) {
      std::cerr << "CLK-JMP: " << "\t";
   }

   std::cerr << AS_UID_toInteger(ulkID) << "\t"
         //what is status?
         << static_cast<char>(type) << "\n";

   return true;
}

bool DBTextOutput::storeJMPList2DB(int jmpType, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID) {
   if (jmpType == ULK_TYPE) {
      std::cerr << "JMP-LIST: " << "\t";
   }
   else if (jmpType == CLK_TYPE) {
      std::cerr << "CLK-JMP-LIST: " << "\t";
   }

   std::cerr << AS_UID_toInteger(jmpID) << "\t"
         << AS_UID_toInteger(fragID) << "\n";

   return true;
}

bool DBTextOutput::storeCCO2DB (
                  AS_UID eaccession,
                  IntFragment_ID iaccession,
                  ContigPlacementStatusType placed,
                  CDS_COORD_t length,
                  const char * consensus,
                  const char * quality,
                  int32 forced,
                  int32 num_pieces,
                  int32 num_unitigs,
                  int32 num_vars) {

   std::cerr << assemblyID << "\t"
               << AS_UID_toInteger(eaccession) << "\t"
               << iaccession << "\t"
               << static_cast<char>(placed) << "\t"
               << length << "\t"
               << "" << "\t"
               << "" << "\t"
               << forced << "\t"
               << num_pieces << "\t"
               << num_unitigs << "\t"
               << num_vars << "\n";

   return true;
}

bool DBTextOutput::storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {

   //TODO: warning using 0 as ciid for CCOMPS
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(ccoMpsID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(ccoID) << "\t"
                  << AS_UID_toInteger(fragID) << "\t"
                  << static_cast<char>(type) << "\t"
                  << source << "\t"
                  << bgn << "\t"
                  << end << "\t"
                  << delta.substr(0,MAX_DELTA).c_str() << "\n";

   return true;
}

bool DBTextOutput::storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {
   //TODO: warning truncating delta in UPS to 1000
   //	   using 0 as ciid for UPS
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(upsID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(ccoID) << "\t"
                  << AS_UID_toInteger(unitigID) << "\t"
                  << static_cast<char>(type) << "\t"
                  << bgn << "\t"
                  << end << "\t"
                  << delta.substr(0,MAX_DELTA).c_str() << "\n";

   return true;
}

bool DBTextOutput::storeVAR2DB(
                  AS_UID varID,
                  AS_UID ccoID,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  uint32 num_reads,
                  uint32 num_conf_alleles,
                  uint32 anchor_size,
                  CDS_COORD_t var_length,
                  int32 curr_var_id,
                  int32 phased_var_id) {
   //TODO: warning using 0 as ciid for VAR
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(varID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(ccoID) << "\t"
                  << bgn << "\t"
                  << end << "\t"
                  << num_reads << "\t"
                  << num_conf_alleles << "\t"
                  << anchor_size << "\t"
                  << var_length << "\t"
                  << curr_var_id << "\t"
                  << phased_var_id << "\n";

   return true;
}

bool DBTextOutput::storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq) {
   //TODO: warning using 0 as ciid for VAR_ALLELE
   std::cerr << assemblyID << "\t"
                     << AS_UID_toInteger(varAlleleID) << "\t"
                     << 0 << "\t"
                     << AS_UID_toInteger(varID) << "\t"
                     << nra << "\t"
                     << wgt << "\t"
                     << seq.c_str() << "\n";

   return true;
}

bool DBTextOutput::storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID) {
   //TODO: warning using 0 as ciid for VAR_AFG
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(varAfgID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(varID) << "\t"
                  << readID << "\n";

   return true;
}

bool DBTextOutput::storeCLK2DB(
                  AS_UID euid,
                  CDS_CID_t ciid,
                  ChunkOrientationType orientation,
                  UnitigOverlapType overlap_type,
                  uint32 is_possible_chimera,
                  uint32 includes_guide,
                  float mean_distance,
                  float std_deviation,
                  uint32 num_contributing,
                  PlacementStatusType status) {
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(euid) << "\t"
                  << ciid << "\t"
                  << static_cast<char>(orientation) << "\t"
                  << static_cast<char>(overlap_type) << "\t"
                  << is_possible_chimera << "\t"
                  << includes_guide << "\t"
                  << mean_distance << "\t"
                  << std_deviation << "\t"
                  << num_contributing << "\t"
                  << static_cast<char>(status) << "\n";

   return true;
}

bool DBTextOutput::storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs) {
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(eaccession) << "\t"
                  << iaccession << "\t"
                  << num_contig_pairs << "\n";

   return true;
}

bool DBTextOutput::storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient) {
   //TODO: warning using 0 as ciid for CTP
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(ctpID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(scfID) << "\t"
                  << mean << "\t"
                  << stddev << "\t"
                  << static_cast<char>(orient) << "\n";

   return true;
}

bool DBTextOutput::storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID) {
   std::cerr << assemblyID << "\t"
                  << AS_UID_toInteger(ctpListID) << "\t"
                  << (uint32)0 << "\t"
                  << AS_UID_toInteger(ctpID) << "\t"
                  << AS_UID_toInteger(ccoID) << "\n";

   return true;
}

bool DBTextOutput::storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd) {
   std::cerr << assemblyID << "\t"
               << AS_UID_toInteger(cpsID) << "\t"
               << (uint32)0 << "\t"
               << AS_UID_toInteger(ctpID) << "\t"
               << AS_UID_toInteger(ccoID) << "\t"
               << ctgStart << "\t"
               << ctgEnd << "\n";

   return true;
}

