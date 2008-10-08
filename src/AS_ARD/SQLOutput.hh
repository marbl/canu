
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
#ifdef SYBASE

#ifndef SQLOutput_H
#define SQLOutput_H

static const char *rcsid_SQLOutput_H = "$Id: SQLOutput.hh,v 1.5 2008-10-08 22:02:54 brianwalenz Exp $";

#include <iostream>
#include "Sybase.hh"
#include "IDBOutput.hh"
#include "IDBConnection.hh"

extern "C" {
   #include <ctpublic.h>
   #include "AS_UTL_Hash.h"
}

namespace AS_ARD {
   class SQLOutput : public IDBOutput {
      private:
         // disable copy constructor
         SQLOutput(SQLOutput &);

         IDBConnection * dbConnection;
protected:
         HashTable_AS *AFG_UID_to_MSGID;
         HashTable_AS *AFG_IID_to_MSGID;
         HashTable_AS *UTG_UID_to_MSGID;
         HashTable_AS *ULK_UID_to_MSGID;
         HashTable_AS *JMP_UID_to_MSGID;
         HashTable_AS *CLK_UID_to_MSGID;
         HashTable_AS *CLK_JMP_UID_to_MSGID;
         HashTable_AS *CCO_UID_to_MSGID;
         HashTable_AS *VAR_UID_to_MSGID;
         HashTable_AS *SCF_UID_to_MSGID;
         HashTable_AS *CTP_UID_to_MSGID;

      public:
         SQLOutput(IDBConnection * connection);
         ~SQLOutput();

         IDBConnection * getConnection() {
            return dbConnection;
         }

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

         virtual bool storeMDI2DB (
                  AS_UID erefines,
                  IntDist_ID irefines,
                  float mean,
                  float stddev,
                  int32 min,
                  int32 max);
         virtual bool storeAFG2DB (
                  AS_UID erefines,
                  IntFragment_ID irefines,
                  MateStatType mate_status,
                  int32 chaff,
                  CDS_COORD_t begin,
                  CDS_COORD_t end);
         virtual bool storeUTG2DB (
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
         virtual bool storeMPS2DB (
                  AS_UID unitigID,
                  AS_UID afgID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         virtual bool storeULK2DB (
                  AS_UID euid,
                  CDS_CID_t ciid,
                  ChunkOrientationType orientation,
                  UnitigOverlapType overlap_type,
                  int32 is_possible_chimera,
                  float mean_distance,
                  float std_deviation,
                  int32 num_contributing,
                  PlacementStatusType status);
         virtual bool storeLKList2DB(int jmpType, AS_UID utgID, AS_UID ulkID);
         virtual bool storeJMP2DB(int jmpType, AS_UID jmpID, AS_UID ulkID, LinkType type);
         virtual bool storeJMPList2DB(int jmpType, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID);
         virtual bool storeCCO2DB (
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
         virtual bool storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         virtual bool storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta);
         virtual bool storeVAR2DB(
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
         virtual bool storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq);
         virtual bool storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID);
         virtual bool storeCLK2DB(
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
         virtual bool storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs);
         virtual bool storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient);
         virtual bool storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID);
         virtual bool storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd);

         virtual bool commitMDI2DB();
         virtual bool commitAFG2DB();
         virtual bool commitUTG2DB();
         virtual bool commitMPS2DB();
         virtual bool commitULK2DB();
         virtual bool commitULKList2DB();
         virtual bool commitJMP2DB();
         virtual bool commitJMPList2DB();

         virtual bool commitCCO2DB();

         virtual bool commitCCOMPS2DB();
         virtual bool commitUPS2DB();
         virtual bool commitVAR2DB();
         virtual bool commitVARAllele2DB();
         virtual bool commitVARAFG2DB();
         virtual bool commitCLK2DB();
         virtual bool commitCLKList2DB();
         virtual bool commitCLKJMP2DB();
         virtual bool commitCLKJMPList2DB();

         virtual bool commitSCF2DB();
         virtual bool commitCTP2DB();
         virtual bool commitCTPList2DB();
         virtual bool commitCPS2DB();
   };
};

#endif // BCPOutput_H
#endif //SYBASE
