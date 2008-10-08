
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
#ifndef IDBOutput_HH
#define IDBOutput_HH

static const char *rcsid_IDBOutput_HH = "$Id: IDBOutput.hh,v 1.6 2008-10-08 22:02:54 brianwalenz Exp $";

extern "C" {
   #include "AS_global.h"
   #include "AS_MSG_pmesg.h"
}

namespace AS_ARD {
   //#define ULK_TYPE 0
   //#define CLK_TYPE 1

   class IDBOutput {
      private:
         // disable copy constructor
         IDBOutput(IDBOutput &);

      public:
         static const int MAX_DELTA = 1000;
         static const int ULK_TYPE = 0;
         static const int CLK_TYPE = 1;

         IDBOutput() {};
         virtual ~IDBOutput() {};

         virtual bool storeMDI2DB (
                  AS_UID erefines,
                  IntDist_ID irefines,
                  float mean,
                  float stddev,
                  int32 min,
                  int32 max) = 0;

         virtual bool storeAFG2DB (
                  AS_UID erefines,
                  IntFragment_ID irefines,
                  MateStatType mate_status,
                  int32 chaff,
                  CDS_COORD_t begin,
                  CDS_COORD_t end) = 0;

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
                  int32 num_frags) = 0;

         virtual bool storeMPS2DB (
                  AS_UID unitigID,
                  AS_UID eident,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) = 0;

         virtual bool storeULK2DB (
                  AS_UID euid,
                  CDS_CID_t ciid,
                  ChunkOrientationType orientation,
                  UnitigOverlapType overlap_type,
                  int32 is_possible_chimera,
                  float mean_distance,
                  float std_deviation,
                  int32 num_contributing,
                  PlacementStatusType status) = 0;

         virtual bool storeLKList2DB(int type, AS_UID utgID, AS_UID ulkID) = 0;

         virtual bool storeJMP2DB(int type, AS_UID jmpID, AS_UID ulkID, LinkType type) = 0;


         virtual bool storeJMPList2DB(int type, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID) = 0;

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
                  int32 num_vars) = 0;

         virtual bool storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) = 0;

         virtual bool storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) = 0;

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
                  int32 phased_var_id) = 0;

         virtual bool storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq) = 0;

         virtual bool storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID)= 0;

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
                  PlacementStatusType status) = 0;

         virtual bool storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs) = 0;

         virtual bool storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient) = 0;
         virtual bool storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID) = 0;
         virtual bool storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd) = 0;

         virtual uint64 storeGenome(
                  const char * study,
                  const char * project,
                  const char * taxon) = 0;
         virtual uint64 storeAssembly(
                  AS_UID assemblyEUID,
                  const char * date,
                  AS_UID genomeIID,
                  const char * op,
                  const char * genProg,
                  const char * ver,
                  const char status,
                  const char * notes) = 0;

         virtual bool commitMDI2DB() = 0;
         virtual bool commitAFG2DB() = 0;
         virtual bool commitUTG2DB() = 0;
         virtual bool commitMPS2DB() = 0;
         virtual bool commitULK2DB() = 0;
         virtual bool commitULKList2DB() = 0;
         virtual bool commitJMP2DB() = 0;
         virtual bool commitJMPList2DB() = 0;
         virtual bool commitCCO2DB() = 0;
         virtual bool commitCCOMPS2DB() = 0;
         virtual bool commitUPS2DB() = 0;
         virtual bool commitVAR2DB() = 0;
         virtual bool commitVARAllele2DB() = 0;
         virtual bool commitVARAFG2DB() = 0;
         virtual bool commitCLK2DB() = 0;
         virtual bool commitCLKList2DB() = 0;
         virtual bool commitCLKJMP2DB() = 0;
         virtual bool commitCLKJMPList2DB() = 0;
         virtual bool commitSCF2DB() = 0;
         virtual bool commitCTP2DB() = 0;
         virtual bool commitCTPList2DB() = 0;
         virtual bool commitCPS2DB() = 0;

         uint64 assemblyID;
   };
};

#endif // IDBOutput_HH
