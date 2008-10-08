
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

#ifndef BCPOutput_H
#define BCPOutput_H

static const char *rcsid_BCPOutput_H = "$Id: BCPOutput.hh,v 1.5 2008-10-08 22:02:54 brianwalenz Exp $";

#include <iostream>
#include "SQLOutput.hh"

namespace AS_ARD {
   class BCPOutput : public SQLOutput {
      private:
         // disable copy constructor
         BCPOutput(BCPOutput &);

         static const int MAX_PATH_LEN = 128;
         static const int MAX_FILE_LEN = 20;
         static const int MAX_STR_LEN  = 1024;

         static const char AFG_FILENAME[MAX_FILE_LEN];
         static const char MDI_FILENAME[MAX_FILE_LEN];
         static const char UTG_FILENAME[MAX_FILE_LEN];
         static const char MPS_FILENAME[MAX_FILE_LEN];
         static const char ULK_FILENAME[MAX_FILE_LEN];
         static const char ULK_LIST_FILENAME[MAX_FILE_LEN];
         static const char JMP_FILENAME[MAX_FILE_LEN];
         static const char JMP_LIST_FILENAME[MAX_FILE_LEN];
         static const char CLK_FILENAME[MAX_FILE_LEN];
         static const char CLK_LIST_FILENAME[MAX_FILE_LEN];
         static const char CLK_JMP_FILENAME[MAX_FILE_LEN];
         static const char CLK_JMP_LIST_FILENAME[MAX_FILE_LEN];
         static const char CCO_FILENAME[MAX_FILE_LEN];
         static const char CCO_MPS_FILENAME[MAX_FILE_LEN];
         static const char UPS_FILENAME[MAX_FILE_LEN];
         static const char VAR_FILENAME[MAX_FILE_LEN];
         static const char VAR_ALLELE_FILENAME[MAX_FILE_LEN];
         static const char VAR_AFG_FILENAME[MAX_FILE_LEN];
         static const char SCF_FILENAME[MAX_FILE_LEN];
         static const char CTP_FILENAME[MAX_FILE_LEN];
         static const char CTP_LIST_FILENAME[MAX_FILE_LEN];
         static const char CPS_FILENAME[MAX_FILE_LEN];


         static const char DEFAULT_BCP[MAX_STR_LEN];

         char * getFileName(const char * fileName);
         std::fstream * openFile(const char * fileName, std::_Ios_Openmode mode = (std::ios::in | std::ios::out | std::ios::trunc));
         bool closeFile(std::fstream **file);
         bool updateFile(
                  const char * fileName,
                  HashTable_AS * hash,
                  int32 position);
         bool runBCP(const char * fileName);
         bool runBCP(const char * fileName, bool eraseFile);

         // file output name prefix
         char *prefix;

         // connection parameters
         char *server;
         char *database;
         char *user;
         char *pass;
         char *bcp;

         std::fstream * mdiBCP;
         std::fstream * afgBCP;
         std::fstream * utgBCP;
         std::fstream * mpsBCP;
         std::fstream * ulkBCP;
         std::fstream * ulkListBCP;
         std::fstream * jmpBCP;
         std::fstream * jmpListBCP;
         std::fstream * clkBCP;
         std::fstream * clkListBCP;
         std::fstream * clkJmpBCP;
         std::fstream * clkJmpListBCP;
         std::fstream * ccoBCP;
         std::fstream * ccoMpsBCP;
         std::fstream * upsBCP;
         std::fstream * varBCP;
         std::fstream * varAlleleBCP;
         std::fstream * varAFGBCP;
         std::fstream * scfBCP;
         std::fstream * ctpBCP;
         std::fstream * ctpListBCP;
         std::fstream * cpsBCP;

      public:
         BCPOutput(
            const char * _prefix,
            const char * _server,
            const char * _database,
            const char * _user,
            const char * _password,
            const char * _bcp = NULL);
         ~BCPOutput();

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
                  AS_UID eident,
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
         bool storeLKList2DB(int jmpType, AS_UID utgID, AS_UID ulkID);
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


         bool commitMDI2DB();
         bool commitAFG2DB();
         bool commitUTG2DB();
         bool commitMPS2DB();
         bool commitULK2DB();
         bool commitULKList2DB();
         bool commitJMP2DB();
         bool commitJMPList2DB();
         bool commitCCO2DB();
         bool commitCCOMPS2DB();
         bool commitUPS2DB();
         bool commitVAR2DB();
         bool commitVARAllele2DB();
         bool commitVARAFG2DB();
         bool commitCLK2DB();
         bool commitCLKList2DB();
         bool commitCLKJMP2DB();
         bool commitCLKJMPList2DB();
         bool commitSCF2DB();
         bool commitCTP2DB();
         bool commitCTPList2DB();
         bool commitCPS2DB();
   };
};

#endif // BCPOutput_H
#endif //SYBASE
