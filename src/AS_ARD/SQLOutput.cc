
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

static const char *rcsid = "$Id: SQLOutput.cc,v 1.8 2008-10-29 06:34:30 brianwalenz Exp $";

#ifdef SYBASE

#include <iostream>
#include <string>

#include "SQLOutput.hh"

using AS_ARD::SQLOutput;
using AS_ARD::Sybase;

SQLOutput::SQLOutput(IDBConnection * connection) {
   dbConnection = connection;

   AFG_UID_to_MSGID     = NULL;
   AFG_IID_to_MSGID     = NULL;
   UTG_UID_to_MSGID     = NULL;
   ULK_UID_to_MSGID     = NULL;
   JMP_UID_to_MSGID     = NULL;
   CLK_UID_to_MSGID     = NULL;
   CLK_JMP_UID_to_MSGID = NULL;
   CCO_UID_to_MSGID     = NULL;
   VAR_UID_to_MSGID     = NULL;
   SCF_UID_to_MSGID     = NULL;
   CTP_UID_to_MSGID     = NULL;

}

SQLOutput::~SQLOutput() {
   delete dbConnection;
}

uint64 SQLOutput::storeGenome(
         const char * study,
         const char * project,
         const char * taxon) {
   assert(study != NULL);
   assert(project != NULL);
   assert(taxon != NULL);

   char cmd[Sybase::MAX_STR_LEN];

   // enforce that only one genome can be loaded into a database
   if (getConnection()->getCount("Genome") != 0) { return 0; }

   sprintf(cmd, "INSERT INTO Genome (Study, Project, Taxon) VALUES ('%s', '%s', '%s')", study, project, taxon);
   if (getConnection()->sqlCommand(cmd) == false) { return 0; }

   return getConnection()->getLast("ID", "Genome");
}

uint64 SQLOutput::storeAssembly(
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

   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO Assembly " \
            "(Creation, Genome_ID, Operator, GeneratingProgram, ProgramVersion, Status, Notes) " \
            "VALUES ('%s', "F_U64", '%s', '%s', '%s', '%c', '%s')",
                     date,
                     AS_UID_toInteger(genomeIID),
                     op,
                     genProg,
                     ver,
                     status,
                     notes
            );
   if (getConnection()->sqlCommand(cmd) == false) { return 0; }

   assemblyID = getConnection()->getLast("ID", "Assembly");
   return assemblyID;
}

bool SQLOutput::storeMDI2DB (
         AS_UID erefines,
         IntDist_ID irefines,
         float mean,
         float stddev,
         int32 min,
         int32 max) {

   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO MDI " \
            "(mdi_AssemblyID, mdi_EUID, mdi_CIID, mdi_mea, mdi_std, mdi_min, mdi_max) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", %f, %f, "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(erefines),
                     irefines,
                     mean,
                     stddev,
                     min,
                     max
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeAFG2DB (
         AS_UID eaccession,
         IntFragment_ID iaccession,
         MateStatType mate_status,
         int32 chaff,
         CDS_COORD_t bgn,
         CDS_COORD_t end) {

   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO AFG " \
            "(afg_AssemblyID, afg_EUID, afg_CIID, afg_mst, afg_cha, afg_clr1, afg_clr2) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", '%c', "F_S32", "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(eaccession),
                     iaccession,
                     static_cast<char>(mate_status),
                     chaff,
                     bgn,
                     end
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (AFG_UID_to_MSGID == NULL) {
      AFG_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   if (AFG_IID_to_MSGID == NULL) {
      AFG_IID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }

   uint64 afg = (uint64)getConnection()->getLast("afg_MSG_ID", "AFG");
   InsertInHashTable_AS(AFG_UID_to_MSGID, AS_UID_toInteger(eaccession), 0, static_cast<uint64>(afg), 0);
   InsertInHashTable_AS(AFG_IID_to_MSGID, static_cast<uint64>(iaccession), 0, AS_UID_toInteger(eaccession), 0);

   return true;
}

bool SQLOutput::storeUTG2DB (
         AS_UID eaccession,
         IntFragment_ID iaccession,
         const char * source,
         float microhet_prob,
         float coverage_stat,
         UnitigStatus status,
         CDS_COORD_t length,
         const char * consensus,
         const char * quality,
         int32 forced,
         int32 num_frags) {

   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO UTG " \
            "(utg_AssemblyID, utg_EUID, utg_CIID, utg_src, utg_mhp, utg_cov, " \
            " utg_sta, utg_abp, utg_bbp, utg_len, utg_cns, utg_qlt, utg_for, utg_nfr) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", '%s', %f, %f, '%c', "F_S32", " \
                     F_S32", "F_S32", '%s', '%s', "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(eaccession),
                     iaccession,
                     source,
                     microhet_prob,
                     coverage_stat,
                     static_cast<char>(status),
                     0,
                     0,
                     length,
                     "",
                     "",
                     forced,
                     num_frags
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (UTG_UID_to_MSGID == NULL) {
      UTG_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(UTG_UID_to_MSGID, AS_UID_toInteger(eaccession), 0, static_cast<uint64>(getConnection()->getLast("utg_MSG_ID", "UTG")), 0);

   return true;
}

bool SQLOutput::storeMPS2DB (
         AS_UID unitigID,
         AS_UID afgID,
         FragType type,
         const char * source,
         CDS_COORD_t bgn,
         CDS_COORD_t end,
         int32 delta_length,
         std::string delta) {
   char cmd[Sybase::MAX_STR_LEN];
   uint64 utg = LookupValueInHashTable_AS(UTG_UID_to_MSGID, AS_UID_toInteger(unitigID), 0);
   uint64 afg = LookupValueInHashTable_AS(AFG_UID_to_MSGID, AS_UID_toInteger(afgID), 0);

   // should be mps_mid not mps_afg_MSG_ID for consistency with asm
   sprintf(cmd,
            "INSERT INTO MPS " \
            "(mps_AssemblyID, mps_utg_MSG_ID, mps_afg_MSG_ID, mps_type, mps_src, mps_pos1, " \
            " mps_pos2, mps_del) " \
            "VALUES ("F_U64", "F_U64", "F_U64", '%c', '%s', "F_S32", "F_S32", '%s')",
                     assemblyID,
                     utg,
                     afg,
                     static_cast<char>(type),
                     source,
                     bgn,
                     end,
                     delta.c_str()
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeULK2DB (
         AS_UID euid,
         CDS_CID_t ciid,
         ChunkOrientationType orientation,
         UnitigOverlapType overlap_type,
         int32 is_possible_chimera,
         float mean_distance,
         float std_deviation,
         int32 num_contributing,
         PlacementStatusType status) {
   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO ULK " \
            "(ulk_assemblyID, ulk_EUID, ulk_CIID, ulk_ori, ulk_ovt, ulk_ipc, " \
            " ulk_mea, ulk_std, ulk_num, ulk_sta) " \
            "VALUES ("F_CID", '"F_U64"', "F_CID", '%c', '%c', "F_S32", %f, %f, "F_S32", '%c')",
                     assemblyID,
                     AS_UID_toInteger(euid),
                     ciid,
                     static_cast<char>(orientation),
                     static_cast<char>(overlap_type),
                     is_possible_chimera,
                     mean_distance,
                     std_deviation,
                     num_contributing,
                     static_cast<char>(status)
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (ULK_UID_to_MSGID == NULL) {
      ULK_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(ULK_UID_to_MSGID, AS_UID_toInteger(euid), 0, static_cast<uint64>(getConnection()->getLast("ulk_MSG_ID", "ULK")), 0);

   return true;
}

bool SQLOutput::storeLKList2DB(int jmpType, AS_UID utgID, AS_UID ulkID) {
   char cmd[Sybase::MAX_STR_LEN];
   uint64 utg;
   uint64 ulk;

std::cerr << "STORING ULK " << AS_UID_toInteger(utgID) << " and " << AS_UID_toInteger(ulkID) << std::endl;
   if (jmpType == ULK_TYPE) {
      sprintf(cmd, "INSERT INTO ULK_LIST (ulk_list_AssemblyID, ulk_list_UTG_ID, ulk_list_ULK_ID) ");
      utg = static_cast<uint64>(LookupValueInHashTable_AS(UTG_UID_to_MSGID, AS_UID_toInteger(utgID), 0));
      ulk = static_cast<uint64>(LookupValueInHashTable_AS(ULK_UID_to_MSGID, AS_UID_toInteger(ulkID), 0));
   }
   else if (jmpType == CLK_TYPE) {
      sprintf(cmd, "INSERT INTO CLK_LIST (clk_list_AssemblyID, clk_list_cco_MSG_ID, clk_list_clk_MSG_ID) ");
      utg = static_cast<uint64>(LookupValueInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(utgID), 0));
      ulk = static_cast<uint64>(LookupValueInHashTable_AS(CLK_UID_to_MSGID, AS_UID_toInteger(ulkID), 0));
   }
   else {
      assert(0);
   }
   sprintf(cmd,
            "%s " \
            "VALUES ("F_U64", "F_U64", "F_U64")",
                     cmd,
                     assemblyID,
                     utg,
                     ulk
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeJMP2DB(int jmpType, AS_UID jmpID, AS_UID ulkID, LinkType type) {
   char cmd[Sybase::MAX_STR_LEN];
   uint64 ulk = 0;

std::cerr << "STORING JMP " << AS_UID_toInteger(jmpID) << " and " << AS_UID_toInteger(ulkID) << std::endl;

   if (jmpType == ULK_TYPE) {
      // should change name of jmp_unitig to jmp_utg_MSG_ID
      sprintf(cmd, "INSERT INTO JMP (jmp_assemblyID, jmp_EUID, jmp_CIID, jmp_utg_MSG_ID, jmp_status, jmp_type) ");
      ulk = static_cast<uint64>(LookupValueInHashTable_AS(ULK_UID_to_MSGID, AS_UID_toInteger(ulkID), 0));

      sprintf(cmd,
               "%s VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", '%c', '%c')",
                        cmd,
                        assemblyID,
                        AS_UID_toInteger(jmpID),
                        0,
                        ulk,
                        'X',
                        static_cast<char>(type)
               );
   }
   else if (jmpType == CLK_TYPE) {
      sprintf(cmd, "INSERT INTO CLK_JMP (clk_jmp_assemblyID, clk_jmp_EUID, clk_jmp_CIID, clk_jmp_cco_MSG_ID, clk_jmp_status) ");
      ulk = static_cast<uint64>(LookupValueInHashTable_AS(CLK_UID_to_MSGID, AS_UID_toInteger(ulkID), 0));

      sprintf(cmd,
               "%s VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", '%c')",
                        cmd,
                        assemblyID,
                        AS_UID_toInteger(jmpID),
                        0,
                        ulk,
                        static_cast<char>(type)
               );
   }
   else {
      assert(0);
   }

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (JMP_UID_to_MSGID == NULL) {
      JMP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   if (CLK_JMP_UID_to_MSGID == NULL) {
      CLK_JMP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }

   if (jmpType == ULK_TYPE) {
      // should be jmp_MSG_ID not jmp_id
      InsertInHashTable_AS(JMP_UID_to_MSGID, AS_UID_toInteger(jmpID), 0, static_cast<uint64>(getConnection()->getLast("jmp_MSG_ID", "JMP")), 0);
   }
   else if (jmpType == CLK_TYPE) {
      InsertInHashTable_AS(CLK_JMP_UID_to_MSGID, AS_UID_toInteger(jmpID), 0, static_cast<uint64>(getConnection()->getLast("clk_jmp_MSG_ID", "CLK_JMP")), 0);
   }
   else {
      assert(0);
   }

   return true;
}

bool SQLOutput::storeJMPList2DB(int jmpType, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID) {
   char cmd[Sybase::MAX_STR_LEN];
   uint64 jmp = 0;
   uint64 afg = 0;

std::cerr << "STORING JMP LIST" << AS_UID_toInteger(jmpID) << " and " << AS_UID_toInteger(fragID) << std::endl;

   if (jmpType == ULK_TYPE) {
      // should be jmp_list_jmp_MSG_ID instead of jmp_list_jmp_id
      // should be jmp_list_afg_MSG_ID instead of jmp_list_frag
      sprintf(cmd, "INSERT INTO JMP_LIST (jmp_list_AssemblyID, jmp_list_EUID, jmp_list_CIID, jmp_list_jmp_id, jmp_list_MSG_ID, jmp_list_frag)");
      jmp = static_cast<uint64>(LookupValueInHashTable_AS(JMP_UID_to_MSGID, AS_UID_toInteger(jmpID), 0));
      afg = static_cast<uint64>(LookupValueInHashTable_AS(AFG_UID_to_MSGID, AS_UID_toInteger(fragID), 0));
   }
   else if (jmpType == CLK_TYPE) {
      sprintf(cmd, "INSERT INTO CLK_JMP_LIST (clk_jmp_list_AssemblyID, clk_jmp_list_EUID, clk_jmp_list_CIID, clk_jmp_list_clk_jmp_MSG_ID, clk_jmp_list_MSG_ID, clk_jmp_list_afg_MSG_ID)");
      jmp = static_cast<uint64>(LookupValueInHashTable_AS(CLK_JMP_UID_to_MSGID, AS_UID_toInteger(jmpID), 0));
      afg = static_cast<uint64>(LookupValueInHashTable_AS(AFG_UID_to_MSGID, AS_UID_toInteger(fragID), 0));
   }
   else {
      assert(0);
   }

   sprintf(cmd,
            "%s " \
            "VALUES ("F_U64", "F_U64", "F_CID", "F_U64", "F_U64")",
                     cmd,
                     assemblyID,
                     AS_UID_toInteger(jmpListID),
                     0,
                     jmp,
                     afg
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeCCO2DB (
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
   char cmd[Sybase::MAX_STR_LEN];

std::cerr << "Storing CCO " << std::endl;
   sprintf(cmd,
            "INSERT INTO CCO " \
            "(cco_AssemblyID, cco_EUID, cco_CIID, cco_pla, cco_len, cco_cns, cco_qlt, " \
            " cco_for, cco_npc, cco_nou, cco_nvr) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", '%c', "F_S32", '%s', '%s', " \
                     F_S32", "F_S32", "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(eaccession),
                     iaccession,
                     static_cast<char>(placed),
                     length,
                     "",
                     "",
                     forced,
                     num_pieces,
                     num_unitigs,
                     num_vars
            );
   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (CCO_UID_to_MSGID == NULL) {
      CCO_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(eaccession), 0, static_cast<uint64>(getConnection()->getLast("cco_MSG_ID", "CCO")), 0);

   return true;
}

bool SQLOutput::storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {
   char cmd[Sybase::MAX_STR_LEN];

std::cerr << "Storing CCOMPS " << std::endl;
   uint64 cco = static_cast<uint64>(LookupValueInHashTable_AS(UTG_UID_to_MSGID, AS_UID_toInteger(ccoID), 0));
   uint64 afg = static_cast<uint64>(LookupValueInHashTable_AS(AFG_UID_to_MSGID, AS_UID_toInteger(fragID), 0));

   sprintf(cmd,
            "INSERT INTO CCO_MPS " \
            "(cco_mps_AssemblyID, cco_mps_EUID, cco_mps_CIID, cco_mps_cco_MSG_ID, cco_mps_mid, cco_mps_type, cco_mps_src, cco_mps_pos1, " \
            " cco_mps_pos2, cco_mps_del) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_U64", '%c', '%s', "F_S32", "F_S32", '%s')",
                     assemblyID,
                     AS_UID_toInteger(ccoMpsID),
                     0,
                     cco,
                     afg,
                     static_cast<char>(type),
                     source,
                     bgn,
                     end,
                     delta.c_str()
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {
   char cmd[Sybase::MAX_STR_LEN];

std::cerr << "Storing UPS " << std::endl;
   uint64 cco = static_cast<uint64>(LookupValueInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(ccoID), 0));
   uint64 utg = static_cast<uint64>(LookupValueInHashTable_AS(UTG_UID_to_MSGID, AS_UID_toInteger(unitigID), 0));

   //TODO: warning truncating delta in UPS to 1000
   sprintf(cmd,
            "INSERT INTO UPS " \
            "(ups_AssemblyID, ups_EUID, ups_CIID, ups_cco_MSG_ID, ups_utg_MSG_ID, ups_type, ups_pos1, " \
            " ups_pos2, ups_del) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_U64", '%c', "F_S32", "F_S32", '%s')",
                     assemblyID,
                     AS_UID_toInteger(upsID),
                     0,
                     cco,
                     utg,
                     static_cast<char>(type),
                     bgn,
                     end,
                     delta.substr(0,MAX_DELTA).c_str()
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeVAR2DB(
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
   char cmd[Sybase::MAX_STR_LEN];
std::cerr << "Storing VAR " << std::endl;

   uint64 cco = static_cast<uint64>(LookupValueInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(ccoID), 0));

   sprintf(cmd,
            "INSERT INTO VAR " \
            "(var_AssemblyID, var_EUID, var_CIID, var_cco_MSG_ID, var_pos1, var_pos2, " \
            " var_nrd, var_nca, var_anc, var_len, var_vid, var_pid) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_S32", "F_S32", " \
                     F_S32", "F_S32", "F_S32", "F_S32", "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(varID),
                     0,
                     cco,
                     bgn,
                     end,
                     num_reads,
                     num_conf_alleles,
                     anchor_size,
                     var_length,
                     curr_var_id,
                     phased_var_id
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (VAR_UID_to_MSGID == NULL) {
      VAR_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(VAR_UID_to_MSGID, AS_UID_toInteger(varID), 0, (uint64)getConnection()->getLast("var_MSG_ID", "VAR"), 0);

   return true;
}

bool SQLOutput::storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq) {
   char cmd[Sybase::MAX_STR_LEN];
std::cerr << "Storing VARALL " << std::endl;

   uint64 var = static_cast<uint64>(LookupValueInHashTable_AS(VAR_UID_to_MSGID, AS_UID_toInteger(varID), 0));

   sprintf(cmd,
            "INSERT INTO VAR_ALLELE " \
            "(var_allele_AssemblyID, var_allele_EUID, var_allele_CIID, var_allele_var_MSG_ID, " \
            " var_allele_nra, var_allele_wgt, var_ellele_seq) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", " \
                     F_S32", "F_S32", '%s')",
                     assemblyID,
                     AS_UID_toInteger(varAlleleID),
                     0,
                     var,
                     nra,
                     wgt,
                     seq.c_str()
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID) {
   char cmd[Sybase::MAX_STR_LEN];

std::cerr << "Storing VAR AFG " << std::endl;

   uint64 var = static_cast<uint64>(LookupValueInHashTable_AS(VAR_UID_to_MSGID, AS_UID_toInteger(varID), 0));
   uint64 euid = static_cast<uint64>(LookupValueInHashTable_AS(AFG_IID_to_MSGID, static_cast<uint64>(readID), 0));
   uint64 afg = static_cast<uint64>(LookupValueInHashTable_AS(AFG_UID_to_MSGID, euid, 0));

   sprintf(cmd,
            "INSERT INTO VAR_AFG " \
            "(var_afg_AssemblyID, var_afg_EUID, var_afg_CIID, var_afg_var_MSG_ID, var_afg_afg_MSG_ID) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_U64")",
                     assemblyID,
                     AS_UID_toInteger(varAfgID),
                     0,
                     var,
                     afg
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeCLK2DB(
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
   char cmd[Sybase::MAX_STR_LEN];
std::cerr << "Storing CLK " << std::endl;
std::cerr << "The value os mean and std is " << mean_distance << " " << std_deviation << std::endl;
   sprintf(cmd,
            "INSERT INTO CLK " \
            "(clk_AssemblyID, clk_EUID, clk_CIID, clk_ori, clk_ovt, clk_ipc, " \
            " clk_gui, clk_mea, clk_std, clk_num, clk_sta) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", '%c', '%c', "F_S32", "F_S32", %f, %f, "F_S32", '%c')",
                     assemblyID,
                     AS_UID_toInteger(euid),
                     ciid,
                     static_cast<char>(orientation),
                     static_cast<char>(overlap_type),
                     is_possible_chimera,
                     includes_guide,
                     mean_distance,
                     std_deviation,
                     num_contributing,
                     static_cast<char>(status)
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (CLK_UID_to_MSGID == NULL) {
      CLK_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(CLK_UID_to_MSGID, AS_UID_toInteger(euid), 0, static_cast<uint64>(getConnection()->getLast("clk_MSG_ID", "CLK")), 0);

   return true;
}

bool SQLOutput::storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs) {
   char cmd[Sybase::MAX_STR_LEN];

   sprintf(cmd,
            "INSERT INTO SCF " \
            "(scf_AssemblyID, scf_EUID, scf_CIID, scf_noc) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID" ,"F_U64")",
                     assemblyID,
                     AS_UID_toInteger(eaccession),
                     iaccession,
                     num_contig_pairs
            );

   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (SCF_UID_to_MSGID == NULL) {
      SCF_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(SCF_UID_to_MSGID, AS_UID_toInteger(eaccession), 0, static_cast<uint64>(getConnection()->getLast("scf_MSG_ID", "SCF")), 0);

   return true;
}

bool SQLOutput::storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient) {
   char cmd[Sybase::MAX_STR_LEN];

   uint64 scf = static_cast<uint64>(LookupValueInHashTable_AS(SCF_UID_to_MSGID, AS_UID_toInteger(scfID), 0));
   sprintf(cmd,
            "INSERT INTO CTP " \
            "(ctp_AssemblyID, ctp_EUID, ctp_CIID, ctp_scf_MSG_ID, ctp_mea, ctp_std, ctp_ori) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", %f, %f, '%c')",
                     assemblyID,
                     AS_UID_toInteger(ctpID),
                     0,
                     scf,
                     mean,
                     stddev,
                     static_cast<char>(orient)
            );
   if (getConnection()->sqlCommand(cmd) == false) { return false; }

   if (CTP_UID_to_MSGID == NULL) {
      CTP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(CTP_UID_to_MSGID, AS_UID_toInteger(ctpID), 0, (uint64)getConnection()->getLast("ctp_MSG_ID", "CTP"), 0);

   return true;
}

bool SQLOutput::storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID) {
   char cmd[Sybase::MAX_STR_LEN];

   uint64 ctp = static_cast<uint64>(LookupValueInHashTable_AS(CTP_UID_to_MSGID, AS_UID_toInteger(ctpID), 0));
   uint64 cco = static_cast<uint64>(LookupValueInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(ccoID), 0));

   sprintf(cmd,
            "INSERT INTO CTP_LIST " \
            "(ctp_list_AssemblyID, ctp_list_EUID, ctp_list_CIID, ctp_list_ctp_MSG_ID, ctp_list_cco_MSG_ID) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_U64")",
                     assemblyID,
                     AS_UID_toInteger(ctpListID),
                     (uint32)0,
                     ctp,
                     cco
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd) {
   char cmd[Sybase::MAX_STR_LEN];

   uint64 ctp = static_cast<uint64>(LookupValueInHashTable_AS(CTP_UID_to_MSGID, AS_UID_toInteger(ctpID), 0));
   uint64 cco = static_cast<uint64>(LookupValueInHashTable_AS(CCO_UID_to_MSGID, AS_UID_toInteger(ccoID), 0));

   sprintf(cmd,
            "INSERT INTO CPS " \
            "(cps_AssemblyID, cps_EUID, cps_CIID, cps_ctp_MSG_ID, cps_cco_MSG_ID, cps_pos1, cps_pos2) " \
            "VALUES ("F_U64", '"F_U64"', "F_CID", "F_U64", "F_U64", "F_S32", "F_S32")",
                     assemblyID,
                     AS_UID_toInteger(cpsID),
                     (uint32)0,
                     ctp,
                     cco,
                     ctgStart,
                     ctgEnd
            );

   return getConnection()->sqlCommand(cmd);
}

bool SQLOutput::commitMDI2DB() {
   return true;
}

bool SQLOutput::commitAFG2DB() {
   return true;
}

bool SQLOutput::commitUTG2DB() {
   return true;
}

bool SQLOutput::commitMPS2DB() {
   return true;
}

bool SQLOutput::commitULK2DB() {
   return true;
}

bool SQLOutput::commitULKList2DB() {
   return true;
}

bool SQLOutput::commitJMP2DB() {
   // we no longer need the ULK structure, dump it
   if (ULK_UID_to_MSGID != NULL) {
      DeleteHashTable_AS(ULK_UID_to_MSGID);
      ULK_UID_to_MSGID = NULL;
   }

   return true;
}

bool SQLOutput::commitJMPList2DB() {
   if (JMP_UID_to_MSGID != NULL) {
      // we no longer need the JMP structure, dump it
      DeleteHashTable_AS(JMP_UID_to_MSGID);
      JMP_UID_to_MSGID = NULL;
   }

   return true;
}

bool SQLOutput::commitCCO2DB() {
   return true;
}

bool SQLOutput::commitCCOMPS2DB() {
   return true;
}

bool SQLOutput::commitUPS2DB() {
   if (UTG_UID_to_MSGID != NULL) {
      // we no longer need the UTG structure, dump it
      DeleteHashTable_AS(UTG_UID_to_MSGID);
      UTG_UID_to_MSGID = NULL;
   }

   return true;
}

bool SQLOutput::commitVAR2DB() {
   return true;
}

bool SQLOutput::commitVARAllele2DB() {
   return true;
}

bool SQLOutput::commitVARAFG2DB() {
   if (AFG_IID_to_MSGID != NULL) {
      //we no longer need the AFG structure, dump it
      DeleteHashTable_AS(AFG_IID_to_MSGID);
      AFG_IID_to_MSGID = NULL;
   }

   if (VAR_UID_to_MSGID != NULL) {
      //we no longer need the VAR structure, dump it
      DeleteHashTable_AS(VAR_UID_to_MSGID);
      VAR_UID_to_MSGID = NULL;
   }

   return true;
}

bool SQLOutput::commitCLK2DB() {
   return true;
}

bool SQLOutput::commitCLKList2DB() {
   return true;
}

bool SQLOutput::commitCLKJMP2DB() {
   return true;
}

bool SQLOutput::commitCLKJMPList2DB() {
   if (CLK_JMP_UID_to_MSGID != NULL) {
      // we no longer need the JMP structure, dump it
      DeleteHashTable_AS(CLK_JMP_UID_to_MSGID);
      CLK_JMP_UID_to_MSGID = NULL;
   }

   if (AFG_UID_to_MSGID != NULL) {
      //we no longer need the AFG structure, dump it
      DeleteHashTable_AS(AFG_UID_to_MSGID);
      AFG_UID_to_MSGID = NULL;
   }

   return true;
}

bool SQLOutput::commitSCF2DB()  {
   return true;
}

bool SQLOutput::commitCTP2DB()  {
   if (SCF_UID_to_MSGID != NULL) {
      // dump the SCF struct
      DeleteHashTable_AS(SCF_UID_to_MSGID);
   }

   return true;
}

bool SQLOutput::commitCTPList2DB() {
   return true;
}

bool SQLOutput::commitCPS2DB() {
   if (CCO_UID_to_MSGID != NULL) {
      // dump the CCO struct
      DeleteHashTable_AS(CCO_UID_to_MSGID);
   }

   if (CTP_UID_to_MSGID != NULL) {
      // dump the CTP struct
      DeleteHashTable_AS(CTP_UID_to_MSGID);
   }

   return true;
}

#endif //SYBASE
