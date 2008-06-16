
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

#include <iostream>
#include <fstream>
#include <assert.h>

#include "BCPOutput.hh"
#include "IDBConnection.hh"
 
using AS_ARD::BCPOutput;
using AS_ARD::SQLOutput;

const char BCPOutput::AFG_FILENAME[MAX_FILE_LEN] = "AFG";
const char BCPOutput::MDI_FILENAME[MAX_FILE_LEN] = "MDI";
const char BCPOutput::UTG_FILENAME[MAX_FILE_LEN] = "UTG";
const char BCPOutput::MPS_FILENAME[MAX_FILE_LEN] = "MPS";
const char BCPOutput::ULK_FILENAME[MAX_FILE_LEN] = "ULK";
const char BCPOutput::ULK_LIST_FILENAME[MAX_FILE_LEN] = "ULK_LIST";
const char BCPOutput::JMP_FILENAME[MAX_FILE_LEN] = "JMP";
const char BCPOutput::JMP_LIST_FILENAME[MAX_FILE_LEN] = "JMP_LIST";
const char BCPOutput::CLK_FILENAME[MAX_FILE_LEN] = "CLK";
const char BCPOutput::CLK_LIST_FILENAME[MAX_FILE_LEN] = "CLK_LIST";
const char BCPOutput::CLK_JMP_FILENAME[MAX_FILE_LEN] = "CLK_JMP";
const char BCPOutput::CLK_JMP_LIST_FILENAME[MAX_FILE_LEN] = "CLK_JMP_LIST";
const char BCPOutput::CCO_FILENAME[MAX_FILE_LEN] = "CCO";
const char BCPOutput::CCO_MPS_FILENAME[MAX_FILE_LEN] = "CCO_MPS";
const char BCPOutput::UPS_FILENAME[MAX_FILE_LEN] = "UPS";
const char BCPOutput::VAR_FILENAME[MAX_FILE_LEN] = "VAR";
const char BCPOutput::VAR_ALLELE_FILENAME[MAX_FILE_LEN] = "VAR_ALLELE";
const char BCPOutput::VAR_AFG_FILENAME[MAX_FILE_LEN] = "VAR_AFG";
const char BCPOutput::SCF_FILENAME[MAX_FILE_LEN] = "SCF";
const char BCPOutput::CTP_FILENAME[MAX_FILE_LEN] = "CTP";
const char BCPOutput::CTP_LIST_FILENAME[MAX_FILE_LEN] = "CTP_LIST";
const char BCPOutput::CPS_FILENAME[MAX_FILE_LEN] = "CPS";

const char BCPOutput::DEFAULT_BCP[MAX_STR_LEN] = "bcp";

BCPOutput::BCPOutput(
            const char * _prefix, 
            const char * _server,
            const char * _database,
            const char * _user,
            const char * _password,
            const char * _bcp) :
      SQLOutput(new Sybase(_server, _user, _password, _database)), 
      mdiBCP(NULL), afgBCP(NULL), utgBCP(NULL), mpsBCP(NULL),
      ulkBCP(NULL), ulkListBCP(NULL), jmpBCP(NULL), jmpListBCP(NULL),
      clkBCP(NULL),clkListBCP(NULL),clkJmpBCP(NULL),clkJmpListBCP(NULL),
      ccoBCP(NULL),ccoMpsBCP(NULL),upsBCP(NULL),
      varBCP(NULL),varAlleleBCP(NULL),varAFGBCP(NULL),
      scfBCP(NULL),ctpBCP(NULL),ctpListBCP(NULL),cpsBCP(NULL) {

   int len = 0;
   assert (_prefix != NULL);
   assert (_server != NULL);
   assert (_user != NULL);
   assert (_password != NULL);
   
   len = strlen(_prefix);
   prefix = new char[len + 1];
   strncpy(prefix, _prefix, len+1);

   len = strlen(_server);
   server = new char[len + 1];
   strncpy(server, _server, len+1);

   len = strlen(_database);
   database = new char[len + 1];
   strncpy(database, _database, len+1);

   len = strlen(_user);
   user = new char[len + 1];
   strncpy(user, _user, len+1);

   len = strlen(_password);
   pass = new char[len + 1];
   strncpy(pass, _password, len+1);
   
   if (_bcp != NULL) { 
      len = strlen(_bcp);
      bcp = new char[len + 1];
      strncpy(bcp, _bcp, len+1);
   }
   else {
      len = strlen(DEFAULT_BCP);
      bcp = new char[len + 1];
      strncpy(bcp, DEFAULT_BCP, len+1);
   }
}

BCPOutput::~BCPOutput() {
   delete[] prefix;
   delete[] server;
   delete[] database;
   delete[] user;
   delete[] pass;
   delete[] bcp;
}

char * BCPOutput::getFileName(const char * fileName) {
   assert(fileName != NULL);
   
   char * path = new char[MAX_PATH_LEN];
   sprintf(path, "%s_%s.bcp", prefix, fileName);

   return path;
}

std::fstream * BCPOutput::openFile(const char * fileName, std::_Ios_Openmode mode) {
   char * path = getFileName(fileName);
   std::fstream * myfile = new std::fstream();   
      
   myfile->open(path, mode);
   assert(myfile->is_open());
    
   delete[] path;
   return myfile;   
}

bool BCPOutput::closeFile(std::fstream **file) {
   (*file)->close();   
   
   delete (*file);
   (*file) = NULL;
   
   return true;
}

bool BCPOutput::updateFile(
         const char * fileName, 
         HashTable_AS * hash, 
         int32 position) {
   std::string line;   
   
   std::ofstream fout;
   std::ifstream fin;
   char * fname = getFileName("tempUpdate");
   char * currName = getFileName(fileName);
   char replace[MAX_STR_LEN];
      
   fin.open(currName, std::ios::in);
   assert(fin.is_open());
   
   fout.open(fname, std::ios::out | std::ios::trunc);
   assert(fout.is_open());
   
   fin.seekg(0);
   while (fin) {
      std::getline(fin, line);
      
      if (!fin.eof()) {
         std::string::size_type pos = 0;
         int counter = 0;
         while (counter < position && (pos = line.find("\t", pos)) != std::string::npos) {
            counter++;
            pos++;
         }
         // do one last find so we get the range to replace
         std::string::size_type end = line.find("\t", pos+1) - pos;
         AS_UID uid = AS_UID_fromInteger(strtoll(line.substr(pos, end).c_str(), NULL, 10));
         sprintf(replace, F_U64, LookupValueInHashTable_AS(hash, AS_UID_toInteger(uid), 0));
         line.replace(pos,end,replace);

         fout << line << std::endl;
      }
   }
   fin.close();
   
   fout.flush();
   fout.close();
   
   rename(fname, currName);
   
   delete[] fname;
   delete[] currName;
   
   return true;
}

bool BCPOutput::runBCP(const char * fileName) {
   return runBCP(fileName, true);
}

bool BCPOutput::runBCP(const char * fileName, bool eraseFile) {
   char * path = getFileName(fileName);
   char command[MAX_STR_LEN];
   char bzip[MAX_STR_LEN];
   bool result = false; 
  
   sprintf(command, "%s %s..%s in %s -c -N -S %s -U %s -P %s -b 10000", bcp, database, fileName, path, server, user, pass);
   result = (system(command) == 0);
   
   if (eraseFile) {
      sprintf(bzip, "rm %s", path);
   }
   else {
      // do harmless command like ls just for filler
      sprintf(bzip, "ls -lh %s", path);
   }
   
   delete[] path;
   
   return result && (system(bzip) == 0);
   //return result;
}

uint64 BCPOutput::storeGenome(
         const char * study,
         const char * project,
         const char * taxon) {
   assert(study != NULL);
   assert(project != NULL);
   assert(taxon != NULL);
   
   return SQLOutput::storeGenome(study, project, taxon);
}
         
uint64 BCPOutput::storeAssembly(
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
   
   return SQLOutput::storeAssembly(assemblyEUID, date, genomeIID, op, genProg, ver, status, notes);
}

bool BCPOutput::storeMDI2DB (
         AS_UID erefines,  
         IntDist_ID irefines,
         float mean,
         float stddev,
         int32 min,
         int32 max) {   
   if (mdiBCP == NULL) {      
      mdiBCP = openFile(MDI_FILENAME);
   }
   
   (*mdiBCP) << assemblyID << "\t" 
         << AS_UID_toInteger(erefines) << "\t"
         << irefines << "\t"
         << mean << "\t"
         << stddev << "\t"
         << min << "\t"
         << max  << "\n";
   
   return true;
}

bool BCPOutput::storeAFG2DB (
         AS_UID eaccession,  
         IntFragment_ID iaccession,
         MateStatType mate_status,
         int32 chimeric,
         int32 chaff,
         CDS_COORD_t bgn,
         CDS_COORD_t end) {
   if (afgBCP == NULL) {
      afgBCP = openFile(AFG_FILENAME);      
   }
   (*afgBCP) << assemblyID << "\t"
         << AS_UID_toInteger(eaccession) << "\t"
         << iaccession << "\t"
         << static_cast<char>(mate_status) << "\t"
         << chimeric << "\t"
         << chaff << "\t"
         << bgn  << "\t"
         << end  << "\n";

   if (AFG_IID_to_MSGID == NULL) {
      AFG_IID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   InsertInHashTable_AS(AFG_IID_to_MSGID, static_cast<uint64>(iaccession), 0, AS_UID_toInteger(eaccession), 0);

   return true;
}

bool BCPOutput::storeUTG2DB (
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
   if (utgBCP == NULL) {
      utgBCP = openFile(UTG_FILENAME);
   }
   
   (*utgBCP) << assemblyID << "\t" 
         << AS_UID_toInteger(eaccession) << "\t"
         << iaccession << "\t"
         << source << "\t"
         << mhp << "\t"
         << coverage_stat << "\t"
         << static_cast<char>(status) << "\t"         
         << "0 \t 0 \t" // what are abp and bbp in CARD spec?
         << length << "\t"                  
         //<< consensus << "\t"
         //<< quality << "\t"
         << " NULL \t NULL \t" // no consensus or quality for now
         << forced << "\t"
         << num_frags << "\n";

   return true;
}

bool BCPOutput::storeMPS2DB (
         AS_UID unitigID,         
         AS_UID afgID,
         FragType type,
         const char * source,
         CDS_COORD_t bgn,
         CDS_COORD_t end,
         int32 delta_length,
         std::string delta) {
   if (mpsBCP == NULL) {
      mpsBCP = openFile(MPS_FILENAME);
   }
     
   (*mpsBCP) << assemblyID << "\t" 
         << AS_UID_toInteger(unitigID) << "\t"         
         << AS_UID_toInteger(afgID) << "\t"
         << static_cast<char>(type) << "\t"
         << source << "\t"
         << bgn << "\t"
         << end << "\t"
         << delta.substr(0,MAX_DELTA).c_str() << "\n";

   return true;
}

bool BCPOutput::storeULK2DB (
         AS_UID euid,
         CDS_CID_t ciid,
         ChunkOrientationType orientation,
         UnitigOverlapType overlap_type,
         int32 is_possible_chimera,
         float mean_distance,
         float std_deviation,
         int32 num_contributing,
         PlacementStatusType status) {
   if (ulkBCP == NULL) {
      ulkBCP = openFile(ULK_FILENAME);
   }
   
   (*ulkBCP) 
         << assemblyID << "\t"   
         << AS_UID_toInteger(euid) << "\t"
         << ciid << "\t"
         << static_cast<char>(orientation) << "\t"
         << static_cast<char>(overlap_type) << "\t"
         << is_possible_chimera << "\t"
         //smp->includes_guide, - why no guide included?
         << mean_distance << "\t"
         << std_deviation << "\t"
         << num_contributing << "\t"
         << static_cast<char>(status) << "\n";

   return true;
}

bool BCPOutput::storeLKList2DB(int jmpType, AS_UID utgID, AS_UID ulkID) {
   std::fstream * stream = NULL;
      
   if (jmpType == ULK_TYPE) {
      if (ulkListBCP == NULL) {
         ulkListBCP = openFile(ULK_LIST_FILENAME);
      }

      stream = ulkListBCP;
   }
   else if (jmpType == CLK_TYPE) {
      if (clkListBCP == NULL) {
         clkListBCP = openFile(CLK_LIST_FILENAME);
      }

      stream = clkListBCP;
   }
   else {
      assert(0);
   }

   (*stream) << assemblyID << "\t" 
         << AS_UID_toInteger(utgID) << "\t"
         << AS_UID_toInteger(ulkID) << "\n";

   return true;
}

bool BCPOutput::storeJMP2DB(int jmpType, AS_UID jmpID, AS_UID ulkID, LinkType type) {
   std::fstream * stream = NULL;
   
   if (jmpType == ULK_TYPE) {
      if (jmpBCP == NULL) {
         jmpBCP = openFile(JMP_FILENAME);
      }

      stream = jmpBCP;
   }
   else if (jmpType == CLK_TYPE) {
      if (clkJmpBCP == NULL) {
         clkJmpBCP = openFile(CLK_JMP_FILENAME);
      }

      stream = clkJmpBCP;
   }
   else {
      assert(0);
   }
   
   //TODO: warning using 0 as ciid for JMP
   (*stream) 
         << assemblyID  << "\t"
         << AS_UID_toInteger(jmpID) << "\t"
         << 0 << "\t"
         << AS_UID_toInteger(ulkID) << "\t";

   if (jmpType == ULK_TYPE) {
      (*stream) << "X" << "\t"; //what is status?
   }
   (*stream) << static_cast<char>(type) << "\n";
   
   return true;
}

bool BCPOutput::storeJMPList2DB(int jmpType, AS_UID jmpListID, AS_UID jmpID, AS_UID fragID) {
   std::fstream * stream = NULL;
   
   if (jmpType == ULK_TYPE) {
      if (jmpListBCP == NULL) {
         jmpListBCP = openFile(JMP_LIST_FILENAME);
      }

      stream = jmpListBCP;
   }
   else if (jmpType == CLK_TYPE) {
      if (clkJmpListBCP == NULL) {
         clkJmpListBCP = openFile(CLK_JMP_LIST_FILENAME);
      }

      stream = clkJmpListBCP;
   }
   else {
      assert(0);
   }

   //TODO: warning using 0 as ciid for JMP_LIST
   (*stream) << assemblyID << "\t" 
         << AS_UID_toInteger(jmpListID) << "\t"
         << 0 << "\t"
         << AS_UID_toInteger(jmpID) << "\t"
         << AS_UID_toInteger(fragID) << "\n";

   return true;
}

bool BCPOutput::storeCCO2DB (
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

   if (ccoBCP == NULL) {
      ccoBCP = openFile(CCO_FILENAME);
   }
   
   (*ccoBCP) << assemblyID << "\t" 
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

bool BCPOutput::storeCCOMPS2DB(
                  AS_UID ccoMpsID,
                  AS_UID ccoID,            
                  AS_UID fragID,
                  FragType type,
                  const char * source,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {
   if (ccoMpsBCP == NULL) {
      ccoMpsBCP = openFile(CCO_MPS_FILENAME);
   }

   //TODO: warning using 0 as ciid for CCOMPS
   (*ccoMpsBCP) << assemblyID << "\t"
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

bool BCPOutput::storeUPS2DB(
                  AS_UID upsID,
                  AS_UID ccoID,            
                  AS_UID unitigID,
                  UnitigType type,
                  CDS_COORD_t bgn,
                  CDS_COORD_t end,
                  int32 delta_length,
                  std::string delta) {
   if (upsBCP == NULL) {
      upsBCP = openFile(UPS_FILENAME);
   }
      
   //TODO: warning truncating delta in UPS to 1000
   //using 0 as ciid for UPS
   (*upsBCP) << assemblyID << "\t" 
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

bool BCPOutput::storeVAR2DB(
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
   if (varBCP == NULL) {
      varBCP = openFile(VAR_FILENAME);
   }

   //TODO: warning using 0 as ciid for VAR
   (*varBCP) << assemblyID << "\t"
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

bool BCPOutput::storeVARAllele2DB(AS_UID varAlleleID, AS_UID varID, uint32 nra, uint32 wgt, std::string seq) {
   if (varAlleleBCP == NULL) {
      varAlleleBCP = openFile(VAR_ALLELE_FILENAME);
   }
   
   //TODO: warning using 0 as ciid for VAR_ALLELE
   (*varAlleleBCP) << assemblyID << "\t"
                     << AS_UID_toInteger(varAlleleID) << "\t"
                     << 0 << "\t" 
                     << AS_UID_toInteger(varID) << "\t"
                     << nra << "\t"
                     << wgt << "\t"
                     << seq.c_str() << "\n";

   return true;
}

bool BCPOutput::storeVARAFG2DB(AS_UID varAfgID, AS_UID varID, CDS_CID_t readID) {
   if (varAFGBCP == NULL) {
      varAFGBCP = openFile(VAR_AFG_FILENAME);
   }
   
   //TODO: warning using 0 as ciid for VAR_AFG
   (*varAFGBCP) << assemblyID << "\t"
                  << AS_UID_toInteger(varAfgID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(varID) << "\t"
                  << readID << "\n";

   return true;
}

bool BCPOutput::storeCLK2DB(
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
   if (clkBCP == NULL) {
      clkBCP = openFile(CLK_FILENAME);
   }

   (*clkBCP) << assemblyID << "\t"
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

bool BCPOutput::storeSCF2DB(AS_UID eaccession, CDS_CID_t iaccession, uint32 num_contig_pairs) {
   if (scfBCP == NULL) {
      scfBCP = openFile(SCF_FILENAME);
   }
   
   (*scfBCP) << assemblyID << "\t"
                  << AS_UID_toInteger(eaccession) << "\t"
                  << iaccession << "\t"
                  << num_contig_pairs << "\n";
   
   return true;
}

bool BCPOutput::storeCTP2DB(AS_UID ctpID, AS_UID scfID, float mean, float stddev, ChunkOrientationType orient) {
   if (ctpBCP == NULL) {
      ctpBCP = openFile(CTP_FILENAME);
   }

   //TODO: warning using 0 as ciid for CTP
   (*ctpBCP) << assemblyID << "\t"
                  << AS_UID_toInteger(ctpID) << "\t"
                  << 0 << "\t"
                  << AS_UID_toInteger(scfID) << "\t"
                  << mean << "\t"
                  << stddev << "\t"
                  << static_cast<char>(orient) << "\n";
   
   return true;
}
         
bool BCPOutput::storeCTPList2DB(AS_UID ctpListID, AS_UID ctpID, AS_UID ccoID) {
   if (ctpListBCP == NULL) {
      ctpListBCP = openFile(CTP_LIST_FILENAME);
   }
   
   (*ctpListBCP) << assemblyID << "\t"
                  << AS_UID_toInteger(ctpListID) << "\t"
                  << (uint32)0 << "\t"
                  << AS_UID_toInteger(ctpID) << "\t"
                  << AS_UID_toInteger(ccoID) << "\n";

   return true;
}

bool BCPOutput::storeCPS2DB(AS_UID cpsID, AS_UID ctpID, AS_UID ccoID, CDS_COORD_t ctgStart, CDS_COORD_t ctgEnd) {
   if (cpsBCP == NULL) {
      cpsBCP = openFile(CPS_FILENAME);
   }
   
   (*cpsBCP) << assemblyID << "\t"
               << AS_UID_toInteger(cpsID) << "\t"
               << (uint32)0 << "\t"
               << AS_UID_toInteger(ctpID) << "\t"
               << AS_UID_toInteger(ccoID) << "\t"
               << ctgStart << "\t"
               << ctgEnd << "\n";

   return true;
}

bool BCPOutput::commitMDI2DB() {
   if (mdiBCP == NULL) { return true; }
      
   SQLOutput::commitMDI2DB();
   bool result = closeFile(&mdiBCP) && runBCP(MDI_FILENAME);
   
   return result;
}

bool BCPOutput::commitAFG2DB() {
   if (afgBCP == NULL) { return true; }
   
   SQLOutput::commitAFG2DB();
   if (AFG_UID_to_MSGID == NULL) {
      AFG_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }
   bool result = closeFile(&afgBCP) && runBCP(AFG_FILENAME);
   result = result && getConnection()->populateHash(AFG_UID_to_MSGID, "afg_EUID", "afg_MSG_ID", "AFG", assemblyID);
      
   return result;
}

bool BCPOutput::commitUTG2DB() {
   if (utgBCP == NULL) { return true; }

   SQLOutput::commitUTG2DB();
   
   if (UTG_UID_to_MSGID == NULL) {
      UTG_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   bool result = closeFile(&utgBCP) && runBCP(UTG_FILENAME);
   result = result && getConnection()->populateHash(UTG_UID_to_MSGID, "utg_EUID", "utg_MSG_ID", "UTG", assemblyID);
      
   return result;
}

bool BCPOutput::commitMPS2DB() {
   if (mpsBCP == NULL) { return true; }
   
   SQLOutput::commitMPS2DB();
   
   bool result = closeFile(&mpsBCP);
   updateFile(MPS_FILENAME, UTG_UID_to_MSGID, 1);
   updateFile(MPS_FILENAME, AFG_UID_to_MSGID, 2);
      
   return result && runBCP(MPS_FILENAME);
}

bool BCPOutput::commitULK2DB() {
   if (ulkBCP == NULL) { return true; }

   SQLOutput::commitULK2DB();
   
   if (ULK_UID_to_MSGID == NULL) {
      ULK_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   bool result = closeFile(&ulkBCP) && runBCP(ULK_FILENAME);   
   result = result && getConnection()->populateHash(ULK_UID_to_MSGID, "ulk_EUID", "ulk_MSG_ID", "ULK", assemblyID);   
   
   return result;   
}

bool BCPOutput::commitULKList2DB() {
   if (ulkListBCP == NULL) { return true; }

   SQLOutput::commitULKList2DB();

   bool result = closeFile(&ulkListBCP);
   updateFile(ULK_LIST_FILENAME, UTG_UID_to_MSGID, 1);
   updateFile(ULK_LIST_FILENAME, ULK_UID_to_MSGID, 2);
      
   return result && runBCP(ULK_LIST_FILENAME);
}

bool BCPOutput::commitJMP2DB() {
   if (jmpBCP == NULL) { return true; }
   
   if (JMP_UID_to_MSGID == NULL) {
      JMP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   bool result = closeFile(&jmpBCP);
   updateFile(JMP_FILENAME, ULK_UID_to_MSGID, 3);

   SQLOutput::commitJMP2DB();
      
   result = result && runBCP(JMP_FILENAME);
   return result && getConnection()->populateHash(JMP_UID_to_MSGID, "jmp_EUID", "jmp_MSG_ID", "JMP", assemblyID);
}
        
bool BCPOutput::commitJMPList2DB() {
   if (jmpListBCP == NULL) { return true; }

   bool result = closeFile(&jmpListBCP);
   updateFile(JMP_LIST_FILENAME, JMP_UID_to_MSGID, 3);
   updateFile(JMP_LIST_FILENAME, AFG_UID_to_MSGID, 4);

   SQLOutput::commitJMPList2DB();
   return result && runBCP(JMP_LIST_FILENAME);
}

bool BCPOutput::commitCCO2DB() {
   if (ccoBCP == NULL) { return true; }

   SQLOutput::commitCCO2DB();

   if (CCO_UID_to_MSGID == NULL) {
      CCO_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   bool result = closeFile(&ccoBCP) && runBCP(CCO_FILENAME);
   result = result && getConnection()->populateHash(CCO_UID_to_MSGID, "cco_EUID", "cco_MSG_ID", "CCO", assemblyID);
      
   return result;
   
}
         
bool BCPOutput::commitCCOMPS2DB() {
   if (ccoMpsBCP == NULL) { return true; }

   SQLOutput::commitCCOMPS2DB();
   
   bool result = closeFile(&ccoMpsBCP);
   updateFile(CCO_MPS_FILENAME, CCO_UID_to_MSGID, 3);
   updateFile(CCO_MPS_FILENAME, AFG_UID_to_MSGID, 4);
      
   result = result && runBCP(CCO_MPS_FILENAME);
}

bool BCPOutput::commitUPS2DB() { 
   if (upsBCP == NULL) { return true; }

   bool result = closeFile(&upsBCP);
   updateFile(UPS_FILENAME, CCO_UID_to_MSGID, 3);
   updateFile(UPS_FILENAME, UTG_UID_to_MSGID, 4);
      
   SQLOutput::commitUPS2DB();
   
   return result && runBCP(UPS_FILENAME);
}

bool BCPOutput::commitVAR2DB() { 
   if (varBCP == NULL) { return true; }
   
   SQLOutput::commitVAR2DB();

   if (VAR_UID_to_MSGID == NULL) {
      VAR_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   bool result = closeFile(&varBCP);
   updateFile(VAR_FILENAME, CCO_UID_to_MSGID, 3);
   result = result && runBCP(VAR_FILENAME);
   result = result && getConnection()->populateHash(VAR_UID_to_MSGID, "var_EUID", "var_MSG_ID", "VAR", assemblyID);
      
   return result; 
}

bool BCPOutput::commitVARAllele2DB() { 
   if (varAlleleBCP == NULL) { return true; }
      
   bool result = closeFile(&varAlleleBCP);
   updateFile(VAR_ALLELE_FILENAME, VAR_UID_to_MSGID, 3);
   SQLOutput::commitVARAllele2DB();
   
   return result && runBCP(VAR_ALLELE_FILENAME);
}

bool BCPOutput::commitVARAFG2DB() {
   if (varAFGBCP == NULL) { return true; }
      
   bool result = closeFile(&varAFGBCP);
   updateFile(VAR_AFG_FILENAME, VAR_UID_to_MSGID, 3);
   updateFile(VAR_AFG_FILENAME, AFG_IID_to_MSGID, 4);
   updateFile(VAR_AFG_FILENAME, AFG_UID_to_MSGID, 4);
   SQLOutput::commitVARAFG2DB();
   
   return result && runBCP(VAR_AFG_FILENAME);
}

bool BCPOutput::commitCLK2DB() {
   if (clkBCP == NULL) { return true; }

   if (CLK_UID_to_MSGID == NULL) {
      CLK_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   SQLOutput::commitCLK2DB();
   bool result = closeFile(&clkBCP) && runBCP(CLK_FILENAME);
   result = result && getConnection()->populateHash(CLK_UID_to_MSGID, "clk_EUID", "clk_MSG_ID", "CLK", assemblyID);
   
   return result;
}

bool BCPOutput::commitCLKList2DB() {
   if (clkListBCP == NULL) { return true; }

   bool result = closeFile(&clkListBCP);
   updateFile(CLK_LIST_FILENAME, CCO_UID_to_MSGID, 1);
   updateFile(CLK_LIST_FILENAME, CLK_UID_to_MSGID, 2);
   
   SQLOutput::commitCLKList2DB();
   return result && runBCP(CLK_LIST_FILENAME);
}

bool BCPOutput::commitCLKJMP2DB() {
   if (clkJmpBCP == NULL) { return true; }

   if (CLK_JMP_UID_to_MSGID == NULL) {
      CLK_JMP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   bool result = closeFile(&clkJmpBCP);
   updateFile(CLK_JMP_FILENAME, CLK_UID_to_MSGID, 3);
   
   SQLOutput::commitCLKJMP2DB();
   result = result && runBCP(CLK_JMP_FILENAME);
   return result && getConnection()->populateHash(CLK_JMP_UID_to_MSGID, "clk_jmp_EUID", "clk_jmp_MSG_ID", "CLK_JMP", assemblyID);
}
        
bool BCPOutput::commitCLKJMPList2DB() {
   if (clkJmpListBCP == NULL) { return true; }

   bool result = closeFile(&clkJmpListBCP);
   updateFile(CLK_JMP_LIST_FILENAME, CLK_JMP_UID_to_MSGID, 3);
   updateFile(CLK_JMP_LIST_FILENAME, AFG_UID_to_MSGID, 4);

   SQLOutput::commitCLKJMPList2DB();
   return result && runBCP(CLK_JMP_LIST_FILENAME);
}

bool BCPOutput::commitSCF2DB()  { 
   if (scfBCP == NULL) { return true; }

   SQLOutput::commitSCF2DB();
   if (SCF_UID_to_MSGID == NULL) {
      SCF_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }   
   
   bool result = closeFile(&scfBCP) && runBCP(SCF_FILENAME, false);
   result = result && getConnection()->populateHash(SCF_UID_to_MSGID, "scf_EUID", "scf_MSG_ID", "SCF", assemblyID);
   
   return result;
}

bool BCPOutput::commitCTP2DB()  {
   if (ctpBCP == NULL) { return true; }

   if (CTP_UID_to_MSGID == NULL) {
      CTP_UID_to_MSGID = CreateScalarHashTable_AS(32 * 1024);
   }

   bool result = closeFile(&ctpBCP);
   updateFile(CTP_FILENAME, SCF_UID_to_MSGID, 3);
   
   SQLOutput::commitCTP2DB();
   result = result && runBCP(CTP_FILENAME, false);
   return result && getConnection()->populateHash(CTP_UID_to_MSGID, "ctp_EUID", "ctp_MSG_ID", "CTP", assemblyID);
}

bool BCPOutput::commitCTPList2DB() {
   if (ctpListBCP == NULL) { return true; }

   bool result = closeFile(&ctpListBCP);
   updateFile(CTP_LIST_FILENAME, CTP_UID_to_MSGID, 3);
   updateFile(CTP_LIST_FILENAME, CCO_UID_to_MSGID, 4);

   SQLOutput::commitCTPList2DB();
   return result && runBCP(CTP_LIST_FILENAME, false);
}

bool BCPOutput::commitCPS2DB() {
   if (cpsBCP == NULL) { return true; }

   bool result = closeFile(&cpsBCP);
   updateFile(CPS_FILENAME, CTP_UID_to_MSGID, 3);
   updateFile(CPS_FILENAME, CCO_UID_to_MSGID, 4);

   SQLOutput::commitCPS2DB();
   return result && runBCP(CPS_FILENAME, false);
}

#endif //SYBASE
