
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
#ifndef AS_ARD_database_HH
#define AS_ARD_database_HH

#include <string>

#include "IAssemblyDB.hh"

extern "C" {
   #include "AS_global.h"
   #include "SYS_UIDclient.h"
   #include "AS_MSG_pmesg.h"
   #include "AS_UTL_Hash.h"
}

namespace AS_ARD {
   class AS_ARD_database : public IAssemblyDB {
      private:
         // disable copy constructor
         AS_ARD_database(AS_ARD_database &);

         bool addMDI2DB     (SnapMateDistMesg * smdm);
         bool addAFG2DB     (AugFragMesg * afg);
         bool addUTG2DB     (SnapUnitigMesg * sum);
         bool addMPS2DB     (AS_UID unitigID, SnapMultiPos * smp);
         bool addULK2DB     (SnapUnitigLinkMesg * ulk);
         bool addLKList2DB  (int type, AS_UID ulkID, AS_UID utgID);
         bool addJMP2DB     (int type, AS_UID ulkID, SnapMate_Pairs * mp);
         bool addJMPList2DB (int type, AS_UID jmpID, AS_UID fragID);
         bool addCCO2DB     (SnapConConMesg * cco);
         bool addCCOMPS2DB  (AS_UID ccoID, SnapMultiPos * piece);
         bool addUPS2DB     (AS_UID ccoID, UnitigPos * unitig);
         bool addVAR2DB     (AS_UID ccoID, IntMultiVar * var);
         bool addVARAllele2DB(AS_UID varID, uint64 nra, uint64 wgt, std::string seq);
         bool addVARAFG2DB  (AS_UID varID, CDS_CID_t readID);

         bool addCLK2DB     (SnapContigLinkMesg *);
         bool addSCF2SDB    (SnapScaffoldMesg *);
         bool addCTP2DB     (AS_UID scfID, SnapContigPairs * ctp, CDS_COORD_t &scfLen);
         bool addCTPList2DB (AS_UID ctpID, AS_UID ccoID);
         bool addCPS2DB     (AS_UID ctpID, AS_UID ccoID, bool isReversed, CDS_COORD_t &scfLen);

         bool addGenericMesg2DB(GenericMesg * gen);

         void fixSource(std::string &source);

         HashTable_AS * contigLens;
      public:
         AS_ARD_database(IDBOutput * output);
         ~AS_ARD_database();
         bool LoadDatabaseFromASMFile(FILE * fi, UIDserver * uid);
   };
};

#endif // AS_ARD_database_HH
