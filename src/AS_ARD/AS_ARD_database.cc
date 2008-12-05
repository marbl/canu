
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

static const char *rcsid = "$Id: AS_ARD_database.cc,v 1.9 2008-12-05 19:06:11 brianwalenz Exp $";

#include <iostream>
#include <string>

extern "C" {
   #include <math.h>
   #include <assert.h>
}

#include "AS_ARD_database.hh"

//#define DEBUG_AS_ARD
#define MAX_DELTA_CHARS 10

using AS_ARD::AS_ARD_database;

AS_ARD_database::AS_ARD_database(IDBOutput * output) :
      IAssemblyDB(output) {
   contigLens = CreateScalarHashTable_AS();
}

AS_ARD_database::~AS_ARD_database() {
   std::cerr << "Deleting db" << std::endl;
   DeleteHashTable_AS(contigLens);
}

void AS_ARD_database::fixSource(std::string &source) {
   std::string::size_type pos = 0;
   while ((pos = source.find("\n", pos)) != std::string::npos) {
      source.replace(pos, 1, " ");
   }
}

bool AS_ARD_database::addMDI2DB(SnapMateDistMesg * smdm) {
   return db->storeMDI2DB(
       smdm->erefines,
       smdm->irefines,
       smdm->mean,
       smdm->stddev,
       smdm->min,
       smdm->max);
}


bool AS_ARD_database::addAFG2DB(AugFragMesg * afg) {
   return db->storeAFG2DB(
         afg->eaccession,
         afg->iaccession,
         afg->mate_status,
         afg->chaff,
         afg->clear_rng.bgn,
         afg->clear_rng.end);
}

bool AS_ARD_database::addUTG2DB(SnapUnitigMesg * utg) {
   int i = 0;
   std::string source;

   bool result = db->storeUTG2DB(
            utg->eaccession,
            utg->iaccession,
            source.c_str(),
            utg->microhet_prob,
            utg->coverage_stat,
            utg->status,
            //what are abp and bbp in CARD spec?
            utg->length,
            utg->consensus,
            utg->quality,
            utg->forced,
            utg->num_frags);

   for (i = 0; i < utg->num_frags; i++) {
      result = result && addMPS2DB(utg->eaccession, &utg->f_list[i]);
   }

   return result;
}

bool AS_ARD_database::addMPS2DB(AS_UID unitigID, SnapMultiPos * mps) {
  int i, pos;
  std::string delta;
  char currDelta[MAX_DELTA_CHARS];
  std::string source;

#warning source is unused in MPS
  source[0] = '\0';

   for (i = 0; i < mps->delta_length; i++) {
      sprintf(currDelta, "%d ", mps->delta[i]);
      delta.append(currDelta);
   }

   return db->storeMPS2DB(
            unitigID,
            mps->eident,
            mps->type,
            source.c_str(),
            mps->position.bgn,
            mps->position.end,
            mps->delta_length,
            delta);
}

bool AS_ARD_database::addULK2DB(SnapUnitigLinkMesg * ulk) {
   int i;
   AS_UID euid = AS_UID_fromInteger(getUID(this->uids));
   CDS_CID_t ciid = 0; // how can we get IID?

   bool result = db->storeULK2DB(
            euid,
            ciid,
            ulk->orientation,
            ulk->overlap_type,
            ulk->is_possible_chimera,
            //smp->includes_guide, - why no guide included?
            ulk->mean_distance,
            ulk->std_deviation,
            ulk->num_contributing,
            ulk->status);

   int32 num = ulk->num_contributing;
   if (ulk->overlap_type != AS_NO_OVERLAP) {
      num--;
   }

   for (i = 0; i < num; i++) {
      result = result && addJMP2DB(IDBOutput::ULK_TYPE, euid, &ulk->jump_list[i]);
   }

   result = result && addLKList2DB(IDBOutput::ULK_TYPE, ulk->eunitig1, euid);
   result = result && addLKList2DB(IDBOutput::ULK_TYPE, ulk->eunitig2, euid);

   return result;
}

bool AS_ARD_database::addLKList2DB(int type, AS_UID childID, AS_UID parentID) {
   return db->storeLKList2DB(type, childID, parentID);
}

bool AS_ARD_database::addJMP2DB(int type, AS_UID lkID, SnapMate_Pairs * mp) {
   AS_UID jmpID = AS_UID_fromInteger(getUID(this->uids));

   bool result = db->storeJMP2DB(type, jmpID, lkID, mp->type);

   result = result && addJMPList2DB(type, jmpID, mp->in1);
   result = result && addJMPList2DB(type, jmpID, mp->in2);

   return result;
}

bool AS_ARD_database::addJMPList2DB(int type, AS_UID jmpID, AS_UID fragID) {
   return db->storeJMPList2DB(type, AS_UID_fromInteger(getUID(this->uids)), jmpID, fragID);
}

bool AS_ARD_database::addCCO2DB(SnapConConMesg * cco)
{
   int i = 0;
   std::string::size_type pos = 0;
   std::string cons = cco->consensus;
   uint64 length = cco->length;

   while ((pos = cons.find("-", pos)) != std::string::npos) {
      pos++;
      length--;
   }
   InsertInHashTable_AS(contigLens, AS_UID_toInteger(cco->eaccession), 0, length, 0);

   bool result = db->storeCCO2DB(
            cco->eaccession,
            cco->iaccession,
            cco->placed,
            cco->length,
            cco->consensus,
            cco->quality,
            cco->forced,
            cco->num_pieces,
            cco->num_unitigs,
            cco->num_vars);

   for (i = 0; i < cco->num_pieces; i++) {
      result = result && addCCOMPS2DB(cco->eaccession, &cco->pieces[i]);
   }

   for (i = 0; i < cco->num_unitigs; i++) {
      result = result && addUPS2DB(cco->eaccession, &cco->unitigs[i]);
   }

   for (i = 0; i < cco->num_vars; i++) {
      result = result && addVAR2DB(cco->eaccession, &cco->vars[i]);
   }

   return result;
}

bool AS_ARD_database::addCCOMPS2DB(AS_UID ccoID, SnapMultiPos * piece) {
  int i;
  std::string delta;
  char currDelta[MAX_DELTA_CHARS];
  std::string source;


#warning source is unused in mps
  source[0] = '\0';

   for (i = 0; i < piece->delta_length; i++) {
      sprintf(currDelta, "%d ", piece->delta[i]);
      delta.append(currDelta);
   }

   return db->storeCCOMPS2DB(
            AS_UID_fromInteger(getUID(this->uids)),
            ccoID,
            piece->eident,
            piece->type,
            source.c_str(),
            piece->position.bgn,
            piece->position.end,
            piece->delta_length,
            delta);
}


bool AS_ARD_database::addUPS2DB(AS_UID ccoID, UnitigPos * unitig) {
   int i, pos;
   std::string delta;
   char currDelta[MAX_DELTA_CHARS];

   for (i = 0; i < unitig->delta_length; i++) {
      sprintf(currDelta, "%d ", unitig->delta[i]);
      delta.append(currDelta);
   }

   return db->storeUPS2DB(
            AS_UID_fromInteger(getUID(this->uids)),
            ccoID,
            unitig->eident,
            unitig->type,
            unitig->position.bgn,
            unitig->position.end,
            unitig->delta_length,
            delta);
}

bool AS_ARD_database::addVAR2DB(AS_UID ccoID, IntMultiVar * var) {
   bool result = false;
   int i = 0;
   AS_UID varID = AS_UID_fromInteger(getUID(this->uids));

   result = db->storeVAR2DB(
            varID,
            ccoID,
            var->position.bgn,
            var->position.end,
            var->num_reads,
            var->num_conf_alleles,
            var->min_anchor_size,
            var->var_length,
            var->curr_var_id,
            var->phased_var_id
            );

   std::string nra = var->nr_conf_alleles;
   std::string wgs = var->weights;
   std::string seq = var->var_seq;

   // trim the \ns off the string
   nra = nra.substr(0, nra.length() - 1);
   wgs = wgs.substr(0, wgs.length() - 1);
   seq = seq.substr(0, seq.length() - 1);

   std::size_t posInNra = 0;
   std::size_t posInWgt = 0;
   std::size_t posInSeq = 0;

   for (i = 0; i < var->num_conf_alleles; i++) {
      std::size_t nraEnd = nra.find("/", posInNra) - posInNra;
      std::size_t wgtEnd = wgs.find("/", posInWgt) - posInWgt;
      std::size_t seqEnd = seq.find("/", posInSeq) - posInSeq;

      result = result & addVARAllele2DB(
                           varID,
                           atoi(nra.substr(posInNra, nraEnd).c_str()),
                           atoi(wgs.substr(posInWgt, wgtEnd).c_str()),
                           seq.substr(posInSeq, seqEnd));

      posInNra += nraEnd+1;
      posInWgt += wgtEnd+1;
      posInSeq += seqEnd+1;
   }

   std::string readIIDs = var->conf_read_iids;
   std::size_t pos = 0;
   for (i = 0; i < var->num_reads; i++) {
      std::size_t end = readIIDs.find("/", pos) - pos;

      result = result & addVARAFG2DB(varID, AS_IID_fromString((char *)readIIDs.substr(pos, end).c_str(), NULL));
      pos += end + 1;
   }

   return result;
}

bool AS_ARD_database::addVARAllele2DB(AS_UID varID, uint64 nra, uint64 wgt, std::string seq) {
   return db->storeVARAllele2DB(AS_UID_fromInteger(getUID(this->uids)), varID, nra, wgt, seq);
}

bool AS_ARD_database::addVARAFG2DB(AS_UID varID, CDS_CID_t readID) {
   return db->storeVARAFG2DB(AS_UID_fromInteger(getUID(this->uids)), varID, readID);
}

bool AS_ARD_database::addCLK2DB(SnapContigLinkMesg * clk) {
   int i;
   AS_UID euid = AS_UID_fromInteger(getUID(this->uids));
   CDS_CID_t ciid = 0; // how can we get IID?

   // should generated EUID and IID (how do we generate IID?)
   bool result = db->storeCLK2DB(
            euid,
            ciid,
            clk->orientation,
            clk->overlap_type,
            clk->is_possible_chimera,
            clk->includes_guide,
            clk->mean_distance,
            clk->std_deviation,
            clk->num_contributing,
            clk->status);

   int32 num = clk->num_contributing;
   if (clk->overlap_type != AS_NO_OVERLAP) {
      num--;
   }

   for (i = 0; i < num; i++) {
      result = result && addJMP2DB(IDBOutput::CLK_TYPE, euid, &clk->jump_list[i]);
   }

   result = result && addLKList2DB(IDBOutput::CLK_TYPE, clk->econtig1, euid);
   result = result && addLKList2DB(IDBOutput::CLK_TYPE, clk->econtig2, euid);

   return result;
}

bool AS_ARD_database::addSCF2SDB(SnapScaffoldMesg * scf) {
   bool result = false;
   int i = 0;
   uint32 contigPairs = 0;
   CDS_COORD_t scfLen = 0;

   result = db->storeSCF2DB(scf->eaccession, scf->iaccession, scf->num_contig_pairs);

   // we still want to add contigs for those scaffolds that have "0" pairs so increment to 1
   contigPairs = (scf->num_contig_pairs == 0 ? 1 : scf->num_contig_pairs);

   for (i = 0; i < contigPairs; i++) {
      result = result & addCTP2DB(scf->eaccession, &scf->contig_pairs[i], scfLen);
   }


   return result;
}

bool AS_ARD_database::addCTP2DB(AS_UID scfID, SnapContigPairs * ctp, CDS_COORD_t &scfLen) {
   bool result = false;
   AS_UID ctpID = AS_UID_fromInteger(getUID(this->uids));
   float stddev = ctp->stddev;
   if (ctp->stddev != ctp->stddev) {
      stddev = 0;
      std::cerr << "WARNING: stddev was nan, replaced with 0" << std::endl;
   }

   result = db->storeCTP2DB(ctpID, scfID, ctp->mean, stddev, ctp->orient);
   result = result & addCTPList2DB(ctpID, ctp->econtig1);

   // only add the first ctp if it is the only one or the first one
   // otherwise, the second is always a repeat of the previous one
   if (scfLen == 0 || AS_UID_compare(ctp->econtig1, ctp->econtig2) == 0) {
      result = result & addCPS2DB(ctpID, ctp->econtig1, (ctp->orient == AS_OUTTIE || ctp->orient == AS_ANTI) ? true : false, scfLen);
   }

   // if there is a gap and this is not a single-contig pair
   if (AS_UID_compare(ctp->econtig1, ctp->econtig2) != 0) {
      scfLen += (CDS_COORD_t)rintf(ctp->mean);

      result = result & addCTPList2DB(ctpID, ctp->econtig2);
      result = result & addCPS2DB(ctpID, ctp->econtig2, (ctp->orient == AS_INNIE || ctp->orient == AS_ANTI) ? true : false, scfLen);
   }

   return result;
}

bool AS_ARD_database::addCTPList2DB(AS_UID ctpID, AS_UID ccoID) {
   return db->storeCTPList2DB(AS_UID_fromInteger(getUID(this->uids)), ctpID, ccoID);
}

bool AS_ARD_database::addCPS2DB(AS_UID ctpID, AS_UID ccoID, bool isReversed, CDS_COORD_t &scfLen) {
   bool result = false;
   CDS_COORD_t ctgStart = scfLen;
   CDS_COORD_t ctgEnd = scfLen + LookupValueInHashTable_AS(contigLens, AS_UID_toInteger(ccoID), 0);

   scfLen = ctgEnd;

   if (isReversed) {
      result = db->storeCPS2DB(AS_UID_fromInteger(getUID(this->uids)), ctpID, ccoID, ctgEnd, ctgStart);
   } else {
      result = db->storeCPS2DB(AS_UID_fromInteger(getUID(this->uids)), ctpID, ccoID, ctgStart, ctgEnd);
   }
   return result;
}

bool AS_ARD_database::addGenericMesg2DB(GenericMesg * gen) {
  bool result = false;

  // by convention, in the ASM file all messages of one type are together so we commit when we get to the next type
  switch(gen->t)
  {
    case MESG_MDI:
      result = addMDI2DB(static_cast<SnapMateDistMesg *>(gen->m));
      numMDI++;
      break;
    case MESG_AFG:
      db->commitMDI2DB();

      result = addAFG2DB(static_cast<AugFragMesg *> (gen->m));
      numAFG++;
      break;
    case MESG_UTG:
      db->commitAFG2DB();

      result = addUTG2DB(static_cast<SnapUnitigMesg *>(gen->m));
      numUTG++;
      break;
    case MESG_ULK:
      db->commitUTG2DB();
      db->commitMPS2DB();

      result = addULK2DB(static_cast<SnapUnitigLinkMesg *>(gen->m));
      numULK++;
      break;
    case MESG_CCO:
      db->commitULK2DB();
      db->commitULKList2DB();
      db->commitJMP2DB();
      db->commitJMPList2DB();

      result = addCCO2DB(static_cast<SnapConConMesg *>(gen->m));
      numCCO++;
      break;
    case MESG_CLK:
      db->commitCCO2DB();

      result = addCLK2DB(static_cast<SnapContigLinkMesg *>(gen->m));
      numCLK++;
      break;
    case MESG_SCF:
      db->commitCLK2DB();
      db->commitCLKList2DB();
      db->commitCLKJMP2DB();
      db->commitCLKJMPList2DB();

      result = addSCF2SDB(static_cast<SnapScaffoldMesg *>(gen->m));
      numSCF++;
      break;
    case MESG_SLK:
      // SLK not processed right now
      std::cerr << "WARNING: skipping SLK message" << std::endl;
      result = true;
      break;
    default:
      std::cerr << "WARNING: ***********************UNKONWN MESSAGE: " << gen->t << std::endl;
      result = true;
      break;
  }

  return result;
}

bool AS_ARD_database::LoadDatabaseFromASMFile(FILE * fi, UIDserver * uids)
{
   GenericMesg * gen;
   unsigned long mesgCount = 0;
   this->uids = uids;

   assert(fi != NULL);

   std::cerr.setf(std::ios::fixed, std::ios::floatfield);
   std::cerr.setf(std::ios::showpoint);

   std::cerr << "Reading .asm file\n";
   // store the information on the assembly
   // where do we get the info?
   uint64 genomeID = db->storeGenome("", "igs", "Ixodes scapularis");
   if (genomeID == 0) {
      return false;
   }
   if (db->storeAssembly(AS_UID_fromInteger(getUID(uids)), "Aug 08 2007 12:00AM", AS_UID_fromInteger(genomeID), "bwalenz", "CA", "IXODES", 'P', "Brian Walenz i3 assembly of igs using CA-IXODES tag version")== 0) {
      return false;
   }

   while(ReadProtoMesg_AS(fi, &gen) != EOF)
   {
#ifndef DEBUG_AS_ARD
      if(++mesgCount % 10000 == 0)
         fprintf(stderr, "\r%20lu", mesgCount);
#endif

      if(addGenericMesg2DB(gen) != true)
      {
         return false;
      }
   }
   // commit the last messages
   db->commitSCF2DB();
   db->commitCTP2DB();
   db->commitCTPList2DB();
   db->commitCPS2DB();

   return true;
}
