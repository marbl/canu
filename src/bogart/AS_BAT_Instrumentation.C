
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Instrumentation.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-27
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2014-DEC-23
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Unitig.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"

#include "intervalList.H"

//  Will fail if a read is in unitig 0, or if a read isn't in a unitig.
//  Used to also compute some simple size statistics.

void
checkUnitigMembership(UnitigVector &unitigs) {
  uint32 *inUnitig = new uint32 [FI->numFragments()+1];

  //  All reads start of not placed in a unitig.

  for (uint32 i=0; i<FI->numFragments()+1; i++)
    inUnitig[i] = noUnitig;

  //  Over all unitigs, remember where each read is.

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *tig = unitigs[ti];
    int32    len = 0;

    if (tig == NULL)
      continue;

    for (uint32 fi=0; fi<tig->ufpath.size(); fi++) {
      ufNode  *frg = &tig->ufpath[fi];

      assert(frg->ident <= FI->numFragments());

      inUnitig[frg->ident] = ti;
    }
  }

  //  Find any read not placed in a unitig.

  for (uint32 i=0; i<FI->numFragments()+1; i++) {
    if (FI->fragmentLength(i) == 0)  //  Deleted read.
      continue;

    assert(inUnitig[i] != 0);         //  There shouldn't be a unitig 0.
    assert(inUnitig[i] != noUnitig);  //  The read should be in a unitig.
  }

  delete [] inUnitig;
}


//  For every unitig, report the best overlaps contained in the
//  unitig, and all overlaps contained in the unitig.
void
reportOverlapsUsed(UnitigVector &unitigs, const char *prefix, const char *name) {

  if (logFileFlagSet(LOG_OVERLAPS_USED) == 0)
    return;

  char ovlPath[FILENAME_MAX];
  sprintf(ovlPath, "%s.%03u.%s.overlaps", prefix, logFileOrder, name);

  FILE *F = fopen(ovlPath, "w");

  if (F == NULL)
    return;

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    for (uint32 fi=0; fi<utg->ufpath.size(); fi++) {
      ufNode  *frg = &utg->ufpath[fi];

      //  Where is our best overlap?  Contained or dovetail?

      BestEdgeOverlap *bestedge5 = OG->getBestEdgeOverlap(frg->ident, false);
      BestEdgeOverlap *bestedge3 = OG->getBestEdgeOverlap(frg->ident, true);

      uint32           bestident5 = 0;
      uint32           bestident3 = 0;

      if (bestedge5)
        bestident5 = bestedge5->fragId();

      if (bestedge3)
        bestident3 = bestedge3->fragId();

      //  Now search ahead, reporting any overlap to any fragment.
      //
      for (uint32 oi=fi+1; oi<utg->ufpath.size(); oi++) {
        ufNode  *ooo = &utg->ufpath[oi];

        int frgbgn = MIN(frg->position.bgn, frg->position.end);
        int frgend = MAX(frg->position.bgn, frg->position.end);

        int ooobgn = MIN(ooo->position.bgn, ooo->position.end);
        int oooend = MAX(ooo->position.bgn, ooo->position.end);

        if ((frgbgn <= ooobgn) && (ooobgn + 40 < frgend)) {
          BestContainment *bestcont  = OG->getBestContainer(ooo->ident);

          uint32           bestident  = 0;
          if (bestcont->isContained)
            bestident = bestcont->container;

          bool isBest = ((frg->ident == bestident) ||
                         (ooo->ident == bestident5) ||
                         (ooo->ident == bestident3));

          fprintf(F, "%d\t%d%s\n", frg->ident, ooo->ident, (isBest) ? ((bestident) ? "\tbc" : "\tbe") : "");
        }

        if (frgend < ooobgn)
          break;
      }
    }
  }

  fclose(F);
}


//  Decides if a unitig is unassembled.  The other classifications (isBubble, isCircular, isRepeat)
//  are made when the type is processed (e.g., when bubbles are popped).
//
//  A unitig is unassembled if:
//    1) it has fewer than R reads (R=2)
//    2) it is shorter than S bases (S=1000)
//    3) a single read spans at least fraction F of the lenth (F=1.0)
//    4) at least fraction F of the unitig is below read depth D (F=1.0, D=2)
//
void
classifyUnitigsAsUnassembled(UnitigVector &unitigs,
                             uint32        fewReadsNumber,
                             uint32        tooShortLength,
                             double        spanFraction,
                             double        lowcovFraction,   uint32  lowcovDepth) {
  uint32  nTooFew   = 0;
  uint32  nShort    = 0;
  uint32  nSingle   = 0;
  uint32  nCoverage = 0;
  uint32  nContig   = 0;

  uint64  bTooFew   = 0;
  uint64  bShort    = 0;
  uint64  bSingle   = 0;
  uint64  bCoverage = 0;
  uint64  bContig   = 0;

  writeLog("==> FILTERING UNASSEMBLED CRUD\n");

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    utg->_isUnassembled = false;

    //  Rule 1.  Too few reads.

    if (utg->ufpath.size() < fewReadsNumber) {
      writeLog("unitig %u unassembled - too few reads (%u < %u)\n", ti, utg->ufpath.size(), fewReadsNumber);
      utg->_isUnassembled = true;
      nTooFew += 1;
      bTooFew += utg->getLength();
      continue;
    }

    //  Rule 2.  Short.

    if (utg->getLength() < tooShortLength) {
      writeLog("unitig %u unassembled - too short (%u < %u)\n", ti, utg->getLength(), tooShortLength);
      utg->_isUnassembled = true;
      nShort += 1;
      bShort += utg->getLength();
      continue;
    }

    //  Rule 3.  Single read spans large fraction of tig.

    for (uint32 oi=0; oi<utg->ufpath.size(); oi++) {
      ufNode  *frg = &utg->ufpath[oi];

      int frgbgn = MIN(frg->position.bgn, frg->position.end);
      int frgend = MAX(frg->position.bgn, frg->position.end);

      if (frgend - frgbgn > utg->getLength() * spanFraction) {
        writeLog("unitig %u unassembled - single read spans unitig (read %u %u-%u spans fraction %f > %f\n",
                 ti, frg->ident, frg->position.bgn, frg->position.end, (double)(frgend - frgbgn) / utg->getLength(), spanFraction);
        utg->_isUnassembled = true;
        nSingle += 1;
        bSingle += utg->getLength();
        break;
      }
    }
    if (utg->_isUnassembled)
      continue;
        
    //  Rule 4.  Low coverage.

    intervalList<int32>  IL;

    for (uint32 oi=0; oi<utg->ufpath.size(); oi++) {
      ufNode  *frg = &utg->ufpath[oi];

      int frgbgn = MIN(frg->position.bgn, frg->position.end);
      int frgend = MAX(frg->position.bgn, frg->position.end);

      IL.add(frgbgn, frgend - frgbgn);
    }
    
    intervalList<int32>  ID(IL);

    uint32  basesLow  = 0;
    uint32  basesHigh = 0;
    
    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
      if (ID.depth(ii) < lowcovDepth)
        basesLow  += ID.hi(ii) - ID.lo(ii) + 1;
      else
        basesHigh += ID.hi(ii) - ID.lo(ii) + 1;

    double  lowcov = (double)basesLow / (basesLow + basesHigh);

    if (lowcov >= lowcovFraction) {
      writeLog("Unitig %u unassembled - low coverage (%.4f > %.4f at < %ux coverage)\n",
               ti, lowcov, lowcovFraction, lowcovDepth);
      utg->_isUnassembled = true;
      nCoverage += 1;
      bCoverage += utg->getLength();
      continue;
    }

    //  Otherwise, unitig is assembled!

    nContig += 1;
    bContig += utg->getLength();
  }

  writeLog("unassembled filter:  %6u unitigs %11lu bases -- too few reads\n", nTooFew, bTooFew);
  writeLog("unassembled filter:  %6u unitigs %11lu bases -- too short\n", nShort, bShort);
  writeLog("unassembled filter:  %6u unitigs %11lu bases -- single spanning read\n", nSingle, bSingle);
  writeLog("unassembled filter:  %6u unitigs %11lu bases -- low coverage\n", nCoverage, bCoverage);
  writeLog("unassembled filter:  %6u unitigs %11lu bases -- acceptable contigs\n", nContig, bContig);
}


void
reportN50(vector<uint32> &data, char const *label, uint64 genomeSize) {
  uint64  cnt = data.size();
  uint64  sum = 0;
  uint64  tot = 0;
  uint64  nnn = 10;
  uint64  siz = 0;

  if (cnt == 0)
    return;

  //  Duplicates tgTigSizeAnalysis::printSummary()

  sort(data.begin(), data.end(), greater<uint32>());

  for (uint64 i=0; i<cnt; i++)
    tot += data[i];

  writeLog("%s (%u tigs) (%u length) (%u average) (%.2fx coverage)\n",
           label, cnt, tot, tot / cnt, (double)tot / genomeSize);

  if (genomeSize > 0)
    siz = genomeSize;
  else
    siz = tot;

  for (uint64 i=0; i<cnt; i++) {
    sum += data[i];

    while (siz * nnn / 100 < sum) {
      writeLog("ng%03"F_U64P"   lg%03d %6"F_U64P"   sum %10"F_U64P"\n",
               nnn, data[i],
               nnn, i+1,
               sum);

      nnn += 10;
    }
  }

}


void
reportUnitigs(UnitigVector &unitigs, const char *prefix, const char *name, uint64 genomeSize) {

  //  Generate n50.  Assumes unitigs have been 'classified' already.

  vector<uint32>   unassembledLength;
  vector<uint32>   bubbleLength;
  vector<uint32>   repeatLength;
  vector<uint32>   circularLength;
  vector<uint32>   contigLength;

  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    if (utg->_isUnassembled) {
      unassembledLength.push_back(utg->getLength());
    }

    else if (utg->_isBubble) {
      bubbleLength.push_back(utg->getLength());
    }

    else if (utg->_isRepeat) {
      repeatLength.push_back(utg->getLength());
    }

    else if (utg->_isCircular) {
      circularLength.push_back(utg->getLength());
    }

    else {
      contigLength.push_back(utg->getLength());
    }
  }

  reportN50(unassembledLength, "UNASSEMBLED", genomeSize);
  reportN50(bubbleLength,      "BUBBLE",      genomeSize);
  reportN50(repeatLength,      "REPEAT",      genomeSize);
  reportN50(circularLength,    "CIRCULAR",    genomeSize);
  reportN50(contigLength,      "CONTIGS",     genomeSize);

  //  Dump to an intermediate store.

  if (logFileFlagSet(LOG_INTERMEDIATE_UNITIGS) == 0)
    return;

  uint32  numFragsT  = 0;
  uint32  numFragsP  = 0;
  uint64  utgLen     = 0;

  //  Compute average frags per partition.
  for (uint32  ti=0; ti<unitigs.size(); ti++) {
    Unitig  *utg = unitigs[ti];

    if (utg == NULL)
      continue;

    numFragsT += utg->ufpath.size();

    if (utg->ufpath.size() > 2)
      utgLen    += utg->getLength();
  }

  if      (utgLen < 16 * 1024 * 1024)
    numFragsP = numFragsT / 7;
  else if (utgLen < 64 * 1024 * 1024)
    numFragsP = numFragsT / 63;
  else
    numFragsP = numFragsT / 127;

  //  Dump the unitigs to an intermediate store.

  char tigStorePath[FILENAME_MAX];
  sprintf(tigStorePath, "%s.%03u.%s.tigStore", prefix, logFileOrder, name);

  setParentAndHang(unitigs);

  writeUnitigsToStore(unitigs, tigStorePath, tigStorePath, numFragsP, false);
}

