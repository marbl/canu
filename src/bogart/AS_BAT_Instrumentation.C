
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
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_SetParentAndHang.H"
#include "AS_BAT_Outputs.H"

#include "intervalList.H"

//  Will fail if a read is in unitig 0, or if a read isn't in a unitig.

void
checkUnitigMembership(UnitigVector &unitigs) {
  uint32 *inUnitig = new uint32 [FI->numFragments()+1];
  uint32  noUnitig = 0xffffffff;

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

      if (frg->ident > FI->numFragments())
        fprintf(stderr, "tig %u ufpath[%d] ident %u more than number of reads %u\n",
                tig->id(), fi, frg->ident, FI->numFragments());

      if (inUnitig[frg->ident] != noUnitig)
        fprintf(stderr, "tig %u ufpath[%d] ident %u placed multiple times\n",
                tig->id(), fi, frg->ident);

      assert(frg->ident <= FI->numFragments());   //  Can't be out of range.
      assert(inUnitig[frg->ident] == noUnitig);   //  Read must be not placed yet.

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
      writeLog("ng%03"F_U64P" %9"F_U64P"   lg%03"F_U64P" %8"F_U64P"   sum %11"F_U64P"  (%s)\n",
               nnn, data[i],
               nnn, i+1,
               sum,
               label);

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

  if (logFileFlagSet(LOG_INTERMEDIATE_UNITIGS) == 0)
    return;

  //  Dump to an intermediate store.

  char tigStorePath[FILENAME_MAX];
  sprintf(tigStorePath, "%s.%03u.%s.tigStore", prefix, logFileOrder, name);

  fprintf(stderr, "Creating intermediate tigStore '%s'\n", tigStorePath);

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

  setParentAndHang(unitigs);

  writeUnitigsToStore(unitigs, tigStorePath, tigStorePath, numFragsP, false);
}




#define tCTG  0  //  To a read in a normal tig
#define tRPT  1  //  To a read in a repeat tig
#define tBUB  2  //  To a read in a bubble tig
#define tUNA  3  //  To a read in an 'unassembled' leftover tig
#define tUNU  4  //  To a read not placed in a tig
#define tNOP  5  //  To no read (for best edges)

struct olapsUsed {

  uint64    total;
  //  By definition, satisfied overlaps are in the same tig.

  uint64    doveSatSame[6];
  uint64    contSatSame[6];

  //  Unsatisfied overlaps can be in the same tig...
  uint64    doveUnsatSame[6];
  uint64    contUnsatSame[6];

  //  ...or can be between tigs.

  uint64    doveUnsatDiff[6][6];
  uint64    contUnsatDiff[6][6];
};



uint32
getTigType(Unitig *tg) {
  if (tg == NULL)          return(tUNU);
  if (tg->_isUnassembled)  return(tUNA);
  if (tg->_isBubble)       return(tBUB);
  if (tg->_isRepeat)       return(tRPT);
  if (1)                   return(tCTG);
}


bool
satisfiedOverlap(uint32 rdAlo, uint32 rdAhi, bool rdAfwd, uint32 rdBlo, uint32 rdBhi, bool rdBfwd, bool flipped) {
  return(((rdAhi < rdBlo) || (rdBhi < rdBlo)) ||          //  Not satisfied, no overlap
         ((rdAfwd == rdBfwd) && (flipped == true)) ||     //  Not satisfied, same orient, but flipped overlap
         ((rdAfwd != rdBfwd) && (flipped == false)));     //  Not satisfied, diff orient, but normal overlap
}


//  Iterate over all overlaps (but the only interface we have is by iterating
//  over all reads), and count the number of overlaps satisfied in unitigs.
void
reportOverlaps(UnitigVector &unitigs, const char *prefix, const char *name) {
  olapsUsed   *dd       = new olapsUsed;  //  Dovetail overlaps to non-contained reads
  olapsUsed   *dc       = new olapsUsed;  //  Dovetail overlaps to contained reads
  olapsUsed   *cc       = new olapsUsed;  //  Containment overlaps
  olapsUsed   *bb       = new olapsUsed;  //  Best overlaps

  memset(dd, 0, sizeof(olapsUsed));
  memset(dc, 0, sizeof(olapsUsed));
  memset(cc, 0, sizeof(olapsUsed));
  memset(bb, 0, sizeof(olapsUsed));


  for (uint32 fi=0; fi<FI->numFragments()+1; fi++) {
    if (FI->fragmentLength(fi) == 0)
      continue;

    uint32           rdAid   = fi;
    uint32           tgAid   = Unitig::fragIn(rdAid);
    Unitig          *tgA     = unitigs[tgAid];
    uint32           tgAtype = getTigType(tgA);

    //  Best overlaps exist if the read isn't contained.

    if (OG->isContained(rdAid) == false) {
      BestEdgeOverlap *b5      = OG->getBestEdgeOverlap(fi, false);
      uint32           rd5id   = b5->fragId();
      uint32           tg5id   = Unitig::fragIn(rd5id);
      Unitig          *tg5     = unitigs[tg5id];
      uint32           tg5type = getTigType(tg5);

      BestEdgeOverlap *b3      = OG->getBestEdgeOverlap(fi, true);
      uint32           rd3id   = b3->fragId();
      uint32           tg3id   = Unitig::fragIn(rd3id);
      Unitig          *tg3     = unitigs[tg3id];
      uint32           tg3type = getTigType(tg3);

      bb->total += 2;

      //  If this read isn't even in a tig, add to the unused categories.

      if (tgAid == 0) {
        if (rd5id == 0)
          bb->doveUnsatDiff[tUNU][tNOP]++;
        else
          bb->doveUnsatDiff[tUNU][tg5type]++;

        if (rd3id == 0)
          bb->doveUnsatDiff[tUNU][tNOP]++;
        else
          bb->doveUnsatDiff[tUNU][tg3type]++;
      }

      //  Otherwise, its in a tig, and we need to compare positions.

      else {
        uint32           rdApos  = unitigs[tgAid]->pathPosition(rdAid);
        ufNode          *rdA     = &unitigs[tgAid]->ufpath[rdApos];
        bool             rdAfwd  = (rdA->position.bgn < rdA->position.end);
        int32            rdAlo   = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
        int32            rdAhi   = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

        //  Different tigs?  Unsatisfied.  Same tig?  Grab the reads and check for overlap.

        if (tgA != tg5) {
          bb->doveUnsatDiff[tgAtype][tg5type]++;

        } else if (rd5id == 0) {
          bb->doveUnsatDiff[tgAtype][tNOP]++;

        } else {
          uint32     rd5pos   = unitigs[tg5id]->pathPosition(rd5id);
          ufNode    *rd5      = &unitigs[tg5id]->ufpath[rd5pos];
          bool       rd5fwd   = (rd5->position.bgn < rd5->position.end);
          int32      rd5lo    = (rd5fwd) ? rd5->position.bgn : rd5->position.end;
          int32      rd5hi    = (rd5fwd) ? rd5->position.end : rd5->position.bgn;

          if (satisfiedOverlap(rdAlo, rdAhi, rdAfwd, rd5lo, rd5hi, rd5fwd, (b5->frag3p() == true))) {
            bb->doveSatSame[tgAtype]++;
          } else {
            bb->doveUnsatSame[tgAtype]++;
          }
        }


        if (tgA != tg3) {
          bb->doveUnsatDiff[tgAtype][tg3type]++;

        } else if (rd3id == 0) {
          bb->doveUnsatDiff[tgAtype][tNOP]++;

        } else {
          uint32     rd3pos   = unitigs[tg3id]->pathPosition(rd3id);
          ufNode    *rd3      = &unitigs[tg3id]->ufpath[rd3pos];
          bool       rd3fwd   = (rd3->position.bgn < rd3->position.end);
          int32      rd3lo    = (rd3fwd) ? rd3->position.bgn : rd3->position.end;
          int32      rd3hi    = (rd3fwd) ? rd3->position.end : rd3->position.bgn;

          if (satisfiedOverlap(rdAlo, rdAhi, rdAfwd, rd3lo, rd3hi, rd3fwd, (b3->frag3p() == false))) {
            bb->doveSatSame[tgAtype]++;
          } else {
            bb->doveUnsatSame[tgAtype]++;
          }
        }
      }
    }


    //  For all overlaps.

    uint32        ovlLen = 0;
    BAToverlap   *ovl    = OC->getOverlaps(fi, AS_MAX_ERATE, ovlLen);


    for (uint32 oi=0; oi<ovlLen; oi++) {
      uint32     rdAid     = ovl[oi].a_iid;
      uint32     tgAid     = Unitig::fragIn(rdAid);
      Unitig    *tgA       = unitigs[tgAid];
      uint32     tgAtype   = getTigType(tgA);

      uint32     rdBid     = ovl[oi].b_iid;
      uint32     tgBid     = Unitig::fragIn(rdBid);
      Unitig    *tgB       = unitigs[tgBid];
      uint32     tgBtype   = getTigType(tgB);

      bool       isDove    = ovl[oi].isDovetail();
      bool       contReads = OG->isContained(rdAid) || OG->isContained(rdBid);

      //  Figure out what class of overlap we're counting.

      olapsUsed  *used = NULL;

      if (isDove == false)
        used = cc;
      else
        if (contReads == true)
          used = dc;
        else
          used = dd;

      used->total++;

      //  If to reads not in a tig, unsatisfied.

      if ((tgAid == 0) || (tgBid == 0)) {
        if (isDove)
          used->doveUnsatDiff[tgAtype][tgBtype]++;
        else
          used->contUnsatDiff[tgAtype][tgBtype]++;
        continue;
      }

      //  If in different tigs, unsatisfied.

      if (tgAid != tgBid) {
        if (isDove)
          used->doveUnsatDiff[tgAtype][tgBtype]++;
        else
          used->contUnsatDiff[tgAtype][tgBtype]++;
        continue;
      }

      //  Else, possibly satisfied.  We need to check positions.

      uint32     rdApos   = unitigs[tgAid]->pathPosition(rdAid);
      ufNode    *rdA      = &unitigs[tgAid]->ufpath[rdApos];
      bool       rdAfwd   = (rdA->position.bgn < rdA->position.end);
      int32      rdAlo    = (rdAfwd) ? rdA->position.bgn : rdA->position.end;
      int32      rdAhi    = (rdAfwd) ? rdA->position.end : rdA->position.bgn;

      uint32     rdBpos   = unitigs[tgBid]->pathPosition(rdBid);
      ufNode    *rdB      = &unitigs[tgBid]->ufpath[rdBpos];
      bool       rdBfwd   = (rdB->position.bgn < rdB->position.end);
      int32      rdBlo    = (rdBfwd) ? rdB->position.bgn : rdB->position.end;
      int32      rdBhi    = (rdBfwd) ? rdB->position.end : rdB->position.bgn;

      //  If overlapping and correctly oriented, good enough for now.  Do we want to care about
      //  overlap length?  Nah, there's enough fudging (still, I think) in placement that it'd be
      //  tough to get that usefully precise.

      if (satisfiedOverlap(rdAlo, rdAhi, rdAfwd, rdBlo, rdBhi, rdBfwd, ovl[oi].flipped)) {
        if (isDove)
          used->doveUnsatSame[tgAtype]++;
        else
          used->contUnsatSame[tgAtype]++;

      } else {
        if (isDove)
          used->doveSatSame[tgAtype]++;
        else
          used->contSatSame[tgAtype]++;
      }
    }
  }

  //  Merge the symmetrical counts

  for (uint32 ii=0; ii<6; ii++) {
    for (uint32 jj=ii+1; jj<6; jj++) {
      bb->doveUnsatDiff[ii][jj] += bb->doveUnsatDiff[jj][ii];      bb->doveUnsatDiff[jj][ii] = UINT64_MAX;
      dd->doveUnsatDiff[ii][jj] += dd->doveUnsatDiff[jj][ii];      dd->doveUnsatDiff[jj][ii] = UINT64_MAX;
      dc->doveUnsatDiff[ii][jj] += dc->doveUnsatDiff[jj][ii];      dc->doveUnsatDiff[jj][ii] = UINT64_MAX;
      cc->doveUnsatDiff[ii][jj] += cc->doveUnsatDiff[jj][ii];      cc->doveUnsatDiff[jj][ii] = UINT64_MAX;

      bb->contUnsatDiff[ii][jj] += bb->contUnsatDiff[jj][ii];      bb->contUnsatDiff[jj][ii] = UINT64_MAX;
      dd->contUnsatDiff[ii][jj] += dd->contUnsatDiff[jj][ii];      dd->contUnsatDiff[jj][ii] = UINT64_MAX;
      dc->contUnsatDiff[ii][jj] += dc->contUnsatDiff[jj][ii];      dc->contUnsatDiff[jj][ii] = UINT64_MAX;
      cc->contUnsatDiff[ii][jj] += cc->contUnsatDiff[jj][ii];      cc->contUnsatDiff[jj][ii] = UINT64_MAX;
    }
  }


  //  Emit a nicely formatted report.

#define B(X)  (100.0 * (X) / (bb->total))
#define P(X)  (100.0 * (X) / (dd->total))
#define Q(X)  (100.0 * (X) / (dc->total))
#define R(X)  (100.0 * (X) / (cc->total))

  writeLog("--------------------------------------------------------------------------------\n");
  writeLog("OVERLAP FATE\n");
  writeLog("\n");
  writeLog("dovetail overlaps (best)              "F_U64"\n", bb->total);
  writeLog("dovetail overlaps                     "F_U64"\n", dd->total);
  writeLog("dovetail overlaps to contained reads  "F_U64"\n", dc->total);
  writeLog("containment overlaps                  "F_U64"\n", cc->total);
  writeLog("\n");
  writeLog("SATISFIED best edges         DOVETAIL\n");
  writeLog("---------                ------------ -------\n");
  writeLog("contig                   %12"F_U64P" %6.2f%%\n", bb->doveSatSame[tCTG], B(bb->doveSatSame[tCTG]));
  writeLog("repeat contig            %12"F_U64P" %6.2f%%\n", bb->doveSatSame[tRPT], B(bb->doveSatSame[tRPT]));
  writeLog("bubble                   %12"F_U64P" %6.2f%%\n", bb->doveSatSame[tBUB], B(bb->doveSatSame[tBUB]));
  writeLog("\n");
  writeLog("UNSATISFIED best edges       DOVETAIL\n");
  writeLog("-----------              ------------ -------\n");
  writeLog("contig                   %12"F_U64P" %6.2f%%\n", bb->doveUnsatSame[tCTG], B(bb->doveUnsatSame[tCTG]));
  writeLog("repeat                   %12"F_U64P" %6.2f%%\n", bb->doveUnsatSame[tRPT], B(bb->doveUnsatSame[tRPT]));
  writeLog("bubble                   %12"F_U64P" %6.2f%%\n", bb->doveUnsatSame[tBUB], B(bb->doveUnsatSame[tBUB]));
  writeLog("unassembled              %12"F_U64P" %6.2f%%\n", bb->doveUnsatSame[tUNA], B(bb->doveUnsatSame[tUNA]));
  writeLog("unused                   %12"F_U64P" %6.2f%%\n", bb->doveUnsatSame[tUNU], B(bb->doveUnsatSame[tUNU]));
  writeLog("\n");
  writeLog("UNSATISFIED best edges       DOVETAIL\n");
  writeLog("-----------              ------------ -------\n");
  writeLog("contig-contig            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tCTG], B(bb->doveUnsatDiff[tCTG][tCTG]));
  writeLog("contig-repeat            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tRPT], B(bb->doveUnsatDiff[tCTG][tRPT]));
  writeLog("contig-bubble            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tBUB], B(bb->doveUnsatDiff[tCTG][tBUB]));
  writeLog("contig-unassembled       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tUNA], B(bb->doveUnsatDiff[tCTG][tUNA]));
  writeLog("contig-unused            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tUNU], B(bb->doveUnsatDiff[tCTG][tUNU]));
  writeLog("contig-none              %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tCTG][tNOP], B(bb->doveUnsatDiff[tCTG][tNOP]));
  writeLog("\n");
//writeLog("repeat-contig            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tCTG], B(bb->doveUnsatDiff[tRPT][tCTG]));
  writeLog("repeat-repeat            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tRPT], B(bb->doveUnsatDiff[tRPT][tRPT]));
  writeLog("repeat-bubble            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tBUB], B(bb->doveUnsatDiff[tRPT][tBUB]));
  writeLog("repeat-unassembled       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tUNA], B(bb->doveUnsatDiff[tRPT][tUNA]));
  writeLog("repeat-unused            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tUNU], B(bb->doveUnsatDiff[tRPT][tUNU]));
  writeLog("repeat-none              %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tRPT][tNOP], B(bb->doveUnsatDiff[tRPT][tNOP]));
  writeLog("\n");
//writeLog("bubble-contig            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tCTG], B(bb->doveUnsatDiff[tBUB][tCTG]));
//writeLog("bubble-repeat            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tRPT], B(bb->doveUnsatDiff[tBUB][tRPT]));
  writeLog("bubble-bubble            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tBUB], B(bb->doveUnsatDiff[tBUB][tBUB]));
  writeLog("bubble-unassembled       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tUNA], B(bb->doveUnsatDiff[tBUB][tUNA]));
  writeLog("bubble-unused            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tUNU], B(bb->doveUnsatDiff[tBUB][tUNU]));
  writeLog("bubble-none              %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tBUB][tNOP], B(bb->doveUnsatDiff[tBUB][tNOP]));
  writeLog("\n");
//writeLog("unassembled-contig       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tCTG], B(bb->doveUnsatDiff[tUNA][tCTG]));
//writeLog("unassembled-repeat       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tRPT], B(bb->doveUnsatDiff[tUNA][tRPT]));
//writeLog("unassembled-bubble       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tBUB], B(bb->doveUnsatDiff[tUNA][tBUB]));
  writeLog("unassembled-unassembled  %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tUNA], B(bb->doveUnsatDiff[tUNA][tUNA]));
  writeLog("unassembled-unused       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tUNU], B(bb->doveUnsatDiff[tUNA][tUNU]));
  writeLog("unassembled-none         %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNA][tNOP], B(bb->doveUnsatDiff[tUNA][tNOP]));
  writeLog("\n");
//writeLog("unused-contig            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tCTG], B(bb->doveUnsatDiff[tUNU][tCTG]))
//writeLog("unused-repeat            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tRPT], B(bb->doveUnsatDiff[tUNU][tRPT]));
//writeLog("unused-bubble            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tBUB], B(bb->doveUnsatDiff[tUNU][tBUB]));
//writeLog("unused-unassembled       %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tUNA], B(bb->doveUnsatDiff[tUNU][tUNA]));
  writeLog("unused-unused            %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tUNU], B(bb->doveUnsatDiff[tUNU][tUNU]));
  writeLog("unused-none              %12"F_U64P" %6.2f%%\n", bb->doveUnsatDiff[tUNU][tNOP], B(bb->doveUnsatDiff[tUNU][tNOP]));
  writeLog("\n");
  writeLog("\n");
  writeLog("\n");
  writeLog("SATISFIED all overlaps       DOVETAIL              DOVECONT           CONTAINMENT\n");
  writeLog("---------                ------------ -------  ------------ -------  ------------ -------\n");
  writeLog("contig                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveSatSame[tCTG], P(dd->doveSatSame[tCTG]), dc->doveSatSame[tCTG], Q(dc->doveSatSame[tCTG]), cc->contSatSame[tCTG], R(cc->contSatSame[tCTG]));
  writeLog("repeat contig            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveSatSame[tRPT], P(dd->doveSatSame[tRPT]), dc->doveSatSame[tRPT], Q(dc->doveSatSame[tRPT]), cc->contSatSame[tRPT], R(cc->contSatSame[tRPT]));
  writeLog("bubble                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveSatSame[tBUB], P(dd->doveSatSame[tBUB]), dc->doveSatSame[tBUB], Q(dc->doveSatSame[tBUB]), cc->contSatSame[tBUB], R(cc->contSatSame[tBUB]));
  writeLog("\n");
  writeLog("UNSATISFIED all overlaps     DOVETAIL              DOVECONT           CONTAINMENT\n");
  writeLog("-----------              ------------ -------  ------------ -------  ------------ -------\n");
  writeLog("contig                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatSame[tCTG], P(dd->doveUnsatSame[tCTG]), dc->doveUnsatSame[tCTG], Q(dc->doveUnsatSame[tCTG]), cc->contUnsatSame[tCTG], R(cc->contUnsatSame[tCTG]));
  writeLog("repeat                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatSame[tRPT], P(dd->doveUnsatSame[tRPT]), dc->doveUnsatSame[tRPT], Q(dc->doveUnsatSame[tRPT]), cc->contUnsatSame[tRPT], R(cc->contUnsatSame[tRPT]));
  writeLog("bubble                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatSame[tBUB], P(dd->doveUnsatSame[tBUB]), dc->doveUnsatSame[tBUB], Q(dc->doveUnsatSame[tBUB]), cc->contUnsatSame[tBUB], R(cc->contUnsatSame[tBUB]));
  writeLog("unassembled              %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatSame[tUNA], P(dd->doveUnsatSame[tUNA]), dc->doveUnsatSame[tUNA], Q(dc->doveUnsatSame[tUNA]), cc->contUnsatSame[tUNA], R(cc->contUnsatSame[tUNA]));
  writeLog("unused                   %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatSame[tUNU], P(dd->doveUnsatSame[tUNU]), dc->doveUnsatSame[tUNU], Q(dc->doveUnsatSame[tUNU]), cc->contUnsatSame[tUNU], R(cc->contUnsatSame[tUNU]));
  writeLog("\n");
  writeLog("UNSATISFIED all overlaps     DOVETAIL              DOVECONT           CONTAINMENT\n");
  writeLog("-----------              ------------ -------  ------------ -------  ------------ -------\n");
  writeLog("contig-contig            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tCTG][tCTG], P(dd->doveUnsatDiff[tCTG][tCTG]), dc->doveUnsatDiff[tCTG][tCTG], Q(dc->doveUnsatDiff[tCTG][tCTG]), cc->contUnsatDiff[tCTG][tCTG], R(cc->contUnsatDiff[tCTG][tCTG]));
  writeLog("contig-repeat            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tCTG][tRPT], P(dd->doveUnsatDiff[tCTG][tRPT]), dc->doveUnsatDiff[tCTG][tRPT], Q(dc->doveUnsatDiff[tCTG][tRPT]), cc->contUnsatDiff[tCTG][tRPT], R(cc->contUnsatDiff[tCTG][tRPT]));
  writeLog("contig-bubble            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tCTG][tBUB], P(dd->doveUnsatDiff[tCTG][tBUB]), dc->doveUnsatDiff[tCTG][tBUB], Q(dc->doveUnsatDiff[tCTG][tBUB]), cc->contUnsatDiff[tCTG][tBUB], R(cc->contUnsatDiff[tCTG][tBUB]));
  writeLog("contig-unassembled       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tCTG][tUNA], P(dd->doveUnsatDiff[tCTG][tUNA]), dc->doveUnsatDiff[tCTG][tUNA], Q(dc->doveUnsatDiff[tCTG][tUNA]), cc->contUnsatDiff[tCTG][tUNA], R(cc->contUnsatDiff[tCTG][tUNA]));
  writeLog("contig-unused            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tCTG][tUNU], P(dd->doveUnsatDiff[tCTG][tUNU]), dc->doveUnsatDiff[tCTG][tUNU], Q(dc->doveUnsatDiff[tCTG][tUNU]), cc->contUnsatDiff[tCTG][tUNU], R(cc->contUnsatDiff[tCTG][tUNU]));
  writeLog("\n");
//writeLog("repeat-contig            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tRPT][tCTG], P(dd->doveUnsatDiff[tRPT][tCTG]), dc->doveUnsatDiff[tRPT][tCTG], Q(dc->doveUnsatDiff[tRPT][tCTG]), cc->contUnsatDiff[tRPT][tCTG], R(cc->contUnsatDiff[tRPT][tCTG]));
  writeLog("repeat-repeat            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tRPT][tRPT], P(dd->doveUnsatDiff[tRPT][tRPT]), dc->doveUnsatDiff[tRPT][tRPT], Q(dc->doveUnsatDiff[tRPT][tRPT]), cc->contUnsatDiff[tRPT][tRPT], R(cc->contUnsatDiff[tRPT][tRPT]));
  writeLog("repeat-bubble            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tRPT][tBUB], P(dd->doveUnsatDiff[tRPT][tBUB]), dc->doveUnsatDiff[tRPT][tBUB], Q(dc->doveUnsatDiff[tRPT][tBUB]), cc->contUnsatDiff[tRPT][tBUB], R(cc->contUnsatDiff[tRPT][tBUB]));
  writeLog("repeat-unassembled       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tRPT][tUNA], P(dd->doveUnsatDiff[tRPT][tUNA]), dc->doveUnsatDiff[tRPT][tUNA], Q(dc->doveUnsatDiff[tRPT][tUNA]), cc->contUnsatDiff[tRPT][tUNA], R(cc->contUnsatDiff[tRPT][tUNA]));
  writeLog("repeat-unused            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tRPT][tUNU], P(dd->doveUnsatDiff[tRPT][tUNU]), dc->doveUnsatDiff[tRPT][tUNU], Q(dc->doveUnsatDiff[tRPT][tUNU]), cc->contUnsatDiff[tRPT][tUNU], R(cc->contUnsatDiff[tRPT][tUNU]));
  writeLog("\n");
//writeLog("bubble-contig            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tBUB][tCTG], P(dd->doveUnsatDiff[tBUB][tCTG]), dc->doveUnsatDiff[tBUB][tCTG], Q(dc->doveUnsatDiff[tBUB][tCTG]), cc->contUnsatDiff[tBUB][tCTG], R(cc->contUnsatDiff[tBUB][tCTG]));
//writeLog("bubble-repeat            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tBUB][tRPT], P(dd->doveUnsatDiff[tBUB][tRPT]), dc->doveUnsatDiff[tBUB][tRPT], Q(dc->doveUnsatDiff[tBUB][tRPT]), cc->contUnsatDiff[tBUB][tRPT], R(cc->contUnsatDiff[tBUB][tRPT]));
  writeLog("bubble-bubble            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tBUB][tBUB], P(dd->doveUnsatDiff[tBUB][tBUB]), dc->doveUnsatDiff[tBUB][tBUB], Q(dc->doveUnsatDiff[tBUB][tBUB]), cc->contUnsatDiff[tBUB][tBUB], R(cc->contUnsatDiff[tBUB][tBUB]));
  writeLog("bubble-unassembled       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tBUB][tUNA], P(dd->doveUnsatDiff[tBUB][tUNA]), dc->doveUnsatDiff[tBUB][tUNA], Q(dc->doveUnsatDiff[tBUB][tUNA]), cc->contUnsatDiff[tBUB][tUNA], R(cc->contUnsatDiff[tBUB][tUNA]));
  writeLog("bubble-unused            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tBUB][tUNU], P(dd->doveUnsatDiff[tBUB][tUNU]), dc->doveUnsatDiff[tBUB][tUNU], Q(dc->doveUnsatDiff[tBUB][tUNU]), cc->contUnsatDiff[tBUB][tUNU], R(cc->contUnsatDiff[tBUB][tUNU]));
  writeLog("\n");
//writeLog("unassembled-contig       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNA][tCTG], P(dd->doveUnsatDiff[tUNA][tCTG]), dc->doveUnsatDiff[tUNA][tCTG], Q(dc->doveUnsatDiff[tUNA][tCTG]), cc->contUnsatDiff[tUNA][tCTG], R(cc->contUnsatDiff[tUNA][tCTG]));
//writeLog("unassembled-repeat       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNA][tRPT], P(dd->doveUnsatDiff[tUNA][tRPT]), dc->doveUnsatDiff[tUNA][tRPT], Q(dc->doveUnsatDiff[tUNA][tRPT]), cc->contUnsatDiff[tUNA][tRPT], R(cc->contUnsatDiff[tUNA][tRPT]));
//writeLog("unassembled-bubble       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNA][tBUB], P(dd->doveUnsatDiff[tUNA][tBUB]), dc->doveUnsatDiff[tUNA][tBUB], Q(dc->doveUnsatDiff[tUNA][tBUB]), cc->contUnsatDiff[tUNA][tBUB], R(cc->contUnsatDiff[tUNA][tBUB]));
  writeLog("unassembled-unassembled  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNA][tUNA], P(dd->doveUnsatDiff[tUNA][tUNA]), dc->doveUnsatDiff[tUNA][tUNA], Q(dc->doveUnsatDiff[tUNA][tUNA]), cc->contUnsatDiff[tUNA][tUNA], R(cc->contUnsatDiff[tUNA][tUNA]));
  writeLog("unassembled-unused       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNA][tUNU], P(dd->doveUnsatDiff[tUNA][tUNU]), dc->doveUnsatDiff[tUNA][tUNU], Q(dc->doveUnsatDiff[tUNA][tUNU]), cc->contUnsatDiff[tUNA][tUNU], R(cc->contUnsatDiff[tUNA][tUNU]));
  writeLog("\n");
//writeLog("unused-contig            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNU][tCTG], P(dd->doveUnsatDiff[tUNU][tCTG]), dc->doveUnsatDiff[tUNU][tCTG], Q(dc->doveUnsatDiff[tUNU][tCTG]), cc->contUnsatDiff[tUNU][tCTG], R(cc->contUnsatDiff[tUNU][tCTG]));
//writeLog("unused-repeat            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNU][tRPT], P(dd->doveUnsatDiff[tUNU][tRPT]), dc->doveUnsatDiff[tUNU][tRPT], Q(dc->doveUnsatDiff[tUNU][tRPT]), cc->contUnsatDiff[tUNU][tRPT], R(cc->contUnsatDiff[tUNU][tRPT]));
//writeLog("unused-bubble            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNU][tBUB], P(dd->doveUnsatDiff[tUNU][tBUB]), dc->doveUnsatDiff[tUNU][tBUB], Q(dc->doveUnsatDiff[tUNU][tBUB]), cc->contUnsatDiff[tUNU][tBUB], R(cc->contUnsatDiff[tUNU][tBUB]));
//writeLog("unused-unassembled       %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNU][tUNA], P(dd->doveUnsatDiff[tUNU][tUNA]), dc->doveUnsatDiff[tUNU][tUNA], Q(dc->doveUnsatDiff[tUNU][tUNA]), cc->contUnsatDiff[tUNU][tUNA], R(cc->contUnsatDiff[tUNU][tUNA]));
  writeLog("unused-unused            %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%  %12"F_U64P" %6.2f%%\n", dd->doveUnsatDiff[tUNU][tUNU], P(dd->doveUnsatDiff[tUNU][tUNU]), dc->doveUnsatDiff[tUNU][tUNU], Q(dc->doveUnsatDiff[tUNU][tUNU]), cc->contUnsatDiff[tUNU][tUNU], R(cc->contUnsatDiff[tUNU][tUNU]));
  writeLog("\n");
  writeLog("\n");

  delete dd;
  delete dc;
  delete cc;
}

