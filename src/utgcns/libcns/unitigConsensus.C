
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
 *    src/AS_CNS/MultiAlignUnitig.C
 *    src/AS_CNS/MultiAlignUnitig.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MultiAlignUnitig.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-OCT-04
 *      are Copyright 2005-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-JUN-05
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2011-OCT-27
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-AUG-11
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "unitigConsensus.H"

// for pbdagcon
#include "Alignment.H"
#include "AlnGraphBoost.H"
#include "edlib.H"

#include "NDalign.H"

#include <set>

using namespace std;



//  Define this.  Use the faster aligner from overlapper.  If not defined,
//  a full O(n^2) DP is computed.
//
#undef  WITH_NDALIGN
#define WITH_NDALIGN



unitigConsensus::unitigConsensus(gkStore  *gkpStore_,
                                 double    errorRate_,
                                 double    errorRateMax_,
                                 uint32    minOverlap_) {

  gkpStore        = gkpStore_;

  tig             = NULL;
  numfrags        = 0;
  trace           = NULL;
  abacus          = NULL;
  utgpos          = NULL;
  cnspos          = NULL;
  tiid            = 0;
  piid            = -1;

  minOverlap      = minOverlap_;
  errorRate       = errorRate_;
  errorRateMax    = errorRateMax_;

  oaPartial       = NULL;
  oaFull          = NULL;
}


unitigConsensus::~unitigConsensus() {
  delete [] trace;
  delete    abacus;

  delete [] utgpos;
  delete [] cnspos;

  delete    oaPartial;
  delete    oaFull;
}




void
unitigConsensus::reportStartingWork(void) {
  if (showProgress())
    fprintf(stderr, "unitigConsensus()-- processing read %u/%u id %d pos %d,%d anchor %d,%d,%d -- length %u\n",
            tiid+1, numfrags,
            utgpos[tiid].ident(),
            utgpos[tiid].min(),
            utgpos[tiid].max(),
            utgpos[tiid].anchor(),
            utgpos[tiid].aHang(),
            utgpos[tiid].bHang(),
            abacus->numberOfColumns());

  if (showPlacementBefore())
    for (int32 x=0; x<=tiid; x++)
      fprintf(stderr, "unitigConsensus()-- mid %10d  utgpos %7d,%7d  cnspos %7d,%7d  anchor %10d,%6d,%6d\n",
              utgpos[x].ident(),
              utgpos[x].min(), utgpos[x].max(),
              cnspos[x].min(), cnspos[x].max(),
              utgpos[x].anchor(), utgpos[x].aHang(), utgpos[x].bHang());
}


void
unitigConsensus::reportFailure(void) {
  fprintf(stderr, "unitigConsensus()-- failed to align fragment %d in unitig %d.\n",
          utgpos[tiid].ident(), tig->tigID());
}


void
unitigConsensus::reportSuccess() {
  //fprintf(stderr, "unitigConsensus()-- fragment %d aligned in unitig %d.\n",
  //        utgpos[tiid].ident(), tig->tigID());
}



//  Dump the unitig and reads to a single file.  We should also probably save any parameters,
//  but then it's not clear what to do with them.
//
bool
unitigConsensus::savePackage(FILE   *outPackageFile,
                             tgTig  *tig) {

  //  Saving the tig is easy, just use the standard dump.

  tig->saveToStream(outPackageFile);

  //  Saving the reads is also easy, but it's a non-standard dump.

  for (uint32 ii=0; ii<tig->numberOfChildren(); ii++)
    gkpStore->gkStore_saveReadToStream(outPackageFile, tig->getChild(ii)->ident());

  return(true);
}



bool
unitigConsensus::generate(tgTig                     *tig_,
                          map<uint32, gkRead *>     *inPackageRead_,
                          map<uint32, gkReadData *> *inPackageReadData_) {

  tig      = tig_;
  numfrags = tig->numberOfChildren();

  if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
    fprintf(stderr, "generate()--  Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    goto returnFailure;
  }

  while (moreFragments()) {
    reportStartingWork();

    //  First attempt, all default parameters

    if (computePositionFromAnchor()    && alignFragment())  goto applyAlignment;
    if (computePositionFromLayout()    && alignFragment())  goto applyAlignment;
    if (computePositionFromAlignment() && alignFragment())  goto applyAlignment;

    //  Second attempt, default parameters after recomputing consensus sequence.

    if (showAlgorithm())
      fprintf(stderr, "generate()-- recompute full consensus\n");

    recomputeConsensus(showMultiAlignments());

    if (computePositionFromAnchor()    && alignFragment())  goto applyAlignment;
    if (computePositionFromLayout()    && alignFragment())  goto applyAlignment;
    if (computePositionFromAlignment() && alignFragment())  goto applyAlignment;

    //  Third attempot, use whatever aligns.  (alignFragment(true) forced it to align, but that's breaking the consensus with garbage alignments)

    if (computePositionFromAlignment() && alignFragment(true))  goto applyAlignment;

    //  Nope, failed to align.

    reportFailure();
    continue;

  applyAlignment:
    setErrorRate(errorRate);
    setMinOverlap(minOverlap);

    reportSuccess();

    abacus->applyAlignment(tiid, traceABgn, traceBBgn, trace, traceLen);

    refreshPositions();
  }

  generateConsensus(tig);

  return(true);

 returnFailure:
  fprintf(stderr, "generate()-- unitig %d FAILED.\n", tig->tigID());

  //  tgTig should have no changes.

  return(false);
}



void
generateTemplateMosaic(abAbacus    *abacus,
                       tgPosition  *utgpos,
                       uint32       numfrags,
                       uint32      &tiglen,
                       char        *tigseq) {

  for (uint32 i=0; i<numfrags; i++) {
    abSequence  *seq      = abacus->getSequence(i);
    char        *fragment = seq->getBases();
    uint32       readLen  = seq->length();

    uint32       start    = utgpos[i].min();
    uint32       end      = utgpos[i].max();

    if (start > tiglen) {
      fprintf(stderr, "WARNING: reset start  from " F_U32 " to " F_U32 "\n", start, tiglen-1);
      start = tiglen - 1;
    }

    if (end - start > readLen) {
      fprintf(stderr, "WARNING: reset end    from " F_U32 " to " F_U32 "\n", end, start+readLen);
      end = start + readLen;
    }

    if (end > tiglen) {
      fprintf(stderr, "WARNING: truncate end from " F_U32 " to " F_U32 "\n", end, tiglen-1);
      end = tiglen - 1;
    }

    //  Read aligns from position start to end.  Skip ahead until we find unset bases.

    uint32 cur = start;
    while ((cur < end) && (tigseq[cur] != 'N'))
      cur++;

#if 1
    if (cur < end)
      fprintf(stderr, "generatePBDAG()-- template from %7d to %7d comes from read %3d id %6d bases (%5d %5d) nominally %6d %6d)\n",
              cur, end, i, seq->gkpIdent(),
              cur - start,
              end - start,
              utgpos[i].min(),
              utgpos[i].max());
#endif

    for (uint32 j=cur; j<end; j++)
      tigseq[j] = fragment[j - start];

    tigseq[end] = 0;
    tiglen = end;
  }
}



void
generateTemplateStitch(abAbacus    *abacus,
                       tgPosition  *utgpos,
                       uint32       numfrags,
                       uint32      &tiglen,
                       char        *tigseq) {
  int32   minOlap  = 500;
  bool    verbose  = false;

  //  Initialize, copy the first read.

  uint32       rid      = 0;

  abSequence  *seq      = abacus->getSequence(rid);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();
  int32        olapLen  = 0;

  if (verbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "COPY READ read #%d %d (len=%d to %d-%d)\n",
            0, utgpos[0].ident(), readLen, utgpos[0].min(), utgpos[0].max());
  }

  for (uint32 ii=0; ii<readLen; ii++)
    tigseq[tiglen++] = fragment[ii];

  tigseq[tiglen] = 0;

  uint32       ePos = utgpos[0].max();   //  Expected end of template, from bogart supplied positions.


  //  Find the next read that has some minimum overlap and a large extension, copy that into the template.

    //  Align read to template.  Expecting to find alignment:
    //
    //        template ---------------
    //        read             ---------------
    //                               ^
    //
    //  All we need to find is where the template ends on the read, the ^ above.  We know the
    //  expected size of the overlap, and can extract those bases from the template and look for a
    //  full alignment to the read.
    //
    //  We'll align 80% of the expected overlap to the read, allowing a 20% buffer on either end.
    //
    //                        |  +-80% expected overlap size
    //                        |  |     +-fPos
    //                        v  v     v
    //        template ----------(-----)
    //        read            (-----------)------
    //

  while (rid < numfrags) {
    uint32 nr = 0;  //  Next read
    uint32 nm = 0;  //  Next read maximum position

    //  Pick the next read as the one with the longest extension from all with some minimum overlap
    //  to the template

    for (uint32 ii=rid+1; (ii < numfrags) && (utgpos[ii].min() + minOlap < ePos); ii++) {

      if (utgpos[ii].max() < ePos) {
        if (verbose)
          fprintf(stderr, "SKIP read #%d/%d %d %d-%d contained\n", ii, numfrags, utgpos[ii].ident(), utgpos[ii].min(), utgpos[ii].max());
        continue;
      }

      if (verbose)
        fprintf(stderr, "TEST read #%d/%d %d %d-%d\n", ii, numfrags, utgpos[ii].ident(), utgpos[ii].min(), utgpos[ii].max());

      if ((nm < utgpos[ii].max()) && (ePos < utgpos[ii].max())) {
        nr = ii;
        nm = utgpos[ii].max();
      }
    }

    if (nr == 0) {
      if (verbose)
        fprintf(stderr, "NO MORE READS TO ALIGN\n");
      break;
    }

    assert(nr != 0);

    rid      = nr;        //  We'll place read 'nr' in the template.

    seq      = abacus->getSequence(rid);
    fragment = seq->getBases();
    readLen  = seq->length();
    olapLen  = ePos - utgpos[nr].min();

    assert(olapLen > 0);

    int32  readBgn;
    int32  readEnd;

    EdlibAlignResult result;
    bool             aligned       = false;

    double           templateSize  = 0.80;
    double           extensionSize = 0.20;

  alignAgain:
    int32            templateLen  = (int32)ceil(olapLen * templateSize);    //  Extract 80% of the expected overlap size
    int32            extensionLen = (int32)ceil(olapLen * extensionSize);   //  Extend read by 20% of the expected overlap size

    readBgn = 0;
    readEnd = olapLen + extensionLen;

    if (readEnd > readLen)
      readEnd = readLen;

    if (verbose) {
      fprintf(stderr, "\n");
      fprintf(stderr, "TRY ALIGN template %d-%d (len=%d) to read #%d %d %d-%d (len=%d actual=%d at %d-%d)\n",
              tiglen - templateLen, tiglen, templateLen,
              nr, utgpos[nr].ident(), readBgn, readEnd, readEnd - readBgn, readLen, utgpos[nr].min(), utgpos[nr].max());
    }

    result = edlibAlign(tigseq + tiglen - templateLen, templateLen,
                        fragment, readEnd - readBgn,
                        edlibNewAlignConfig(olapLen * 0.06, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    //  We're expecting the template to align inside the read.
    //  If the alignment bumps up against the start of the read, we have too much template.

    bool   tryAgain = false;

    bool   noResult      = (result.numLocations == 0);
    bool   gotResult     = (result.numLocations  > 0);

    bool   hitTheStart   = (gotResult) && (result.startLocations[0] == 0);

    bool   hitTheEnd     = (gotResult) && (result.endLocations[0] + 1 == readEnd - readBgn);
    bool   moreToExtend  = (readEnd < readLen);

    if (noResult)
      if (verbose)
        fprintf(stderr, "FAILED to align\n");

    if ((noResult) || (hitTheStart)) {
      if ((verbose) && (gotResult))
        fprintf(stderr, "FAILED to align - startLocation = %d\n", result.startLocations[0]);
      tryAgain = true;
      templateSize -= 0.10;
    }

    //  If the alignment bumps up against the end of the read, we don't have enough read.

    if ((noResult) || (hitTheEnd && moreToExtend)) {
      if ((verbose) && (gotResult))
        fprintf(stderr, "FAILED to align - endLocation = %d (readEnd = %d readBgn = %d)\n", result.endLocations[0], readEnd, readBgn);
      tryAgain = true;
      extensionSize += 0.10;
    }

    if (tryAgain) {
      edlibFreeAlignResult(result);
      goto alignAgain;
    }

    readBgn = result.startLocations[0];     //  Expected to be zero
    readEnd = result.endLocations[0] + 1;   //  Where we need to start copying the read

    edlibFreeAlignResult(result);

    if (verbose)
      fprintf(stderr, "Aligned template %d-%d to read %u %d-%d; copy read %d-%d to template.\n", tiglen - templateLen, tiglen, nr, readBgn, readEnd, readEnd, readLen);

    for (uint32 ii=readEnd; ii<readLen; ii++)
      tigseq[tiglen++] = fragment[ii];

    tigseq[tiglen] = 0;

    if (verbose) {
      fprintf(stderr, "Template now length %d\n", tiglen);
      fprintf(stderr, "Reset ePos to %d\n", utgpos[rid].max());
    }

    ePos = utgpos[rid].max();
  }
}



bool
alignEdLib(dagAlignment      &aln,
           tgPosition        &utgpos,
           char              *fragment,
           uint32             fragmentLength,
           char              *tigseq,
           uint32             tiglen,
           double             lengthScale,
           double             errorRate,
           bool               normalize,
           bool               verbose) {

  EdlibAlignResult align;

  int32   padding        = (int32)ceil(fragmentLength * 0.10);
  double  bandErrRate    = errorRate / 2;
  bool    aligned        = false;
  double  alignedErrRate = 0.0;

  //  Decide on where to align this read.

  //  But, the utgpos positions are largely bogus, especially at the end of the tig.  utgcns (the
  //  original) used to track positions of previously placed reads, find an overlap beterrn this
  //  read and the last read, and use that info to find the coordinates for the new read.  That was
  //  very complicated.  Here, we just linearly scale.

  int32  tigbgn = max((int32)0,      (int32)floor(lengthScale * utgpos.min() - padding));
  int32  tigend = min((int32)tiglen, (int32)floor(lengthScale * utgpos.max() + padding));

  if (verbose)
    fprintf(stderr, "alignEdLib()-- align read %7u eRate %.4f at %9d-%-9d", utgpos.ident(), bandErrRate, tigbgn, tigend);

  //  This occurs if we don't lengthScale the positions.

  if (tigend < tigbgn)
    fprintf(stderr, "alignEdLib()-- ERROR: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
            tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
  assert(tigend > tigbgn);

  //  Align!  If there is an alignment, compute error rate and declare success if acceptable.

  align = edlibAlign(fragment, fragmentLength,
                     tigseq + tigbgn, tigend - tigbgn,
                     edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

  if (align.alignmentLength > 0) {
    alignedErrRate = (double)align.editDistance / align.alignmentLength;
    aligned        = (alignedErrRate <= errorRate);
    if (verbose)
      fprintf(stderr, " - ALIGNED %.4f at %9d-%-9d\n", alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
  } else {
    if (verbose)
      fprintf(stderr, "\n");
  }

  for (uint32 ii=0; ((ii < 4) && (aligned == false)); ii++) {
    tigbgn = max((int32)0,      tigbgn - 2 * padding);
    tigend = min((int32)tiglen, tigend + 2 * padding);

    bandErrRate += errorRate / 2;

    edlibFreeAlignResult(align);

    if (verbose)
      fprintf(stderr, "alignEdLib()--                    eRate %.4f at %9d-%-9d", bandErrRate, tigbgn, tigend);

    align = edlibAlign(fragment, strlen(fragment),
                       tigseq + tigbgn, tigend - tigbgn,
                       edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (align.alignmentLength > 0) {
      alignedErrRate = (double)align.editDistance / align.alignmentLength;
      aligned        = (alignedErrRate <= errorRate);
      if (verbose)
        fprintf(stderr, " - ALIGNED %.4f at %9d-%-9d\n", alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
    } else {
      if (verbose)
        fprintf(stderr, "\n");
    }
  }

  if (aligned == false) {
    edlibFreeAlignResult(align);
    return(false);
  }

  char *tgtaln = new char [align.alignmentLength+1];
  char *qryaln = new char [align.alignmentLength+1];

  memset(tgtaln, 0, sizeof(char) * (align.alignmentLength+1));
  memset(qryaln, 0, sizeof(char) * (align.alignmentLength+1));

  edlibAlignmentToStrings(align.alignment,               //  Alignment
                          align.alignmentLength,         //    and length
                          align.startLocations[0],       //  tgtStart
                          align.endLocations[0]+1,       //  tgtEnd
                          0,                             //  qryStart
                          fragmentLength,                //  qryEnd
                          tigseq + tigbgn,               //  tgt sequence
                          fragment,                      //  qry sequence
                          tgtaln,                   //  output tgt alignment string
                          qryaln);                  //  output qry alignment string

  //  Populate the output.  AlnGraphBoost does not handle mismatch alignments, at all, so convert
  //  them to a pair of indel.

  uint32 nMatch = 0;

  for (uint32 ii=0; ii<align.alignmentLength; ii++)   //  Edlib guarantees aln[alignmentLength] == 0.
    if ((tgtaln[ii] != '-') &&
        (qryaln[ii] != '-') &&
        (tgtaln[ii] != qryaln[ii]))
      nMatch++;

  aln.start  = tigbgn + align.startLocations[0] + 1;   //  AlnGraphBoost expects 1-based positions.
  aln.end    = tigbgn + align.endLocations[0] + 1;     //  EdLib returns 0-based positions.

  aln.qstr   = new char [align.alignmentLength + nMatch + 1];
  aln.tstr   = new char [align.alignmentLength + nMatch + 1];

  for (uint32 ii=0, jj=0; ii<align.alignmentLength; ii++) {
    char  tc = tgtaln[ii];
    char  qc = qryaln[ii];

    if ((tc != '-') &&
        (qc != '-') &&
        (tc != qc)) {
      aln.tstr[jj] = '-';   aln.qstr[jj] = qc;    jj++;
      aln.tstr[jj] = tc;    aln.qstr[jj] = '-';   jj++;
    } else {
      aln.tstr[jj] = tc;    aln.qstr[jj] = qc;    jj++;
    }

    aln.length = jj;
  }

  aln.qstr[aln.length] = 0;
  aln.tstr[aln.length] = 0;

  delete [] tgtaln;
  delete [] qryaln;

  edlibFreeAlignResult(align);

  if (aln.end > tiglen)
    fprintf(stderr, "ERROR:  alignment from %d to %d, but tiglen is only %d\n", aln.start, aln.end, tiglen);
  assert(aln.end <= tiglen);

  return(true);
}



void
realignReads() {

#ifdef REALIGN
  // update positions, this requires remapping but this time to the final consensus, turned off for now
  uint32 minPos = cns.size();
  uint32 maxPos = 0;

#pragma omp parallel for schedule(dynamic)
  for (uint32 i=0; i<numfrags; i++) {
    abSequence  *seq     = abacus->getSequence(i);

    uint32 bandTolerance = (int32)round((double)(seq->length() * errorRate)) * 2;
    uint32 maxExtend     = (int32)round((double)seq->length() * 0.01) + 1;
    int32  padding       = bandTolerance;
    uint32 start         = max((int32)0, (int32)utgpos[i].min() - padding);
    uint32 end           = min((int32)cns.size(), (int32)utgpos[i].max() + padding);

    EdlibAlignResult align = edlibAlign(seq->getBases(), seq->length()-1, cns.c_str()+start, end-start+1,  edlibNewAlignConfig(bandTolerance, EDLIB_MODE_HW, EDLIB_TASK_LOC));
    if (align.numLocations > 0) {
      cnspos[i].setMinMax(align.startLocations[0]+start, align.endLocations[0]+start+1);
      // when we are very close to end extend
      if (cnspos[i].max() < cns.size() && cns.size() - cnspos[i].max() <= maxExtend && (align.editDistance + cns.size() - cnspos[i].max()) < bandTolerance) {
        cnspos[i].setMinMax(cnspos[i].min(), cns.size());
      }
#pragma omp critical (trackMin)
      if (cnspos[i].min() < minPos) minPos = cnspos[i].min();
#pragma omp critical (trackMax)
      if (cnspos[i].max() > maxPos) maxPos = cnspos[i].max();
    } else {
    }
    edlibFreeAlignResult(align);
  }
  memcpy(tig->getChild(0), cnspos, sizeof(tgPosition) * numfrags);

  // trim consensus if needed
  if (maxPos < cns.size())
    cns = cns.substr(0, maxPos);

  assert(minPos == 0);
  assert(maxPos == cns.size());
#endif
}



bool
unitigConsensus::generatePBDAG(char                       aligner,
                               bool                       normalize,
                               tgTig                     *tig_,
                               map<uint32, gkRead *>     *inPackageRead_,
                               map<uint32, gkReadData *> *inPackageReadData_) {

  bool  verbose = false;

  tig      = tig_;
  numfrags = tig->numberOfChildren();

  if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
    fprintf(stderr, "generatePBDAG()-- Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    return(false);
  }

  //  First we need to load into Unitig data structure the quick cns

  uint32  tiglen = 0;
  char   *tigseq = new char [2 * tig->_layoutLen + 1];

  memset(tigseq, 'N', sizeof(char) * 2 * tig->_layoutLen);

  tigseq[2 * tig->_layoutLen] = 0;

  //  Build a quick consensus to align to.

  fprintf(stderr, "Generating template.\n");

  //generateTemplateMosaic(abacus, utgpos, numfrags, tiglen, tigseq);
  generateTemplateStitch(abacus, utgpos, numfrags, tiglen, tigseq);

  uint32  pass = 0;
  uint32  fail = 0;

  for (uint32 jj=0; jj<tiglen; jj++)
    if (tigseq[jj] == 'N')
      fprintf(stdout, "generatePBDAG()-- WARNING: template position %u not defined.\n", jj);

  assert(tigseq[tiglen] == 0);

  fprintf(stderr, "Generated template of length %d\n", tiglen);

  //  Compute alignments of each sequence in parallel

  fprintf(stderr, "Aligning reads.\n");

  dagAlignment *aligns = new dagAlignment [numfrags];

#pragma omp parallel for schedule(dynamic)
  for (uint32 ii=0; ii<numfrags; ii++) {
    abSequence  *seq      = abacus->getSequence(ii);
    bool         aligned  = false;

    assert(aligner == 'E');  //  Maybe later we'll have more than one aligner again.

    aligned = alignEdLib(aligns[ii],
                         utgpos[ii],
                         seq->getBases(), seq->length(),
                         tigseq, tiglen,
                         (double)tiglen / tig->_layoutLen,
                         errorRate,
                         normalize,
                         verbose);

    if (aligned == false) {
      if (verbose)
        fprintf(stderr, "generatePBDAG()--    read %7u FAILED\n", utgpos[ii].ident());

      fail++;

      continue;
    }

    pass++;
  }

  fprintf(stderr, "Finished aligning reads.  %d failed, %d passed.\n", fail, pass);

  //  Construct the graph from the alignments.  This is not thread safe.

  fprintf(stderr, "Constructing graph\n");

  AlnGraphBoost ag(string(tigseq, tiglen));

  for (uint32 ii=0; ii<numfrags; ii++) {
    cnspos[ii].setMinMax(aligns[ii].start, aligns[ii].end);

    if ((aligns[ii].start == 0) &&
        (aligns[ii].end   == 0))
      continue;

    ag.addAln(aligns[ii]);

    aligns[ii].clear();
  }

  delete [] aligns;

  fprintf(stderr, "Merging graph\n");

  //  Merge the nodes and call consensus
  ag.mergeNodes();

  fprintf(stderr, "Calling consensus\n");

  std::string cns = ag.consensus(1);

  //  Realign reads to get precise endpoints

  realignReads();

  //  Save consensus

  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, 0, tig->_gappedMax, (uint32) cns.length() + 1, resizeArray_doNothing);

  std::string::size_type len = 0;

  for (len=0; len<cns.size(); len++) {
    tig->_gappedBases[len] = cns[len];
    tig->_gappedQuals[len] = CNS_MIN_QV;
  }

  //  Terminate the string.

  tig->_gappedBases[len] = 0;
  tig->_gappedQuals[len] = 0;
  tig->_gappedLen        = len;
  tig->_layoutLen        = len;

  assert(len < tig->_gappedMax);

  return(true);
}



bool
unitigConsensus::generateQuick(tgTig                     *tig_,
                               map<uint32, gkRead *>     *inPackageRead_,
                               map<uint32, gkReadData *> *inPackageReadData_) {
  tig      = tig_;
  numfrags = tig->numberOfChildren();

  if (initialize(inPackageRead_, inPackageReadData_) == FALSE) {
    fprintf(stderr, "generatePBDAG()-- Failed to initialize for tig %u with %u children\n", tig->tigID(), tig->numberOfChildren());
    return(false);
  }

  //  First we need to load into Unitig data structure the quick cns

  uint32  tiglen = 0;
  char   *tigseq = new char [2 * tig->_layoutLen + 1];

  memset(tigseq, 'N', sizeof(char) * 2 * tig->_layoutLen);

  tigseq[2 * tig->_layoutLen] = 0;

  //  Build a quick consensus to align to.

  fprintf(stderr, "Generating template.\n");

  //generateTemplateMosaic(abacus, utgpos, numfrags, tiglen, tigseq);
  generateTemplateStitch(abacus, utgpos, numfrags, tiglen, tigseq);

  //
  //  The above and below came from generatePBDAG(), which should be modified to handle 'quick'.
  //  generagePBDAG() has a bunch of other stuff here.
  //

  //  Save consensus

  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, 0, tig->_gappedMax, tiglen + 1, resizeArray_doNothing);

  for (uint32 ii=0; ii<tiglen; ii++) {
    tig->_gappedBases[ii] = tigseq[ii];
    tig->_gappedQuals[ii] = CNS_MIN_QV;
  }

  //  Terminate the string.

  tig->_gappedBases[tiglen] = 0;
  tig->_gappedQuals[tiglen] = 0;
  tig->_gappedLen           = tiglen;
  tig->_layoutLen           = tiglen;

  return(true);
}



int
unitigConsensus::initialize(map<uint32, gkRead *>     *inPackageRead,
                            map<uint32, gkReadData *> *inPackageReadData) {

  int32 num_columns = 0;
  //int32 num_bases   = 0;

  if (numfrags == 0) {
    fprintf(stderr, "utgCns::initialize()-- unitig has no children.\n");
    return(false);
  }

  utgpos = new tgPosition [numfrags];
  cnspos = new tgPosition [numfrags];

  memcpy(utgpos, tig->getChild(0), sizeof(tgPosition) * numfrags);
  memcpy(cnspos, tig->getChild(0), sizeof(tgPosition) * numfrags);

  traceLen   = 0;
  trace      = new int32 [2 * AS_MAX_READLEN];

  traceABgn  = 0;
  traceBBgn  = 0;

  memset(trace, 0, sizeof(int32) * 2 * AS_MAX_READLEN);

  abacus     = new abAbacus();

  //  Clear the cnspos position.  We use this to show it's been placed by consensus.
  //  Guess the number of columns we'll end up with.
  //  Initialize abacus with the reads.

  for (int32 i=0; i<numfrags; i++) {
    cnspos[i].setMinMax(0, 0);

    num_columns  = (utgpos[i].min() > num_columns) ? utgpos[i].min() : num_columns;
    num_columns  = (utgpos[i].max() > num_columns) ? utgpos[i].max() : num_columns;

    abacus->addRead(gkpStore,
                    utgpos[i].ident(),
                    utgpos[i]._askip, utgpos[i]._bskip,
                    utgpos[i].isReverse(),
                    inPackageRead,
                    inPackageReadData);
  }

  //  Check for duplicate reads

  {
    set<uint32>  dupFrag;

    for (uint32 i=0; i<numfrags; i++) {
      if (utgpos[i].isRead() == false) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is not a read.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      if (dupFrag.find(utgpos[i].ident()) != dupFrag.end()) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is a duplicate.\n",
                tig->tigID(), utgpos[i].ident());
        return(false);
      }

      dupFrag.insert(utgpos[i].ident());
    }
  }

  //  Initialize with the first read.

  abacus->applyAlignment(0, 0, 0, NULL, 0);

  //  And set the placement of the first read.

  cnspos[0].setMinMax(0, abacus->numberOfColumns());

  return(true);
}



int
unitigConsensus::computePositionFromAnchor(void) {

  assert(piid == -1);

  uint32 anchor = utgpos[tiid].anchor();

  if (anchor == 0)
    //  No anchor?!  Damn.
    goto computePositionFromAnchorFail;

  for (piid = tiid-1; piid >= 0; piid--) {
    abSequence *aseq = abacus->getSequence(piid);

    if (anchor != aseq->gkpIdent())
      //  Not the anchor.
      continue;

    if ((cnspos[piid].min() == 0) &&
        (cnspos[piid].max() == 0))
      //  Is the anchor, but that isn't placed.
      goto computePositionFromAnchorFail;

    if ((utgpos[piid].max() < utgpos[tiid].min()) ||
        (utgpos[tiid].max() < utgpos[piid].min())) {
      //  Is the anchor, and anchor is placed, but the anchor doesn't agree with the placement.
      if (showPlacement())
        fprintf(stderr, "computePositionFromAnchor()-- anchor %d at utg %d,%d doesn't agree with my utg %d,%d.  FAIL\n",
                anchor,
                utgpos[piid].min(), utgpos[piid].max(),
                utgpos[tiid].min(), utgpos[tiid].max());
      goto computePositionFromAnchorFail;
    }

    //  Scale the hangs by the change in the anchor size between bogart and consensus.

#if 0
    double   anchorScale = (double)(cnspos[piid].max() - cnspos[piid].min()) / (double)(utgpos[piid].max() - utgpos[piid].min());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- hangs %d,%d -- scale %f -- final hangs %.0f,%.0f\n",
              utgpos[tiid].ident(),
              utgpos[piid].ident(),
              utgpos[tiid].aHang(),
              utgpos[tiid].bHang(),
              anchorScale,
              utgpos[tiid].aHang() * anchorScale,
              utgpos[tiid].bHang() * anchorScale);

    cnspos[tiid].setMinMax(cnspos[piid].min() + utgpos[tiid].aHang() * anchorScale,
                           cnspos[piid].max() + utgpos[tiid].bHang() * anchorScale);

    //  Hmmm, but if we shrank the read too much, add back in some of the length.  We want to end up
    //  with the read scaled by anchorScale, and centered on the hangs.

    int32   fragmentLength = utgpos[tiid].max() - utgpos[tiid].min();

    if ((cnspos[tiid].min() >= cnspos[tiid].max()) ||
        (cnspos[tiid].max() - cnspos[tiid].min() < 0.75 * fragmentLength)) {
      int32  center = (cnspos[tiid].min() + cnspos[tiid].max()) / 2;

      if (showPlacement()) {
        fprintf(stderr, "computePositionFromAnchor()--  frag %u in anchor %u -- too short.  reposition around center %d with adjusted length %.0f\n",
                utgpos[tiid].ident(),
                utgpos[piid].ident(),
                center, fragmentLength * anchorScale);
      }

      cnspos[tiid].setMinMax(center - fragmentLength * anchorScale / 2,
                             center + fragmentLength * anchorScale / 2);

      //  We seem immune to having a negative position.  We only use this to pull out a region from
      //  the partial consensus to align to.
      //
      //if (cnspos[tiid].min() < 0) {
      //  cnspos[tiid].min() = 0;
      //  cnspos[tiid].max() = fragmentLength * anchorScale;
      //}
    }
#else
    assert(0 <= utgpos[tiid].aHang());

    uint32  bgn = abacus->getColumn(piid, cnspos[piid].min() + utgpos[tiid].aHang() - cnspos[piid].min());
    uint32  end = abacus->getColumn(piid, cnspos[piid].max() + utgpos[tiid].bHang() - cnspos[piid].min());

    cnspos[tiid].setMinMax(bgn, end);
#endif

    assert(cnspos[tiid].min() < cnspos[tiid].max());

    if (showPlacement())
      fprintf(stderr, "computePositionFromAnchor()-- anchor %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
              anchor,
              cnspos[piid].min(), cnspos[piid].max(),
              cnspos[tiid].min(), cnspos[tiid].max(),
              abacus->numberOfColumns());
    return(true);
  }

 computePositionFromAnchorFail:
  cnspos[tiid].setMinMax(0, 0);

  piid = -1;

  return(false);
}



int
unitigConsensus::computePositionFromLayout(void) {
  int32   thickestLen = 0;

  assert(piid == -1);

  //  Find the thickest qiid overlap to any cnspos fragment
  for (int32 qiid = tiid-1; qiid >= 0; qiid--) {
    if ((utgpos[tiid].min() < utgpos[qiid].max()) &&
        (utgpos[tiid].max() > utgpos[qiid].min()) &&
        ((cnspos[qiid].min() != 0) ||
         (cnspos[qiid].max() != 0))) {
      cnspos[tiid].setMinMax(cnspos[qiid].min() + utgpos[tiid].min() - utgpos[qiid].min(),
                             cnspos[qiid].max() + utgpos[tiid].max() - utgpos[qiid].max());

      //  This assert triggers.  It results in 'ooo' below being negative, and we
      //  discard this overlap anyway.
      //
      //assert(cnspos[tiid].min() < cnspos[tiid].max());

      int32 ooo = MIN(cnspos[tiid].max(), abacus->numberOfColumns()) - cnspos[tiid].min();

#if 1
      if (showPlacement())
        fprintf(stderr, "computePositionFromLayout()-- layout %d at utg %d,%d cns %d,%d --> utg %d,%d cns %d,%d -- overlap %d\n",
                utgpos[qiid].ident(),
                utgpos[qiid].min(), utgpos[qiid].max(), cnspos[qiid].min(), cnspos[qiid].max(),
                utgpos[tiid].min(), utgpos[tiid].max(), cnspos[tiid].min(), cnspos[tiid].max(),
                ooo);
#endif

      //  Occasionally we see an overlap in the original placement (utgpos overlap) by after
      //  adjusting our fragment to the consensus position, we no longer have an overlap.  This
      //  seems to be caused by a bad original placement.
      //
      //  Example:
      //  utgpos[a] = 13480,14239    cnspos[a] = 13622,14279
      //  utgpos[b] = 14180,15062
      //
      //  Our placement is 200bp different at the start, but close at the end.  When we compute the
      //  new start placement, it starts after the end of the A read -- the utgpos say the B read
      //  starts 700bp after the A read, which is position 13622 + 700 = 14322....50bp after A ends.

      if ((cnspos[tiid].min() < abacus->numberOfColumns()) &&
          (thickestLen < ooo)) {
        thickestLen = ooo;

        assert(cnspos[tiid].min() < cnspos[tiid].max());  //  But we'll still assert cnspos is ordered correctly.

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].min();
        int32 bhang = cnspos[tiid].max() - abacus->numberOfColumns();

        piid  = qiid;
      }
    }
  }

  //  If we have a VALID thickest placement, use that (recompute the placement that is likely
  //  overwritten -- ahang, bhang and piid are still correct).

  if (thickestLen >= minOverlap) {
    assert(piid != -1);

    cnspos[tiid].setMinMax(cnspos[piid].min() + utgpos[tiid].min() - utgpos[piid].min(),
                           cnspos[piid].max() + utgpos[tiid].max() - utgpos[piid].max());

    assert(cnspos[tiid].min() < cnspos[tiid].max());

    if (showPlacement())
      fprintf(stderr, "computePositionFromLayout()-- layout %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
              utgpos[piid].ident(),
              cnspos[piid].min(), cnspos[piid].max(),
              cnspos[tiid].min(), cnspos[tiid].max(),
              abacus->numberOfColumns());

    return(true);
  }

  cnspos[tiid].setMinMax(0, 0);

  piid = -1;

  return(false);
}



//  Occasionally we get a fragment that just refuses to go in the correct spot.  Search for the
//  correct placement in all of consensus, update ahang,bhang and retry.
//
//  We don't expect to have big negative ahangs, and so we don't allow them.  To unlimit this, use
//  "-fragmentLen" instead of the arbitrary cutoff below.
int
unitigConsensus::computePositionFromAlignment(void) {

  assert(piid == -1);

  int32        minlen      = minOverlap;
  int32        ahanglimit  = -10;

  abSequence  *seq         = abacus->getSequence(tiid);
  char        *fragment    = seq->getBases();
  int32        fragmentLen = seq->length();

  bool         foundAlign  = false;

  //
  //  Try NDalign.
  //

  if (foundAlign == false) {

    if (oaPartial == NULL)
      oaPartial = new NDalign(pedLocal, errorRate, 17);  //  partial allowed!

    oaPartial->initialize(0, abacus->bases(), abacus->numberOfColumns(), 0, abacus->numberOfColumns(),
                          1, fragment,        fragmentLen,               0, fragmentLen,
                          false);

    if ((oaPartial->findMinMaxDiagonal(minOverlap) == true) &&
        (oaPartial->findSeeds(false)               == true) &&
        (oaPartial->findHits()                     == true) &&
        (oaPartial->chainHits()                    == true) &&
        (oaPartial->processHits()                  == true)) {

      cnspos[tiid].setMinMax(oaPartial->abgn(), oaPartial->aend());

      //fprintf(stderr, "computePositionFromAlignment()-- cnspos[%3d] mid %d %d,%d (from NDalign)\n", tiid, utgpos[tiid].ident(), cnspos[tiid].min(), cnspos[tiid].max());

      foundAlign = true;
    }
  }

  //
  //  Fail.
  //

  if (foundAlign == false) {
    cnspos[tiid].setMinMax(0, 0);
    piid = -1;

    if (showAlgorithm())
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail (no alignment).\n");
    return(false);
  }

  //  From the overlap and existing placements, find the thickest overlap, to set the piid and
  //  hangs, then reset the original placement based on that anchors original placement.
  //
  //  To work with fixFailures(), we need to scan the entire fragment list.  This isn't so bad,
  //  really, since before we were scanning (on average) half of it.

  assert(cnspos[tiid].min() < cnspos[tiid].max());

  int32   thickestLen = 0;

  for (int32 qiid = numfrags-1; qiid >= 0; qiid--) {
    if ((tiid != qiid) &&
        (cnspos[tiid].min() < cnspos[qiid].max()) &&
        (cnspos[tiid].max() > cnspos[qiid].min())) {
      int32 ooo = (MIN(cnspos[tiid].max(), cnspos[qiid].max()) -
                   MAX(cnspos[tiid].min(), cnspos[qiid].min()));

      if (thickestLen < ooo) {
        thickestLen = ooo;

        int32 ovl   = ooo;
        int32 ahang = cnspos[tiid].min();
        int32 bhang = cnspos[tiid].max() - abacus->numberOfColumns();

        piid  = qiid;
      }
    }
  }

  //  No thickest?  Dang.

  if (thickestLen == 0) {
    cnspos[tiid].setMinMax(0, 0);
    piid = -1;
    if (showAlgorithm())
      fprintf(stderr, "computePositionFromAlignment()-- Returns fail (no thickest).\n");
    return(false);
  }

  //  Success, yay!

  assert(piid != -1);

  if (showPlacement())
    fprintf(stderr, "computePositionFromAlignment()-- layout %d at %d,%d --> beg,end %d,%d (tigLen %d)\n",
            utgpos[piid].ident(),
            cnspos[piid].min(), cnspos[piid].max(),
            cnspos[tiid].min(), cnspos[tiid].max(),
            abacus->numberOfColumns());

  return(true);
}


void
unitigConsensus::generateConsensus(tgTig *tig) {

  abacus->recallBases(true);  //  Do one last base call, using the full works.

  abacus->refine(abAbacus_Smooth);
  abacus->mergeColumns(true);

  abacus->refine(abAbacus_Poly_X);
  abacus->mergeColumns(true);

  abacus->refine(abAbacus_Indel);
  abacus->mergeColumns(true);

  abacus->recallBases(true);  //  The bases are possibly all recalled, depending on the above refinements keeping things consistent.
  //abacus->refreshColumns();    //  Definitely needed, this copies base calls into _cnsBases and _cnsQuals.

  //  Copy the consensus and positions into the tig.

  abacus->getConsensus(tig);
  abacus->getPositions(tig);

  //  While we have fragments in memory, compute the microhet probability.  Ideally, this would be
  //  done in CGW when loading unitigs (the only place the probability is used) but the code wants
  //  to load sequence and quality for every fragment, and that's too expensive.
}



//  Update the position of each fragment in the consensus sequence.
//  Update the anchor/hang of the fragment we just placed.
void
unitigConsensus::refreshPositions(void) {

  for (int32 i=0; i<=tiid; i++) {
    if ((cnspos[i].min() == 0) &&
        (cnspos[i].max() == 0))
      //  Uh oh, not placed originally.
      continue;

    abColumn *fcol = abacus->readTofBead[i].column;
    abColumn *lcol = abacus->readTolBead[i].column;

    cnspos[i].setMinMax(fcol->position(),
                        lcol->position() + 1);

    assert(cnspos[i].min() >= 0);
    assert(cnspos[i].max() > cnspos[i].min());
  }

  if (piid >= 0)
    utgpos[tiid].setAnchor(utgpos[piid].ident(),
                           cnspos[tiid].min() - cnspos[piid].min(),
                           cnspos[tiid].max() - cnspos[piid].max());

  piid = -1;
}



//  Run abacus to rebuild the consensus sequence.  VERY expensive.
void
unitigConsensus::recomputeConsensus(bool display) {

  //abacus->recallBases(false);  //  Needed?  We should be up to date.

  abacus->refine(abAbacus_Smooth);
  abacus->mergeColumns(false);

  abacus->refine(abAbacus_Poly_X);
  abacus->mergeColumns(false);

  abacus->refine(abAbacus_Indel);
  abacus->mergeColumns(false);

  abacus->recallBases(false);  //  Possibly not needed.  If this is removed, the following refresh is definitely needed.
  //abacus->refreshColumns();    //  Definitely needed, this copies base calls into _cnsBases and _cnsQuals.

  refreshPositions();

  if (display)
    abacus->display(stderr);
}



//  This stub lets alignFragmnet() cleanup and return on alignment failures.  The original
//  implementation did the same thing with a goto to the end of the function.  Opening up the if
//  statements exposed variable declarations that prevented the goto from compiling.
//
bool
unitigConsensus::alignFragmentFailure(void) {
  cnspos[tiid].setMinMax(0, 0);
  piid = -1;

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- No alignment found.\n");

  return(false);
}


//  Generates an alignment of the current read to the partial consensus.
//  The primary output is a trace stored in the object data.

bool
unitigConsensus::alignFragment(bool forceAlignment) {

  assert((cnspos[tiid].min() != 0) || (cnspos[tiid].max() != 0));
  assert(piid != -1);

  assert(cnspos[tiid].min() < cnspos[tiid].max());

  abSequence *bSEQ    = abacus->getSequence(tiid);
  char       *fragSeq = bSEQ->getBases();
  int32       fragLen = bSEQ->length();

  //  Decide on how much to align.  Pick too little of consensus, and we leave some of the read
  //  unaligned.  Pick too much, and the read aligns poorly.
  //
  //  endTrim is trimmed from the 3' of the read.  This is the stuff we don't expect to align to consensus.
  //  bgnExtra is trimmed from the 5' of consensus.  Same idea as endTrim.
  //  endExtra is trimmed from the 3' of consensus.  Only for contained reads.
  //
  //  These values are adjusted later, in trimStep increments, based on the alignments returned.
  //
  //  Of the two choices, making the two Extra's small at the start is probably safer.  That case is
  //  easy to detect, and easy to fix.
  //

  //  The expectedAlignLen is almost always an underestimate.  Any gaps inserted will make the real
  //  alignment length longer.  This used to be multiplied by the error rate.
  int32  expectedAlignLen = cnspos[tiid].max() - cnspos[tiid].min();

  //  If the read is contained, the full read is aligned.
  //  Otherwise, an extra 1/32 of the align length is added for padding.
  int32  fragBgn = 0;
  int32  fragEnd = (cnspos[tiid].max() < abacus->numberOfColumns()) ? (fragLen) : (33 * expectedAlignLen / 32);

  if (fragEnd > fragLen)
    fragEnd = fragLen;

  //  Given the usual case of using an actual overlap to a read in the multialign to find the region
  //  to align to, we expect that region to be nearly perfect.  Thus, we shouldn't need to extend it
  //  much.  If anything, we'll need to extend the 3' end of the read.

  int32  bgnExtra = 10;    //  Start with a small 'extra' allowance, easy to make bigger.
  int32  endExtra = 10;    //

  int32  trimStep = max(10, expectedAlignLen / 50);  //  Step by 10 bases or do at most 50 steps.

  //  Find an alignment!

  bool  allowedToTrim = true;

  assert(abacus->bases()[abacus->numberOfColumns()] == 0);  //  Consensus must be NUL terminated
  assert(fragSeq[fragLen]                           == 0);  //  The read must be NUL terminated

 alignFragmentAgain:

  //  Truncate consensus and the read to prevent false alignments.

  if (cnspos[tiid].max() + endExtra > abacus->numberOfColumns())
    endExtra = abacus->numberOfColumns() - cnspos[tiid].max();

  int32 cnsBgn     = MAX(0, cnspos[tiid].min() - bgnExtra);   //  Start position in consensus
  int32 cnsEnd     = cnspos[tiid].max() + endExtra;           //  Truncation of consensus
  int32 cnsEndBase = abacus->bases()[cnsEnd];                 //  Saved base (if not truncated, it's the NUL byte at the end)

  char *aseq         = abacus->bases() + cnsBgn;
  char *bseq         = fragSeq;

  char  fragEndBase  = bseq[fragEnd];                         //  Saved base

  abacus->bases()[cnsEnd] = 0;                                //  Do the truncations.
  bseq[fragEnd]           = 0;

  //  Report!

  if (showAlgorithm())
    fprintf(stderr, "alignFragment()-- Allow bgnExtra=%d and endExtra=%d (cnsBgn=%d cnsEnd=%d cnsLen=%d) (fragBgn=0 fragEnd=%d fragLen=%d)\n",
            bgnExtra, endExtra, cnsBgn, cnsEnd, abacus->numberOfColumns(), fragEnd, fragLen);

  //  Create new aligner object.  'Global' in this case just means to not stop early, not a true global alignment.

  if (oaFull == NULL)
    oaFull = new NDalign(pedGlobal, errorRate, 17);

  oaFull->initialize(0, aseq, cnsEnd  - cnsBgn,   0, cnsEnd  - cnsBgn,
                     1, bseq, fragEnd - fragBgn,  0, fragEnd - fragBgn,
                     false);

  //  Generate a null hit, then align it and then realign, from both endpoints, and save the better
  //  of the two.

  if ((oaFull->makeNullHit() == true) &&
      (oaFull->processHits() == true)) {
    if (showAlignments())
      oaFull->display("utgCns::alignFragment()--", true);

    oaFull->realignBackward(showAlgorithm(), showAlignments());
    oaFull->realignForward (showAlgorithm(), showAlignments());
  }

  //  Restore the bases we removed to end the strings early.

  if (cnsEndBase)   abacus->bases()[cnsEnd] = cnsEndBase;
  if (fragEndBase)  bseq[fragEnd]           = fragEndBase;

  //  If no alignment, bail.

  if (oaFull->length() == 0)
    return(alignFragmentFailure());

  //
  //  Check quality and fail if it sucks.
  //

  bool  isBad = oaFull->scanDeltaForBadness(showAlgorithm(), showAlignments());

  //  Check for bad (under) trimming of input sequences.
  //
  //  If the alignment is bad, and we hit the start of the consensus sequence (or the end of
  //  same), chances are good that the aligner returned a (higher scoring) global alignment instead
  //  of a (lower scoring) local alignment.  Trim off some of the extension and try again.

  if ((allowedToTrim == true) && (isBad == true) && (oaFull->ahg5() == 0) && (bgnExtra > 0)) {
    int32  adj = (bgnExtra < trimStep) ? 0 : bgnExtra - trimStep;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment is bad, hit the trimmed start of consensus, decrease bgnExtra from %u to %u\n", bgnExtra, adj);

    bgnExtra = adj;
    goto alignFragmentAgain;
  }

  if ((allowedToTrim == true) && (isBad == true) && (oaFull->ahg3() == 0) && (endExtra > 0)) {
    int32  adj = (endExtra < trimStep) ? 0 : endExtra - trimStep;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment is bad, hit the trimmed end of consensus, decrease endExtra from %u to %u\n", endExtra, adj);

    endExtra = adj;
    goto alignFragmentAgain;
  }

  //  Check for bad (over) trimming of input sequences.  Bad if:
  //     we don't hit the start of the read, and we chopped the start of consensus.
  //     we don't hit the end   of the read, and we chopped the end   of consensus.
  //     we do    hit the end   of the read, and we chopped the end   of the read.
  //     (we don't chop the start of the read, so the fourth possible case never happens)

  allowedToTrim = false;  //  No longer allowed to reduce bgnExtra or endExtra.  We'd hit infinite loops otherwise.

  if ((oaFull->bhg5() > 0) && (cnsBgn > 0)) {
    int32  adj = bgnExtra + 2 * oaFull->bhg5();

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed start of consensus, increase bgnExtra from %u to %u\n", bgnExtra, adj);

    bgnExtra = adj;
    goto alignFragmentAgain;
  }

  if ((oaFull->bhg3() > 0) && (cnsEnd < abacus->numberOfColumns())) {
    int32  adj = endExtra + 2 * oaFull->bhg3();

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed end of consensus, increase endExtra from %u to %u\n", endExtra, adj);

    endExtra = adj;
    goto alignFragmentAgain;
  }

  if ((oaFull->bhg3() == 0) && (fragEnd < fragLen)) {
    int32  adj = (fragEnd + trimStep < fragLen) ? fragEnd + trimStep : fragLen;

    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- hit the trimmed end of the read, increase fragEnd from %d to %d\n", fragEnd, adj);

    fragEnd = adj;
    goto alignFragmentAgain;
  }

  //  If we get here, and it's still bad, well, not much we can do.

  if ((forceAlignment == false) && (isBad == true)) {
    if (showAlgorithm())
      fprintf(stderr, "utgCns::alignFragment()-- alignment bad after realigning\n");
    return(alignFragmentFailure());
  }

  if ((forceAlignment == false) && (oaFull->erate() > errorRate)) {
    if (showAlgorithm()) {
      fprintf(stderr, "utgCns::alignFragment()-- alignment is low quality: %f > %f\n",
              oaFull->erate(), errorRate);
      oaFull->display("utgCns::alignFragment()-- ", true);
    }
    return(alignFragmentFailure());
  }

  //  Otherwise, its a good alignment.  Process the trace to the 'consensus-format' and return true.
  //
  //  Set the begin points of the trace.  We probably need at least one of abgn and bbgn to be
  //  zero.  If both are nonzero, then we have a branch in the alignment.
  //
  //  If traceABgn is negative, we insert gaps into A before it starts, but I'm not sure how that works.
  //
  //  Overlap encoding:
  //    Add N matches or mismatches.  If negative, insert a base in the first sequence.  If
  //    positive, delete a base in the first sequence.
  //
  //  Consensus encoding:
  //    If negative, align (-trace - apos) bases, then add a gap in A.
  //    If positive, align ( trace - bpos) bases, then add a gap in B.
  //

  if (oaFull->abgn() > 0)   assert(oaFull->bbgn() == 0);  //  read aligned fully if consensus isn't
  if (oaFull->bbgn() > 0)   assert(oaFull->abgn() == 0);  //  read extends past the begin, consensus aligned fully
  if (oaFull->bbgn() > 0)   assert(cnsBgn == 0);          //  read extends past the begin, consensus not trimmed at begin


  traceABgn = cnsBgn + oaFull->abgn() - oaFull->bbgn();
  traceBBgn =            oaFull->bbgn();

  int32   apos = oaFull->abgn();
  int32   bpos = oaFull->bbgn();

  traceLen = 0;

  for (uint32 ii=0; ii<oaFull->deltaLen(); ii++, traceLen++) {
    if (oaFull->delta()[ii] < 0) {
      apos += -oaFull->delta()[ii] - 1;
      bpos += -oaFull->delta()[ii];

      trace[traceLen] = -apos - cnsBgn - 1;

    } else {
      apos +=  oaFull->delta()[ii];
      bpos +=  oaFull->delta()[ii] - 1;  // critical

      trace[traceLen] = bpos + 1;
    }
  }

  trace[traceLen] = 0;

  return(true);
}
