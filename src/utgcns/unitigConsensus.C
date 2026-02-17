
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "unitigConsensus.H"

#include "bits.H"

#include "Alignment.H"
#include "AlnGraphBoost.H"
#include "align.H"

#include "htslib/hts/sam.h"
//nclude <htslib/bgzf.h>
//nclude "kmers.H"
//nclude "merlin-globals.H"

#include <string>
#include <vector>
#include <algorithm>

// alignment parameter controlling retries and error rate scaling
#define ERROR_RATE_FACTOR 4
#define NUM_BANDS         2
#define MAX_RETRIES       ERROR_RATE_FACTOR * NUM_BANDS
#define TRIM_BP           500
#define MAX_PADDING       250
#define MAX_PADDING_MULT  25

// alignment parameters controlling when lower quality reads get promoted to be used
#define MIN_COV           10
#define MIN_COV_SIZE      5000

// alignment decrease for building template
// we are conservative when building the template but then allow a bit more freedom because we want to include as many reads as we can in the consensus
#define ALIGN_FACTOR      1.5

abSequence::abSequence(uint32  readID,
                       uint32  length,
                       char   *seq,
                       bool    isReverse) {
  _iid              = readID;
  _length           = length;
  _bases            = nullptr;

  if (length == 0)
    return;

  _bases            = new char  [_length + 1];

  //  Make a complement table

  char inv[256] = {0};

  inv['a'] = 't';  inv['A'] = 'T';
  inv['c'] = 'g';  inv['C'] = 'G';
  inv['g'] = 'c';  inv['G'] = 'C';
  inv['t'] = 'a';  inv['T'] = 'A';
  inv['n'] = 'n';  inv['N'] = 'N';
  inv['-'] = '-';

  //  Stash the bases/quals

  for (uint32 ii=0; ii<_length; ii++)
    assert((seq[ii] == 'A') ||
           (seq[ii] == 'C') ||
           (seq[ii] == 'G') ||
           (seq[ii] == 'T') ||
           (seq[ii] == 'N'));

  if (isReverse == false)
    for (uint32 ii=0, pp=0; ii<_length; ii++, pp++)
      _bases[pp] = seq[ii];
  else
    for (uint32 ii=_length, pp=0; ii-->0; pp++)
      _bases[pp] = inv[ seq[ii] ];

  _bases[_length] = 0;  //  NUL terminate the strings so we can use them in aligners.
}


abSequence::~abSequence() {
  delete [] _bases;
}



unitigConsensus::unitigConsensus(sqStore  *seqStore_,
                                 double    errorRate_,
                                 double    errorRateMax_,
                                 uint32    errorRateMaxID_,
                                 uint32    minOverlap_,
                                 uint32    minCoverage_) {
  _seqStore        = seqStore_;
  _minOverlap      = minOverlap_;
  _errorRate       = errorRate_;
  _errorRateMax    = errorRateMax_;
  _errorRateMaxID  = errorRateMaxID_;
  _minCoverage     = minCoverage_;
}


unitigConsensus::~unitigConsensus() {

  for (uint32 ss=0; ss<_sequencesLen; ss++)
    delete _sequences[ss];

  delete [] _sequences;
  delete [] _utgpos;
  delete [] _utghpc;
  delete [] _cnspos;
  delete [] _adjpos;
  delete [] _templateToCNS;
}



void
unitigConsensus::addRead(uint32     readID,
                         uint32     askip,
                         uint32     bskip,
                         bool       complemented,
                         u32toRead &reads) {

  //  Grab the read.  If there is no package, load the read from the store.  Otherwise, load the
  //  read from the package.  This REQUIRES that the package be in-sync with the unitig.  We fail
  //  otherwise.  Hey, it's used for debugging only...

  sqRead      *readToDelete = NULL;
  sqRead      *read         = NULL;

  if (reads.size() == 0) {
    readToDelete = new sqRead;
    read         = _seqStore->sqStore_getRead(readID, readToDelete);
  }

  else {
    read         = reads[readID];
  }

  if (read == NULL)
    fprintf(stderr, "Failed to load read %u\n", readID);
  assert(read != NULL);

  //  Grab seq/qlt from the read, offset to the proper begin and length.

  uint32  seqLen = read->sqRead_length() - askip - bskip;
  char   *seq    = read->sqRead_sequence()  + ((complemented == false) ? askip : bskip);

  //  Add it to our list.

  increaseArray(_sequences, _sequencesLen, _sequencesMax, 1);

  _sequences[_sequencesLen++] = new abSequence(readID, seqLen, seq, complemented);

  delete readToDelete;
}



void unitigConsensus::promoteLowQualUniques(uint32 tiglen, uint8_t minCov, uint32 minLen) {
  // now go through, compute coverage and flip unused reads to be promoted to consensus in regions w/low unique coverage
  uint32_t* cov = new uint32_t[tiglen + 1];
  memset(cov, 0, sizeof(uint32_t) * (tiglen + 1));
  #pragma omp parallel for schedule(dynamic)
  for (uint32 child = 0; child < _numReads; child++) {
     if (_utgpos[child]._numPlacement == 1 && !_utgpos[child].skipConsensus()) {
        for (uint32 j = _utgpos[child].min(); j < _utgpos[child].max(); j++) {
           if (j < tiglen) {
              #pragma omp atomic
              cov[j]++;
           }
        }
     }
  }

  #pragma omp parallel for schedule(dynamic)
  for (uint32 child = 0; child < _numReads; child++) {
     uint32 basesBelow = 0;
     for (uint32 j = _utgpos[child].min(); j < _utgpos[child].max() && j < tiglen; j++) {
        assert(j < tiglen);
        if (cov[j] <= minCov)
           basesBelow++;
     }
     if (basesBelow >= minLen && _utgpos[child]._numPlacement == 1 && _utgpos[child].skipConsensus()) {
        if (showAlignments()) {
            #pragma omp critical  // Protect stderr output
            {
              fprintf(stderr, "promoteLowQualUniques()-- PROMOTE READ read %d/%d ident %d at coordinates %d-%d to be used in consensus because it's unique and bases it covers below %d is %d\n", child, _numReads, _utgpos[child].ident(), _utgpos[child].min(), _utgpos[child].max(), minCov, basesBelow);
            }
        }
        _utgpos[child].skipConsensus(false);
     }
  }
  delete[] cov;
}



void unitigConsensus::switchToUncompressedCoordinates(void) {
  // update coordinates of the tig when needed (when it was assembled in compressed space, in normal space this will be skiped)
  // we do this by tracking the read reaching furthest to the right and keeping its offset + homopolymer coordinate translation
  // the read that overlaps it is then updated to start at that reads uncompressed offset + uncompressed bases based on the overlapping coordinate positions
  //

  // just rely on first non-loq qual read
  uint32 checkID = 0;
  while (checkID < _numReads && (_utgpos[checkID].isLowQuality() || _utgpos[checkID].skipConsensus()))
     checkID++;
  if (checkID >= _numReads) { // if we have hit numreads and we had no reads where we didn't skip consensus/low qual, try again ignoring any skip reads, may be a gap-fill read only tig
     checkID = 0;
     while (checkID < _numReads && (_utgpos[checkID].skipConsensus()))
        checkID++;
  }
  assert(checkID < _numReads); // if we have hit numreads and we had no reads where we didn't skip consensus/low qual, what is this tig?
  if (showPlacement()) fprintf(stderr, "switchToUncompressedCoordinates()-- INFO: Checking read %d which had coordinates %d - %d aka %d bp vs raw %d or %f vs 1.2 threshold\n", _utgpos[checkID].ident(), _utgpos[checkID].min(), _utgpos[checkID].max(), (_utgpos[checkID].max()-_utgpos[checkID].min()), getSequence(checkID)->length(), (double)getSequence(checkID)->length() / (_utgpos[checkID].max()-_utgpos[checkID].min()));
  if ((double)getSequence(checkID)->length() / (_utgpos[checkID].max()-_utgpos[checkID].min()) <= 1.2)
    return;

  uint32 compressedOffset   = 0;
  uint32 uncompressedOffset = 0;
  uint32 currentEnd         = _utgpos[0].max();
  uint32 layoutLen          = 0;
  uint32 nlen               = 0;
  uint32* ntoc              = NULL;

  for (uint32 child = 0; child < _numReads; child++) {

    if (compressedOffset > _utgpos[child].min())
      fprintf(stderr, "switchToUncompressedCoordinates()-- ERROR1 in gap in positioning, last read ends at %d; next read starts at %d\n",
              compressedOffset, _utgpos[child].min());
    assert(_utgpos[child].min() >= compressedOffset);

    uint32 readCompressedPosition =  _utgpos[child].min() - compressedOffset;
    uint32 compressedEnd          =  _utgpos[child].max();
    uint32 compressedStart        =  _utgpos[child].min();
    bool   isHQ                   = !_utgpos[child].skipConsensus();

    // find the start position in normal read based on position in compressed read
    uint32 i = 0;
    while (i < nlen && (ntoc[i] < readCompressedPosition))
      i++;

    _utghpc[child].setMinMax(_utgpos[child].min(), _utgpos[child].max());
    _utgpos[child].setMinMax(i+uncompressedOffset, i+uncompressedOffset+getSequence(child)->length());

    // update best end if needed
    if ((ntoc == NULL || compressedEnd > currentEnd) && isHQ) {
      nlen  = getSequence(child)->length();
      delete[] ntoc;
      ntoc  = new uint32 [ nlen + 1 ];
      uint32 clen  = homopolyCompress(getSequence(child)->getBases(), nlen, NULL, ntoc);

      currentEnd         = compressedEnd;
      compressedOffset   = compressedStart;
      uncompressedOffset = _utgpos[child].min();
    }

    if (_utgpos[child].max() > layoutLen)
      layoutLen = _utgpos[child].max();
  }
  delete[] ntoc;
  _tig->_layoutLen = layoutLen;
}


// we assume some reads have good positions (recorded in cnspos) but most do
// not use those reads with good positions as anchors to re-compute everyone
// else's positions too this is very similar to the strategy in the above
// function to convert between compressed and uncompressed coordinates and it
// would be nice to unify the logic
//
void unitigConsensus::updateReadPositions(void) {
  uint32 newOffset = 0;
  uint32 oldOffset = 0;

  for (uint32 child = 0; child < _numReads; child++) {
    if (oldOffset > _utgpos[child].min())
      fprintf(stderr, "switchToUncompressedCoordinates()-- ERROR1 in gap in positioning, last read ends at %d; next read starts at %d\n",
              oldOffset, _utgpos[child].min());
    assert(_utgpos[child].min() >= oldOffset);

    // update best end if needed
    if (_cnspos[child].max() != 0) {
      newOffset = _cnspos[child].min();
      oldOffset = _utgpos[child].min();

      //if (showAlgorithm())
      // fprintf(stderr, "updatePostTemplate()-- Updating guide read to be %d which ends at %d. Now my guide originally was at %d now is at %d\n", _utgpos[child].ident(), _utgpos[child].max(), oldOffset, newOffset);
    }

    uint32 readPosition = _utgpos[child].min() - oldOffset;

    //if (showAlgorithm())
    //  fprintf(stderr, "updatePostTemplate()-- I'm trying to find start of child %d (dist from guide read is %d) currently from %d-%d\n", _utgpos[child].ident(), readPosition, _utgpos[child].min(), _utgpos[child].max());

    _utgpos[child].setMinMax(readPosition+newOffset, readPosition+newOffset+getSequence(child)->length());

    //if (showAlgorithm())
    //  fprintf(stderr, "updatePostTemplate() --Updated read %d which has length %d to be from %d - %d\n", _utgpos[child].ident(), getSequence(child)->length(), _utgpos[child].min(), _utgpos[child].max());
  }
}

// this function checks if we have a template gap if we only use high-quality data
// if we do, allow low-quality reads to be selected as the next to be integrated into the template
// otherwise, only use high-quality reads to avoid introducing errors into the template
bool unitigConsensus::needLowQualReads(uint32 rid, uint32 ePos, std::set<uint32>& badToAdd) {
    uint32 currMax       = ePos;
    uint32 lowQualMax    = ePos;
    uint32 lowQualReadID = 0;
    uint32 nonCNSQualMax = ePos;
    uint32 nonCNSReadID  = 0;

    // back up a bit from the last added read and look for a suitable extension
    for (uint32 ii=rid+1; ii < _numReads; ii++) {
      if (_cnspos[ii].max() != 0 || _utgpos[ii].max() < ePos || badToAdd.count(ii) != 0)    //  If read is contained in the template
        continue;                                                                           //  skip the read.

      // if we ran out of large overlaps, check if this read would still be OK to use given our selected extension, if not, allow low quality reads
      if (_utgpos[ii].min() + _minOverlap >= ePos) {
          if (ePos + _minOverlap >= currMax) {
             // we have an in consensus low quality read so we will promote it and tell template to use it
             if (ePos + _minOverlap < lowQualMax) {
                if (showPlacement()) fprintf(stderr, "generateTemplateStitch()-- read %d/%d ident %d will allow low qual read %d aka %d at template %d\n", ii, _numReads, _utgpos[ii].ident(), lowQualReadID, _utgpos[lowQualReadID].ident(), ePos);
             // we have a read but it's not currently in connsesus, promote it
             } else if (ePos + _minOverlap < nonCNSQualMax) {
                if (showPlacement()) fprintf(stderr, "generateTemplateStitch()-- read %d/%d ident %d will allow non-cns read %d aka %d at template %d\n", ii, _numReads, _utgpos[ii].ident(), nonCNSReadID, _utgpos[nonCNSReadID].ident(), ePos);
                _utgpos[nonCNSReadID].skipConsensus(false);
             } else {
                 if (showPlacement()) fprintf(stderr, "generateTemplateStitch()-- read %d/%d ident %d no read to close a gap in consenus at template %d\n", ii, _numReads, _utgpos[ii].ident(), ePos);
             }
             return true;
         }
         return false;
      } else {
         if (currMax < _utgpos[ii].max() && !_utgpos[ii].isLowQuality()) {
            currMax = _utgpos[ii].max();
         }
         if (lowQualMax < _utgpos[ii].max() && !_utgpos[ii].skipConsensus()) {
            lowQualMax = _utgpos[ii].max();
            lowQualReadID = ii;
         }
         if (nonCNSQualMax < _utgpos[ii].max() && _utgpos[ii].skipConsensus()) {
            nonCNSQualMax = _utgpos[ii].max();
            nonCNSReadID = ii;
         }
      }
   }

   return false;
}

bool
unitigConsensus::templateIsLowQual(uint32 rid) {
   return templateIsLowQual(_utgpos[rid].min(), _utgpos[rid].max());
}

bool
unitigConsensus::templateIsLowQual(uint32 min, uint32 max) {
   for (uint32 i = 0; i < _templateLowQual.size(); i++) {
      if (_templateLowQual.bgn(i) >= max)
         break;
      if (_templateLowQual.end(i) > min)
         return true;
   }
   return false;
}

char *
unitigConsensus::generateTemplateStitch(void) {
  uint32       minOlap  = _minOverlap;
  int32        padding = MAX_PADDING;
  double       savedErrorRate = (_errorRate / ALIGN_FACTOR);
  bool         allowLargerShift  = false;
  bool         allowLowQV        = false;
  int32        lastAddedBP       = 0; // track how much we are adding to the template each time


  //  Find the first non-omitted read, copy that to the template.
retryTig:
  uint32       rid      = 0;

  while ((rid < _numReads) && (_utgpos[rid].skipConsensus() == true))
    rid++;

  abSequence  *seq      = getSequence(rid);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();

  if (showAlgorithm())
    fprintf(stderr, "generateTemplate on tig %u begins at %f seconds.\n", _tig->tigID(), getProcessTime());

  if (showPlacement()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "generateTemplateStitch()-- COPY READ read #%d %d (len=%d to %d-%d)\n",
            rid, _utgpos[rid].ident(), readLen, _utgpos[rid].min(), _utgpos[rid].max());
  }

  uint32       tigmax = std::max(readLen, AS_MAX_READLEN);  //  Must be at least AS_MAX_READLEN, else resizeArray() could fail
  uint32       tiglen = 0;
  char        *tigseq = nullptr;

  allocateArray(tigseq, tigmax, _raAct::clearNew);

  for (uint32 ii=0; ii<readLen; ii++)
    tigseq[tiglen++] = fragment[ii];

  tigseq[tiglen] = 0;

  uint32       ePos = _utgpos[rid].max();   //  Expected end of template, from bogart supplied positions.


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

  while (rid < _numReads) {
    uint32 nr             = 0;   //  Next read
    uint32 nm             = 0;   //  Next read maximum position

    uint32 lastStart      = rid; // Track where we start looking for the next read from, we'll come back and look again if first choice fails
    std::set<uint32> badToAdd;   // Track the list of bad reads, we'll add these at the end
    uint32 firstCandidate = rid; // Track the first read we can use so we know when we have to give up

    _errorRate = savedErrorRate;
    allowLargerShift = false;
    padding = MAX_PADDING;

  retryCandidate:
    //  Reset the candidate back to the beginning and scan forward again
    //  if this is the second or later time we hit this, then we will not pick the failed read and will pick prior good
    rid = lastStart;
    nm  = 0;
    allowLowQV = needLowQualReads(rid, ePos, badToAdd);

    //  Pick the next read as the one with the longest extension from all with some minimum overlap
    //  to the template

    if (showPlacement())
      fprintf(stderr, "\n");

    uint32 ii = rid+1;
    for (ii=rid+1; ii < _numReads; ii++) {
      if ((_utgpos[ii].skipConsensus() == true) ||   //  If  told to not use this read, or
          _cnspos[ii].max() != 0                ||   //  read is already used, or
          (_utgpos[ii].max() <= ePos))               //  read is contained in the template
        continue;                                    //  skip the read.

      //  If a bigger end position, save the overlap.  One quirk: if we've already saved an overlap, and this
      //  overlap is thin or a low-quality read, don't save the overlap.

      bool   thick = (_utgpos[ii].min() + minOlap < ePos);
      bool   first = (nm == 0);
      bool   lowQual = (!allowLowQV && _utgpos[ii].isLowQuality());
      bool   save  = false;

      if ((nm < _utgpos[ii].max() && !lowQual) && (thick || (first && badToAdd.size() == 0)) && badToAdd.count(ii) == 0) {
        save = true;
        nr   = ii;
        nm   = _utgpos[ii].max();
        if ((first && firstCandidate == 0) || ii < firstCandidate)
          firstCandidate = ii;
      }

      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- read #%d/%d ident %d position %d-%d templateEnd %d%s%s%s%s%s\n",
                ii, _numReads, _utgpos[ii].ident(), _utgpos[ii].min(), _utgpos[ii].max(),
                ePos,
                (save  == true)             ? " SAVE"  : "",
                (thick == false)            ? " THIN"  : "",
                (first == true)             ? " FIRST" : "",
                (_utgpos[ii].isLowQuality() ? " LOW" : ""),
                (badToAdd.count(ii) > 0)    ? " BAD" : "");


      //  If this read has an overlap smaller than we want, stop searching.
      if (thick == false) {
        // when we ended on a gap and this is the first time we looked at this position, not just a small overlap, disallow this read and try again to select one that won't leave a gap
        if (badToAdd.size() == 0 && nr != 0 && allowLowQV && ((first && _utgpos[nr].min() > ePos+minOlap) || _utgpos[nr].max() < _utgpos[ii].min()+minOlap)) {
           badToAdd.insert(nr);
           if (showPlacement()) fprintf(stderr, "ERROR: ended up with a gap at %d aka %d from %d-%d and last added was %d going to try to go backwards to %d and allow earlier longer reads to be selected\n", ii, _utgpos[ii].ident(), _utgpos[ii].min(), _utgpos[ii].max(), rid, 0);
           lastStart = 0;
           goto retryCandidate;
        }
        break;
      }
    }

    if (nr == 0) {
      if (ii < _numReads) { // we failed to get to the end but didn't have a read to align? restart from scratch with higher error rate
	     if (_errorRateMax < savedErrorRate) assert (ii >= _numReads); // we can't increase the error rate, really give up
		 fprintf(stderr, "generateTemplateStitch()-- Warning: failed to build tig %d template at error rate %f, trying %f\n",  _tig->tigID(), savedErrorRate, (savedErrorRate*1.25));
		 delete[] tigseq;
		 savedErrorRate *= 1.25;
		 goto retryTig;
	  }
      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- NO MORE READS TO ALIGN\n");
      break;
    }

    assert(nr != 0);

    if ((rid != 0 && _utgpos[rid].isLowQuality()) || _utgpos[nr].isLowQuality()) {
      if (showPlacement()) fprintf(stderr, "generateTemplateStitch()-- Increasing threshold because either %d (%d) or %d (%d) is below requested ID %d\n", rid, _utgpos[rid].ident(), nr, _utgpos[nr].ident(), _errorRateMaxID);
      _errorRate = _errorRateMax / ALIGN_FACTOR;
       padding = MAX_PADDING;
    }
    else {
       _errorRate = savedErrorRate;
       padding = MAX_PADDING;
    }
    int32 olapDiff = (rid == 0 || _utghpc[rid].max() == _utgpos[rid].max()) ? 0 : abs((_utghpc[rid].max()-_utghpc[nr].min())*ALIGN_FACTOR - (_utgpos[rid].max() - _utgpos[nr].min()));
    if (showPlacement()) fprintf(stderr, "generateTemplateStitch()-- Trying to align read %d to template and hpc difference is %d last added read iid %d and new is %d and their coordinates in hpc were %d-%d and %d-%d\n", _utgpos[nr].ident(), olapDiff, rid, nr, _utghpc[rid].min(), _utghpc[rid].max(), _utghpc[nr].min(), _utghpc[nr].max());


    rid      = nr;        //  We'll place read 'nr' in the template.

    seq      = getSequence(rid);
    fragment = seq->getBases();
    readLen  = seq->length();

    int32  readBgn;
    int32  readEnd;

    EdlibAlignResult result;
    bool             aligned       = false;

    double           templateSize  = 0.90;
    double           extensionSize = 0.10;
    double           bandErrRate   = _errorRate / ERROR_RATE_FACTOR;

    int32            olapLen       = ePos - _utgpos[nr].min();  //  The expected size of the overlap
    int32            origLen       = (olapLen < 0 ? 0 : olapLen);

    // compare as ints to ensure that <0 overlap sizes are caught
    if ((int32)olapLen < (int32)minOlap) {
      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- WARNING: increasing min overlap from %d to %u for read %u (%d - %d)\n",
                olapLen, std::min(ePos, minOlap), nr, _utgpos[nr].min(), _utgpos[nr].max());
      if (olapLen > 50 || olapLen < 0) {
        olapLen = std::min(ePos, minOlap);
      } else {
        // hack for mikkos consensus to prevent overlaps
        olapLen=10;
        _errorRate=0;
      }
    }

    int32            templateLen      = 0;
    int32            extensionLen     = 0;

  alignAgain:
    templateLen  = (int32)ceil(olapLen * templateSize);    //  Extract 80% of the expected overlap size
    extensionLen = (int32)ceil(olapLen * extensionSize);   //  Extend read by 20% of the expected overlap size

    readBgn = 0;
    readEnd = olapLen + extensionLen;

    if (readEnd > readLen)
      readEnd = readLen;
    // enforce minimum template length
    if (templateLen <= 1)
      templateLen ++;
    if (templateLen > tiglen)
      templateLen = tiglen-1;

    if (showPlacement()) {
      fprintf(stderr, "\n");
      fprintf(stderr, "generateTemplateStitch()-- ALIGN template %d-%d (len=%d) to read #%d %d %d-%d (len=%d actual=%d at %d-%d)  expecting olap of %d templateLen %d\n",
              tiglen - templateLen, tiglen, templateLen,
              nr, _utgpos[nr].ident(), readBgn, readEnd, readEnd - readBgn, readLen,
              _utgpos[nr].min(), _utgpos[nr].max(),
              olapLen, templateLen);
    }

    result = edlibAlign(tigseq + tiglen - templateLen, templateLen,
                        fragment, readEnd - readBgn,
                        edlibNewAlignConfig(olapLen * bandErrRate, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    //  We're expecting the template to align inside the read.
    //
    //                                                        v- always the end
    //    TEMPLATE  --------------------------[---------------]
    //    READ                          [------------------------------]---------
    //                always the start -^
    //
    //  If we don't find an alignment at all, we move the template start point to the right (making
    //  the template smaller) and also move the read end point to the right (making the read
    //  bigger).

    bool   tryAgain = false;

    bool   noResult      = (result.numLocations == 0);
    bool   gotResult     = (result.numLocations  > 0);

    bool   hitTheStart   = (gotResult) && (result.startLocations[0] == 0);

    bool   hitTheEnd     = (gotResult) && (result.endLocations[0] + 1 == readEnd - readBgn);
    bool   moreToExtend  = (readEnd < readLen);

    // try to find an overlap matching our expectation
    // if we're having issues finding an extension and are at the maximum error rate allowed, give some more tolerance
    int32 maxDifference = std::min(padding, (int32)ceil(0.30*olapLen));
    maxDifference = std::max(maxDifference, olapDiff);

    // pad the difference by the difference between our overlap in compressed space (estimated) and current overlap, this is in rare cases where our uncompression caused a shift
    //  Reset if the edit distance is waay more than our error rate allows or it's very short and we haven't topped out on error.  This seems to be a quirk with
    //  edlib when aligning to N's - I got startLocation = endLocation = 0 and editDistance = alignmentLength.
    if ((double)result.editDistance / result.alignmentLength > bandErrRate || (double)result.alignmentLength / olapLen < 0.51 ||
        (result.numLocations > 0 && abs(result.endLocations[0]-olapLen) > maxDifference && bandErrRate < _errorRate && !allowLargerShift)) {
      noResult    = true;
      gotResult   = false;
      hitTheStart = false;
      hitTheEnd   = false;
    }

    //  HOWEVER, if we get a result and it's near perfect, declare success even if we hit the start.
    //  These are simple repeats that will align with any overlap.  The one BPW debugged was 99+% A.

    if ((gotResult == true) &&
        (hitTheStart == true) &&
        ((double)result.editDistance / result.alignmentLength < 0.01)) {
      hitTheStart = false;
    }

    //  NOTE that if we hit the end with the same conditions, we should try again, unless there
    //  isn't anything left.  In that case, we don't extend the template.

    if ((gotResult == true) &&
        (hitTheEnd == true) &&
        (moreToExtend == false) &&
        ((double)result.editDistance / result.alignmentLength < 0.01)) {
      hitTheEnd = false;
    }

    //  Now, report what happened, and maybe try again.

    if ((showPlacement()) && (noResult == true))
      fprintf(stderr, "generateTemplateStitch()-- FAILED to align - no result\n");

    if ((showPlacement()) && (noResult == false))
      fprintf(stderr, "generateTemplateStitch()-- FOUND alignment at %d-%d editDist %d alignLen %d %.f%% expected %d\n",
              result.startLocations[0], result.endLocations[0]+1,
              result.editDistance,
              result.alignmentLength,
              100.0 * result.editDistance / result.alignmentLength, olapLen);

    if ((noResult) || (hitTheStart)) {
      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - %s - decrease template size by 10%%\n",
                (noResult == true) ? "no result" : "hit the start");
      tryAgain = true;
      templateSize -= 0.10;
    }

    if ((noResult) || (hitTheEnd && moreToExtend)) {
      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - %s - increase read size by 10%%\n",
                (noResult == true) ? "no result" : "hit the end");
      tryAgain = true;
      extensionSize += 0.10;
    }

    if (templateSize < 0.01) {
      if (showPlacement())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - no more template to remove!\n");

      if (bandErrRate + _errorRate / ERROR_RATE_FACTOR > _errorRate || _errorRate - bandErrRate < 1e-9) {
        if (showPlacement()) fprintf(stderr, "  Fail!\n");
        tryAgain = false;
        olapLen = origLen;
      }
      else {
        if (showPlacement())
          fprintf(stderr, "generateTemplateStitch()-- FAILED to align at %.2f error rate, increasing to %.2f\n", bandErrRate, bandErrRate + _errorRate/ERROR_RATE_FACTOR);
        tryAgain = true;
        templateSize  = 0.90;
        extensionSize = 0.10;
        bandErrRate  += _errorRate / ERROR_RATE_FACTOR;
      }
    }

    if (tryAgain) {
      edlibFreeAlignResult(result);
      goto alignAgain;
    }

    //  Use the alignment (or the overlap) to figure out what bases in the read
    //  need to be appended to the template.

    if (noResult == false) {
      readBgn = result.startLocations[0];     //  Expected to be zero
      readEnd = result.endLocations[0] + 1;   //  Where we need to start copying the read

      // record updated read coordinates which we will use later as anchors to re-computed everyone else
      _cnspos[nr].setMinMax(tiglen-readEnd, tiglen + readLen - readEnd);

      if (showPlacement()) {
        fprintf(stderr, "generateTemplateStitch()--\n");
        fprintf(stderr, "generateTemplateStitch()-- Aligned template %d-%d to read %u %d-%d; copy read %d-%d to template.\n",
                tiglen - templateLen, tiglen, nr, readBgn, readEnd, readEnd, readLen);
        fprintf(stderr, "generateTemplateStitch()-- New position for read %d is  %d-%d\n", _utgpos[nr].ident(), _cnspos[nr].min(), _cnspos[nr].max());
      }
    } else {
      // try to go back and pick a different read to extend with if that's an option
      // it's not an option if we're already the next read
      // if we're out of reads but didn't try contains, try again anyway, after this we give up and trim the template below
      // note that we check if this read was already tried because when we go back and re-try the contained reads, we may end up not at the firstCandidate (since it was uncontained) but if we run out of new reads to try we need to stop
      // when the read was a low quality read, we may have promoted it to consensus so demote it back, we can always try it again in the larger shift phase, if even that fails then we don't want it
      bool isAlreadyBad = badToAdd.count(nr) != 0;
      badToAdd.insert(nr);
      if (_utgpos[nr].isLowQuality()) _utgpos[nr].skipConsensus(true);
      if ((rid != firstCandidate && !isAlreadyBad) || allowLargerShift == false) {
        allowLargerShift = (allowLargerShift || rid == firstCandidate || isAlreadyBad);  // we hit the end so now we can try to allow larger shifts for the same reads
        // when we've switched to allowLargerShift, we reset our list of bad reads so we can try them again allowing larger shift
        // we also make sure to reset the error rate back to default since we'll be processing a new read
        if (rid == firstCandidate || isAlreadyBad) {
           badToAdd.clear();
           badToAdd.insert(rid);
        }
        _errorRate = savedErrorRate;
        olapLen = origLen;
        if (showPlacement()) {
          fprintf(stderr, "generateTemplateStitch()--\n");
          fprintf(stderr, "generateTemplateStitch()-- FAILED to align read #%d ident %d to last added %d ident %d, will go back and re-try an earlier read until we hit %d and allowing shift is %d\n", rid, _utgpos[nr].ident(), lastStart, _utgpos[lastStart].ident(), firstCandidate, allowLargerShift);
          fprintf(stderr, "generateTemplateStitch()--\n");
        }
        edlibFreeAlignResult(result);
        goto retryCandidate;
      }
      else if (tiglen > minOlap && origLen > (minOlap + TRIM_BP) && lastAddedBP > (minOlap + TRIM_BP) && lastStart > 0) {
        int32 trimbp = std::min((_utgpos[lastStart].max()-_utgpos[lastStart].min()) - minOlap, lastAddedBP - minOlap - TRIM_BP);
        if (trimbp > tiglen) trimbp = tiglen - minOlap;
        assert(trimbp > 0);
        assert(lastAddedBP > 0);
        assert(ePos > trimbp);

        if (showPlacement()) {
          fprintf(stderr, "generateTemplateStitch()--\n");
          fprintf(stderr, "generateTemplateStitch()-- FAILED to align read #%d ident %d to last added %d ident %d, will try to trim template of %d bp by %d bases\n", rid, _utgpos[nr].ident(), lastStart, _utgpos[lastStart].ident(), tiglen, trimbp);
        }

        // trim the template and also reset error rate and the list of bad reads since we're going to be retrying them all
        ePos -= trimbp;
        lastAddedBP = 0;
        tigseq[tiglen-trimbp] = 0;
        tiglen=strlen(tigseq);
        firstCandidate = 0;
        badToAdd.clear();
        _errorRate = savedErrorRate;
        allowLargerShift = false;

        if (showPlacement()) {
          fprintf(stderr, "generateTemplateStitch()-- Trimmed template to %d bp\n", tiglen);
          fprintf(stderr, "generateTemplateStitch()--\n");
        }
        edlibFreeAlignResult(result);
        goto retryCandidate;
      }

      readBgn = 0;
      readEnd = olapLen;

      // avoid pathological case where we try to extend the template by nothing
      if (readEnd-readLen == 0) {
         if (showPlacement()) fprintf(stderr, "Nothing to append going to try again\n");
         lastStart=rid;
         edlibFreeAlignResult(result);
         goto retryCandidate;
      }
      if (showPlacement()) {
        fprintf(stderr, "generateTemplateStitch()--\n");
        fprintf(stderr, "generateTemplateStitch()-- Alignment failed, use original overlap; copy read %d-%d to template.\n",
                readEnd, readLen);
      }
    }

    edlibFreeAlignResult(result);

    resizeArray(tigseq, tiglen, tigmax, tiglen + readLen - readEnd + 1);

    //  Append the read bases to the template.
    //
    //  When the read fails to align, we used to go back in the overlap region and fill in any N
    //  bases, expected to be large blocks of Ns from shredded scaffolds, with bases from the read.
    //  As of March 2019, this isn't really needed, for two reasons: (1) sqStore is trimming N's
    //  from the end of a read, and (2) edlib is allowing matches to Ns.

    for (uint32 ii=readEnd; ii<readLen; ii++)
      tigseq[tiglen++] = fragment[ii];
    tigseq[tiglen] = 0;
    assert(tiglen < tigmax);

    ePos = _utgpos[rid].max();
    lastAddedBP = readLen-readEnd;

    // record if this alignment used a low quality read of we were forced to join a read
    // then we can allow more error tolerance when aligning to this part of the template
    if ((rid != 0 && _utgpos[rid].isLowQuality()) || _utgpos[nr].isLowQuality() || readEnd == 0)
       _templateLowQual.add_position(std::max((uint32)0, tiglen-readLen), tiglen);

    if (showAlgorithm())
      fprintf(stderr, "generateTemplateStitch()-- Append read #%-8d ident %-8d %9u-%-9u to tig %9u-%-9u expectedolap %-9u expected %-9u difference %7.4f%%\n",
              rid, _utgpos[nr].ident(),
              readEnd, readLen,
              tiglen-(readLen-readEnd), tiglen,
              olapLen,
              ePos, 200.0 * ((int32)tiglen - (int32)ePos) / ((int32)tiglen + (int32)ePos));
  }

  //  Report the expected and final size.  We used to guard against long tigs getting chopped, but
  //  the new template construction (early 2019) will simply append reads that fail to extend the
  //  template - so long tigs never should get chopped.

  double  pd = 200.0 * ((int32)tiglen - (int32)ePos) / ((int32)tiglen + (int32)ePos);

  if (showAlgorithm()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "generateTemplateStitch()-- generated template of length %d, expected length %d, %7.4f%% difference.\n",
            tiglen, ePos, pd);
  }

  // restore saved error rate
  _errorRate = savedErrorRate * ALIGN_FACTOR;
  updateReadPositions();

  return(tigseq);
}



bool
alignEdLib(dagAlignment      &aln,
           tgPosition        &utgpos,
           char              *fragment,
           uint32             fragmentLength,
           char              *tigseq,
           uint32             tiglen,
           double             errorRate,
           int32              maxpad,
           bool               verbose) {

  EdlibAlignResult align;

  int32   padding        = std::min(maxpad, (int32)ceil(fragmentLength * 0.05));
  double  bandErrRate    = errorRate / ERROR_RATE_FACTOR;
  bool    aligned        = false;
  double  alignedErrRate = 0.0;

  //  Decide on where to align this read.

  int32  tigbgn = std::max((int32)0,      (int32)floor(utgpos.min() - padding));
  int32  tigend = std::min((int32)tiglen, (int32)floor(utgpos.max() + padding));

  if (verbose)
    fprintf(stderr, "alignEdLib()-- align read %7u padding %7u maxpad %7u eRate %.4f at %9d-%-9d\n", utgpos.ident(), padding, maxpad, bandErrRate, tigbgn, tigend);

  if (tigend <= tigbgn) {
    fprintf(stderr, "alignEdLib()-- WARNING: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
            tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
    // try to align it to full
    tigbgn = 0;
    tigend = utgpos.max();
    tigend = (tigend < 0)      ? tiglen : tigend;   //  if utgpos.max is invalid,
    tigend = (tigend > tiglen) ? tiglen : tigend;   //  use tiglen.

    fprintf(stderr, "alignEdLib()-- WARNING: updated tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
            tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
  }
  assert(tigend > tigbgn);

  //  Align!  If there is an alignment, compute error rate and declare success if acceptable.

  align = edlibAlign(fragment, fragmentLength,
                     tigseq + tigbgn, tigend - tigbgn,
                     edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

  if (align.alignmentLength > 0) {
    alignedErrRate = (double)align.editDistance / align.alignmentLength;
    aligned        = (alignedErrRate <= bandErrRate) && align.alignmentLength >= fragmentLength;
    if (verbose)
      fprintf(stderr, "alignEdLib()-- - ALIGNED read %7u %.4f at %9d-%-9d\n", utgpos.ident(), alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
  } else {
    if (verbose)
      fprintf(stderr, "\n");
  }

  for (uint32 ii=1 /*the 0ths was above*/; ((ii < MAX_RETRIES) && (aligned == false)); ii++) {
    tigbgn = std::max((int32)0,      tigbgn - 5 * padding);
    tigend = std::min((int32)tiglen, tigend + 5 * padding);
    // last attempt make a very wide band
    if ((ii+1) % NUM_BANDS == 0) {
      tigbgn = std::max((int32)0,      tigbgn - MAX_PADDING_MULT * padding);
      tigend = std::min((int32)tiglen, tigend + MAX_PADDING_MULT * padding);
    }

    // let the band increase without increasing error rate for a while, if we give up increase error rate
    // and reset bnad
    if (ii != 0 && ii % NUM_BANDS == 0) {
      bandErrRate += errorRate / ERROR_RATE_FACTOR;
      tigbgn = std::max((int32)0,      (int32)floor(utgpos.min() - padding));
      tigend = std::min((int32)tiglen, (int32)floor(utgpos.max() + padding));
    }

    if (tigend <= tigbgn) {
      fprintf(stderr, "alignEdLib()-- WARNING: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
              tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
      // try to align it to full
      tigbgn = 0;
      tigend = utgpos.max();
      tigend = (tigend < 0)      ? tiglen : tigend;   //  if utgpos.max is invalid,
      tigend = (tigend > tiglen) ? tiglen : tigend;   //  use tiglen.

      fprintf(stderr, "alignEdLib()-- WARNING: updated tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
              tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
    }
    assert(tigend > tigbgn);

    edlibFreeAlignResult(align);

    if (verbose)
      fprintf(stderr, "alignEdLib()--                    read %7u eRate %.4f at %9d-%-9d\n", utgpos.ident(), bandErrRate, tigbgn, tigend);

    align = edlibAlign(fragment, strlen(fragment),
                       tigseq + tigbgn, tigend - tigbgn,
                       edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (align.alignmentLength > 0) {
      alignedErrRate = (double)align.editDistance / align.alignmentLength;
      aligned        = (alignedErrRate <= bandErrRate && align.alignmentLength >= fragmentLength);
      if (verbose)
        fprintf(stderr, "alignEdLib()--  - ALIGNED read %7u %.4f at %9d-%-9d\n", utgpos.ident(), alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
    } else {
      if (verbose)
        fprintf(stderr, "\n");
    }
  }

  if (aligned == false) {
    edlibFreeAlignResult(align);
    return(false);
  }

  aln.start  = tigbgn + align.startLocations[0] + 1;   //  AlnGraphBoost expects 1-based positions.
  aln.end    = tigbgn + align.endLocations[0] + 1;     //  EdLib returns 0-based positions.

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
                          tgtaln,                        //  output tgt alignment string
                          qryaln);                       //  output qry alignment string

  //  Populate the output.  AlnGraphBoost does not handle mismatch alignments, at all, so convert
  //  them to a pair of indel.

  uint32 nMatch = 0;

  for (uint32 ii=0; ii<align.alignmentLength; ii++)   //  Edlib guarantees aln[alignmentLength] == 0.
    if ((tgtaln[ii] != '-') &&
        (qryaln[ii] != '-') &&
        (tgtaln[ii] != qryaln[ii]))
      nMatch++;

  aln.qstr   = new char [align.alignmentLength + nMatch + 1];
  aln.tstr   = new char [align.alignmentLength + nMatch + 1];
  memset(aln.tstr, 0, sizeof(char)*(align.alignmentLength + nMatch + 1));
  memset(aln.qstr, 0, sizeof(char)*(align.alignmentLength + nMatch + 1));

  for (uint32 ii=0, jj=0; ii<align.alignmentLength; ii++) {
    assert(jj < align.alignmentLength + nMatch + 1);

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

    // push gaps to the right, but not past the end
    for (uint32 ii=0; ii < aln.length-1; ii++) {
        // pushing target gaps
        if (aln.tstr[ii] == '-') {
            uint32 jj = ii;
            while (jj < aln.length - 1) {
                char c = aln.tstr[++jj];
                if (c != '-') {
                    assert(jj < aln.length);
                    if (c == aln.qstr[ii]) {
                        aln.tstr[ii] = c;
                        aln.tstr[jj] = '-';
                    }
                    break;
                }
            }
        }

        // pushing query gaps
        if (aln.qstr[ii] == '-') {
            uint32 jj = ii;
            while (jj < aln.length -1) {
                char c = aln.qstr[++jj];
                if (c != '-') {
                    assert(jj < aln.length);
                    if (c == aln.tstr[ii]) {
                        aln.qstr[ii] = c;
                        aln.qstr[jj] = '-';
                    }
                    break;
                }
            }
        }
    }

    // generate the final, normalized alignment strings
    uint32 write_ii = 0;
    for (uint32 ii = 0; ii < aln.length; ii++) {
        if (aln.tstr[ii] != '-' || aln.qstr[ii] != '-') {
            aln.tstr[write_ii] = aln.tstr[ii];
            aln.qstr[write_ii] = aln.qstr[ii];
            write_ii++;
        }
    }
    aln.length = write_ii;


  aln.qstr[aln.length] = 0;
  aln.tstr[aln.length] = 0;

  delete [] tgtaln;
  delete [] qryaln;

  edlibFreeAlignResult(align);
  if (verbose)
    fprintf(stderr, "alignEdLib()--  - ALIGNED read %7u strings:\n%s\n%s\n", utgpos.ident(), aln.tstr, aln.qstr);

  if (aln.end > tiglen)
    fprintf(stderr, "ERROR:  alignment from %d to %d, but tiglen is only %d\n", aln.start, aln.end, tiglen);
  assert(aln.end <= tiglen);

  return(true);
}



bool
unitigConsensus::generatePBDAG(char aligner_, uint32 numIterations_, u32toRead &reads_) {

  //  Build a quick consensus to align to.

  if (showAlgorithm())
    fprintf(stderr, "generatePBDAG on tig %u begins at %f seconds.\n", _tig->tigID(), getProcessTime());

  char   *tigseq = generateTemplateStitch();
  uint32  tiglen = strlen(tigseq);

  if (showAlgorithm()) {
    fprintf(stderr, "Generated template of length %d\n", tiglen);
    fprintf(stderr, "Aligning reads at %f seconds.\n", getProcessTime());
  }
  promoteLowQualUniques(tiglen, MIN_COV, MIN_COV_SIZE);

  //  Compute alignments of each sequence in parallel
  std::string cns;
  cns.reserve(tiglen*1.1);
  for (uint32 iteration = 1; iteration <= numIterations_; iteration++) {
    dagAlignment *aligns = new dagAlignment [_numReads];
    uint32        pass = 0;
    uint32        fail = 0;
    uint32        skip = 0;

    #pragma omp parallel for schedule(dynamic)
    for (uint32 ii=0; ii<_numReads; ii++) {
      abSequence  *seq      = getSequence(ii);
      bool         aligned  = false;

      assert(aligner_ == 'E');  //  Maybe later we'll have more than one aligner again.

      if (_utgpos[ii].skipConsensus() == false) {
         aligned = alignEdLib(aligns[ii],
                             _utgpos[ii],
                             seq->getBases(), seq->length(),
                            tigseq, tiglen,
                            (_utgpos[ii].isLowQuality() == true || templateIsLowQual(ii)) ? (_errorRateMax) : (_errorRate),   //  Allow higher error for low-quality reads or low quality template
                            (_utgpos[ii].isLowQuality() == true || templateIsLowQual(ii)) ? MAX_PADDING*5 : MAX_PADDING,      //  Pad low-quality reads by more
                            showAlignments());

         if (aligned == false) {
           if (showAlgorithm())
             fprintf(stderr, "generatePBDAG()--    read %7u FAILED\n", _utgpos[ii].ident());

           fail++;

           continue;
         }

        pass++;
      } else
        skip++;
    }

    if (showAlgorithm())
      fprintf(stderr, "generatePBDAG()--    read alignment: %d failed, %d skipped, %d passed.\n", fail, skip, pass);

    //  Construct the graph from the alignments.  This is not thread safe.

    if (showAlgorithm())
      fprintf(stderr, "Constructing graph at %f seconds.\n", getProcessTime());

    AlnGraphBoost ag(std::string(tigseq, tiglen));

    for (uint32 ii=0; ii<_numReads; ii++) {
      _cnspos[ii].setMinMax(aligns[ii].start, aligns[ii].end);

      if ((aligns[ii].start == 0) &&
          (aligns[ii].end   == 0))
        continue;

      if (_utgpos[ii].skipConsensus() == false)
        ag.addAln(aligns[ii]);

      aligns[ii].clear();
    }

    delete [] aligns;

    if (showAlgorithm())
      fprintf(stderr, "Merging graph iteration %d at %f seconds.\n", iteration, getProcessTime());

    //  Merge the nodes and call consensus
    ag.mergeNodes();

    if (showAlgorithm())
      fprintf(stderr, "Calling consensus iteration %d at %f seconds.\n", iteration, getProcessTime());

    //FIXME why do we have 0weight nodes (template seq w/o support even from the read that generated them)?

    assert(_templateToCNS == nullptr);

    _templateToCNS  = new uint32 [tiglen + 1];
    _templateLength = tiglen;

    delete [] tigseq;
    tigseq=nullptr;
    cns = ag.consensusNoSplit(((_tig->_suggestNoTrim == 0 || templateIsLowQual(0, std::min((uint32)MIN_COV_SIZE*10, tiglen)) || templateIsLowQual(std::max(0, (int32)tiglen-MIN_COV_SIZE*10), tiglen)) && iteration == numIterations_ ? _minCoverage : 0), _templateToCNS, _templateLength);

    if (iteration == numIterations_ || cns.size() == 0)
      break;

    tiglen = (uint32) cns.size();
    allocateArray(tigseq, tiglen+1, _raAct::clearNew);
    memcpy(tigseq, cns.data(), cns.size());
    tigseq[tiglen] = 0;

    // update positions for second round
    for (uint32 jj=0; jj<_numReads; jj++)
      adjustPosition(_utgpos[jj], _cnspos[jj], _utgpos[jj], false);
    delete[] _templateToCNS;
    _templateToCNS=nullptr;
  }

  //  Save consensus

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, (uint32) cns.length() + 1, _raAct::doNothing);

  std::string::size_type len = 0;

  for (len=0; len<cns.size(); len++) {
    _tig->_bases[len] = cns[len];
    _tig->_quals[len] = CNS_MIN_QV;
  }

  //  Terminate the string.

  _tig->_bases[len] = 0;
  _tig->_quals[len] = 0;
  _tig->_basesLen   = len;
  _tig->_layoutLen  = len;

  assert(len < _tig->_basesMax);

  if (showAlgorithm())
    fprintf(stderr, "generatePBDAG on tig %u finishes at %f seconds.\n", _tig->tigID(), getProcessTime());

  return(true);
}



bool
unitigConsensus::generateQuick(u32toRead &reads_) {

  //  Quick is just the template sequence, so one and done!

  char   *tigseq = generateTemplateStitch();
  uint32  tiglen = strlen(tigseq);

  //  Save consensus

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, tiglen + 1, _raAct::doNothing);

  for (uint32 ii=0; ii<tiglen; ii++) {
    _tig->_bases[ii] = tigseq[ii];
    _tig->_quals[ii] = CNS_MIN_QV;
  }

  //  Set positions of all the reads.  We don't know anything and default to the incoming positions.

  for (uint32 ii=0; ii<_numReads; ii++)
    _cnspos[ii] = _utgpos[ii];

  //  Terminate the string.

  _tig->_bases[tiglen] = 0;
  _tig->_quals[tiglen] = 0;
  _tig->_basesLen      = tiglen;
  _tig->_layoutLen     = tiglen;

  delete [] tigseq;

  return(true);
}



bool
unitigConsensus::generateSingleton(u32toRead &reads_) {

  assert(_numReadsUsable == 1);

  //  Find the first usable read.

  uint32 fid=0;
  while ((fid < _numReads) && (_utgpos[fid].skipConsensus() == true))
    fid++;

  assert(fid < _numReads);

  //  Copy the single read to the tig sequence.

  abSequence  *seq      = getSequence(fid);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, readLen + 1, _raAct::doNothing);

  for (uint32 ii=0; ii<readLen; ii++) {
    _tig->_bases[ii] = fragment[ii];
    _tig->_quals[ii] = CNS_MIN_QV;
  }

  //  Set positions of all the reads.

  _cnspos[0].setMinMax(0, readLen);

  //  Terminate the string.

  _tig->_bases[readLen] = 0;
  _tig->_quals[readLen] = 0;
  _tig->_basesLen       = readLen;
  _tig->_layoutLen      = readLen;

  return(true);
}



void
edlibToDelta(tgTig             *tig,
             tgPosition        *child,
             uint32             readlen,
             EdlibAlignResult  &align) {
  const unsigned char *alignment       = align.alignment;
  const int32          alignmentLength = align.alignmentLength;

        int32 tgtPos = align.startLocations[0];
  const int32 tgtEnd = align.endLocations[0] + 1;
        int32 qryPos = 0;
  const int32 qryEnd = readlen;

  stuffedBits *bits = tig->_childDeltaBits;

  int32  nMatch  = 0;   //  Counting the number of match/mismatch in a delta block.
  int32  nBlocks = 0;   //  Number of blocks in the delta encoding.

#pragma omp critical (tgTigLoadAlign)
  {
    child->_deltaOffset = bits->getPosition();

    for (int32 a=0; a<alignmentLength; a++) {
      assert(tgtPos <= tgtEnd);   //  We don't run off the end of the tig
      assert(qryPos <= qryEnd);   //  We don't run off the end of the read

      if ((alignment[a] == EDLIB_EDOP_MATCH) ||      //  Match or mismatch.
          (alignment[a] == EDLIB_EDOP_MISMATCH)) {   //  Just count how many we see.
        nMatch++;      qryPos++;   tgtPos++;
      }
      else if (alignment[a] == EDLIB_EDOP_INSERT) {  //  Insertion in target.
        bits->setEliasDelta(nMatch + 1);   bits->setBit(0);
        nMatch = 0;   nBlocks++;   qryPos++;
      }
      else if (alignment[a] == EDLIB_EDOP_DELETE) {  //  Insertion in query.
        bits->setEliasDelta(nMatch + 1);   bits->setBit(1);
        nMatch = 0;   nBlocks++;   tgtPos++;
      }
    }
  }

  child->_deltaLen = nBlocks;
}



//  Find the position of the read on the final consensus.
//  If the read was placed by consensus, cnspos is valid, and we use that
//  to decide position.  If not, cnspos is not valid and we fall back
//  to the slightly incorrect utgpos.
//
void
unitigConsensus::adjustPosition(tgPosition   utgpos,
                                tgPosition   cnspos,
                                tgPosition  &adjusted,
                                bool         isSingle) {
  tgPosition  original = ((cnspos.min() == 0) && (cnspos.max() == 0)) ? utgpos : cnspos;

  //  Adjust the 'original' position so it is between 0 and templateLength.

  int32 omin = std::min(std::max(0, original.min()), (int32)_templateLength-1);
  int32 omax = std::min(            original.max() , (int32)_templateLength-1);

  //  But if omin == omax
  assert(omin < _templateLength);

  //  Map that position to the final consensus.

  int32 amin = _templateToCNS[omin];
  int32 amax = _templateToCNS[omax];

  //  If that template position did not make it to final consensus,
  //  search left/right for the first one that did make it.

  while ((amin == UINT32_MAX) && (omin > 0))
    amin = _templateToCNS[--omin];

  while ((amax == UINT32_MAX) && (omax < _templateLength))
    amax = _templateToCNS[++omax];

  //  If still invalid, or the tig is a singleton, reset to the whole tig.

  if ((isSingle == true) || (amin == UINT32_MAX))  amin = 0;
  if ((isSingle == true) || (amax == UINT32_MAX))  amax = _templateToCNS[_templateLength-1];

  adjusted.setMinMax(amin, amax);
}



//  Align each read to the region of the consensus sequence the read claims
//  to be from, extended by 5% of the read length on either end.  If it fails
//  to align full length, make the extensions larger and/or error rate more
//  tolerant until it aligns.
//
//  Reads that fail to align have their cnspos set to 0.
//
void
unitigConsensus::findCoordinates(char algorithm_, u32toRead  &reads_) {
  if (algorithm_ != 'P')      return;   //  -norealign or -quick
  if (_tig->length() == 0)    return;   //  Failed consensus.
  if (_templateLength == 0)   return;   //  Singleton.

  if (showAlgorithm())
    fprintf(stderr, "findCoordinates on tig %u starts at %f seconds.\n", _tig->tigID(), getProcessTime());

  _tig->_childDeltaBitsLen = 1;
  _tig->_childDeltaBits    = new stuffedBits();

  _tig->_childCIGAR        = new char * [_numReads];

  for (uint32 ii=0; ii<_numReads; ii++)
    _tig->_childCIGAR[ii] = nullptr;

  uint32     pass        = 0;
  uint32     fail        = 0;

  char      *tigseq      = _tig->bases();
  int32      tiglen      = strlen(tigseq);

  bool       isSingleton = (_numReadsUsable == 1);

#pragma omp parallel for schedule(dynamic)
  for (uint32 ii=0; ii<_numReads; ii++) {
    char const *readseq = getSequence(ii)->getBases();
    uint32      readlen = getSequence(ii)->length();

    adjustPosition(_utgpos[ii], _cnspos[ii], _adjpos[ii], isSingleton);

    //  Decide on where to align this read and the expected quality.

    int32  padding  = std::min(MAX_PADDING, (int32)ceil(readlen * 0.05));
    double erate    = _utgpos[ii].isLowQuality() ? _errorRateMax : _errorRate;

    int32  bandpad[9]  = {  1 * padding,  5 * padding, MAX_PADDING_MULT * padding,
                            1 * padding,  5 * padding, MAX_PADDING_MULT * padding,
                            1 * padding,  5 * padding, MAX_PADDING_MULT * padding };

    double bandqual[9] = { 0.25 * erate,  0.20 * erate, 0.15 * erate,
                           0.50 * erate,  0.45 * erate, 0.40 * erate,
                           0.75 * erate,  0.90 * erate, 1.00 * erate };

    EdlibAlignResult align;

    int32   alignbgn  = 0;       //  How evil is it to have variables
    int32   alignend  = 0;       //       alignend
    double  alignqual = 1.0;     //  and  |||||| |    (they even align!)
    bool    aligned   = false;   //       aligne d

    //  Make several attempts to align.
    //
    //  The general idea is to make several attempts at high-, medium- and
    //  low-quality alignments (bandqual), each attempt lowering the
    //  precision of where the read is expected to align (bandpad).
    //
    //  A further refinement is to require slightly higher identity as the
    //  location precision is reduced - most alignments succeed on the first
    //  attempt.

    for (uint32 attempt=0; (attempt < 9) && (aligned == false); attempt++) {
      int32 aspanbgn  = std::max(0, _adjpos[ii].min() - bandpad[attempt]);
      int32 aspanend  = std::min(   _adjpos[ii].max() + bandpad[attempt], tiglen);

      assert(aspanbgn < aspanend);

      if (showPlacement()) {
        if (attempt == 0)
          fprintf(stderr, "ALIGN read #%-6d %8u length %6u - utgpos %8d,%-8d cnspos %8d,%-8d -> consensus %8d,%-8d - attempt %d:",
                  ii, _utgpos[ii].ident(), readlen,
                  _utgpos[ii].min(), _utgpos[ii].max(),
                  _cnspos[ii].min(), _cnspos[ii].max(),
                  aspanbgn,          aspanend, attempt+1);
        else
          fprintf(stderr, "                                                                                              -> consensus %8d,%-8d - attempt %d:",
                  aspanbgn,          aspanend, attempt+1);
      }

      align = edlibAlign(readseq, readlen,
                         _tig->bases() + aspanbgn, aspanend - aspanbgn,
                         edlibNewAlignConfig(bandqual[attempt] * readlen, EDLIB_MODE_HW, EDLIB_TASK_PATH));

      if (align.alignmentLength == 0) {
        alignbgn  = 0;
        alignend  = 0;
        alignqual = 1.0;
        aligned   = false;
      }
      else {
        alignbgn  = aspanbgn + align.startLocations[0];
        alignend  = aspanbgn + align.endLocations[0] + 1;
        alignqual = (double)align.editDistance / align.alignmentLength;
        aligned   = (double)align.editDistance / align.alignmentLength <= bandqual[attempt] && align.alignmentLength >= readlen;
      }

      if (showPlacement()) {
        if (aligned)
          fprintf(stderr, " - SUCCESS %8d,%-8d %6u/%-6u = %6.3f%% < %6.3f%%\n",
                  alignbgn, alignend, align.editDistance, align.alignmentLength,
                  100.0 * alignqual, 100.0 * bandqual[attempt]);
        else if (alignqual < 1.0)
          fprintf(stderr, " - FAILURE %8d,%-8d %6u/%-6u = %6.3f%% > %6.3f%%\n",
                  alignbgn, alignend, align.editDistance, align.alignmentLength,
                  100.0 * alignqual, 100.0 * bandqual[attempt]);
        else
          fprintf(stderr, " - FAILURE\n");
      }

      if (aligned == false)
        edlibFreeAlignResult(align);
    }

    if (aligned == false) {
      _cnspos[ii].        setMinMax(0, 0);                  //  NO final position
      _tig->getChild(ii)->setMinMax(0, 0);                  //  in consensus.

      _tig->_childCIGAR[ii]    = new char [1];              //  NO final alignment.
      _tig->_childCIGAR[ii][0] = 0;

      fail++;
    }
    else {
      _cnspos[ii].        setMinMax(alignbgn, alignend);    //  Final position in consensus.
      _tig->getChild(ii)->setMinMax(alignbgn, alignend);    //  Save in tig.

      _tig->_childCIGAR[ii] = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_EXTENDED);

      edlibToDelta(_tig, _tig->getChild(ii), readlen, align);

      pass++;
    }

    edlibFreeAlignResult(align);
  }    //  Looping over reads.

  if (showAlgorithm())
    fprintf(stderr, "findCoordinates on tig %u finishes at %f seconds. %u pass %u fail.\n", _tig->tigID(), getProcessTime(), pass, fail);
}




void
unitigConsensus::trimCircular(void) {
  char    *bases  = _tig->bases();
  uint32   length = _tig->length();

  if (_tig->length() == 0)              return;
  if (_tig->_suggestCircular == false)  return;
  if (_tig->_circularLength == 0)       return;

  if (showAlgorithm()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Circularizing tig %u with expected overlap %u bases.\n", _tig->tigID(), _tig->_circularLength);
  }

  //  Try alignments of various lengths until we find a proper overlap
  //  between the end and the start at acceptable quality, or until we get
  //  longer than expected.
  //
  //  Stop searching if we found a proper overlap.
  //
  //                   bgnA -v        v- end of the tig
  //            A:  ....---------------
  //            B:           ----------------......
  //   beginning of the tig -^        ^- endB

  bool          found = false;
  parasailLib   align(1, -4, 5, 3);

  for (double scale = 1.05; (found == false) && (scale < 2.10); scale += 0.1) {
    uint32 trylen = std::min(length, (uint32)(_tig->_circularLength * scale));

    if (showPlacement())
      fprintf(stderr, "  Align   A %6u-%6u against B %6u-%6u\n",
              length-trylen, length,
              0, trylen);

    align.alignDovetail(bases, length, length-trylen, length,           //  A:  the end of the tig
                        bases, length, 0,             trylen, false);   //  B:  the start of the tig

    if ((align.errorRate() < 0.05) &&                //  "Found" if the error rate is reasonable, and
        (align.alignmentLength() >= _minOverlap) &&  //  and length is reasonable
        (align.bgnA() > length - trylen) &&          //  A has unaligned stuff at the start, and
        (align.endB() < trylen))                     //  B has unaligned stuff at the end.
      found = true;

    if ((found == false) && (showPlacement()))
      fprintf(stderr, "  Aligned A %6u-%6u against B %6u-%6u at %7.3f%%\n",
              align.bgnA(), align.endA(), align.bgnB(), align.endB(), align.percentIdentity());
  }

  if (found == false)   //  No overlap found, probably not circular.
    return;

  uint32   headBgn = align.bgnB();   //  Should be zero!
  uint32   headEnd = align.endB();   //

  uint32   tailBgn = align.bgnA();   //
  uint32   tailEnd = align.endA();   //  Should be 'length'.

  uint32   keepBgn = 0;              //  Default to trimming the
  uint32   keepEnd = tailBgn;        //  end part off.

  if (showAlgorithm())
    fprintf(stderr, "  Aligned A %6u-%6u against B %6u-%6u at %7.3f%%\n",
            tailBgn, tailEnd, headBgn, headEnd, align.percentIdentity());

  if (showAlignments())
    align.display(20);

  //  Find the biggest match block, and use the middle of that for the trim point.

  uint32  CC = UINT32_MAX;

  for (uint32 cc=0; cc<align.cigarLength(); cc++) {
    if (align.cigarCode(cc) != '=')
      continue;

    if ((CC == UINT32_MAX) ||                           //  If not a valid biggest block,
        (align.cigarValu(CC) <  align.cigarValu(cc)))   //  or the current block is bigger,
      CC = cc;                                          //  reset the biggest block.
  }

  if (CC != UINT32_MAX) {
    uint32  mid = align.cigarToMapBgn(CC) + (align.cigarToMapEnd(CC) - align.cigarToMapBgn(CC)) / 2;

    keepBgn = align.alignMapB(mid);   //  B is the end of the tig, remember?
    keepEnd = align.alignMapA(mid);   //  A is the start of the tig.
  }

  if (showAlgorithm())
    fprintf(stderr, "  Trim to %u-%u\n", keepBgn, keepEnd);

  //  Test the trimming.
  //
  //  Align the begin 0..keepBgn   to the end   of the tig buffer+keepEnd..len
  //  Align the end   keepEnd..len to the begin of the tig            len..keepBgn+buffer

  bool  headTest = false;
  bool  tailTest = false;

  //  Align trimmed-off start to end of tig.  EDLIB_MODE_HW allows free gaps
  //  at the start/end of the second sequence.
  {
    int32  aBgn = 0;
    int32  aEnd = keepBgn;
    int32  aLen = aEnd - aBgn;

    int32  bBgn = (tailBgn >= 222) ? (tailBgn - 222) : (0);
    int32  bEnd = tailEnd;
    int32  bLen = bEnd - bBgn;

    if (showAlgorithm())
      fprintf(stderr, "  Test   trimmed head %6d-%6d against   kept tail %6d-%6d\n",
              aBgn, aEnd, bBgn, bEnd);

    EdlibAlignResult  result;
    result = edlibAlign(bases + aBgn, aLen,   //  Piece of the start we trim off
                        bases + bBgn, bLen,   //  Should align into the tail we keep
                        edlibNewAlignConfig(0.5 * aLen, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (bBgn + result.endLocations[0] + 1 == keepEnd)
      headTest = true;

    if ((headTest == false) && (showAlgorithm()))
      fprintf(stderr, "  Tested trimmed head %6d-%6d aligns to kept tail %6d-%6d%s\n",
              aBgn, aEnd, bBgn + result.startLocations[0], bBgn + result.endLocations[0] + 1, (headTest) ? "" : " - FAIL");

    edlibFreeAlignResult(result);
  }

  //  Align trimmed-off end to start of tig.
  {
    int32  aBgn = keepEnd;
    int32  aEnd = length;
    int32  aLen = aEnd - aBgn;

    int32  bBgn = headBgn;
    int32  bEnd = (headEnd + 222 <= length) ? (headEnd + 222) : (length);
    int32  bLen = bEnd - bBgn;

    if (showAlgorithm())
      fprintf(stderr, "  Test   trimmed tail %6d-%6d against   kept head %6d-%6d\n",
              aBgn, aEnd, bBgn, bEnd);

    EdlibAlignResult  result;
    result = edlibAlign(bases + aBgn, aLen,   //  Piece of the end we trim off
                        bases + bBgn, bLen,   //  Should align into the head we keep
                        edlibNewAlignConfig(0.5 * aLen, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (bBgn + result.startLocations[0] == keepBgn)
      tailTest = true;

    if (showAlgorithm())
      fprintf(stderr, "  Tested trimmed head %6d-%6d aligns to kept tail %6d-%6d%s\n",
              aBgn, aEnd, bBgn + result.startLocations[0], bBgn + result.endLocations[0] + 1, (tailTest) ? "" : " - FAIL");

    edlibFreeAlignResult(result);
  }

  if ((headTest == true) &&
      (tailTest == true)) {
    if (showAlgorithm())
      fprintf(stderr, "  Tests pass.\n");

    _tig->_trimBgn = keepBgn;
    _tig->_trimEnd = keepEnd;
  }
}



bool
unitigConsensus::generate(tgTig      *tig_,
                          char        algorithm_,
                          char        aligner_,
                          uint32      numIterations_,
                          u32toRead  &reads_) {
  bool  success = false;

  _tig            = tig_;
  _numReads       = _tig->numberOfChildren();
  _numReadsUsable = 0;                                 //  Counted below.

  if (_numReads == 0) {
    fprintf(stderr, "unitigConsensus::initializeGenerate()-- unitig %u has no children.\n", _tig->tigID());
    return false;
  }

  memcpy(_utghpc = new tgPosition [_numReads], _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_utgpos = new tgPosition [_numReads], _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_cnspos = new tgPosition [_numReads], _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_adjpos = new tgPosition [_numReads], _tig->getChild(0), sizeof(tgPosition) * _numReads);

  for (int32 ii=0; ii<_numReads; ii++) {
    if (_tig->getChild(ii)->skipConsensus() == false)   //  Count the number usable
      _numReadsUsable++;                               //  for consensus.
    // initialize low quality flag based on supplied ID cutoff
    _utgpos[ii].isLowQuality(_utgpos[ii].ident() <= _errorRateMaxID);
    _cnspos[ii].isLowQuality(_utgpos[ii].ident() <= _errorRateMaxID);
    _adjpos[ii].isLowQuality(_utgpos[ii].ident() <= _errorRateMaxID);

    addRead(_utgpos[ii].ident(),                       //  Copy read metadata and seq
            _utgpos[ii]._askip, _utgpos[ii]._bskip,    //  from utgpos and seqStore
            _utgpos[ii].isReverse(),                   //  to out _sequences array.
            reads_);

    //tgpos[ii].setMinMax(0, );                        //  Leave input positions alone!
    _cnspos[ii].setMinMax(0, 0);                       //  Clear output positions.
    _adjpos[ii].setMinMax(0, 0);
  }
  assert(_numReads       > 0);
  assert(_numReadsUsable > 0);

  switchToUncompressedCoordinates();

  if      (_numReadsUsable == 1)                       //  Boring, just a single read.
    success = generateSingleton(reads_);
  else if (algorithm_ == 'Q')                          //  Quick consensus, just the template.
    success = generateQuick(reads_);
  else if ((algorithm_ == 'P') ||                      //  Normal utgcns.
           (algorithm_ == 'p'))                        //  'norealign' variant of normal utgcns.
    success = generatePBDAG(aligner_, numIterations_, reads_);

  _tig->_trimBgn = 0;
  _tig->_trimEnd = _tig->length();

  findCoordinates(algorithm_, reads_);
  trimCircular();

  return success;
}
