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

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "kmers.H"
#include "merlin-globals.H"

#include <string>
#include <vector>
#include <algorithm>

#define ERROR_RATE_FACTOR 4
#define NUM_BANDS         2
#define MAX_RETRIES       ERROR_RATE_FACTOR * NUM_BANDS
#define TRIM_BP           500

abSequence::abSequence(uint32  readID,
                       uint32  length,
                       char   *seq,
                       uint32  complemented,
                       uint32  isIgnored,
                       uint32  isONT) {
  _iid              = readID;
  _length           = length;
  _complement       = complemented;
  _isIgnored        = isIgnored;
  _isONT            = isONT;

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

  if (complemented == false)
    for (uint32 ii=0, pp=0; ii<_length; ii++, pp++)
      _bases[pp] = seq[ii];

  else
    for (uint32 ii=_length, pp=0; ii-->0; pp++)
      _bases[pp] = inv[ seq[ii] ];

  _bases[_length] = 0;  //  NUL terminate the strings so we can use them in aligners.
};



unitigConsensus::unitigConsensus(sqStore  *seqStore_,
                                 double    errorRate_,
                                 double    errorRateMax_,
                                 uint32    errorRateMaxID_,
                                 uint32    minOverlap_,
                                 uint32    minCoverage_,
                                 merlinGlobal *merlinGlobal_ = NULL) {

  _seqStore        = seqStore_;

  _tig             = NULL;
  _numReads        = 0;

  _sequencesMax   = 0;
  _sequencesLen   = 0;
  _sequences      = NULL;

  _utgpos          = NULL;
  _cnspos          = NULL;
  _adjpos          = NULL;

  _minOverlap      = minOverlap_;
  _errorRate       = errorRate_;
  _errorRateMax    = errorRateMax_;
  _errorRateMaxID  = errorRateMaxID_;
  _minCoverage     = minCoverage_;
  _merlinGlobal    = merlinGlobal_;
}


unitigConsensus::~unitigConsensus() {

  for (uint32 ss=0; ss<_sequencesLen; ss++)
    delete _sequences[ss];

  delete [] _sequences;
  delete [] _utgpos;
  delete [] _cnspos;
  delete [] _adjpos;
}



void
unitigConsensus::addRead(uint32   readID,
                         uint32   askip, uint32 bskip,
                         bool     complemented,
                         bool     isIgnored, bool isONT,
                         std::map<uint32, sqRead *>     *inPackageRead) {

  //  Grab the read.  If there is no package, load the read from the store.  Otherwise, load the
  //  read from the package.  This REQUIRES that the package be in-sync with the unitig.  We fail
  //  otherwise.  Hey, it's used for debugging only...

  sqRead      *readToDelete = NULL;
  sqRead      *read         = NULL;

  if (inPackageRead == NULL) {
    readToDelete = new sqRead;
    read         = _seqStore->sqStore_getRead(readID, readToDelete);
  }

  else {
    read         = (*inPackageRead)[readID];
  }

  if (read == NULL)
    fprintf(stderr, "Failed to load read %u\n", readID);
  assert(read != NULL);

  //  Grab seq/qlt from the read, offset to the proper begin and length.

  uint32  seqLen = read->sqRead_length() - askip - bskip;
  char   *seq    = read->sqRead_sequence()  + ((complemented == false) ? askip : bskip);

  //  Add it to our list.

  increaseArray(_sequences, _sequencesLen, _sequencesMax, 1);

  _sequences[_sequencesLen++] = new abSequence(readID, seqLen, seq, complemented, isIgnored, isONT);

  delete readToDelete;
}



bool
unitigConsensus::initialize(std::map<uint32, sqRead *>     *reads) {

  if (_numReads == 0) {
    fprintf(stderr, "utgCns::initialize()-- unitig has no children.\n");
    return(false);
  }

  _utgpos = new tgPosition [_numReads];
  _cnspos = new tgPosition [_numReads];
  _adjpos = new tgPosition [_numReads];

  memcpy(_utgpos, _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_cnspos, _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_adjpos, _tig->getChild(0), sizeof(tgPosition) * _numReads);

  //  Clear the cnspos position.  We use this to show it's been placed by consensus.
  //  Guess the number of columns we'll end up with.
  //  Initialize abacus with the reads.

  for (int32 i=0; i<_numReads; i++) {
    _cnspos[i].setMinMax(0, 0);
    _adjpos[i].setMinMax(0, 0);

    addRead(_utgpos[i].ident(),
            _utgpos[i]._askip, _utgpos[i]._bskip,
            _utgpos[i].isReverse(), _utgpos[i].isIgnored(), _utgpos[i].isONT(), 
            reads);
  }

  return(true);
}

void unitigConsensus::switchToUncompressedCoordinates(std::map<uint32, sqRead *>  *reads_) {
    // update coordinates of the tig when needed (when it was assembled in compressed space, in normal space this will be skiped)
    // we do this by tracking the read reaching furthest to the right and keeping its offset + homopolymer coordinate translation
    // the read that overlaps it is then updated to start at that reads uncompressed offset + uncompressed bases based on the overlapping coordinate positions
    //

    // check that we need to do something first
    // just rely on first read
    if ((double)getSequence(0)->length() / (_utgpos[0].max()-_utgpos[0].min()) <= 1.2)
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

      uint32 readCompressedPosition = _utgpos[child].min() - compressedOffset;
      uint32 compressedEnd   = _utgpos[child].max();
      uint32 compressedStart = _utgpos[child].min();

      // find the start position in normal read based on position in compressed read
      uint32 i = 0;
      while (i < nlen && (ntoc[i] < readCompressedPosition))
         i++;

      //if (showAlgorithm())
       //fprintf(stderr, "switchToUncompressedCoordinates()-- I'm trying to find start of child %d %s compressed %d (dist from guide read is %d and in uncompressed in becomes %d)\n", _utgpos[child].ident(), (*reads_)[ _utgpos[child].ident() ]->sqRead_name(), _utgpos[child].min(), readCompressedPosition, i);

      _utgpos[child].setMinMax(i+uncompressedOffset, i+uncompressedOffset+getSequence(child)->length());

      //if (showAlgorithm())
       //fprintf(stderr, "switchToUncompressedCoordinates()-- Updated read %d %s which has length %d to be from %d - %d\n", _utgpos[child].ident(), (*reads_)[ _utgpos[child].ident() ]->sqRead_name(), getSequence(child)->length(), _utgpos[child].min(), _utgpos[child].max());

      // update best end if needed
      if ((ntoc == NULL || compressedEnd > currentEnd) && getSequence(child)->isHiFi()) { 
        // && getSequence(child)->isHiFi()
        nlen  = getSequence(child)->length();
        delete[] ntoc;
        ntoc  = new uint32 [ nlen + 1 ];
        uint32 clen  = homopolyCompress(getSequence(child)->getBases(), nlen, NULL, ntoc);

        currentEnd         = compressedEnd;
        compressedOffset   = compressedStart;
        uncompressedOffset = _utgpos[child].min();

        //if (showAlgorithm())
         //fprintf(stderr, "switchToUncompressedCoordinates()-- Updating guide read to be %d %s which ends at %d. Best before ended at %d. Now my guide is at %d (%d uncompressed)\n", _utgpos[child].ident(), (*reads_)[ _utgpos[child].ident() ]->sqRead_name(), compressedEnd, currentEnd, compressedOffset, uncompressedOffset);
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

char *
unitigConsensus::generateTemplateStitch(std::map<uint32, sqRead *>  *reads_) {
  uint32   minOlap  = _minOverlap;

  //  Initialize, copy the first read.

  uint32       rid      = 0;
  // set rid to be the first non-ignored read idx
  fprintf(stderr, "generateTemplateStitch()-- Starting with read #%d %d %s isIgnored %u %u\n", rid, _utgpos[rid].ident(), (*reads_)[ _utgpos[rid].ident() ]->sqRead_name(), _utgpos[rid].isIgnored(), getSequence(rid)->isIgnored());
  while (getSequence(rid)->isIgnored()) {
    fprintf(stderr, "generateTemplateStitch()-- Skipping read #%d %d %s isIgnored %u %u\n", rid, _utgpos[rid].ident(), (*reads_)[ _utgpos[rid].ident() ]->sqRead_name(), _utgpos[rid].isIgnored(), getSequence(rid)->isIgnored());
    rid++;
  }
  fprintf(stderr, "generateTemplateStitch()-- Using read #%d %d %s isIgnored %u %u\n", rid, _utgpos[rid].ident(), (*reads_)[ _utgpos[rid].ident() ]->sqRead_name(), _utgpos[rid].isIgnored(), getSequence(rid)->isIgnored());

  abSequence  *seq      = getSequence(rid);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();

  uint32       tigmax = std::max(readLen, AS_MAX_READLEN);  //  Must be at least AS_MAX_READLEN, else resizeArray() could fail
  uint32       tiglen = 0;
  char        *tigseq = NULL;

  double	savedErrorRate = _errorRate;

  bool          allowContains  = false;

  allocateArray(tigseq, tigmax, _raAct::clearNew);
  char *readName = (*reads_)[ _utgpos[rid].ident() ]->sqRead_name();
  if (showAlgorithm()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "generateTemplateStitch()-- COPY READ read #%d %d %s (len=%d to %d-%d)\n",
            rid, _utgpos[rid].ident(), readName, readLen, _utgpos[rid].min(), _utgpos[rid].max());
    // fprintf(stderr, "generateTemplateStitch()-- bases: %s\n", getSequence(rid)->getBases());
  }

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
    uint32 firstCandidate = rid;   // Track the first read we can use so we know when we have to give up

    _errorRate = savedErrorRate;
    allowContains = false;

  retryCandidate:
    //  Reset the candidate back to the beginning and scan forward again
    //  if this is the second or later time we hit this, then we will not pick the failed read and will pick prior good
    rid = lastStart;
    nm  = 0;

    //  Pick the next read as the one with the longest extension from all with some minimum overlap
    //  to the template

    if (showAlgorithm())
      fprintf(stderr, "\n");

    for (uint32 ii=rid+1; ii < _numReads; ii++) {
      char *readName = (*reads_)[ _tig->getChild(ii)->ident() ]->sqRead_name();

      //  If contained, move to the next read.  (Not terribly useful to log, so we don't)

      if (_utgpos[ii].isIgnored() || (_utgpos[ii].max() < ePos && allowContains == false))
        continue;

      //  If a bigger end position, save the overlap.  One quirk: if we've already saved an overlap, and this
      //  overlap is thin, don't save the thin overlap.

      bool   thick = (_utgpos[ii].min() + minOlap < ePos);
      bool   first = (nm == 0);
      bool   save  = false;

      if ((nm < _utgpos[ii].max()) && (thick || (first && badToAdd.size() == 0)) && badToAdd.count(ii) == 0) {
        save = true;
        nr   = ii;
        nm   = _utgpos[ii].max();
        if ((first && firstCandidate == 0) || ii < firstCandidate)
           firstCandidate = ii;
      }

      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- read #%d/%d ident %d %s position %d-%d%s%s%s\n",
                ii, _numReads, _utgpos[ii].ident(), readName, _utgpos[ii].min(), _utgpos[ii].max(),
                (save  == true)  ? " SAVE"  : "",
                (thick == false) ? " THIN"  : "",
                (first == true)  ? " FIRST" : "");


      //  If this read has an overlap smaller than we want, stop searching.

      if (thick == false)
        break;
    }

    if (nr == 0) {
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- NO MORE READS TO ALIGN\n");
      break;
    }

    assert(nr != 0);

   if (_utgpos[rid].ident() <= _errorRateMaxID || _utgpos[nr].ident() <= _errorRateMaxID) {
      if (showAlgorithm()) fprintf(stderr, "generateTemplateStitch()-- Increasing threshold because either %d (%d) or %d (%d) name %s is below requested ID %d\n", rid, _utgpos[rid].ident(), nr, _utgpos[nr].ident(), readName, _errorRateMaxID);
      _errorRate = _errorRateMax;
   }
  
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
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- WARNING, increasing min overlap from %d to %u for read %u %s (%d - %d)\n",
                olapLen, std::min(ePos, minOlap), nr, (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), _utgpos[nr].min(), _utgpos[nr].max());
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

    if (showAlgorithm()) {
      fprintf(stderr, "\n");
      fprintf(stderr, "generateTemplateStitch()-- ALIGN template %d-%d (len=%d) to read #%d %d %s %d-%d (len=%d actual=%d at %d-%d)  expecting olap of %d templateLen %d\n",
              tiglen - templateLen, tiglen, templateLen,
              nr, _utgpos[nr].ident(), (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), readBgn, readEnd, readEnd - readBgn, readLen,
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


    int32 maxDifference = std::min(2500, (int32)ceil(0.30*olapLen));
    //  Reset if the edit distance is waay more than our error rate allows or it's very short and we haven't topped out on error.  This seems to be a quirk with
    //  edlib when aligning to N's - I got startLocation = endLocation = 0 and editDistance = alignmentLength.
    if ((double)result.editDistance / result.alignmentLength > bandErrRate || 
        (abs(olapLen-result.alignmentLength) > maxDifference && bandErrRate < _errorRate)) {
      noResult    = true;
      gotResult   = false;
      hitTheStart = false;
      hitTheEnd   = false;
    }

    //  HOWEVER, if we get a result and it's near perfect, declare success even if we hit the start.
    //  These are simple repeats that will align with any overlap.  The one BPW debugged was 99+% A.

    if ((gotResult == true) &&
        (hitTheStart == true) &&
        ((double)result.editDistance / result.alignmentLength < 0.1)) {
      hitTheStart = false;
    }

    //  NOTE that if we hit the end with the same conditions, we should try again, unless there
    //  isn't anything left.  In that case, we don't extend the template.

    if ((gotResult == true) &&
        (hitTheEnd == true) &&
        (moreToExtend == false) &&
        ((double)result.editDistance / result.alignmentLength < 0.1)) {
      hitTheEnd = false;
    }

    //  Now, report what happened, and maybe try again.

    if ((showAlgorithm()) && (noResult == true))
      fprintf(stderr, "generateTemplateStitch()-- FAILED to align - no result\n");

    if ((showAlgorithm()) && (noResult == false))
      fprintf(stderr, "generateTemplateStitch()-- FOUND alignment at %d-%d editDist %d alignLen %d %.f%% expected %d\n",
              result.startLocations[0], result.endLocations[0]+1,
              result.editDistance,
              result.alignmentLength,
              100.0 * result.editDistance / result.alignmentLength, olapLen);

    if ((noResult) || (hitTheStart)) {
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - %s - decrease template size by 10%%\n",
                (noResult == true) ? "no result" : "hit the start");
      tryAgain = true;
      templateSize -= 0.10;
    }

    if ((noResult) || (hitTheEnd && moreToExtend)) {
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - %s - increase read size by 10%%\n",
                (noResult == true) ? "no result" : "hit the end");
      tryAgain = true;
      extensionSize += 0.10;
    }

    if (templateSize < 0.01) {
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - no more template to remove!");

      if (bandErrRate + _errorRate / ERROR_RATE_FACTOR > _errorRate || _errorRate - bandErrRate < 1e-9) {
        if (showAlgorithm()) fprintf(stderr, "  Fail!\n");
        tryAgain = false;
        olapLen = origLen;
      }
      else {
        if (showAlgorithm())
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

      if (showAlgorithm()) {
        fprintf(stderr, "generateTemplateStitch()--\n");
        fprintf(stderr, "generateTemplateStitch()-- Aligned template %d-%d to read %u %s %d-%d; copy read %d-%d to template.\n",
                tiglen - templateLen, tiglen, nr, (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), readBgn, readEnd, readEnd, readLen);
        fprintf(stderr, "generateTemplateStitch()-- New position for read %d %s is %d-%d\n", _utgpos[nr].ident(), (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), _cnspos[nr].min(), _cnspos[nr].max());
      }
    } else {
      // try to go back and pick a different read to extend with if that's an option
      // it's not an option if we're already the next read
      // if we're out of reads but didn't try contains, try again anyway, after this we give up and trim the template below
      // note that we check if this read was already tried because when we go back and re-try the contained reads, we may end up not at the firstCandidate (since it was uncontained) but if we run out of new reads to try we need to stop
      bool isAlreadyBad = badToAdd.count(nr) != 0;
      badToAdd.insert(nr);
      if ((rid != firstCandidate && !isAlreadyBad) || allowContains == false) {
         allowContains = (allowContains || rid == firstCandidate || isAlreadyBad);  // we hit the end so now we can try to allow contains
         olapLen = origLen;
         if (showAlgorithm()) {
            fprintf(stderr, "generateTemplateStitch()--\n");
            fprintf(stderr, "generateTemplateStitch()-- FAILED to align read #%d ident %d %s to last added %d ident %d, will go back and re-try an earlier read until we hit %d and allowing contains is %d\n", rid, _utgpos[nr].ident(), (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), lastStart, _utgpos[lastStart].ident(), firstCandidate, allowContains);
            fprintf(stderr, "generateTemplateStitch()--\n");
         }
         edlibFreeAlignResult(result);
         goto retryCandidate;
      }
      else if (tiglen > minOlap && origLen > minOlap + TRIM_BP ) {
         int32 trimbp = std::min(tiglen - minOlap,  origLen - minOlap - TRIM_BP);
         assert(trimbp > 0);
         assert(origLen > minOlap + trimbp);
         assert(ePos > trimbp);

         if (showAlgorithm()) {
            fprintf(stderr, "generateTemplateStitch()--\n");
            fprintf(stderr, "generateTemplateStitch()-- FAILED to align read #%d ident %d %s to last added %d ident %d, will try to trim template of %d bp by %d bases\n", rid, _utgpos[nr].ident(), (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), lastStart, _utgpos[lastStart].ident(), tiglen, trimbp);
         }

         ePos -= trimbp; 
         tigseq[tiglen-trimbp] = 0;
         tiglen=strlen(tigseq);
         firstCandidate = 0;
         badToAdd.clear();

         if (showAlgorithm()) {
            fprintf(stderr, "generateTemplateStitch()-- Trimmed template to %d bp\n", tiglen);
            fprintf(stderr, "generateTemplateStitch()--\n");
         }
         edlibFreeAlignResult(result);
         goto retryCandidate;
      }

      readBgn = 0;
      readEnd = olapLen;

      if (showAlgorithm()) {
        fprintf(stderr, "generateTemplateStitch()--\n");
        fprintf(stderr, "generateTemplateStitch()-- Alignment failed, use original overlap; copy read %s %d-%d to template.\n",
                (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), readEnd, readLen);
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

    if (showAlgorithm())
      fprintf(stderr, "generateTemplateStitch()-- Append read %s from %u to %u, starting at tiglen %u\n", (*reads_)[ _utgpos[nr].ident() ]->sqRead_name(), readEnd, readLen, tiglen);

    for (uint32 ii=readEnd; ii<readLen; ii++)
      tigseq[tiglen++] = fragment[ii];

    tigseq[tiglen] = 0;

    assert(tiglen < tigmax);

    ePos = _utgpos[rid].max();

    if (showAlgorithm())
      fprintf(stderr, "generateTemplateStitch()-- Template now length %d, expected %d, difference %7.4f%%\n",
              tiglen, ePos, 200.0 * ((int32)tiglen - (int32)ePos) / ((int32)tiglen + (int32)ePos));
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

  updateReadPositions();
  // char tigNameTemplateFa[sizeof "tig00000000_template.fa"];
  // char tigName[sizeof "tig00000000"];
  // sprintf(tigNameTemplateFa, "tig%08d_template.fa", _tig->tigID());
  // sprintf(tigName, "tig%08d", _tig->tigID());
  // FILE *fp = fopen(tigNameTemplateFa, "w");
  // fprintf(fp, ">%s\n%s\n", tigName, tigseq);
  // fclose(fp);
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
           bool               verbose) {

  EdlibAlignResult align;

  int32   padding        = std::min(250, (int32)ceil(fragmentLength * 0.05)); 
  double  bandErrRate    = errorRate / ERROR_RATE_FACTOR;
  bool    aligned        = false;
  double  alignedErrRate = 0.0;

  //  Decide on where to align this read.

  int32  tigbgn = std::max((int32)0,      (int32)floor(utgpos.min() - padding));
  int32  tigend = std::min((int32)tiglen, (int32)floor(utgpos.max() + padding));

  if (verbose)
    fprintf(stderr, "alignEdLib()-- align read %7u eRate %.4f at %9d-%-9d\n", utgpos.ident(), bandErrRate, tigbgn, tigend);

  if (tigend < tigbgn) {
    fprintf(stderr, "alignEdLib()-- WARNING: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
            tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
    // try to align it to full
    tigbgn = 0;
    tigend = (utgpos.max() < 0 ? tiglen : std::min((int32)tiglen, (int32)utgpos.max()));
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
    aligned        = (alignedErrRate <= bandErrRate);
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
       tigbgn = std::max((int32)0,      tigbgn - 25 * padding);
       tigend = std::min((int32)tiglen, tigend + 25 * padding);
    }

    // let the band increase without increasing error rate for a while, if we give up increase error rate
    // and reset bnad
    if (ii != 0 && ii % NUM_BANDS == 0) {
       bandErrRate += errorRate / ERROR_RATE_FACTOR;
       tigbgn = std::max((int32)0,      (int32)floor(utgpos.min() - padding));
       tigend = std::min((int32)tiglen, (int32)floor(utgpos.max() + padding));
    }

    if (tigend < tigbgn) {
       fprintf(stderr, "alignEdLib()-- WARNING: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
               tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
       // try to align it to full
       tigbgn = 0;
       tigend = (utgpos.max() < 0 ? tiglen : std::min((int32)tiglen, (int32)utgpos.max()));
       fprintf(stderr, "alignEdLib()-- WARNING: updated tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
               tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
    }
    assert(tigend > tigbgn);

    edlibFreeAlignResult(align);

    if (verbose)
      fprintf(stderr, "alignEdLib()--                  read %7u eRate %.4f at %9d-%-9d\n", utgpos.ident(), bandErrRate, tigbgn, tigend);

    align = edlibAlign(fragment, strlen(fragment),
                       tigseq + tigbgn, tigend - tigbgn,
                       edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (align.alignmentLength > 0) {
      alignedErrRate = (double)align.editDistance / align.alignmentLength;
      aligned        = (alignedErrRate <= bandErrRate);
      if (verbose)
        fprintf(stderr, "alignEdLib()-- - ALIGNED read %7u %.4f at %9d-%-9d\n", utgpos.ident(), alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
    } else {
      if (verbose)
        fprintf(stderr, "\n");
    }
  }

  if (aligned == false) {
    fprintf(stderr, "alignEdLib()-- failed to align read %7u.\n", utgpos.ident());
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

EdlibAlignResult
alignEdLibFindCoordinates(dagAlignment    &aln,
            bool            &success,
            tgPosition      &cnspos,
            uint32          isONT,
            char            *fragment,
            uint32          fragmentLength,
            char            *tigseq,
            uint32          tiglen,
            double          errorRate,
            bool            verbose,
            char            *readName) {

  EdlibAlignResult align;
  int32 padding;
  if (isONT) {
    padding        = std::min(2000, (int32)ceil(fragmentLength * 0.05));
  } else {
    padding        = std::min(250, (int32)ceil(fragmentLength * 0.05));
  }

  double  bandErrRate    = errorRate / ERROR_RATE_FACTOR;
  bool    aligned        = false;
  double  alignedErrRate = 0.0;

  //  Decide on where to align this read.

  int32  tigbgn = std::max((int32)0,      (int32)floor(cnspos.min() - padding));
  int32  tigend = std::min((int32)tiglen, (int32)floor(cnspos.max() + padding));

  if (verbose)
    fprintf(stderr, "alignEdLibFindCoordinates()-- align read %7u %s eRate %.4f adjusted pos %d-%d at %d-%d\n", cnspos.ident(), readName, bandErrRate, cnspos.min(), cnspos.max(), tigbgn, tigend);

  if (tigend < tigbgn) {
    fprintf(stderr, "alignEdLibFindCoordinates()-- WARNING: read %s tigbgn %d > tigend %d - tiglen %d cnspos %d-%d padding %d\n",
            readName, tigbgn, tigend, tiglen, cnspos.min(), cnspos.max(), padding);
    // try to align it to full
    tigbgn = 0;
    tigend = (cnspos.max() < 0 ? tiglen : std::min((int32)tiglen, (int32)cnspos.max()));
    fprintf(stderr, "alignEdLibFindCoordinates()-- WARNING: read %s updated tigbgn %d > tigend %d - tiglen %d cnspos %d-%d padding %d\n",
            readName, tigbgn, tigend, tiglen, cnspos.min(), cnspos.max(), padding);
  }
  assert(tigend > tigbgn);

  //  Align!  If there is an alignment, compute error rate and declare success if acceptable.

  align = edlibAlign(fragment, fragmentLength,
                     tigseq + tigbgn, tigend - tigbgn,
                     edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

  if (align.alignmentLength > 0) {
    alignedErrRate = (double)align.editDistance / align.alignmentLength;
    aligned        = (alignedErrRate <= bandErrRate);
    if (verbose)
      fprintf(stderr, "alignEdLibFindCoordinates()-- - ALIGNED read %s, %.4f at %9d-%9d\n", readName, alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
  } else {
    if (verbose)
      fprintf(stderr, "\n");
  }

  for (uint32 ii=1 /*the 0ths was above*/; ((ii < MAX_RETRIES) && (aligned == false)); ii++) {
    tigbgn = std::max((int32)0,      tigbgn - 5 * padding);
    tigend = std::min((int32)tiglen, tigend + 5 * padding);
    // last attempt make a very wide band
    if ((ii+1) % NUM_BANDS == 0) {
       tigbgn = std::max((int32)0,      tigbgn - 25 * padding);
       tigend = std::min((int32)tiglen, tigend + 25 * padding);
    }

    // let the band increase without increasing error rate for a while, if we give up increase error rate
    // and reset band
    if (ii != 0 && ii % NUM_BANDS == 0) {
       bandErrRate += errorRate / ERROR_RATE_FACTOR;
       tigbgn = std::max((int32)0,      (int32)floor(cnspos.min() - padding));
       tigend = std::min((int32)tiglen, (int32)floor(cnspos.max() + padding));
    }

    if (tigend < tigbgn) {
       fprintf(stderr, "alignEdLibFindCoordinates()-- WARNING: read %s tigbgn %d > tigend %d - tiglen %d cnspos %d-%d padding %d\n",
               readName, tigbgn, tigend, tiglen, cnspos.min(), cnspos.max(), padding);
       // try to align it to full
       tigbgn = 0;
       tigend = (cnspos.max() < 0 ? tiglen : std::min((int32)tiglen, (int32)cnspos.max()));
       fprintf(stderr, "alignEdLibFindCoordinates()-- WARNING: read %s updated tigbgn %d > tigend %d - tiglen %d cnspos %d-%d padding %d\n",
               readName, tigbgn, tigend, tiglen, cnspos.min(), cnspos.max(), padding);
    }
    assert(tigend > tigbgn);

    edlibFreeAlignResult(align);

    if (verbose)
      fprintf(stderr, "alignEdLibFindCoordinates()--                 read %s   eRate %.4f at %9d-%-9d\n", readName, bandErrRate, tigbgn, tigend);

    align = edlibAlign(fragment, strlen(fragment),
                       tigseq + tigbgn, tigend - tigbgn,
                       edlibNewAlignConfig(bandErrRate * fragmentLength, EDLIB_MODE_HW, EDLIB_TASK_PATH));

    if (align.alignmentLength > 0) {
      alignedErrRate = (double)align.editDistance / align.alignmentLength;
      aligned        = (alignedErrRate <= bandErrRate);
      if (verbose)
        fprintf(stderr, "alignEdLibFindCoordinates()-- - ALIGNED read %s, %.4f at %9d-%-9d\n", readName, alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
    } else {
      if (verbose)
        fprintf(stderr, "\n");
    }
  }

  if (aligned == false) {
    success = aligned;
    return align;
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

  if (tigbgn + align.endLocations[0] + 1 > tiglen)
    fprintf(stderr, "ERROR:  alignment from %d to %d, but tiglen is only %d\n", tigbgn + align.startLocations[0] + 1, tigbgn + align.endLocations[0] + 1, tiglen);
  assert(tigbgn + align.endLocations[0] + 1 <= tiglen);

  success = aligned;
  return align;
}



bool
unitigConsensus::initializeGenerate(tgTig                       *tig_,
                                    std::map<uint32, sqRead *>  *reads_) {

  _tig      = tig_;
  _numReads = _tig->numberOfChildren();

  if (initialize(reads_) == false) {
    fprintf(stderr, "Failed to initialize for tig %u with %u children\n", _tig->tigID(), _tig->numberOfChildren());
    return(false);
  }

  switchToUncompressedCoordinates(reads_);

  return(true);
}


bool
unitigConsensus::generatePBDAG(tgTig                       *tig_,
                               char                         aligner_,
                               std::map<uint32, sqRead *>  *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  //  Build a quick consensus to align to.

  char   *tigseq = generateTemplateStitch(reads_);
  uint32  tiglen = strlen(tigseq);

  if (showAlgorithm())
    fprintf(stderr, "Generated template of length %d\n", tiglen);

  //  Compute alignments of each sequence in parallel

  if (showAlgorithm())
    fprintf(stderr, "Aligning reads.\n");

  dagAlignment *aligns = new dagAlignment [_numReads];
  uint32        pass = 0;
  uint32        fail = 0;

#pragma omp parallel for schedule(dynamic)
  for (uint32 ii=0; ii<_numReads; ii++) {
    abSequence  *seq      = getSequence(ii);
    bool         aligned  = false;
    char *rn = (*reads_)[ _tig->getChild(ii)->ident() ]->sqRead_name();

    assert(aligner_ == 'E');  //  Maybe later we'll have more than one aligner again.

    aligned = alignEdLib(aligns[ii],
                         _utgpos[ii],
                         seq->getBases(), seq->length(),
                         tigseq, tiglen,
                         _errorRateMax,
                         showAlgorithm());

    if (aligned == false) {
      if (showAlgorithm())
        fprintf(stderr, "generatePBDAG()--    read %7u %s FAILED\n", _utgpos[ii].ident(), rn);

      fail++;

      continue;
    }

    pass++;
  }

  if (showAlgorithm())
    fprintf(stderr, "generatePBDAG()--    TIG %u read alignment: %d failed, %d passed.\n", _tig->tigID(), fail, pass);

  //  Construct the graph from the alignments.  This is not thread safe.

  if (showAlgorithm())
    fprintf(stderr, "Constructing graph\n");

  fprintf(stderr, "TIG %u tiglen %u\n", _tig->tigID(), tiglen);
  AlnGraphBoost ag(std::string(tigseq, tiglen));

  for (uint32 ii=0; ii<_numReads; ii++) {
    _cnspos[ii].setMinMax(aligns[ii].start, aligns[ii].end);

    if ((aligns[ii].start == 0) &&
        (aligns[ii].end   == 0))
      continue;

    // here we ignore flagged reads, not added to the DAG
    if (!_utgpos[ii].isIgnored())
      ag.addAln(aligns[ii]);

    aligns[ii].clear();
  }

  delete [] aligns;

  if (showAlgorithm())
    fprintf(stderr, "Merging graph\n");

  //  Merge the nodes and call consensus
  ag.mergeNodes();

  if (showAlgorithm())
    fprintf(stderr, "Calling consensus\n");

  //FIXME why do we have 0weight nodes (template seq w/o support even from the read that generated them)?
  _templateToCns = new uint32_t[tiglen + 1];
  _templateLength = tiglen;
  for (uint32 i=0; i < tiglen; i++) {
    _templateToCns[i] = UINT32_MAX;
  }
  _startTrim = 0;
  std::string cns = ag.consensusNoSplitMap(_minCoverage, _templateToCns, &_startTrim);
  delete [] tigseq;

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

  return(true);
}



bool
unitigConsensus::generateQuick(tgTig                       *tig_,
                               std::map<uint32, sqRead *>  *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  //  Quick is just the template sequence, so one and done!

  char   *tigseq = generateTemplateStitch(reads_);
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
unitigConsensus::generateSingleton(tgTig                       *tig_,
                                   std::map<uint32, sqRead *>  *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  assert((getNumHiFiReads(tig_, reads_) == 1));
  // assert(_numReads == 1);

  //  Copy the single read to the tig sequence.
  uint32 rid = 0;
  while (rid < tig_->numberOfChildren()) {
    if (getSequence(rid)->isHiFi()) {
      break;
    }
    rid++;
  }
  
  fprintf(stderr, "generateSingleton()-- read rid %u id %u name %s isHiFi %u\n", rid, tig_->getChild(rid)->ident(), (*reads_)[tig_->getChild(rid)->ident()]->sqRead_name(), tig_->getChild(rid)->isHiFi());

  abSequence  *seq      = getSequence(rid);
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




//  This has more canu dependency than edlib dependency, so it's here instead of in edlib.c
//
uint32
edlibAlignmentToCanu(stuffedBits *align,
                     const unsigned char* alignment,
                     int alignmentLength,
                     int tgtStart, int tgtEnd,
                     int qryStart, int qryEnd) {
  int32  qryPos = qryStart;
  int32  tgtPos = tgtStart;

  //  Count the length of the alignment.

  int32  nMatch  = 0;   //  Counting the number of match/mismatch in a delta block.
  int32  nBlocks = 0;   //  Number of blocks in the delta encoding.

  //  Code the alignment.

  for (int32 a=0; a<alignmentLength; a++) {
    assert(qryPos <= qryEnd);
    assert(tgtPos <= tgtEnd);

    //  Match or mismatch.
    if ((alignment[a] == EDLIB_EDOP_MATCH) ||
        (alignment[a] == EDLIB_EDOP_MISMATCH)) {
      nMatch++;

      qryPos++;
      tgtPos++;
    }

    //  Insertion in target.
    else if (alignment[a] == EDLIB_EDOP_INSERT) {
      align->setEliasDelta(nMatch + 1);
      align->setBit(0);

      nMatch = 0;
      nBlocks++;

      qryPos++;
    }

    //  Insertion in query.
    else if (alignment[a] == EDLIB_EDOP_DELETE) {
      align->setEliasDelta(nMatch + 1);
      align->setBit(1);

      nMatch = 0;
      nBlocks++;

      tgtPos++;
    }
  }

  //  Don't forget the last match/mismatch block.
  //align->setEliasDelta(nMatch + 1);
  //align->setBit(1);

  //nMatch = 0;
  //nBlocks++;

  //fprintf(stderr, "Align of length %d with %d gaps stored in %lu bits.\n", alignmentLength, nBlocks, align->getLength());

  return(nBlocks);
}


tgPosition
unitigConsensus::adjustPosition(tgPosition original, char* readname) {
  int32 original_min = original.min();
  int32 original_max = original.max();
  int32 adjusted_min = 0;
  int32 adjusted_max = 0;

  fprintf(stderr, "adjustPosition()-- read %s original min %u original max %u _templateLength = %u _tig->length() = %u _startTrim %u\n", readname, original_min, original_max, _templateLength, _tig->length(), _startTrim);

  if (original_min > _templateLength) {
    //fprintf(stderr, "adjustPosition()-- read %s original min %u > _templateLength %u - return (-1,-1) don't align\n", readname, original_min, _templateLength);
    tgPosition pos = tgPosition();
    pos.initialize();
    pos.setMinMax(-1, -1);
    return pos;
  }

  bool endOffTheTig = false;

  if (original_max > _templateLength) {
    endOffTheTig = true;
  }

  adjusted_max = endOffTheTig ? _templateLength : _templateToCns[original_max];
  if (adjusted_max > _startTrim && adjusted_max - _startTrim > _tig->length()) {
    endOffTheTig = true;
  }

  while ( (adjusted_max == UINT32_MAX || adjusted_max == 0) && original_max < _templateLength) {
    adjusted_max = _templateToCns[++original_max];
  }

  if (original_max == _templateLength) {
    endOffTheTig = true;
    adjusted_max = _tig->length();
  } else if (_startTrim > 0 && adjusted_max > _startTrim && !endOffTheTig) {
    adjusted_max -= _startTrim;
  }

  adjusted_min = _templateToCns[original_min];


  while (adjusted_min == UINT32_MAX && original_min > 0) {
    adjusted_min = _templateToCns[--original_min];
  }

  if (original_min == 0) {
    adjusted_min = 0;
  } else if (_startTrim > 0) {
    if (adjusted_min > _startTrim) {
      adjusted_min = adjusted_min - _startTrim;
    } else {
      adjusted_min = 0;
    }
  }

  tgPosition pos = tgPosition();
  pos.initialize();
  pos.setMinMax(adjusted_min, adjusted_max);
  fprintf(stderr, "adjustPosition()-- read %s return adjust position min %u max %u\n", readname, pos.min(), pos.max());

  return pos;
}

//  Align each read to the region of the consensus sequence the read claims
//  to be from, adjusted by the map we created while constructing the cns.
//  Output these alignments to SAM files, along with some tags we create.
void
unitigConsensus::findCoordinates(uint32_t                     oneCopyPeak,
                                 char                         aligner_,
                                 std::map<uint32, sqRead *>   *reads_) {
  #warning Kmer normalization assumes a genome where the 2-copy peak is taller than the 1-copy peak to estimate the haploid coverage, if this does not hold true this may return the 1-copy peak / 2.
  fprintf(stderr, "findCoordinates()-- One copy peak estimate used to normalize kmer counts: " F_U32 "x. WARNING: assumes a genome where 2-copy peak is taller than 1-copy peak! \n", oneCopyPeak);
  char filenameSam[sizeof "tig00000000.sam"];
  char tigName[sizeof "tig00000000"];
  sprintf(tigName, "tig%08d", _tig->tigID());
  sprintf(filenameSam, "tig%08d.sam", _tig->tigID());
  FILE *fp = fopen(filenameSam, "w");
  fprintf(fp, "@SQ\tSN:%s\tLN:%u\n", tigName, _tig->length());
  fclose(fp);
  samFile *samFp = sam_open(filenameSam, "a+");
  if (samFp == NULL) {
    fprintf(stderr, "Failed to open output file!\n");
    exit(1);
  }

  sam_hdr_t *samHeader = sam_hdr_init();
  samHeader->n_targets = 1;
  samHeader->target_len = (uint32_t*)malloc(samHeader->n_targets * sizeof(uint32_t));
  samHeader->target_len[0] = _tig->length();
  samHeader->target_name = (char**)malloc(samHeader->n_targets * sizeof(char*));
  char targetName[100];
  sprintf(targetName, "tig%08d", _tig->tigID());
  samHeader->target_name[0] = strdup(targetName);
  // int ret = bam_hdr_write(bamFp->fp.bgzf, bamHeader);
  int ret = sam_hdr_write(samFp, samHeader);
  if (ret < 0) {
    fprintf(stderr, "Failed to write header to SAM file!\n");
    exit(1);
  }


  if (showPlacement()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "TIG %u length %u\n", _tig->tigID(), _tig->length());
    fprintf(stderr, "\n");
  }
  

  if ( _tig->length() == 0) {
    fprintf(stderr, "TIG %u length is %u, exit findCoordinates()\n", _tig->tigID(), _tig->length());
    sam_close(samFp);
    sam_hdr_destroy(samHeader);
    delete[] _templateToCns;
    return;
  }

  _tig->_childDeltaBitsLen  = 1;
  _tig->_childDeltaBits     = new stuffedBits();
  dagAlignment *alns        = new dagAlignment [_numReads];
  EdlibAlignResult *aligns  = new EdlibAlignResult [_numReads];
  uint32           pass     = 0;
  uint32           fail     = 0;
  char             *tigseq  = _tig->bases();
  uint32           tiglen   = strlen(tigseq);


  #pragma omp parallel for schedule(dynamic)
  for (uint32 ii=0; ii<_numReads; ii++) {
    abSequence  *seq      = getSequence(ii);
    bool        aligned   = false;
    int32       tigbgn    = 0;

    assert(aligner_ == 'E');  //  Maybe later we'll have more than one aligner again.
    char *rn = (*reads_)[ _tig->getChild(ii)->ident() ]->sqRead_name();

    // adjust coordinates of reads based on map we constructed during consensus
    if (getNumHiFiReads(_tig, reads_) > 1) {
      if (_cnspos[ii].min() == 0 && _cnspos[ii].max() == 0) {
        _adjpos[ii] = adjustPosition(_utgpos[ii], rn);
      } else {
        _adjpos[ii] = adjustPosition(_cnspos[ii], rn);
      }
    } else { // if singleton, just align to whole tig
      if (_cnspos[ii].min() == 0 && _cnspos[ii].max() == 0) {
        _adjpos[ii].setMinMax(0, tiglen);
      } else {
        _adjpos[ii].setMinMax(0, tiglen);
      }
    }

    if (_adjpos[ii].min() == -1 && _adjpos[ii].max() == -1) {
        fail++;
        continue;
    }

    EdlibAlignResult align = alignEdLibFindCoordinates(alns[ii],
                                         aligned,
                                         _adjpos[ii],
                                         _tig->getChild(ii)->isONT(),
                                         seq->getBases(),
                                         seq->length(),
                                         tigseq,
                                         tiglen,
                                         _errorRateMax,
                                         showAlgorithm(),
                                         rn);
    if (aligned == false) {
      if (showAlgorithm())
        fprintf(stderr, "findCoordinates()--    read %7u %s FAILED\n", _utgpos[ii].ident(), rn);

      fail++;
      continue;
    }

    aligns[ii] = align;
    pass++;
    if (showAlgorithm())
      fprintf(stderr, "findCoordinates()-- read %7u %s PASSED - PBDAG placed at %u %u realigned to %u %u, final pos is %u %u\n", _utgpos[ii].ident(), rn, _cnspos[ii].min(), _cnspos[ii].max(), _adjpos[ii].min(), _adjpos[ii].max(), alns[ii].start - 1, alns[ii].end);
  }

  if (getNumHiFiReads(_tig, reads_) > 1) {
    delete[] _templateToCns;
  }

  fprintf(stderr, "findCoordinates()--    TIG %u read alignment: %d failed, %d passed.\n", _tig->tigID(), fail, pass);

  for (uint32 ii=0; ii<_numReads; ii++) {
    abSequence   *seq      = getSequence(ii);
    char         *readSeq = seq->getBases();
    uint32       readLen = seq->length();
    EdlibAlignResult align = aligns[ii];
  
    if ((alns[ii].start == 0) &&
        (alns[ii].end   == 0))
      continue;
    
    char *cigar = edlibAlignmentToCigar(align.alignment,
                                        align.alignmentLength, 
                                        EDLIB_CIGAR_EXTENDED);
    bam1_t *bamRecord = bam_init1();
	  char *cigarCopy = strdup(cigar);
    bam1_t *cigarRecord = bam_init1();
    ssize_t numCigarOpsBam = bam_parse_cigar(cigarCopy,
                              NULL,
                              cigarRecord);
    if (numCigarOpsBam < 0) {
      fprintf(stderr, "Failed to process cigar string to numCigarOps!\n");
      exit(1);
    }
    uint32_t *cigarData = (uint32_t *) malloc((numCigarOpsBam + 1) * sizeof(uint32_t));
    size_t numOpsCast = static_cast<size_t>(numCigarOpsBam);
    ssize_t numCigarOpsSam = sam_parse_cigar(
      cigarCopy,
      NULL,
      &cigarData,
      &numOpsCast);
        
	  char *readName = (*reads_)[ _tig->getChild(ii)->ident() ]->sqRead_name();
    char qual[readLen + 1];
    for (int i = 0; i < readLen; i++) {
      qual[i] = 'F';
    }
    qual[readLen] = '\0';
    const char* qualPtr = (const char*) qual;

    // create bam record with alignment
    int ret = bam_set1(
      bamRecord,
      strlen(readName),
      readName,
      seq->isForward() ? 0 : 16, 
      0,
      alns[ii].start - 1,
      255,
      numCigarOpsSam,
      cigarData,
      0,
      0,
      _tig->length(),
      readLen,
      readSeq,
      qualPtr,
      0);
    if (ret < 0) { 
      fprintf(stderr, "Failed to create bam record!\n");
      exit(1);
    }

    free(cigarCopy);

    // get clipping information for IH tag
    char leftClipStr [11] = "0";
    char rightClipStr [11] = "0";
    uint32_t leftClipLen = 0;
    uint32_t rightClipLen = 0;
    char leftOp = '\0';
    char rightOp = '\0';
    uint32_t first_pos = 0;
    uint32_t totClipping = 0;
    while (first_pos < strlen(cigar)) {
      if (!isDecDigit(cigar[first_pos])) {
        leftOp = cigar[first_pos];
        if (leftOp != '=' && leftOp != 'M') {
          strncpy(leftClipStr, cigar, first_pos);
          leftClipLen = std::atol(leftClipStr);
          leftOp = cigar[first_pos];
          totClipping += leftClipLen;
        }
        break;
      }
      first_pos++;
    }

    uint32_t last_pos = strlen(cigar) - 1;
    rightOp = cigar[last_pos--];
    if (rightOp != '=' && rightOp != 'M') {
        while (last_pos > 0) {
          if (!isDecDigit(cigar[last_pos])) {
            strncpy(rightClipStr, cigar + last_pos + 1, strlen(cigar) - last_pos - 2);
            rightClipLen = std::atol(rightClipStr);
            totClipping += rightClipLen;
            break;
          }
          last_pos--;
        }
    }

    std::string totClippingStr = std::to_string(totClipping);

    // get normalized kmer counts for CO tag
    if (_merlinGlobal) {
      std::string merylStr;
      uint32_t count;
      float normCount;
      std::string countStr;
      kmerIterator  kiter(readSeq, readLen);
      bool first = true;
      while (kiter.nextMer()) {
        count = _merlinGlobal->lookupMarker(kiter.fmer(), kiter.rmer());
        normCount = count / float(oneCopyPeak);
        countStr = std::to_string(normCount);
        if (first) {
          merylStr += countStr;
          first = false;
        } else {
          merylStr += ',' + countStr;
        }
      }

      // write tags to bam record

      uint8_t* auxMerylData = (uint8_t*) malloc(sizeof(uint8_t) * (merylStr.length() + 1));
      memcpy(auxMerylData, merylStr.c_str(), merylStr.length() + 1);    
      const char* COtag = "CO";
      ret = bam_aux_append(
        bamRecord,
        COtag,
        'Z',
        merylStr.length() + 1,
        auxMerylData);
      if (ret < 0) { 
        fprintf(stderr, "Failed to add meryl data tag!\n");
        exit(1);
      }
      free(auxMerylData);
    }

    uint8_t* auxClippingData = (uint8_t*) malloc(sizeof(uint8_t) * (totClippingStr.length() + 1));
    memcpy(
      auxClippingData,
      totClippingStr.c_str(),
      totClippingStr.length() + 1);
    const char* IHtag = "IH";
    ret = bam_aux_append(bamRecord,
      IHtag,
      'Z',
      totClippingStr.length() + 1,
      auxClippingData);
    if (ret < 0) { 
      fprintf(stderr, "Failed to add clipping tag!\n");
      exit(1);
    }
    free(auxClippingData);

    // write bam record to file
    ret = sam_write1(
      samFp,
      samHeader,
      bamRecord);
    if (ret < 0) {
      fprintf(stderr, "Failed to write sam record!\n");
      exit(1);
    }
    bam_destroy1(cigarRecord);
    bam_destroy1(bamRecord);
  	free(cigarData);
    delete[] cigar;
    edlibFreeAlignResult(align);
  }
  delete[] alns;
  delete[] aligns;
  sam_close(samFp);
  sam_hdr_destroy(samHeader);
}


void
unitigConsensus::findRawAlignments(void) {

#if 0
  for (uint32 ii=0; ii<_numReads; ii++) {
    fprintf(stderr, "read %4u original %9u-%9u aligned %9u-%9u\n",
            ii,
            _utgpos[ii].min(), _utgpos[ii].max(),
            _cnspos[ii].min(), _cnspos[ii].max());
  }
#endif
}



void
unitigConsensus::trimCircular(void) {
  char    *bases  = _tig->bases();
  uint32   length = _tig->length();

  if (_tig->_suggestCircular == false)
    return;

  if (_tig->_circularLength == 0)
    return;

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

uint32
unitigConsensus::getNumHiFiReads(tgTig                      *tig_, 
                                 std::map<uint32, sqRead *> *reads_) {
  int32 rid = 0;
  int32 numHiFiReads = 0;
  while (rid < tig_->numberOfChildren()) {
    //fprintf(stderr, "getNumHiFiReads()-- read %u %s isHiFi: %u\n", tig_->getChild(rid)->ident(), (*reads_)[tig_->getChild(rid)->ident()]->sqRead_name(), tig_->getChild(rid)->isHiFi());
    if (tig_->getChild(rid)->isHiFi()) {
      numHiFiReads++;
    }
    rid++;
  }
  return numHiFiReads;
}


bool
unitigConsensus::generate(tgTig                       *tig_,
                          char                         algorithm_,
                          char                         aligner_,
                          double                       maxCov,
                          tgTigStashed                 S,
                          uint32_t                     oneCopyPeak,
                          std::map<uint32, sqRead *>  *reads_) {
  bool  success = false;

  //  Flag excess coverage to be ignored.  Singletons report no logging.
  tig_->flagContains(maxCov, S, reads_);
  if (getNumHiFiReads(tig_, reads_) == 1) {
      success = generateSingleton(tig_, reads_);
  }
  else if (tig_->numberOfChildren() == 1) {
    fprintf(stderr, "generate()-- TIG %u has only one ONT read, don't generate ONT singleton\n", tig_->tigID());
  }


  else if (algorithm_ == 'Q') {
    success = generateQuick(tig_, reads_);
    // findCoordinates(oneCopyPeak, aligner_, reads_);
  }

  else if ((algorithm_ == 'P') ||   //  Normal utgcns.
           (algorithm_ == 'p')) {   //  'norealign'
    success = generatePBDAG(tig_, aligner_, reads_);
  }

  if (success) {
    _tig->_trimBgn = 0;
    _tig->_trimEnd = _tig->length();

    if (algorithm_ == 'P') {
      findCoordinates(oneCopyPeak, aligner_, reads_);
      // findRawAlignments();
    }

    // trimCircular();
  }

  return(success);
}