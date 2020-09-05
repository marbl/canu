
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

// for pbdagcon
#include "Alignment.H"
#include "AlnGraphBoost.H"
#include "edlib.H"
#include "align-ssw.H"
#include "align-ssw-driver.H"

#include <set>

using namespace std;


abSequence::abSequence(uint32  readID,
                       uint32  length,
                       char   *seq,
                       uint32  complemented) {
  _iid              = readID;
  _length           = length;
  _complement       = complemented;

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
                                 uint32    minOverlap_) {

  _seqStore        = seqStore_;

  _tig             = NULL;
  _numReads        = 0;

  _sequencesMax   = 0;
  _sequencesLen   = 0;
  _sequences      = NULL;

  _utgpos          = NULL;
  _cnspos          = NULL;

  _minOverlap      = minOverlap_;
  _errorRate       = errorRate_;
  _errorRateMax    = errorRateMax_;
}


unitigConsensus::~unitigConsensus() {

  for (uint32 ss=0; ss<_sequencesLen; ss++)
    delete _sequences[ss];

  delete [] _sequences;
  delete [] _utgpos;
  delete [] _cnspos;
}



void
unitigConsensus::addRead(uint32   readID,
                         uint32   askip, uint32 bskip,
                         bool     complemented,
                         map<uint32, sqRead *>     *inPackageRead) {

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

  _sequences[_sequencesLen++] = new abSequence(readID, seqLen, seq, complemented);

  delete readToDelete;
}



bool
unitigConsensus::initialize(map<uint32, sqRead *>     *reads) {

  if (_numReads == 0) {
    fprintf(stderr, "utgCns::initialize()-- unitig has no children.\n");
    return(false);
  }

  _utgpos = new tgPosition [_numReads];
  _cnspos = new tgPosition [_numReads];

  memcpy(_utgpos, _tig->getChild(0), sizeof(tgPosition) * _numReads);
  memcpy(_cnspos, _tig->getChild(0), sizeof(tgPosition) * _numReads);

  //  Clear the cnspos position.  We use this to show it's been placed by consensus.
  //  Guess the number of columns we'll end up with.
  //  Initialize abacus with the reads.

  for (int32 i=0; i<_numReads; i++) {
    _cnspos[i].setMinMax(0, 0);

    addRead(_utgpos[i].ident(),
            _utgpos[i]._askip, _utgpos[i]._bskip,
            _utgpos[i].isReverse(),
            reads);
  }

  //  Check for duplicate reads

  {
    set<uint32>  dupFrag;

    for (uint32 i=0; i<_numReads; i++) {
      if (_utgpos[i].isRead() == false) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is not a read.\n",
                _tig->tigID(), _utgpos[i].ident());
        return(false);
      }

      if (dupFrag.find(_utgpos[i].ident()) != dupFrag.end()) {
        fprintf(stderr, "unitigConsensus()-- Unitig %d FAILED.  Child %d is a duplicate.\n",
                _tig->tigID(), _utgpos[i].ident());
        return(false);
      }

      dupFrag.insert(_utgpos[i].ident());
    }
  }

  return(true);
}

void unitigConsensus::switchToUncompressedCoordinates(void) {
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
         fprintf(stderr, "switchToUncompressedCoordinates()-- Error: there is a gap in positioning, the last read I have ended at %d and the next starts at %d\n", compressedOffset, _utgpos[child].min());
      assert(_utgpos[child].min() >= compressedOffset);

      uint32 readCompressedPosition = _utgpos[child].min() - compressedOffset;
      uint32 compressedEnd   = _utgpos[child].max();
      uint32 compressedStart = _utgpos[child].min();

      // find the start position in normal read based on position in compressed read
      uint32 i = 0;
      while (i < nlen && (ntoc[i] < readCompressedPosition))
         i++;

      if (showAlgorithm())
          fprintf(stderr, "switchToUncompressedCoordinates()-- I'm trying to find start of child %d compressed %d (dist from guide read is %d and in uncompressed in becomes %d)\n", _utgpos[child].ident(), _utgpos[child].min(), readCompressedPosition, i);

      _utgpos[child].setMinMax(i+uncompressedOffset, i+uncompressedOffset+getSequence(child)->length());
      if (showAlgorithm())
           fprintf(stderr, "switchToUncompressedCoordinates() --Updated read %d which has length %d to be from %d - %d\n", _utgpos[child].ident(), getSequence(child)->length(), _utgpos[child].min(), _utgpos[child].max());

       // update best end if needed
       if (ntoc == NULL || compressedEnd > currentEnd) {
          nlen  = getSequence(child)->length();
          delete[] ntoc;
          ntoc  = new uint32 [ nlen + 1 ];
          uint32 clen  = homopolyCompress(getSequence(child)->getBases(), nlen, NULL, ntoc);

          currentEnd         = compressedEnd;
          compressedOffset   = compressedStart;
          uncompressedOffset = _utgpos[child].min();

          if (showAlgorithm())
             fprintf(stderr, "switchToUncompressedCoordinates()-- Updating guide read to be %d which ends at %d. Best before ended at %d. Now my guide is at %d (%d uncompressed)\n", _utgpos[child].ident(), compressedEnd, currentEnd, compressedOffset, uncompressedOffset);
       }

      if (_utgpos[child].max() > layoutLen)
        layoutLen = _utgpos[child].max();
    }
    delete[] ntoc;
    _tig->_layoutLen = layoutLen;
}

char *
unitigConsensus::generateTemplateStitch(void) {
  int32   minOlap  = _minOverlap;

  //  Initialize, copy the first read.

  uint32       rid      = 0;

  abSequence  *seq      = getSequence(rid);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();

  uint32       tigmax = AS_MAX_READLEN;  //  Must be at least AS_MAX_READLEN, else resizeArray() could fail
  uint32       tiglen = 0;
  char        *tigseq = NULL;

  allocateArray(tigseq, tigmax, resizeArray_clearNew);

  if (showAlgorithm()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "generateTemplateStitch()-- COPY READ read #%d %d (len=%d to %d-%d)\n",
            0, _utgpos[0].ident(), readLen, _utgpos[0].min(), _utgpos[0].max());
  }

  for (uint32 ii=0; ii<readLen; ii++)
    tigseq[tiglen++] = fragment[ii];

  tigseq[tiglen] = 0;

  uint32       ePos = _utgpos[0].max();   //  Expected end of template, from bogart supplied positions.


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
    uint32 nr = 0;  //  Next read
    uint32 nm = 0;  //  Next read maximum position

    //  Pick the next read as the one with the longest extension from all with some minimum overlap
    //  to the template

    if (showAlgorithm())
      fprintf(stderr, "\n");

    for (uint32 ii=rid+1; ii < _numReads; ii++) {

      //  If contained, move to the next read.  (Not terribly useful to log, so we don't)

      if (_utgpos[ii].max() < ePos)
        continue;

      //  If a bigger end position, save the overlap.  One quirk: if we've already saved an overlap, and this
      //  overlap is thin, don't save the thin overlap.

      bool   thick = (_utgpos[ii].min() + minOlap < ePos);
      bool   first = (nm == 0);
      bool   save  = false;

      if ((nm < _utgpos[ii].max()) && (thick || first)) {
        save = true;
        nr   = ii;
        nm   = _utgpos[ii].max();
      }

      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- read #%d/%d ident %d position %d-%d%s%s%s\n",
                ii, _numReads, _utgpos[ii].ident(), _utgpos[ii].min(), _utgpos[ii].max(),
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
    double           bandErrRate   = _errorRate / 4;

    int32            olapLen       = ePos - _utgpos[nr].min();  //  The expected size of the overlap

    if (olapLen < minOlap) {
      if (showAlgorithm())
        fprintf(stderr, "generateTemplateStitch()-- WARNING, increasing min overlap from %d to %u for read %u (%d - %d\n",
                olapLen, readLen, nr, _utgpos[nr].min(), _utgpos[nr].max());
      olapLen = readLen; 
    }

    int32            templateLen  = 0;
    int32            extensionLen = 0;

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
      fprintf(stderr, "generateTemplateStitch()-- ALIGN template %d-%d (len=%d) to read #%d %d %d-%d (len=%d actual=%d at %d-%d)  expecting olap of %d\n",
              tiglen - templateLen, tiglen, templateLen,
              nr, _utgpos[nr].ident(), readBgn, readEnd, readEnd - readBgn, readLen,
              _utgpos[nr].min(), _utgpos[nr].max(),
              olapLen);
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

    //  Reset if the edit distance is waay more than our error rate allows.  This seems to be a quirk with
    //  edlib when aligning to N's - I got startLocation = endLocation = 0 and editDistance = alignmentLength.

    if ((double)result.editDistance / result.alignmentLength > _errorRate) {
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
      fprintf(stderr, "generateTemplateStitch()-- FOUND alignment at %d-%d editDist %d alignLen %d %.f%%\n",
              result.startLocations[0], result.endLocations[0]+1,
              result.editDistance,
              result.alignmentLength,
              100.0 * result.editDistance / result.alignmentLength);

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
        fprintf(stderr, "generateTemplateStitch()-- FAILED to align - no more template to remove!  Fail!\n");
      if (bandErrRate + _errorRate / 4 >= _errorRate)
         tryAgain = false;
      else {
        if (showAlgorithm()) 
          fprintf(stderr, "generateTemplateStitch()-- FAILED to align at %.2f error rate, increasing to %.2f\n", bandErrRate, bandErrRate + _errorRate/4);
        tryAgain = true;
        templateSize  = 0.90;
        extensionSize = 0.10;
        bandErrRate  += _errorRate / 4;
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

      if (showAlgorithm()) {
        fprintf(stderr, "generateTemplateStitch()--\n");
        fprintf(stderr, "generateTemplateStitch()-- Aligned template %d-%d to read %u %d-%d; copy read %d-%d to template.\n",
                tiglen - templateLen, tiglen, nr, readBgn, readEnd, readEnd, readLen);
      }
    } else {
      readBgn = 0;
      readEnd = olapLen;

      if (showAlgorithm()) {
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

    if (showAlgorithm())
      fprintf(stderr, "generateTemplateStitch()-- Append read from %u to %u, starting at tiglen %u\n", readEnd, readLen, tiglen);

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

  return(tigseq);
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
           bool               verbose) {

  EdlibAlignResult align;

  int32   padding        = (int32)ceil(fragmentLength * 0.15);
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

  if (tigend < tigbgn) {
    fprintf(stderr, "alignEdLib()-- WARNING: tigbgn %d > tigend %d - tiglen %d utgpos %d-%d padding %d\n",
            tigbgn, tigend, tiglen, utgpos.min(), utgpos.max(), padding);
    // try to align it to full
    tigbgn = 0;
    tigend = utgpos.max();
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
    aligned        = (alignedErrRate <= errorRate);
    if (verbose)
      fprintf(stderr, " - ALIGNED %.4f at %9d-%-9d\n", alignedErrRate, tigbgn + align.startLocations[0], tigbgn + align.endLocations[0]+1);
  } else {
    if (verbose)
      fprintf(stderr, "\n");
  }

  uint32 MAX_RETRIES=5;

  for (uint32 ii=0; ((ii < MAX_RETRIES) && (aligned == false)); ii++) {
    tigbgn = max((int32)0,      tigbgn - 3 * padding);
    tigend = min((int32)tiglen, tigend + 3 * padding);
    // last attempt make a very wide band
    if (ii == (MAX_RETRIES - 1)) {
       tigbgn = max((int32)0,      tigbgn - 20 * padding);
       tigend = min((int32)tiglen, tigend + 20 * padding);
    }

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



bool
unitigConsensus::initializeGenerate(tgTig                     *tig_,
                                    map<uint32, sqRead *>     *reads_) {

  _tig      = tig_;
  _numReads = _tig->numberOfChildren();

  if (initialize(reads_) == false) {
    fprintf(stderr, "Failed to initialize for tig %u with %u children\n", _tig->tigID(), _tig->numberOfChildren());
    return(false);
  }

  switchToUncompressedCoordinates();

  return(true);
}



bool
unitigConsensus::generatePBDAG(tgTig                     *tig_,
                               char                       aligner_,
                               map<uint32, sqRead *>     *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  //  Build a quick consensus to align to.

  char   *tigseq = generateTemplateStitch();
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

    assert(aligner_ == 'E');  //  Maybe later we'll have more than one aligner again.

    aligned = alignEdLib(aligns[ii],
                         _utgpos[ii],
                         seq->getBases(), seq->length(),
                         tigseq, tiglen,
                         (double)tiglen / _tig->_layoutLen,
                         _errorRate,
                         showAlgorithm());

    if (aligned == false) {
      if (showAlgorithm())
        fprintf(stderr, "generatePBDAG()--    read %7u FAILED\n", _utgpos[ii].ident());

      fail++;

      continue;
    }

    pass++;
  }

  if (showAlgorithm())
    fprintf(stderr, "generatePBDAG()--    read alignment: %d failed, %d passed.\n", fail, pass);

  //  Construct the graph from the alignments.  This is not thread safe.

  if (showAlgorithm())
    fprintf(stderr, "Constructing graph\n");

  AlnGraphBoost ag(string(tigseq, tiglen));

  for (uint32 ii=0; ii<_numReads; ii++) {
    _cnspos[ii].setMinMax(aligns[ii].start, aligns[ii].end);

    if ((aligns[ii].start == 0) &&
        (aligns[ii].end   == 0))
      continue;

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
  std::string cns = ag.consensus(0);

  delete [] tigseq;

  //  Save consensus

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, (uint32) cns.length() + 1, resizeArray_doNothing);

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
unitigConsensus::generateQuick(tgTig                     *tig_,
                               map<uint32, sqRead *>     *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  //  Quick is just the template sequence, so one and done!

  char   *tigseq = generateTemplateStitch();
  uint32  tiglen = strlen(tigseq);

  //  Save consensus

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, tiglen + 1, resizeArray_doNothing);

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
unitigConsensus::generateSingleton(tgTig                     *tig_,
                                   map<uint32, sqRead *>     *reads_) {

  if (initializeGenerate(tig_, reads_) == false)
    return(false);

  assert(_numReads == 1);

  //  Copy the single read to the tig sequence.

  abSequence  *seq      = getSequence(0);
  char        *fragment = seq->getBases();
  uint32       readLen  = seq->length();

  resizeArrayPair(_tig->_bases, _tig->_quals, 0, _tig->_basesMax, readLen + 1, resizeArray_doNothing);

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



//  Align each read to the region of the consensus sequence the read claims
//  to be from, extended by 5% of the read length on either end.  If it fails
//  to align full length, make the extensions larger and/or error rate more
//  tolerant until it aligns.
//
//  Reads that fail to align have their cnspos set to 0.
//
void
unitigConsensus::findCoordinates(void) {

  if (showPlacement()) {
    fprintf(stderr, "\n");
    fprintf(stderr, "TIG %u length %u\n", _tig->tigID(), _tig->length());
    fprintf(stderr, "\n");
  }

  _tig->_childDeltaBitsLen = 1;
  _tig->_childDeltaBits    = new stuffedBits();

  int32         alignShift = 0;

  #pragma omp parallel for schedule(dynamic)
  for (uint32 ii=0; ii<_numReads; ii++) {
    abSequence   *read    = getSequence(ii);
    char         *readSeq = read->getBases();
    uint32        readLen = read->length();

    int32         origbgn = _cnspos[ii].min();
    int32         origend = _cnspos[ii].max();
    int32         origlen = origend - origbgn;

    int32         ext5    = readLen * 0.01;  //  Try an initial alignment nice and tight.  This is
    int32         ext3    = readLen * 0.01;  //  important for the first/last reads so we get them
    double        era     = _errorRate;      //  aligned at the end of the tig.

    if (showPlacement()) {
      fprintf(stderr, "\n");
      fprintf(stderr, "ALIGN read #%d %u length %u cnspos %d %d\n",
              ii, read->seqIdent(), readLen, _cnspos[ii].min(), _cnspos[ii].max());
    }

    _cnspos[ii].setMinMax(0, 0);  //  When the read aligns, we set the true position.

    if ((origbgn == 0) &&
        (origend == 0)) {
      if (showPlacement())
        fprintf(stderr, "skip, failed to align to template initially.\n");
      continue;
    }

    while ((ext5 < readLen * 1.5) &&
           (ext3 < readLen * 1.5) &&
           (era  < 4 * _errorRate)) {
      int32 bgn = max(0, origbgn + alignShift - ext5);                          //  First base in tig we'll align to.
      int32 end = min(   origend + alignShift + ext3, (int32)_tig->length());   //  Last base
      int32 len = end - bgn;   //  WAS: +1

      //  If the length of the region is (wildly) less than expected, reset
      //  bgn to give the alignment a fighting chance.  BPW has only seen
      //  this when running consensus on 10% error reads, and only on reads
      //  near the end of the tig -- that is, when end (above) is reset to
      //  tig->length().

      if ((bgn == 0) &&
          (len < ext5 + origlen + ext3)) {
        int32 endnew = bgn + (origend - origbgn) + ext5 + ext3;

        if (showPlacement())
          fprintf(stderr, "reset position from %d-%d to %d-%d based on len=%d ext5=%d origlen=%d ext3=%d\n",
                  bgn, end,
                  bgn, endnew,
                  len, ext5, origlen, ext3);

        end = (endnew < (int32)_tig->length()) ? endnew : _tig->length();
        len = end - bgn;
      }

      if ((end == _tig->length()) &&
          (len < ext5 + origlen + ext3)) {
        int32 bgnnew = end - (origend - origbgn) - ext5 - ext3;

        if (showPlacement())
          fprintf(stderr, "reset position from %d-%d to %d-%d based on len=%d ext5=%d origlen=%d ext3=%d\n",
                  bgn,    end,
                  bgnnew, end,
                  len, ext5, origlen, ext3);

        bgn = (bgnnew < 0) ? 0 : bgnnew;
        len = end - bgn;
      }

      //  Some logging.

      if (showPlacement())
        fprintf(stderr, "align read #%u length %u to %d-%d - shift %d extension %d %d error rate %.3f\n",
                ii, readLen, bgn, end, alignShift, ext5, ext3, era);

      //  More logging.

      if (showAlignments()) {
        char  save = _tig->bases()[bgn + len];

        _tig->bases()[bgn + len] = 0;

        fprintf(stderr, "READ: %s\n", readSeq);
        fprintf(stderr, "TIG:  %s\n", _tig->bases() + bgn);

        _tig->bases()[bgn + len] = save;
      }

      // Align the entire read into a subsequence of the tig.

      assert(bgn < end);

      EdlibAlignResult align = edlibAlign(readSeq, readLen,
                                          _tig->bases() + bgn, len,
                                          edlibNewAlignConfig(readLen * era * 2, EDLIB_MODE_HW, EDLIB_TASK_PATH));

      //  If nothing aligned, make the tig subsequence bigger and allow more errors.

      if (align.numLocations == 0) {
        ext5 += readLen * 0.05;
        ext3 += readLen * 0.05;
        era  += 0.025;

        if (showPlacement())
          fprintf(stderr, "  NO ALIGNMENT - Increase extension to %d / %d and error rate to %.3f\n", ext5, ext3, era);

        edlibFreeAlignResult(align);

        continue;
      }

      //  Something aligned.  Find how much of the tig subsequence was not used in the alignment.
      //  If we hit either end - the unaligned bit is length 0 - increase that extension and retry.

      int32 unaligned5 = align.startLocations[0];
      int32 unaligned3 = len - (align.endLocations[0] + 1);   //  0-based position of last character in alignment.

      if (showPlacement())
        fprintf(stderr, "               - read %4u original %9u-%9u claimed %9u-%9u aligned %9u-%9u unaligned %d %d\n",
                ii,
                _utgpos[ii].min(), _utgpos[ii].max(),
                origbgn, origend,
                bgn, end,
                unaligned5, unaligned3);

      //  If bump the start, make bigger.

      if ((bgn > 0) && (unaligned5 <= 0)) {
        ext5 += readLen * 0.05;

        if (showPlacement())
          fprintf(stderr, "  BUMPED START - unaligned hangs %d %d - increase 5' extension to %d\n",
                  unaligned5, unaligned3, ext5);

        edlibFreeAlignResult(align);

        continue;
      }

      //  If bump the end, make bigger.

      if ((end < _tig->length()) && (unaligned3 <= 0)) {
        ext3 += readLen * 0.05;

        if (showPlacement())
          fprintf(stderr, "  BUMPED END   - unaligned hangs %d %d - increase 3' extension to %d\n",
                  unaligned5, unaligned3, ext3);

        edlibFreeAlignResult(align);

        continue;
      }

      //  Otherwise, SUCCESS!
      //
      //  Save, for the next alignment, the difference between where the read though it was placed,
      //  and where we actually aligned it to.  This will adjust the next subsequence to (hopefully)
      //  account for any expansion/collapses in the consensus sequence.

      int32   abgn = bgn + align.startLocations[0];
      int32   aend = bgn + align.endLocations[0] + 1;

      alignShift = ((abgn - origbgn) +
                    (aend - origend)) / 2;

      _cnspos[ii].setMinMax(abgn, aend);

      _tig->getChild(ii)->setMinMax(abgn, aend);

      if (showPlacement())
        fprintf(stderr, "  SUCCESS aligned to %d %d\n", abgn, aend);

#pragma omp critical (tgTigLoadAlign)
      {
        _tig->getChild(ii)->_deltaOffset = _tig->_childDeltaBits->getPosition();
        _tig->getChild(ii)->_deltaLen    = edlibAlignmentToCanu(_tig->_childDeltaBits,
                                                                align.alignment,
                                                                align.alignmentLength,
                                                                align.startLocations[0],
                                                                align.endLocations[0] + 1,
                                                                0,
                                                                readLen);
      }

      edlibFreeAlignResult(align);

      break;   //  Stop looping over extension and error rate.
    }  //  Looping over extension and error rate.
  }    //  Looping over reads.
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

  if (showAlgorithm())
    fprintf(stderr, "Circularizing tig %u with expected overlap %u bases.\n", _tig->tigID(), _tig->_circularLength);

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

  bool     found = false;
  sswLib   align(1, -4, -5, -5);

  for (double scale = 1.05; (found == false) && (scale < 2.10); scale += 0.1) {
    uint32 trylen = min(length, (uint32)(_tig->_circularLength * scale));

    if (showPlacement())
      fprintf(stderr, "  Align   A %6u-%6u against B %6u-%6u\n",
              length-trylen, length,
              0, trylen);

    align.align(bases, length, length-trylen, length,           //  A:  the end of the tig
                bases, length, 0,             trylen, false);   //  B:  the start of the tig

    if ((align.errorRate() < 0.05) &&         //  "Found" if the error rate is reasonable, and
        (align.bgnA() > length - trylen) &&   //  A has unaligned stuff at the start, and
        (align.endB() < trylen))              //  B has unaligned stuff at the end.
      found = true;
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
    fprintf(stderr, "  Aligned A %6u-%6u against B %6u-%6u at %.2f%%\n",
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

  //  Align trimmed-off start to end of tig.
  {
    sswLib testBgn;

    if (showAlgorithm())
      fprintf(stderr, "  Test   head %6d-%6d against   tail %6d-%6d\n",
              0, keepBgn, tailBgn-222, tailEnd);

    testBgn.align(bases, length, 0,             keepBgn,
                  bases, length, tailBgn - 222, tailEnd);

    if (showAlgorithm())
      fprintf(stderr, "  Tested head %6d-%6d aligns to tail %6d-%6d\n",
              testBgn.bgnA(), testBgn.endA(), testBgn.bgnB(), testBgn.endB());

    if (testBgn.endB() == keepEnd)
      headTest = true;
  }

  //  Align trimmed-off end to start of tig.
  {
    sswLib testEnd;

    if (showAlgorithm())
      fprintf(stderr, "  Test   tail %6d-%6d against   head %6d-%6d\n",
              keepEnd, length, headBgn, headEnd + 222);

    testEnd.align(bases, length, keepEnd, length,
                  bases, length, headBgn, headEnd + 222);

    if (showAlgorithm())
      fprintf(stderr, "  Tested tail %6d-%6d aligns to head %6d-%6d\n",
              testEnd.bgnA(), testEnd.endA(), testEnd.bgnB(), testEnd.endB());

    if (testEnd.bgnB() == keepBgn)
      tailTest = true;
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
unitigConsensus::generate(tgTig                     *tig_,
                          char                       algorithm_,
                          char                       aligner_,
                          map<uint32, sqRead *>     *reads_) {
  bool  success = false;

  if      (tig_->numberOfChildren() == 1) {
    success = generateSingleton(tig_, reads_);
  }

  else if (algorithm_ == 'Q') {
    success = generateQuick(tig_, reads_);
  }

  else if ((algorithm_ == 'P') ||   //  Normal utgcns.
           (algorithm_ == 'p')) {   //  'norealign'
    success = generatePBDAG(tig_, aligner_, reads_);
  }

  if (success) {
    _tig->_trimBgn = 0;
    _tig->_trimEnd = _tig->length();

    if (algorithm_ == 'P') {
      findCoordinates();
      findRawAlignments();
    }

    trimCircular();
  }

  return(success);
}
