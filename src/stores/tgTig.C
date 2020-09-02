
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

#include "tgTig.H"
#include "sqStore.H"

#include "files.H"
#include "sequence.H"

#include "strings.H"
#include "intervalList.H"



void
tgTig::saveToRecord(tgTigRecord &tr) {
  tr._tigID               = _tigID;

  tr._unused32            = _unused32;
  tr._trimBgn             = _trimBgn;
  tr._trimEnd             = _trimEnd;

  tr._class               = _class;
  tr._suggestRepeat       = _suggestRepeat;
  tr._suggestBubble       = _suggestBubble;
  tr._suggestCircular     = _suggestCircular;
  tr._circularLength      = _circularLength;
  tr._spare               = _spare;

  tr._layoutLen           = _layoutLen;
  tr._basesLen            = _basesLen;
  tr._childrenLen         = _childrenLen;
  tr._childDeltaBitsLen   = _childDeltaBitsLen;
}


void
tgTig::restoreFromRecord(tgTigRecord &tr) {
  _tigID               = tr._tigID;

  _unused32            = tr._unused32;
  _trimBgn             = tr._trimBgn;
  _trimEnd             = tr._trimEnd;

  _class               = tr._class;
  _suggestRepeat       = tr._suggestRepeat;
  _suggestBubble       = tr._suggestBubble;
  _suggestCircular     = tr._suggestCircular;
  _circularLength      = tr._circularLength;
  _spare               = tr._spare;

  _layoutLen           = tr._layoutLen;
  _basesLen            = tr._basesLen;
  _childrenLen         = tr._childrenLen;
  _childDeltaBitsLen   = tr._childDeltaBitsLen;
}



void
tgPosition::initialize(void) {
  _objID       = UINT32_MAX;

  _isRead      = true;   //  Bogus values.
  _isUnitig    = true;
  _isContig    = true;
  _isReverse   = false;

  _spare       = 0;

  _anchor      = UINT32_MAX;
  _ahang       = INT32_MAX;
  _bhang       = INT32_MAX;

  _askip       = 0;
  _bskip       = 0;

  _min         = INT32_MIN;
  _max         = INT32_MAX;

  _deltaOffset = UINT32_MAX;
  _deltaLen    = 0;
}




tgTig::tgTig() {
  _tigID                = UINT32_MAX;

  _unused32             = 0;
  _trimBgn              = 0;
  _trimEnd              = 0;

  _utgcns_verboseLevel  = 0;

  _class                = tgTig_noclass;
  _suggestRepeat        = 0;
  _suggestBubble        = 0;
  _suggestCircular      = 0;
  _circularLength       = 0;
  _spare                = 0;

  _layoutLen            = 0;

  _basesLen             = 0;
  _basesMax             = 0;
  _bases                = NULL;
  _quals                = NULL;

  _children             = NULL;
  _childrenLen          = 0;
  _childrenMax          = 0;

  _stashed              = nullptr;
  _stashedLen           = 0;
  _stashedMax           = 0;

  _childDeltaBitsLen    = 0;
  _childDeltaBits       = NULL;
}

tgTig::~tgTig() {
  delete [] _bases;
  delete [] _quals;
  delete [] _children;

  assert(_stashed == nullptr);

  delete    _childDeltaBits;
}




//  Deep copy the tig.
tgTig &
tgTig::operator=(tgTig & tg) {
  _tigID               = tg._tigID;

  _unused32            = tg._unused32;
  _trimBgn             = tg._trimBgn;
  _trimEnd             = tg._trimEnd;

  _class               = tg._class;
  _suggestRepeat       = tg._suggestRepeat;
  _suggestBubble       = tg._suggestBubble;
  _suggestCircular     = tg._suggestCircular;
  _circularLength      = tg._circularLength;
  _spare               = tg._spare;

  _layoutLen           = tg._layoutLen;

  _basesLen            = tg._basesLen;
  duplicateArray(_bases, _basesLen, _basesMax, tg._bases, tg._basesLen, tg._basesMax);
  duplicateArray(_quals, _basesLen, _basesMax, tg._quals, tg._basesLen, tg._basesMax, true);

  if (_basesLen > 0) {
    assert(_basesMax > _basesLen);
    _bases[_basesLen] = 0;
    _quals[_basesLen] = 0;
  }

  _childrenLen = tg._childrenLen;
  duplicateArray(_children, _childrenLen, _childrenMax, tg._children, tg._childrenLen, tg._childrenMax);

  return(*this);
}



double
tgTig::computeCoverage(void) {
  intervalList<int32>  allL;

  for (uint32 ci=0; ci<numberOfChildren(); ci++) {
    tgPosition *read = getChild(ci);
    uint32      bgn  = read->min();
    uint32      end  = read->max();

    allL.add(bgn, end - bgn);
  }

  intervalDepth<int32>  ID(allL);

  double  aveDepth    = 0;

  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
    aveDepth += (ID.hi(ii) - ID.lo(ii) + 1) * ID.depth(ii);

  if (length() == 0)
    return(0);

  return(aveDepth / length());
}



//  Clears the data but doesn't release memory.  The only way to do that is to delete it.
void
tgTig::clear(void) {
  _tigID                = UINT32_MAX;

  _unused32             = 0;
  _trimBgn              = 0;
  _trimEnd              = 0;

  _class                = tgTig_noclass;
  _suggestRepeat        = 0;
  _suggestBubble        = 0;
  _suggestCircular      = 0;
  _circularLength       = 0;
  _spare                = 0;

  _layoutLen            = 0;
  _basesLen             = 0;
  _childrenLen          = 0;

  _stashed              = nullptr;
  _stashedLen           = 0;
  _stashedMax           = 0;

  delete _childDeltaBits;

  _childDeltaBitsLen    = 0;
  _childDeltaBits       = NULL;
}



bool
tgTig::loadFromStreamOrLayout(FILE *F) {

  //  Decide if the file contains an ASCII layout or a binary stream.  It's probably rather fragile,
  //  testing if the first byte is 't' (from 'tig') or 'T' (from 'TIGR').

  int ch = getc(F);

  ungetc(ch, F);

  if (ch == 't')
    return(loadLayout(F));

  else if (ch == 'T')
    return(loadFromStream(F));

  else
    return(false);
}









void
tgTig::saveToBuffer(writeBuffer *B) {
  char         tag[4] = {'T', 'I', 'G', 'R', };  //  That's tigRecord, not TIGR
  tgTigRecord  tr;

  saveToRecord(tr);

  B->write( tag, 4);
  B->write(&tr,  sizeof(tgTigRecord));

  //  We could save the null byte too, but don't.  It's explicitly added during the load.

  if (_basesLen > 0) {
    B->write(_bases, _basesLen);
    B->write(_quals, _basesLen);
  }

  if (_childrenLen > 0)
    B->write(_children, sizeof(tgPosition) * _childrenLen);

  if (_childDeltaBitsLen > 0)
    _childDeltaBits->dumpToBuffer(B);
}



bool
tgTig::loadFromBuffer(readBuffer *B) {
  char         tag[4];
  tgTigRecord  tr;

  //  Read the tgTigRecord from disk and copy it into our tgTig.

  if (4 != B->read(tag, 4)) {
    fprintf(stderr, "tgTig::loadFromStream()-- failed to read four byte code: %s\n", strerror(errno));
    return(false);
  }

  if ((tag[0] != 'T') ||
      (tag[1] != 'I') ||
      (tag[2] != 'G') ||
      (tag[3] != 'R')) {
    fprintf(stderr, "tgTig::loadFromStream()-- not at a tigRecord, got bytes '%c%c%c%c' (0x%02x%02x%02x%02x).\n",
            tag[0], tag[1], tag[2], tag[3],
            tag[0], tag[1], tag[2], tag[3]);
    return(false);
  }

  if (sizeof(tgTigRecord) != B->read(&tr, sizeof(tgTigRecord))) {
    fprintf(stderr, "tgTig::loadFromStream()-- failed to read tgTigRecord: %s\n", strerror(errno));
    return(false);
  }

  clear();
  restoreFromRecord(tr);

  //  Allocate space for bases/quals and load them.  Be sure to terminate them, too.

  resizeArrayPair(_bases, _quals, 0, _basesMax, _basesLen + 1, resizeArray_doNothing);

  if (_basesLen > 0) {
    B->read(_bases, _basesLen);
    B->read(_quals, _basesLen);

    _bases[_basesLen] = 0;
    _quals[_basesLen] = 0;
  }

  //  Allocate space for reads and alignments, and load them.

  resizeArray(_children,    0, _childrenMax,    _childrenLen,    resizeArray_doNothing);

  if (_childrenLen > 0)
    B->read(_children, sizeof(tgPosition) * _childrenLen);

  if (_childDeltaBitsLen > 0)
    _childDeltaBits = new stuffedBits(B);

  //  Return success.

  return(true);
}



void
tgTig::saveToStream(FILE *F) {
  char         tag[4] = {'T', 'I', 'G', 'R', };  //  That's tigRecord, not TIGR
  tgTigRecord  tr;

  saveToRecord(tr);

  writeToFile(tag, "tgTig::saveToStream::tigr", 4, F);
  writeToFile(tr,  "tgTig::saveToStream::tr",      F);

  //  We could save the null byte too, but don't.  It's explicitly added during the load.

  if (_basesLen > 0) {
    writeToFile(_bases, "tgTig::saveToStream::bases", _basesLen, F);
    writeToFile(_quals, "tgTig::saveToStream::quals", _basesLen, F);
  }

  if (_childrenLen > 0)
    writeToFile(_children, "tgTig::saveToStream::children", _childrenLen, F);

  if (_childDeltaBitsLen > 0)
    _childDeltaBits->dumpToFile(F);
}



bool
tgTig::loadFromStream(FILE *F) {
  char         tag[4];
  tgTigRecord  tr;

  clear();

  //  Read the tgTigRecord from disk and copy it into our tgTig.

  if (4 != loadFromFile(tag, "tgTig::saveToStream::tigr", 4, F, false)) {
    fprintf(stderr, "tgTig::loadFromStream()-- failed to read four byte code: %s\n", strerror(errno));
    return(false);
  }

  if ((tag[0] != 'T') ||
      (tag[1] != 'I') ||
      (tag[2] != 'G') ||
      (tag[3] != 'R')) {
    fprintf(stderr, "tgTig::loadFromStream()-- not at a tigRecord, got bytes '%c%c%c%c' (0x%02x%02x%02x%02x).\n",
            tag[0], tag[1], tag[2], tag[3],
            tag[0], tag[1], tag[2], tag[3]);
    return(false);
  }

  if (0 == loadFromFile(tr, "tgTig::loadFromStream::tr", F, false)) {
    fprintf(stderr, "tgTig::loadFromStream()-- failed to read tgTigRecord: %s\n", strerror(errno));
    return(false);
  }

  restoreFromRecord(tr);

  //  Allocate space for bases/quals and load them.  Be sure to terminate them, too.

  if (_basesLen > 0) {
    resizeArrayPair(_bases, _quals, 0, _basesMax, _basesLen + 1, resizeArray_doNothing);
    loadFromFile(_bases, "tgTig::loadFromStream::bases", _basesLen, F);
    loadFromFile(_quals, "tgTig::loadFromStream::quals", _basesLen, F);

    _bases[_basesLen] = 0;
    _quals[_basesLen] = 0;
  }

  //  Allocate space for reads and alignments, and load them.

  resizeArray(_children,    0, _childrenMax,    _childrenLen,    resizeArray_doNothing);

  if (_childrenLen > 0)
    loadFromFile(_children, "tgTig::savetoStream::children", _childrenLen, F);

  if (_childDeltaBitsLen > 0)
    _childDeltaBits = new stuffedBits(F);

  //  Return success.

  return(true);
};







void
tgTig::dumpLayout(FILE *F, bool withSequence) {
  char  deltaString[128] = {0};
  char  trimString[128]  = {0};

  fprintf(F, "tig " F_U32 "\n", _tigID);
  fprintf(F, "len %d\n",      _layoutLen);

  //  Dump the sequence and quality

  if ((_basesLen == 0) || (withSequence == false)) {
    fputs("cns\n", F);
    fputs("qlt\n", F);

  } else {
    char  *qvString = new char [_basesLen + 1];

    for (uint32 ii=0; ii<_basesLen; ii++)     //  Adjust QV's to Sanger encoding (and make it
      qvString[ii] = _quals[ii] + '!';   //  a character string so we can actually print it).

    qvString[_basesLen] = 0;

    fputs("cns ", F);  fputs(_bases, F);        fputs("\n", F);
    fputs("qlt ", F);  fputs(qvString,     F);  fputs("\n", F);

    delete [] qvString;
  }

  //  Properties.

  //fprintf(F, "unused32        %u\n", _unused32);
  fprintf(F, "trimBgn         %u\n", _trimBgn);
  fprintf(F, "trimEnd         %u\n", _trimEnd);
  fprintf(F, "class           %s\n", toString(_class));
  fprintf(F, "suggestRepeat   %c\n", _suggestRepeat   ? 'T' : 'F');
  fprintf(F, "suggestBubble   %c\n", _suggestBubble   ? 'T' : 'F');
  fprintf(F, "suggestCircular %c\n", _suggestCircular ? 'T' : 'F');
  fprintf(F, "circularLength  %u\n", _circularLength);
  fprintf(F, "numChildren     " F_U32 "\n", _childrenLen);

  //  And the reads.

  for (uint32 i=0; i<_childrenLen; i++) {
    tgPosition *imp = _children + i;

    trimString[0]  = 0;
    deltaString[0] = 0;

    if (imp->_askip + imp->_bskip > 0)
      snprintf(trimString,  128, " trim %6u %6u", imp->_askip, imp->_bskip);

    if (imp->_deltaLen > 0)
      snprintf(deltaString, 128, " delta %5u at %u", imp->_deltaLen, imp->_deltaOffset);


    if (imp->_isRead)
      fprintf(F, "read   %9" F_U32P " anchor %9" F_U32P " hang %7" F_S32P " %7" F_S32P " position %9" F_U32P " %9" F_U32P "%s%s\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), trimString, deltaString);

    if (imp->_isUnitig)
      fprintf(F, "unitig %9" F_U32P " anchor %9" F_U32P " hang %7" F_S32P " %7" F_S32P " position %9" F_U32P " %9" F_U32P "%s%s\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), trimString, deltaString);

    if (imp->_isContig)
      fprintf(F, "contig %9" F_U32P " anchor %9" F_U32P " hang %7" F_S32P " %7" F_S32P " position %9" F_U32P " %9" F_U32P "%s%s\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), trimString, deltaString);
  }

  fprintf(F, "tigend\n");

  //  Fail if it looks weird.

  if (_basesLen > 0) {
    if (_basesLen != _layoutLen) {
      fprintf(stderr, "ERROR: basesLen %u differs from layoutLen %u\n", _basesLen, _layoutLen);
      fprintf(stderr, "ERROR:   _bases length %ld\n", strlen(_bases));
    }
    assert(_basesLen == _layoutLen);
  }
}


bool
tgTig::loadLayout(FILE *F) {
  uint64   LINEnum = 0;
  uint32   LINElen = 0;
  uint32   LINEmax = 1 * 1024 * 1024;
  char    *LINE    = new char [LINEmax];

  uint32   nChildren = 0;

  clear();

  fgets(LINE, LINEmax, F);  LINEnum++;

  if (feof(F)) {
    delete [] LINE;
    return(false);
  }

  while (!feof(F)) {
    splitToWords  W(LINE);

    if        ((W.numWords() == 0) ||
               (W[0][0] == '#') ||
               (W[0][0] == '!')) {
      //  Comment, ignore.

    } else if (strcmp(W[0], "tig") == 0) {
      _tigID = strtouint32(W[1]);

    } else if (strcmp(W[0], "len") == 0) {
      _layoutLen = strtouint32(W[1]);
      resizeArray(LINE, LINElen, LINEmax, _layoutLen + 1, resizeArray_doNothing);

    } else if (((strcmp(W[0], "cns") == 0) || (strcmp(W[0], "qlt") == 0)) && (W.numWords() == 1)) {
      _basesLen = 0;

    } else if (((strcmp(W[0], "cns") == 0) || (strcmp(W[0], "qlt") == 0)) && (W.numWords() == 2)) {
      _basesLen  = strlen(W[1]);
      _layoutLen = _basesLen;    //  Must be enforced, probably should be an explicit error.

      resizeArrayPair(_bases, _quals, 0, _basesMax, _basesLen + 1, resizeArray_doNothing);

      if (W[0][0] == 'c')
        memcpy(_bases, W[1], sizeof(char) * (_basesLen + 1));  //  W[1] is null terminated, and we just copy it in
      else
        memcpy(_quals, W[1], sizeof(char) * (_basesLen + 1));

    //} else if (strcmp(W[0], "unused32") == 0) {
    //  _unused32 = strtouint32(W[1]);
    } else if (strcmp(W[0], "trimBgn") == 0) {
      _trimBgn = strtouint32(W[1]);
    } else if (strcmp(W[0], "trimEnd") == 0) {
      _trimEnd = strtouint32(W[1]);

    } else if (strcmp(W[0], "class") == 0) {
      if      (strcmp(W[1], "unassembled") == 0)
        _class = tgTig_unassembled;
      else if (strcmp(W[1], "contig") == 0)
        _class = tgTig_contig;
      else
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line " F_U64 " invalid: '%s'\n", W[0], LINEnum, LINE), exit(1);

    } else if (strcmp(W[0], "suggestRepeat") == 0) {
      _suggestRepeat = strtouint32(W[1]);

    } else if (strcmp(W[0], "suggestBubble") == 0) {
      _suggestBubble = strtouint32(W[1]);

    } else if (strcmp(W[0], "suggestCircular") == 0) {
      _suggestCircular = strtouint32(W[1]);

    } else if (strcmp(W[0], "circularLength") == 0) {
      _circularLength = strtouint32(W[1]);

    } else if (strcmp(W[0], "numChildren") == 0) {
      //_numChildren = strtouint32(W[1]);
      //resizeArray(_children, 0, _childrenMax, _childrenLen, resizeArray_doNothing);

    } else if ((strcmp(W[0], "read")   == 0) ||
               (strcmp(W[0], "unitig") == 0) ||
               (strcmp(W[0], "contig") == 0)) {

      if (W.numWords() < 10)
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line " F_U64 " invalid: '%s'\n", W[0], LINEnum, LINE), exit(1);

      if (nChildren >= _childrenLen) {
        resizeArray(_children, _childrenLen, _childrenMax, _childrenLen + 1, resizeArray_copyData);
        _childrenLen++;
      }

      _children[nChildren]._objID       = strtouint32(W[1]);
      _children[nChildren]._isRead      = (strcmp(W[0], "read")   == 0);
      _children[nChildren]._isUnitig    = (strcmp(W[0], "unitig") == 0);
      _children[nChildren]._isContig    = (strcmp(W[0], "contig") == 0);
      _children[nChildren]._isReverse   = false;
      _children[nChildren]._spare       = 0;
      _children[nChildren]._anchor      = strtouint32(W[3]);
      _children[nChildren]._ahang       = strtouint32(W[5]);
      _children[nChildren]._bhang       = strtouint32(W[6]);
      _children[nChildren]._askip       = 0;
      _children[nChildren]._bskip       = 0;
      _children[nChildren]._min         = strtouint32(W[8]);
      _children[nChildren]._max         = strtouint32(W[9]);
      _children[nChildren]._deltaOffset = 0;
      _children[nChildren]._deltaLen    = 0;

      if (_children[nChildren]._max < _children[nChildren]._min) {
        _children[nChildren]._min       = strtouint32(W[9]);
        _children[nChildren]._max       = strtouint32(W[8]);
        _children[nChildren]._isReverse = true;
      }


      for (uint32 pos=10; (pos < W.numWords()); pos++) {
        if (strcmp(W[pos], "delta") == 0) {
          _children[nChildren]._deltaLen    = strtouint32(W[++pos]);
          pos++;  //  "at"
          _children[nChildren]._deltaOffset = strtouint32(W[++pos]);
        }

        if (strcmp(W[pos], "trim") == 0) {
          _children[nChildren]._askip = strtouint32(W[++pos]);
          _children[nChildren]._bskip = strtouint32(W[++pos]);
        }

        pos++;
      }

      nChildren++;

    } else if (strcmp(W[0], "tigend") == 0) {
      //  All done, get out of the reading loop.
      break;

    } else {
      //  LINE is probably munged by splitToWords.
      fprintf(stderr, "tgTig::loadLayout()-- unknown line '%s'\n", LINE);
    }

    fgets(LINE, LINEmax, F);  LINEnum++;
  }

  delete [] LINE;

  return(true);
}



//  Dump the tig and all data referenced to a file.
//  For correction, we also need to dump the read this tig is representing.
//
void
tgTig::exportData(writeBuffer  *exportDataFile,
                  sqStore      *seqStore,
                  bool          isForCorrection=false) {
  sqRead           *rd = new sqRead;
  sqReadDataWriter *wr = new sqReadDataWriter;

  //  Export the metadata.

  saveToBuffer(exportDataFile);

  //  Export a read.
  //
  //  This is either the read we're trying to correct, or the first read in
  //  the tig (redundantly stored).

  if (isForCorrection)
    seqStore->sqStore_saveReadToBuffer(exportDataFile, tigID(), rd, wr);

  else
    seqStore->sqStore_saveReadToBuffer(exportDataFile, getChild(0)->ident(), rd, wr);

  //  Now export all the reads in the layot.

  for (uint32 ii=0; ii<numberOfChildren(); ii++)
    seqStore->sqStore_saveReadToBuffer(exportDataFile, getChild(ii)->ident(), rd, wr);

  delete wr;
  delete rd;
}



//  Undo the dump.  This tig is populated with the data from disk,
//  and the reads are loaded into a pair of map<>s.
//
//  Returns true if data was loaded, but minimal checking is done.
//
bool
tgTig::importData(readBuffer                 *importDataFile,
                  map<uint32, sqRead     *>  &reads,
                  FILE                       *layoutOutput,
                  FILE                       *sequenceOutput) {

  //  Try to load the metadata.  If nothing there, we're done.

  if (loadFromBuffer(importDataFile) == false)
    return(false);

  if (layoutOutput)
    dumpLayout(layoutOutput);

  //  We stored numberOfChildren() + 1 reads in the export.  The first read is either a redundant
  //  copy of the first read in the layout, or the read we're trying to correct.  If it's redundant,
  //  we'll just ignore the next copy.

  for (int32 ii=0; ii<numberOfChildren() + 1; ii++) {
    sqRead     *read = new sqRead;

    sqStore::sqStore_loadReadFromBuffer(importDataFile, read);

    if (reads[read->sqRead_readID()] != NULL)     //  If we already have data, just nuke it.  We've
      delete read;                                //  got to read the data from disk regardless.

    else {
      if (sequenceOutput)
        fprintf(sequenceOutput, ">read%u\n%s\b", read->sqRead_readID(), read->sqRead_sequence());

      reads[read->sqRead_readID()] = read;
    }
  }

  return(true);
}



void
tgTig::reverseComplement(void) {

  //  Note that _anchor and the hangs are now invalid.

  ::reverseComplement(_bases, _quals, _basesLen);

  for (uint32 ii=0; ii<_childrenLen; ii++) {
    int32  bgn = _basesLen - _children[ii].bgn();
    int32  end = _basesLen - _children[ii].end();

    _children[ii].set(_children[ii].ident(), 0, 0, 0, bgn, end);
  }

  //  _childDeltas are also invalid.
}


void
tgTig::dumpFASTA(FILE *F) {

  if (consensusExists() == false)
    return;

  AS_UTL_writeFastA(F,
                    bases(), length(), 100,
                    ">tig%08u len=" F_U32 " reads=" F_U32 " class=%s suggestRepeat=%s suggestBubble=%s suggestCircular=%s trim=%u-%u\n",
                    tigID(),
                    length(),
                    numberOfChildren(),
                    toString(_class),
                    _suggestRepeat ? "yes" : "no",
                    _suggestBubble ? "yes" : "no",
                    _suggestCircular ? "yes" : "no",
                    _trimBgn, _trimEnd);
}


void
tgTig::dumpFASTQ(FILE *F) {

  if (consensusExists() == false)
    return;

  AS_UTL_writeFastQ(F,
                    bases(), length(),
                    quals(), length(),
                    "@tig%08u len=" F_U32 " reads=" F_U32 " class=%s suggestRepeat=%s suggestBubble=%s suggestCircular=%s trim=%u-%u\n",
                    tigID(),
                    length(),
                    numberOfChildren(),
                    toString(_class),
                    _suggestRepeat ? "yes" : "no",
                    _suggestBubble ? "yes" : "no",
                    _suggestCircular ? "yes" : "no",
                    _trimBgn, _trimEnd);
}
