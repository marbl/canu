
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2014-DEC-22
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "tgTig.H"

#include "AS_UTL_fileIO.H"
#include "AS_UTL_fasta.C"

#include "splitToWords.H"


tgPosition::tgPosition() {
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


tgTigRecord::tgTigRecord() {
  _tigID           = UINT32_MAX;

  _coverageStat    = 0.0;
  _microhetProb    = 0.0;

  _suggestRepeat   = false;
  _suggestUnique   = false;
  _suggestCircular = false;
  _suggestHaploid  = false;
  _spare           = 0;

  _layoutLen       = 0;
  _gappedLen       = 0;
  _childrenLen     = 0;
  _childDeltasLen  = 0;
}




tgTig::tgTig() {
  _tigID                = UINT32_MAX;

  _coverageStat         = 0;
  _microhetProb         = 0;

  _utgcns_verboseLevel  = 0;

  _suggestRepeat        = 0;
  _suggestUnique        = 0;
  _suggestCircular      = 0;
  _suggestHaploid       = 0;
  _spare                = 0;

  _layoutLen            = 0;

  _gappedBases          = NULL;
  _gappedQuals          = NULL;
  _gappedLen            = 0;
  _gappedMax            = 0;

  _ungappedBases        = NULL;
  _ungappedQuals        = NULL;
  _ungappedLen          = 0;
  _ungappedMax          = 0;

  _children             = NULL;
  _childrenLen          = 0;
  _childrenMax          = 0;

  _childDeltas          = NULL;
  _childDeltasLen       = 0;
  _childDeltasMax       = 0;
}

tgTig::~tgTig() {
  delete [] _gappedBases;
  delete [] _gappedQuals;
  delete [] _ungappedBases;
  delete [] _ungappedQuals;
  delete [] _children;
  delete [] _childDeltas;
}




//  Copy data from an in-core tgTig to an on-disk tgTigRecord.
tgTigRecord &
tgTigRecord::operator=(tgTig & tg) {
  _tigID               = tg._tigID;

  _coverageStat        = tg._coverageStat;
  _microhetProb        = tg._microhetProb;

  _suggestRepeat       = tg._suggestRepeat;
  _suggestUnique       = tg._suggestUnique;
  _suggestCircular     = tg._suggestCircular;
  _suggestHaploid      = tg._suggestHaploid;
  _spare               = tg._spare;

  _layoutLen           = tg._layoutLen;

  _gappedLen           = tg._gappedLen;
  _childrenLen         = tg._childrenLen;
  _childDeltasLen      = tg._childDeltasLen;

  return(*this);
}



//  Copy data from an on-disk tgTigRecord to an in-core tgTig.
tgTig &
tgTig::operator=(tgTigRecord & tr) {
  _tigID               = tr._tigID;

  _coverageStat        = tr._coverageStat;
  _microhetProb        = tr._microhetProb;

  _suggestRepeat       = tr._suggestRepeat;
  _suggestUnique       = tr._suggestUnique;
  _suggestCircular     = tr._suggestCircular;
  _suggestHaploid      = tr._suggestHaploid;
  _spare               = tr._spare;

  _layoutLen           = tr._layoutLen;
  _gappedLen           = tr._gappedLen;
  _childrenLen         = tr._childrenLen;
  _childDeltasLen      = tr._childDeltasLen;

  return(*this);
}





//  Deep copy the tig.
tgTig &
tgTig::operator=(tgTig & tg) {
  _tigID               = tg._tigID;

  _coverageStat        = tg._coverageStat;
  _microhetProb        = tg._microhetProb;

  _suggestRepeat       = tg._suggestRepeat;
  _suggestUnique       = tg._suggestUnique;
  _suggestCircular     = tg._suggestCircular;
  _suggestHaploid      = tg._suggestHaploid;
  _spare               = tg._spare;

  _layoutLen = tg._layoutLen;

  _gappedLen = tg._gappedLen;
  duplicateArray(_gappedBases, _gappedLen, _gappedMax, tg._gappedBases, tg._gappedLen, tg._gappedMax);
  duplicateArray(_gappedQuals, _gappedLen, _gappedMax, tg._gappedQuals, tg._gappedLen, tg._gappedMax, true);

  if (_gappedLen > 0) {
    assert(_gappedMax > _gappedLen);
    _gappedBases[_gappedLen] = 0;
    _gappedQuals[_gappedLen] = 0;
  }

  _ungappedLen = tg._ungappedLen;
  duplicateArray(_ungappedBases, _ungappedLen, _ungappedMax, tg._ungappedBases, tg._ungappedLen, tg._ungappedMax);
  duplicateArray(_ungappedQuals, _ungappedLen, _ungappedMax, tg._ungappedQuals, tg._ungappedLen, tg._ungappedMax, true);

  if (_ungappedLen > 0) {
    assert(_ungappedMax > _ungappedLen);
    _ungappedBases[_ungappedLen] = 0;
    _ungappedQuals[_ungappedLen] = 0;
  }

  _childrenLen = tg._childrenLen;
  duplicateArray(_children, _childrenLen, _childrenMax, tg._children, tg._childrenLen, tg._childrenMax);

  _childDeltasLen = tg._childDeltasLen;
  duplicateArray(_childDeltas, _childDeltasLen, _childDeltasMax, tg._childDeltas, tg._childDeltasLen, tg._childDeltasMax);

  return(*this);
}




void
tgTig::buildUngapped(void) {

  if (_ungappedLen > 0)
    //  Already computed.  Return what is here.
    return;

  if (_gappedLen == 0)
    //  No gapped sequence to convert to ungapped.
    return;

  //  Allocate more space, if needed.  We'll need no more than gappedMax.

  resizeArrayPair(_ungappedBases, _ungappedQuals, 0, _ungappedMax, _gappedMax, resizeArray_doNothing);

  //  gappedLen doesn't include the terminating null, but gappedMax does.
  //  See abMultiAlign.C, among other places.

  if (_gappedLen >= _gappedMax)
    fprintf(stderr, "ERROR: gappedLen = %u >= gappedMax = %u\n",
            _gappedLen+1, _gappedMax);
  assert(_gappedLen < _gappedMax);

  //  Copy all but the gaps.

  _ungappedLen = 0;

  for (uint32 gp=0; gp<_gappedLen; gp++) {
    if (_gappedBases[gp] == '-')
      continue;

    _ungappedBases[_ungappedLen] = _gappedBases[gp];
    _ungappedQuals[_ungappedLen] = _gappedQuals[gp];

    _ungappedLen++;
  }

  assert(_ungappedLen < _ungappedMax);

  //  Terminate it.  Lots of work just for printf.

  _ungappedBases[_ungappedLen] = 0;
  _ungappedQuals[_ungappedLen] = 0;
}





//  Clears the data but doesn't release memory.  The only way to do that is to delete it.
void
tgTig::clear(void) {
  _tigID                = UINT32_MAX;

  _coverageStat         = 0;
  _microhetProb         = 0;

  _suggestRepeat        = 0;
  _suggestUnique        = 0;
  _suggestCircular      = 0;
  _suggestHaploid       = 0;
  _spare                = 0;

  _layoutLen            = 0;
  _gappedLen            = 0;
  _ungappedLen          = 0;
  _childrenLen          = 0;
  _childDeltasLen       = 0;
}



bool
tgTig::loadFromStreamOrLayout(FILE *F) {

  //  Decide if the file contains an ASCII layout or a binary stream.  It's probably rather fragile,
  //  testing if the first byte is 't' (from 'tig') or 'T' (from 'TIGR').

  char ch = getc(F);

  ungetc(ch, F);

  if (ch == 't')
    return(loadLayout(F));

  else if (ch == 'T')
    return(loadFromStream(F));

  else
    return(false);
}



void
tgTig::saveToStream(FILE *F) {
  tgTigRecord  tr = *this;
  char         tag[4] = {'T', 'I', 'G', 'R', };  //  That's tigRecord, not TIGR

  AS_UTL_safeWrite(F,  tag, "tgTig::saveToStream::tigr", sizeof(char), 4);
  AS_UTL_safeWrite(F, &tr,  "tgTig::saveToStream::tr",   sizeof(tgTigRecord), 1);

  //  We could save the null byte too, but don't.  It's explicitly added during the load.

  if (_gappedLen > 0) {
    AS_UTL_safeWrite(F, _gappedBases, "tgTig::saveToStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeWrite(F, _gappedQuals, "tgTig::saveToStream::gappedQuals", sizeof(char), _gappedLen);
  }

  if (_childrenLen > 0)
    AS_UTL_safeWrite(F, _children, "tgTig::saveToStream::children", sizeof(tgPosition), _childrenLen);

  if (_childDeltasLen > 0)
    AS_UTL_safeWrite(F, _childDeltas, "tgTig::saveToStream::childDeltas", sizeof(int32), _childDeltasLen);
}





bool
tgTig::loadFromStream(FILE *F) {
  char    tag[4];

  clear();

  //  Read the tgTigRecord from disk and copy it into our tgTig.

  tgTigRecord  tr;

  if (4 != AS_UTL_safeRead(F, tag, "tgTig::saveToStream::tigr", sizeof(char), 4)) {
    return(false);
  }

  if ((tag[0] != 'T') ||
      (tag[1] != 'I') ||
      (tag[2] != 'G') ||
      (tag[3] != 'R')) {
    return(false);
  }

  if (0 == AS_UTL_safeRead(F, &tr, "tgTig::loadFromStream::tr", sizeof(tgTigRecord), 1)) {
    //  Nothing loaded, end of file.
    return(false);
  }

  *this = tr;

  //  Allocate space for bases/quals and load them.  Be sure to terminate them, too.

  resizeArrayPair(_gappedBases, _gappedQuals, 0, _gappedMax, _gappedLen + 1, resizeArray_doNothing);

  if (_gappedLen > 0) {
    AS_UTL_safeRead(F, _gappedBases, "tgTig::loadFromStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeRead(F, _gappedQuals, "tgTig::loadFromStream::gappedQuals", sizeof(char), _gappedLen);

    _gappedBases[_gappedLen] = 0;
    _gappedQuals[_gappedLen] = 0;
  }

  //  Allocate space for reads and alignments, and load them.

  resizeArray(_children,    0, _childrenMax,    _childrenLen,    resizeArray_doNothing);
  resizeArray(_childDeltas, 0, _childDeltasMax, _childDeltasLen, resizeArray_doNothing);

  if (_childrenLen > 0)
    AS_UTL_safeRead(F, _children, "tgTig::savetoStream::children", sizeof(tgPosition), _childrenLen);

  if (_childDeltasLen > 0)
    AS_UTL_safeRead(F, _childDeltas, "tgTig::loadFromStream::childDeltas", sizeof(int32), _childDeltasLen);

  //  Return success.

  return(true);
};







void
tgTig::dumpLayout(FILE *F) {

  if (_gappedLen > 0)
    assert(_gappedLen == _layoutLen);

  fprintf(F, "tig "F_U32"\n", _tigID);
  fprintf(F, "len %d\n",      _layoutLen);

  if (_gappedLen == 0) {
    fputs("cns\n", F);
    fputs("qlt\n", F);
  } else {
    //fputs("cns ", F);  fwrite(_gappedBases, sizeof(char), _gappedLen, F);  fputs("\n", F);
    //fputs("qlt ", F);  fwrite(_gappedQuals, sizeof(char), _gappedLen, F);  fputs("\n", F);
    fputs("cns ", F);  fputs(_gappedBases, F);  fputs("\n", F);  //  strings are null terminated now, but expected to be long.
    fputs("qlt ", F);  fputs(_gappedQuals, F);  fputs("\n", F);
  }

  fprintf(F, "coverageStat    %f\n", _coverageStat);
  fprintf(F, "microhetProb    %f\n", _microhetProb);
  fprintf(F, "suggestRepeat   %c\n", _suggestRepeat   ? 'T' : 'F');
  fprintf(F, "suggestUnique   %c\n", _suggestUnique   ? 'T' : 'F');
  fprintf(F, "suggestCircular %c\n", _suggestCircular ? 'T' : 'F');
  fprintf(F, "suggestHaploid  %c\n", _suggestHaploid  ? 'T' : 'F');
  fprintf(F, "numChildren     "F_U32"\n", _childrenLen);

  for (uint32 i=0; i<_childrenLen; i++) {
    tgPosition *imp = _children + i;

    if (imp->_isRead)
      fprintf(F, "read   %9"F_U32P" anchor %9"F_U32P" hang %6"F_S32P" %6"F_S32P" position %6"F_U32P" %6"F_U32P" deltas %u at %u\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), imp->_deltaLen, imp->_deltaOffset);

    if (imp->_isUnitig)
      fprintf(F, "unitig %9"F_U32P" anchor %9"F_U32P" hang %6"F_S32P" %6"F_S32P" position %6"F_U32P" %6"F_U32P" deltas %u at %u\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), imp->_deltaLen, imp->_deltaOffset);

    if (imp->_isContig)
      fprintf(F, "contig %9"F_U32P" anchor %9"F_U32P" hang %6"F_S32P" %6"F_S32P" position %6"F_U32P" %6"F_U32P" deltas %u at %u\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end(), imp->_deltaLen, imp->_deltaOffset);
  }

  fprintf(F, "tigend\n");
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

  if (feof(F))
    return(false);

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
      _gappedLen = 0;

    } else if (((strcmp(W[0], "cns") == 0) || (strcmp(W[0], "qlt") == 0)) && (W.numWords() == 2)) {
      _gappedLen = strlen(W[1]);
      _layoutLen = _gappedLen;    //  Must be enforced, probably should be an explicit error.

      resizeArrayPair(_gappedBases, _gappedQuals, 0, _gappedMax, _gappedLen+1, resizeArray_doNothing);

      if (W[0][0] == 'c')
        memcpy(_gappedBases, W[1], sizeof(char) * (_gappedLen + 1));  //  W[1] is null terminated, and we just copy it in
      else
        memcpy(_gappedQuals, W[1], sizeof(char) * (_gappedLen + 1));

    } else if (strcmp(W[0], "coverageStat") == 0) {
      _coverageStat = strtodouble(W[1]);

    } else if (strcmp(W[0], "microhetProb") == 0) {
      _microhetProb = strtodouble(W[1]);

    } else if (strcmp(W[0], "suggestRepeat") == 0) {
      _suggestRepeat = strtouint32(W[1]);

    } else if (strcmp(W[0], "suggestUnique") == 0) {
      _suggestUnique = strtouint32(W[1]);

    } else if (strcmp(W[0], "suggestCircular") == 0) {
      _suggestCircular = strtouint32(W[1]);

    } else if (strcmp(W[0], "suggestHaploid") == 0) {
      _suggestHaploid = strtouint32(W[1]);

    } else if (strcmp(W[0], "numChildren") == 0) {
      _childrenLen = strtouint32(W[1]);
      resizeArray(_children, 0, _childrenMax, _childrenLen, resizeArray_doNothing);

    } else if ((strcmp(W[0], "read")   == 0) ||
               (strcmp(W[0], "unitig") == 0) ||
               (strcmp(W[0], "contig") == 0)) {

      if (W.numWords() < 10)
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line "F_U64" invalid: '%s'\n", W[0], LINEnum, LINE), exit(1);

      if (nChildren >= _childrenLen) {
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line "F_U64" invalid: more children than claimed, adjusting\n", W[0], LINEnum);
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
      _children[nChildren]._min         = strtouint32(W[8]);
      _children[nChildren]._max         = strtouint32(W[9]);
      _children[nChildren]._deltaOffset = 0;
      _children[nChildren]._deltaLen    = 0;

      if (_children[nChildren]._max < _children[nChildren]._min) {
        _children[nChildren]._min       = strtouint32(W[9]);
        _children[nChildren]._max       = strtouint32(W[8]);
        _children[nChildren]._isReverse = true;
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



void
tgTig::dumpFASTQ(FILE *F) {
  AS_UTL_writeFastQ(F,
                    ungappedBases(), ungappedLength(),
                    ungappedQuals(), ungappedLength(), "@utg%08u\n", tigID());
}
