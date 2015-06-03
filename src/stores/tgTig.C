
#include "tgTig.H"
#include "AS_UTL_fileIO.H"
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

  _utgcns_verboseLevel  = 0;      //  extern uint32 VERBOSE_MULTIALIGN_OUTPUT;
  _utgcns_smoothWindow  = 11;     //  #define CNS_OPTIONS_MIN_ANCHOR_DEFAULT    11
  _utgcns_splitAlleles  = true;   //  #define CNS_OPTIONS_SPLIT_ALLELES_DEFAULT  1
  _utgcns_doPhasing     = true;   //  #define CNS_OPTIONS_DO_PHASING_DEFAULT     1

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

  _ungappedLen = tg._ungappedLen;
  duplicateArray(_ungappedBases, _ungappedLen, _ungappedMax, tg._ungappedBases, tg._ungappedLen, tg._ungappedMax);
  duplicateArray(_ungappedQuals, _ungappedLen, _ungappedMax, tg._ungappedQuals, tg._ungappedLen, tg._ungappedMax, true);

  _childrenLen = tg._childrenLen;
  duplicateArray(_children, _childrenLen, _childrenMax, tg._children, tg._childrenLen, tg._childrenMax);

  _childDeltasLen = tg._childDeltasLen;
  duplicateArray(_childDeltas, _childDeltasLen, _childDeltasMax, tg._childDeltas, tg._childDeltasLen, tg._childDeltasMax);
}




void
tgTig::buildUngapped(void) {

  if (_ungappedLen > 0)
    //  Already computed.  Return what is here.
    return;

  if (_gappedLen == 0)
    //  No gapped sequence to convert to ungapped.
    return;

  //  Allocate more space, if needed.

  if (_ungappedMax < _gappedMax) {
    delete [] _ungappedBases;
    delete [] _ungappedQuals;

    _ungappedMax = _gappedMax;

    _ungappedBases = new char [_ungappedMax];
    _ungappedQuals = new char [_ungappedMax];
  }

  //  gappedMax could be as small as gappedLen + 1 nul byte.  See abMultiAlign:getConsensus().

  assert(_gappedLen + 1 <= _gappedMax);

  //  Copy all but the gaps.

  _ungappedLen = 0;

  for (uint32 gp=0; gp<_gappedLen; gp++) {
    if (_gappedBases[gp] == '-')
      continue;

    _ungappedBases[_ungappedLen] = _gappedBases[gp];
    _ungappedQuals[_ungappedLen] = _gappedQuals[gp];

    _ungappedLen++;
  }

  assert(_ungappedLen + 1 <= _ungappedMax);

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




void
tgTig::saveToStream(FILE *F) {
  tgTigRecord  tr = *this;

  //fprintf(stderr, "tgTig::saveToStream()-- at "F_U64" - start\n", AS_UTL_ftell(F));

  AS_UTL_safeWrite(F, &tr, "tgTig::saveToStream::tr", sizeof(tgTigRecord), 1);

  if (_gappedLen > 0) {
    AS_UTL_safeWrite(F, _gappedBases, "tgTig::saveToStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeWrite(F, _gappedQuals, "tgTig::saveToStream::gappedQuals", sizeof(char), _gappedLen);
  }

  //fprintf(stderr, "tgTig::saveToStream()-- at "F_U64" - before children, saving "F_U32"\n", AS_UTL_ftell(F), _childrenLen);

  if (_childrenLen > 0) {
    AS_UTL_safeWrite(F, _children, "tgTig::saveToStream::children", sizeof(tgPosition), _childrenLen);
  }

  //fprintf(stderr, "tgTig::saveToStream()-- at "F_U64" - before deltas, saving "F_U32"\n", AS_UTL_ftell(F), _childDeltasLen);

  if (_childDeltasLen > 0) {
    AS_UTL_safeWrite(F, _childDeltas, "tgTig::saveToStream::childDeltas", sizeof(int32), _childDeltasLen);
  }

  //fprintf(stderr, "tgTig::saveToStream()-- at "F_U64" - finsihed\n", AS_UTL_ftell(F));
}





void
tgTig::loadFromStream(FILE *F) {
  tgTigRecord  tr;

  clear();

  //fprintf(stderr, "tgTig::loadFromStream()-- Loading at position "F_U64" - start\n", AS_UTL_ftell(F));

  AS_UTL_safeRead(F, &tr, "tgTig::loadFromStream::tr", sizeof(tgTigRecord), 1);

  *this = tr;

  //  After that copy, the various Len fields are bigger than the Max fields.  Quick!  Allocate arrays!

  resizeArrayPair(_gappedBases,   _gappedQuals,   0, _gappedMax,   _gappedLen,   resizeArray_doNothing);
  //resizeArrayPair(_ungappedBases, _ungappedQuals, 0, _ungappedMax, _ungappedLen, resizeArray_doNothing);

  resizeArray(_children, 0, _childrenMax, _childrenLen, resizeArray_doNothing);

  resizeArray(_childDeltas, 0, _childDeltasMax, _childDeltasLen, resizeArray_doNothing);

  if (_gappedLen > 0) {
    //fprintf(stderr, "tgTig::loadFromStream()-- loading %u gapped bases\n", _gappedLen);
    AS_UTL_safeRead(F, _gappedBases, "tgTig::loadFromStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeRead(F, _gappedQuals, "tgTig::loadFromStream::gappedQuals", sizeof(char), _gappedLen);
  }

  //if (_ungappedLen > 0) {
  //  //fprintf(stderr, "tgTig::loadFromStream()-- loading %u ungapped bases\n", _ungappedLen);
  //  AS_UTL_safeRead(F, _ungappedBases, "tgTig::loadFromStream::ungappedBases", sizeof(char), _ungappedLen);
  //  AS_UTL_safeRead(F, _ungappedQuals, "tgTig::loadFromStream::ungappedQuals", sizeof(char), _ungappedLen);
  //}

  //fprintf(stderr, "tgTig::loadFromStream()-- Loading at position "F_U64" - before children, need to load "F_U32" at "F_U64" bytes each.\n",
  //        AS_UTL_ftell(F), _childrenLen, sizeof(tgPosition));

  if (_childrenLen > 0) {
    //fprintf(stderr, "tgTig::loadFromStream()-- loading %u children\n", _childrenLen);
    AS_UTL_safeRead(F, _children, "tgTig::savetoStream::children", sizeof(tgPosition), _childrenLen);
  }

  //fprintf(stderr, "tgTig::loadFromStream()-- Loading at position "F_U64" - before deltas, need to load "F_U32" at "F_U64" bytes each.\n",
  //        AS_UTL_ftell(F), _childDeltasLen, sizeof(tgPosition));

  if (_childDeltasLen > 0) {
    AS_UTL_safeRead(F, _childDeltas, "tgTig::loadFromStream::childDeltas", sizeof(int32), _childDeltasLen);
  }

  //fprintf(stderr, "tgTig::loadFromStream()-- Loading at position "F_U64" - finished\n", AS_UTL_ftell(F));
};







void
tgTig::dumpLayout(FILE *F) {

  fprintf(F, "tig "F_U32"\n", _tigID);
  fprintf(F, "len %d\n",      _layoutLen);
  fprintf(F, "cns %s\n",     (_gappedLen > 0) ? _gappedBases : "");
  fprintf(F, "qlt %s\n",     (_gappedLen > 0) ? _gappedQuals : "");
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

    } else if (strcmp(W[0], "cns") == 0) {
      if (W.numWords() > 1) {
        resizeArrayPair(_gappedBases, _gappedQuals, _gappedLen, _gappedMax, _layoutLen, resizeArray_doNothing);
        _gappedLen = _layoutLen;
        memcpy(_gappedBases, W[1], sizeof(char) * (_gappedLen + 1));
      }

    } else if (strcmp(W[0], "qlt") == 0) {
      if (W.numWords() > 1) {
        memcpy(_gappedQuals, W[1], sizeof(char) * (_gappedLen + 1));
      }

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

      if (W.numWords() != 10)
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line "F_U64" invalid: '%s'\n", W[0], LINEnum, LINE), exit(1);

      if (nChildren >= _childrenLen) {
        fprintf(stderr, "tgTig::loadLayout()-- '%s' line "F_U64" invalid: more children than claimed, adjusting\n", W[0], LINEnum, LINE);
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

