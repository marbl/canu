
#include "tgTig.H"
#include "AS_UTL_fileIO.H"
#include "splitToWords.H"

tgPosition::tgPosition() {
  _objID       = UINT32_MAX;

  _isRead      = true;   //  Bogus values.
  _isUnitig    = true;
  _isContig    = true;

  _spare       = 0;

  _anchor      = UINT32_MAX;
  _ahang       = INT32_MAX;
  _bhang       = INT32_MAX;

  _bgn         = UINT32_MAX;
  _end         = UINT32_MAX;

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
  _ungappedLen     = 0;
  _childrenLen     = 0;
  _childDeltasLen  = 0;
  _variantsLen     = 0;
  _variantDataLen  = 0;
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

  _variants             = NULL;
  _variantsLen          = 0;
  _variantsMax          = 0;

  _variantData          = NULL;
  _variantDataLen       = 0;
  _variantDataMax       = 0;
}

tgTig::~tgTig() {
  delete [] _gappedBases;
  delete [] _gappedQuals;
  delete [] _ungappedBases;
  delete [] _ungappedQuals;
  delete [] _children;
  delete [] _childDeltas;
  delete [] _variants;
  delete [] _variantData;
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
  _ungappedLen         = tg._ungappedLen;
  _childrenLen         = tg._childrenLen;
  _childDeltasLen      = tg._childDeltasLen;
  _variantsLen         = tg._variantsLen;
  _variantDataLen      = tg._variantDataLen;
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
  _ungappedLen         = tr._ungappedLen;
  _childrenLen         = tr._childrenLen;
  _childDeltasLen      = tr._childDeltasLen;
  _variantsLen         = tr._variantsLen;
  _variantDataLen      = tr._variantDataLen;
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

  _variantsLen = tg._variantsLen;
  duplicateArray(_variants, _variantsLen, _variantsMax, tg._variants, tg._variantsLen, tg._variantsMax);

  _variantDataLen = tg._variantDataLen;
  duplicateArray(_variantData, _variantDataLen, _variantDataMax, tg._variantData, tg._variantDataLen, tg._variantDataMax);
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
  _variantsLen          = 0;
};




void
tgTig::saveToStream(FILE *F) {
  tgTigRecord  tr = *this;

  AS_UTL_safeWrite(F, &tr, "tgTig::saveToStream::tr", sizeof(tgTigRecord), 1);

  if (_gappedLen > 0) {
    AS_UTL_safeWrite(F, _gappedBases, "tgTig::saveToStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeWrite(F, _gappedQuals, "tgTig::saveToStream::gappedQuals", sizeof(char), _gappedLen);
  }

  if (_ungappedLen > 0) {
    AS_UTL_safeWrite(F, _ungappedBases, "tgTig::saveToStream::ungappedBases", sizeof(char), _ungappedLen);
    AS_UTL_safeWrite(F, _ungappedQuals, "tgTig::saveToStream::ungappedQuals", sizeof(char), _ungappedLen);
  }

  if (_childrenLen > 0) {
    AS_UTL_safeWrite(F, _children, "tgTig::saveToStream::children", sizeof(char), _childrenLen);
  }

  if (_childDeltasLen > 0) {
    AS_UTL_safeWrite(F, _childDeltas, "tgTig::saveToStream::childDeltas", sizeof(int32), _childDeltasLen);
  }

  if (_variantsLen > 0) {
    AS_UTL_safeWrite(F, _variants, "tgTig::saveToStream::variants", sizeof(char), _variantsLen);
  }
};

void
tgPosition::saveToStream(FILE *F) {
  fprintf(stderr, "tgPosition::saveToStream unimplemented\n");
}

void
tgVariantPosition::saveToStream(FILE *F) {
  fprintf(stderr, "tgVariantPosition::saveToStream unimplemented\n");
}






void
tgTig::loadFromStream(FILE *F) {
  tgTigRecord  tr;

  clear();

  AS_UTL_safeRead(F, &tr, "tgTig::loadFromStream::tr", sizeof(tgTigRecord), 1);

  resizeArrayPair(_gappedBases,   _gappedQuals,   0, _gappedMax,   tr._gappedLen,   resizeArray_doNothing);
  resizeArrayPair(_ungappedBases, _ungappedQuals, 0, _ungappedMax, tr._ungappedLen, resizeArray_doNothing);

  resizeArray(_children, 0, _childrenMax, tr._childrenLen, resizeArray_doNothing);

  resizeArray(_variants, 0, _variantsMax, tr._variantsLen, resizeArray_doNothing);

  if (_gappedLen > 0) {
    AS_UTL_safeRead(F, _gappedBases, "tgTig::loadFromStream::gappedBases", sizeof(char), _gappedLen);
    AS_UTL_safeRead(F, _gappedQuals, "tgTig::loadFromStream::gappedQuals", sizeof(char), _gappedLen);
  }

  if (_ungappedLen > 0) {
    AS_UTL_safeRead(F, _ungappedBases, "tgTig::loadFromStream::ungappedBases", sizeof(char), _ungappedLen);
    AS_UTL_safeRead(F, _ungappedQuals, "tgTig::loadFromStream::ungappedQuals", sizeof(char), _ungappedLen);
  }

  if (_childrenLen > 0) {
    AS_UTL_safeRead(F, _children, "tgTig::savetoStream::children", sizeof(tgPosition), _childrenLen);

    //for (uint32 cc=0; cc<_childrenLen; cc++)
    //  _children[cc].loadFromStream(F);
  }

  if (_variantsLen > 0) {
    AS_UTL_safeRead(F, _variants, "tgTig::loadFromStream::variants", sizeof(tgVariantPosition), _variantsLen);

    //for (uint32 vv=0; vv<_variantsLen; vv++)
    //  _variants[vv].loadFromStream(F);
  }
};

void
tgPosition::loadFromStream(FILE *F) {
  fprintf(stderr, "tgPosition::loadFromStream unimplemented\n");
}

void
tgVariantPosition::loadFromStream(FILE *F) {
  fprintf(stderr, "tgVariantPosition::loadFromStream unimplemented\n");
}








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
      fprintf(F, "read   %9d anchor %9d hang %6d %6d position %6d %6d\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end());

    if (imp->_isUnitig)
      fprintf(F, "unitig %9d anchor %9d hang %6d %6d position %6d %6d\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end());

    if (imp->_isContig)
      fprintf(F, "contig %9d anchor %9d hang %6d %6d position %6d %6d\n",
              imp->ident(), imp->anchor(), imp->aHang(), imp->bHang(), imp->bgn(), imp->end());
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
      _children[nChildren]._spare       = 0;
      _children[nChildren]._anchor      = strtouint32(W[3]);
      _children[nChildren]._ahang       = strtouint32(W[5]);
      _children[nChildren]._bhang       = strtouint32(W[6]);
      _children[nChildren]._bgn         = strtouint32(W[8]);
      _children[nChildren]._end         = strtouint32(W[9]);
      _children[nChildren]._deltaOffset = 0;
      _children[nChildren]._deltaLen    = 0;

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

