
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
 *    Brian P. Walenz beginning on 2016-OCT-25
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStoreHistogram.H"

#include <algorithm>

using namespace std;


void
ovStoreHistogram::initialize(gkStore *gkp) {
  _gkp = gkp;

  _maxOlength = 0;
  _maxEvalue  = 0;

  _epb = 1;     //  Evalues per bucket
  _bpb = 250;   //  Bases per bucket

  _opelLen = 0;
  _opel    = NULL;

  _oprLen  = 0;
  _oprMax  = 0;
  _opr     = NULL;

  _scoresListLen = 0;
  _scoresListMax = 0;
  _scoresList    = NULL;
  _scoresListAid = 0;

  _scoresBgn     = 0;
  _scoresLen     = 0;
  _scoresMax     = 0;
  _scores        = NULL;
}



void
ovStoreHistogram::clear(void) {

  if (_opel)
    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
      delete [] _opel[ii];

  delete [] _opel;
  delete [] _opr;
  delete [] _scoresList;
  delete [] _scores;

  initialize(_gkp);
}



ovStoreHistogram::ovStoreHistogram() {
  initialize();
}



ovStoreHistogram::ovStoreHistogram(char *path) {
  initialize();
  loadData(path);
}



ovStoreHistogram::ovStoreHistogram(gkStore *gkp, ovFileType type) {
  initialize(gkp);

  //  When writing full overlaps out of an overlapper (ovFileFullWrite) we want
  //  to keep track of the number of overlaps per read.  We could pre-allocate
  //  the array based on the size of gkpStore, but if we don't have that, it's
  //  easy enough to grow the array.
  //
  //  _opr is notably skipped if ovFileFullWriteNoCounts is used.  That symbol
  //  isn't actually used anywhere except in this comment (and when some ovFile
  //  is created) so we mention it here for grep.

  if (type == ovFileFullWrite) {
    _oprLen = 0;
    _oprMax = (_gkp == NULL) ? (256 * 1024) : (_gkp->gkStore_getNumReads() + 1);
    _opr    = new uint32   [_oprMax];

    memset(_opr, 0, sizeof(uint32) * _oprMax);
  }

  //  When writing store overlaps (ovFileNormalWrite) we want to keep track of
  //  how many overlaps for each evalue X length.
  //
  //  A gkpStore is required here so we can allocate the correct amount of
  //  space and compute the length of an overlap.
  //
  //  The histogram always allocates one pointer for each eValue (there's only 4096 of them),
  //  but defers allocating the vector until needed.

  if (type == ovFileNormalWrite) {
    if (_gkp == NULL)
      fprintf(stderr, "ovStoreHistogram()-- ERROR: I need a valid gkpStore.\n"), exit(1);

    for (uint32 ii=1; ii<_gkp->gkStore_getNumReads(); ii++)
      if (_opelLen < _gkp->gkStore_getRead(ii)->gkRead_sequenceLength())
        _opelLen = _gkp->gkStore_getRead(ii)->gkRead_sequenceLength();
    _opelLen = _opelLen * 1.40 / _bpb + 1;  //  the overlap could have 40% insertions.

    _opel = new uint32 * [AS_MAX_EVALUE + 1];

    memset(_opel, 0, sizeof(uint32 *) * (AS_MAX_EVALUE + 1));
  }

  //  When writing store overlaps (ovFileNormalWrite) _for_correction_ (which we don't know here) we
  //  want to keep a profile of the overlap scores for later filtering.  We don't technically need
  //  to allocate stuff here, but if we don't, we never collect these stats because _scores isn't
  //  allocated.  Oh, the quandry!

  if (type == ovFileNormalWrite) {
    if (_gkp == NULL)
      fprintf(stderr, "ovStoreHistogram()-- ERROR: I need a valid gkpStore.\n"), exit(1);

    _scoresListLen = 0;
    _scoresListMax = 16384;                          //  Enough for 16k overlaps per read.
    _scoresList    = new uint16 [_scoresListMax];
    _scoresListAid = 0;

    _scoresBgn     = 0;
    _scoresLen     = 0;
    _scoresMax     = 65535;                          //  Enough for 64k reads.
    _scores        = new oSH_ovlSco [_scoresMax];
  }
}



ovStoreHistogram::~ovStoreHistogram() {
  clear();
}



void
ovStoreHistogram::processScores(uint32 Aid) {

  if (_scoresBgn == 0) {    //  Bootstrap.  On the first overlap, _scoresListAid == 0 and we call this function
    _scoresListAid = Aid;   //  (because the overlap can't have a_iid == 0).  We just need to initialize by
    _scoresBgn     = Aid;   //  setting _scoresListAid (the overlaps being loaded) and _scoreBgn (the first
    return;                 //  stored in these stats) to Aid.
  }

  if (_scoresListLen == 0)  //  Likewise, if we haven't added any data, there's nothing for us to do.
    return;                 //  This happens when we just add data, e.g., the master stats for a store.

  assert(Aid > _scoresListAid);  //  Next overlap (Aid) must be strictly larger than what we've loaded.

  uint32  scoff = _scoresListAid - _scoresBgn;   //  Index into _scores array for this read.

  if (scoff >= _scoresMax)
    resizeArray(_scores, _scoresLen, _scoresMax, scoff + 65536, resizeArray_copyData | resizeArray_clearNew);

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#endif
  sort(_scoresList, _scoresList + _scoresListLen, greater<uint16>());


  //  Decide on a set of points to save.  Eventually, maybe, we'll analyze the graph to find inflection points.
  //  For now, we just sample evenly-ish.

  double  step  = (double)_scoresListLen / N_OVL_SCORE;
  double  point = 0;

  if (step < 1)
    step = 1;

  if (step > 10)   //  With current N_OVL_SCORE=16, this gets us to 150x coverage.
    step = 10;

  //  Then just save the points.  We first fill the array with the last point in case there are fewer
  //  than space for.

  for (uint32 ii=0; ii<N_OVL_SCORE; ii++)
    _scores[scoff].points[ii] = _scoresListLen-1;

  for (uint32 ii=0; ii<N_OVL_SCORE-1 && point<_scoresListLen; ii++, point += step)
    _scores[scoff].points[ii] = (uint16)floor(point + 0.5);

  //  Now just fill in the scores.

  for (uint32 ii=0; ii<N_OVL_SCORE; ii++)
    _scores[scoff].scores[ii] = _scoresList[ _scores[scoff].points[ii] ];


  _scoresListLen = 0;           //  Reset the scores list for the overlaps for the next read.
  _scoresListAid = Aid;         //  Remember that next read ID.

  _scoresLen = scoff + 1;       //  Remember the last valid data in _scores.
}



void
ovStoreHistogram::addOverlap(ovOverlap *overlap) {

  //  For overlaps out of overlapper, track the number of overlaps per read.

  if (_opr) {
    uint32   maxID = max(overlap->a_iid, overlap->b_iid);

    if (_oprMax < maxID)
      resizeArray(_opr, _oprLen, _oprMax, maxID + maxID/2, resizeArray_copyData | resizeArray_clearNew);

    if (_oprLen < maxID + 1)
      _oprLen = maxID + 1;

    _opr[overlap->a_iid]++;
    _opr[overlap->b_iid]++;
  }

  //  For overlaps in the store, track the number of overlaps per evalue-length

  if (_opel) {
    uint32 ev  = overlap->evalue();
    uint32 len = (_gkp->gkStore_getRead(overlap->a_iid)->gkRead_sequenceLength() - overlap->dat.ovl.ahg5 - overlap->dat.ovl.ahg3 +
                  _gkp->gkStore_getRead(overlap->b_iid)->gkRead_sequenceLength() - overlap->dat.ovl.bhg5 - overlap->dat.ovl.bhg3) / 2;

    if (_maxEvalue  < ev)    _maxEvalue  = ev;
    if (_maxOlength < len)   _maxOlength = len;

    ev  /= _epb;
    len /= _bpb;

    if (_opel[ev] == NULL) {
      _opel[ev] = new uint32 [_opelLen];
      memset(_opel[ev], 0, sizeof(uint32) * _opelLen);
    }

    int32  alen = _gkp->gkStore_getRead(overlap->a_iid)->gkRead_sequenceLength();
    int32  blen = _gkp->gkStore_getRead(overlap->b_iid)->gkRead_sequenceLength();

    if (len < _opelLen) {
      //fprintf(stderr, "overlap %8u (len %6d) %8u (len %6d) hangs %6" F_U64P " %6d %6" F_U64P " - %6" F_U64P " %6d %6" F_U64P " flip " F_U64 "\n",
      //        overlap->a_iid, alen,
      //        overlap->b_iid, blen,
      //        overlap->dat.ovl.ahg5, (int32)alen - (int32)overlap->dat.ovl.ahg5 - (int32)overlap->dat.ovl.ahg3, overlap->dat.ovl.ahg3,
      //        overlap->dat.ovl.bhg5, (int32)blen - (int32)overlap->dat.ovl.bhg5 - (int32)overlap->dat.ovl.bhg3, overlap->dat.ovl.bhg3,
      //        overlap->dat.ovl.flipped);
      _opel[ev][len]++;
    } else {
      fprintf(stderr, "overlap %8u (len %6d) %8u (len %6d) hangs %6" F_U64P " %6d %6" F_U64P " - %6" F_U64P " %6d %6" F_U64P " flip " F_U64 " -- BOGUS\n",
              overlap->a_iid, alen,
              overlap->b_iid, blen,
              overlap->dat.ovl.ahg5, (int32)alen - (int32)overlap->dat.ovl.ahg5 - (int32)overlap->dat.ovl.ahg3, overlap->dat.ovl.ahg3,
              overlap->dat.ovl.bhg5, (int32)blen - (int32)overlap->dat.ovl.bhg5 - (int32)overlap->dat.ovl.bhg3, overlap->dat.ovl.bhg3,
              overlap->dat.ovl.flipped);
    }
  }

  //  For overlap scoring, process an existing scoresList if the ID changed, then and add the new overlap to the scoresList.

  if (_scores) {
    if (_scoresListAid != overlap->a_iid)
      processScores(overlap->a_iid);

    if (_scoresListLen >= _scoresListMax)
      resizeArray(_scoresList, _scoresListLen, _scoresListMax, _scoresListMax + 32768);

    _scoresList[_scoresListLen++] = overlap->overlapScore();

    //fprintf(stderr, "ADD OLAP #%u for aid %u  a %u b %u - score %u - aid\n",
    //        _scoresListLen-1, _scoresListAid, overlap->a_iid, overlap->b_iid, _scoresList[_scoresListLen-1]);
  }
}




//  Build an output file name from a prefix and a suffix based
//  on if the prefix is a directory or a file.  If a directory,
//  the new name will be a file in the directory, otherwise,
//  it will be an extension to the origianl name.
void
createDataName(char *name, char *prefix, char *suffix) {

  if (AS_UTL_fileExists(prefix, true, false)) {
    snprintf(name, FILENAME_MAX, "%s/%s", prefix, suffix);
  }

  else {
    AS_UTL_findBaseFileName(name, prefix);
    strcat(name, ".");
    strcat(name, suffix);
  }
}



void
ovStoreHistogram::saveData(char *prefix) {
  char  name[FILENAME_MAX];

  //  If we have overlaps-per-read data, dump it.  Just a simple array.

  if (_opr) {
    createDataName(name, prefix, "counts");

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to open counts file '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &_oprLen, "ovStoreHistogram::nr",  sizeof(uint32), 1);
    AS_UTL_safeWrite(F,  _opr,    "ovStoreHistogram::opr", sizeof(uint32), _oprLen);

    fclose(F);
  }

  //  If we have overlaps-per-evalue-length, dump it.  This is a bit more complicated, as it has
  //  holes in the array.

  if (_opel) {
    createDataName(name, prefix, "evalueLen");

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to open evalueLen file '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    uint32  nArr = 0;

    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
      if (_opel[ii])
        nArr++;

    AS_UTL_safeWrite(F, &_opelLen,    "ovStoreHistogram::opelLen",    sizeof(uint32), 1);
    AS_UTL_safeWrite(F, &_maxOlength, "ovStoreHistogram::maxOlength", sizeof(uint32), 1);
    AS_UTL_safeWrite(F, &_maxEvalue,  "ovStoreHistogram::maxEvalue",  sizeof(uint32), 1);

    AS_UTL_safeWrite(F, &_epb, "ovStoreHistogram::epb", sizeof(uint32), 1);
    AS_UTL_safeWrite(F, &_bpb, "ovStoreHistogram::bpb", sizeof(uint32), 1);

    AS_UTL_safeWrite(F, &nArr, "ovStoreHistogram::nArr", sizeof(uint32), 1);

    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++) {
      if (_opel[ii] == NULL)
        continue;

      AS_UTL_safeWrite(F, &ii,       "ovStoreHistogram::evalue",    sizeof(uint32),  1);
      AS_UTL_safeWrite(F, _opel[ii], "ovStoreHistogram::evalueLen", sizeof(uint32), _opelLen);
    }

    fclose(F);
  }

  //  If we have overlap scores, process the last one and then dump it.

  if (_scores) {
    createDataName(name, prefix, "overlapScores");

    processScores();

    errno = 0;
    FILE *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "failed to open overlapScores file '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &_scoresBgn, "ovStoreHistogram::scoresBgn",   sizeof(uint32), 1);
    AS_UTL_safeWrite(F, &_scoresLen, "ovStoreHistogram::scoresLen",   sizeof(uint32), 1);
    AS_UTL_safeWrite(F,  _scores,    "ovStoreHistogram::scores",      sizeof(oSH_ovlSco), _scoresLen);

    fclose(F);
  }
}



void
ovStoreHistogram::loadData(char *prefix, uint32 maxIID) {
  char    name[FILENAME_MAX];

  //  Add in any overlaps-per-read data.

  createDataName(name, prefix, "counts");

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "failed to open counts file '%s' for reading: %s\n", name, strerror(errno)), exit(1);

    uint32  inLen = 0;

    AS_UTL_safeRead(F, &inLen, "ovStoreHistogram::opr", sizeof(uint32), 1);         //  How many values on disk?

    if (_oprMax < inLen)                                                            //  Resize to fit those values
      resizeArray(_opr, _oprLen, _oprMax, inLen + inLen/2, resizeArray_copyData | resizeArray_clearNew);

    if (maxIID < inLen)
      fprintf(stderr, "WARNING: histogram file '%s' has data for %u reads, but only %u reads known.\n", name, inLen, maxIID);

    if (_oprLen < inLen)                                                            //  Remember the new length
      _oprLen = inLen;

    uint32 *in = new uint32 [inLen];                                                //  Allocate temp space for new values

    AS_UTL_safeRead(F, in, "ovStoreHistogram::opr", sizeof(uint32), inLen);         //  Load new values

    for (uint32 ii=0; ii<inLen; ii++)                                               //  Add in new values
      _opr[ii] += in[ii];

    delete [] in;

    fclose(F);
  }

  //  Add in any overlaps-per-evalue-length data.

  createDataName(name, prefix, "evalueLen");

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "failed to open evalueLen file '%s' for reading: %s\n", name, strerror(errno)), exit(1);

    uint32  nArr = 0;

    AS_UTL_safeRead(F, &_opelLen,    "ovStoreHistogram::opelLen",    sizeof(uint32), 1);  //  Load parameters of the data
    AS_UTL_safeRead(F, &_maxOlength, "ovStoreHistogram::maxOlength", sizeof(uint32), 1);
    AS_UTL_safeRead(F, &_maxEvalue,  "ovStoreHistogram::maxEvalue",  sizeof(uint32), 1);
    AS_UTL_safeRead(F, &_epb,        "ovStoreHistogram::epb",        sizeof(uint32), 1);
    AS_UTL_safeRead(F, &_bpb,        "ovStoreHistogram::bpb",        sizeof(uint32), 1);
    AS_UTL_safeRead(F, &nArr,        "ovStoreHistogram::nArr",       sizeof(uint32), 1);

    if (_opel == NULL)
      allocateArray(_opel, AS_MAX_EVALUE+1, resizeArray_clearNew);                        //  Abuse resizeArray() to alloc new

    uint32 *in = new uint32 [_opelLen];                                                   //  Allocate space for a single vector

    for (uint32 ev=0; nArr-- > 0; ) {                                                     //  For each saved vector:
      AS_UTL_safeRead(F, &ev, "ovStoreHistogram::evalue", sizeof(uint32), 1);             //    Load the evalue it is for
      AS_UTL_safeRead(F, in, "ovStoreHistogram::evalueLen", sizeof(uint32), _opelLen);    //    Load the data.

      if (_opel[ev] == NULL)                                                              //    More abuse, if needed
        allocateArray(_opel[ev], _opelLen, resizeArray_clearNew);

      for (uint32 kk=0; kk<_opelLen; kk++)                                                //    Add new data to old data
        _opel[ev][kk] += in[kk];
    }

    delete [] in;

    fclose(F);
  }

  //  Add in any overlap score data.

  createDataName(name, prefix, "overlapScores");

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "failed to open evalueLen file '%s' for reading: %s\n", name, strerror(errno)), exit(1);

    uint32  nArr = 0;

    AS_UTL_safeRead(F, &_scoresBgn,  "ovStoreHistogram::scoresBgn",  sizeof(uint32), 1);
    AS_UTL_safeRead(F, &_scoresLen,  "ovStoreHistogram::scoresLen",  sizeof(uint32), 1);

    _scoresMax = _scoresLen;
    _scores    = new oSH_ovlSco [_scoresMax];

    AS_UTL_safeRead(F,  _scores,     "ovStoreHistogram::scores",     sizeof(oSH_ovlSco), _scoresLen);
  }
}



void
ovStoreHistogram::removeData(char *prefix) {
  char    name[FILENAME_MAX];

  createDataName(name, prefix, "counts");          AS_UTL_unlink(name);
  createDataName(name, prefix, "evalueLen");       AS_UTL_unlink(name);
  createDataName(name, prefix, "overlapScores");   AS_UTL_unlink(name);
}



void
ovStoreHistogram::add(ovStoreHistogram *input) {

  //  Add in any overlaps-per-read data.

  if (input->_opr) {
    resizeArray(_opr, _oprLen, _oprMax, input->_oprMax, resizeArray_copyData | resizeArray_clearNew);

    for (uint32 ii=0; ii<input->_oprMax; ii++)
      _opr[ii] += input->_opr[ii];

    _oprLen = max(_oprLen, input->_oprLen);
  }

  //  Add in any overlaps-per-evalue-length data.

  if (input->_opel) {
    if (_opel == NULL) {
      allocateArray(_opel, AS_MAX_EVALUE+1, resizeArray_clearNew);
      _opelLen    = input->_opelLen;
      _maxOlength = input->_maxOlength;
      _maxEvalue  = input->_maxEvalue;
      _epb        = input->_epb;
      _bpb        = input->_bpb;
    }

    if ((_opelLen != input->_opelLen) ||
        (_epb     != input->_epb) ||
        (_bpb     != input->_bpb)) {
      fprintf(stderr, "ERROR: can't merge histogram; parameters differ.\n");
      fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _opelLen, input->_opelLen);
      fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _epb,     input->_epb);
      fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _bpb,     input->_bpb);
      exit(1);
    }

    _maxOlength = max(_maxOlength, input->_maxOlength);
    _maxEvalue  = max(_maxEvalue,  input->_maxEvalue);

    for (uint32 ev=0; ev<AS_MAX_EVALUE+1; ev++) {
      if (input->_opel[ev] == NULL)
        continue;

      if (_opel[ev] == NULL)
        allocateArray(_opel[ev], _opelLen, resizeArray_clearNew);

      for (uint32 kk=0; kk<_opelLen; kk++)
        _opel[ev][kk] += input->_opel[ev][kk];
    }
  }

  //  Add in any overlap score data.

  if (input->_scores) {
    uint32   oBgn = _scoresBgn;
    uint32   oEnd = _scoresBgn + _scoresLen;
    uint32   oLen = _scoresLen;

    input->processScores();  //  Make sure the input is all up-to-date.

    uint32   iBgn = input->_scoresBgn;
    uint32   iEnd = input->_scoresBgn + input->_scoresLen;
    uint32   iLen = input->_scoresLen;

    //  Copy new scores to an empty allocation
    if (oLen == 0) {
      oSH_ovlSco *sccopy = new oSH_ovlSco [iLen];

      memcpy(sccopy, input->_scores, sizeof(oSH_ovlSco) * iLen);

      _scores    = sccopy;
      _scoresBgn = iBgn;
      _scoresMax = iLen;
      _scoresLen = iLen;
    }

    //  Copy new scores to middle of (existing) scores.
    else if ((oBgn <= iBgn) &&
             (iEnd <= oEnd)) {
      memcpy(_scores + iBgn - oBgn, input->_scores, sizeof(oSH_ovlSco) * iLen);
    }

    //  Copy new scores to end of (reallocated) scores.
    else if (oEnd < iEnd) {
      oSH_ovlSco *sccopy = new oSH_ovlSco [iEnd - oBgn];

      memset(sccopy, 0, sizeof(oSH_ovlSco) * (iEnd - oBgn));

      memcpy(sccopy,                      _scores, sizeof(oSH_ovlSco) * oLen);
      memcpy(sccopy + iBgn - oBgn, input->_scores, sizeof(oSH_ovlSco) * iLen);

      delete [] _scores;
      _scores    = sccopy;
      _scoresBgn = oBgn;
      _scoresMax = iEnd - oBgn;
      _scoresLen = iEnd - oBgn;
    }

    //  Copy new scores to start of (reallocated) scores.
    else if (iBgn < oBgn) {
      oSH_ovlSco *sccopy = new oSH_ovlSco [oEnd - iBgn];

      memset(sccopy, 0, sizeof(oSH_ovlSco) * (oEnd - iBgn));

      memcpy(sccopy + oBgn - iBgn,        _scores, sizeof(oSH_ovlSco) * oLen);
      memcpy(sccopy,               input->_scores, sizeof(oSH_ovlSco) * iLen);

      delete [] _scores;
      _scores    = sccopy;
      _scoresBgn = iBgn;
      _scoresMax = oEnd - iBgn;
      _scoresLen = oEnd - iBgn;
    }
  }
}



uint64
ovStoreHistogram::getOverlapsPerRead(uint32 *oprOut, uint32 oprOutLen) {
  uint64 tot = 0;

  if (oprOutLen < _oprLen)
    fprintf(stderr, "ERROR: more reads in histogram than available for output?  oprOutLen=%u  _oprLen=%u\n", oprOutLen, _oprLen), exit(1);

  for (uint32 ii=0; ii<_oprLen; ii++) {
    oprOut[ii] += _opr[ii];
    tot        += _opr[ii];
  }

  return(tot);
}



void
ovStoreHistogram::dumpEvalueLength(FILE *out) {
  uint32  maxEvalue = 0;
  uint32  maxLength = 0;

  //  Find the largest Evalue and length with values in the histogram

  for (uint32 ee=0; ee<AS_MAX_EVALUE + 1; ee++) {
    if (_opel[ee] == NULL)
      continue;

    maxEvalue = ee;

    for (uint32 ll=maxLength; ll<_opelLen; ll++)
      if (_opel[ee][ll] > 0)
        maxLength = ll;
  }

  //  Dump those values

  for (uint32 ee=0; ee<=maxEvalue; ee++) {
    for (uint32 ll=0; ll<=maxLength; ll++)
      fprintf(out, "%u\t%.4f\t%u\n",
              ll * _bpb,
              AS_OVS_decodeEvalue(ee),
              (_opel[ee] == NULL) ? 0 : _opel[ee][ll]);

    fprintf(out, "\n");
  }

  fprintf(stderr, "MAX Evalue  %.4f\n", AS_OVS_decodeEvalue(maxEvalue));
  fprintf(stderr, "MAX Length  %u\n",   maxLength * _bpb);
}



uint16
ovStoreHistogram::overlapScoreEstimate(uint32 id, uint32 coverage) {

  if ((id < _scoresBgn) ||                             //  Return the highest score
      (_scoresBgn + _scoresLen <= id))                 //  if the read is out of range.
    return(UINT16_MAX);

  id -= _scoresBgn;                                    //  Offset the id into the array, and check.

  if (coverage == 0)                                   //  Return the highest score if the coverage is zero.
    return(UINT16_MAX);

  if (_scores[id].points[N_OVL_SCORE-1] < coverage)    //  Return the lowest score if the coverage is higher
    return(0);                                         //  than the number of overlaps.

  uint32   cp = 1;

  for (; cp<N_OVL_SCORE; cp++)                         //  Search the list of data points for the pair surrounding 'coverage'.
    if (coverage <= _scores[id].points[cp])
      break;

  assert(_scores[id].points[cp-1] <  coverage);        //  'coverage' is now between cp-1 and cp.
  assert(coverage <= _scores[id].points[cp]);          //  Linearly interpolate to find the score.

  double  x     = _scores[id].points[cp] - _scores[id].points[cp-1];
  double  y     = _scores[id].scores[cp] - _scores[id].scores[cp-1];          
  double  score = _scores[id].scores[cp-1] + y / x * (coverage - _scores[id].points[cp-1]);

  if (score < 0)      score = 0;
  if (score > 65535)  score = 65535;

#if 0
  fprintf(stdout, "%8u scores %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u - %f\n",
          id + _scoresBgn,
          _scores[id].points[0], _scores[id].scores[0],
          _scores[id].points[1], _scores[id].scores[1],
          _scores[id].points[2], _scores[id].scores[2],
          _scores[id].points[3], _scores[id].scores[3],
          _scores[id].points[4], _scores[id].scores[4],
          _scores[id].points[5], _scores[id].scores[5],
          _scores[id].points[6], _scores[id].scores[6],
          _scores[id].points[7], _scores[id].scores[7],
          _scores[id].points[8], _scores[id].scores[8],
          _scores[id].points[9], _scores[id].scores[9],
          _scores[id].points[10], _scores[id].scores[10],
          _scores[id].points[11], _scores[id].scores[11],
          _scores[id].points[12], _scores[id].scores[12],
          _scores[id].points[13], _scores[id].scores[13],
          _scores[id].points[14], _scores[id].scores[14],
          _scores[id].points[15], _scores[id].scores[15],
          score);
#endif

  return((uint16)floor(score));
}
