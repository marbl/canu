
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
 *    Sergey Koren beginning on 2017-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "ovStoreHistogram.H"

#include <set>
#include <algorithm>

using namespace std;



void
ovStoreHistogram::clear(gkStore *gkp) {

  _gkp           = gkp;

  _bgnID         = UINT32_MAX;
  _endID         = 0;
  _maxID         = (gkp == NULL) ? 0 : gkp->gkStore_getNumReads();

  _epb           = 1;     //  Evalues per bucket
  _bpb           = 250;   //  Bases per bucket

  _opelLen       = 0;
  _opel          = NULL;

  _scoresListLen = 0;
  _scoresListMax = 0;
  _scoresList    = NULL;
  _scoresListAid = 0;

  _scoresLen     = 0;
  _scores        = NULL;
};



ovStoreHistogram::~ovStoreHistogram() {

  if (_opel)
    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
      delete [] _opel[ii];

  delete [] _opel;
  delete [] _scoresList;
  delete [] _scores;
}



//  For use in ovStoreIndexer (ovStoreSliceWriter::mergeHistogram()), allocate
//  enough space for all the data, so we can stop reallocating over and over.
ovStoreHistogram::ovStoreHistogram() {
  clear();
}



//  For use when writing ovStore files.  Data allocated as needed.
ovStoreHistogram::ovStoreHistogram(gkStore *gkp) {

  clear(gkp);

  if (_gkp == NULL)
    fprintf(stderr, "ovStoreHistogram()-- ERROR: I need a valid gkpStore.\n"), exit(1);
}



//  Read only access to existing data.
ovStoreHistogram::ovStoreHistogram(const char *path) {

  clear();

  char    name[FILENAME_MAX+1];
  createDataName(name, path);

  if (AS_UTL_fileExists(name, false, false) == false)    //  If no data, nothing to
    return;                                              //  load, so leave it empty.

  //  Load!

  FILE *F = AS_UTL_openInputFile(name);

  AS_UTL_safeRead(F, &_bgnID, "ovStoreHistogram::bgnID", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &_endID, "ovStoreHistogram::endID", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &_maxID, "ovStoreHistogram::maxID", sizeof(uint32), 1);

  //  Data for overlapsPerEvalueLength

  uint32  nArr = 0;
  uint32 *aArr = new uint32 [AS_MAX_EVALUE + 1];

  AS_UTL_safeRead(F, &_opelLen, "ovStoreHistogram::opelLen", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &_epb,     "ovStoreHistogram::epb",     sizeof(uint32), 1);
  AS_UTL_safeRead(F, &_bpb,     "ovStoreHistogram::bpb",     sizeof(uint32), 1);
  AS_UTL_safeRead(F, &nArr,     "ovStoreHistogram::nArr",    sizeof(uint32), 1);
  AS_UTL_safeRead(F,  aArr,     "ovStoreHistogram::nArr",    sizeof(uint32), nArr);

  allocateArray(_opel, AS_MAX_EVALUE+1, resizeArray_clearNew);

  for (uint32 ii=0; ii<nArr; ii++)
    _opel[aArr[ii]] = new uint32 [_opelLen];

  for (uint32 ii=0; ii<nArr; ii++ )
    AS_UTL_safeRead(F, _opel[aArr[ii]], "ovStoreHistogram::evalueLen", sizeof(uint32), _opelLen);

  delete [] aArr;

  //  Data for overlapsScores

  _scoresLen = _endID - _bgnID + 1;
  _scores    = new oSH_ovlSco [_scoresLen];

  AS_UTL_safeRead(F,  _scores,    "ovStoreHistogram::scores",     sizeof(oSH_ovlSco), _scoresLen);

  AS_UTL_closeFile(F);
}



//  If 'prefix' refers to a directory, the new name will be a file in the directory.
//  Otherwise, it will be an extension to the origianl name.
//
char *
ovStoreHistogram::createDataName(char *name, const char *prefix) {

  if (AS_UTL_fileExists(prefix, true, false)) {              //  If name is a directory,
    snprintf(name, FILENAME_MAX, "%s/statistics", prefix);   //  return a file in the directory.
  }

  else {                                                     //  Otherwise, replace the suffix
    AS_UTL_findBaseFileName(name, prefix);                   //  with '.statistics'.
    strcat(name, ".statistics");
  }

  return(name);
}



void
ovStoreHistogram::saveHistogram(char *prefix) {
  char  name[FILENAME_MAX+1];

  //  If no data, don't make any file.

  if ((_opel   == NULL) &&
      (_scores == NULL))
    return;

  //  Otherwise, make an output file.

  createDataName(name, prefix);

  FILE   *F = AS_UTL_openOutputFile(name);

  fprintf(stderr, "SAVE HISTOGRAM DATA to '%s' for ID %u to %u\n", name, _bgnID, _endID);

  //  Save all the parameters.

  AS_UTL_safeWrite(F, &_bgnID, "ovStoreHistogram::bgnID", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &_endID, "ovStoreHistogram::endID", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &_maxID, "ovStoreHistogram::maxID", sizeof(uint32), 1);

  //  Data for overlapsPerEvalueLength

  uint32  nArr = 0;
  uint32 *aArr = new uint32 [AS_MAX_EVALUE + 1];

  for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
    if ((_opel != NULL) && (_opel[ii] != NULL))
      aArr[nArr++] = ii;

  AS_UTL_safeWrite(F, &_opelLen, "ovStoreHistogram::opelLen", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &_epb,     "ovStoreHistogram::epb",     sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &_bpb,     "ovStoreHistogram::bpb",     sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &nArr,     "ovStoreHistogram::nArr",    sizeof(uint32), 1);
  AS_UTL_safeWrite(F,  aArr,     "ovStoreHistogram::nArr",    sizeof(uint32), nArr);

  for (uint32 ii=0; ii<nArr; ii++)
    AS_UTL_safeWrite(F, _opel[aArr[ii]], "ovStoreHistogram::evalueLen", sizeof(uint32), _opelLen);

  delete [] aArr;

  //  Data for overlapsScores

  if (_scores)          //  Process the data for the last read added!
    processScores();

  AS_UTL_safeWrite(F,  _scores, "ovStoreHistogram::scores",      sizeof(oSH_ovlSco), _endID - _bgnID + 1);

  //  That's it!

  AS_UTL_closeFile(F, name);
}



void
ovStoreHistogram::mergeOPEL(ovStoreHistogram *other) {

  if (other->_opel == NULL)
    return;

  if (_opel == NULL) {
    _epb        = other->_epb;
    _bpb        = other->_bpb;
    _opelLen    = other->_opelLen;
    allocateArray(_opel, AS_MAX_EVALUE+1, resizeArray_clearNew);
  }

  if ((_epb     != other->_epb) ||
      (_bpb     != other->_bpb) ||
      (_opelLen != other->_opelLen)) {
    fprintf(stderr, "ERROR: can't merge histogram; parameters differ.\n");
    fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _opelLen, other->_opelLen);
    fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _epb,     other->_epb);
    fprintf(stderr, "ERROR:   opelLen = %7u vs %7u\n", _bpb,     other->_bpb);
    exit(1);
  }

  for (uint32 ev=0; ev<AS_MAX_EVALUE+1; ev++) {
    if (other->_opel[ev] == NULL)
      continue;

    if (_opel[ev] == NULL)
      allocateArray(_opel[ev], _opelLen, resizeArray_clearNew);

    for (uint32 kk=0; kk<_opelLen; kk++)
      _opel[ev][kk] += other->_opel[ev][kk];
  }
}



void
ovStoreHistogram::mergeScores(ovStoreHistogram *other) {

  if (other->_scores == NULL)
    return;

  if (_scores == NULL) {
    _maxID         = other->_maxID;
    _scoresLen     = _maxID + 1;
    allocateArray(_scores, _scoresLen, resizeArray_clearNew);
  }

  if (_maxID != other->_maxID) {
    fprintf(stderr, "ERROR: can't merge histogram; parameters differ.\n");
    fprintf(stderr, "ERROR:   maxID = %9u vs %9u\n", _maxID, other->_maxID);
    exit(1);
  }

  //  Make sure all the data in 'other' is processed.  This occurs in the sequential store
  //  build when a file gets full; usually this last processScores() is handled in the
  //  destructor, just before the data is dumped to disk, but we need to force it here.

  other->processScores();

  //  Copy other scores to our array.  No checking of overwriting data is performed.

  fprintf(stderr, "COPYING scores for ID %u to %u\n", other->_bgnID, other->_endID);

  memcpy(_scores + other->_bgnID,
         other->_scores,
         sizeof(oSH_ovlSco) * (other->_endID - other->_bgnID + 1));
}



void
ovStoreHistogram::processScores(uint32 Aid) {
  uint32  scoff = _scoresListAid - _bgnID;

  if (_scoresListLen == 0)       //  If we haven't added any data, there's nothing for us to do.
    return;                      //  This happens when we've just merged in existing data.

  //  Make space for new scores.

  while (scoff >= _scoresLen)
    resizeArray(_scores, _scoresLen, _scoresLen, scoff + 65536, resizeArray_copyData | resizeArray_clearNew);

  //  Sort the scores in decreasing order.

#ifdef _GLIBCXX_PARALLEL
  __gnu_sequential::
#endif
  sort(_scoresList, _scoresList + _scoresListLen, greater<uint16>());

  //  Decide on a set of points to save.  Eventually, maybe, we'll analyze the plot
  //  to find inflection points.  For now, we just sample evenly-ish.

  double  step  = (double)_scoresListLen / N_OVL_SCORE;
  double  point = 0;

  if (step < 1)
    step = 1;

  if (step > 10)   //  With current N_OVL_SCORE=16, this gets us to 150x coverage.
    step = 10;

  //  Then just save the points.  We first fill the array with the
  //  last point in case there are fewer than space for.

  for (uint32 ii=0; ii<N_OVL_SCORE; ii++)                                   //  Initialize all points with
    _scores[scoff].points[ii] = _scoresListLen-1;                           //  the last data value.

  for (uint32 ii=0; ((ii<N_OVL_SCORE-1) &&                                  //  Fill valid data.
                     (point<_scoresListLen)); ii++, point += step)
    _scores[scoff].points[ii] = (uint16)floor(point + 0.5);

  for (uint32 ii=0; ii<N_OVL_SCORE; ii++)                                   //  And add the scores.
    _scores[scoff].scores[ii] = _scoresList[ _scores[scoff].points[ii] ];

  //  Reset for the next overlap.  The next overlap must be larger than what we just processed.

  assert(Aid > _scoresListAid);

  _scoresListLen = 0;
  _scoresListAid = Aid;
}



void
ovStoreHistogram::addOverlap(ovOverlap *overlap) {

  assert(_gkp != NULL);                  //  Must have a valid gkpStore so we can get read lengths.

  if (_bgnID != UINT32_MAX)              //  If we've seen an overlap, all remaining overlaps
    assert(_bgnID <= overlap->a_iid);    //  must be larger than the first ID seen.

  _bgnID = min(_bgnID, overlap->a_iid);  //  Save the min/max ID of the overlaps we've seen.
  _endID = max(_endID, overlap->a_iid);

  //

  if (_opelLen == 0) {
    for (uint32 ii=1; ii<_gkp->gkStore_getNumReads(); ii++)
      _opelLen = max(_opelLen, _gkp->gkStore_getRead(ii)->gkRead_sequenceLength());

    _opelLen = _opelLen * 1.40 / _bpb + 1;  //  the overlap could have 40% insertions.
  }

  if (_opel == NULL) {
    allocateArray(_opel, AS_MAX_EVALUE + 1);
  }

  //  For overlaps in the store, track the number of overlaps per evalue-length

  int32  alen = _gkp->gkStore_getRead(overlap->a_iid)->gkRead_sequenceLength();
  int32  blen = _gkp->gkStore_getRead(overlap->b_iid)->gkRead_sequenceLength();

  uint32 ev   = overlap->evalue();
  uint32 len  = (alen - overlap->dat.ovl.ahg5 - overlap->dat.ovl.ahg3 +
                 blen - overlap->dat.ovl.bhg5 - overlap->dat.ovl.bhg3) / 2;

  ev  /= _epb;
  len /= _bpb;

  if (_opel[ev] == NULL) {
    _opel[ev] = new uint32 [_opelLen];
    memset(_opel[ev], 0, sizeof(uint32) * _opelLen);
  }

  if (len < _opelLen) {
    _opel[ev][len]++;
  }

  else {
    fprintf(stderr, "overlap %8u (len %6d) %8u (len %6d) hangs %6" F_OVP " %6d %6" F_OVP " - %6" F_OVP " %6d %6" F_OVP " flip " F_OV " -- BOGUS\n",
            overlap->a_iid, alen,
            overlap->b_iid, blen,
            overlap->dat.ovl.ahg5, alen - (int32)overlap->dat.ovl.ahg5 - (int32)overlap->dat.ovl.ahg3, overlap->dat.ovl.ahg3,
            overlap->dat.ovl.bhg5, blen - (int32)overlap->dat.ovl.bhg5 - (int32)overlap->dat.ovl.bhg3, overlap->dat.ovl.bhg3,
            overlap->dat.ovl.flipped);
  }

  //

  if (_scores == NULL) {
    _scoresListLen = 0;
    _scoresListMax = 16384;
    _scoresListAid = overlap->a_iid;

    allocateArray(_scoresList, _scoresListMax);

    _scoresLen     = 65535;

    allocateArray(_scores, _scoresLen);
  }

  if (_scoresListAid != overlap->a_iid)   //  Process existing overlaps if we
    processScores(overlap->a_iid);        //  have an overlap for a different ID.

  increaseArray(_scoresList,              //  Ensure there is space for
                _scoresListLen,           //  one more overlap.
                _scoresListMax, 32768);

  _scoresList[_scoresListLen++] = overlap->overlapScore();
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
ovStoreHistogram::overlapScoreEstimate(uint32 id, uint32 coverage, bool showLog) {

  if ((id < _bgnID) ||                                 //  Return the highest score possible
      (_endID < id))                                   //  if the read is out of range.
    return(UINT16_MAX);

  id -= _bgnID    ;                                    //  Offset the id into the array, and check.

  if (coverage == 0)                                   //  Return the highest score if the coverage is zero.
    return(UINT16_MAX);

  if (_scores[id].points[N_OVL_SCORE-1] < coverage)    //  Return the lowest score if the coverage is higher than the
    return(0);                                         //  number of overlaps.  This also handles the no overlaps case.

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

  if (showLog == true)
    fprintf(stdout, "%8u scores %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u %3u/%5u - %f\n",
            id + _bgnID,
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

  return((uint16)floor(score));
}
