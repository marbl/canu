
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

#include "ovStoreHistogram.H"

#include <set>
#include <algorithm>

using namespace std;



ovErateLengthHistogram::~ovErateLengthHistogram() {
  if (_opel)
    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
      delete [] _opel[ii];

  delete [] _opel;
}


ovStoreHistogram::~ovStoreHistogram() {
  delete [] _scoresList;
  delete [] _scores;
}



//  For use in ovStoreDump, computing a length-x-erate histogram.
ovErateLengthHistogram::ovErateLengthHistogram(sqStore *seq) {
  _seq           = seq;
  _maxID         = seq->sqStore_lastReadID();

  _epb           = 1;     //  Evalues per bucket
  _bpb           = 250;   //  Bases per bucket

  _opelLen       = 0;
  _opel          = NULL;
}


//  For use when writing ovStore files.  Data allocated as needed.
ovStoreHistogram::ovStoreHistogram(sqStore *seq) {

  if (seq == NULL)
    fprintf(stderr, "ovStoreHistogram()-- ERROR: I need a valid seqStore.\n"), exit(1);

  _seq           = seq;
  _maxID         = seq->sqStore_lastReadID();

  _scoresListLen = 0;
  _scoresListMax = 0;
  _scoresList    = NULL;
  _scoresListAid = 0;

  _scoresBaseID  = UINT32_MAX;
  _scoresLastID  = 0;
  _scoresAlloc   = 0;
  _scores        = NULL;
}



//  Read only access to existing data.
ovStoreHistogram::ovStoreHistogram(const char *path) {

  _seq           = NULL;
  _maxID         = 0;

  _scoresListLen = 0;
  _scoresListMax = 0;
  _scoresList    = NULL;
  _scoresListAid = 0;

  _scoresBaseID  = UINT32_MAX;
  _scoresLastID  = 0;
  _scoresAlloc   = 0;
  _scores        = NULL;

  char    name[FILENAME_MAX+1];

  createDataName(name, path);

  if (fileExists(name) == false)    //  If no data, nothing to
    return;                         //  load, so leave it empty.

  //  Load!

  FILE *F = AS_UTL_openInputFile(name);

  loadFromFile(_maxID, "ovStoreHistogram::maxID", F);

  loadFromFile(_scoresBaseID, "ovStoreHistogram::scoresBaseID", F);
  loadFromFile(_scoresLastID, "ovStoreHistogram::scoresBaseID", F);

  _scoresAlloc = _scoresLastID - _scoresBaseID + 1;
  _scores      = new oSH_ovlSco [_scoresAlloc];

  loadFromFile(_scores,       "ovStoreHistogram::scores",       _scoresAlloc, F);

  AS_UTL_closeFile(F, name);
}



//  If 'prefix' refers to a directory, the new name will be a file in the directory.
//  Otherwise, it will be an extension to the origianl name.
//
char *
ovStoreHistogram::createDataName(char *name, const char *prefix) {

  if (directoryExists(prefix))
    snprintf(name, FILENAME_MAX, "%s/statistics", prefix);
  else
    snprintf(name, FILENAME_MAX, "%s.statistics", prefix);

  return(name);
}



void
ovStoreHistogram::saveHistogram(char *prefix) {
  char  name[FILENAME_MAX+1];

  //  If no data, don't make any file.

  if (_scores == NULL)
    return;

  //  Otherwise, make an output file.

  createDataName(name, prefix);

  FILE   *F = AS_UTL_openOutputFile(name);

  //  Save all the parameters.

  writeToFile(_maxID, "ovStoreHistogram::maxID", F);

  //  And the data.  There's no apparent guard against getting here with
  //  _scores == NULL, and, if so, we write one element from the NULL
  //  pointer.  Are _scores always set?  Do we just not saveHistogram() when
  //  scores don't exist?  Not sure.

  if (_scores)          //  Process the data for the last read added!
    processScores();

  assert(_scores != NULL);

  writeToFile(_scoresBaseID, "ovStoreHistogram::scoresBaseID", F);
  writeToFile(_scoresLastID, "ovStoreHistogram::scoresLastID", F);
  writeToFile(_scores,       "ovStoreHistogram::scores",       _scoresLastID - _scoresBaseID + 1, F);

  //  That's it!

  AS_UTL_closeFile(F, name);
}



void
ovStoreHistogram::mergeScores(ovStoreHistogram *other) {

  if (other->_scores == NULL)
    return;

  if (_scores == NULL) {
    _maxID         = other->_maxID;

    _scoresBaseID  = 0;
    _scoresLastID  = _maxID;
    _scoresAlloc   = _maxID + 1;

    allocateArray(_scores, _scoresAlloc, resizeArray_clearNew);
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

  assert(_scoresBaseID == 0);  //  Can't copy into a histogram used for counting overlaps.

  memcpy(_scores + other->_scoresBaseID,
         other->_scores,
         sizeof(oSH_ovlSco) * (other->_scoresLastID - other->_scoresBaseID + 1));
}



void
ovStoreHistogram::processScores(uint32 Aid) {
  uint32  scoff = _scoresListAid - _scoresBaseID;

  if (_scoresListLen == 0)       //  If we haven't added any data, there's nothing for us to do.
    return;                      //  This happens when we've just merged in existing data.

  //  Make space for new scores.

  while (scoff >= _scoresAlloc)
    resizeArray(_scores, _scoresAlloc, _scoresAlloc, scoff + 65536, resizeArray_copyData | resizeArray_clearNew);

  //  Sort the scores in decreasing order.

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

  assert(_seq != NULL);                  //  Must have a valid seqStore so we can get read lengths.

  //  Allocate space for the scores data.

  if (_scores == NULL) {
    _scoresListLen = 0;
    _scoresListMax = 16384;
    _scoresListAid = overlap->a_iid;

    allocateArray(_scoresList, _scoresListMax);

    _scoresAlloc = 65535;

    allocateArray(_scores, _scoresAlloc);
  }

  //  And save the overlap, maybe processing the last batch.

  if (_scoresBaseID != UINT32_MAX)                     //  If we've seen an overlap, all remaining overlaps
    assert(_scoresBaseID <= overlap->a_iid);           //  must be larger than the first ID seen.

  _scoresBaseID = min(_scoresBaseID, overlap->a_iid);  //  Save the min/max ID of the overlaps we've seen.
  _scoresLastID = max(_scoresLastID, overlap->a_iid);

  if (_scoresListAid != overlap->a_iid)                //  Process existing overlaps if we
    processScores(overlap->a_iid);                     //  have an overlap for a different ID.

  increaseArray(_scoresList,                           //  Ensure there is space for
                _scoresListLen,                        //  one more overlap.
                _scoresListMax, 32768);

  _scoresList[_scoresListLen++] = overlap->overlapScore();
}



void
ovErateLengthHistogram::addOverlap(ovOverlap *overlap) {

  //  Allocate space for the overlaps-per-evalue-len data.

  if (_opelLen == 0) {
    for (uint32 ii=1; ii<_seq->sqStore_lastReadID(); ii++)
      _opelLen = max(_opelLen, _seq->sqStore_getReadLength(ii));

    _opelLen = _opelLen * 1.40 / _bpb + 1;  //  the overlap could have 40% insertions.
  }

  if (_opel == NULL) {
    allocateArray(_opel, AS_MAX_EVALUE + 1);
  }

  //  Add one to the appropriate entry.

  int32  alen = _seq->sqStore_getReadLength(overlap->a_iid);
  int32  blen = _seq->sqStore_getReadLength(overlap->b_iid);

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
    fprintf(stderr, "BOGUS overlap - id %8u (len %6d) hangs %6" F_OVP " %6" F_OVP " -- id %8u (len %6d) hangs %6" F_OVP " %6" F_OVP "%s\n",
            overlap->a_iid, alen, overlap->dat.ovl.ahg5, overlap->dat.ovl.ahg3,
            overlap->b_iid, blen, overlap->dat.ovl.bhg5, overlap->dat.ovl.bhg3,
            overlap->dat.ovl.flipped ? " flipped" : "");
  }
}



uint32
ovErateLengthHistogram::maxEvalue(void) {
  uint32  maxE = 0;

  for (uint32 ee=0; ee<AS_MAX_EVALUE + 1; ee++) {
    if (_opel[ee] == NULL)
      continue;

    maxE = ee;
  }

  return(maxE);
}



double
ovErateLengthHistogram::maxErate(void) {
  return(AS_OVS_decodeEvalue(maxEvalue()));
}



uint32
ovErateLengthHistogram::maxLength(void) {
  uint32  maxL = 0;

  for (uint32 ee=0; ee<AS_MAX_EVALUE + 1; ee++) {
    if (_opel[ee] == NULL)
      continue;

    for (uint32 ll=maxL; ll<_opelLen; ll++)
      if (_opel[ee][ll] > 0)
        maxL = ll;
  }

  return(maxL * _bpb);
}



void
ovErateLengthHistogram::dumpEvalueLength(FILE *out) {
  uint32  maxE = maxEvalue();
  uint32  maxL = maxLength() / _bpb;

  fprintf(out, "# MAX Evalue  %.4f\n", AS_OVS_decodeEvalue(maxE));
  fprintf(out, "# MAX Length  %u\n",   maxL * _bpb);
  fprintf(out, "$\n");

  for (uint32 ee=0; ee<=maxE; ee++) {
    for (uint32 ll=0; ll<=maxL; ll++)
      fprintf(out, "%u\t%.4f\t%u\n",
              ll * _bpb,
              AS_OVS_decodeEvalue(ee),
              (_opel[ee] == NULL) ? 0 : _opel[ee][ll]);

    fprintf(out, "\n");
  }
}



uint16
ovStoreHistogram::overlapScoreEstimate(uint32 id, uint32 coverage, FILE *scoreDumpFile) {

  if ((id < _scoresBaseID) ||                          //  Return the highest score possible
      (_scoresLastID < id))                            //  if the read is out of range.
    return(UINT16_MAX);                                //  (should never hit this, since we should have all scores when this is used)

  if (coverage == 0)                                   //  Return the highest score if the coverage is zero.
    return(UINT16_MAX);

  id -= _scoresBaseID;                                 //  Offset the id into the array, and check.  (_scoresBaseID should be zero though)

  //  If the coverage requested is within our range, estimate the score.  Otherwise,
  //  return the minimum score - no overlaps will be filtered as there aren't enough.

  double  score = 0.0;

  if (coverage <= _scores[id].points[N_OVL_SCORE-1]) {
    uint32   cp = 1;

    for (; cp<N_OVL_SCORE; cp++)                         //  Search the list of data points for the pair surrounding 'coverage'.
      if (coverage <= _scores[id].points[cp])
        break;

    assert(cp < N_OVL_SCORE);

    assert(_scores[id].points[cp-1] <  coverage);        //  'coverage' is now between cp-1 and cp.
    assert(coverage <= _scores[id].points[cp]);          //  Linearly interpolate to find the score.

    double  x = _scores[id].points[cp] - _scores[id].points[cp-1];
    double  y = _scores[id].scores[cp] - _scores[id].scores[cp-1];

    score     = _scores[id].scores[cp-1] + y / x * (coverage - _scores[id].points[cp-1]);
  }

  if (score < 0)      score = 0;
  if (score > 65535)  score = 65535;

  if (scoreDumpFile != NULL)
    fprintf(scoreDumpFile, "%8u scores %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u %4u/%5u - %f\n",
            id + _scoresBaseID,
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
