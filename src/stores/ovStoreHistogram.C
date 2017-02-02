
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



ovStoreHistogram::ovStoreHistogram() {

  _gkp = NULL;

  _maxOlength = 0;
  _maxEvalue  = 0;

  _epb = 0;
  _bpb = 0;

  _opelLen = 0;
  _opel    = NULL;

  _oprLen  = 0;
  _oprMax  = 0;
  _opr     = NULL;
}


ovStoreHistogram::ovStoreHistogram(char *path) {

  _gkp = NULL;

  _maxOlength = 0;
  _maxEvalue  = 0;

  _epb = 0;
  _bpb = 0;

  _opelLen = 0;
  _opel    = NULL;

  _oprLen  = 0;
  _oprMax  = 0;
  _opr     = NULL;

  loadData(path);
}


ovStoreHistogram::ovStoreHistogram(gkStore *gkp, ovFileType type) {

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
}



ovStoreHistogram::~ovStoreHistogram() {

  if (_opel)
    for (uint32 ii=0; ii<AS_MAX_EVALUE + 1; ii++)
      delete [] _opel[ii];

  delete [] _opel;
  delete [] _opr;
}



void
ovStoreHistogram::addOverlap(ovOverlap *overlap) {

  if (_opr) {
    uint32   maxID = max(overlap->a_iid, overlap->b_iid);

    if (_oprMax < maxID)
      resizeArray(_opr, _oprLen, _oprMax, maxID + maxID/2, resizeArray_copyData | resizeArray_clearNew);

    if (_oprLen < maxID + 1)
      _oprLen = maxID + 1;

    _opr[overlap->a_iid]++;
    _opr[overlap->b_iid]++;
  }

  if (_opel) {
    uint32 ev  = overlap->evalue();
    uint32 len = (_gkp->gkStore_getRead(overlap->a_iid)->gkRead_sequenceLength() - overlap->dat.ovl.ahg5 - overlap->dat.ovl.ahg3 +
                  _gkp->gkStore_getRead(overlap->b_iid)->gkRead_sequenceLength() - overlap->dat.ovl.bhg5 - overlap->dat.ovl.bhg3) / 2;

    ev  /= _epb;
    len /= _bpb;

    if (_opel[ev] == NULL) {
      _opel[ev] = new uint32 [_opelLen];
      memset(_opel[ev], 0, sizeof(uint32) * _opelLen);
    }

    int32  alen = _gkp->gkStore_getRead(overlap->a_iid)->gkRead_sequenceLength();
    int32  blen = _gkp->gkStore_getRead(overlap->b_iid)->gkRead_sequenceLength();

    if (len < _opelLen) {
      //fprintf(stderr, "overlap %8u (len %6d) %8u (len %6d) hangs %6lu %6d %6lu - %6lu %6d %6lu flip %lu\n",
      //        overlap->a_iid, alen,
      //        overlap->b_iid, blen,
      //        overlap->dat.ovl.ahg5, (int32)alen - (int32)overlap->dat.ovl.ahg5 - (int32)overlap->dat.ovl.ahg3, overlap->dat.ovl.ahg3,
      //        overlap->dat.ovl.bhg5, (int32)blen - (int32)overlap->dat.ovl.bhg5 - (int32)overlap->dat.ovl.bhg3, overlap->dat.ovl.bhg3,
      //        overlap->dat.ovl.flipped);
      _opel[ev][len]++;
    } else {
      fprintf(stderr, "overlap %8u (len %6d) %8u (len %6d) hangs %6lu %6d %6lu - %6lu %6d %6lu flip %lu -- BOGUS\n",
              overlap->a_iid, alen,
              overlap->b_iid, blen,
              overlap->dat.ovl.ahg5, (int32)alen - (int32)overlap->dat.ovl.ahg5 - (int32)overlap->dat.ovl.ahg3, overlap->dat.ovl.ahg3,
              overlap->dat.ovl.bhg5, (int32)blen - (int32)overlap->dat.ovl.bhg5 - (int32)overlap->dat.ovl.bhg3, overlap->dat.ovl.bhg3,
              overlap->dat.ovl.flipped);
    }
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
}



void
ovStoreHistogram::removeData(char *prefix) {
  char    name[FILENAME_MAX];

  createDataName(name, prefix, "counts");      AS_UTL_unlink(name);
  createDataName(name, prefix, "evalueLen");   AS_UTL_unlink(name);
}



void
ovStoreHistogram::add(ovStoreHistogram *input) {

  if (input->_opr) {
    resizeArray(_opr, _oprLen, _oprMax, input->_oprMax, resizeArray_copyData | resizeArray_clearNew);

    for (uint32 ii=0; ii<input->_oprMax; ii++)
      _opr[ii] += input->_opr[ii];

    _oprLen = max(_oprLen, input->_oprLen);
  }

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
