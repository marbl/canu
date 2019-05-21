
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
 *  This file is derived from:
 *
 *    src/sequence/sequence.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2019APR-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

class sampledDistribution {
public:
  sampledDistribution() {
    _dataLen = 0;
    _dataMax = 0;
    _data    = NULL;

    _dataSum = 0;
  };

  ~sampledDistribution() {
    delete [] _data;
  };


  void     loadDistribution(char *path) {

    if ((path == NULL) || (path[0] == 0))
      return;

    _dataLen = 0;
    _dataMax = 1048576;
    _data    = new uint32 [_dataMax];

    memset(_data, 0, sizeof(uint32) * _dataMax);

    splitToWords  S;

    uint32  Lnum = 0;
    uint32  Llen = 0;
    uint32  Lmax = 1024;
    char   *L    = new char [Lmax];

    FILE   *D    = AS_UTL_openInputFile(path);

    while (AS_UTL_readLine(L, Llen, Lmax, D) == true) {
      S.split(L);

      uint32  val = Lnum++;
      uint32  cnt = 0;

      if      (S.numWords() == 1) {
        val = S.touint32(0);
        cnt = 1;
      }

      else if (S.numWords() == 2) {
        val = S.touint32(0);
        cnt = S.touint32(1);
      }

      else {
        fprintf(stderr, "too many words on line '%s'\n", L);
        exit(1);
      }

      while (_dataMax <= val)
        resizeArray(_data, _dataLen, _dataMax, 2 * _dataMax, resizeArray_copyData | resizeArray_clearNew);

      _data[val] += cnt;
      _dataSum   += cnt;

      _dataLen = max(_dataLen, val + 1);
    }

    AS_UTL_closeFile(D);
  };


  uint32   getValue(double d) {

    if (d < 0.0)  d = 0.0;
    if (d > 1.0)  d = 1.0;

    uint64  pos = 0;
    uint64  lim = (uint64)floor(_dataSum * d);

    if (pos >= _dataSum)   pos = _dataSum;

    while (_data[pos] < lim) {
      lim -= _data[pos];
      pos += 1;
    }

    assert(pos < _dataLen);

    return(pos);
  };


public:
  uint32   _dataLen;   //  Highest valid point
  uint32   _dataMax;   //  Number allocated
  uint32  *_data;

  uint64   _dataSum;
};

