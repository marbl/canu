
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
 *    Brian P. Walenz beginning on 2019-APR-18
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

    delete [] L;
  };


  uint32   getValue(double d) {

    //  Imagine our input histogram (value occurences) is expanded into an
    //  array of (_dataSum) values, where each 'value' in the input histogram
    //  is listed 'occurences' times (and that it's sorted).
    //
    //  For input 0 <= 'd' < 1, we want to return the 'value' that is at that
    //  position in the array.

    if (d < 0.0)  d = 0.0;   //  Limit our input (random) number
    if (d > 1.0)  d = 1.0;   //  to valid scaling range.

    uint64  lim = (uint64)floor(_dataSum * d);

    //  Scan the _data, looking for the _data element that contains the 'lim'th
    //  entry in the (imagined) array.

    uint64  value = 0;              //  value we're testing.

    while (_data[value] <= lim) {   //  If _data[value] is more than the current
      lim -= _data[value];          //  limit, we're found the value, otherwise,
      value += 1;                   //  decrement the limit by the occurrences
    }                               //  for this value and move to the next.

    assert(value < _dataLen);

    return(value);
  };


public:
  uint32   _dataLen;   //  Length of the valid data in _data: _data[_dataLen-1] is the last valid data point.
  uint32   _dataMax;   //  Allocated size of _data.
  uint32  *_data;      //  Number of data points of value x:  _data[x] == v <-> 'x' was present 'v' times in the input.

  uint64   _dataSum;   //  Number of points in the input.  It's the sum of all _data[x].
};

