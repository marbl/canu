
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
 *    Brian P. Walenz beginning on 2018-MAR-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef OVSTORECONFIG_H
#define OVSTORECONFIG_H


#include "AS_global.H"

#include <vector>
using namespace std;



class ovStoreConfig {
public:
  ovStoreConfig() {
    _maxID         = 0;

    _numBuckets    = 0;
    _numSlices     = 0;
    _sortMemory    = 0;

    _numInputs     = 0;
    _inputNames    = NULL;

    _inputToBucket = NULL;
    _readToSlice   = NULL;
  };

  ovStoreConfig(vector<char *> &names, uint32 maxID) {
    _maxID         = maxID;

    _numBuckets    = 0;
    _numSlices     = 0;
    _sortMemory    = 0;

    _numInputs     = names.size();
    _inputNames    = new char * [_numInputs];

    for (uint32 ii=0; ii<names.size(); ii++)
      _inputNames[ii] = duplicateString(names[ii]);

    _inputToBucket = new uint32 [_numInputs];
    _readToSlice   = new uint16 [_maxID+1];
  };

  ovStoreConfig(const char *configName) {
    _maxID         = 0;

    _numBuckets    = 0;
    _numSlices     = 0;
    _sortMemory    = 0;

    _numInputs     = 0;
    _inputNames    = NULL;

    _inputToBucket = NULL;
    _readToSlice   = NULL;

    loadConfig(configName);
  };

  ~ovStoreConfig() {
    for (uint32 ii=0; ii<_numInputs; ii++)
      delete [] _inputNames[ii];

    delete [] _inputNames;

    delete [] _inputToBucket;
    delete [] _readToSlice;
  };

  void    loadConfig(const char *configName) {
    FILE *C = AS_UTL_openInputFile(configName);

    AS_UTL_safeRead(C, &_maxID,      "maxID",      sizeof(uint32), 1);
    AS_UTL_safeRead(C, &_numBuckets, "numBuckets", sizeof(uint32), 1);
    AS_UTL_safeRead(C, &_numSlices,  "numSlices",  sizeof(uint32), 1);
    AS_UTL_safeRead(C, &_sortMemory, "sortMemory", sizeof(double), 1);
    AS_UTL_safeRead(C, &_numInputs,  "numInputs",  sizeof(uint32), 1);

    _inputNames = new char * [_numInputs];

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32  nl = 0;

      AS_UTL_safeRead(C, &nl, "nameLen", sizeof(uint32), 1);

      _inputNames[ii] = new char [nl+1];

      AS_UTL_safeRead(C, _inputNames[ii], "name", sizeof(char), nl+1);
    }

    _inputToBucket = new uint32 [_numInputs];
    _readToSlice   = new uint16 [_maxID+1];

    AS_UTL_safeRead(C, _inputToBucket, "inputToBucket", sizeof(uint32), _numInputs);
    AS_UTL_safeRead(C, _readToSlice,   "readToSlice",   sizeof(uint16), _maxID+1);

    AS_UTL_closeFile(C, configName);
  };

  void    writeConfig(const char *configName) {
    FILE *C = AS_UTL_openOutputFile(configName);

    AS_UTL_safeWrite(C, &_maxID,      "maxID",      sizeof(uint32), 1);
    AS_UTL_safeWrite(C, &_numBuckets, "numBuckets", sizeof(uint32), 1);
    AS_UTL_safeWrite(C, &_numSlices,  "numSlices",  sizeof(uint32), 1);
    AS_UTL_safeWrite(C, &_sortMemory, "sortMemory", sizeof(double), 1);
    AS_UTL_safeWrite(C, &_numInputs,  "numInputs",  sizeof(uint32), 1);

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32 nl = strlen(_inputNames[ii]);

      AS_UTL_safeWrite(C, &nl,             "nameLen", sizeof(uint32), 1);
      AS_UTL_safeWrite(C, _inputNames[ii], "name",    sizeof(char),   nl+1);
    }

    AS_UTL_safeWrite(C, _inputToBucket, "inputToBucket", sizeof(uint32), _numInputs);
    AS_UTL_safeWrite(C, _readToSlice,   "readToSlice",   sizeof(uint16), _maxID+1);

    AS_UTL_closeFile(C, configName);

    fprintf(stderr, "\n");
    fprintf(stderr, "Saved configuration to '%s'.\n", configName);
  };

  uint32  numBuckets(void) { return(_numBuckets); };
  uint32  numSlices(void)  { return(_numSlices);  };
  double  sortMemory(void) { return(_sortMemory); };


  uint32  numInputs(uint32 bucketNumber) {
    uint32 ni = 0;

    bucketNumber--;  //  Internally starting at 0, externally at 1.

    for (uint32 ii=0; ii<_numInputs; ii++)
      if (_inputToBucket[ii] == bucketNumber)
        ni++;

    return(ni);
  };

  char   *getInput(uint32 bucketNumber, uint32 fileNumber) {
    uint32 ni = 0;

    bucketNumber--;  //  Internally starting at 0, externally at 1.

    for (uint32 ii=0; ii<_numInputs; ii++)
      if (_inputToBucket[ii] == bucketNumber)
        if (ni++ == fileNumber)
          return(_inputNames[ii]);

    return(NULL);
  }


  uint32  getAssignedSlice(uint32 id) {
    return(_readToSlice[id] + 1);
  };


  void    assignReadsToSlices(sqStore *seq,
                              uint64   minMemory,
                              uint64   maxMemory);

private:
  uint32     _maxID;

  uint32     _numBuckets;
  uint32     _numSlices;
  double     _sortMemory;      //  Expected maximum memory usage in GB (for sorting).

  uint32     _numInputs;       //  Number of input ovb files.
  char     **_inputNames;      //  Input ovb files.

  uint32    *_inputToBucket;   //  Maps an input name to a bucket.
  uint16    *_readToSlice;      //  Map each read ID to a slice.
};



#endif  //  OVSTORECONFIG_H