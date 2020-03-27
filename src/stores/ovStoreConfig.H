
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

#ifndef OVSTORECONFIG_H
#define OVSTORECONFIG_H


#include "runtime.H"

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

  ovStoreConfig(vector<char const *> &names, uint32 maxID) {
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

  ovStoreConfig(char const *configName) {
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

  void    loadConfig(char const *configName) {
    FILE *C = AS_UTL_openInputFile(configName);

    loadFromFile(_maxID,      "maxID",      C);
    loadFromFile(_numBuckets, "numBuckets", C);
    loadFromFile(_numSlices,  "numSlices",  C);
    loadFromFile(_sortMemory, "sortMemory", C);
    loadFromFile(_numInputs,  "numInputs",  C);

    _inputNames = new char * [_numInputs];

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32  nl = 0;

      loadFromFile(nl, "nameLen", C);

      _inputNames[ii] = new char [nl+1];

      loadFromFile(_inputNames[ii], "name", nl+1, C);
    }

    _inputToBucket = new uint32 [_numInputs];
    _readToSlice   = new uint16 [_maxID+1];

    loadFromFile(_inputToBucket, "inputToBucket", _numInputs, C);
    loadFromFile(_readToSlice,   "readToSlice",   _maxID+1,   C);

    AS_UTL_closeFile(C, configName);
  };

  void    writeConfig(char const *configName) {
    FILE *C = AS_UTL_openOutputFile(configName);

    writeToFile(_maxID,      "maxID",      C);
    writeToFile(_numBuckets, "numBuckets", C);
    writeToFile(_numSlices,  "numSlices",  C);
    writeToFile(_sortMemory, "sortMemory", C);
    writeToFile(_numInputs,  "numInputs",  C);

    for (uint32 ii=0; ii<_numInputs; ii++) {
      uint32 nl = strlen(_inputNames[ii]);

      writeToFile(nl,              "nameLen",       C);
      writeToFile(_inputNames[ii], "name",    nl+1, C);
    }

    writeToFile(_inputToBucket, "inputToBucket", _numInputs, C);
    writeToFile(_readToSlice,   "readToSlice",   _maxID + 1, C);

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

  char const *getInput(uint32 bucketNumber, uint32 fileNumber) {
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
