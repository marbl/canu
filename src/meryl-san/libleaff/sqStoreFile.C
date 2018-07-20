
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStoreFile.H"
//#include "AS_UTL_fileIO.H"


sqStoreFile::sqStoreFile() {
  clear();
  seq = NULL;
}



sqStoreFile::sqStoreFile(const char *name) {

  clear();

  strncpy(_filename, name, FILENAME_MAX-1);

  seq = sqStore::sqStore_open(_filename);

  _numberOfSequences = seq->sqStore_getNumReads();
  //fprintf(stderr, "Opened '%s' with %u reads\n", _filename, _numberOfSequences);
}

sqStoreFile::~sqStoreFile() {
  seq->sqStore_close();
}



seqFile *
sqStoreFile::openFile(const char *name) {
  struct stat  st;

  //  Assume it's a sqStore if it is a directory, and the info / reads / blobs files exist.
  //
  //  Well, we used to check for the blobs file, but it doesn't exist if
  //  an object store is being used.

  char  infoName[FILENAME_MAX];
  char  readName[FILENAME_MAX];

  sprintf(infoName, "%s/info",  name);
  sprintf(readName, "%s/reads", name);

  if ((AS_UTL_fileExists(name, true) == false) ||
      (AS_UTL_fileExists(infoName) == false) ||
      (AS_UTL_fileExists(readName) == false))
    return(0L);

  //  Yup, probably a sqStore.  If it isn't, the sqStore() constructor blows up.

  return(new sqStoreFile(name));
}



bool
sqStoreFile::getSequence(uint32 iid,
                         char *&h, uint32 &hLen, uint32 &hMax,
                         char *&s, uint32 &sLen, uint32 &sMax) {

  if (iid > _numberOfSequences) {
    fprintf(stderr, "sqStoreFile::getSequence()-- iid %u exceeds number in store %u\n", iid, _numberOfSequences);
    return(false);
  }

  iid++;

  uint32  rLength = seq->sqStore_getRead(iid)->sqRead_sequenceLength();

  if (hMax < 32) {
    delete h;
    h    = new char [32];
    hMax = 32;
  }

  if (sMax < rLength) {
    delete s;
    s    = new char [rLength + 1];
    sMax = rLength;
  }

  hLen = sprintf(h, F_U32, iid);
  sLen = rLength;

  seq->sqStore_loadReadData(iid, &readData);
  memcpy(s, readData.sqReadData_getSequence(), sizeof(char) * rLength);

  s[sLen] = 0;

  return(true);
}



bool
sqStoreFile::getSequence(uint32 iid,
                         uint32 bgn, uint32 end, char *s) {

  if (iid > _numberOfSequences) {
    fprintf(stderr, "sqStoreFile::getSequence()-- iid %u exceeds number in store %u\n", iid, _numberOfSequences);
    return(false);
  }

  iid++;

  uint32  rLength = seq->sqStore_getRead(iid)->sqRead_sequenceLength();

  //fprintf(stderr, "return canu iid %u of length %u\n", iid, rLength);

  assert(bgn <  end);
  assert(bgn <= rLength);
  assert(end <= rLength);

  seq->sqStore_loadReadData(iid, &readData);
  memcpy(s, readData.sqReadData_getSequence() + bgn, sizeof(char) * (end - bgn));

  s[end-bgn] = 0;

  return(true);
}



void
sqStoreFile::clear(void) {

  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "SQSTORE");

  _numberOfSequences = 0;
}
