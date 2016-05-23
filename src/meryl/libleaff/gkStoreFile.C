
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
 *    Brian P. Walenz from 2015-FEB-04 to 2015-AUG-14
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "gkStoreFile.H"
//#include "AS_UTL_fileIO.H"


gkStoreFile::gkStoreFile() {
  clear();
}



gkStoreFile::gkStoreFile(const char *name) {

  clear();

  strcpy(_filename, name);

  gkp = gkStore::gkStore_open(_filename);

  _numberOfSequences = gkp->gkStore_getNumReads();
  //fprintf(stderr, "Opened '%s' with %u reads\n", _filename, _numberOfSequences);
}

gkStoreFile::~gkStoreFile() {
  gkp->gkStore_close();
}



seqFile *
gkStoreFile::openFile(const char *name) {
  struct stat  st;

  //  Assume it's a gkStore if it is a directory, and the info / reads / blobs files exist.

  char  infoName[FILENAME_MAX];
  char  readName[FILENAME_MAX];
  char  blobName[FILENAME_MAX];

  sprintf(infoName, "%s/info",  name);
  sprintf(readName, "%s/reads", name);
  sprintf(blobName, "%s/blobs", name);

  if ((AS_UTL_fileExists(name, true) == false) ||
      (AS_UTL_fileExists(infoName) == false) ||
      (AS_UTL_fileExists(readName) == false) ||
      (AS_UTL_fileExists(blobName) == false))
    return(0L);

  //  Yup, probably a gkStore.  If it isn't, the gkStore() constructor blows up.

  return(new gkStoreFile(name));
}



bool
gkStoreFile::getSequence(uint32 iid,
                         char *&h, uint32 &hLen, uint32 &hMax,
                         char *&s, uint32 &sLen, uint32 &sMax) {
  if (iid > _numberOfSequences) {
    fprintf(stderr, "gkStoreFile::getSequence()-- iid %u exceeds number in store %u\n", iid, _numberOfSequences);
    return(false);
  }

  iid++;

  uint32  rLength = gkp->gkStore_getRead(iid)->gkRead_sequenceLength();

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

  gkp->gkStore_loadReadData(iid, &readData);
  memcpy(s, readData.gkReadData_getSequence(), sizeof(char) * rLength);

  s[sLen] = 0;

  return(true);
}



bool
gkStoreFile::getSequence(uint32 iid,
                         uint32 bgn, uint32 end, char *s) {

  if (iid > _numberOfSequences) {
    fprintf(stderr, "gkStoreFile::getSequence()-- iid %u exceeds number in store %u\n", iid, _numberOfSequences);
    return(false);
  }

  iid++;

  uint32  rLength = gkp->gkStore_getRead(iid)->gkRead_sequenceLength();

  //fprintf(stderr, "return canu iid %u of length %u\n", iid, rLength);

  assert(bgn <  end);
  assert(bgn <= rLength);
  assert(end <= rLength);

  gkp->gkStore_loadReadData(iid, &readData);
  memcpy(s, readData.gkReadData_getSequence() + bgn, sizeof(char) * (end - bgn));

  s[end-bgn] = 0;

  return(true);
}



void
gkStoreFile::clear(void) {

  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "GKSTORE");

  _numberOfSequences = 0;
}
