#include "gkStoreFile.H"
//#include "AS_UTL_fileIO.H"


gkStoreFile::gkStoreFile() {
  clear();
}



gkStoreFile::gkStoreFile(const char *name) {

  clear();

  strcpy(_filename, name);

  gkp = new gkStore(_filename);

  _numberOfSequences = gkp->gkStore_getNumReads();
  fprintf(stderr, "Opened '%s' with %u reads\n", _filename, _numberOfSequences);
}

gkStoreFile::~gkStoreFile() {
  delete gkp;
}



seqFile *
gkStoreFile::openFile(const char *name) {
  struct stat  st;

  //  Assume it's a gkStore if it is a directory, and the files 'name/info' and 'name/reads' exist.

  char  infoName[FILENAME_MAX];
  char  readName[FILENAME_MAX];

  sprintf(infoName, "%s/info",  name);
  sprintf(readName, "%s/reads", name);
  
  if ((AS_UTL_fileExists(name, true) == false) &&
      (AS_UTL_fileExists(infoName) == false) &&
      (AS_UTL_fileExists(readName) == false))
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

  //fprintf(stderr, "return ca3g iid %u of length %u\n", iid, rLength);

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
