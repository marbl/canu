#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "libbri.H"
#include "fasta.H"

//  Fri Jul 12 11:01:13 EDT 2002: Improved error and status messages.


//  XXX:  error handling of fread(), fwrite() and rewind() (probably all done, but need to check)
//  XXX:  This would probably benefit from being threaded.


//  Magic number; the last number is the version
//
static char const magicNumber[8] = { 0x45, 0x52, 0x50, 0x6c, 0x61, 0x6e, 0x64, 0x01 };


u32bit
readBuffer(int file, char *fastaname, char *buffer, u32bit bufferMax, u32bit &bufferPos) {
  bufferPos = 0;
  errno = 0;
  int len = (u32bit)read(file, buffer, bufferMax * sizeof(char));
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() couldn't read %d bytes from the FastA file.\n%s\n",
            fastaname, bufferMax * sizeof(char), strerror(errno));
    exit(1);
  }
  return(len);
}


void
buildIndex(char *fastaname, char *indexname) {
  u32bit                     numSeqs = 0;
  u32bit                     maxSeqs = 16;
  _idxfa_desc               *theSeqs = new _idxfa_desc [maxSeqs];

  u32bit                     bufferPos = 0;
  u32bit                     bufferLen = 0;
  u32bit                     bufferMax = 16 * 1024 * 1024;
  char                      *buffer  = new char [bufferMax];

  _idxfa_pos                 filePos   = 0;
  _idxfa_len                 defLen    = 0;
  _idxfa_pos                 defStart  = 0;
  _idxfa_len                 seqLen    = 0;
  _idxfa_pos                 seqStart  = 0;

  errno = 0;
  int file = open(fastaname, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() couldn't open the FastA file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  bufferLen = readBuffer(file, fastaname, buffer, bufferMax, bufferPos);

  //  Generate data for the current sequence.
  //
  while (bufferLen > 0) {

    //  Skip any whitespace
    //
    while ((bufferLen > 0) &&
           isspace(buffer[bufferPos])) {
      filePos++;
      bufferPos++;

      if (bufferPos >= bufferLen)
        bufferLen = readBuffer(file, fastaname, buffer, bufferMax, bufferPos);
    }

    //  Save the start of the defline
    //
    defStart = filePos;
    defLen   = 0;

    //  Count the length of the defline
    //
    while ((bufferLen > 0) &&
           (buffer[bufferPos] != '\n')) {
      defLen++;
      filePos++;
      bufferPos++;

      if (bufferPos >= bufferLen)
        bufferLen = readBuffer(file, fastaname, buffer, bufferMax, bufferPos);
    }

    //  Skip any whitespace
    //
    while ((bufferLen > 0) &&
           isspace(buffer[bufferPos])) {
      filePos++;
      bufferPos++;

      if (bufferPos >= bufferLen)
        bufferLen = readBuffer(file, fastaname, buffer, bufferMax, bufferPos);
    }

    //  Save the start of the sequence
    //
    seqStart = filePos;
    seqLen   = 0;

    //  Count the length of the sequence
    //
    while ((bufferLen > 0) &&
           (buffer[bufferPos] != '>')) {
      if (!isspace(buffer[bufferPos]))
        seqLen++;

      filePos++;
      bufferPos++;

      if (bufferPos >= bufferLen)
        bufferLen = readBuffer(file, fastaname, buffer, bufferMax, bufferPos);
    }

    //  Make the index bigger
    //
    if (numSeqs >= maxSeqs) {
      if (maxSeqs == 0)
        maxSeqs = 16;
      maxSeqs *= 2;
      
      _idxfa_desc *sd = new _idxfa_desc [ maxSeqs ];

      for (u32bit i=0; i<numSeqs; i++) {
        sd[i]._headerStart = theSeqs[i]._headerStart;
        sd[i]._headerLen   = theSeqs[i]._headerLen;
        sd[i]._seqStart    = theSeqs[i]._seqStart;
        sd[i]._seqLen      = theSeqs[i]._seqLen;
      }
      
      delete [] theSeqs;
      theSeqs = sd;
    }


    //  Add the new sequence description to the list.
    //
    theSeqs[numSeqs]._headerStart = defStart;
    theSeqs[numSeqs]._headerLen   = defLen;
    theSeqs[numSeqs]._seqStart    = seqStart;
    theSeqs[numSeqs]._seqLen      = seqLen;

    numSeqs++;
  }

  errno = 0;
  close(file);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() couldn't close the FastA file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }



  ///////////////////////////////////////
  //
  //  Write the data file; version, number of sequences and sequence
  //  descriptions.
  //
  errno = 0;
  int info = open(indexname,
                  O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                  S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() can't open FastA index file for writing.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  errno = 0;
  write(info, magicNumber, sizeof(char) * 8);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() can't write FastA index file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  errno = 0;
  write(info, &numSeqs,  sizeof(u32bit));
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() can't write FastA index file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  errno = 0;
  write(info, theSeqs,   sizeof(_idxfa_desc) * numSeqs);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() can't write FastA index file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  errno = 0;
  close(info);
  if (errno) {
    fprintf(stderr, "%s: FastA::buildIndex() can't close FastA index file '%s'.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  delete [] theSeqs;
}




//  Open the index, building/rebuilding when necessary
//
int
openIndex(char *fastaname, bool beVerbose) {
  struct stat  fStat;
  struct stat  iStat;


  //  Build the index filename.  It's just the prefix passed in with
  //  '.fastaidx' tacked on the end.  If the prefix has '.fasta' at
  //  the end, strip it off so that we don't get '*.fasta.fastaidx'.
  //
  int   fastanamelen = (int)strlen(fastaname);
  char *indexname    = new char [fastanamelen + 16];

  strcpy(indexname, fastaname);
  if (strcmp(fastaname+fastanamelen-6, ".fasta") == 0)
    indexname[fastanamelen-6] = 0;
  strcat(indexname, ".fastaidx");


  //  Make sure that the fasta file actually exists.
  //
  errno = 0;
  if (stat(fastaname, &fStat) != 0) {
    fprintf(stderr, "%s: FastA::openIndex() couldn't stat the FastA file.\n%s\n", fastaname, strerror(errno));
    exit(1);
  }

  int idx = 0;

  //  Determine if the index exists, if it is up to date, building/rebuilding
  //  as necessary.
  //
  if (stat(indexname, &iStat) == 0) {
    if (fStat.st_mtime <= iStat.st_mtime) {
      errno = 0;
      idx = open(indexname, O_RDONLY | O_LARGEFILE);
      if (errno) {
        fprintf(stderr, "%s: FastA::openIndex() couldn't open the FastA index file.\n%s\n", fastaname, strerror(errno));
        exit(1);
      }
    } else {
      if (beVerbose)
        fprintf(stderr, "%s: Index is out of date.  Rebuilding.\n", fastaname);
    }
  } else {
    if (beVerbose)
      fprintf(stderr, "%s: Index not found.  Building.\n", fastaname);
  }

 again:

  //  Build/rebuild the index
  //
  if (idx == 0) {
    buildIndex(fastaname, indexname);
    errno = 0;
    idx = open(indexname, O_RDONLY | O_LARGEFILE);
    if (errno) {
      fprintf(stderr, "%s: FastA::openIndex() couldn't open the FastA index file.\n%s\n", fastaname, strerror(errno));
      exit(1);
    }
  }

  //  Verify that the index is the correct version.  If it is not,
  //  rebuild it.
  //
  char   magicNumberInFile[8];
  char   version = 0;
  read(idx, magicNumberInFile, sizeof(char) * 8);

  if ((magicNumberInFile[0] == magicNumber[0]) &&
      (magicNumberInFile[1] == magicNumber[1]) &&
      (magicNumberInFile[2] == magicNumber[2]) &&
      (magicNumberInFile[3] == magicNumber[3]) &&
      (magicNumberInFile[4] == magicNumber[4]) &&
      (magicNumberInFile[5] == magicNumber[5]) &&
      (magicNumberInFile[6] == magicNumber[6])) {
    version = magicNumberInFile[7];
  }

  if (version != 1) {
    if (beVerbose)
      fprintf(stderr, "%s: Index is version %d; Building a version %d index.\n", fastaname, version, magicNumber[7]);
    errno = 0;
    close(idx);
    if (errno) {
      fprintf(stderr, "%s: FastA::openIndex() couldn't close the FastA index file.\n%s\n", fastaname, strerror(errno));
      exit(1);
    }
    idx = 0L;
    goto again;
  }

  delete [] indexname;

  return(idx);
}




FastA::FastA(char *fastaname, bool randomAccess, bool beVerbose) {

  //  Initialize things to their default values.
  //
  _file        = 0;
  _eof         = false;

  _numSeqs     = 0;
  _maxSeqs     = 0;
  _descr       = 0L;

  //
  //  This buffers the input file.
  //
  _bufferLen   = 0;
  _bufferPos   = 0;
  _bufferMax   = 512 * 1024;
  _buffer      = new unsigned char [_bufferMax];

  //  Try to open the fasta file 'fastaname'.  This is done for both
  //  sequential and random access.
  //
  //  If the fastaname is 0L, use stdin, and force sequential access.
  //
  if (fastaname == 0L) {
    _file        = fileno(stdin);
    randomAccess = false;
    _fastaname   = new char [6];
    strcpy(_fastaname, "stdin");
  } else {
    errno = 0;
    _file = open(fastaname, O_RDONLY | O_LARGEFILE);
    if (errno) {
      fprintf(stderr, "%s: FastA::FastA() couldn't open FastA file.\n%s\n", fastaname, strerror(errno));
      exit(1);
    }
    _fastaname   = new char [strlen(fastaname)+1];
    strcpy(_fastaname, fastaname);
  }
 
  //  Random Access??  Read the index, building one if necessary.
  //
  if (randomAccess) {
    int idx = openIndex(_fastaname, beVerbose);

    if (idx) {
      errno = 0;
      read(idx, &_numSeqs, sizeof(u32bit));
      if (errno) {
        fprintf(stderr, "%s: FastA::FastA() couldn't read the FastA index file.\n%s\n", _fastaname, strerror(errno));
        exit(1);
      }

      if (beVerbose)
#ifdef TRUE64BIT
        fprintf(stderr, "%s: Index found, reading %u sequence descriptions.\n", _fastaname, _numSeqs);
#else
        fprintf(stderr, "%s: Index found, reading %lu sequence descriptions.\n", _fastaname, _numSeqs);
#endif
      _descr = new _idxfa_desc [_numSeqs];


      errno = 0;
      read(idx, _descr, sizeof(_idxfa_desc) * _numSeqs);
      if (errno) {
        fprintf(stderr, "%s: FastA::FastA() couldn't read the FastA index file.\n%s\n", _fastaname, strerror(errno));
        exit(1);
      }

      errno = 0;
      close(idx);
      if (errno) {
        fprintf(stderr, "%s: FastA::FastA() couldn't close the FastA index file.\n%s\n", _fastaname, strerror(errno));
        exit(1);
      }
    }
  }
}


FastA::~FastA() {
  delete [] _buffer;
  delete [] _descr;
  delete [] _fastaname;
  close(_file);
}


void
FastA::first(FastABuffer &b) {
  _eof = false;
  errno = 0;
  lseek(_file, 0, SEEK_SET);
  if (errno) {
    fprintf(stderr, "%s: FastA::first() failed.\n%s\n", _fastaname, strerror(errno));
    exit(1);
  }

  _bufferPos = 0;
  _bufferLen = 0;

  next(b);

  b.setIndex(0);
}



bool
FastA::next(FastABuffer &b) {

  u32bit  idx = b.index() + 1;
  b.clear();
  b.setIndex(idx);

  //  Fill the buffer
  //
  fillBuffer();

  //  Move to the next fasta file.  If we hit eof here, we return
  //  false (no more sequences).
  //
  while ((eof() == false) && (theCharacter() != '>'))
    nextCharacter();

  //  If we hit the end of file, there are no more sequences,
  //  return false.
  //
  if (eof())
    return(false);

  //  Copy the description header
  //
  while ((eof() == false) && (theCharacter() != '\n')) {
    b.pushHeader(theCharacter());
    nextCharacter();
  }
  b.pushHeader(0);

  //  Bump over the '\n'
  //
  nextCharacter();

  //  Copy the sequence into _seq, up until we hit another '>'
  //
  while ((eof() == false) && (theCharacter() != '>')) {
    if (!isspace(theCharacter()))
      b.pushSequence(theCharacter());
    nextCharacter();
  }
  b.terminateSequence();

  //  If we hit the end of file, we have just read the last sequence.
  //  Reset _eof to be false.  Otherwise, theCharacter() is a
  //  '>', and we don't need to do anything.
  //
  _eof = false;

  return(true);
}


//  Fail if the id is invalid, otherwise load the specified sequence.
//
bool
FastA::seek(FastABuffer &b, u32bit id) {

  if (id >= _numSeqs) {
    fprintf(stderr, "%s: FastA::seek() ID out of range: %u >= %u\n", _fastaname, id, _numSeqs);
    return(false);
  }

  if (id == b.index()) {
    //  Not a very useful message.
    //fprintf(stderr, "%s: FastA::seek() ID %u already loaded.\n", fastaname, id);
    return(true);
  }

#ifdef DEBUG_SEEK
  fprintf(stderr, "%s: FastA::seek() to pos %llu for seq %lu\n",
          _fastaname, _descr[id]._headerStart, id);
#endif

 again:
  errno = 0;
  off_t lsr = lseek(_file, _descr[id]._headerStart, SEEK_SET);
  if (errno) {
    if (errno == EAGAIN) {
      fprintf(stderr, "%s: FastA::seek() got EAGAIN.  Attempting again.\n", _fastaname);
      goto again;
    }
    fprintf(stderr, "%s: FastA::seek() failed!\n%s\n", _fastaname, strerror(errno));
    exit(1);
  }
  if (lsr != _descr[id]._headerStart) {
    fprintf(stderr, "%s: FastA::seek() wanted pos %llu for seq %lu but got %llu??\n",
            _fastaname, _descr[id]._headerStart, id, lsr);
  }

  _bufferPos = 0;
  _bufferLen = 0;

  next(b);

  b.setIndex(id);

  return(true);
}
