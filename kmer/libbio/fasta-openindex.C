#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "bio++.H"

//  If the index exists, check the global information in the index
//  against the version of this code, and the fasta file it is
//  indexing.  If everything checks out OK, just return.  Otherwise,
//  continue and build a new index.
//
//  Returns false if
//    The index doesn't exist
//    The index exists, but is empty
//    The index doesn't match the fasta file and would need to be rebuilt
//
bool
FastAFile::isIndexValid(u32bit indextype, bool beVerbose) {
  errno = 0;

  //  No indexname?  Must be stdin.
  if (_indexname == 0L)
    return(false);

  //  Index exists?
  struct stat    st;
  stat(_indexname, &st);
  if (errno)
    return(false);

  //  Index is empty?
  if (st.st_size == 0)
    return(false);

  //  Stat the -> fasta <- file to get the modification time
  // 
  stat(_filename, &st);
  if (errno) {
    fprintf(stderr, "ERROR: Can't stat '%s': %s\n", _filename, strerror(errno));
    exit(1);
  }

  if (st.st_size == 0) {
    //  file is empty
    fprintf(stderr, "ERROR: '%s' is empty?\n", _filename);
    exit(1);
  }

  int indexfile = open(_indexname, O_RDONLY | O_LARGEFILE);
  if (errno) {
    //  Index doesn't exist??
    fprintf(stderr, "ERROR: Index file exists, but can't be opened for reading: %s\n", strerror(errno));
    exit(1);
  }

  //  Opened an index, so read the description.
  //
  _idxfa_global  theGlobalDesc;
  read(indexfile, &theGlobalDesc, sizeof(_idxfa_global));
  if (errno) {
    //  can't read
    fprintf(stderr, "ERROR: Index file exists, but can't read description: %s\n", strerror(errno));
    exit(1);
  }

  //  Close the index
  //
  errno = 0;
  close(indexfile);
  if (errno == 0) {
    //  can't close
  }

  //  Fix endianess issues in the index.
  //
#if FASTA_VERSIONNUMBER > 4
  if (theGlobalDesc._magic == u64bitSwap(FASTA_MAGICNUMBER)) {
    //fprintf(stderr, "WARNING:  Fixing index for endianess swap!\n");

    theGlobalDesc._magic                 = u64bitSwap(theGlobalDesc._magic);
    theGlobalDesc._version               = u32bitSwap(theGlobalDesc._version);
    theGlobalDesc._indexType             = u32bitSwap(theGlobalDesc._indexType);
    theGlobalDesc._fastaType             = u32bitSwap(theGlobalDesc._fastaType);
    theGlobalDesc._numberOfSequences     = u32bitSwap(theGlobalDesc._numberOfSequences);
    theGlobalDesc._fastaFileSize         = idxfaPosSwap(theGlobalDesc._fastaFileSize);
    theGlobalDesc._fastaModificationTime = idxfaTimeSwap(theGlobalDesc._fastaModificationTime);
    theGlobalDesc._fastaCreationTime     = idxfaTimeSwap(theGlobalDesc._fastaCreationTime);
    theGlobalDesc._seqlineLength         = idxfaLenSwap(theGlobalDesc._seqlineLength);
    theGlobalDesc._seqendlLength         = idxfaLenSwap(theGlobalDesc._seqendlLength);
    theGlobalDesc._fixedWidth            = u32bitSwap(theGlobalDesc._fixedWidth);
    theGlobalDesc._squeezedSequences     = u32bitSwap(theGlobalDesc._squeezedSequences);
  }
#endif


  //  Check that the magic number and version are correct and that
  //  the modification time of the fasta file is the same as stored
  //  in the index.
  //
  if ((theGlobalDesc._magic                          == FASTA_MAGICNUMBER) &&
      (theGlobalDesc._version                        == FASTA_VERSIONNUMBER) &&
      (theGlobalDesc._fastaModificationTime          == st.st_mtime) &&
      ((theGlobalDesc._indexType & FASTA_INDEX_MASK) >= (indextype & FASTA_INDEX_MASK))) {
    return(true);
  } else {
    if (beVerbose) {
      fprintf(stderr, "WARNING: Index found, but stale or wrong type!\n");

      if (theGlobalDesc._magic != FASTA_MAGICNUMBER)
        fprintf(stderr, "           Magic number incorrect; perhaps not an index file?  got="u64bitHEX" expected="u64bitHEX"\n",
                theGlobalDesc._magic, (u64bit)FASTA_MAGICNUMBER);

      if (theGlobalDesc._version != FASTA_VERSIONNUMBER)
        fprintf(stderr, "           Version number incorrect; got "u32bitFMT", expected "u32bitFMT".\n",
                theGlobalDesc._version, (u32bit)FASTA_VERSIONNUMBER);

      if (theGlobalDesc._fastaModificationTime != st.st_mtime)
        fprintf(stderr, "           File age not correct; got "s64bitFMT", expected %ld.\n",
                theGlobalDesc._fastaModificationTime, st.st_mtime);

      if ((theGlobalDesc._indexType & FASTA_INDEX_MASK) < (indextype & FASTA_INDEX_MASK))
        fprintf(stderr, "           Type of index insufficient; got %s, need %s.\n",
                indexTypeNames(theGlobalDesc._indexType),
                indexTypeNames(indextype));
    }
  }

  return(false);
}



void
FastAFile::openIndex(u32bit indextypetoload) {


  //  If we've been told to open an index, but we're not random
  //  access, complain.
  //
  if (_isStreamInput) {
    fprintf(stderr, "FastAFile()--  '%s' is a stream and not valid for indexing!\n", _filename);
    return;
  }

  //  If the index is already open, and it's the same type as desired, we
  //  can just return.  Downgrading an index is not allowed -- if you open
  //  a FASTA_INDEX_PLUS_DEFLINES, you're stuck with it.
  //
  if ((_isRandomAccess) &&
      (indextypetoload & FASTA_INDEX_MASK <= _theGlobalDesc._indexType & FASTA_INDEX_MASK))
    return;


  //  If no index, or it's not current, build a new one.
  //
  if (isIndexValid(indextypetoload, true) == false)
    createIndex(indextypetoload);

  //  Index exists, or we wouldn't have made it here.
  //
  _isRandomAccess = true;

  //  Open the index, read it in.
  //
  errno = 0;
  int indexfile = open(_indexname, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "FastAFile()-- couldn't open the index '%s': %s\n", _indexname, strerror(errno));
    exit(1);
  }

  //
  //  Read the index
  //
  read(indexfile, &_theGlobalDesc, sizeof(_idxfa_global));
  if (errno) {
    fprintf(stderr, "FastAFile()-- couldn't read description from the index '%s': %s\n", _indexname, strerror(errno));
    exit(1);
  }


  //  Fix endianess issues in the index.
  //
#if FASTA_VERSIONNUMBER > 4
  bool swapIndex = _theGlobalDesc._magic == u64bitSwap(FASTA_MAGICNUMBER);

  if (swapIndex) {
    //fprintf(stderr, "WARNING:  Fixing index for endianess swap!\n");

    _theGlobalDesc._magic                 = u64bitSwap(_theGlobalDesc._magic);
    _theGlobalDesc._version               = u32bitSwap(_theGlobalDesc._version);
    _theGlobalDesc._indexType             = u32bitSwap(_theGlobalDesc._indexType);
    _theGlobalDesc._fastaType             = u32bitSwap(_theGlobalDesc._fastaType);
    _theGlobalDesc._numberOfSequences     = u32bitSwap(_theGlobalDesc._numberOfSequences);
    _theGlobalDesc._fastaFileSize         = idxfaPosSwap(_theGlobalDesc._fastaFileSize);
    _theGlobalDesc._fastaModificationTime = idxfaTimeSwap(_theGlobalDesc._fastaModificationTime);
    _theGlobalDesc._fastaCreationTime     = idxfaTimeSwap(_theGlobalDesc._fastaCreationTime);
    _theGlobalDesc._seqlineLength         = idxfaLenSwap(_theGlobalDesc._seqlineLength);
    _theGlobalDesc._seqendlLength         = idxfaLenSwap(_theGlobalDesc._seqendlLength);
    _theGlobalDesc._fixedWidth            = u32bitSwap(_theGlobalDesc._fixedWidth);
    _theGlobalDesc._squeezedSequences     = u32bitSwap(_theGlobalDesc._squeezedSequences);
  }
#endif


  _theSeqs = new _idxfa_desc [ _theGlobalDesc._numberOfSequences ];

  read(indexfile, _theSeqs, sizeof(_idxfa_desc) * _theGlobalDesc._numberOfSequences);
  if (errno) {
    fprintf(stderr, "FastAFile()-- couldn't read sequence descriptions from the index '%s': %s\n", _indexname, strerror(errno));
    exit(1);
  }

#if FASTA_VERSIONNUMBER > 4
  if (swapIndex) {
    for (u32bit i=0; i<_theGlobalDesc._numberOfSequences; i++) {
      _theSeqs[i]._headerStart = idxfaPosSwap(_theSeqs[i]._headerStart);
      _theSeqs[i]._seqStart    = idxfaPosSwap(_theSeqs[i]._seqStart);
      _theSeqs[i]._headerLen   = idxfaLenSwap(_theSeqs[i]._headerLen);
      _theSeqs[i]._seqLen      = idxfaLenSwap(_theSeqs[i]._seqLen);
    }
  }
#endif


  //  If indextype is FASTA_INDEX_ANY, open the index as reported by
  //  the file.  Otherwise, open the index as specified.
  //
  if ((indextypetoload & FASTA_INDEX_MASK) == FASTA_INDEX_ANY)
    indextypetoload = _theGlobalDesc._indexType;
  else
    _theGlobalDesc._indexType = indextypetoload;

  if (((_theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_IDS) ||
      ((_theGlobalDesc._indexType & FASTA_INDEX_MASK) == FASTA_INDEX_PLUS_DEFLINES)) {
    read(indexfile, &_theNamesLen, sizeof(u32bit));
    if (errno) {
      fprintf(stderr, "FastAFile()-- couldn't read lengths of names from the index '%s': %s\n", _indexname, strerror(errno));
      exit(1);
    }

#if FASTA_VERSIONNUMBER > 4
    if (swapIndex)
      _theNamesLen = u32bitSwap(_theNamesLen);
#endif

    _theNames = new char [sizeof(char) * _theNamesLen];

    read(indexfile, _theNames, sizeof(char) * _theNamesLen);
    if (errno) {
      fprintf(stderr, "FastAFile()-- couldn't read names from the index '%s': %s\n", _indexname, strerror(errno));
      exit(1);
    }
  }

  errno = 0;
  close(indexfile);
  if (errno)
    fprintf(stderr, "FastAFile()-- couldn't close the index '%s': %s\n", _indexname, strerror(errno));
}


//  Analyze the index, reset the bufferSize of the readBuffer if the
//  sequences are short (e.g., ESTs, SNPs).  This will make a huge
//  difference in the completely random access case -- instead of
//  always reading 32k (default buffer size) for each sequence, we
//  read a little more than average.  But, it sucks for sequential
//  -- we'll be doing a read on almost every sequence.
//
void
FastAFile::optimizeRandomAccess(void) {

  if (_isRandomAccessOpt)
    return;

  if (_isRandomAccess == false)
    openIndex();

  if (_theGlobalDesc._numberOfSequences < 10)
    return;

  u64bit  aveLen = 0;
  for (u32bit i=0; i<_theGlobalDesc._numberOfSequences; i++)
    aveLen += _theSeqs[i]._seqLen + _theSeqs[i]._headerLen;
  aveLen /= _theGlobalDesc._numberOfSequences;

  u64bit stdDev = 0;
  for (u32bit i=0; i<_theGlobalDesc._numberOfSequences; i++)
    stdDev += ((_theSeqs[i]._seqLen + _theSeqs[i]._headerLen - aveLen) *
               (_theSeqs[i]._seqLen + _theSeqs[i]._headerLen - aveLen));
  stdDev /= _theGlobalDesc._numberOfSequences - 1;

  stdDev = (u64bit)ceil(sqrt((double)stdDev));

  fprintf(stderr, "FastAFile::optimizeRandomAccess()-- For "u32bitFMT" seqs, ave="u64bitFMT" stddev="u64bitFMT", reset buffer to "u64bitFMT"\n",
          _theGlobalDesc._numberOfSequences, aveLen, stdDev, aveLen + stdDev);

  if (aveLen + stdDev < 32768) {

    //  Make it a nice power of two
    aveLen += stdDev;
    aveLen |= aveLen >> 1;
    aveLen |= aveLen >> 2;
    aveLen |= aveLen >> 4;
    aveLen |= aveLen >> 8;
    aveLen |= aveLen >> 16;
    aveLen |= aveLen >> 32;
    aveLen++;

    fprintf(stderr, "FastAFile::optimizeRandomAccess()-- Make new filebuffer of size "u64bitFMT".\n", aveLen);

    delete _filebuffer;
    _filebuffer    = new readBuffer(_filename, aveLen);
  }

  _isRandomAccessOpt = true;
}



//  seq == unknown sequence
//  chr == chromosomes
//  scf == scaffold
//  ctg == contig
//
void
FastAFile::printDescription(FILE *out, char *name) {

  fprintf(out, "!format ata 1.0\n");

  //  Remember the alphabet.
  //
  char    alpha[257] = {0};
  u32bit  alphalen = 0;
  for (u32bit i=0; i<256; i++)
    if (_theGlobalDesc._alphabet[i])
      alpha[alphalen++] = (char)i;
  alpha[alphalen] = 0;


  //  Print the description of these sequences as comments
  //
  fprintf(out, "!filename              = %s\n",          _filename);
  fprintf(out, "!numberOfSequences     = "u32bitFMT"\n", _theGlobalDesc._numberOfSequences);
  fprintf(out, "!fastaFileSize         = "u64bitFMT"\n", _theGlobalDesc._fastaFileSize);
  fprintf(out, "!fastaModificationTime = "s64bitFMT"\n", _theGlobalDesc._fastaModificationTime);
  fprintf(out, "!fastaCreationTime     = "s64bitFMT"\n", _theGlobalDesc._fastaCreationTime);
  fprintf(out, "!seqlineLength         = "u32bitFMT"\n", _theGlobalDesc._seqlineLength);
  fprintf(out, "!seqendlLength         = "u32bitFMT"\n", _theGlobalDesc._seqendlLength);
  fprintf(out, "!fixedWidth            = %s\n",          _theGlobalDesc._fixedWidth ? "yes" : "no");
  fprintf(out, "!squeezedSequences     = %s\n",          _theGlobalDesc._squeezedSequences ? "yes" : "no");
  fprintf(out, "!alphabet              = %s\n",          alpha);

  //  Print the same stuff on a single line
  //
  fprintf(out, "S %s %s FASTA DNA %s %s "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
          name,
          _filename,
          alpha,
          (_theGlobalDesc._squeezedSequences) ? "SQUEEZED" : 
          (_theGlobalDesc._fixedWidth) ? "FIXED" : "VARIABLE",
          _theGlobalDesc._seqlineLength,
          _theGlobalDesc._seqendlLength,
          _theGlobalDesc._numberOfSequences);

  char *names = _theNames;

  for (u32bit iid=0; iid<_theGlobalDesc._numberOfSequences; iid++) {
    fprintf(out, "G seq %s:%u . 0 1 %s 0 "u32bitFMT" "u64bitFMT" "u32bitFMT" "u32bitFMT" "u64bitFMT" %s\n",
            name,
            iid,
            name,
            _theSeqs[iid]._seqLen,
            _theSeqs[iid]._seqStart,
            iid,
            _theSeqs[iid]._headerLen,
            _theSeqs[iid]._headerStart,
            names ? names : ".");

    if (names)
      names = moveToNextName(names, iid);
  }
}
