#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>

#include "bio++.H"


fastaFile::fastaFile() {
  memset(&_theDesc,  0, sizeof(_idxfa_global));
  memset(_filename,  0, sizeof(char) * FILENAME_MAX);
  memset(_indexname, 0, sizeof(char) * FILENAME_MAX);

  _theSeqs               = 0L;
  _theNamesLen           = 0;
  _theNames              = 0L;
  _filebuffer            = 0L;
  _curIID                = 0;
  _isStreamInput         = false;
  _isIndexed             = false;
  _isRandomAccessOpt     = false;
}


fastaFile::fastaFile(char const *filename) {

  memset(&_theDesc,  0, sizeof(_idxfa_global));
  memset(_filename,  0, sizeof(char) * FILENAME_MAX);
  memset(_indexname, 0, sizeof(char) * FILENAME_MAX);

  _theSeqs               = 0L;
  _theNamesLen           = 0;
  _theNames              = 0L;
  _filebuffer            = 0L;
  _curIID                = 0;
  _isStreamInput         = false;
  _isIndexed             = false;
  _isRandomAccessOpt     = false;

  _fileTimeStamp         = time(NULL);

  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    strcpy(_filename, "stdin");

    _filebuffer    = new readBuffer(fileno(stdin), _filename);
    _isStreamInput = true;

    return;
  }

  strcpy(_filename,  filename);
  strcpy(_indexname, filename);

  //  If the filename ends in '.fasta' then append a 'idx',
  //  otherwise, append '.fastaidx'.
  //
  u32bit l = strlen(filename);
  if ((l > 5) && (strcmp(filename + l - 6, ".fasta") == 0))
    strcat(_indexname, "idx");
  else
    strcat(_indexname, ".fastaidx");

  _filebuffer    = new readBuffer(_filename);

  _fileTimeStamp = timeOfFile(_filename);
}



fastaFile::~fastaFile() {
  delete    _filebuffer;
  delete [] _theSeqs;
  delete [] _theNames;
}



seqFile *
fastaFile::openFile(const char *filename) {
  struct stat    st;

  //  We assume that if it's stdin, it's a fasta.
  //
  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    fprintf(stderr, "stdin supplied; assuming FastA format.\n");
    return(new fastaFile(filename));
  }

  errno = 0;
  stat(filename, &st);
  if (errno)
    return(0L);
  if ((st.st_mode & S_IFREG) == 0)
    return(0L);

  //  Otherwise, open and see if we can get the first sequence.  We
  //  assume it's fasta if we find a '>' denoting a defline the first
  //  thing in the file.
  //
  fastaFile *f = new fastaFile(filename);

  while ((!f->_filebuffer->eof()) && whitespaceSymbol[f->_filebuffer->get()])
    f->_filebuffer->next();

  if (f->_filebuffer->get() == '>') {
    f->_filebuffer->seek(0);
    return(f);
  }

  delete f;

  return(0L);
}



void
fastaFile::openIndex(u32bit indextypetoload) {
  struct stat    st;

  //  If we've been told to open an index, but we're not random
  //  access, complain.
  //
  if (_isStreamInput) {
    fprintf(stderr, "fastaFile()--  '%s' is a stream and not valid for indexing!\n", _filename);
    return;
  }

  //  If the index is already open, and it's the same type as desired, we
  //  can just return.  Downgrading an index is not allowed -- if you open
  //  a FASTA_INDEX_PLUS_DEFLINES, you're stuck with it.
  //
  if ((_isIndexed) &&
      (indextypetoload <= _theDesc._indexType))
    return;

  assert(_indexname[0] != 0);

  //  Stat the -> fasta <- file to get the modification time
  // 
  errno = 0;
  stat(_filename, &st);
  if (errno) {
    fprintf(stderr, "ERROR: Can't stat '%s': %s\n", _filename, strerror(errno));
    exit(1);
  }
  if (st.st_size == 0) {
    fprintf(stderr, "ERROR: source fasta file '%s' is empty!\n", _filename);
    exit(1);
  }

  int indexfile = open(_indexname, O_RDONLY | O_LARGEFILE);
  if (errno == ENOENT) {
    createIndex(indextypetoload);
    return;
  } else if (errno) {
    fprintf(stderr, "ERROR: Index file '%s' exists, but can't be opened for reading: %s\n",
            _indexname, strerror(errno));
    exit(1);
  }

  //  Opened an index, so read the description.
  //
  read(indexfile, &_theDesc, sizeof(_idxfa_global));
  if (errno) {
    fprintf(stderr, "ERROR: Index file '%s' exists, but can't read description: %s\n",
            _indexname, strerror(errno));
    exit(1);
  }

  //  Fix endianess issues in the index.
  //
  if (_theDesc._magic == u64bitSwap(FASTA_MAGICNUMBER)) {
    _theDesc._magic                 = u64bitSwap(_theDesc._magic);
    _theDesc._version               = u32bitSwap(_theDesc._version);
    _theDesc._indexType             = u32bitSwap(_theDesc._indexType);
    _theDesc._numberOfSequences     = u32bitSwap(_theDesc._numberOfSequences);
    _theDesc._fastaFileSize         = u64bitSwap(_theDesc._fastaFileSize);
    _theDesc._fastaModificationTime = u64bitSwap(_theDesc._fastaModificationTime);
    _theDesc._fastaCreationTime     = u64bitSwap(_theDesc._fastaCreationTime);
    _theDesc._seqlineLength         = u32bitSwap(_theDesc._seqlineLength);
    _theDesc._seqendlLength         = u32bitSwap(_theDesc._seqendlLength);
    _theDesc._fixedWidth            = u32bitSwap(_theDesc._fixedWidth);
    _theDesc._squeezedSequences     = u32bitSwap(_theDesc._squeezedSequences);
  }

  //  Check that the magic number and version are correct and that
  //  the modification time of the fasta file is the same as stored
  //  in the index.
  //
  if ((_theDesc._magic                 != FASTA_MAGICNUMBER) ||
      (_theDesc._version               != FASTA_VERSIONNUMBER) ||
      (_theDesc._fastaModificationTime != (u64bit)st.st_mtime) ||
      (_theDesc._indexType              < indextypetoload)) {
    const char *indexTypeNames[4] = {
      "any-index",
      "index-only",
      "index-plus-names",
      "index-plus-deflines"
    };

    fprintf(stderr, "WARNING: Index found, but stale or wrong type; I'll recreate it.\n");

    if (_theDesc._magic != FASTA_MAGICNUMBER)
      fprintf(stderr, "           Magic number incorrect; perhaps not an index file?  got="u64bitHEX" expected="u64bitHEX"\n",
              _theDesc._magic, (u64bit)FASTA_MAGICNUMBER);
    if (_theDesc._version != FASTA_VERSIONNUMBER)
      fprintf(stderr, "           Version number incorrect; got "u32bitFMT", expected "u32bitFMT".\n",
              _theDesc._version, (u32bit)FASTA_VERSIONNUMBER);
    if (_theDesc._fastaModificationTime != (u64bit)st.st_mtime)
      fprintf(stderr, "           File age not correct; got "s64bitFMT", expected %ld.\n",
              _theDesc._fastaModificationTime, st.st_mtime);
    if (_theDesc._indexType < indextypetoload)
      fprintf(stderr, "           Type of index insufficient; got %s, need %s.\n",
              indexTypeNames[_theDesc._indexType],
              indexTypeNames[indextypetoload]);

    createIndex(indextypetoload);
    return;
  }

  //  Continue reading the index

  _theSeqs = new _idxfa_desc [ _theDesc._numberOfSequences ];

  read(indexfile, _theSeqs, sizeof(_idxfa_desc) * _theDesc._numberOfSequences);
  if (errno) {
    fprintf(stderr, "fastaFile()-- couldn't read sequence descriptions from the index '%s': %s\n", _indexname, strerror(errno));
    exit(1);
  }

  if (_theDesc._magic == u64bitSwap(FASTA_MAGICNUMBER)) {
    for (u32bit i=0; i<_theDesc._numberOfSequences; i++) {
      _theSeqs[i]._headerStart = u64bitSwap(_theSeqs[i]._headerStart);
      _theSeqs[i]._seqStart    = u64bitSwap(_theSeqs[i]._seqStart);
      _theSeqs[i]._headerLen   = u32bitSwap(_theSeqs[i]._headerLen);
      _theSeqs[i]._seqLen      = u32bitSwap(_theSeqs[i]._seqLen);
    }
  }


  //  If indextype is FASTA_INDEX_ANY, open the index as reported by
  //  the file.  Otherwise, open the index as specified.
  //
  if (indextypetoload == FASTA_INDEX_ANY)
    indextypetoload = _theDesc._indexType;
  else
    _theDesc._indexType = indextypetoload;


  if ((_theDesc._indexType == FASTA_INDEX_PLUS_IDS) ||
      (_theDesc._indexType == FASTA_INDEX_PLUS_DEFLINES)) {
    read(indexfile, &_theNamesLen, sizeof(u32bit));
    if (errno) {
      fprintf(stderr, "fastaFile()-- couldn't read lengths of names from the index '%s': %s\n", _indexname, strerror(errno));
      exit(1);
    }

    if (_theDesc._magic == u64bitSwap(FASTA_MAGICNUMBER))
      _theNamesLen = u32bitSwap(_theNamesLen);

    _theNames = new char [sizeof(char) * _theNamesLen];

    read(indexfile, _theNames, sizeof(char) * _theNamesLen);
    if (errno) {
      fprintf(stderr, "fastaFile()-- couldn't read names from the index '%s': %s\n", _indexname, strerror(errno));
      exit(1);
    }
  }

  close(indexfile);

  _isIndexed = true;
}


seqInCore*
fastaFile::getSequenceInCore(void) {
  seqInCore *c = 0L;

  //  Skip whitespace at the start of the sequence.
  //
  while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
    _filebuffer->next();

  if (_filebuffer->eof())
    return(0L);

  //  We should be at a '>' character now.  Fail if not.
  //
  if (_filebuffer->get() != '>') {
    fprintf(stderr, "getSequenceInCore()-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, _filebuffer->get());
    exit(1);
  }

  //
  //  If we have an index, we can be quicker about reading things.
  //    If squeezed, suck in the whole thing in two reads
  //    If fixed width, suck in whole lines
  //    Else do character by character
  //

  if        ((isIndexed() && _theDesc._squeezedSequences) ||
             (isIndexed() && _theDesc._fixedWidth)) {

    //  Yea!  We have an index!  Preallocate the space, and do some quick
    //  block reads!
    //
    char *h = new char [ _theSeqs[_curIID]._headerLen + 1 ];
    char *s = new char [ _theSeqs[_curIID]._seqLen + 1 ];

    _filebuffer->read(h, _theSeqs[_curIID]._headerLen);

    //  Skip any whitespace between the defline and the start of the sequence
    //  (Rather than a seek, we just skip the spaces)
    //
    while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
      _filebuffer->next();

    if (_theDesc._squeezedSequences) {
      _filebuffer->read(s, _theSeqs[_curIID]._seqLen);
    } else {
      u32bit  bytesRead = 0;
      while (bytesRead < _theSeqs[_curIID]._seqLen) {

        u32bit  lengthToRead = _theDesc._seqlineLength;

        if (bytesRead + lengthToRead > _theSeqs[_curIID]._seqLen)
          lengthToRead = _theSeqs[_curIID]._seqLen - bytesRead;

        bytesRead += (u32bit)_filebuffer->read(s + bytesRead, lengthToRead);

        //  Skip any whitespace between the defline and the start of the sequence
        //  (Rather than a seek, we just skip the spaces)
        //
        while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
          _filebuffer->next();
      }
    }

    h[_theSeqs[_curIID]._headerLen] = 0;
    s[_theSeqs[_curIID]._seqLen] = 0;

    c = new seqInCore(_curIID,
                      h, _theSeqs[_curIID]._headerLen,
                      s, _theSeqs[_curIID]._seqLen);
  } else {
    //  Shucks!  We know nothing about the file.  Go character by
    //  character.
    //
    u32bit hLen = 0;
    u32bit hMax = 128;
    u32bit sLen = 0;
    u32bit sMax = 16 * 1024 * 1024;

    //  If we have an index, we know the sizes.  The file must not
    //  be fixed width (or squeezed).
    //
    if (isIndexed()) {
      hMax = _theSeqs[_curIID]._headerLen + 1;
      sMax = _theSeqs[_curIID]._seqLen + 1;
    }

    char   *h   = new char [ hMax + 1 ];
    char   *s   = new char [ sMax + 1 ];

    while ((!_filebuffer->eof()) && (_filebuffer->get() != '\r') && (_filebuffer->get() != '\n')) {
      if (hLen >= hMax) {
        hMax += 128;
        char *htmp = new char [hMax + 1];
        memcpy(htmp, h, sizeof(char) * hLen);
        delete [] h;
        h = htmp;
      }
      h[hLen++] = _filebuffer->get();

      _filebuffer->next();
    }
    h[hLen] = 0;

    //  Skip any whitespace between the defline and the start of the sequence
    //
    while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
      _filebuffer->next();

    while ((!_filebuffer->eof()) &&
           (_filebuffer->get() != '>')) {
      if (!whitespaceSymbol[_filebuffer->get()]) {
        if (sLen >= sMax) {
          sMax += 32 * 1024 * 1024;
          char *stmp = new char [sMax + 1];
          memcpy(stmp, s, sizeof(char) * sLen);
          delete [] s;
          s = stmp;
        }
        s[sLen++] = _filebuffer->get();
      }

      _filebuffer->next();
    }
    s[sLen] = 0;

    c = new seqInCore(_curIID,
                      h, hLen, s, sLen);
  }


  //  Skip whitespace at the end of the sequence.
  //
  while ((!_filebuffer->eof()) && whitespaceSymbol[_filebuffer->get()])
    _filebuffer->next();


  //  By reading the sequence, we have implicitly moved to the next
  //  sequence
  //
  _curIID++;

  return(c);
}


seqOnDisk*
fastaFile::getSequenceOnDisk(void) {

  if (_curIID >= _theDesc._numberOfSequences)
    return(0L);

  //  XXX: A find has already been performed.  This is probably wasted
  //  (it loads the readbuffer).  We should pass the read buffer into
  //  the SeqOnDisk, and allocate a new one for the wrapper.

  seqOnDisk  *s = new seqOnDisk(getSourceName(),
                                _theSeqs[_curIID]._headerStart,
                                _theSeqs[_curIID]._headerLen,
                                _theSeqs[_curIID]._seqStart,
                                _theSeqs[_curIID]._seqLen,
                                _curIID,
                                _theDesc._squeezedSequences,
                                _theDesc._fixedWidth,
                                _theDesc._seqlineLength,
                                _theDesc._seqendlLength);

  //  To be compatible with the SequenceInCore, we need to advance to
  //  the next sequence.
  //
  _curIID++;

  return(s);
}



bool
fastaFile::find(seqIID  iid) {

  if (!_isIndexed) {
    fprintf(stderr, "fastaFile::find(IID)-- ERROR: '%s' is not open for random access; find() failed.\n", _filename);
    return(false);
  }

  if (iid >= _theDesc._numberOfSequences) {
    fprintf(stderr, "fastaFile::find(IID)-- ERROR: index of "u32bitFMT" too large for '%s' (only "u32bitFMT" sequences).\n",
            iid, _filename, _theDesc._numberOfSequences);
    return(false);
  }

  if (iid != _curIID) {
    _filebuffer->seek(_theSeqs[iid]._headerStart);
    _curIID = iid;
  }

  return(true);
}



bool
fastaFile::find(char *id) {
  u32bit         iid = 0;
  char          *ptr = _theNames;
  char          *res = 0L;

  if (id == 0L)
    return(false);

  //  XXX: probably shoujld trim off the whitespace from the front/end
  //  of the id

  while (iid < _theDesc._numberOfSequences) {
    res = strstr(ptr, id);
    if (res) {
      _filebuffer->seek(_theSeqs[iid]._headerStart);
      _curIID = iid;
      return(true);
    }

    //  Move to the next name.  Sadly, if our index has saved only the
    //  id, we don't know the length of this field.  We must blindly
    //  search for the start of the next id.
    //
    if (_theDesc._indexType == FASTA_INDEX_PLUS_IDS) {
      while (*ptr)
        ptr++;
      ptr++;
    } else {
      ptr += _theSeqs[iid]._headerLen + 1;
    }

    iid++;
  }
    
  return(false);
}



//  seq == unknown sequence
//  chr == chromosomes
//  scf == scaffold
//  ctg == contig
//
void
fastaFile::printDescription(FILE *out, char *name) {

  fprintf(out, "!format ata 1.0\n");

  //  Remember the alphabet.
  //
  char    alpha[257] = {0};
  u32bit  alphalen = 0;
  for (u32bit i=0; i<256; i++)
    if (_theDesc._alphabet[i])
      alpha[alphalen++] = (char)i;
  alpha[alphalen] = 0;


  //  Print the description of these sequences as comments
  //
  fprintf(out, "!filename              = %s\n",          _filename);
  fprintf(out, "!numberOfSequences     = "u32bitFMT"\n", _theDesc._numberOfSequences);
  fprintf(out, "!fastaFileSize         = "u64bitFMT"\n", _theDesc._fastaFileSize);
  fprintf(out, "!fastaModificationTime = "s64bitFMT"\n", _theDesc._fastaModificationTime);
  fprintf(out, "!fastaCreationTime     = "s64bitFMT"\n", _theDesc._fastaCreationTime);
  fprintf(out, "!seqlineLength         = "u32bitFMT"\n", _theDesc._seqlineLength);
  fprintf(out, "!seqendlLength         = "u32bitFMT"\n", _theDesc._seqendlLength);
  fprintf(out, "!fixedWidth            = %s\n",          _theDesc._fixedWidth ? "yes" : "no");
  fprintf(out, "!squeezedSequences     = %s\n",          _theDesc._squeezedSequences ? "yes" : "no");
  fprintf(out, "!alphabet              = %s\n",          alpha);

  //  Print the same stuff on a single line
  //
  fprintf(out, "S %s %s FASTA DNA %s %s "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
          name,
          _filename,
          alpha,
          (_theDesc._squeezedSequences) ? "SQUEEZED" : 
          (_theDesc._fixedWidth) ? "FIXED" : "VARIABLE",
          _theDesc._seqlineLength,
          _theDesc._seqendlLength,
          _theDesc._numberOfSequences);

  char *names = _theNames;

  for (u32bit iid=0; iid<_theDesc._numberOfSequences; iid++) {
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

    if (names) {
      //  Move to the next name.  Sadly, if our index has saved only the
      //  id, we don't know the length of this field.  We must blindly
      //  search for the start of the next id.
      //
      if (_theDesc._indexType == FASTA_INDEX_PLUS_IDS) {
        while (*names)
          names++;
        names++;
      } else {
        names += _theSeqs[iid]._headerLen + 1;
      }
    }
  }
}




void
fastaFile::createIndex(u32bit indextype) {

  //  If the user didn't specify which type of index to build, build
  //  just the basic one.
  //
  if (indextype == FASTA_INDEX_ANY)
    indextype = FASTA_INDEX_ONLY;

  //  Clear the description structure
  //
  memset(&_theDesc, 0, sizeof(_idxfa_global));

  _theDesc._magic                 = FASTA_MAGICNUMBER;
  _theDesc._version               = FASTA_VERSIONNUMBER;
  _theDesc._indexType             = indextype;
  _theDesc._numberOfSequences     = 0;
  _theDesc._fastaFileSize         = 0;
  _theDesc._fastaModificationTime = 0;
  _theDesc._fastaCreationTime     = 0;
  _theDesc._seqlineLength         = 0;
  _theDesc._seqendlLength         = 0;
  _theDesc._fixedWidth            = true;
  _theDesc._squeezedSequences     = true;

  //  Remove any existing index -- this will reset permissions if they
  //  are hosed, and generally clean the slate for us.
  //
  errno = 0;
  unlink(_indexname);
  if ((errno) && (errno != ENOENT))
    fprintf(stderr, "fastaFile::buildIndex()-- Can't remove old index '%s' (but I'll rebuild it anyway): %s\n",
            _indexname, strerror(errno));


  //  Before we compute the index, try to open the output file.  This
  //  saves us the compute if we cannot actually open the file
  //  (read-only directory, usually).
  //
  errno = 0;
  int indexfile = open(_indexname,
                       O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                       S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno) {
    fprintf(stderr, "fastaFile::buildIndex()-- Can't open index file '%s' for writing: %s\n",
            _indexname, strerror(errno));
    exit(1);
  }


  //  Stat the file to get size and times.  Also checks if the file exists.
  //
  struct stat  fastastat;
  errno = 0;
  if (stat(_filename, &fastastat)) {
    fprintf(stderr, "FastA::buildIndex()-- stat() of file '%s' failed.\n%s\n",
            _filename, strerror(errno));
    exit(1);
  }
  _theDesc._fastaFileSize         = fastastat.st_size;
  _theDesc._fastaModificationTime = fastastat.st_mtime;
  _theDesc._fastaCreationTime     = fastastat.st_ctime;


  //  Describes the sequences themselves
  //
  _numSeqs = 0;
  _maxSeqs = 1024;
  _theSeqs = new _idxfa_desc [_maxSeqs];

  //  Copies of the deflines or names
  //
  _theNamesLen = 0;
  _theNamesMax = 16 * 1024;
  _theNames    = new char [_theNamesMax];

  //  Info about the sequence we just read
  //
  u32bit    defLen    = 0;
  u64bit    defStart  = 0;
  u32bit    seqLen    = 0;
  u64bit    seqStart  = 0;

  //  Copy of the defline for the sequence we just read
  //
  u32bit    theDefLineLen    = 0;
  u32bit    theDefLineMax    = 2 * 1024;
  char     *theDefLine       = new char [theDefLineMax];

  //  Open a buffered input sequence
  //
  readBuffer    B(_filename);

  //  We need to remember the length of the first sequence line
  //  so we can test if all sequence lines are the same length.
  //  This could be extended to be per-sequence, but that makes
  //  the index 25% larger, and doesn't (IMHO) gain much.
  //
  _theDesc._seqlineLength = ~u32bitZERO;

  while (!B.eof()) {

    //  Skip any whitespace before the defline
    //
    while ((!B.eof()) && whitespaceSymbol[B.get()])
      B.next();

    //  We should be at a '>' character now.  Fail if not.
    //
    if (B.get() != '>') {
      fprintf(stderr, "FastA::buildIndex()-- In file %s, expected '>' at beginning of defline, got '%c' instead.\n",
              _filename, B.get());
      fprintf(stderr, "FastA::buildIndex()-- File position is "s64bitFMT", sequence number "u32bitFMT"\n", (s64bit)B.tell(), _numSeqs);

      exit(1);
    }

    //  Save the start of the defline
    //
    defStart = B.tell();
    defLen   = 0;

    //  Count the length of the defline -- counts until the first \r or \n
    //
    //  Also saves a copy of the characters if 'saveDefLine' is enabled.  This
    //  duals as a method of saving just the first word.
    //
    //  theDefLine is NOT zero terminated!
    //
    bool saveDefLine = true;
    if (_theDesc._indexType == FASTA_INDEX_ONLY)
      saveDefLine = false;

    theDefLineLen = 0;

    while ((!B.eof()) && (B.get() != '\r') && (B.get() != '\n')) {

      if ((_theDesc._indexType == FASTA_INDEX_PLUS_IDS) &&
          (whitespaceSymbol[B.get()]))
        saveDefLine = false;

      if (saveDefLine) {
        if (theDefLineLen > theDefLineMax) {
          theDefLineMax *= 2;
          char *newdef = new char [theDefLineMax];
          memcpy(newdef, theDefLine, theDefLineLen);
          delete [] theDefLine;
          theDefLine = newdef;
        }

        theDefLine[theDefLineLen++] = B.get();
      }

      defLen++;

      B.next();
    }

    //  Skip any whitespace between the defline and the start of the sequence
    //
    while ((!B.eof()) && whitespaceSymbol[B.get()])
      B.next();

    //  Save the start of the sequence
    //
    seqStart = B.tell();
    seqLen   = 0;

    //  We'll compare thisLineLen to firstLineLen.  If they differ,
    //  then the global _fixedWidth flag is cleared.
    //
    u32bit   thisLineLen      = 0;      //  length of the line we've read in, including whitespace
    u32bit   thisLineSpaces   = 0;      //  How many spaces before hitting end of line?

    bool     lastLineDisagree = false;
    bool     multipleLines    = false;

    while ((!B.eof()) &&
           (B.get() != '>')) {

      if (!whitespaceSymbol[B.get()]) {
        seqLen++;
        thisLineLen++;
        _theDesc._alphabet[B.get()] = 1;

        //  If we've seen space already, then we have embedded space,
        //  and we're not fixed width or squeezed.
        //
        if (thisLineSpaces) {
          _theDesc._fixedWidth        = false;
          _theDesc._squeezedSequences = false;
        }

        if (multipleLines)
          _theDesc._squeezedSequences = false;

        B.next();

      } else {

        //  Count the number of spaces we see.  If we get more
        //  sequence before the end of line, that's bad.  See above.
        //
        thisLineSpaces++;

        //  If our space is a newline, then set the multipleLines
        //  flag.  If we see more sequence in this sequence, then we
        //  are not squeezed.
        //
        if ((B.get() == '\r') ||
            (B.get() == '\n')) {
          multipleLines = true;

          //  First, count any remaining whitespace on this "line"
          //  (extra linebreak characters, etc).  This is used much
          //  the same as seqlineLength.
          //
          B.next();
          while ((!B.eof()) && (whitespaceSymbol[B.get()])) {
            thisLineSpaces++;
            B.next();
          }

          //  Check the line lengths.
          //
          //  0) if the length of the first line in the file isn't
          //  set, set it.
          //
          //  1) if the *last* line we read in was of a different
          //  length than the first, then set the non-fixed width flag
          //
          //  2) if this line length differs with the first line
          //  length, set the lastLineDisagree flag.  We'll update the
          //  global flag if this turns out to be not the last line.
          //
          if (_theDesc._seqlineLength == ~u32bitZERO) {
            _theDesc._seqlineLength = thisLineLen;
            _theDesc._seqendlLength = thisLineSpaces;
          }

          if (lastLineDisagree)
            _theDesc._fixedWidth = false;

          if ((thisLineLen != _theDesc._seqlineLength) ||
              (thisLineSpaces != _theDesc._seqendlLength))
            lastLineDisagree = true;

          //  Reset this line length to zero; we're done with it.
          //
          thisLineLen      = 0;
          thisLineSpaces   = 0;
        } else {
          B.next();
        }
      }
    }

    //  Make the index bigger, if needed
    //
    if (_numSeqs >= _maxSeqs) {
      if (_maxSeqs == 0)
        _maxSeqs = 16;
      _maxSeqs *= 2;

      _idxfa_desc *sa = new _idxfa_desc [ _maxSeqs ];
      memcpy(sa, _theSeqs, sizeof(_idxfa_desc) * _numSeqs);
      delete [] _theSeqs;
      _theSeqs = sa;
    }

    //  Add the new sequence description to the list.
    //
    _theSeqs[_numSeqs]._headerStart  = defStart;
    _theSeqs[_numSeqs]._headerLen    = defLen;
    _theSeqs[_numSeqs]._seqStart     = seqStart;
    _theSeqs[_numSeqs]._seqLen       = seqLen;

    //  Add the description of the sequence to the list of descriptions
    //
    if ((_theDesc._indexType == FASTA_INDEX_PLUS_IDS) ||
        (_theDesc._indexType == FASTA_INDEX_PLUS_DEFLINES)) {
      if (_theNamesLen + defLen >= _theNamesMax) {
        _theNamesMax *= 2;
        char *nd = new char [_theNamesMax];
        memcpy(nd, _theNames, sizeof(char) * _theNamesLen);
        delete [] _theNames;
        _theNames = nd;
      }
      memcpy(_theNames + _theNamesLen, theDefLine, theDefLineLen);
      _theNamesLen += theDefLineLen;
      _theNames[_theNamesLen++] = 0;
    }

    _numSeqs++;
  }


  //  All done
  //
  _theDesc._numberOfSequences = _numSeqs;


  //  XXX:
  //
  //  If the sequences are all on one line (e.g., squeezed) then the
  //  method used to check for lines of the same length fails.
  //  Ideally, we should now check all the _seqLen's in the index to
  //  decide if they are really fixed length or not......but why?
  //
  if (_theDesc._squeezedSequences)
    _theDesc._fixedWidth = false;



  //  Write the data file; version, number of sequences and sequence
  //  descriptions.
  //
  errno = 0;
  write(indexfile, &_theDesc, sizeof(_idxfa_global));
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't write header to index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  errno = 0;
  write(indexfile, _theSeqs, sizeof(_idxfa_desc) * _numSeqs);
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't write index to index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  if ((_theDesc._indexType == FASTA_INDEX_PLUS_IDS) ||
      (_theDesc._indexType == FASTA_INDEX_PLUS_DEFLINES)) {
    errno = 0;
    write(indexfile, &_theNamesLen,  sizeof(u32bit));
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() can't write nameslen to index file '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }

    errno = 0;
    write(indexfile, _theNames,  sizeof(char) * _theNamesLen);
    if (errno) {
      fprintf(stderr, "FastA::buildIndex() can't write names to index file '%s'.\n%s\n", _filename, strerror(errno));
      exit(1);
    }
  }


  //  Close the index file.
  //
  errno = 0;
  close(indexfile);
  if (errno) {
    fprintf(stderr, "FastA::buildIndex() can't close index file '%s'.\n%s\n", _filename, strerror(errno));
    exit(1);
  }

  _isIndexed = true;
}



//  Analyze the index, reset the bufferSize of the readBuffer if the
//  sequences are short (e.g., ESTs, SNPs).  This will make a huge
//  difference in the completely random access case -- instead of
//  always reading 32k (default buffer size) for each sequence, we
//  read a little more than average.  But, it sucks for sequential
//  -- we'll be doing a read on almost every sequence.
//
void
fastaFile::optimizeRandomAccess(void) {

  if (_isRandomAccessOpt)
    return;

  if (_isIndexed == false)
    openIndex();

  if (_theDesc._numberOfSequences < 10)
    return;

  u64bit  aveLen = 0;
  for (u32bit i=0; i<_theDesc._numberOfSequences; i++)
    aveLen += _theSeqs[i]._seqLen + _theSeqs[i]._headerLen;
  aveLen /= _theDesc._numberOfSequences;

  u64bit stdDev = 0;
  for (u32bit i=0; i<_theDesc._numberOfSequences; i++)
    stdDev += ((_theSeqs[i]._seqLen + _theSeqs[i]._headerLen - aveLen) *
               (_theSeqs[i]._seqLen + _theSeqs[i]._headerLen - aveLen));
  stdDev /= _theDesc._numberOfSequences - 1;

  stdDev = (u64bit)ceil(sqrt((double)stdDev));

  fprintf(stderr, "fastaFile::optimizeRandomAccess()-- For "u32bitFMT" seqs, ave="u64bitFMT" stddev="u64bitFMT", reset buffer to "u64bitFMT"\n",
          _theDesc._numberOfSequences, aveLen, stdDev, aveLen + stdDev);

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

    fprintf(stderr, "fastaFile::optimizeRandomAccess()-- Make new filebuffer of size "u64bitFMT".\n", aveLen);

    delete _filebuffer;
    _filebuffer    = new readBuffer(_filename, aveLen);
  }

  _isRandomAccessOpt = true;
}
