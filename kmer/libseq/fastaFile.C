#include "fastaFile.H"
#include "alphabet.h"


#undef DEBUG
#undef DEBUGINDEX

//  Says 'kmerFastaFileIdx'
#define FASTA_MAGICNUMBER1  0x7473614672656d6bULL
#define FASTA_MAGICNUMBER2  0x786449656c694661ULL


fastaFile::fastaFile(const char *filename) {
  clear();

#ifdef DEBUG
  fprintf(stderr, "fastaFile::fastaFile()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  strcpy(_filename, filename);

  constructIndex();

  _rb    = new readBuffer(_filename);

  _numberOfSequences = _header._numberOfSequences;
}



fastaFile::fastaFile() {
  clear();
}



fastaFile::~fastaFile() {
  delete    _rb;
  delete [] _index;
  delete [] _names;
}





seqFile *
fastaFile::openFile(const char *filename) {
  struct stat    st;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::openFile()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  if (((filename == 0L) && (isatty(fileno(stdin)) == 0)) ||
      ((filename != 0L) && (filename[0] == '-') && (filename[1] == 0)))
    return(0L);

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
  //  Use of a readBuffer here is a bit heavyweight, but it's safe and
  //  easy.  Opening a fastaFile isn't, after all, lightweight anyway.
  //
  fastaFile   *f = 0L;
  readBuffer  *r = new readBuffer(filename);
  char         x = r->read();

  while ((r->eof() == false) && (whitespaceSymbol[x] == true))
    x = r->read();

  //  If we get a fasta record separator assume it's a fasta file.  If
  //  it's eof, the file is empty, and we might as well return this
  //  fasta file and let the client deal with the lack of sequence.
  //
  if ((x == '>') || (r->eof() == true))
    f = new fastaFile(filename);

  delete r;

  return(f);
}



u32bit
fastaFile::find(const char *sequencename) {
  char   *ptr = _names;

  //  If this proves far too slow, rewrite the _names string to
  //  separate IDs with 0xff, then use strstr on the whole thing.  To
  //  find the ID, scan down the string counting the number of 0xff's.
  //
  //  Similar code is used for seqStore::find()

  for (u32bit iid=0; iid < _header._numberOfSequences; iid++) {
    //fprintf(stderr, "fastaFile::find()-- '%s' vs '%s'\n", sequencename, ptr);
    if (strcmp(sequencename, ptr) == 0)
      return(iid);

    while (*ptr)
      ptr++;
    ptr++;
  }

  return(~u32bitZERO);
}



u32bit
fastaFile::getSequenceLength(u32bit iid) {

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequenceLength()-- "u32bitFMT"\n", iid);
#endif

  return((iid < _numberOfSequences) ? _index[iid]._seqLength : 0);
}



bool
fastaFile::getSequence(u32bit iid,
                       char *&h, u32bit &hLen, u32bit &hMax,
                       char *&s, u32bit &sLen, u32bit &sMax) {

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequence(full)-- "u32bitFMT"\n", iid);
#endif

  //  Assume there is no index.  Without being horribly complicated
  //  (as in the previous versions of this codebase) all we'd get from
  //  having an index around is the length of the sequence.
  //
  //  Previous versions used to use the index to tell if the sequence
  //  was squeezed (and so a direct copy to the output), if it was
  //  fixed width (mostly direct copies) or unknown.  Now we just
  //  assume it's unknown and go byte by byte.  If speed is a concern,
  //  use the seqFile instead.

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "fastaFile::getSequence(full)--  iid "u32bitFMT" more than number of sequences "u32bitFMT"\n",
      iid, _header._numberOfSequences);
    return(false);
  }

  if (sMax == 0) {
    sMax = 2048;
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 2048;
    h    = new char [hMax];
  }

  if ((_index) && (sMax < _index[iid]._seqLength)) {
    sMax = _index[iid]._seqLength;
    delete [] s;
    s = new char [sMax];
  }

  hLen = 0;
  sLen = 0;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequence(full)-- seek to iid="u32bitFMT" at pos="u32bitFMT"\n",
          iid, _index[iid]._seqPosition);
#endif
  _rb->seek(_index[iid]._seqPosition);

  char x = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof())
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaFile::getSequence(full)-- ERROR1: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

  //  Skip the '>' in the defline
  x = _rb->read();

  //  Skip whitespace between the '>' and the defline
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true) && (x != '\r') && (x != '\n'))
    x = _rb->read();

  //  Copy the defline, until the first newline.
  while ((_rb->eof() == false) && (x != '\r') && (x != '\n')) {
    h[hLen++] = x;
    if (hLen >= hMax) {
      hMax += 2048;
      char *H = new char [hMax];
      memcpy(H, h, hLen);
      delete [] h;
      h = H;
    }
    x = _rb->read();
  }
  h[hLen] = 0;

  //  Skip whitespace between the defline and the sequence.
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  Copy the sequence, until EOF or the next '>'.
  while ((_rb->eof() == false) && (_rb->peek() != '>')) {
    if (whitespaceSymbol[x] == false) {
      s[sLen++] = x;
      if (sLen >= sMax) {
        if (sMax == 4294967295)  //  4G - 1
          fprintf(stderr, "fastaFile::getSequence()-- ERROR: sequence is too long; must be less than 4 Gbp.\n"), exit(1);
        if (sMax >= 2147483648)  //  2G
          sMax = 4294967295;
        else
          sMax *= 2;
        char *S = new char [sMax];
        memcpy(S, s, sLen);
        delete [] s;
        s = S;
      }
    }
    x = _rb->read();
  }
  s[sLen] = 0;

  _nextID++;

  return(true);
}


// slow
bool
fastaFile::getSequence(u32bit iid,
                       u32bit bgn, u32bit end, char *s) {

  if (iid >= _header._numberOfSequences) {
    fprintf(stderr, "fastaFile::getSequence(part)--  iid "u32bitFMT" more than number of sequences "u32bitFMT"\n",
      iid, _header._numberOfSequences);
    return(false);
  }

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequence(part)-- "u32bitFMT"\n", iid);
#endif

  //  It is impossible to be efficient here; see the big comment in
  //  the other getSequence() above.
  //
  //  We can't even guess where to start scanning the sequence; we
  //  just don't have any information about how much whitespace is in
  //  the sequence.

  _rb->seek(_index[iid]._seqPosition);

  u32bit pos = 0;
  char   x   = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof())
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaFile::getSequence(part)-- ERROR2: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

  //  Skip the defline.
  while ((_rb->eof() == false) && (x != '\r') && (x != '\n'))
    x = _rb->read();

  //  Skip whitespace between the defline and the sequence.
  while ((_rb->eof() == false) && (whitespaceSymbol[x] == true))
    x = _rb->read();

  //  Skip sequence up until bgn.
  while ((_rb->eof() == false) && (pos < bgn)) {
    if (whitespaceSymbol[x] == false)
      pos++;

    x = _rb->read();
  }

  //  Copy sequence
  while ((_rb->eof() == false) && (pos < end)) {
    if (whitespaceSymbol[x] == false)
      s[pos++ - bgn] = x;

    x = _rb->read();
  }
  s[pos - bgn] = 0;

  //  Fail if we didn't copy enough stuff.
  return((pos == end) ? true : false);
}




void
fastaFile::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "FastA");

  _numberOfSequences = 0;

  _rb                = 0L;
  memset(&_header, 0, sizeof(fastaFileHeader));
  _index = 0L;
  _names = 0L;
  _nextID = 0;
}



void
fastaFile::loadIndex(char *indexname) {
  struct stat  fastastat;

  if (fileExists(indexname) == false)
    return;

  errno = 0;
  if (stat(_filename, &fastastat)) {
    fprintf(stderr, "fastaFile::constructIndex()-- stat of file '%s' failed: %s\n",
            _filename, strerror(errno));
    return;
  }

  FILE *I = fopen(indexname, "r");
  if (errno) {
    fprintf(stderr, "fastaFile::constructIndex()-- open of file '%s' failed: %s\n",
            indexname, strerror(errno));
    return;
  }

  fread(&_header, sizeof(fastaFileHeader), 1, I);

  if ((_header._magic[0] != FASTA_MAGICNUMBER1) &&
      (_header._magic[1] != FASTA_MAGICNUMBER2)) {
    fprintf(stderr, "fastaFile::constructIndex()-- magic mismatch.\n");
    fclose(I);
    return;
  }

  if ((_header._fastaFileSize         != (u64bit)fastastat.st_size) ||
      (_header._fastaModificationTime != (u64bit)fastastat.st_mtime) ||
      (_header._fastaCreationTime     != (u64bit)fastastat.st_ctime)) {
    fprintf(stderr, "fastaFile::constructIndex()-- stat mismatch.\n");
    fclose(I);
    return;
  }

  _index = new fastaFileIndex [_header._numberOfSequences];
  _names = new char           [_header._namesLength];

  fread(_index, sizeof(fastaFileIndex), _header._numberOfSequences, I);
  fread(_names, sizeof(char),           _header._namesLength,       I);

#ifdef DEBUG
  fprintf(stderr, "fastaFile::constructIndex()-- '%s' LOADED\n", _filename);
#endif

  fclose(I);
  return;
}


void
fastaFile::constructIndex(void) {

  if (_index)
    return;

  //  If the filename ends in '.fasta' then append a 'idx',
  //  otherwise, append '.fastaidx'.

  char  indexname[FILENAME_MAX];

  strcpy(indexname, _filename);
  u32bit l = strlen(_filename);
  if ((l > 5) && (strcmp(_filename + l - 6, ".fasta") == 0))
    strcat(indexname, "idx");
  else
    strcat(indexname, ".fastaidx");

  //  If the index exists, suck it in and return.

  loadIndex(indexname);

  if (_index)
    return;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::constructIndex()-- '%s' BUILDING\n", _filename);
#endif

  //  Allocate some space for the index structures.

  u32bit  indexMax = 64 * 1024 * 1024 / sizeof(fastaFileIndex);
  u32bit  indexLen = 0;

  _index = new fastaFileIndex [indexMax];

  u32bit  namesMax = 32 * 1024 * 1024;
  u32bit  namesLen = 0;

  _names = new char [namesMax];

  //  Some local storage

  u64bit       seqStart;
  u32bit       seqLen;
  u32bit       seqLenMax = ~u32bitZERO;
  u32bit       namePos;

  readBuffer   ib(_filename);
  char         x = ib.read();

#ifdef DEBUGINDEX
  fprintf(stderr, "readBuffer '%s' eof=%d x=%c %d\n", _filename, ib.eof(), x, x);
#endif

  //  Build it.

  //  Skip whitespace at the start of the sequence.
  while ((ib.eof() == false) && (whitespaceSymbol[x] == true)) {
#ifdef DEBUGINDEX
    fprintf(stderr, "skip '%c' %d\n", x, x);
#endif
    x = ib.read();
  }

  while (ib.eof() == false) {
#ifdef DEBUGINDEX
    fprintf(stderr, "index\n");
#endif

    //  We should be at a '>' character now.  Fail if not.
    if (x != '>')
      fprintf(stderr, "fastaFile::constructIndex()-- ERROR3: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
              _filename, x), exit(1);

    //  Save info - ib's position is correctly at the first letter in
    //  the defline (which might be whitespace), but the reader
    //  expects our position to be at the '>' -- hence the -1.
    seqStart = ib.tell() - 1;
    seqLen   = 0;
    namePos  = namesLen;

    //  Read that first letter
    x = ib.read();

    //  Copy the name to the names
    while ((ib.eof() == false) && (whitespaceSymbol[x] == false)) {
      if (namesLen + 1 >= namesMax) {
        namesMax += 32 * 1024 * 1024;
        char *nt = new char [namesMax];
        memcpy(nt, _names, namesLen);
        delete [] _names;
        _names = nt;
      }

      _names[namesLen++] = x;
#ifdef DEBUGINDEX
      fprintf(stderr, "name += %c\n", x);
#endif
      x = ib.read();
    }

    if (namesLen + 1 >= namesMax) {
      namesMax += 32 * 1024 * 1024;
      char *nt = new char [namesMax];
      memcpy(nt, _names, namesLen);
      delete [] _names;
      _names = nt;
    }
    _names[namesLen++] = 0;

    //  Skip the rest of the defline
    while ((ib.eof() == false) && (x != '\r') && (x != '\n')) {
#ifdef DEBUGINDEX
      fprintf(stderr, "skip let %c\n", x);
#endif
      x = ib.read();
    }

    //  Skip whitespace between the defline and the sequence.
    while ((ib.eof() == false) && (whitespaceSymbol[x] == true)) {
#ifdef DEBUGINDEX
      fprintf(stderr, "skip num %d\n", x);
#endif
      x = ib.read();
    }

#ifdef DEBUGINDEX
    fprintf(stderr, "x=%c peek=%c\n", x, ib.peek());
#endif

    //  Count sequence length
    while ((ib.eof() == false) && (ib.peek() != '>')) {
#ifdef DEBUGINDEX
      fprintf(stderr, "seqlen %s %c\n", (whitespaceSymbol[x] == false) ? "save" : "skip", x);
#endif
      if (whitespaceSymbol[x] == false)
        seqLen++;
      if (seqLen >= seqLenMax)
        fprintf(stderr, "fastaFile::constructIndex()-- ERROR: In %s, sequence '%s' is too long.  Maximum length is %u bases.\n",
                _filename, _names + namePos, seqLenMax), exit(1);
      x = ib.read();
    }

    //  Save to the index.

    if (indexLen >= indexMax) {
      indexMax *= 2;
      fastaFileIndex *et = new fastaFileIndex[indexMax];
      memcpy(et, _index, sizeof(fastaFileIndex) * indexLen);
      delete [] _index;
      _index = et;
    }

    _index[indexLen]._seqPosition = seqStart;
    _index[indexLen]._seqLength   = seqLen;

#ifdef DEBUG
    fprintf(stderr, "INDEX iid="u32bitFMT" len="u32bitFMT" pos="u64bitFMT"\n",
            indexLen, seqLen, seqStart);
#endif

    indexLen++;

    //  Load the '>' for the next iteration.
    x = ib.read();
  }

  //  Fill out the index meta data

  struct stat  fastastat;
  errno = 0;
  if (stat(_filename, &fastastat))
    fprintf(stderr, "fastaFile::constructIndex()-- stat() of file '%s' failed: %s\n",
            _filename, strerror(errno)), exit(1);

  _header._magic[0]              = FASTA_MAGICNUMBER1;
  _header._magic[1]              = FASTA_MAGICNUMBER2;
  _header._numberOfSequences     = indexLen;
  _header._namesLength           = namesLen;
  _header._fastaFileSize         = fastastat.st_size;
  _header._fastaModificationTime = fastastat.st_mtime;
  _header._fastaCreationTime     = fastastat.st_ctime;

  //  Dump the index, if possible.

  errno = 0;
  FILE *I = fopen(indexname, "w");
  if (errno)
    return;

  fwrite(&_header,  sizeof(fastaFileHeader), 1,                          I);
  fwrite( _index,   sizeof(fastaFileIndex),  _header._numberOfSequences, I);
  fwrite( _names,   sizeof(char),            _header._namesLength,       I);

  fclose(I);
}
