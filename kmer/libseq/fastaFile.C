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

  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    strcpy(_filename, "stdin");

    _rb    = new readBuffer("-");
    _isStreamInput = true;

    //  Unknown number of sequeuces.
    _numberOfSequences = ~u32bitZERO;
    return;
  }

  strcpy(_filename, filename);

  fprintf(stderr, "fastaFile::fastaFile()-- building index for '%s'\n", (filename) ? filename : "NULLPOINTER");
  constructIndex();

  _rb    = new readBuffer(_filename);

  _numberOfSequences = _index._numberOfSequences;
}



fastaFile::fastaFile() {
  clear();
}



fastaFile::~fastaFile() {
  delete    _rb;
  delete [] _entry;
  delete [] _names;
}





seqFile *
fastaFile::openFile(const char *filename) {
  struct stat    st;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::openFile()-- '%s'\n", (filename) ? filename : "NULLPOINTER");
#endif

  //  We assume that if it's stdin, it's a fasta.
  //
  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    fprintf(stderr, "seqFile()-- stdin supplied; assuming FastA format.\n");
    return(new fastaFile("-"));
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
  fastaFile   *f = 0L;
  readBuffer  *r = new readBuffer(filename);
  char         x = r->read();

  while ((!r->eof()) && whitespaceSymbol[x])
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

  for (u32bit iid=0; iid < _index._numberOfSequences; iid++) {
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

  return((iid < _numberOfSequences) ? _entry[iid]._seqLen : 0);
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

  if (sMax == 0) {
    sMax = 2048;
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 2048;
    h    = new char [hMax];
  }

  if ((_entry) && (sMax < _entry[iid]._seqLen)) {
    sMax = _entry[iid]._seqLen;
    delete [] s;
    s = new char [sMax];
  }

  hLen = 0;
  sLen = 0;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequence(full)-- seek to iid="u32bitFMT" at pos="u32bitFMT"\n",
          iid, _entry[iid]._position);
#endif
  _rb->seek(_entry[iid]._position);

  char x = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((!_rb->eof()) && whitespaceSymbol[x])
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
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  Copy the defline, until the first newline.
  while ((!_rb->eof()) && (x != '\r') && (x != '\n')) {
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
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  Copy the sequence, until EOF or the next '>'.
  while ((!_rb->eof()) && (_rb->peek() != '>')) {
    if (!whitespaceSymbol[x]) {
      s[sLen++] = x;
      if (sLen >= sMax) {
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



bool
fastaFile::getSequence(u32bit iid,
                       u32bit bgn, u32bit end, char *s) {
  u32bit pos = 0;

#ifdef DEBUG
  fprintf(stderr, "fastaFile::getSequence(part)-- "u32bitFMT"\n", iid);
#endif

  //  It is impossible to be efficient here; see the big comment in
  //  the other getSequence() above.
  //
  //  We can't even guess where to start scanning the sequence; we
  //  just don't have any information about how much whitespace is in
  //  the sequence.

  _rb->seek(_entry[iid]._position);

  char x = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof())
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaFile::getSequence(part)-- ERROR2: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

  //  Skip the defline.
  while ((!_rb->eof()) && (x != '\r') && (x != '\n'))
    x = _rb->read();

  //  Skip whitespace between the defline and the sequence.
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  Skip sequence up until bgn.
  while ((!_rb->eof()) && (pos < bgn)) {
    if (!whitespaceSymbol[x])
      pos++;

    x = _rb->read();
  }

  //  Copy sequence
  while ((!_rb->eof()) && (pos < end)) {
    if (!whitespaceSymbol[x])
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
  memset(&_index, 0, sizeof(fastaFileIndex));
  _entry = 0L;
  _names = 0L;
  _nextID = 0;
  _isStreamInput = false;
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

  fread(&_index, sizeof(fastaFileIndex), 1, I);

  if ((_index._magic[0] != FASTA_MAGICNUMBER1) &&
      (_index._magic[1] != FASTA_MAGICNUMBER2)) {
    fprintf(stderr, "fastaFile::constructIndex()-- magic mismatch.\n");
    fclose(I);
    return;
  }

  if ((_index._fastaFileSize         != (u64bit)fastastat.st_size) ||
      (_index._fastaModificationTime != (u64bit)fastastat.st_mtime) ||
      (_index._fastaCreationTime     != (u64bit)fastastat.st_ctime)) {
    fprintf(stderr, "fastaFile::constructIndex()-- stat mismatch.\n");
    fclose(I);
    return;
  }

  _entry = new fastaFileEntry [_index._numberOfSequences];
  _names = new char           [_index._namesLength];

  fread(_entry, sizeof(fastaFileEntry), _index._numberOfSequences, I);
  fread(_names, sizeof(char),           _index._namesLength,       I);

  fprintf(stderr, "fastaFile::constructIndex()-- '%s' LOADED\n", _filename);

  fclose(I);
  return;
}


void
fastaFile::constructIndex(void) {

  if (_entry)
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

  if (_entry)
    return;

  fprintf(stderr, "fastaFile::constructIndex()-- '%s' BUILDING\n", _filename);

  //  Allocate some space for the index structures.

  u32bit  entryMax = 64 * 1024 * 1024 / sizeof(fastaFileEntry);
  u32bit  entryLen = 0;

  _entry = new fastaFileEntry [entryMax];

  u32bit  namesMax = 32 * 1024 * 1024;
  u32bit  namesLen = 0;

  _names = new char [namesMax];

  //  Some local storage

  u64bit       seqStart;
  u32bit       seqLen;

  readBuffer   ib(_filename);
  char         x = ib.read();

#ifdef DEBUGINDEX
  fprintf(stderr, "readBuffer '%s' eof=%d x=%c %d\n", _filename, ib.eof(), x, x);
#endif

  //  Build it.

  //  Skip whitespace at the start of the sequence.
  while ((!ib.eof()) && whitespaceSymbol[x]) {
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

    //  Read that first letter
    x = ib.read();

    //  Copy the name to the names
    while ((!ib.eof()) && !whitespaceSymbol[x]) {
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

    _names[namesLen++] = 0;

    //  Skip the rest of the defline
    while ((!ib.eof()) && (x != '\r') && (x != '\n')) {
#ifdef DEBUGINDEX
      fprintf(stderr, "skip let %c\n", x);
#endif
      x = ib.read();
    }

    //  Skip whitespace between the defline and the sequence.
    while ((!ib.eof()) && whitespaceSymbol[x]) {
#ifdef DEBUGINDEX
      fprintf(stderr, "skip num %d\n", x);
#endif
      x = ib.read();
    }

#ifdef DEBUGINDEX
    fprintf(stderr, "x=%c peek=%c\n", x, ib.peek());
#endif

    //  Count sequence length
    while ((!ib.eof()) && (ib.peek() != '>')) {
#ifdef DEBUGINDEX
      fprintf(stderr, "seqlen %s %c\n", (!whitespaceSymbol[x]) ? "save" : "skip", x);
#endif
      if (!whitespaceSymbol[x])
        seqLen++;
      x = ib.read();
    }

    //  Save to the index.

    if (entryLen >= entryMax) {
      entryMax *= 2;
      fastaFileEntry *et = new fastaFileEntry[entryMax];
      memcpy(et, _entry, sizeof(fastaFileEntry) * entryLen);
      delete [] _entry;
      _entry = et;
    }

    _entry[entryLen]._position = seqStart;
    _entry[entryLen]._seqLen   = seqLen;

#ifdef DEBUG
    fprintf(stderr, "INDEX iid="u32bitFMT" len="u32bitFMT" pos="u64bitFMT"\n",
            entryLen, seqLen, seqStart);
#endif

    entryLen++;

    //  Load the '>' for the next iteration.
    x = ib.read();
  }

  //  Fill out the index meta data

  struct stat  fastastat;
  errno = 0;
  if (stat(_filename, &fastastat))
    fprintf(stderr, "fastaFile::constructIndex()-- stat() of file '%s' failed: %s\n",
            _filename, strerror(errno)), exit(1);

  _index._magic[0]              = FASTA_MAGICNUMBER1;
  _index._magic[1]              = FASTA_MAGICNUMBER2;
  _index._numberOfSequences     = entryLen;
  _index._namesLength           = namesLen;
  _index._fastaFileSize         = fastastat.st_size;
  _index._fastaModificationTime = fastastat.st_mtime;
  _index._fastaCreationTime     = fastastat.st_ctime;

  //  Dump the index, if possible.

  if (_isStreamInput == false) {
    errno = 0;
    FILE *I = fopen(indexname, "w");
    if (errno)
      return;

    fwrite(&_index,   sizeof(fastaFileIndex), 1,                         I);
    fwrite( _entry,   sizeof(fastaFileEntry), _index._numberOfSequences, I);
    fwrite( _names,   sizeof(char),           _index._namesLength,       I);

    fclose(I);
  }
}
