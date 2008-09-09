#include "fastaFile.H"

#include "alphabet.h"

//  Says 'kmerFastaFileIdx'
#define FASTA_MAGICNUMBER1  0x7473614672656d6bULL
#define FASTA_MAGICNUMBER2  0x786449656c694661ULL


fastaFile::fastaFile(const char *filename) {
  clear();

  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    strcpy(_filename, "stdin");
    strcpy(_typename, "fasta");

    _rb    = new readBuffer(fileno(stdin), _filename);
    _isStreamInput = true;

    return;
  }

  strcpy(_filename, filename);
  strcpy(_typename, "fasta");

  _rb    = new readBuffer(_filename);
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

  //  We assume that if it's stdin, it's a fasta.
  //
  if ((filename == 0L) || (strcmp(filename, "-") == 0)) {
    fprintf(stderr, "seqFile()-- stdin supplied; assuming FastA format.\n");
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
  char       x = f->_rb->read();

  while ((!f->_rb->eof()) && whitespaceSymbol[x])
    x = f->_rb->read();

  //  If we get a fasta record separator assume it's a fasta file.  If
  //  it's eof, the file is empty, and we might as well return this
  //  fasta file and let the client deal with the lack of sequence.
  //
  if ((x == '>') || (f->_rb->eof() == true)) {
    f->_rb->seek(0);
    return(f);
  }

  delete f;

  return(0L);
}



u32bit
fastaFile::find(const char *sequencename) {
  constructIndex();
  fprintf(stderr, "fastaFile::find()--  Not implemented.\n");
  exit(1);
}



u32bit
fastaFile::getSequenceLength(u32bit iid) {
  constructIndex();
  return((iid < _numberOfSequences) ? _entry[iid]._seqLen : 0);
}



bool
fastaFile::getSequence(u32bit iid,
                       char *&h, u32bit &hLen, u32bit &hMax,
                       char *&s, u32bit &sLen, u32bit &sMax) {
  if (iid != _nextID) {
    constructIndex();
    _rb->seek(_entry[iid]._position);
  }

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
  
  char x = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof())
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaFile::getSequence()-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
            _filename, x), exit(1);

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

  _nextID++;

  return(true);
}



bool
fastaFile::getSequence(u32bit iid,
                       u32bit bgn, u32bit end, char *s) {
  constructIndex();

  _rb->seek(_entry[iid]._position);

  //  It is impossible to be efficient here; see the big comment in
  //  the other getSequence() above.
  //
  //  We can't even guess where to start scanning the sequence; we
  //  just don't have any information about how much whitespace is in
  //  the sequence.

  u32bit pos = 0;

  char x = _rb->read();

  //  Skip whitespace at the start of the sequence.
  while ((!_rb->eof()) && whitespaceSymbol[x])
    x = _rb->read();

  //  We should be at a '>' character now.  Fail if not.
  if (_rb->eof())
    return(false);
  if (x != '>')
    fprintf(stderr, "fastaFile::getSequence()-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
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

  //  Fail if we didn't copy enough stuff.
  return((pos == end) ? true : false);
}




void
fastaFile::clear(void) {
  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);
  _numberOfSequences = 0;

  _rb                = 0L;
  memset(&_index, 0, sizeof(fastaFileIndex));
  _entry = 0L;
  _names = 0L;
  _nextID = 0;
  _isStreamInput = false;
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

  if (fileExists(indexname)) {

  }

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

  //  Build it.

  //  Skip whitespace at the start of the sequence.
  while ((!ib.eof()) && whitespaceSymbol[x])
    x = ib.read();

  while (ib.eof() == false) {

    //  We should be at a '>' character now.  Fail if not.
    if (x != '>')
      fprintf(stderr, "fastaFile::constructIndex()-- ERROR: In %s, expected '>' at beginning of defline, got '%c' instead.\n",
              _filename, x), exit(1);

    //  Save info.
    seqStart = ib.tell();
    seqLen   = 0;

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
      x = ib.read();
    }

    _names[namesLen++] = 0;

    //  Skip the rest of the defline
    while ((!ib.eof()) && (x != '\r') && (x != '\n'))
      x = ib.read();

    //  Skip whitespace between the defline and the sequence.
    while ((!ib.eof()) && whitespaceSymbol[x])
      x = ib.read();

    //  Count sequence length
    while ((!_rb->eof()) && (_rb->peek() != '>')) {
      if (!whitespaceSymbol[x])
        seqLen++;
      x = _rb->read();
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

    entryLen++;
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

    fwrite(&_index,   sizeof(fastaFileIndex), 1,        I);
    fwrite( _entry,   sizeof(fastaFileEntry), entryLen, I);
    fwrite( _names,   sizeof(char),           namesLen, I);

    fclose(I);
  }
}
