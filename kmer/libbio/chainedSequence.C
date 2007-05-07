#include "bio++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>


static
char     magic[16] = { 'c','h','a','i','n','e','d','S','e','q','.','v','1','.','.','.' };
static
char     faild[16] = { 'c','h','a','i','n','e','d','S','e','q','F','A','I','L','E','D' };


chainedSequence::chainedSequence() {
  _filename           = 0L;
  _file               = 0L;
  _fileToDelete       = 0L;
  _sequence           = 0L;

  _useListLen         = 0;
  _useListMax         = 32;
  _useList            = new use_s [_useListMax + 1];

  _currentSeq         = 0;

  _positionInSequence = 0;
  _positionInStream   = 0;
  _lengthOfSequences  = 0;

  _eof                = false;

  _separator          = '.';
  _separatorLength    = 9;  //  Default of 10 .'s
  _separatorPosition  = 0;
  _separatorDone      = false;
}


chainedSequence::~chainedSequence() {
  delete [] _filename;
  delete    _fileToDelete;
  delete [] _useList;
  delete    _sequence;
}


void
chainedSequence::setSeparator(char sep) {
  _separator = sep;
}


void
chainedSequence::setSeparatorLength(u32bit len) {
  if (len == 0) {
    fprintf(stderr, "chainedSequence()-- ERROR: the separator length must be at least one.\n");
    exit(1);
  }
  _separatorLength = len - 1;
}


void
chainedSequence::setSource(char const *filename) {
  setSource(_fileToDelete = openSeqFile(filename));
}


void
chainedSequence::setSource(seqFile *file) {
  _file = file;
  _file->openIndex();

  if (_file->isSqueezed() == false) {
    fprintf(stderr, "ERROR:  chainedSequence::setSource()-- source file '%s' isn't squeezed,\n", _file->getSourceName());
    fprintf(stderr, "ERROR:      and this causes errors in indexing.  Bri should really fix this.\n");
    fprintf(stderr, "ERROR:      until then, squeeze it with:\n");
    fprintf(stderr, "ERROR:  leaff -f %s -W > %s.squeezed\n", _file->getSourceName(), _file->getSourceName());
    exit(1);
  }

  _filename = new char [strlen(_file->getSourceName()) + 1];
  strcpy(_filename, _file->getSourceName());
}


void
chainedSequence::parse(char *line) {
  u32bit  v=0, u=0;
  char   *rest = line;

  //  line can be either a list of numbers, or a file.  See if "line" opens
  //  as a file.  If it does, read in the numbers.
  //
  FILE *F = fopen(line, "r");
  if (F) {
    while (!feof(F)) {
      if (fscanf(F, " "u32bitFMT" ", &v))
        add(v);
    }
    fclose(F);
  } else {
    while (*line) {

      //  We are at a number.  Get it.
      //
      v = (u32bit)strtoul(line, &rest, 10);
      line = rest;
      
      //  If we get ',' or 0, add the single number to the list.
      //  If we get '-', decode the rest of the range.
      //
      switch (*line) {
        case 0:
        case ',':
          add(v);
          break;
        case '-':
          line++;
          u = (u32bit)strtoul(line, &rest, 10);
          line = rest;

          if (v > u)
            for (; u <= v; u++)
              add(u);
          else
            for (; v <= u; v++)
              add(v);
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- no separator found -- format should be, e.g., \"1,2,4-9,10\".\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }

      //  We should be at a ',' (or end of line).  Anything else is
      //  an error.  If ',', move to the next number.
      //
      switch (*line) {
        case 0:
          break;
        case ',':
          line++;
          break;
        default:
          fprintf(stderr, "Invalid -use specification -- no separator found -- format should be, e.g., \"1,2,4-9,10\".\n");
          fprintf(stderr, "Trouble starts at '%s'\n", line);
          exit(1);
          break;
      }
    }
  }
}


void
chainedSequence::add(u32bit v) {

  if (_useListLen >= _useListMax) {
    _useListMax <<= 1;
    use_s *u = new use_s [_useListMax + 1];
    memcpy(u, _useList, sizeof(use_s) * _useListLen);
    delete [] _useList;
    _useList = u;
  }

  _useList[_useListLen].iid    = v;
  _useList[_useListLen].length = 0;
  _useList[_useListLen].start  = 0;
  _useListLen++;
}


u64bit
chainedSequence::lengthOf(u32bit s) {
  if (s >= _useListLen)
    return(~u32bitZERO);
  return(_useList[s].length);
}


u64bit
chainedSequence::startOf(u32bit s) {
  if (s >= _useListLen)
    return(~u32bitZERO);
  return(_useList[s].start);
}


u64bit
chainedSequence::IIDOf(u32bit s) {
  if (s >= _useListLen)
    return(~u32bitZERO);
  return(_useList[s].iid);
}


u32bit
chainedSequence::sequenceNumberOfPosition(u64bit p) {
  u32bit   sret = ~u32bitZERO;

  //  binary search on our list of start positions, to find the
  //  sequence that p is in.

  if (_lengthOfSequences < p) {
    fprintf(stderr, "chainedSequence::sequenceNumberOfPosition()-- WARNING! Position p="u64bitFMT" too big; only "u64bitFMT" positions.\n",
            p, _lengthOfSequences);
    return(~u32bitZERO);
  }

  if (_useListLen < 16) {
    for (u32bit ss=0; ss<_useListLen; ss++) {
      if ((_useList[ss].start <= p) && (p < _useList[ss+1].start)) {
        sret = ss;
        break;
      }
    }
  } else {
    u32bit  lo = 0;
    u32bit  hi = _useListLen;
    u32bit  md = 0;

    while (lo <= hi) {
      md = (lo + hi) / 2;

      if        (p < _useList[md].start) {
        //  This block starts after the one we're looking for.  
        hi = md;

      } else if ((_useList[md].start <= p) && (p < _useList[md+1].start)) {
        //  Got it!
        lo           = md + 1;
        hi           = md;
        sret         = md;

      } else {
        //  By default, then, the block is too low.
        lo = md;
      }
    }
  }

  //fprintf(stderr, "chainedSequence::sequenceNumberOfPosition()-- p="u64bitFMT" --> s="u32bitFMT"\n", p, sret);

  return(sret);
}


bool
chainedSequence::rewind(void) {
  _currentSeq         = 0;
  _positionInSequence = 0;
  _positionInStream   = 0;
  _eof                = false;
  _separatorDone      = false;
  _separatorPosition  = 0;

  delete _sequence;

  _file->find(_useList[0].iid);
  _sequence = _file->getSequenceOnDisk();

  return(true);
}


void
chainedSequence::finish(void) {
  bool    useAll   = false;
  bool   *seen     = 0L;
  u64bit  seenLen  = 0;
  u64bit  startPos = 0;

  //  If silly user didn't add any sequences, use all the sequences.
  //  Otherwise, fix the existing one (remove duplicates, sort).

  if (_useListLen == 0) {
    useAll  = true;
    seen    = 0L;
    seenLen = _file->getNumberOfSequences();

    delete [] _useList;

    _useListMax = seenLen;
    _useList    = new use_s [_useListMax + 1];
  } else {
    useAll  = false;
    seen    = new bool [_file->getNumberOfSequences()];
    seenLen = 0;

    for (u32bit i=_file->getNumberOfSequences(); i--; )
      seen[i] = false;

    for (u32bit i=0; i<_useListLen; i++) {
      if (_useList[i].iid >= _file->getNumberOfSequences()) {
        fprintf(stderr, "WARNING: chainedSequence requested sequence iid "u32bitFMT" which isn't in '%s'\n",
                i, _filename);
      } else {
        if (seen[_useList[i].iid] == false) {
          seen[_useList[i].iid] = true;
          seenLen++;
        }
      }
    }
  }

  _useListLen = 0;

  for (u32bit i=0; i<_file->getNumberOfSequences(); i++) {
    if ((useAll || seen[i]) &&
        (_file->sequenceLength(i) > 0)) {
      _useList[_useListLen].iid    = i;
      _useList[_useListLen].length = _file->sequenceLength(i);
      _useList[_useListLen].start  = startPos;

#if 0
      fprintf(stderr, "chainedSequence::finish()-- added iid="u32bitFMTW(6)" len="u32bitFMTW(9)" start="u32bitFMTW(9)"\n",
              _useList[_useListLen].iid,
              _useList[_useListLen].length,
              _useList[_useListLen].start);
#endif

      startPos += _useList[_useListLen].length + _separatorLength + 1;

      _useListLen++;
    }
  }

  delete [] seen;

  _lengthOfSequences = startPos;

  _useList[_useListLen].iid      = ~u32bitZERO;
  _useList[_useListLen].length   = 0;
  _useList[_useListLen].start    = _lengthOfSequences;

  //  Initialize -- grab the first sequence
  //
  _file->find(_useList[0].iid);
  _sequence = _file->getSequenceOnDisk();
}


void
chainedSequence::saveState(char const *filename) {

  errno = 0;
  int F = open(filename, O_RDWR | O_CREAT | O_LARGEFILE,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno)
    fprintf(stderr, "Can't open '%s' for writing chainedSequence.\n%s\n", filename, strerror(errno)), exit(1);

  write(F, magic, sizeof(char) * 16);
  if (errno)
    fprintf(stderr, "chainedSequence::saveState()-- Write failure on magic.\n%s\n", strerror(errno)), exit(1);

  u32bit len = strlen(_filename) + 1;

  write(F, &len,                 sizeof(u32bit));
  write(F,  _filename,           sizeof(char) * len);
  write(F, &_useListLen,         sizeof(u32bit));
  write(F,  _useList,            sizeof(use_s) * _useListLen);
  write(F, &_positionInSequence, sizeof(u64bit));
  write(F, &_positionInStream,   sizeof(u64bit));
  write(F, &_lengthOfSequences,  sizeof(u64bit));
  write(F, &_eof,                sizeof(bool));
  write(F, &_separatorDone,      sizeof(bool));
  write(F, &_separatorLength,    sizeof(u32bit));
  write(F, &_separatorPosition,  sizeof(u32bit));
  write(F, &_separator,          sizeof(char));

  write(F, magic, sizeof(char) * 16);
  if (errno)
    fprintf(stderr, "chainedSequence::saveState()-- Write failure on final magic.\n%s\n", strerror(errno)), exit(1);

  close(F);
}


bool
chainedSequence::loadState(char const *filename, bool beNoisy, bool loadData) {
  char  cigam[16];

  errno = 0;
  int F = open(filename, O_RDONLY | O_LARGEFILE);
  if (errno)
    fprintf(stderr, "Can't open '%s' for reading chainedSequence.\n%s\n", filename, strerror(errno)), exit(1);

  read(F, cigam, sizeof(char) * 16);
  if (errno)
    fprintf(stderr, "Can't read magic from chainedSequence '%s'.\n%s\n", filename, strerror(errno)), exit(1);
  if (strncmp(cigam, magic, 16) != 0)
    fprintf(stderr, "Invalid magic in chainedSequence '%s'.\n%s\n", filename, strerror(errno)), exit(1);

  u32bit len = 0;
  read(F, &len,                 sizeof(u32bit));
  _filename = new char [len];
  read(F,  _filename,           sizeof(char) * len);
  read(F, &_useListLen,         sizeof(u32bit));
  _useList = new use_s [_useListLen];
  read(F,  _useList,            sizeof(use_s) * _useListLen);
  read(F, &_positionInSequence, sizeof(u64bit));
  read(F, &_positionInStream,   sizeof(u64bit));
  read(F, &_lengthOfSequences,  sizeof(u64bit));
  read(F, &_eof,                sizeof(bool));
  read(F, &_separatorDone,      sizeof(bool));
  read(F, &_separatorLength,    sizeof(u32bit));
  read(F, &_separatorPosition,  sizeof(u32bit));
  read(F, &_separator,          sizeof(char));
  if (errno)
    fprintf(stderr, "Can't read header from chainedSequence '%s'.\n%s\n", filename, strerror(errno)), exit(1);

  read(F, cigam, sizeof(char) * 16);
  if (errno)
    fprintf(stderr, "Can't read final magic from chainedSequence '%s'.\n%s\n", filename, strerror(errno)), exit(1);
  if (strncmp(cigam, magic, 16) != 0)
    fprintf(stderr, "Invalid final magic in chainedSequence '%s'.\n%s\n", filename, strerror(errno)), exit(1);

  setSource(_filename);
  finish();

  return(true);
}
