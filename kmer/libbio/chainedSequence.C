#include "bri++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>


chainedSequence::chainedSequence() {
  _filename           = 0L;
  _file               = 0L;
  _fileToDelete       = 0L;
  _sequence           = 0L;

  _useListLen         = 0;
  _useListMax         = 32;
  _useList            = new use_s [_useListMax];

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
chainedSequence::setSource(char *filename) {
  setSource(_fileToDelete = new FastAWrapper(filename));
}


void
chainedSequence::setSource(FastAWrapper *file) {
  _file     = file;
  _file->openIndex();

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
    use_s *u = new use_s [_useListMax];
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
    _useList    = new use_s [_useListMax];
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
    if (useAll || seen[i]) {
      _useList[_useListLen].iid    = i;
      _useList[_useListLen].length = _file->sequenceLength(i);
      _useList[_useListLen].start  = startPos;

      startPos += _useList[_useListLen].length + _separatorLength;

      _useListLen++;
    }
  }

  delete [] seen;

  _lengthOfSequences = startPos;

  //  Initialize -- grab the first sequence
  //
  _file->find(_useList[0].iid);
  _sequence = _file->getSequenceOnDisk();
}
