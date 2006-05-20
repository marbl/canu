#include "hitReader.H"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>


static
int
hitCompareCoverage(const void *a, const void *b) {
  const hit_s  *A = (const hit_s *)a;
  const hit_s  *B = (const hit_s *)b;

  if (A->coverage > B->coverage)
    return(-1);
  return(A->coverage < B->coverage);
}

static
int
hitCompareGenPos(const void *a, const void *b) {
  const hit_s  *A = (const hit_s *)a;
  const hit_s  *B = (const hit_s *)b;

  if (A->a._dsIdx < B->a._dsIdx) return(-1);
  if (A->a._dsIdx > B->a._dsIdx) return(1);

  if (A->a._forward < B->a._forward) return(-1);
  if (A->a._forward > B->a._forward) return(1);

  if (A->a._dsLo  < B->a._dsLo) return(-1);
  if (A->a._dsLo  > B->a._dsLo) return(1);

  if (A->a._dsHi  < B->a._dsHi) return(-1);
  if (A->a._dsHi  > B->a._dsHi) return(1);

  return(0);
}







hitReader::hitReader(int m) {
  _filesMax   = m;
  _filesLen   = 0;
  _files      = new hitFile_s [_filesMax];

  _listLen    = 0;
  _listMax    = 1024 * 1024;
  _list       = new hit_s [_listMax];

  _iid        = u32bitZERO;
  _bestScore  = 0.0;
  _worstScore = 1.0;
}


hitReader::~hitReader() {
  for (u32bit i=0; i<_filesLen; i++)
    //fclose(_files[i].file);
    delete _files[i].buff;

  delete [] _files;
  delete [] _list;
}

void
hitReader::addInputFile(char *filename) {
  errno = 0;

  _files[_filesLen].stillMore = true;

  if (strcmp(filename, "-") == 0) {
    fprintf(stderr, "hitReader::addInputFile()-- stdin not supported!\n"), exit(1);
  } else {
    _files[_filesLen].buff = new readBuffer(filename, 16 * 1048576);
  }

  //  Binary or ASCII input?
  //
  char x = _files[_filesLen].buff->get();

  _files[_filesLen].isBINARY = (x != '-');

  //  Load the first hit
  loadHit(_files+_filesLen);

  _filesLen++;
}


void
hitReader::loadHit(hitFile_s *HF) {

  if (HF->isBINARY) {
    ahit_readBinary(&HF->a, HF->buff);
  } else {
    fprintf(stderr, "ERROR:  hitReader::loadHit() ascii not supported right now.\n");
    exit(1);
    //fgets(HF->b, 1024, HF->file);
    //ahit_parseString(&HF->a, HF->b);
  }

  if (HF->buff->eof())
    HF->stillMore = false;
};


bool
hitReader::loadHits(void) {

  _listLen    = 0;
  _iid        = u32bitZERO;
  _bestScore  = 0.0;
  _worstScore = 1.0;

  //  See if there are more hits to process.
  //
  bool  keepGoing = false;
  for (u32bit i=0; i<_filesLen; i++)
    keepGoing |= _files[i].stillMore;

  if (keepGoing == false)
    return(false);

  //  Find the lowest ESTid
  //
  _iid = 1 << 30;
  for (u32bit i=0; i<_filesLen; i++)
    if ((_files[i].stillMore) && (_iid > _files[i].a._qsIdx))
      _iid = _files[i].a._qsIdx;


  //  For each file, load the next hit if it's the est
  //  we're looking at
  //
  for (u32bit i=0; i<_filesLen; i++) {
    while ((_files[i].stillMore) && (_files[i].a._qsIdx == _iid)) {
      if (_listLen >= _listMax) {
        _listMax *= 2;
        hit_s *new_list = new hit_s [_listMax];
        memcpy(new_list, _list, _listLen * sizeof(hit_s));
        delete [] _list;
        _list = new_list;
      }

      memcpy(&_list[_listLen].a, &_files[i].a, sizeof(aHit));

      _list[_listLen].coverage     = (double)_files[i].a._covered / (double)_files[i].a._numMers;
      _list[_listLen].multiplicity = (double)_files[i].a._matched / (double)_files[i].a._covered;

      //  aHit->_covered is in bases, but aHit->_numMers is the
      //  number of mers.  Possible for coverage to be > 1.0.
      //
      if (_list[_listLen].coverage > 1.0)
        _list[_listLen].coverage = 1.0;

      if (_list[_listLen].coverage > _bestScore)
        _bestScore = _list[_listLen].coverage;

      if (_list[_listLen].coverage < _worstScore)
        _worstScore = _list[_listLen].coverage;

#ifdef WITH_ANSWERS
      //  Look for the answer string.  If not found, set to zero.
      //
      _list[_listLen].mappedIdentity = 0;
      _list[_listLen].mappedCoverage = 0;

      for (int p=0; _files[i].b[p]; p++) {
        if ((_files[i].b[p] == 'Y') || (_files[i].b[p] == 'N')) {
          char *c = _files[i].b+p+1;
          _list[_listLen].mappedIdentity = (u32bit)strtoul(c, &c, 10);
          _list[_listLen].mappedCoverage = (u32bit)strtoul(c, &c, 10);
        }
      }
#endif

      _listLen++;

      loadHit(_files+i);
    }
  }

  mergeOverlappingHits();

  return(true);
}



void
hitReader::sortByCoverage(void) {
  qsort(_list, _listLen, sizeof(hit_s), hitCompareCoverage);
};





//  scan the list of hits (for a single EST, remember) and merge
//  any that are overlapping
//
void
hitReader::mergeOverlappingHits(void) {

  //  Sort by the genomic position
  //
  qsort(_list, _listLen, sizeof(hit_s), hitCompareGenPos);

  //  Scan through the list, merging.
  //
  u32bit   cur = 0;  //  Currently active entry
  u32bit   exa = 1;  //  Entry we examine for merging
  while (exa < _listLen) {

    //  Do they overlap?
    if ((_list[cur].a._dsIdx     == _list[exa].a._dsIdx) &&
        (_list[cur].a._forward   == _list[exa].a._forward) &&
        (_list[cur].a._dsHi      >= _list[exa].a._dsLo)) {

      //  Yup, merge.  Extend the current hit if it is smaller.

      if ((_list[cur].a._dsLo == _list[exa].a._dsLo) &&
          (_list[cur].a._dsHi == _list[exa].a._dsHi)) {
        //  Nop, they're the same.
      } else if (_list[cur].a._dsHi >= _list[exa].a._dsHi) {
        //  Nop, exa is contained in cur.
      } else {

        //  exa extends cur!

        //  If cur is contained in exa, just get rid of cur.
        //  Otherwise, we need to fudge up new scores -- but we
        //  instead just mark them as merged, and don't filter them.
        //
        if (_list[cur].a._dsLo == _list[exa].a._dsLo) {
          memcpy(_list+cur, _list+exa, sizeof(hit_s));
        } else {
#ifdef DEBUG_HITREADER
          fprintf(stderr,
                  "MERGE:   ("u32bitFMT","u32bitFMT") -e "u32bitFMT" "
                  u32bitFMT":"u32bitFMT"-"u32bitFMT"%c("u32bitFMT"-"u32bitFMT"-"u32bitFMT") "
                  u32bitFMT":"u32bitFMT"-"u32bitFMT"%c("u32bitFMT"-"u32bitFMT"-"u32bitFMT")\n",
                  cur, exa,
                  _list[cur].a._qsIdx,
                  _list[cur].a._dsIdx,
                  _list[cur].a._dsLo,    _list[cur].a._dsHi,    _list[cur].a._forward ? 'f' : 'r',
                  _list[cur].a._covered, _list[cur].a._matched, _list[cur].a._numMers,
                  _list[exa].a._dsIdx,
                  _list[exa].a._dsLo,    _list[exa].a._dsHi,    _list[exa].a._forward ? 'f' : 'r',
                  _list[exa].a._covered, _list[exa].a._matched, _list[exa].a._numMers);
#endif
          _list[cur].a._merged  = true;
          _list[cur].a._covered = 0;
          _list[cur].a._matched = 0;
          _list[cur].a._dsHi    = _list[exa].a._dsHi;
        }
      }

      //  By now, we've updated cur to include all that exa did.  exa is junk.

    } else {
      //  Nope, copy exa to the next spot (unless they're the same)
      //  and move there.
      //
      cur++;
      if (cur != exa)
        memcpy(_list+cur, _list+exa, sizeof(hit_s));
    }

    //  Move to the next examination!
    exa++;
  }

  _listLen = cur + 1;
}
