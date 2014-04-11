#ifndef HITREADER_H
#define HITREADER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bio.h"
#include "aHit.H"

//  Define this if your hits have answers and you're curious about how
//  well the filter is performing.
//
//  This is used in filterEST.C also.
//
//#define WITH_ANSWERS


//  XXX:  Lots of this stuff can be moved into hitReader as temporary variables
//
typedef struct {
  bool         stillMore;
  readBuffer  *buff;
  char         b[1024];
  aHit         a;
  bool         isBINARY;
} hitFile_s;


typedef struct {
  aHit     a;
  double   coverage;
  double   multiplicity;
  uint32   estid;
#ifdef WITH_ANSWERS
  uint32   mappedIdentity;
  uint32   mappedCoverage;
#endif
} hit_s;



class hitReader {
public:
  hitReader(int m);
  ~hitReader();

  void   addInputFile(char *filename);

  void   loadHit(hitFile_s *HF);
  bool   loadHits(void);

  double bestScore(void)  { return(_bestScore);  };
  double worstScore(void) { return(_worstScore); };
  uint32 iid(void)        { return(_iid);        };
  uint32 numHits(void)    { return(_listLen);    };

  void   sortByCoverage(void);
  void   mergeOverlappingHits(void);

  hit_s &operator[](uint32 x) {
    if (x >= _listLen) {
      fprintf(stderr, "hitReader::operator[]()-- ERROR: asked for hit "uint32FMT" out of "uint32FMT".\n", x, _listLen);
      exit(1);
    }

    return(_list[x]);
  };
private:
  uint32        _filesMax;
  uint32        _filesLen;
  hitFile_s    *_files;

  uint32        _listLen;
  uint32        _listMax;
  hit_s        *_list;

  double        _bestScore;
  double        _worstScore;
  uint32        _iid;
};

#endif  //  HITREADER_H
