#ifndef HITREADER_H
#define HITREADER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libbritypes.h"
#include "aHit.H"


//  XXX:  Lots of this stuff can be moved into hitReader as temporary variables
//
typedef struct {
  bool     stillMore;
  FILE    *file;
  char     b[1024];
  aHit     a;
  bool     isBINARY;
} hitFile_s;


typedef struct {
  aHit     a;
  u32bit   estid;
  double   coverage;
  double   multiplicity;
#ifdef WITH_ANSWERS
  char     mappedIdentity;
  char     mappedCoverage;
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
  u32bit iid(void)        { return(_iid);        };
  u32bit numHits(void)    { return(_listLen);    };

  void   sortByCoverage(void);

  hit_s &operator[](u32bit x) {
    if (x >= _listLen) {
      fprintf(stderr, "hitReader::operator[]()-- ERROR: asked for hit "u32bitFMT" out of "u32bitFMT".\n", x, _listLen);
      exit(1);
    }

    return(_list[x]);
  };
private:
  void   grow_List(void) {
    if (_listLen >= _listMax) {
      _listMax *= 2;
      hit_s *new_list = new hit_s [_listMax];
      memcpy(new_list, _list, _listLen * sizeof(hit_s));
      delete [] _list;
      _list = new_list;
    }
  };

  u32bit        _filesMax;
  u32bit        _filesLen;
  hitFile_s    *_files;

  u32bit        _listLen;
  u32bit        _listMax;
  hit_s        *_list;

  u32bit        _iid;
  double        _bestScore;
  double        _worstScore;
};





#endif  //  HITREADER_H
