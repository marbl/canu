#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//  A very simple seatac filter.  It reports the single longest match for each pair.
//
//  Also shows how to use a C++ object as a filter.  C is pretty much the same thing.

#include "bio.h"
#include "util++.H"

extern "C" {
  void    *construct(char *options);
  void     destruct(void *handle);
  void     addHit(void *handle,
                  char    orientation,
                  uint32  id1,
                  uint32  pos1,
                  uint32  len1,
                  uint32  id2,
                  uint32  pos2,
                  uint32  len2,
                  uint32  filled);
  void     filter(void *handle);
  uint64   output(void *handle, FILE *file, uint64 matchid);

  void    *constructStats(char *options);
  void     destructStats(void *handle);
  void     addStats(void *handle, void *filterhandle);
  void     showStats(void *handle, FILE *file);
}



class filterLongest {
public:
  filterLongest(char *n1, char *n2) {
    fprintf(stderr, "Creating a filterLongest\n");
    strncpy(name1, n1, 31);
    strncpy(name2, n2, 31);
  };

  ~filterLongest() {
    fprintf(stderr, "Destroyed a filterLongest\n");
  };

  void addHit(char    orientation,
              uint32  id1,
              uint32  pos1,
              uint32  len1,
              uint32  id2,
              uint32  pos2,
              uint32  len2,
              uint32  filled) {

    if (maxfilled < filled) {
      fprintf(stderr, "filterNOP-- addHit\n");

      maxfilled = filled;
#if 0
      sprintf(outstring,
              "-%c -e "uint32FMT" "uint32FMT" "uint32FMT" -D "uint32FMT" "uint32FMT" "uint32FMT" -F "uint32FMT"\n",
              orientation, id1, pos1, len1, id2, pos2, len2, filled);
#endif

      sprintf(outstring,
              "M x . . %s:"uint32FMT" "uint32FMT" "uint32FMT" 1 %s:"uint32FMT" "uint32FMT" "uint32FMT" %s "uint32FMT"\n",
              name1, id1, pos1, len1, name2, id2, pos2, len2, (orientation == 'f') ? "1" : "-1", filled);
    }
  };

  void filter(void) {
    fprintf(stderr, "filterNOP-- filter\n");
  };

  uint64 output(FILE *file, uint64 matchid) {
    fprintf(stderr, "filterNOP-- output (ignoring matchid)\n");
    fprintf(file, "%s", outstring);
    return(matchid);
  };
private:
  char     outstring[512];
  char     name1[32], name2[32];
  uint32   maxfilled;
};



class statLongest {
public:
  statLongest() {
    num = 0;
  }
  ~statLongest() {
  }

  void add(filterLongest *F) {
    num++;
  }

  void show(FILE *file) {
    fprintf(file, "/statObjNum=%d\n", num);
  }
  
private:
  int    num;
};






void*
construct(char *opts) {
  char *seq1 = "UNK";
  char *seq2 = "UNK";

  //  Parse the options to find the parameters
  //
  splitToWords  W(opts);

  uint32 arg = 0;
  while (arg < W.numWords()) {
    if        (strcmp(W.getWord(arg), "-1") == 0) {
      seq1 = W.getWord(++arg);
    } else if (strcmp(W.getWord(arg), "-2") == 0) {
      seq2 = W.getWord(++arg);
    }

    arg++;
  }

  return(new filterLongest(seq1, seq2));
}

void
destruct(void *handle) {
  delete (filterLongest *)handle;
}

void
addHit(void   *handle,
       char    orientation,
       uint32  id1,
       uint32  pos1,
       uint32  len1,
       uint32  id2,
       uint32  pos2,
       uint32  len2,
       uint32  filled) {
  ((filterLongest *)handle)->addHit(orientation, id1, pos1, len1, id2, pos2, len2, filled);
}

void
filter(void *handle) {
  ((filterLongest *)handle)->filter();
}

uint64
output(void *handle, FILE *file, uint64 matchid) {
  return(((filterLongest *)handle)->output(file, matchid));
}




void*
constructStats(char *options) {
  return(new statLongest);
}

void
destructStats(void *handle) {
  delete (statLongest *)handle;
}

void
addStats(void *handle, void *filterhandle) {
  ((statLongest *)handle)->add((filterLongest *)filterhandle);
}

void
showStats(void *handle, FILE *file) {
  ((statLongest *)handle)->show(file);
}
