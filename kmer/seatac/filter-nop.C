#include <stdio.h>
#include <stdlib.h>

//  A very simple seatac filter.  It reports the single longest match for each pair.
//
//  Also shows how to use a C++ object as a filter.  C is pretty much the same thing.

#include "libbritypes.h"

extern "C" {
  void    *construct(char *options);
  void     destruct(void *handle);
  void     addHit(void *handle,
                  char    orientation,
                  u32bit  id1,
                  u32bit  pos1,
                  u32bit  len1,
                  u32bit  id2,
                  u32bit  pos2,
                  u32bit  len2,
                  u32bit  filled);
  void     filter(void *handle);
  void     output(void *handle, FILE *file);

  void    *constructStats(char *options);
  void     destructStats(void *handle);
  void     addStats(void *handle, void *filterhandle);
  void     showStats(void *handle, FILE *file);
}



class filterLongest {
public:
  filterLongest() {
    fprintf(stderr, "Creating a filterLongest\n");
    id1 = pos1 = len1 = id2 = pos2 = len2 = filled = maxfilled = 0;
  };

  ~filterLongest() {
    fprintf(stderr, "Destroyed a filterLongest\n");
  };

  void addHit(char    orientation,
              u32bit  id1,
              u32bit  pos1,
              u32bit  len1,
              u32bit  id2,
              u32bit  pos2,
              u32bit  len2,
              u32bit  filled) {

    if (maxfilled < filled) {
      fprintf(stderr, "filterNOP-- addHit\n");

      maxfilled = filled;
      sprintf(outstring,
              "-%c -e "u32bitFMT" "u32bitFMT" "u32bitFMT" -D "u32bitFMT" "u32bitFMT" "u32bitFMT" -F "u32bitFMT"\n",
              orientation, id1, pos1, len1, id2, pos2, len2, filled);
    }
  };

  void filter(void) {
    fprintf(stderr, "filterNOP-- filter\n");
  };

  void output(FILE *file) {
    fprintf(stderr, "filterNOP-- output\n");
    fprintf(stderr, "OUTPUT %s", outstring);
    fprintf(file, "%s", outstring);
  };
private:
  char   outstring[512];
  u32bit id1, pos1, len1, id2, pos2, len2, filled;
  u32bit maxfilled;
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
    fprintf(file, "statObj: num = %d\n", num);
  }
  
private:
  int    num;
};






void*
construct(char *opts) {
  return(new filterLongest);
}

void
destruct(void *handle) {
  delete (filterLongest *)handle;
}

void
addHit(void   *handle,
       char    orientation,
       u32bit  id1,
       u32bit  pos1,
       u32bit  len1,
       u32bit  id2,
       u32bit  pos2,
       u32bit  len2,
       u32bit  filled) {
  ((filterLongest *)handle)->addHit(orientation, id1, pos1, len1, id2, pos2, len2, filled);
}

void
filter(void *handle) {
  ((filterLongest *)handle)->filter();
}

void
output(void *handle, FILE *file) {
  ((filterLongest *)handle)->output(file);
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
