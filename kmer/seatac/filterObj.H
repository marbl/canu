#ifndef FILTEROBJ_H
#define FILTEROBJ_H

#include <stdio.h>
#include <stdlib.h>
#include "sharedObj.H"

//
//  Object that will filter and output hits.  If no sharedObj is
//  supplied, the default behavior is output all hits.
//

class filterObj {
public:
  filterObj(sharedObj *so, char *soOpts);
  ~filterObj();

  void   addHit(char   direction,
                u32bit id1,
                u32bit pos1,
                u32bit len1,
                u32bit id2,
                u32bit pos2,
                u32bit len2,
                u32bit filled);
  void   filter(void);
  void   output(FILE *);

private:
  char        *soOpts;
  sharedObj   *so;

  void        *handle;

  void *     (*soconstruct)(char *);
  void       (*sodestruct)(void *);
  void       (*soaddHit)(void *, char, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit);
  void       (*sofilter)(void *);
  void       (*sooutput)(void*,FILE *);

  u32bit       theOutputPos;
  u32bit       theOutputMax;
  char        *theOutput;

  friend class statObj;
};



inline
filterObj::filterObj(sharedObj *so_, char *op_) {
  soOpts       = 0L;
  so           = so_;

  handle       = 0L;

  soconstruct  = 0L;
  sodestruct   = 0L;
  soaddHit     = 0L;
  sofilter     = 0L;
  sooutput     = 0L;

  theOutputPos = 0;
  theOutputMax = 0;
  theOutput    = 0L;

  if (op_) {
    soOpts = new char (strlen(op_) + 1);
    strcpy(soOpts, op_);
  }

  if (so) {
    soconstruct = (void* (*)(char *))so->get("construct");
    sodestruct  = (void  (*)(void*))so->get("destruct");
    soaddHit    = (void  (*)(void *, char, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit))so->get("addHit");
    sofilter    = (void  (*)(void*))so->get("filter");
    sooutput    = (void  (*)(void*,FILE*))so->get("output");

    handle = (*soconstruct)(soOpts);
  }

  if (!so) {
    theOutputPos = 0;
    theOutputMax = 1048576;
    theOutput    = new char [theOutputMax];
  }
}

inline
filterObj::~filterObj() {
  if (sodestruct) {
    (*sodestruct)(handle);
  }
  delete [] soOpts;
  delete [] theOutput;
}


inline
void
filterObj::addHit(char   direction,
                  u32bit id1,
                  u32bit pos1,
                  u32bit len1,
                  u32bit id2,
                  u32bit pos2,
                  u32bit len2,
                  u32bit filled) {

  if (soaddHit) {
    (*soaddHit)(handle, direction, id1, pos1, len1, id2, pos2, len2, filled);
  } else {
    if (theOutputPos + 128 >= theOutputMax) {
      theOutputMax <<= 1;
      char *o = 0L;
      try {
        o = new char [theOutputMax];
      } catch (std::bad_alloc) {
        fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
        fprintf(stderr, "hitMatrix::filter()-- tried to extend output string from %lu to %lu bytes.\n", theOutputPos, theOutputMax);
        exit(1);
      }
      memcpy(o, theOutput, theOutputPos);
      delete [] theOutput;
      theOutput = o;
    }

    sprintf(theOutput + theOutputPos,
            "-%c -e "u32bitFMT" "u32bitFMT" "u32bitFMT" -D "u32bitFMT" "u32bitFMT" "u32bitFMT" -F "u32bitFMT"\n",
            direction, id1, pos1, len1, id2, pos2, len2, filled);

    while (theOutput[theOutputPos])
      theOutputPos++;
  }
}



inline
void
filterObj::filter(void) {

  if (sofilter) {
    (*sofilter)(handle);
  }
}

inline
void
filterObj::output(FILE *F) {
  if (sooutput) {
    (*sooutput)(handle, F);
  } else {
    fwrite(theOutput, sizeof(char), theOutputPos, F);
  }  
}


#endif  //  FILTEROBJ_H
