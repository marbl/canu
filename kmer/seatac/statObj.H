#ifndef STATOBJ_H
#define STATOBJ_H

#include <stdio.h>
#include <stdlib.h>
#include "sharedObj.H"


class statObj {
public:
  statObj(sharedObj *so, char *soOpts);
  ~statObj();

  void   add(filterObj *FO) {
    if (soadd)
      (*soadd)(handle, FO->handle);
  };

  void   show(FILE *file) {
    if (soshow)
      (*soshow)(handle, file);
  };

private:
  char        *soOpts;
  sharedObj   *so;

  void        *handle;

  void *     (*soconstruct)(char *);
  void       (*sodestruct)(void *);
  void       (*soadd)(void *, void *);
  void       (*soshow)(void *, FILE *);
};


inline
statObj::statObj(sharedObj *so_, char *op_) {
  soOpts       = 0L;
  so           = so_;

  handle       = 0L;

  soconstruct  = 0L;
  soadd        = 0L;
  soshow       = 0L;
  sodestruct   = 0L;

  if (op_) {
    soOpts = new char [strlen(op_) + 1];
    strcpy(soOpts, op_);
  }

  if (so) {
    soconstruct = (void* (*)(char*))so->get("constructStats");
    sodestruct  = (void  (*)(void*))so->get("destructStats");
    soadd       = (void  (*)(void*, void*))so->get("addStats");
    soshow      = (void  (*)(void*, FILE*))so->get("showStats");

    handle = (*soconstruct)(soOpts);
  }
}


inline
statObj::~statObj() {
  if (sodestruct)
    (*sodestruct)(handle);
  delete [] soOpts;
}


#endif  //  STATOBJ_H
