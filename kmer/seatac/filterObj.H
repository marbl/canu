#ifndef FILTEROBJ_H
#define FILTEROBJ_H

#include <stdio.h>
#include <stdlib.h>
#include "sharedObj.H"
#include "util++.H"

//
//  Object that will filter and output hits.  If no sharedObj is
//  supplied, the default behavior is output all hits.
//

//  The default filter here inserts lots of '#'s into the output string,
//  then replaces those with the real match id on output.  An alternative
//  (and probably better idea) is to build a vector of structs of things
//  to output.

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
  u64bit output(FILE *, u64bit);

private:
  char        *soOpts;
  sharedObj   *so;

  void        *handle;

  void *     (*soconstruct)(char *);
  void       (*sodestruct)(void *);
  void       (*soaddHit)(void *, char, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit);
  void       (*sofilter)(void *);
  u64bit     (*sooutput)(void *, FILE *, u64bit);

  u32bit       theOutputPos;
  u32bit       theOutputMax;
  char        *theOutput;

  char         name1[32];
  char         name2[32];

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

  strcpy(name1, "UNK");
  strcpy(name2, "UNK");

  if (op_) {
    soOpts = new char [strlen(op_) + 1];
    strcpy(soOpts, op_);

    splitToWords  W(soOpts);

    for (u32bit arg = 0; arg < W.numWords(); arg++) {
      if        (strcmp(W.getWord(arg), "-1") == 0) {
        strncpy(name1, W.getWord(++arg), 31);
      } else if (strcmp(W.getWord(arg), "-2") == 0) {
        strncpy(name2, W.getWord(++arg), 31);
      }
    }
  }

  if (so) {
    soconstruct = (void*  (*)(char *))so->get("construct");
    sodestruct  = (void   (*)(void*))so->get("destruct");
    soaddHit    = (void   (*)(void *, char, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit, u32bit))so->get("addHit");
    sofilter    = (void   (*)(void*))so->get("filter");
    sooutput    = (u64bit (*)(void*,FILE*,u64bit))so->get("output");

    if (!soconstruct) fprintf(stderr, "construct not found!\n");
    if (!sodestruct)  fprintf(stderr, "destruct not found!\n");
    if (!soaddHit)    fprintf(stderr, "addHit not found!\n");
    if (!sofilter)    fprintf(stderr, "filter not found!\n");
    if (!sooutput)    fprintf(stderr, "output not found!\n");

    handle = (*soconstruct)(soOpts);
  }

  if (!so) {
    theOutputPos = 0;
    theOutputMax = 1048576;
    theOutput    = new char [theOutputMax];
    theOutput[0] = 0;
  }
}

inline
filterObj::~filterObj() {
  if (sodestruct)
    (*sodestruct)(handle);
  delete [] soOpts;
  delete [] theOutput;
}


inline
void
filterObj::addHit(char   orientation,
                  u32bit id1,
                  u32bit pos1,
                  u32bit len1,
                  u32bit id2,
                  u32bit pos2,
                  u32bit len2,
                  u32bit filled) {

  if (soaddHit) {
    (*soaddHit)(handle, orientation, id1, pos1, len1, id2, pos2, len2, filled);
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
            "M x ############ . %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %s "u32bitFMT"\n",
            name1, id1, pos1, len1, name2, id2, pos2, len2, (orientation == 'f') ? "1" : "-1", filled);

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
u64bit
filterObj::output(FILE *F, u64bit matchid) {
  if (sooutput) {
    matchid = (*sooutput)(handle, F, matchid);
  } else {
    char    matchIDstring[32] = {0};

    //  Insert the match id's for all these matches.  We have to
    //  do this here (not during searches) because we're threaded.
    //
    char *pos = theOutput;
    while (*pos) {

      //  Construct a string holding the text version of the match id.
      //
      matchid++;
      sprintf(matchIDstring, u64bitFMT, matchid);

      //  At the start of an output record.  Skip the row type and
      //  sub type, 'M x ', which should put us at the start of
      //  the match id.
      //
      pos += 4;

      //  Copy the number into the space, removing any extra #
      //  marks, warning if we run out of space.
      //
      char *matchIDiterator = matchIDstring;
      while ((*pos == '#') && (*matchIDiterator != 0))
        *pos++ = *matchIDiterator++;

      while (*pos == '#')
        *pos++ = ' ';

      if (*matchIDiterator != 0)
        fprintf(stderr, "WARNING:  there isn't enough space in the match to insert the match id "u64bitFMT" '%s'!\n",
                matchid, matchIDstring);

      //  Skip to the next record
      //
      while (*pos++ != '\n')
        ;
    }

    fwrite(theOutput, sizeof(char), theOutputPos, F);
  }

  return(matchid);
}


#endif  //  FILTEROBJ_H
