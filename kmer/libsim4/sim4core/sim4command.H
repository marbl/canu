#ifndef SIM4COMMAND_H
#define SIM4COMMAND_H

#include "bio++.H"

//
//  Contains the variable stuff for an execution of sim4
//    access to sequences (via FastAWrappers)
//    genomic iid
//    genomic range
//    est iid (maybe more than one)
//    forward only
//    reverse only
//

class sim4command {
public:

  //  Run a single EST against a genomic range
  //
  sim4command(u32bit        ESTid,
              FastAWrapper *ESTs,
              u32bit        GENid,
              u32bit        GENlo,
              u32bit        GENhi,
              FastAWrapper *GENs,
              bool          doForward,
              bool          doReverse);

  //  Single EST against a genomic range, using alternative
  //  interface.
  //
  sim4command(FastASequenceInCore  *EST,
              FastASequenceInCore  *GEN,
              u32bit                GENlo,
              u32bit                GENhi,
              bool                  doForward,
              bool                  doReverse);

  //  Use two char*'s for sequence sources -- both sequence deflines
  //  and iid's are undefined!
  //
  sim4command(char             *EST,
              u32bit            ESTlen,
              char             *GEN,
              u32bit            GENlen,
              u32bit            GENlo,
              u32bit            GENhi,
              bool              doForward,
              bool              doReverse);

  ~sim4command();

  //  These methods allow the initial seed detection to be
  //  done outside Sim4::run().  If used:
  //    each seed is extended as before.
  //    for interative alignments, seeds are masked out
  //

  //  addSeed() takes coordinates relative to the start of the GEN
  //  sequence supplied to the constructor.
  //
  void   addSeed(u32bit GENpos, u32bit ESTpos, u32bit length);

  void   sortExternalSeeds(void);
  bool   externalSeedsExist(void)          { return(_externalSeedsLen > 0); };
  u32bit numberOfExternalSeeds(void)       { return(_externalSeedsLen); };
  u32bit externalSeedESTPosition(u32bit i) { return(_externalSeeds[i]._ESTposition); };
  u32bit externalSeedGENPosition(u32bit i) { return(_externalSeeds[i]._GENposition); };
  u32bit externalSeedLength(u32bit i)      { return(_externalSeeds[i]._length);      };
  void   maskExternalSeed(u32bit i)        { _externalSeeds[i]._length = 0; };


  //  Load the sequences, make some checks.  This isn't done in the
  //  constructor so that it is possible to make a big list of
  //  commands, then give them to a processor.  If we loaded all the
  //  genomics at creation....
  //
  void   finalize(void);

  void   setForward(bool x) { _doForward = x; };
  void   setReverse(bool x) { _doReverse = x; };

  bool   doForward(void)    { return(_doForward); };
  bool   doReverse(void)    { return(_doReverse); };

  void   setGenomic(u32bit idx, u32bit lo, u32bit hi) {
    _genIdx = idx;
    _genLo  = lo;
    _genHi  = hi;
  };

  u32bit          getESTidx();
  char           *getESTheader();
  char           *getESTsequence();
  u32bit          getESTlength();

  u32bit          getGENidx(void)           { return(_genIdx); };
  u32bit          getGENlo(void)            { return(_genLo); };
  u32bit          getGENhi(void)            { return(_genHi); };
  char           *getGENheader(void);
  char           *getGENsequence(void);
  u32bit          getGENlength(void);
private:
  void            loadEST(void);
  void            loadGEN(void);

  u32bit                _estIdx;

  FastAWrapper         *_ESTs;
  FastASequenceInCore  *_ESTloaded;
  char                 *_ESTsequence;
  u32bit                _ESTsequenceLength;  //  valid only for _ESTsequence

  u32bit                _genIdx;
  u32bit                _genLo;
  u32bit                _genHi;

  FastAWrapper         *_GENs;
  FastASequenceInCore  *_GENloaded;
  char                 *_GENsequence;
  u32bit                _GENsequenceLength;

  bool                  _doForward;
  bool                  _doReverse;

  //  For external seeding
  //
  struct externalSeed {
    u32bit    _GENposition;
    u32bit    _ESTposition;
    u32bit    _length;
  };

  int                   _externalSeedsCompare(const void *a, const void *b) {
    externalSeed *A = *((externalSeed **)a);
    externalSeed *B = *((externalSeed **)b);

    if (A->_GENposition < B->_GENposition)
      return(-1);
    if (A->_GENposition > B->_GENposition)
      return(1);
    return(0);
  };

  u32bit                _externalSeedsLen;
  u32bit                _externalSeedsMax;
  externalSeed         *_externalSeeds;

  friend int sortExternalSeedsCompare(const void *, const void *);
};


#endif  //  SIM4COMMAND_H
