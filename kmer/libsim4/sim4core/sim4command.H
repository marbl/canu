#ifndef SIM4COMMAND_H
#define SIM4COMMAND_H

#include "bio++.H"

#include "seqCache.H"

//
//  Contains the variable stuff for an execution of sim4
//    access to sequences (via seqCache)
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
  sim4command(uint32        ESTid,
              seqCache     *ESTs,
              uint32        GENid,
              uint32        GENlo,
              uint32        GENhi,
              seqCache     *GENs,
              bool          doForward,
              bool          doReverse);

  //  Single EST against a genomic range, using alternative
  //  interface.
  //
  sim4command(seqInCore  *EST,
              seqInCore  *GEN,
              uint32      GENlo,
              uint32      GENhi,
              bool        doForward,
              bool        doReverse);

  //  Use two char*'s for sequence sources -- both sequence deflines
  //  and iid's are undefined!
  //
  sim4command(char             *EST,
              uint32            ESTlen,
              char             *GEN,
              uint32            GENlen,
              uint32            GENlo,
              uint32            GENhi,
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
  void   addSeed(uint32 GENpos, uint32 ESTpos, uint32 length);

  void   sortExternalSeeds(void);
  bool   externalSeedsExist(void)          { return(_externalSeedsLen > 0); };
  uint32 numberOfExternalSeeds(void)       { return(_externalSeedsLen); };
  uint32 externalSeedESTPosition(uint32 i) { return(_externalSeeds[i]._ESTposition); };
  uint32 externalSeedGENPosition(uint32 i) { return(_externalSeeds[i]._GENposition); };
  uint32 externalSeedLength(uint32 i)      { return(_externalSeeds[i]._length);      };
  void   maskExternalSeed(uint32 i)        { _externalSeeds[i]._length = 0; };


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

  void   setGenomic(uint32 idx, uint32 lo, uint32 hi) {
    _genIdx = idx;
    _genLo  = lo;
    _genHi  = hi;
  };

  uint32          getESTidx();
  char           *getESTheader();
  char           *getESTsequence();
  uint32          getESTlength();

  uint32          getGENidx(void)           { return(_genIdx); };
  uint32          getGENlo(void)            { return(_genLo); };
  uint32          getGENhi(void)            { return(_genHi); };
  char           *getGENheader(void);
  char           *getGENsequence(void);
  uint32          getGENlength(void);
private:
  void            loadEST(void);
  void            loadGEN(void);

  uint32                _estIdx;

  seqCache             *_ESTs;
  seqInCore            *_ESTloaded;
  char                 *_ESTsequence;
  uint32                _ESTsequenceLength;  //  valid only for _ESTsequence

  uint32                _genIdx;
  uint32                _genLo;
  uint32                _genHi;

  seqCache             *_GENs;
  seqInCore            *_GENloaded;
  char                 *_GENsequence;
  uint32                _GENsequenceLength;

  bool                  _doForward;
  bool                  _doReverse;

  //  For external seeding
  //
  class externalSeed {
  public:
    uint32    _GENposition;
    uint32    _ESTposition;
    uint32    _length;

    bool       operator<(const externalSeed &that) const {
      return(_GENposition < that._GENposition);
    };
  };

  uint32                _externalSeedsLen;
  uint32                _externalSeedsMax;
  externalSeed         *_externalSeeds;
};


#endif  //  SIM4COMMAND_H
