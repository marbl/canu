#ifndef SIM4COMMAND_H
#define SIM4COMMAND_H

#include "libbri.H"
#include "fasta.H"

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
};


#endif  //  SIM4COMMAND_H
