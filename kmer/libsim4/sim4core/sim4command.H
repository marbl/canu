#ifndef SIM4COMMAND_H
#define SIM4COMMAND_H

#include "libbri.H"
#include "fastafile.H"

//
//  Contains the variable stuff for an execution of sim4
//    access to sequences (via FastAFileWrappers)
//    genomic iid
//    genomic range
//    est iid (maybe more than one)
//    forward only
//    reverse only
//    an optional strand indicator
//

class sim4command {
public:

  //  Run a single EST against a genomic range
  //
  sim4command(u32bit            ESTid,
              FastAFileWrapper *ESTs,
              u32bit            GENid,
              u32bit            GENlo,
              u32bit            GENhi,
              FastAFileWrapper *GENs,
              bool              doForward,
              bool              doReverse);

  //  Single EST against a genomic range, using alternative
  //  interface.
  //
  sim4command(FastASequenceInCore  *EST,
              FastASequenceInCore  *GEN,
              u32bit                GENlo,
              u32bit                GENhi,
              bool                  doForward,
              bool                  doReverse);

  //  Run multiple ESTs against a genomic range
  //
#if 0
  sim4command(FastAFileWrapper *ESTs,
              FastAFileWrapper *GENs);
#endif

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

  void   addESTid(u32bit id);

  u32bit          getNumberOfESTs(void) {
    if (_ESTsequence)
      return(1);
    else
      return(_ESTlistLen);
  };

  u32bit          getESTidx(u32bit i);
  char           *getESTheader(u32bit i);
  char           *getESTsequence(u32bit i);
  u32bit          getESTlength(u32bit i);

  u32bit          getGENidx(void)           { return(_genIdx); };
  u32bit          getGENlo(void)            { return(_genLo); };
  u32bit          getGENhi(void)            { return(_genHi); };
  char           *getGENheader(void);
  char           *getGENsequence(void);
  u32bit          getGENlength(void);

  void            setStrandIndicator(char *si) {
    _strandIndicator = new char [strlen(si) + 1];
    strcpy(_strandIndicator, si);
  }
  char           *getStrandIndicator(void) {
    return(_strandIndicator);
  };
private:
  void            loadEST(u32bit i);
  void            loadGEN(void);

  u32bit                _ESTlistLen;
  u32bit                _ESTlistMax;
  u32bit               *_ESTlist;

  u32bit                _genIdx;
  u32bit                _genLo;
  u32bit                _genHi;

  FastAFileWrapper     *_ESTs;
  FastASequenceInCore  *_ESTloaded;
  char                 *_ESTsequence;
  u32bit                _ESTsequenceLength;  //  valid only for _ESTsequence

  FastAFileWrapper     *_GENs;
  FastASequenceInCore  *_GENloaded;
  char                 *_GENsequence;
  u32bit                _GENsequenceLength;

  bool                  _doForward;
  bool                  _doReverse;

  char                 *_strandIndicator;
};


#endif  //  SIM4COMMAND_H
