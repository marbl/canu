#ifndef SIM4PARAMETERS
#define SIM4PARAMETERS

#include "fasta.H"
#include "sequencemanager.H"

class sim4parameters {
public:
  sim4parameters(char *cmd, SequenceManager *ESTs, SequenceManager *DBs);
  sim4parameters(unsigned int ESTid,
                 unsigned int DBid,
                 SequenceManager *ESTs, SequenceManager *DBs);
  ~sim4parameters();

  //  COMMAND STRING
  //
  char           *_cmd;
  char          **_arg;

  //  EST SEQUENCE SPECIFICATION
  //
  unsigned int      _numESTs;
  unsigned int      _maxESTs;
  unsigned int     *_ESTlist;

  unsigned int      _dbIdx;
  unsigned int      _dbLo;
  unsigned int      _dbHi;

  SequenceManager  *_ESTs;
  SequenceManager  *_DBs;
  bool              _DBsIsPrivate;

  void            loadAll(void) {
    for (unsigned int i=0; i<_numESTs; i++)
      _ESTs->getSequence(_ESTlist[i]);

    _DBs->getSequence(_dbIdx);
  }

  unsigned char  *getESTsequence(unsigned int i) {
    return(_ESTs->getSequence(i));
  };
  unsigned char  *getDBsequence(unsigned int i) {
    return(_DBs->getSequence(i));
  };

  unsigned int    getESTlength(unsigned int i) {
    return(_ESTs->getSequenceLength(i));
  };
  unsigned int    getDBlength(unsigned int i) {
    return(_DBs->getSequenceLength(i));
  };

  unsigned char  *getESTheader(unsigned int i) {
    return(_ESTs->getHeader(i));
  };
  unsigned char  *getDBheader(unsigned int i) {
    return(_DBs->getHeader(i));
  };

  //  OPTIONS FOR SIM4
  //
  bool            _doForward;
  bool            _doReverse;

  char           *_strandIndicator;

  unsigned int    _accurateSequences;
  unsigned int    DRANGE;
  //unsigned int    weight;
  unsigned int    cutoff; 

};


#endif  //  SIM4PARAMETERS
