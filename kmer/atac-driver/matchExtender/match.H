#ifndef MATCH_H
#define MATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "bio++.H"

class match_s {
public:
  char    _matchId[32];

  FastASequenceInCore  *_seq1;
  FastAAccessor  *_acc1;
  char    _id1[32];
  u32bit  _iid1;
  u32bit  _pos1;
  u32bit  _len1;
  u32bit  _ori1;

  FastASequenceInCore  *_seq2;
  FastAAccessor  *_acc2;
  char    _id2[32];
  u32bit  _iid2;
  u32bit  _pos2;
  u32bit  _len2;
  u32bit  _ori2;

  bool    _isDeleted;

  //  Our diagonal never, ever changes
  u32bit  _diagonal;

public:
  match_s(char *matchId,
          FastASequenceInCore *seq1, char *id1, u32bit iid1, u32bit pos1, u32bit len1, u32bit ori1,
          FastASequenceInCore *seq2, char *id2, u32bit iid2, u32bit pos2, u32bit len2, u32bit ori2) {

    strncpy(_matchId, matchId, 32);

    _seq1 = seq1;
    _acc1 = new FastAAccessor(_seq1, false);
    strncpy(_id1, id1, 32);
    _iid1 = iid1;
    _pos1 = pos1;
    _len1 = len1;
    _ori1 = ori1;

    _seq2 = seq2;
    _acc2 = new FastAAccessor(_seq2, (ori1 != ori2));
    strncpy(_id2, id2, 32);
    _iid2 = iid2;
    _pos2 = pos2;
    _len2 = len2;
    _ori2 = ori2;

    _isDeleted = false;

    _acc1->setRange(_pos1, _len1);
    _acc2->setRange(_pos2, _len2);

    _acc1->setPosition(_pos1);
    _acc2->setPosition(_pos2);

    //  the diagonal is....
    if (ori1 == ori2)
      _diagonal = seq1->sequenceLength() - _pos1 + _pos2;
    else
      _diagonal = seq1->sequenceLength() - _pos1 + seq2->sequenceLength() - (_pos2 + _len2);
    //_diagonal = seq2->sequenceLength() - (_pos2 + _len2) + _pos1;
  };

  //  Compare by diagonal, then by position in the first sequence.
  //
  bool operator<(const match_s& r) const {
    if (_diagonal  < r._diagonal)  return(true);
    if (_diagonal == r._diagonal)  return(_pos1 < r._pos1);
    return(false);
  };

  void dump(FILE *out, const char *descr, bool showSeq=false);


  //  For compatibility

  void extendLeft(s32bit num) {
    _acc1->extendLeft(num);
    _acc2->extendLeft(num);
  };

  void extendRight(s32bit num) {
    _acc1->extendRight(num);
    _acc2->extendRight(num);
  };

  bool isDeleted(void)           { return(_isDeleted); };
  void setDeleted(void)          { _isDeleted = true;  };


  //  Since we're on the same diagonal, and ungapped, the choice of
  //  testing sequence 1 or 2 is arbirary.
  //
  bool canMergeWith(match_s *m) {
    return((m != 0L) &&
           (_diagonal == m->_diagonal) &&
           (_acc1->getRangeBegin() + _acc1->getRangeEnd() >= m->_acc1->getRangeBegin()));
  };

  //  Extend us to end where m ends
  //
  void consume(match_s *m) {
    if (m->_acc1->getRangeEnd() > _acc1->getRangeEnd())
      extendRight(m->_acc1->getRangeEnd() - _acc1->getRangeEnd());
  };

  u32bit  len(void) { return(_acc1->getRangeEnd() - _acc1->getRangeBegin()); };
  u32bit  pos1(void) { return(_acc1->getRangeBegin()); };
  u32bit  pos2(void) { return(_acc2->getRangeBegin()); };
  FastASequenceInCore  *s1(void) { return(_seq1); };
  FastASequenceInCore  *s2(void) { return(_seq2); };
};

#endif  //  MATCH_H
