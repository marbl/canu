#ifndef MATCH_H
#define MATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "atac.H"
#include "bio++.H"
#include "seqCache.H"

class match_s {
public:
  char    _matchId[32];

  seqInCore  *_seq1;
  FastAAccessor  *_acc1;
  u32bit  _iid1;
  u32bit  _pos1;
  u32bit  _len1;
  u32bit  _ori1;

  seqInCore  *_seq2;
  FastAAccessor  *_acc2;
  u32bit  _iid2;
  u32bit  _pos2;
  u32bit  _len2;
  u32bit  _ori2;

  bool    _isDeleted;

  //  Our diagonal never, ever changes
  u32bit  _diagonal;

public:
  match_s(char *matchId,
          seqInCore *s1, u32bit i1, u32bit p1, u32bit l1, u32bit o1,
          seqInCore *s2, u32bit i2, u32bit p2, u32bit l2, u32bit o2) {

    strncpy(_matchId, matchId, 32);

    _seq1 = s1;
    _acc1 = new FastAAccessor(_seq1, false);
    _iid1 = i1;
    _pos1 = p1;
    _len1 = l1;
    _ori1 = o1;

    _seq2 = s2;
    _acc2 = new FastAAccessor(_seq2, (o1 != o2));
    _iid2 = i2;
    _pos2 = p2;
    _len2 = l2;
    _ori2 = o2;

    _isDeleted = false;

    _acc1->setRange(_pos1, _len1);
    _acc2->setRange(_pos2, _len2);

    _acc1->setPosition(_pos1);
    _acc2->setPosition(_pos2);

    //  the diagonal is....
    if (_ori1 == _ori2)
      _diagonal = _seq1->sequenceLength() - _pos1 + _pos2;
    else
      _diagonal = _seq1->sequenceLength() - _pos1 + _seq2->sequenceLength() - (_pos2 + _len2);
  };

  ~match_s() {
    delete _acc1;
    delete _acc2;
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
           (_acc1->getRangeEnd() >= m->_acc1->getRangeBegin()));
  };

  //  Extend us to end where m ends
  //
  void consume(match_s *m) {
    if (m->_acc1->getRangeEnd() > _acc1->getRangeEnd())
      extendRight(m->_acc1->getRangeEnd() - _acc1->getRangeEnd());
  };

  u32bit      len(void)  { return(_acc1->getRangeEnd() - _acc1->getRangeBegin()); };
  u32bit      pos1(void) { return(_acc1->getRangeBegin()); };
  u32bit      pos2(void) { return(_acc2->getRangeBegin()); };
  seqInCore  *seq1(void) { return(_seq1); };
  seqInCore  *seq2(void) { return(_seq2); };
};

#endif  //  MATCH_H
