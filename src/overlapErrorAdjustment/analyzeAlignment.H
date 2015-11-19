
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-18 to 2015-JUL-01
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"

#include "correctionOutput.H"


//  Allocate one per base in read under evaluation.
//    5 values per 64 bit word -> 12 bits per value + 4 left over
//    6 values per 64 bit word -> 10 bits per value + 4 left over
//    7 values per 64 bit word ->  9 bits per value + 1 left over
//
//  Since there are only 11 values to store, we use the 10 bit size.
//
//  This is only in-memory, so we don't worry about padding it to maintain
//  proper sizes between runs.

#define VOTETALLY_BITS  10
#define VOTETALLY_MAX   1023

struct Vote_Tally_t {
  uint64  confirmed : VOTETALLY_BITS;
  uint64  deletes   : VOTETALLY_BITS;
  uint64  a_subst   : VOTETALLY_BITS;
  uint64  c_subst   : VOTETALLY_BITS;
  uint64  g_subst   : VOTETALLY_BITS;
  uint64  t_subst   : VOTETALLY_BITS;

  uint64  no_insert : VOTETALLY_BITS;
  uint64  a_insert  : VOTETALLY_BITS;
  uint64  c_insert  : VOTETALLY_BITS;
  uint64  g_insert  : VOTETALLY_BITS;
  uint64  t_insert  : VOTETALLY_BITS;
};


//  Maps an original location to a fixed location.  Supposedly
//  smaller/faster than the obvious array mapping...
struct Adjust_t {
  int32  adjpos;
  int32  adjust;
};


class analyzeAlignment {
public:
  analyzeAlignment() {
    Degree_Threshold   = 2;       //  ??
    Use_Haplo_Ct       = true;    //  Use haplotype counts to correct
    End_Exclude_Len    = 3;       //  ??
    Kmer_Len           = 9;       //  ??
    Vote_Qualify_Len   = 9;       //  ??
    Min_Haplo_Occurs   = 3;       //  This many or more votes for the same base indicates a haplotype

    _readID           = 0;

    _seqLen           = 0;
    _seq              = NULL;

    _corSeqLen        = 0;
    _corSeqMax        = AS_MAX_READLEN;
    _corSeq           = new char [_corSeqMax];

    _voteMax          = 0;
    _vote             = NULL;

    _lDegree          = 0;
    _rDegree          = 0;

    _readSub          = new int32        [AS_MAX_READLEN];
    _algnSub          = new int32        [AS_MAX_READLEN];
    _voteValue        = new Vote_Value_t [AS_MAX_READLEN];

    _corLen           = 0;
    _corMax           = AS_MAX_READLEN / 4;
    _cor              = new Correction_Output_t [_corMax];

    //  The original version converted ACGT to lowercase, and replaced non-ACGT with 'a'.

    for (uint32 i=0; i<256; i++)
      _filter[i] = 'a';

    _filter['A'] = _filter['a'] = 'a';
    _filter['C'] = _filter['c'] = 'c';
    _filter['G'] = _filter['g'] = 'g';
    _filter['T'] = _filter['t'] = 't';
  };


  ~analyzeAlignment() {
    delete [] _corSeq;
    delete [] _vote;
    delete [] _readSub;
    delete [] _algnSub;
    delete [] _voteValue;
    delete [] _cor;
  };

public:
  void   reset(uint32 id, char *seq, uint32 seqLen) {
    fprintf(stderr, "reset() for read id %u of length %u\n", id, seqLen);
    _readID    = id;

    _seqLen    = seqLen;
    _seq       = seq;

    _corSeqLen = 0;

    if (_voteMax < _seqLen) {
      delete [] _vote;

      _voteMax  = _seqLen + 2000;
      _vote     = new Vote_Tally_t [_voteMax];
    }

    //  Seems like overkill, but I don't know where to put it
    memset(_vote, 0, sizeof(Vote_Tally_t) * _voteMax);
  };

  void   analyze(char  *a, int32 aLen,  int32 aOffset,
                 char  *b, int32 bLen,
                 int32  deltaLen,
                 int32 *delta);

  void   outputDetails(uint32 j);
  void   outputDetails(void);

  void   generateCorrections(FILE *corFile=NULL);

  void   generateCorrectedRead(Adjust_t *fadj=NULL, uint32 *fadjLen=NULL,
                               uint64   *changes=NULL);


private:
  void
  castVote(Vote_Value_t  val,
           int32         pos) {
    int32 v = 0;

    switch (val) {
      case DELETE:    if (_vote[pos].deletes  < VOTETALLY_MAX)  v = ++_vote[pos].deletes;   break;
      case A_SUBST:   if (_vote[pos].a_subst  < VOTETALLY_MAX)  v = ++_vote[pos].a_subst;   break;
      case C_SUBST:   if (_vote[pos].c_subst  < VOTETALLY_MAX)  v = ++_vote[pos].c_subst;   break;
      case G_SUBST:   if (_vote[pos].g_subst  < VOTETALLY_MAX)  v = ++_vote[pos].g_subst;   break;
      case T_SUBST:   if (_vote[pos].t_subst  < VOTETALLY_MAX)  v = ++_vote[pos].t_subst;   break;
      case A_INSERT:  if (_vote[pos].a_insert < VOTETALLY_MAX)  v = ++_vote[pos].a_insert;  break;
      case C_INSERT:  if (_vote[pos].c_insert < VOTETALLY_MAX)  v = ++_vote[pos].c_insert;  break;
      case G_INSERT:  if (_vote[pos].g_insert < VOTETALLY_MAX)  v = ++_vote[pos].g_insert;  break;
      case T_INSERT:  if (_vote[pos].t_insert < VOTETALLY_MAX)  v = ++_vote[pos].t_insert;  break;
      case NO_VOTE:
        break;
      default :
        fprintf(stderr, "ERROR:  Illegal vote type\n");
        break;
    }
  };



private:

  //  Parameters

  int32                   Degree_Threshold;
  int32                   Use_Haplo_Ct;
  int32                   End_Exclude_Len;
  int32                   Kmer_Len;
  int32                   Vote_Qualify_Len;
  int32                   Min_Haplo_Occurs;

  //  Per-read data

  uint32                 _readID;

  uint32                 _seqLen;
  char                  *_seq;

public:
  uint32                 _corSeqLen;
  uint32                 _corSeqMax;
  char                  *_corSeq;

private:
  uint32                 _voteMax;
  Vote_Tally_t          *_vote;

  uint32                 _lDegree;
  uint32                 _rDegree;

  int32                 *_readSub;
  int32                 *_algnSub;
  Vote_Value_t          *_voteValue;

  //  Outputs (used to generate the corrected read)

  uint32                 _corLen;
  uint32                 _corMax;
  Correction_Output_t   *_cor;

  char                   _filter[256];
};
