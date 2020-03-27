
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#pragma once

//  Definitions for our exportable data.

enum Vote_Value_t {
  IDENT = 0,    //  Just an iid in this record.
  DELETE,
  A_SUBST,
  C_SUBST,
  G_SUBST,
  T_SUBST,  //  DON'T rearrange this!  Code depends on the ordering.
  A_INSERT,
  C_INSERT,
  G_INSERT,
  T_INSERT,
  NO_VOTE,
  EXTENSION
};

static char
VoteChar(Vote_Value_t val) {
  switch (val) {
    case A_INSERT:
    case A_SUBST:
      return 'a';
    case C_INSERT:
    case C_SUBST:
      return 'c';
    case G_INSERT:
    case G_SUBST:
      return 'g';
    case T_INSERT:
    case T_SUBST:
      return 't';
    default:
      assert(false);
  }
  return 'a';
}

static Vote_Value_t
InsVote(char c) {
  switch (c) {
    case 'a':
      return A_INSERT;
    case 'c':
      return C_INSERT;
    case 'g':
      return G_INSERT;
    case 't':
      return T_INSERT;
    default:
      fprintf(stderr, "ERROR: Bad sequence '%c' 0x%02x)\n", c, c);
      assert(false);
  }
  return NO_VOTE;
}

static Vote_Value_t
SubstVote(char c) {
  switch (c) {
    case 'a':
      return A_SUBST;
    case 'c':
      return C_SUBST;
    case 'g':
      return G_SUBST;
    case 't':
      return T_SUBST;
    default:
      fprintf(stderr, "ERROR: Bad sequence '%c' 0x%02x)\n", c, c);
      assert(false);
  }
  return NO_VOTE;
}


struct Correction_Output_t {
  bool  keep_left   :  1;    //  set true if left overlap degree is low
  bool  keep_right  :  1;    //  set true if right overlap degree is low
  Vote_Value_t  type        :  4;    //  Vote_Value_t
  uint32  pos         :  26;   //
  uint32  readID;
};
