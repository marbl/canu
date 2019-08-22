
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */
#pragma once

#include <string>
#include <iostream>

//  Definitions for our exportable data.

//TODO add #pragma once
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

//TODO fix types
struct Correction_Output_t {
  bool  keep_left;   //:  1;    //  set true if left overlap degree is low
  bool  keep_right;  //:  1;    //  set true if right overlap degree is low
  uint32  type;        //:  4;    //  Vote_Value_t
  uint32  pos;         //:  26;   //
  uint32  readID;
  std::string insertion_str;
};

namespace std {
  static std::ostream & operator<<(std::ostream &os, const Correction_Output_t &corr_out) {
    os << corr_out.keep_left << ' ' << corr_out.keep_right << ' ';
    os << corr_out.type << ' ' << corr_out.pos << ' ';
    os << corr_out.readID;
    if (corr_out.type == EXTENSION) {
      assert(!corr_out.insertion_str.empty());
      os << ' ' << corr_out.insertion_str;
    }
    return os;
  }

  static std::istream & operator>>(std::istream &is, Correction_Output_t &corr_out) {
    is >> corr_out.keep_left >> corr_out.keep_right;
    is >> corr_out.type >> corr_out.pos;
    is >> corr_out.readID;
    if (corr_out.type == EXTENSION) {
      is >> corr_out.insertion_str;
    }
    return is;
  }
}

