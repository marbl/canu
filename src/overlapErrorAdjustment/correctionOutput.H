
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

//  Definitions for our exportable data.

enum Vote_Value_t {
  IDENT,    //  Just an iid in this record.
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


struct Correction_Output_t {
  uint32  keep_left   :  1;    //  set true if left overlap degree is low
  uint32  keep_right  :  1;    //  set true if right overlap degree is low
  uint32  type        :  4;    //  Vote_Value_t
  uint32  pos         :  26;   //
  uint32  readID;
};


