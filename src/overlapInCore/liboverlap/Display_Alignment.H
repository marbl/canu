
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

#ifndef DISPLAY_ALIGNMENT_H
#define DISPLAY_ALIGNMENT_H

#include "types.H"

void
Display_Alignment(char  *a,    int32   a_len,
                  char  *b,    int32   b_len,
                  int32 *delta,
                  int32  delta_ct);

#endif
