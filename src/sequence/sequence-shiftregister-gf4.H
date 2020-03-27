
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

#ifndef SEQUENCE_SHIFTREGISTER_GF4
#define SEQUENCE_SHIFTREGISTER_GF4

#include "types.H"


//  Addition and multiplication tables for GF(4).

const
uint32
gf4add[4][4] =  { { 0, 1, 2, 3 },
                  { 1, 0, 3, 2 },
                  { 2, 3, 0, 1 },
                  { 3, 2, 1, 0 } };

const
uint32
gf4mult[4][4] = { { 0, 0, 0, 0 },
                  { 0, 1, 2, 3 },
                  { 0, 2, 3, 1 },
                  { 0, 3, 1, 2 } };


#endif  //  SEQUENCE_SHIFTREGISTER_GF4
