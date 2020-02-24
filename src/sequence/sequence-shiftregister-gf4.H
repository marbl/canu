
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
 *    Brian P. Walenz beginning on 2020-FEB-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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
