
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
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignment_CNS.C
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/MultiAlignment_CNS.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2009-MAY-14
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-08
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-14
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"

//  Shouldn't be global, but some things -- like abBaseCount -- need it.

bool    DATAINITIALIZED                     = false;

double  EPROB[CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };
double  PROB [CNS_MAX_QV - CNS_MIN_QV + 1]  = { 0 };

uint32  baseToIndex[256]                    = { 0 };
char    indexToBase[CNS_NUM_SYMBOLS]        = { 0 };

void
abAbacus::initializeGlobals(void) {

  for (int32 i=0; i<256; i++)
    baseToIndex[i] = UINT32_MAX;

  indexToBase[ 0] = '-';
  indexToBase[ 1] = 'A';
  indexToBase[ 2] = 'C';
  indexToBase[ 3] = 'G';
  indexToBase[ 4] = 'T';
  indexToBase[ 5] = 'N';
#if 0
  indexToBase[ 6] = 'a';  //  -A
  indexToBase[ 7] = 'c';  //  -C
  indexToBase[ 8] = 'g';  //  -G
  indexToBase[ 9] = 't';  //  -T
  indexToBase[10] = 'M';  //  AC
  indexToBase[11] = 'R';  //  AG
  indexToBase[12] = 'W';  //  AT
  indexToBase[13] = 'S';  //  CG
  indexToBase[14] = 'Y';  //  CT
  indexToBase[15] = 'K';  //  GT
  indexToBase[16] = 'm';  //  -AC
  indexToBase[17] = 'r';  //  -AG
  indexToBase[18] = 'w';  //  -AT
  indexToBase[19] = 's';  //  -CG
  indexToBase[20] = 'y';  //  -CT
  indexToBase[21] = 'k';  //  -GT
  indexToBase[22] = 'V';  //  ACG
  indexToBase[23] = 'H';  //  ACT
  indexToBase[24] = 'D';  //  AGT
  indexToBase[25] = 'B';  //  CGT
  indexToBase[26] = 'v';  //  -ACG
  indexToBase[27] = 'h';  //  -ACT
  indexToBase[28] = 'd';  //  -AGT
  indexToBase[29] = 'b';  //  -CGT
  indexToBase[30] = 'X';  //  ACGT
  indexToBase[31] = 'x';  //  -ACGT
#endif

  for (int32 i=0; i<CNS_NUM_SYMBOLS; i++)
    baseToIndex[indexToBase[i]] = i;

  baseToIndex['n'] = baseToIndex['N'];  //  Used in baseCount


  double TAU_MISMATCH = 1.0 / (5.0 - 1.0);

  for (int32 i=0, qv=CNS_MIN_QV; i<CNS_MAX_QV - CNS_MIN_QV + 1; i++, qv++) {
    EPROB[i]= log(TAU_MISMATCH * pow(10, -qv/10.0));
    PROB[i] = log(1.0 - pow(10, -qv/10.0));
  }

  DATAINITIALIZED = true;
}
