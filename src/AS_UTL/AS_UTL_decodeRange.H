
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
 *    Brian P. Walenz from 2012-FEB-12 to 2013-OCT-11
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2015-MAY-28
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_DECODERANGE_H
#define AS_UTL_DECODERANGE_H

#include "AS_global.H"

#include <set>

using namespace std;

void  AS_UTL_decodeRange(char *range, set<uint64> &ranges);
void  AS_UTL_decodeRange(char *range, set<uint32> &ranges);

void  AS_UTL_decodeRange(char *range, uint64 &lo, uint64 &hi);
void  AS_UTL_decodeRange(char *range,  int64 &lo,  int64 &hi);
void  AS_UTL_decodeRange(char *range, uint32 &lo, uint32 &hi);
void  AS_UTL_decodeRange(char *range,  int32 &lo,  int32 &hi);

void  AS_UTL_decodeRange(char *range, double &lo, double &hi);

#endif  //  AS_UTL_DECODERANGE_H
