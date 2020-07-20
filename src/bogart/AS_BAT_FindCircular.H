
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

#ifndef AS_BAT_FIND_CIRCULAR_H
#define AS_BAT_FIND_CIRCULAR_H

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

bool
isCircularizingEdge(Unitig   *tig,
                    uint32    rdFid,
                    uint32    rdTid);


void
findCircularContigs(TigVector &tigs,
                    const char *prefix);

#endif  //  AS_BAT_FIND_CIRCULAR_H
