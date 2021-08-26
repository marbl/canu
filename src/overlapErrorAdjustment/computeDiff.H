
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

#include <utility>
#include "types.H"

static const uint32 DEFAULT_FLANK_IGNORE = 5;

inline
int32
TrimStartingIndels(char *&s, int32 &s_end, int32 *delta, int32 &deltaLen, int32 val) {
    int32  i = 0;

    //while (i < stop && ped->delta[i] == 1)
    while (i < deltaLen && delta[i] == val)
      i++;

    //fprintf(stderr, "RESET 1 i=%d delta=%d\n", i, ped->delta[i]);
    //assert(i == stop || ped->delta[i] != -1);

    if (i == 0)
      return 0;

    deltaLen -= i;
    memmove(delta, delta + i, deltaLen * sizeof(int32));

    s += i;
    s_end -= i;
    return i;
}

std::pair<std::size_t, std::size_t>
ComputeErrors(const char* const a_part, const char* const b_part,
    int32 delta_len, int32 *deltas,
    int32 a_len, int32 b_len,
    bool check_trivial_dna,
    uint32 ignore_flank = DEFAULT_FLANK_IGNORE);
