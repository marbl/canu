
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
 *    Sergey Koren beginning on 2015-DEC-28
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Brian P. Walenz beginning on 2016-JAN-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef __GCON_ALIGNMENT_HPP__
#define __GCON_ALIGNMENT_HPP__

///
/// Super-simple alignment representation.  Represents an alignment between two
/// PacBio reads, one of which we're trying to correct.  The read to correct
/// may be either the target or the query, depending on how the alignment was
/// done.
///

#include <string>
#include <stdint.h>

namespace dagcon {
class Alignment {
public:
    // conforming offsets are 1-based
    uint32_t start;

    uint32_t end;

    // Fragment ID
    std::string frgid;

    // query and target strings must be equal length
    std::string qstr;
    std::string tstr;

    Alignment();
};
}

struct Unitig {
    uint32_t id;
    std::string seq;
};

/// Simplifies the alignment by normalizing gaps.  Converts mismatches into
/// indels ...
///      query: CAC        query:  C-AC
///             | |  --->          |  |
///     target: CGC       target:  CG-C
///
/// Shifts equivalent gaps to the right in the reference ...
///      query: CAACAT        query: CAACAT
///             | | ||  --->         |||  |
///     target: C-A-AT       target: CAA--T
///
/// Shifts equivalent gaps to the right in the read ...
///      query: -C--CGT       query: CCG--T
///              |  | |  --->        |||  |
///     target: CCGAC-T      target: CCGACT
dagcon::Alignment normalizeGaps(dagcon::Alignment& aln);

#endif // __GCON_ALIGNMENT_HPP__
