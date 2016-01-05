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
