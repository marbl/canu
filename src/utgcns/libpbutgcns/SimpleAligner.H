#ifndef __GCON_SIMPLE_ALIGNER__
#define __GCON_SIMPLE_ALIGNER__

#include "dw.H"

class SimpleAligner {
public:
    SimpleAligner();
    void align(dagcon::Alignment& aln, double errorRate);
};

#endif // __GCON_SIMPLE_ALIGNER__
