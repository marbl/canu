#ifndef FASTA_SEQUENCEINCORE_H
#define FASTA_SEQUENCEINCORE_H

#include "bri.h"
#include "fasta-types.H"

class FastASequenceInCore {
public:
  FastASequenceInCore(IID_t iid);
  FastASequenceInCore(IID_t iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen);
  ~FastASequenceInCore();

  char           *header(void) const         { return(_header); };
  u32bit          headerLength(void) const   { return(_headerLen); };

  char           *sequence(void) const       { return(_seq); };
  u32bit          sequenceLength(void) const { return(_seqLen); };

  u32bit          getIID(void) const         { return(_idx); };
private:
  u32bit         _headerLen;
  u32bit         _headerMax;
  char          *_header;

  u32bit         _seqLen;
  u32bit         _seqMax;
  char          *_seq;

  u32bit         _idx;
};

#endif  //  FASTA_SEQUENCEINCORE_H
