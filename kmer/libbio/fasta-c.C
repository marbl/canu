#include "fasta-c.h"
#include "libbri.H"

//
//  A limited FastA reader, callable from C
//

typedef struct {
  FastA        *F;
  FastABuffer   B;
} _fasta_c_private;

void *createFastA(char *file) {
  _fasta_c_private  *c = new _fasta_c_private;
  c->F    = new FastA(file, true, false);
  return(c);
}

unsigned char  *getFastAsequence(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  C->F->seek(C->B, idx);
  return(C->B.sequence());
}

unsigned char  *getFastAsequenceReverseComplement(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  C->F->seek(C->B, idx);
  C->B.reverseComplementSequence();
  return(C->B.sequence());
}

int getFastAsequenceLength(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  return(C->F->sequenceLength(idx));
}

unsigned char  *getFastAheader(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  C->F->seek(C->B, idx);
  return(C->B.header());
}

void  destroyFastA(void *c) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  delete C->F;
  delete C;
}

