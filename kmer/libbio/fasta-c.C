#include "bio++.H"
#include "fasta-c.h"

typedef struct {
  FastABase           *F;
  FastASequenceInCore *B;
} _fasta_c_private;

void*
createFastA(char *file) {
  _fasta_c_private  *c = new _fasta_c_private;
  c->F = new FastAFile(file);
  c->F->openIndex();
  c->B = 0L;
  return(c);
}

void
destroyFastA(void *c) {
  _fasta_c_private *C = (_fasta_c_private *)c;
  delete C->F;
  delete C->B;
  delete C;
}

char *
getFastAsequence(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;

  if (C->B == 0L) {
    if (C->F->find(idx) == false) {
      fprintf(stderr, "Can't find iid=%u in %s\n", idx, C->F->getSourceName());
      exit(1);
    }
    C->B = C->F->getSequence();
  }

  if (C->B->getIID() != (u32bit)idx) {
    delete C->B;
    if (C->F->find(idx) == false) {
      fprintf(stderr, "Can't find iid=%u in %s\n", idx, C->F->getSourceName());
      exit(1);
    }
    C->B = C->F->getSequence();
  }

  return(C->B->sequence());
}

int
getFastAsequenceLength(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;

  return(C->F->sequenceLength(idx));
}

char *
getFastAheader(void *c, int idx) {
  _fasta_c_private *C = (_fasta_c_private *)c;

  if (C->B == 0L) {
    if (C->F->find(idx) == false) {
      fprintf(stderr, "Can't find iid=%u in %s\n", idx, C->F->getSourceName());
      exit(1);
    }
    C->B = C->F->getSequence();
  }

  if (C->B->getIID() != (u32bit)idx) {
    delete C->B;
    if (C->F->find(idx) == false) {
      fprintf(stderr, "Can't find iid=%u in %s\n", idx, C->F->getSourceName());
      exit(1);
    }
    C->B = C->F->getSequence();
  }

  return(C->B->header());
}
