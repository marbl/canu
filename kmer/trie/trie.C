#include "util++.H"
#include "bio++.H"

class trieNode {
public:
  trieNode(void) {
    next[0] = ~u32bitZERO;
    next[1] = ~u32bitZERO;
    next[2] = ~u32bitZERO;
    next[3] = ~u32bitZERO;
    numseq  = 0;
    seqptr  = ~u32bitZERO;
  };

  u32bit   next[4];  //  next node for A, C, G, T input, ~u32bitZERO if no next
  u32bit   numseq;   //  number of seqs we have in here
  u32bit   seqptr;   //  pointer to seqs
};


class trieSeqPtr {
public:
  trieSeqPtr(void) {
    seqiid   = ~u32bitZERO;
    nodeiid  = ~u32bitZERO;
    defline  = 0L;
    reversed = false;
  };

  u32bit   seqiid;
  u32bit   nodeiid;
  char    *defline;
  bool     reversed;
};




int
trieSeqPtrCompare(const void *a, const void *b) {
  const trieSeqPtr *A = (const trieSeqPtr *)a;
  const trieSeqPtr *B = (const trieSeqPtr *)b;

  if (A->nodeiid < B->nodeiid)
    return(-1);
  if (A->nodeiid > B->nodeiid)
    return(1);
  return(0);
}


u32bit
addSequence(trieNode *nodes,    u32bit &nodesLen,
            trieSeqPtr *seqptr, u32bit &seqptrLen,
            FastASequenceInCore *S,
            bool isReverse) {
  char   *s = 0L;
  u32bit  n = 0;

  if (S->sequenceLength() < 12)
    return(0);

  for (s = S->sequence(); *s; s++)
    if (validSymbol[(int)*s] == 0)
      return(0);

  for (s = S->sequence(); *s; s++) {
    u32bit  v = compressSymbol[(int)*s];

    //  add a new pointer if needed
    if (nodes[n].next[v] == ~u32bitZERO)
      nodes[n].next[v] = nodesLen++;

    //  Go there
    n = nodes[n].next[v];
  }

  //  add this sequence to node i -- after all sequences have been
  //  added, we'll sort this list and build pointers.

  seqptr[seqptrLen].seqiid   = S->getIID();
  seqptr[seqptrLen].nodeiid  = n;
  seqptr[seqptrLen].defline  = strdup(S->header());
  seqptr[seqptrLen].reversed = isReverse;
  seqptrLen++;

  return(1);
}



int
main(int argc, char **argv) {
  char  *queries = 0L;
  char  *genome  = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-q") == 0) {
      queries = argv[++arg];
    } else if (strcmp(argv[arg], "-g") == 0) {
      genome = argv[++arg];
    } else {
      err++;
    }
    arg++;
  }

  if (queries == 0L)
    err = 1;
  if (genome == 0L)
    err = 1;

  if (err) {
    fprintf(stderr, "usage: %s -q queries.fasta -g genome.fasta\n", argv[0]);
    exit(1);
  }

  //  Build the trie from the queries.  We just allocate a large
  //  number of nodes so we don't need to deal with reallocation.
  //
  //  32M * 28 = 896M.
  //
  u32bit     nodesLen = 1;
  u32bit     nodesMax = 16 * 1024 * 1024;
  trieNode  *nodes    = new trieNode [nodesMax];

  FastAWrapper         *F = new FastAWrapper(queries);
  FastASequenceInCore  *S = 0L;

  u32bit      seqptrLen = 0;
  u32bit      seqptrMax = 2 * 1024 * 1024;
  trieSeqPtr *seqptr    = new trieSeqPtr [seqptrMax];

  while ((S = F->getSequence()) != 0L) {
    u32bit success = 0;

    success += addSequence(nodes, nodesLen, seqptr, seqptrLen, S, false);

    reverseComplementSequence(S->sequence(), S->sequenceLength());
    success += addSequence(nodes, nodesLen, seqptr, seqptrLen, S, true);

    if (success != 2) {
      fprintf(stderr, "Failed to add sequence '%s' ('%s').\n", S->header(), S->sequence());
    }

    if (nodesLen >= nodesMax)
      fprintf(stderr, "ERROR: out of node space.\n"), exit(1);
    if (seqptrLen >= seqptrMax)
      fprintf(stderr, "ERROR: out of seqptr space.\n"), exit(1);

    delete S;
  }

  delete F;

  fprintf(stderr, "Used "u32bitFMT" trie nodes.           \n", nodesLen);

  //  Fix up sequence pointers - we could probably do this inplace
  //  with some trickery, but why?

  qsort(seqptr, seqptrLen, sizeof(trieSeqPtr), trieSeqPtrCompare);

  //  Now sorted by node iid, so run through both arrays and set
  //  pointers.  We point to the first thing found, and remember
  //  the number of things found.

  for (u32bit i=0; i<seqptrLen; i++) {
    u32bit  ni = seqptr[i].nodeiid;

    if (nodes[ni].seqptr == ~u32bitZERO)
      nodes[ni].seqptr = i;

    nodes[ni].numseq++;
  }

  //

  F = new FastAWrapper(genome);
  S = 0L;
  while ((S = F->getSequence()) != 0L) {
    char    *s    = S->sequence();
    u32bit   siid = S->getIID();
    u32bit   spos = 0;

    u32bit   n[256] = {0};  //  Pointer into the trie
    u32bit   d[256] = {0};  //  Depth this pointer is at (== sequence length)
    u32bit   nLen = 0;

    while (*s) {
      if (validSymbol[(int)*s] == 0) {

        //  Not a valid symbol, all node pointers are killed, no exact matches
        //  possible!
        nLen = 0;

      } else {

        //  Valid symbol.  Advance all pointers, print out any
        //  matches, kill any pointers, and then finally add a new
        //  one.

        u32bit  v = compressSymbol[(int)*s];
        u32bit  ni;
        u32bit  nj;

        //  Advance pointers.
        //
        for (ni=0; ni<nLen; ni++) {
          n[ni] = nodes[n[ni]].next[v];
          d[ni]++;
          //fprintf(stderr, "pointer "u32bitFMT" moves to "u32bitFMT"\n", ni, n[ni]);
        }

        //  Kill any thing that is now dead - copy nj into ni
        //
        for (ni=0, nj=0; nj<nLen; nj++)
          if (n[nj] != ~u32bitZERO) {
            if ((ni != nj)) {
              n[ni] = n[nj];
              d[ni] = d[nj];
              //fprintf(stderr, "pointer "u32bitFMT" replaces dead pointer "u32bitFMT"\n", nj, ni);
            }
            ni++;
          }
        nLen = ni;

        //  Print any matches
        //
        for (ni=0; ni<nLen; ni++) {
          if (nodes[n[ni]].numseq > 0) {
            //fprintf(stderr, "nLen="u32bitFMT"\n", nLen);

            for (nj=0; nj<nodes[n[ni]].numseq; nj++) {
              u32bit  p = nodes[n[ni]].seqptr + nj;

#if 0
              fprintf(stdout, "match ni="u32bitFMT" n[ni]="u32bitFMT" iid="u32bitFMT" gen="u32bitFMT" pos="u32bitFMT"\n",
                      ni,
                      n[ni],
                      seqptr[p].seqiid,
                      siid,
                      spos - d[ni]);
#endif

              fprintf(stdout, "sim4begin\n");
              fprintf(stdout, u32bitFMT"["u32bitFMT"-0-0] "u32bitFMT"[0-0] <"u32bitFMT"-0-100-%s-unknown>\n",
                      seqptr[p].seqiid,
                      d[ni] + 1,
                      siid,
                      d[ni] + 1,
                      seqptr[p].reversed ? "complement" : "forward");
              fprintf(stdout, "edef=%s\n", seqptr[p].defline);
              fprintf(stdout, "ddef=%s\n", S->header());
              fprintf(stdout, "1-"u32bitFMT" ("u32bitFMT"-"u32bitFMT") <"u32bitFMT"-0-100>\n",
                      d[ni] + 1,
                      spos - d[ni] + 1,
                      spos + 1,
                      d[ni]  + 1);
              fprintf(stdout, "sim4end\n");
            }
          }
        }

        //  Add a new pointer for the just seen letter
        //
        if (nodes[0].next[v] != ~u32bitZERO) {
          d[nLen]   = 0;
          n[nLen++] = nodes[0].next[v];
          //fprintf(stderr, "add new pointer "u32bitFMT" pointing to "u32bitFMT"\n", nLen-1, n[nLen-1]);
        }
      }

      s++;
      spos++;
    }

    delete S;
  }

  delete F;
}
