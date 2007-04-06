#include "util++.H"
#include "bio++.H"

//#define ALPHALEN   4
#define ALPHALEN   20

//  NOTE that our list of letters is not alphabetic.  The DNA letters
//  are first, then the rest of the protein letters.
//
const char  *trieAlpha         = "acgtdefhiklmnpqrsvwy";
u32bit       trieAlphaMap[256] = {0};

class trieNode {
public:
  trieNode(void) {
    for (u32bit i=0; i<ALPHALEN; i++)
      next[i] = ~u32bitZERO;
    numseq  = 0;
    seqptr  = ~u32bitZERO;
  };

  u32bit   next[ALPHALEN];  //  next node for A, C, G, T input, ~u32bitZERO if no next
  u32bit   numseq;          //  number of seqs we have in here
  u32bit   seqptr;         //  pointer to seqs
};


class trieSeqPtr {
public:
  trieSeqPtr(void) {
    seqiid     = ~u32bitZERO;
    nodeiid    = ~u32bitZERO;
    defline    = 0L;
    reversed   = false;
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
            seqInCore *S,
            bool isReverse) {
  char   *s = 0L;
  u32bit  n = 0;

  if (S->sequenceLength() < 12)
    return(0);

  for (s = S->sequence(); *s; s++)
    if (trieAlphaMap[*s] == 0)
      return(0);

  for (s = S->sequence(); *s; s++) {
    u32bit  v = trieAlphaMap[*s] - 1;

    //  add a new pointer if needed
    if (nodes[n].next[v] == ~u32bitZERO)
      nodes[n].next[v] = nodesLen++;

    //  Go there
    n = nodes[n].next[v];
  }

  //  add this sequence to node i -- after all sequences have been
  //  added, we'll sort this list and build pointers.

  seqptr[seqptrLen].seqiid     = S->getIID();
  seqptr[seqptrLen].nodeiid    = n;
  seqptr[seqptrLen].defline    = strdup(S->header());
  seqptr[seqptrLen].reversed   = isReverse;
  seqptrLen++;

  return(1);
}



int
main(int argc, char **argv) {
  char  *queries = 0L;
  char  *genome  = 0L;
  FILE  *logfile = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-q") == 0) {
      queries = argv[++arg];
    } else if (strcmp(argv[arg], "-g") == 0) {
      genome = argv[++arg];
    } else if (strcmp(argv[arg], "-l") == 0) {
      errno = 0;
      logfile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open logfile '%s': %s\n", argv[arg], strerror(errno)), exit(1);
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
    fprintf(stderr, "  -q queries.fasta   -- the input with the short stuff\n");
    fprintf(stderr, "  -q queries.fasta   -- the input with the short stuff\n");
    exit(1);
  }

  for (u32bit i=0; i<ALPHALEN; i++) {
    trieAlphaMap[trieAlpha[i]] = i+1;
    trieAlphaMap[trieAlpha[i] + 'A' - 'a'] = i+1;
  }


  //  Build the trie from the queries.  We just allocate a large
  //  number of nodes so we don't need to deal with reallocation.
  //
  //  32M * 28 = 896M.
  //
  u32bit     nodesLen = 1;
  u32bit     nodesMax = 32 * 1024 * 1024;
  trieNode  *nodes    = new trieNode [nodesMax];

  seqFile              *F = openSeqFile(queries);
  seqInCore            *S = 0L;

  u32bit      seqptrLen = 0;
  u32bit      seqptrMax = 2 * 1024 * 1024;
  trieSeqPtr *seqptr    = new trieSeqPtr [seqptrMax];

  //  Number of matches per IID, not seqptr (which has two entries for
  //  each iid, one forward, one reverse).
  //
  u32bit     *nummatches = new u32bit [seqptrMax];
  for (u32bit i=0; i<seqptrMax; i++)
    nummatches[i] = 0;

  while ((S = F->getSequenceInCore()) != 0L) {
    u32bit success = 0;

    success += addSequence(nodes, nodesLen, seqptr, seqptrLen, S, false);

#if ALPHALEN == 4
    reverseComplementSequence(S->sequence(), S->sequenceLength());
    success += addSequence(nodes, nodesLen, seqptr, seqptrLen, S, true);
#else
    success++;
#endif

    if (success != 2)
      if (logfile)
        fprintf(logfile, "Failed to add sequence '%s' ('%s').\n", S->header(), S->sequence());

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

  F = openSeqFile(genome);
  S = 0L;
  while ((S = F->getSequenceInCore()) != 0L) {
    char    *s    = S->sequence();
    u32bit   siid = S->getIID();
    u32bit   spos = 0;

    u32bit   n[256] = {0};  //  Pointer into the trie
    u32bit   d[256] = {0};  //  Depth this pointer is at (== sequence length)
    u32bit   nLen = 0;

    //fprintf(stderr, "WORKING ON '%s'\n", S->header());

    while (*s) {
      if (trieAlphaMap[*s] == 0) {

        //  Not a valid symbol, all node pointers are killed, no exact matches
        //  possible!
        nLen = 0;

      } else {

        //  Valid symbol.  Advance all pointers, print out any
        //  matches, kill any pointers, and then finally add a new
        //  one.

        u32bit  v = trieAlphaMap[*s] - 1;
        u32bit  ni;
        u32bit  nj;

        //  Advance pointers.
        //
        for (ni=0; ni<nLen; ni++) {
          n[ni] = nodes[n[ni]].next[v];
          d[ni]++;
        }

        //  Kill any thing that is now dead - copy nj into ni
        //
        for (ni=0, nj=0; nj<nLen; nj++)
          if (n[nj] != ~u32bitZERO) {
            if ((ni != nj)) {
              n[ni] = n[nj];
              d[ni] = d[nj];
            }
            ni++;
          }
        nLen = ni;

        //  Print any matches
        //
        for (ni=0; ni<nLen; ni++) {
          if (nodes[n[ni]].numseq > 0) {
            for (nj=0; nj<nodes[n[ni]].numseq; nj++) {
              u32bit  p = nodes[n[ni]].seqptr + nj;

              nummatches[seqptr[p].seqiid]++;
              if (nummatches[seqptr[p].seqiid] == 1000) {
                if (logfile)
                  fprintf(logfile, "sequence "u32bitFMT" '%s' has too many matches, not reporting any more.\n",
                          seqptr[p].seqiid,
                          seqptr[p].defline);
              } else if (nummatches[seqptr[p].seqiid] < 1000) {
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
        }

        //  Add a new pointer for the just seen letter
        //
        if (nodes[0].next[v] != ~u32bitZERO) {
          d[nLen]   = 0;
          n[nLen++] = nodes[0].next[v];
        }
      }

      s++;
      spos++;
    }

    delete S;
  }

  //  We should print out the total number of matches for each
  //  sequence....  Report those with matches first.
  //
  if (logfile) {
    for (u32bit i=0; i<seqptrLen; i++)
      if ((seqptr[i].reversed == false) && (nummatches[seqptr[i].seqiid] > 0))
        fprintf(logfile, "sequence "u32bitFMT" '%s' has "u32bitFMT" matches.\n",
                seqptr[i].seqiid,
                seqptr[i].defline,
                nummatches[seqptr[i].seqiid]);

    for (u32bit i=0; i<seqptrLen; i++)
      if ((seqptr[i].reversed == false) && (nummatches[seqptr[i].seqiid] == 0))
        fprintf(logfile, "sequence "u32bitFMT" '%s' has no matches.\n",
                seqptr[i].seqiid,
                seqptr[i].defline,
                nummatches[seqptr[i].seqiid]);
  }

  delete F;
}
