#include "bio++.H"
#include "libmeryl.H"
#include "kmerlite.H"

//  This tests that all the mers in an input fasta file are counted
//  properly.  It does not test that the meryl output contains exactly
//  those mers, just that those mers are there.
//
//  If you can fit into one batch, then it _will_ verift that the
//  meryl output is exactly correct.
//
//  Reads a meryl-format kmer count in chunks.  Each chunk is stored
//  in a searchable structure (we should be using, say, an extended
//  existDB, but we're using a balanced binary tree).  The entire
//  source fasta file is then streamed against the kmer chunk,
//  decrementing the count for each mer.  When the whole file is
//  streamed, any kmers with positive count are reported.


class kMerLiteCount {
public:
  kMerLite  K;
  s32bit    C;
};


//  NB: My hacked kazlib returns a pointer to whatever we give it.
//  Since we gave it a pointer to an object, it gives us back a
//  pointer to "a pointer to an object".  Hence, this ugliness.
//
int
kMerLiteCountSort(void const *a, void const *b) {
  kMerLiteCount const *A = *((kMerLiteCount const **)a);
  kMerLiteCount const *B = *((kMerLiteCount const **)b);

  if (A->K < B->K) return(-1);
  if (A->K > B->K) return(1);
  return(0);
}



int
main(int argc, char **argv) {

  char   *merylCount = 0L;
  char   *fastaName = 0L;

  int arg=1;
  while (arg < argc) {

    if        (strcmp(argv[arg], "-m") == 0) {
      merylCount = argv[++arg];
    } else if (strcmp(argv[arg], "-f") == 0) {
      fastaName = argv[++arg];
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((merylCount == 0L) || (fastaName == 0L)) {
    fprintf(stderr, "usage: %s -m <meryl-name-prefix> -f <fasta-file>\n", argv[0]);
    exit(1);
  }


  //  Open the count files
  //
  merylStreamReader  *MSR = new merylStreamReader(merylCount);

  fprintf(stderr, "Mers are "u32bitFMT" bases.\n", MSR->merSize());
  fprintf(stderr, "There are "u64bitFMT" unique (copy = 1) mers.\n", MSR->numberOfUniqueMers());
  fprintf(stderr, "There are "u64bitFMT" distinct mers.\n", MSR->numberOfDistinctMers());
  fprintf(stderr, "There are "u64bitFMT" mers total.\n", MSR->numberOfTotalMers());

  //  Guess how many mers we can fit into 512MB, then report how many chunks we need to do.

  u32bit  merSize      = MSR->merSize();
  u64bit  memoryLimit  = 700 * 1024 * 1024;
  u64bit  perMer       = sizeof(kMerLite) + sizeof(dnode_t);
  u64bit  mersPerBatch = memoryLimit / perMer;
  u32bit  numBatches   = MSR->numberOfDistinctMers() / mersPerBatch;
  u32bit  batch        = 0;

  dnode_t        *nodes      = new dnode_t       [mersPerBatch];
  kMerLiteCount  *mers       = new kMerLiteCount [mersPerBatch];

  if (MSR->numberOfDistinctMers() % mersPerBatch)
    numBatches++;

  fprintf(stderr, "perMer:  "u64bitFMT" bytes ("u64bitFMT" for kMerLite, "u64bitFMT" for dnode_t.\n",
          perMer, (u64bit)sizeof(kMerLite), (u64bit)sizeof(dnode_t));
  fprintf(stderr, "We can fit "u64bitFMT" mers into "u64bitFMT"MB.\n", mersPerBatch, memoryLimit >> 20);
  fprintf(stderr, "So we need "u64bitFMT" batches to verify the count.\n", numBatches);

  while (MSR->validMer()) {
    u64bit          mersRemain = mersPerBatch;
    dict_t         *merDict    = dict_create(mersPerBatch, kMerLiteCountSort);

    batch++;

    //  STEP 1:  Insert mersPerBatch into the merDict
    //
    fprintf(stderr, "STEP 1 BATCH "u32bitFMTW(2)":  Insert into merDict\n", batch);
    while (MSR->nextMer() && mersRemain) {
      mersRemain--;

      mers[mersRemain].K = MSR->theFMer();
      mers[mersRemain].C = MSR->theCount();

#if 0
      char str1[1024];
      char str2[1024];
      fprintf(stderr, "insert '%s' ->\n       '%s' (%p)\n",
              MSR->theFMer().merToString(str1),
              mers[mersRemain].K.merToString(merSize, str2), &mers[mersRemain]);
#endif

      //  initialize the node with the value, then insert the node
      //  into the tree using the key
      //
      //  XXX: The node's value (dnode_init) is the same as the key
      //  (dict_insert), we could use the value to store the count.
      //
      dnode_init(&nodes[mersRemain], &mers[mersRemain]);
      dict_insert(merDict, &nodes[mersRemain], &mers[mersRemain]);
    }

    //  STEP 2:  Stream the original file, decrementing the count
    //
    fprintf(stderr, "STEP 2 BATCH "u32bitFMTW(2)":  Stream fasta\n", batch);
    FastAstream  *FS = new FastAstream(fastaName);
    merStream    *MS = new merStream(merSize, FS);

    kMerLiteCount  mer;
    dnode_t       *nod;
    kMerLiteCount *nodmer;

    while (MS->nextMer()) {
      mer.K = MS->theFMer();
      mer.C = 0;

      nod = dict_lookup(merDict, &mer);

      if (nod != 0L) {
        nodmer = (kMerLiteCount *)dnode_get(nod);
        nodmer->C--;
      } else {
        //  Unless the whole meryl file fit into our merDict, we cannot warn if
        //  we don't find mers.
        //
        if (numBatches == 1) {
          char str[1024];
          fprintf(stderr, "Didn't find node for mer '%s'\n", mer.K.merToString(merSize, str));
        }
      }
    }

    delete MS;
    delete FS;

    //  STEP 3:  Check every node in the tree to make sure that the counts
    //  are exactly zero.
    //
    fprintf(stderr, "STEP 3 BATCH "u32bitFMTW(2)":  Check\n", batch);
    nod = dict_first(merDict);
    while (nod) {
      nodmer = (kMerLiteCount *)dnode_get(nod);

      if (nodmer->C != 0) {
        char str[1024];
        fprintf(stderr, "Got count "s32bitFMT" for mer '%s'\n", nodmer->C, nodmer->K.merToString(merSize, str));
      }

      nod = dict_next(merDict, nod);
    }


    //  STEP 4:  Destroy the dictionary.
    //
    fprintf(stderr, "STEP 4 BATCH "u32bitFMTW(2)":  Destroy\n", batch);
    while ((nod = dict_first(merDict)))
      dict_delete(merDict, nod);
    dict_destroy(merDict);
  }
}
