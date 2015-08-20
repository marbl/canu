
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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


//  NB: My hacked kazlib returns a pointer to whatever we give it.
//  Since we gave it a pointer to an object, it gives us back a
//  pointer to "a pointer to an object".  Hence, this ugliness.
//
int
kMerLiteSort(void const *a, void const *b) {
  kMerLite const *A = *((kMerLite * const *)a);
  kMerLite const *B = *((kMerLite * const *)b);

  if (*A < *B) return(-1);
  if (*A > *B) return(1);
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

  fprintf(stderr, "Mers are "uint32FMT" bases.\n", MSR->merSize());
  fprintf(stderr, "There are "uint64FMT" unique (copy = 1) mers.\n", MSR->numberOfUniqueMers());
  fprintf(stderr, "There are "uint64FMT" distinct mers.\n", MSR->numberOfDistinctMers());
  fprintf(stderr, "There are "uint64FMT" mers total.\n", MSR->numberOfTotalMers());

  //  Guess how many mers we can fit into 512MB, then report how many chunks we need to do.

  uint32  merSize      = MSR->merSize();
  uint64  memoryLimit  = 700 * 1024 * 1024;
  uint64  perMer       = sizeof(kMerLite) + sizeof(dnode_t);
  uint64  mersPerBatch = memoryLimit / perMer;
  uint32  numBatches   = MSR->numberOfDistinctMers() / mersPerBatch;
  uint32  batch        = 0;

  dnode_t   *nodes     = new dnode_t  [mersPerBatch];
  kMerLite  *mers      = new kMerLite [mersPerBatch];

  if (MSR->numberOfDistinctMers() % mersPerBatch)
    numBatches++;

  fprintf(stderr, "perMer:  "uint64FMT" bytes ("uint64FMT" for kMerLite, "uint64FMT" for dnode_t.\n",
          perMer, (uint64)sizeof(kMerLite), (uint64)sizeof(dnode_t));
  fprintf(stderr, "We can fit "uint64FMT" mers into "uint64FMT"MB.\n", mersPerBatch, memoryLimit >> 20);
  fprintf(stderr, "So we need "uint32FMT" batches to verify the count.\n", numBatches);

  while (MSR->validMer()) {
    uint64          mersRemain = mersPerBatch;
    dict_t         *merDict    = dict_create(mersPerBatch, kMerLiteSort);

    batch++;

    //  STEP 1:  Insert mersPerBatch into the merDict
    //
    fprintf(stderr, "STEP 1 BATCH "uint32FMTW(2)":  Insert into merDict\n", batch);
    while (MSR->nextMer() && mersRemain) {
      mersRemain--;

      mers[mersRemain] = MSR->theFMer();

      //  initialize the node with the value, then insert the node
      //  into the tree using the key

      int32 val = (int32)MSR->theCount();
      dnode_init(&nodes[mersRemain], (void *)val);
      dict_insert(merDict, &nodes[mersRemain], &mers[mersRemain]);
    }

    //  STEP 2:  Stream the original file, decrementing the count
    //
    fprintf(stderr, "STEP 2 BATCH "uint32FMTW(2)":  Stream fasta\n", batch);
    seqStream    *CS = new seqStream(fastaName, true);
    merStream    *MS = new merStream(new kMerBuilder(merSize), CS);

    kMerLite       mer;
    dnode_t       *nod;

    while (MS->nextMer()) {
      mer = MS->theFMer();

      nod = dict_lookup(merDict, &mer);

      if (nod != 0L) {
        int32 val = (int32)dnode_get(nod);
        val--;
        dnode_put(nod, (void *)val);
      } else {
        //  Unless the whole meryl file fit into our merDict, we cannot warn if
        //  we don't find mers.
        //
        if (numBatches == 1) {
          char str[1024];
          fprintf(stderr, "Didn't find node for mer '%s'\n", mer.merToString(merSize, str));
        }
      }
    }

    delete MS;
    delete CS;

    //  STEP 3:  Check every node in the tree to make sure that the counts
    //  are exactly zero.
    //
    fprintf(stderr, "STEP 3 BATCH "uint32FMTW(2)":  Check\n", batch);
    nod = dict_first(merDict);
    while (nod) {
      int32           val = (int32)dnode_get(nod);
      kMerLite const  *nodmer = (kMerLite const *)dnode_getkey(nod);

      if (val != 0) {
        char str[1024];
        fprintf(stderr, "Got count "int32FMT" for mer '%s'\n",
                val,
                nodmer->merToString(merSize, str));
      }

      nod = dict_next(merDict, nod);
    }


    //  STEP 4:  Destroy the dictionary.
    //
    fprintf(stderr, "STEP 4 BATCH "uint32FMTW(2)":  Destroy\n", batch);
    while ((nod = dict_first(merDict)))
      dict_delete(merDict, nod);
    dict_destroy(merDict);
  }
}
