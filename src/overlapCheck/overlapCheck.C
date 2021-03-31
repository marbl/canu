
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"

#include "files.H"
#include "kmers.H"

#include "sqStore.H"
#include "sqCache.H"

#include "ovStore.H"

#include "edlib.H"

#include <vector>
#include <algorithm>



bool
overlapBefore(ovOverlap &A, ovOverlap &B) {
  if ((A.a_iid     < B.a_iid))
    return(true);

  if ((A.a_iid    == B.a_iid) &&
      (A.b_iid     < B.b_iid))
    return(true);

  if ((A.a_iid    == B.a_iid) &&
      (A.b_iid    == B.b_iid) &&
      (A.flipped() < B.flipped()))
    return(true);

  return(false);
}



void
findReadKmers(char const      *seq,
              uint32           seqLen,
              uint32           olapBgn,
              uint32           olapEnd,
              std::set<kmer>  &inOverlap,
              std::set<kmer>  &allKmer) {

  kmerIterator  it(seq, seqLen);

  inOverlap.clear();
  allKmer.clear();

  if (olapEnd < olapBgn)
    std::swap(olapBgn, olapEnd);

  while (it.nextMer() == true) {
    kmer f = it.fmer();
    kmer r = it.rmer();
    kmer c = (f < r) ? f : r;

    allKmer.insert(c);

    if ((olapBgn < it.bgnPosition()) &&
        (it.endPosition() <= olapEnd))
      inOverlap.insert(c);
  }
}



void
evaluateKmers(ovOverlap ovl,
              const char *Aalign, uint32 AalignLen,
              const char *Balign, uint32 BalignLen,
              uint32 &sharedRepeat,
              uint32 &sharedUnique,
              merylExactLookup *lookup) {

  kmerIterator  Ait(Aalign, AalignLen);
  kmerIterator  Bit(Balign, BalignLen);

  while ((Ait.nextBase() == true) &&
         (Bit.nextBase() == true)) {
    kmer af = Ait.fmer(), ar = Ait.rmer();
    kmer bf = Bit.fmer(), br = Bit.rmer();

    kmer ac = (af < ar) ? af : ar;
    kmer bc = (bf < br) ? bf : br;

    kmvalu av = lookup->value(ac);
    kmvalu bv = lookup->value(bc);


    if ((Ait.isValid() == false) && (Bit.isValid() == false)) {
      continue;
    }

    if ((Ait.isValid() == true) && (Bit.isValid() == false)) {
      continue;
    }
    if ((Ait.isValid() == false) && (Bit.isValid() == true)) {
      continue;
    }

    //  Both valid.

    if (ac != bc) {
      continue;
    }

    //  Both the same.

    assert(av == bv);

    if (av <= 1333) {
      sharedUnique++;
    }

    else {
      sharedRepeat++;
    }
  }
}




void
dumpKmers(ovOverlap   ovl,
          char const *Aalign,  uint32 AalignLen,  uint32 abgn, uint32 aend,
          char const *Balign,  uint32 BalignLen,  uint32 bbgn, uint32 bend,
          merylExactLookup *lookup) {

  char  NA[FILENAME_MAX];
  char  NK[FILENAME_MAX];

  sprintf(NA, "align/%08u-%08u-%s.align", ovl.a_iid, ovl.b_iid, ovl.flipped() ? "flip" : "norm");
  sprintf(NK, "align/%08u-%08u-%s.kmers", ovl.a_iid, ovl.b_iid, ovl.flipped() ? "flip" : "norm");

  FILE *FA = AS_UTL_openOutputFile(NA);
  FILE *FK = AS_UTL_openOutputFile(NK);

  //  Output the alignments.

  fprintf(FA, "<pre>\n");
  fprintf(FA, "%-8u %6u-%-6u %s %s\n", ovl.a_iid, abgn, aend, ovl.flipped() ? "fwd" : "fwd", Aalign);
  fprintf(FA, "%-8u %6u-%-6u %s %s\n", ovl.b_iid, bbgn, bend, ovl.flipped() ? "rev" : "fwd", Balign);
  fprintf(FA, "</pre>\n");

  //  Output the kmers.

  kmerIterator  Ait(Aalign, AalignLen);
  kmerIterator  Bit(Balign, BalignLen);

  while ((Ait.nextBase() == true) &&
         (Bit.nextBase() == true)) {
    kmer af = Ait.fmer(), ar = Ait.rmer();
    kmer bf = Bit.fmer(), br = Bit.rmer();

    kmer ac = (af < ar) ? af : ar;
    kmer bc = (bf < br) ? bf : br;

    if      (Ait.isValid() &&
             Bit.isValid()) {
      uint64  pos = Ait.bgnPosition();

      fprintf(FK, "%lu\t%u\t%u\t%c\t%c\n",
              pos,
              lookup->value(ac), lookup->value(bc),
              Aalign[pos], Balign[pos]);
    }

    else if (Ait.isValid()) {
      //fprintf(FK, "%u\t%u\t%u\n", Ait.bgnPosition(), lookup->value(ac), 0);
    }

    else if (Bit.isValid()) {
      //fprintf(FK, "%u\t%u\t%u\n", Ait.bgnPosition(), 0, lookup->value(bc));
    }

    else {
    }
  }

  //  Cleanup.

  AS_UTL_closeFile(FA);
  AS_UTL_closeFile(FK);
}



void
compareOverlapKmers(sqCache *cache,
                    merylExactLookup *lookup,
                    ovOverlap ovl) {

  std::set<kmer>  Ain, Aall;
  std::set<kmer>  Bin, Ball;

  uint32  AbasesLen = 0,  AbasesMax = 0;    char *Abases = nullptr;
  uint32  BbasesLen = 0,  BbasesMax = 0;    char *Bbases = nullptr;

  uint32  AalignLen = 0;                    char *Aalign = nullptr;
  uint32  BalignLen = 0;                    char *Balign = nullptr;

  cache->sqCache_getSequence(ovl.a_iid, Abases, AbasesLen, AbasesMax);
  cache->sqCache_getSequence(ovl.b_iid, Bbases, BbasesLen, BbasesMax);

  if (ovl.flipped() == true)
    reverseComplementSequence(Bbases, BbasesLen);

  EdlibAlignResult  res;

  int32  abgn = ovl.a_bgn();
  int32  aend = ovl.a_end();
  int32  bbgn = ovl.flipped() ? BbasesLen - ovl.b_bgn() : ovl.b_bgn();
  int32  bend = ovl.flipped() ? BbasesLen - ovl.b_end() : ovl.b_end();

  //fprintf(stderr, "Align %u %u-%u %u  to  %u %u-%u %u\n",
  //        ovl.a_iid, abgn, aend, AbasesLen,
  //        ovl.b_iid, bbgn, bend, BbasesLen);

  res = edlibAlign(Abases + abgn, aend - abgn,
                   Bbases + bbgn, bend - bbgn,
                   edlibNewAlignConfig((aend - abgn) * 0.25, EDLIB_MODE_NW, EDLIB_TASK_PATH));

  AalignLen = res.alignmentLength + 1;    Aalign = new char [AalignLen];
  BalignLen = res.alignmentLength + 1;    Balign = new char [BalignLen];

  if (res.numLocations == 0) {
    fprintf(stderr, "Not aligned.\n");
    assert(0);
  } else {
    //fprintf(stderr, "Aligned with distance %d length %d\n", res.editDistance, res.alignmentLength);
  }

  edlibAlignmentToStrings(res,
                          Abases+abgn, aend-abgn,
                          Bbases+bbgn, bend-bbgn,
                          Aalign,
                          Balign);




  uint32 sharedRepeat = 0;
  uint32 sharedUnique = 0;

  evaluateKmers(ovl,
                Aalign, AalignLen,
                Balign, BalignLen,
                sharedRepeat,
                sharedUnique,
                lookup);


#if 0
  findReadKmers(Aalign, AalignLen, 0, AalignLen, Ain, Aall);
  findReadKmers(Balign, BalignLen, 0, BalignLen, Bin, Ball);

  for (auto it=Ain.begin(); it != Ain.end(); ++it) {
    kmer    k = *it;
    kmvalu  v = lookup->value(k);

    if (Bin.count(k) == 0)
      continue;

    if (v == 0) {
      nocount++;
      continue;
    }

    if (v < 1333) {
      sharedUnique++;
    } else {
      sharedRepeat++;
    }
  }
#endif

  if (sharedUnique > 0) {
    fprintf(stdout, "KMERS unique  %u\n", sharedUnique);
    fprintf(stdout, "      repeat  %u\n", sharedRepeat);
    fprintf(stdout, "OVERLAP A %6u-%6u %6lu/%6lu kmers\n", ovl.a_bgn(), ovl.a_end(), Ain.size(), Aall.size());
    fprintf(stdout, "        B %6u-%6u %6lu/%6lu kmers\n", ovl.b_bgn(), ovl.b_end(), Bin.size(), Ball.size());
    fprintf(stdout, "perl plot-kmer.pl %u %u %u %u  %u %u %u %u\n",
            ovl.a_iid, ovl.a_bgn(), ovl.a_end(), AbasesLen,
            ovl.b_iid, ovl.b_bgn(), ovl.b_end(), BbasesLen);
    fprintf(stdout, "plot [][1:6666] 'align/%08d-%08d-%s.kmers' using 1:2, 'align/%08d-%08d-%s.kmers' using 1:3, 1333\n",
            ovl.a_iid, ovl.b_iid, ovl.flipped() ? "flip" : "norm",
            ovl.a_iid, ovl.b_iid, ovl.flipped() ? "flip" : "norm");
    fprintf(stdout, "\n");

    dumpKmers(ovl,
              Aalign, AalignLen, abgn, aend,
              Balign, BalignLen, bbgn, bend, lookup);
  }

  delete [] Abases;   delete [] Aalign;
  delete [] Bbases;   delete [] Balign;
}



void
compareFiles(sqStore *seq,
             sqCache *cache,
             ovFile *ovfA, ovFile *ovfB,
             merylExactLookup *lookup) {
  uint64      ovlAmax = ovfA->getCounts()->numOverlaps();
  uint64      ovlBmax = ovfB->getCounts()->numOverlaps();

  //  Allocate space for and load overlaps.

  fprintf(stderr, "Allocating 2x %lu A and 2x %lu B overlaps.\n", ovlAmax, ovlBmax);

  ovOverlap  *ovlA    = new ovOverlap [2 * ovlAmax];
  ovOverlap  *ovlB    = new ovOverlap [2 * ovlBmax];

  fprintf(stderr, "Loading overlaps.\n");

  uint64      ovlAlen = ovfA->readOverlaps(ovlA, ovlAmax);
  uint64      ovlBlen = ovfB->readOverlaps(ovlB, ovlBmax);

  //  Duplicate all overlaps into their other twin

  fprintf(stderr, "Duplicating overlaps.\n");

  for (uint32 ii=0; ii<ovlAmax; ii++, ovlAlen++)
    ovlA[ovlAlen].swapIDs(ovlA[ii]);

  for (uint32 ii=0; ii<ovlBmax; ii++, ovlBlen++)
    ovlB[ovlBlen].swapIDs(ovlB[ii]);

  //  Sort.

  fprintf(stderr, "Sorting overlaps.\n");

  std::sort(ovlA, ovlA + ovlAlen);
  std::sort(ovlB, ovlB + ovlBlen);

  uint64  onlyA=0, onlyB=0, same=0;

  fprintf(stderr, "Comparing overlaps.\n");

  for (uint32 aa=0, bb=0; (aa < ovlAlen) || (bb < ovlBlen); ) {

    //  If one runs out of overlaps, examine the other one.

    if      (aa == ovlAlen) {
      fprintf(stdout, "                                     -- B#%-5u %8u-%8u %c %6.3f%%\n",
              bb, ovlB[bb].a_iid, ovlB[bb].b_iid, ovlB[bb].flipped() ? 'f' : 'n', ovlB[bb].erate() * 100.0);

      onlyB++;
      bb++;
      continue;
    }

    if (bb == ovlBlen) {
      fprintf(stdout, "A#%-5u %8u-%8u %c %6.3f%% --\n",
              aa, ovlA[aa].a_iid, ovlA[aa].b_iid, ovlA[aa].flipped() ? 'f' : 'n', ovlA[aa].erate() * 100.0);

      onlyA++;
      aa++;
      continue;
    }

    //  If they're the same reads and the same orientation, count
    //  it as a match.
    if ((ovlA[aa].a_iid     == ovlB[bb].a_iid) &&
        (ovlA[aa].b_iid     == ovlB[bb].b_iid) &&
        (ovlA[aa].flipped() == ovlB[bb].flipped())) {
      fprintf(stdout, "A#%-5u %8u-%8u %c %6.3f%% == B#%-5u %8u-%8u %c %6.3f%%\n",
              aa, ovlA[aa].a_iid, ovlA[aa].b_iid, ovlA[aa].flipped() ? 'f' : 'n', ovlA[aa].erate() * 100.0,
              bb, ovlB[bb].a_iid, ovlB[bb].b_iid, ovlB[bb].flipped() ? 'f' : 'n', ovlB[bb].erate() * 100.0);

      same++;
      aa++;
      bb++;
      continue;
    }

    //  

    if (overlapBefore(ovlA[aa], ovlB[bb]) == true) {
      fprintf(stdout, "A#%-5u %8u-%8u %c %6.3f%% <+\n",
              aa, ovlA[aa].a_iid, ovlA[aa].b_iid, ovlA[aa].flipped() ? 'f' : 'n', ovlA[aa].erate() * 100.0);

      onlyA++;
      aa++;
    } else {
      fprintf(stdout, "                                    +> B#%-5u %8u-%8u %c %6.3f%%\n",
              bb, ovlB[bb].a_iid, ovlB[bb].b_iid, ovlB[bb].flipped() ? 'f' : 'n', ovlB[bb].erate() * 100.0);

      compareOverlapKmers(cache, lookup, ovlB[bb]);

      onlyB++;
      bb++;
    }
  }

  fprintf(stderr, "onlyA:  %lu\n", onlyA);
  fprintf(stderr, "onlyB:  %lu\n", onlyB);
  fprintf(stderr, "same:   %lu\n", same);

  delete [] ovlA;
  delete [] ovlB;
}



void
compareStores(sqStore *seq,
              ovStore *ovfA, ovStore *ovfB,
              merylExactLookup *lookup) {
}



int
main(int argc, char **argv) {
  char const   *seqStoreName  = nullptr;
  char const   *ovlFileNameA  = nullptr;
  char const   *ovlStoreNameA = nullptr;
  char const   *ovlFileNameB  = nullptr;
  char const   *ovlStoreNameB = nullptr;
  char const   *kmerDBname    = nullptr;
  uint64        maxMemory     = 32llu * 1024 * 1024 * 1024;

  argc = AS_configure(argc, argv);

  omp_set_num_threads(16);

  std::vector<char const *>  err;
  for (int arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-S") == 0) {
      seqStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-Af") == 0) {
      ovlFileNameA = argv[++arg];
    } else if (strcmp(argv[arg], "-A") == 0) {
      ovlStoreNameA = argv[++arg];

    } else if (strcmp(argv[arg], "-Bf") == 0) {
      ovlFileNameB = argv[++arg];
    } else if (strcmp(argv[arg], "-B") == 0) {
      ovlStoreNameB = argv[++arg];

    } else if (strcmp(argv[arg], "-K") == 0) {
      kmerDBname = argv[++arg];

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (seqStoreName == nullptr)   err.push_back("No sequence store (-S) supplied.\n");
  if (ovlFileNameA == nullptr)   err.push_back("No overlap file A (-A) supplied.\n");
  if (ovlFileNameB == nullptr)   err.push_back("No overlap file B (-B) supplied.\n");
  if (kmerDBname   == nullptr)   err.push_back("No kmer database (-S) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore -A ovlA -B ovlB -K kmers\n", argv[0]);
    fprintf(stderr, "\n");
    return(1);
  }

  //  Open inputs.

  sqStore  *seq   = new sqStore(seqStoreName);
  sqCache  *cache = new sqCache(seq);

  fprintf(stderr, "Loading reads.\n");
  cache->sqCache_loadReads(true);

  //  Make a kmer lookup table.

  merylFileReader   *reader = new merylFileReader(kmerDBname);
  merylExactLookup  *lookup = new merylExactLookup();

  fprintf(stderr, "Loading kmers.\n");
  lookup->load(reader, maxMemory, true, false, uint32min, uint32max);

  delete reader;

  fprintf(stderr, "Lookup table loaded.\n");

  if (ovlFileNameA && ovlFileNameB) {
    ovFile   *ovfA = new ovFile(seq, ovlFileNameA, ovFileFull);
    ovFile   *ovfB = new ovFile(seq, ovlFileNameB, ovFileFull);

    compareFiles(seq, cache, ovfA, ovfB, lookup);

    delete ovfA;
    delete ovfB;
  }


  if (ovlStoreNameA && ovlStoreNameB) {
    ovStore  *ovsA = new ovStore(ovlStoreNameA, nullptr);
    ovStore  *ovsB = new ovStore(ovlStoreNameB, nullptr);

    compareStores(seq, ovsA, ovsB, lookup);

    delete ovsA;
    delete ovsB;
  }


  //  Done.  Cleanup.

  delete lookup;
  delete seq;

  return(0);
}
