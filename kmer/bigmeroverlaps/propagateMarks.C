#include "util++.H"
//#include "gmp.h"   //  Include this before bio++ so that kMer knows we have GMP!
#include "bio++.H"

mt_s    *mt = 0L;

#include "randomhash-shifting.H"
//#include "randomhash-mt.H"
//#include "randomhash-mod.H"
#include "markFile.H"

#warning NEED TO USE CANONICAL MERS!

int
main(int argc, char **argv) {
  char    *inName        = 0L;
  char    *outName       = 0L;
  u64bit   seed          = time(0L);
  u32bit   iterations    = 16;
  u64bit   partitionBits = 0;
  u64bit   partitionMask = 0L;
  u64bit   partition     = 0;
  u32bit   merSize       = 100;
  u32bit   hashSize      = 33;

  bool     error         = false;

  bool     pass1         = false;
  bool     pass2         = false;
  char    *pass1mask     = 0L;

  bool     init          = false;

  bool     isFasta       = false;
  bool     isMerStream   = false;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      isFasta     = true;
      isMerStream = false;
      inName      = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      isFasta     = false;
      isMerStream = true;
      inName      = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outName = argv[++arg];

    } else if (strcmp(argv[arg], "-seed") == 0) {
      seed = strtou64bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-iters") == 0) {
      iterations = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-hashsize") == 0) {
      hashSize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-partitionbits") == 0) {
      partitionBits = strtou64bit(argv[++arg], 0L);
      partitionMask = u64bitMASK(partitionBits);
    } else if (strcmp(argv[arg], "-partition") == 0) {
      partition  = strtou64bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-mersize") == 0) {
      merSize = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-pass1") == 0) {
      pass1     = true;
      pass2     = false;
    } else if (strcmp(argv[arg], "-pass2") == 0) {
      pass1     = false;
      pass2     = true;
      pass1mask = argv[++arg];

    } else {
      fprintf(stderr, "invalid option '%s'\n", argv[arg]);
      error = true;
    }

    arg++;
  }

  if (inName == 0L)
    fprintf(stderr, "ERROR:  No input file given!\n"), error = true;

  if (outName == 0L)
    fprintf(stderr, "ERROR:  No output prefix given!\n"), error = true;

    fprintf(stderr, "checking partition ("u64bitFMT"); 0 <= p < "u64bitFMT".\n",
            partition, partitionMask+1);

  if ((partitionBits > 0) && (partition > partitionMask))
    fprintf(stderr, "ERROR: invalid partition ("u64bitFMT"); 0 <= p < "u64bitFMT".\n",
            partition, partitionMask+1), error = true;

  if (!pass1 && !pass2 && !init)
    fprintf(stderr, "ERROR: which pass?\n"), error = true;

  if (error) {
    fprintf(stderr, "usage: %s [options] -f in.fasta -o out-prefix\n", argv[0]);
    fprintf(stderr, "  -seed             seed randomness, a big number\n");
    fprintf(stderr, "  -iters            number of iterations\n");
    fprintf(stderr, "  -hashsize         size of the hash, in bits\n");
    fprintf(stderr, "  -partitionbits    hash partitons, in bits\n"); 
    fprintf(stderr, "  -partition        hash partition to compute (0 <= p < 2^P)\n");
    fprintf(stderr, "  -mersize          mersize, default 100\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -pass1\n");
    fprintf(stderr, "  -pass2 pass1file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -f   read sequence from here\n");
    fprintf(stderr, "  -o   prefix output files with this\n");
    exit(1);
  }

  //  If we are initializing, we just convert the fastafile into a
  //  merStreamFile.
  //
  if (init) {
    fprintf(stderr, "Converting '%s' into '%s.merStream'\n", inName, outName);

    merStreamFileBuilder  *B = new merStreamFileBuilder(100, inName, outName);
    B->build();
    delete B;

    exit(0);
  }

  mt = mtInit(seed);

  FastAstream          *FS   = 0L;
  merStream            *MS   = 0L;
  merStreamFileReader  *MSFR = 0L;

  if (isFasta) {
    FS   = new FastAstream(inName);
    MS   = new merStream(merSize, FS);
  } else if (isMerStream) {
    MSFR = new merStreamFileReader(inName, merSize);
    MS   = new merStream(MSFR);
  } else {
    fprintf(stderr, "ERROR!  Isn't fasta and isn't merstream?\n");
    exit(1);
  }

  speedCounter *SC = 0L;

  //  Space to store our annotation of "mer that starts here is
  //  distinct".  This doesn't need to be in core, and it isn't.
  //
  //  Supposedly, the bitPackedFile ensures that bits are all set to
  //  zero.
  //
  markFile   *distinct = new markFile(outName);
  markFile   *ignore   = 0L;

  if (pass2) {
    if (!fileExists(pass1mask)) {
      fprintf(stderr, "ERROR: -pass2 file '%s' doesn't exist!\n", pass1mask);
      exit(1);
    }
    ignore = new markFile(pass1mask);
  }

  //  Early versions used to pick a random hashSize, but why bother?
  //  Just use the largest we can fit in memory.
  //
  randomHash       H(merSize, hashSize);

  for (u32bit i=0; i<iterations; i++) {
    MS->rewind();
    u64bit           merNumber  = 0;
    u64bit           hm;
    u64bit           hf;

    u64bit           markedLastRound = distinct->numberMarked();

    //  We also need space to store the result of our hashing -- we can
    //  use any size hash for each step, they don't need to be the same
    //  size!
    //
    bitArray        *hashed = new bitArray(512);
    hashed->set(u64bitONE << (hashSize - partitionBits));
    hashed->clear();

    //  Build a new random hash for this iteration.
    //
    H.resetRandomness();

    SC = new speedCounter(" %8f Mbp (%8.5f Mbp/sec)\r", 1000000, 0xfffff, true);

    while (MS->nextMer()) {
      SC->tick();

      //  Pass 2?  If so, completely ignore this mer if pass1 told us
      //  to.  Could we get away without a goto?  Yes; if pass1 OR
      //  ignore==0 then proceed with the block.
      //
      if (pass2 && (ignore->get(merNumber) == 1))
        goto skip;

      hm = hf = H.hash(MS->theFMer());

      hm >>= partitionBits;
      hf  &= partitionMask;
      hf  ^= partition;

      //  Partitioned?  If not, both partitionMask and partition are
      //  zero, and so hf will be zero.  If so, hf is zero iff the low
      //  order bits are the same as the partition.
      //
      //  Then, mark that we've seen this hashed value (hm), and if
      //  this is the first time we've seen it, continue into the
      //  test.
      //
      //  If this isn't the first pass through, we can avoid lots of
      //  unnecessary updates if we check that this mer is not already
      //  marked as distinct -- that has been moved into the set()
      //  function.
      //
      if ((hf == u64bitZERO) &&
          (hashed->getAndSet(hm) == 0)) {
        distinct->set(merNumber);
      }

    skip:
      merNumber++;
    }

    //  Update our last mer number
    //
    distinct->setLastMer(merNumber);

    delete SC;
    delete hashed;

    fprintf(stderr, "hs="u32bitFMT" iter="u32bitFMT" mers="u64bitFMT" marked="u64bitFMT" tot="u64bitFMT"\n",
            hashSize, i, merNumber,
            distinct->numberMarked() - markedLastRound,
            distinct->numberMarked());

    if (distinct->numberMarked() - markedLastRound == 0)
      break;
  }

  delete distinct;
  delete ignore;

  delete MS;
  delete MSFR;
  delete FS;

  return(0);
}

