#include "util++.H"
//#include "gmp.h"   //  Include this before bio++ so that kMer knows we have GMP!
#include "bio++.H"

#include "randomhash-shifting.H"
//#include "randomhash-mt.H"
//#include "randomhash-mod.H"
#include "markFile.H"

#warning NEED TO USE CANONICAL MERS!

class configuration {
public:
  configuration() {
    inName        = 0L;
    outName       = 0L;
    seed          = time(0L);
    iterations    = 16;
    partitionBits = 0;
    partitionMask = 0L;
    partition     = 0;
    merSize       = 100;
    hashSize      = 33;

    pass1         = false;
    pass2         = false;
    init          = false;
    propagate     = false;

    existingMask  = 0L;

    isFasta       = false;
    isMerStream   = false;
  };
  ~configuration() {
  };

  void usage(char *name) {
    fprintf(stderr, "usage: %s [options] -f in.fasta -o out-prefix\n", name);
    fprintf(stderr, "  -seed             seed randomness, a big number\n");
    fprintf(stderr, "  -iters            number of iterations\n");
    fprintf(stderr, "  -hashsize         total size of the hash, in bits\n");
    fprintf(stderr, "  -partitionbits    paritioned size of the hash, in bits\n"); 
    fprintf(stderr, "  -partition        hash partition to compute (0 <= p < 2^P)\n");
    fprintf(stderr, "  -mersize          mersize, default 100\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -pass1\n");
    fprintf(stderr, "  -pass2 pass1file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -f   read sequence from here\n");
    fprintf(stderr, "  -o   prefix output files with this\n");
    exit(1);
  };

  void parseArgs(int argc, char **argv) {
    bool     error  = false;
    int      arg    = 1;

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
      } else if (strcmp(argv[arg], "-partition") == 0) {
        partition  = strtou64bit(argv[++arg], 0L);

      } else if (strcmp(argv[arg], "-mersize") == 0) {
        merSize = strtou32bit(argv[++arg], 0L);

      } else if (strcmp(argv[arg], "-pass1") == 0) {
        pass1        = true;
        pass2        = false;
      } else if (strcmp(argv[arg], "-pass2") == 0) {
        pass1        = false;
        pass2        = true;
        existingMask = argv[++arg];

      } else if (strcmp(argv[arg], "-init") == 0) {
        init = true;

      } else if (strcmp(argv[arg], "-propagate") == 0) {
        propagate    = true;
        existingMask = argv[++arg];

      } else {
        fprintf(stderr, "invalid option '%s'\n", argv[arg]);
        error = true;
      }

      arg++;
    }

    //  Update partitionBits and the Mask to be what we need to mask out.
    //
    if (partitionBits > 0) {
      partitionBits = hashSize - partitionBits;
      partitionMask = u64bitMASK(partitionBits);
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

    if (!pass1 && !pass2 && !init && !propagate)
      fprintf(stderr, "ERROR: which pass?\n"), error = true;

    if (error)
      usage(argv[0]);
  };

public:
  char    *inName;
  char    *outName;
  u64bit   seed;
  u32bit   iterations;
  u64bit   partitionBits;
  u64bit   partitionMask;
  u64bit   partition;
  u32bit   merSize;
  u32bit   hashSize;

  bool     pass1;
  bool     pass2;
  bool     init;
  bool     propagate;

  char    *existingMask;

  bool     isFasta;
  bool     isMerStream;
};

configuration config;



void
doMark(merStream *MS) {

  //  Space to store our annotation of "mer that starts here is
  //  distinct".  This doesn't need to be in core, and it isn't.
  //
  //  Supposedly, the bitPackedFile ensures that bits are all set to
  //  zero.
  //
  markFile   *distinct = new markFile(config.outName);
  markFile   *ignore   = 0L;

  if (config.pass2) {
    if (!fileExists(config.existingMask)) {
      fprintf(stderr, "ERROR: -pass2 file '%s' doesn't exist!\n", config.existingMask);
      exit(1);
    }
    ignore = new markFile(config.existingMask);
  }

  //  Early versions used to pick a random hashSize, but why bother?
  //  Just use the largest we can fit in memory.
  //
  randomHash       H(config.seed, config.merSize, config.hashSize);

  for (u32bit i=0; i<config.iterations; i++) {
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
    hashed->set(u64bitONE << (config.hashSize - config.partitionBits));
    hashed->clear();

    //  Build a new random hash for this iteration.
    //
    H.resetRandomness();

    speedCounter *SC = new speedCounter(" %8f Mbp (%8.5f Mbp/sec)\r", 1000000, 0xfffff, true);

    while (MS->nextMer()) {
      SC->tick();

      //  Pass 2?  If so, completely ignore this mer if pass1 told us
      //  to.  Could we get away without a goto?  Yes; if pass1 OR
      //  ignore==0 then proceed with the block.
      //
      if (config.pass2 && (ignore->get(merNumber) == 1))
        goto skip;

      hm = hf = H.hash(MS->theCMer());

      hm >>= config.partitionBits;
      hf  &= config.partitionMask;
      hf  ^= config.partition;

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
            config.hashSize, i, merNumber,
            distinct->numberMarked() - markedLastRound,
            distinct->numberMarked());

    if (distinct->numberMarked() - markedLastRound == 0)
      break;
  }

  delete distinct;
  delete ignore;
}



void
doPropagate(merStream *MS) {

  //  Two step process:
  //
  //  First, we generate a random hash function, and load it
  //  with any mer marked in our input.  A mer should be marked if it's useful
  //  for overlapping.
  //
  //  Then, reading the input again, we add marks for any mer that
  //  hashes to something in the hash table.  Thus, we propagate marks
  //  to duplicate mers in the input.  Unfortunately, the world isn't
  //  perfect, and we'll include some unique mers here too.
  //
  //  If CPU turns out to be cheap, we could optimize the hash
  //  function to minimize the number of unique mers we catch.
  //  Construct several hash functions, load them, count how many new
  //  mers we include, pick the one that includes the smallest number.
  //
  //  It also might be just as cheap to load into core a subset of the
  //  mers we want to propagate, and then binary search each mer in
  //  the full set against those to figure out what to mark.


  markFile   *distinct   = new markFile(config.existingMask);

  randomHash  H(config.seed, config.merSize, config.hashSize);
  H.resetRandomness();

  bitArray        *hashed = new bitArray(512);
  hashed->set(u64bitONE << (config.hashSize - config.partitionBits));
  hashed->clear();

  MS->rewind();
  u64bit           merNumber  = 0;
  u64bit           hm;
  u64bit           hf;

  u64bit  collisions = 0;
  u64bit  marked     = 0;

  speedCounter *SC = new speedCounter(" %8f Mbp (%8.5f Mbp/sec)\r", 1000000, 0xfffff, true);
  while (MS->nextMer()) {
    SC->tick();

    if (distinct->get(merNumber) == 1) {
      hm = hf = H.hash(MS->theCMer());

      hm >>= config.partitionBits;
      hf  &= config.partitionMask;
      hf  ^= config.partition;

      //  See comments in doMark()

      if ((hf == u64bitZERO) &&
          (hashed->getAndSet(hm) == 1))
        collisions++;

      marked++;
    }

    merNumber++;
  }
  delete SC;

  fprintf(stderr, "propagate:  Had "u64bitFMT" collisions (out of "u64bitFMT" marked) when just building the hash.\n",
          collisions, marked);

  //  pass 2

  markFile   *propagated = new markFile(config.outName);

  MS->rewind();
  merNumber  = 0;
  collisions = 0;
  marked     = 0;

  SC = new speedCounter(" %8f Mbp (%8.5f Mbp/sec)\r", 1000000, 0xfffff, true);
  while (MS->nextMer()) {
    SC->tick();

    hm = hf = H.hash(MS->theCMer());

    hm >>= config.partitionBits;
    hf  &= config.partitionMask;
    hf  ^= config.partition;

    //  See comments in doMark()

    if ((hf == u64bitZERO) &&
        (hashed->getAndSet(hm) == 1)) {
      propagated->set(merNumber);
      marked++;
    }

    merNumber++;
  }
  delete SC;

  fprintf(stderr, "propagate:  Maked "u64bitFMT" mers during propagation.\n",
          marked);
  
  propagated->setLastMer(merNumber);

  delete hashed;

  fprintf(stderr, "hs="u32bitFMT" mers="u64bitFMT" marked="u64bitFMT"\n",
          config.hashSize, merNumber,
          propagated->numberMarked());

  delete distinct;
  delete propagated;
}




int
main(int argc, char **argv) {

  config.parseArgs(argc, argv);

  //  If we are initializing, we just convert the fastafile into a
  //  merStreamFile.
  //
  if (config.init) {
    fprintf(stderr, "Converting '%s' into '%s.merStream'\n", config.inName, config.outName);

    merStreamFileBuilder  *B = new merStreamFileBuilder(config.merSize,
                                                        config.inName,
                                                        config.outName);
    B->build();
    delete B;

    exit(0);
  }

  FastAstream          *FS   = 0L;
  merStream            *MS   = 0L;
  merStreamFileReader  *MSFR = 0L;

  if (config.isFasta) {
    FS   = new FastAstream(config.inName);
    MS   = new merStream(config.merSize, FS);
  } else if (config.isMerStream) {
    MSFR = new merStreamFileReader(config.inName, config.merSize);
    MS   = new merStream(MSFR);
  } else {
    fprintf(stderr, "ERROR!  Isn't fasta and isn't merstream?\n");
    exit(1);
  }

  if (config.pass1 || config.pass2)
    doMark(MS);

  if (config.propagate)
    doPropagate(MS);

  delete MS;
  delete MSFR;
  delete FS;

  return(0);
}
