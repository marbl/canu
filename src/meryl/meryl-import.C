
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
 *    Brian P. Walenz beginning on 2018-NOV-02
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "kmers.H"
#include "sequence.H"
#include "strings.H"
#include "bits.H"

#include "merylCountArray.H"


int
main(int argc, char **argv) {
  char   *inputName    = NULL;
  char   *outputDBname = NULL;
  uint32  kLen         = 0;
  uint64  maxValue     = 0;
  bool    doMultiSet   = false;

  bool    useC         = true;
  bool    useF         = false;

  uint32  threads      = 1;
  uint64  memory       = 8;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-kmers") == 0) {
      inputName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-k") == 0) {
      kLen = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-maxvalue") == 0) {
      maxValue = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-multiset") == 0) {
      doMultiSet = true;

    } else if (strcmp(argv[arg], "-forward") == 0) {
      useC = false;
      useF = true;

    } else if (strcmp(argv[arg], "-reverse") == 0) {
      useC = false;
      useF = false;

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {   //  Not implemented.  If implemented, merylCountArray::initializeValues()
      memory = strtouint64(argv[++arg]);              //  needs to return a memory size, etc, etc.

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputName == NULL)
    err.push_back("No input kmer file (-kmers) supplied.\n");
  if (outputDBname == NULL)
    err.push_back("No output database name (-output) supplied.\n");
  if (kLen == 0)
    err.push_back("No kmer size (-k) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [...] -k <kmer-size> -kmers <input-kmers> -output <db.meryl>\n", argv[0]);
    fprintf(stderr, "  Loads the kmers and values listed in <input-kmers> into a meryl kmer database.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS and OUTPUTS\n");
    fprintf(stderr, "  -kmers <input-kmers>  A file consisting of kmers and values, one per line, separated\n");
    fprintf(stderr, "                        by white space ('AGTTGCC 4').  Order of kmers is not important.\n");
    fprintf(stderr, "                        Duplicate kmers will be handled according to the -multiset\n");
    fprintf(stderr, "                        option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                        A persistent value can be specified as '#<value>' (e.g., '#3')\n");
    fprintf(stderr, "                        All kmers with no value after this line will use this value.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -k <size>             The size of a kmer, in bases.  Setting this larger than the\n");
    fprintf(stderr, "                        kmers in the input will probably lead to a crash.  Setting it\n");
    fprintf(stderr, "                        smaller will result in only the left-most bases being used.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -output <db.meryl>    Create (or overwrite) meryl database 'database.meryl'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS\n");
    fprintf(stderr, "  -multiset             Write duplicate kmers in the input to the database as individual\n");
    fprintf(stderr, "                        entries.  A kmer AGTTGCC in the input twice with values 4 and 7\n");
    fprintf(stderr, "                        will be listed in the output database twice, once with value 4,\n");
    fprintf(stderr, "                        and once with value 7.  Without this option, the values will be\n");
    fprintf(stderr, "                        summed: AGTTGCC will be listed once with value 11.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -maxvalue <value>     An optional memory and time optimization, useful if your values\n");
    fprintf(stderr, "                        are randomly distributed and below some known maximum value.\n");
    fprintf(stderr, "                        For data whose values are the counts from actual data, it is\n");
    fprintf(stderr, "                        probably best to not set this option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -forward              By default, the canonical kmer is loaded into the database.  These\n");
    fprintf(stderr, "  -reverse              options force either the forward or reverse-complement kmer to be\n");
    fprintf(stderr, "                        loaded instead.  These options are mutually exclusive.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -threads <t>          Use <t> compute threads when sorting and writing data.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -memory <m>           (accepted, but not implemented, sorry)\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  //  Figure out the kmer size.  We need this to set up encoding parameters.

  kmerTiny::setSize(kLen);

  //  Decide on some parameters.

  uint32  wPrefix   = 10;
  uint32  nPrefix   = 1 << wPrefix;

  uint32  wData     = 2 * kmerTiny::merSize() - wPrefix;
  uint64  wDataMask = uint64MASK(wData);

  //  Open the input kmer file, allocate space for reading kmer lines.

  FILE         *K    = AS_UTL_openInputFile(inputName);
  uint32        Llen = 0;
  uint32        Lmax = 1023;
  char         *L    = new char [Lmax + 1];

  splitToWords  W;

  kmerTiny      kmerF;
  kmerTiny      kmerR;

  uint64        nKmers = 0;

  //  Allocate a bunch of counting arrays.  The naming here follows merylOp-count.C.

  merylCountArray<uint64>  *data  = new merylCountArray<uint64> [nPrefix];

  for (uint32 pp=0; pp<nPrefix; pp++) {
    data[pp].initialize(pp, wData);
    data[pp].initializeValues(maxValue);
    data[pp].enableMultiSet(doMultiSet);
  }

  //  Read each kmer and value, stuff into a merylCountArray, writing when the array is full.

  uint64  persistentValue = 1;

  while (AS_UTL_readLine(L, Llen, Lmax, K) == true) {
    W.split(L);

    if (W.numWords() == 0)
      continue;

    //  Decode the line, make a kmer.

    char   *kstr = W[0];
    uint64  vv   = persistentValue;

    if (kstr[0] == '#') {
      persistentValue = strtouint64(kstr + 1);
      continue;
    }

    if (W.numWords() > 1)
      vv = W.touint64(1);

    for (uint32 ii=0; kstr[ii]; ii++)
      kmerF.addR(kstr[ii]);

    kmerR = kmerF;
    kmerR.reverseComplement();

    //  Decide to use the F or the R kmer.

    if (useC == true) {
      useF = (kmerF < kmerR) ? true : false;
    }

    //  And use it.

    uint64  pp = (useF == true) ? ((uint64)kmerF >> wData)     : ((uint64)kmerR >> wData);
    uint64  mm = (useF == true) ? ((uint64)kmerF  & wDataMask) : ((uint64)kmerR  & wDataMask);

    assert(pp < nPrefix);

    data[pp].add(mm);
    data[pp].addValue(vv);

    nKmers++;
  }

  //  All data loaded, cleanup.

  fprintf(stderr, "Found %lu kmers in the input.\n", nKmers);



  //  And dump to the output.

  kmerCountFileWriter   *output = new kmerCountFileWriter(outputDBname);

  output->initialize(wPrefix);

  kmerCountBlockWriter  *writer = output->getBlockWriter();

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<output->numberOfFiles(); ff++) {
    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, output->firstPrefixInFile(ff), output->lastPrefixInFile(ff));

    for (uint64 pp=output->firstPrefixInFile(ff); pp <= output->lastPrefixInFile(ff); pp++) {
      data[pp].countKmers();                //  Convert the list of kmers into a list of (kmer, count).
      data[pp].dumpCountedKmers(writer);    //  Write that list to disk.
      data[pp].removeCountedKmers();        //  And remove the in-core data.
    }
  }

  writer->finish();

  delete    writer;
  delete    output;
  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}
