
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

#include "sqStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"

#include <vector>
#include <algorithm>

using namespace std;


static
void
writeToDumpFile(sqStore          *seq,
                ovOverlap       *overlap,
                ovFile          **dumpFile,
                uint64           *dumpLength,
                uint32           *iidToBucket,
                char             *ovlName) {
  uint32 df = iidToBucket[overlap->a_iid];

  //  If the dump file isn't opened, open it.

  if (dumpFile[df] == NULL) {
    char name[FILENAME_MAX];
    snprintf(name, FILENAME_MAX, "%s/tmp.sort.%04d", ovlName, df);
    fprintf(stderr, "-- Create bucket '%s'\n", name);
    dumpFile[df]   = new ovFile(seq, name, ovFileFullWriteNoCounts);
    dumpLength[df] = 0;
  }

  //  And write the overlap.

  dumpFile[df]->writeOverlap(overlap);
  dumpLength[df]++;
}



int
main(int argc, char **argv) {
  char const     *ovlName        = NULL;
  char const     *seqName        = NULL;
  char const     *cfgName        = NULL;

  double          maxErrorRate   = 1.0;

  bool            eValues        = false;
  char const     *configOut      = NULL;

  argc = AS_configure(argc, argv);

  vector<char const *>  err;
  int                   arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-S") == 0) {
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      cfgName = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      maxErrorRate = atof(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (ovlName == NULL)
    err.push_back("ERROR: No overlap store (-O) supplied.\n");

  if (seqName == NULL)
    err.push_back("ERROR: No sequence store (-S) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -S asm.seqStore -C ovStoreConfig [opts]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to overlap store to create\n");
    fprintf(stderr, "  -S asm.seqStore       path to a sequence store\n");
    fprintf(stderr, "  -C config             path to ovStoreConfig configuration file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e e                  filter overlaps above e fraction error\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Load the config, open the store, create a filter.

  ovStoreConfig    *config = new ovStoreConfig(cfgName);
  sqStore          *seq    = new sqStore(seqName);
  ovStoreFilter    *filter = new ovStoreFilter(seq, maxErrorRate);

  //  Figure out how many overlaps there are, quit if too many.

  uint32  maxID       = seq->sqStore_lastReadID();
  uint64  ovlsTotal   = 0;  //  Total in inputs.
  uint32  numInputs   = 0;

  fprintf(stderr, "\n");
  fprintf(stderr, "-- SCANNING INPUTS --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   Moverlaps\n");
  fprintf(stderr, "------------ ----------------------------------------\n");

  for (uint32 bb=1; bb<=config->numBuckets(); bb++) {
    for (uint32 ii=0; ii<config->numInputs(bb); ii++) {
      char const        *inputName = config->getInput(bb, ii);
      ovFile            *inputFile = new ovFile(seq, inputName, ovFileFull);

      ovlsTotal += inputFile->getCounts()->numOverlaps() * 2;
      numInputs += 1;

      fprintf(stderr, "%12.3f %40s\n",
              inputFile->getCounts()->numOverlaps() / 1000000.0,
              inputName);

      delete inputFile;
    }
  }

  fprintf(stderr, "------------ ----------------------------------------\n");
  fprintf(stderr, "%12.3f Moverlaps in inputs\n", ovlsTotal / 2 / 1000000.0);
  fprintf(stderr, "%12.3f Moverlaps to sort\n",   ovlsTotal     / 1000000.0);
  fprintf(stderr, "\n");

  if (ovlsTotal == 0)
    fprintf(stderr, "Found no overlaps to sort.\n");

  //  Load overlaps into memory.

  fprintf(stderr, "\n");
  fprintf(stderr, "Allocating space for " F_U64 " overlaps.\n", ovlsTotal);
  fprintf(stderr, "\n");

  ovOverlap      *ovls       = new ovOverlap [ovlsTotal];
  uint64          ovlsInput  = 0;
  uint64          ovlsLoaded = 0;

  fprintf(stderr, "\n");
  fprintf(stderr, "-- LOADING OVERLAPS --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Input       Loaded  Percent  Percent\n");
  fprintf(stderr, "   Moverlaps    Moverlaps   Loaded Complete\n");
  fprintf(stderr, "------------ ------------ -------- -------- ----------------------------------------\n");

  for (uint32 bb=1; bb<=config->numBuckets(); bb++) {
    for (uint32 ii=0; ii<config->numInputs(bb); ii++) {
      char const *inputName = config->getInput(bb, ii);

      fprintf(stderr, "%12.3f %12.3f %7.2f%% %7.2f%% %40s\n",
              ovlsInput   / 1000000.0,
              ovlsLoaded  / 1000000.0,
              100.0 * ovlsInput   / ovlsTotal,
              (ovlsInput == 0) ? (100.0) : (100.0 * ovlsLoaded / ovlsInput),
              inputName);

      ovOverlap foverlap;
      ovOverlap roverlap;

      ovFile   *inputFile = new ovFile(seq, inputName, ovFileFull);

      while (inputFile->readOverlap(&foverlap)) {
        filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r, and checks IDs

        ovlsInput += 2;

        //  Write the overlap if anything requests it.  These can be non-symmetric; e.g., if
        //  we only want to trim reads 1-1000, we'll not output any overlaps for a_iid > 1000.

        if ((foverlap.dat.ovl.forUTG == true) ||
            (foverlap.dat.ovl.forOBT == true) ||
            (foverlap.dat.ovl.forDUP == true))
          ovls[ovlsLoaded++] = foverlap;

        if ((roverlap.dat.ovl.forUTG == true) ||
            (roverlap.dat.ovl.forOBT == true) ||
            (roverlap.dat.ovl.forDUP == true))
          ovls[ovlsLoaded++] = roverlap;

        //  Report every 15.5 million overlaps (it's the millionth prime, why not).

        if ((ovlsLoaded % 15485863) == 0)
          fprintf(stderr, "%12.3f %12.3f %7.2f%% %7.2f%%\n",
                  ovlsInput   / 1000000.0,
                  ovlsLoaded  / 1000000.0,
                  100.0 * ovlsInput   / ovlsTotal,
                  (ovlsInput == 0) ? (100.0) : (100.0 * ovlsLoaded / ovlsInput));

        //  Make sure we didn't blow our space.

        assert(ovlsLoaded <= ovlsTotal);
      }

      delete inputFile;
    }
  }

  fprintf(stderr, "------------ ------------ -------- -------- ----------------------------------------\n");
  fprintf(stderr, "%12.3f %12.3f %7.2f%% %7.2f%%\n",
          ovlsInput   / 1000000.0,
          ovlsLoaded  / 1000000.0,
          100.0 * ovlsInput   / ovlsTotal,
          (ovlsInput == 0) ? (100.0) : (100.0 * ovlsLoaded / ovlsInput));

  //  Report what was filtered and loaded.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- OVERLAP FILTERING --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "TRIMMING OVERLAPS\n");
  fprintf(stderr, "Saved      " F_U64 " trimming overlaps\n", filter->savedTrimming());
  fprintf(stderr, "Discarded  " F_U64 " don't care\n",        filter->filteredNoTrim());
  fprintf(stderr, "\n");
  fprintf(stderr, "UNITIGGING OVERLAPS\n");
  fprintf(stderr, "Saved      " F_U64 " unitigging overlaps\n", filter->savedUnitigging());
  fprintf(stderr, "\n");
  fprintf(stderr, "Discarded  " F_U64 " low quality, more than %.4f fraction error\n", filter->filteredErate(), maxErrorRate);
  fprintf(stderr, "Discarded  " F_U64 " opposite orientation\n", filter->filteredFlipped());
  fprintf(stderr, "\n");

  delete filter;

  //  Sort the assorted overlaps.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- SORT OVERLAPS --\n");
  fprintf(stderr, "\n");

  sort(ovls, ovls + ovlsLoaded);

  //  Write.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- OUTPUT OVERLAPS --\n");
  fprintf(stderr, "\n");

  ovStoreWriter  *writer = new ovStoreWriter(ovlName, seq);

  for (uint64 oo=0; oo<ovlsLoaded; oo++)
    writer->writeOverlap(ovls + oo);

  delete    writer;
  delete [] ovls;

  //  Test.  Open the store and get the number of overlaps per read.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- TEST STORE --\n");
  fprintf(stderr, "\n");

  ovStore *tester = new ovStore(ovlName, seq);
  tester->testStore();
  delete    tester;

  //  And we have a store.  Cleanup and success!

  delete seq;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}
