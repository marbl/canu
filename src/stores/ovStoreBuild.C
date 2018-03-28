
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
 *  This file is derived from:
 *
 *    src/AS_OVS/overlapStoreBuild.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-APR-02 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-JUL-31 to 2015-SEP-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"
#include "ovStoreConfig.H"

#include "AS_UTL_decodeRange.H"

#include <vector>
#include <algorithm>

using namespace std;


static
void
writeToDumpFile(gkStore          *gkp,
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
    dumpFile[df]   = new ovFile(gkp, name, ovFileFullWriteNoCounts);
    dumpLength[df] = 0;
  }

  //  And write the overlap.

  dumpFile[df]->writeOverlap(overlap);
  dumpLength[df]++;
}



int
main(int argc, char **argv) {
  char           *ovlName        = NULL;
  char           *gkpName        = NULL;
  char           *cfgName        = NULL;

  double          maxErrorRate   = 1.0;

  bool            eValues        = false;
  char           *configOut      = NULL;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

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

  if (gkpName == NULL)
    err.push_back("ERROR: No gatekeeper store (-G) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -O asm.ovlStore -G asm.gkpStore -C ovStoreConfig [opts]\n", argv[0]);
    fprintf(stderr, "  -O asm.ovlStore       path to overlap store to create\n");
    fprintf(stderr, "  -G asm.gkpStore       path to gatekeeper store\n");
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
  gkStore          *gkp    = gkStore::gkStore_open(gkpName);
  ovStoreFilter    *filter = new ovStoreFilter(gkp, maxErrorRate);

  //  Figure out how many overlaps there are, quit if too many.

  uint32  maxID       = gkp->gkStore_getNumReads();
  uint64  totOverlaps = 0;  //  Total in inputs.
  uint32  numInputs   = 0;

  fprintf(stderr, "\n");
  fprintf(stderr, "-- SCANNING INPUTS --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "      Molaps\n");
  fprintf(stderr, "------------ ----------------------------------------\n");

  for (uint32 bb=1; bb<=config->numBuckets(); bb++) {
    for (uint32 ii=0; ii<config->numInputs(bb); ii++) {
      char              *inputName = config->getInput(bb, ii);
      ovFile            *inputFile = new ovFile(gkp, inputName, ovFileFull);

      totOverlaps += inputFile->getCounts()->numOverlaps() * 2;
      numInputs   += 1;

      fprintf(stderr, "%12.3f %40s\n", 
              inputFile->getCounts()->numOverlaps() / 1000000.0,
              inputName);

      delete inputFile;
    }
  }

  fprintf(stderr, "------------ ----------------------------------------\n");
  fprintf(stderr, "%12.3f overlaps in inputs\n", totOverlaps / 2 / 1000000.0);
  fprintf(stderr, "%12.3f overlaps to sort\n",   totOverlaps     / 1000000.0);
  fprintf(stderr, "\n");

  if (totOverlaps == 0)
    fprintf(stderr, "Found no overlaps to sort.\n"), exit(1);

  //  Load overlaps into memory.

  fprintf(stderr, "\n");
  fprintf(stderr, "Allocating space for " F_U64 " overlaps.\n", totOverlaps);
  fprintf(stderr, "\n");

  ovOverlap      *ovls    = ovOverlap::allocateOverlaps(gkp, totOverlaps);
  uint64          ovlsLen = 0;

  fprintf(stderr, "\n");
  fprintf(stderr, "-- LOADING OVERLAPS --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Input       Loaded Percent\n");
  fprintf(stderr, "      Molaps       Molaps  Loaded\n");
  fprintf(stderr, "------------ ------------ ------- ----------------------------------------\n");

  for (uint32 bb=1; bb<=config->numBuckets(); bb++) {
    for (uint32 ii=0; ii<config->numInputs(bb); ii++) {
      char     *inputName = config->getInput(bb, ii);

      fprintf(stderr, "%12.3f %12.3f %6.2f%% %40s\n",
              totOverlaps / 1000000.0,
              ovlsLen     / 1000000.0,
              0.0,
              inputName);

      ovOverlap foverlap(gkp);
      ovOverlap roverlap(gkp);

      ovFile   *inputFile = new ovFile(gkp, inputName, ovFileFull);

      while (inputFile->readOverlap(&foverlap)) {
        filter->filterOverlap(foverlap, roverlap);  //  The filter copies f into r, and checks IDs

        //  Write the overlap if anything requests it.  These can be non-symmetric; e.g., if
        //  we only want to trim reads 1-1000, we'll not output any overlaps for a_iid > 1000.

        if ((foverlap.dat.ovl.forUTG == true) ||
            (foverlap.dat.ovl.forOBT == true) ||
            (foverlap.dat.ovl.forDUP == true))
          ovls[ovlsLen++] = foverlap;

        if ((roverlap.dat.ovl.forUTG == true) ||
            (roverlap.dat.ovl.forOBT == true) ||
            (roverlap.dat.ovl.forDUP == true))
          ovls[ovlsLen++] = roverlap;

        //  Report every 15.5 million overlaps (it's the millionth prime, why not).

        if ((ovlsLen % 15485863) == 0)
          fprintf(stderr, "%12.3f %12.3f %6.2f%%\n",
                  totOverlaps / 1000000.0,
                  ovlsLen     / 1000000.0,
                  0.0);

        //  Make sure we didn't blow our space.

        assert(ovlsLen <= totOverlaps);
      }

      delete inputFile;
    }
  }

  fprintf(stderr, "------------ ------------ ------- ----------------------------------------\n");
  fprintf(stderr, "%12.3f %12.3f %6.2f%%\n",
          totOverlaps / 1000000.0,
          ovlsLen     / 1000000.0,
          0.0);

  //  Report what was filtered and loaded.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- OVERLAP FILTERING --\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "DEDUPE OVERLAPS\n");
  fprintf(stderr, "Saved      " F_U64 " dedupe overlaps\n",          filter->savedDedupe());
  fprintf(stderr, "Discarded  " F_U64 " don't care\n",               filter->filteredNoDedupe());
  fprintf(stderr, "Discarded  " F_U64 " different library\n",        filter->filteredNotDupe());
  fprintf(stderr, "Discarded  " F_U64 " obviously not duplicates\n", filter->filteredDiffLib());
  fprintf(stderr, "\n");
  fprintf(stderr, "TRIMMING OVERLAPS\n");
  fprintf(stderr, "Saved      " F_U64 " trimming overlaps\n", filter->savedTrimming());
  fprintf(stderr, "Discarded  " F_U64 " don't care\n",        filter->filteredNoTrim());
  fprintf(stderr, "Discarded  " F_U64 " too similar\n",       filter->filteredBadTrim());
  fprintf(stderr, "Discarded  " F_U64 " too short\n",         filter->filteredShortTrim());
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

#ifdef _GLIBCXX_PARALLEL
  //  If we have the parallel STL, don't use it!  Sort is not inplace!
  __gnu_sequential::
#endif
  sort(ovls, ovls + ovlsLen);

  //  Write.

  fprintf(stderr, "\n");
  fprintf(stderr, "-- OUTPUT OVERLAPS --\n");
  fprintf(stderr, "\n");

  ovStoreWriter  *store = new ovStoreWriter(ovlName, gkp);

  for (uint64 oo=0; oo<ovlsLen; oo++)
    store->writeOverlap(ovls + oo);

  delete    store;
  delete [] ovls;

  gkp->gkStore_close();

  //  And we have a store.

  fprintf(stderr, "Bye.\n");

  exit(0);
}
