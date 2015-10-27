
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
 *    src/AS_OVS/overlapStore.C
 *    src/AS_OVS/overlapStore.c
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2014-OCT-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

const char *mainid = "$Id$";

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "stddev.H"
#include "intervalList.H"


#define PLOT_READ_LENGTH      0x01
#define PLOT_OVERLAP_LENGTH   0x02
#define PLOT_OVERLAP_ERROR    0x04

#define OVL_5                 0x01
#define OVL_3                 0x02
#define OVL_CONTAINED         0x04
#define OVL_CONTAINER         0x08
#define OVL_PARTIAL           0x10

#define READ_FULL_COVERAGE    0x01
#define READ_MISSING_5        0x02
#define READ_MISSING_3        0x04
#define READ_MISSING_MIDDLE   0x08
#define READ_CONTAINED        0x10
#define READ_CONTAINER        0x20
#define READ_PARTIAL          0x40




void
plotCoverageHistogram(intervalList<uint32> &depth, uint32 numReads, char *outPrefix, char *label) {
  char  N[FILENAME_MAX];
  FILE *F;

  sprintf(N, "%s.coverage%s.histogram", outPrefix, label);
  errno = 0;

  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
  
  for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
    for (uint32 xx=depth.lo(ii); xx<depth.hi(ii); xx++)
      fprintf(F, "%u\t%.2f\n", xx, (double)depth.depth(ii) / numReads);

  fclose(F);

  sprintf(N, "%s.coverage%s.gp", outPrefix, label);
  errno = 0;

  F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

  fprintf(F, "set terminal 'png'\n");
  fprintf(F, "set output '%s.coverage%s.png'\n", outPrefix, label);
  fprintf(F, "plot '%s.coverage%s.histogram' using 1:2 with lines\n", outPrefix, label);

  fclose(F);

  sprintf(N, "gnuplot < %s.coverage%s.gp", outPrefix, label);
  system(N);
}





int
main(int argc, char **argv) {
  char           *gkpName      = NULL;
  char           *ovlName      = NULL;

  char           *outPrefix    = NULL;

  uint32          bgnID        = 0;
  uint32          endID        = UINT32_MAX;
  uint32          qryID        = 0;

  uint32          plotType     = 0;
  uint32          ovlSelect    = 0;
  double          ovlAtMost    = AS_OVS_encodeEvalue(1.0);
  double          ovlAtLeast   = AS_OVS_encodeEvalue(0.0);
  uint32          readSelect   = 0;

  bool            verbose      = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {

    if      (strcmp(argv[arg], "-G") == 0)
      gkpName = argv[++arg];

    else if (strcmp(argv[arg], "-O") == 0)
      ovlName = argv[++arg];


    else if (strcmp(argv[arg], "-o") == 0)
      outPrefix = argv[++arg];


    else if (strcmp(argv[arg], "-b") == 0)
      bgnID = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-e") == 0)
      endID = atoi(argv[++arg]);


    else if (strcmp(argv[arg], "-v") == 0)
      verbose = true;


    else if (strcmp(argv[arg], "-plot") == 0) {
      arg++;

      if      (strcmp(argv[arg], "read-length") == 0)
        plotType |= PLOT_READ_LENGTH;

      else if (strcmp(argv[arg], "overlap-length") == 0)
        plotType |= PLOT_OVERLAP_LENGTH;
      
      else if (strcmp(argv[arg], "overlap-error") == 0)
        plotType |= PLOT_OVERLAP_ERROR;

      else {
        fprintf(stderr, "ERROR: unknown -plot '%s'\n", argv[arg]);
        exit(1);
      }
    }


    else if (strcmp(argv[arg], "-overlap") == 0) {
      arg++;

      if      (strcmp(argv[arg], "5") == 0)
        ovlSelect |= OVL_5;

      else if (strcmp(argv[arg], "3") == 0)
        ovlSelect |= OVL_3;
      
      else if (strcmp(argv[arg], "contained") == 0)
        ovlSelect |= OVL_CONTAINED;

      else if (strcmp(argv[arg], "container") == 0)
        ovlSelect |= OVL_CONTAINER;

      else if (strcmp(argv[arg], "partial") == 0)
        ovlSelect |= OVL_PARTIAL;

      else if (strcmp(argv[arg], "atmost") == 0)
        ovlAtMost = atof(argv[++arg]);

      else if (strcmp(argv[arg], "atleast") == 0)
        ovlAtLeast = atof(argv[++arg]);

      else {
        fprintf(stderr, "ERROR: unknown -overlap '%s'\n", argv[arg]);
        exit(1);
      }
    }


    else if (strcmp(argv[arg], "-read") == 0) {
      arg++;

      if      (strcmp(argv[arg], "full-coverage") == 0)
        readSelect |= READ_FULL_COVERAGE;

      else if (strcmp(argv[arg], "missing-5") == 0)
        readSelect |= READ_MISSING_5;
      
      else if (strcmp(argv[arg], "missing-3") == 0)
        readSelect |= READ_MISSING_3;

      else if (strcmp(argv[arg], "missing-middle") == 0)
        readSelect |= READ_MISSING_MIDDLE;

      else if (strcmp(argv[arg], "contained") == 0)
        readSelect |= READ_CONTAINED;

      else if (strcmp(argv[arg], "container") == 0)
        readSelect |= READ_CONTAINER;

      else if (strcmp(argv[arg], "partial") == 0)
        readSelect |= READ_PARTIAL;

      else {
        fprintf(stderr, "ERROR: unknown -read '%s'\n", argv[arg]);
        exit(1);
      }
    }


    else {
      fprintf(stderr, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (ovlName == NULL)
    err++;
  if (outPrefix == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o outPrefix [-b bgnID] [-e endID] ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Generates statistics for an overlap store.  By default all possible classes\n");
    fprintf(stderr, "are generated, options can disable specific classes.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Logging\n");
    fprintf(stderr, "  -v                       For each read included in the statistics, report....stuff.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Data to Plot:\n");
    fprintf(stderr, "  -plot read-length        read length\n");
    fprintf(stderr, "  -plot overlap-length     overlap length\n");
    fprintf(stderr, "  -plot overlap-error      overlap error rate\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  What type of plot to generate.  Multiple types can be supplied.\n");
    fprintf(stderr, "  The plot type is appended to the outPrefix -- 'outPrefix.read-length.png'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Overlap Selection:\n");
    fprintf(stderr, "  -overlap 5               5' overlaps only\n");
    fprintf(stderr, "  -overlap 3               3' overlaps only\n");
    fprintf(stderr, "  -overlap contained       contained overlaps only\n");
    fprintf(stderr, "  -overlap container       container overlaps only\n");
    fprintf(stderr, "  -overlap partial         overlap is not valid for assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  An overlap is classified as exactly one of 5', 3', contained or container.\n");
    fprintf(stderr, "  By default, all overlaps are selected.  Specifying any of these options will\n");
    fprintf(stderr, "  restrict overlaps to just those classifications.  E.g., '-overlap 5 -overlap 3'\n");
    fprintf(stderr, "  will select dovetail overlaps off either end of the read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -overlap atmost x        at most fraction x error  (overlap-erate <= x)\n");
    fprintf(stderr, "  -overlap atleast x       at least fraction x error (x <= overlap-erate)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Overlaps can be further filtered by fraction error.  Usually, this will be an\n");
    fprintf(stderr, "  'atmost' filtering to use only the higher qualtiy overlaps.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Read Selection\n");
    fprintf(stderr, "  -read full-coverage      reads with full coverage in overlaps\n");
    fprintf(stderr, "  -read missing-5          reads with no 5' overlaps\n");
    fprintf(stderr, "  -read missing-3          reads with no 3' overlaps\n");
    fprintf(stderr, "  -read missing-middle     reads with no interior overlaps (implies both ends are covered)\n");
    fprintf(stderr, "                           (the bogart 'suspicious' reads)\n");
    fprintf(stderr, "  -read contained          reads that are contained in at least one other read\n");
    fprintf(stderr, "  -read container          reads that contain at least one other read\n");
    fprintf(stderr, "  -read partial            reads with overlaps that are not valid for assembly\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  By default, all reads are used.  Specifying any of these options will restrict\n");
    fprintf(stderr, "  reads to just those categories.  The decision is made after 'overlap selection' is\n");
    fprintf(stderr, "  performed.  Some of these are nonsense - '-overlap 5 -read full-coverage' will retain\n");
    fprintf(stderr, "  overlaps and reads that have full coverage in overlaps extending off the 5' end of the\n");
    fprintf(stderr, "  read (meaning that one of those overlaps must end exactly at the 3' end of the read).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -read contained          contained reads\n");
    fprintf(stderr, "  -read container          container reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A contained read has at least one container overlap.  Container read    -> ---------------\n");
    fprintf(stderr, "  A container read has at least one contained overlap.  Contained overlap ->      -----\n");
    fprintf(stderr, "\n");

    exit(1);
  }

  //  Set the default to 'all' if nothing set.

  if (plotType == 0)
    plotType = 0xff;

  if (ovlSelect == 0)
    ovlSelect = 0xff;

  if (readSelect == 0)
    readSelect = 0xff;

  //  Open inputs, find limits.

  gkStore    *gkpStore = new gkStore(gkpName);
  ovStore    *ovlStore = new ovStore(ovlName, gkpStore);

  if (endID > gkpStore->gkStore_getNumReads())
    endID = gkpStore->gkStore_getNumReads();

  if (endID < bgnID)
    fprintf(stderr, "ERROR: invalid bgn/end range bgn=%u end=%u; only %u reads in the store\n", bgnID, endID, gkpStore->gkStore_getNumReads()), exit(1);

  ovlStore->setRange(bgnID, endID);

  //  Allocate histograms.

  histogramStatistics   *readLength    = new histogramStatistics;
  histogramStatistics   *overlapLength = new histogramStatistics;
  histogramStatistics   *overlapError  = new histogramStatistics;

  //  Coverage interval lists, of all overlaps selected.

  intervalList<uint32>   coverage5;
  intervalList<uint32>   coverage3;
  intervalList<uint32>   coverageC;

  uint32                 numReads = 0;

  //  Open outputs.

  char N[FILENAME_MAX];
  sprintf(N, "%s.per-read.log", outPrefix);

  FILE  *LOG = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

  //  Compute!

  uint32        overlapsMax = 1024 * 1024;
  uint32        overlapsLen = 0;
  ovOverlap    *overlaps    = ovOverlap::allocateOverlaps(gkpStore, overlapsMax);

  overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);

  while (overlapsLen > 0) {
    intervalList<uint32>   cov;
    uint32                 covID = 0;

    bool    readFullCoverage  = true;
    bool    readMissingMiddle = true;

    bool    readCoverage5     = false;  //  Opposite of the user flag missing5
    bool    readCoverage3     = false;
    bool    readContained     = false;
    bool    readContainer     = false;
    bool    readPartial       = false;

    uint32  readID  = overlaps[0].a_iid;
    uint32  readLen = gkpStore->gkStore_getRead(readID)->gkRead_sequenceLength();

    for (uint32 oo=0; oo<overlapsLen; oo++) {
      bool  is5prime    = (overlaps[oo].overlapAEndIs5prime()  == true) && (ovlSelect & OVL_5);
      bool  is3prime    = (overlaps[oo].overlapAEndIs3prime()  == true) && (ovlSelect & OVL_3);
      bool  isContained = (overlaps[oo].overlapAIsContained()  == true) && (ovlSelect & OVL_CONTAINED);
      bool  isContainer = (overlaps[oo].overlapAIsContainer()  == true) && (ovlSelect & OVL_CONTAINER);
      bool  isPartial   = (overlaps[oo].overlapIsPartial()     == true) && (ovlSelect & OVL_PARTIAL);

      if ((is5prime    == false) &&
          (is3prime    == false) &&
          (isContained == false) &&
          (isContainer == false) &&
          (isPartial   == false))
        continue;

      if (overlaps[oo].evalue() < ovlAtLeast)
        continue;

      if (overlaps[oo].evalue() > ovlAtMost)
        continue;

      readCoverage5    |= is5prime;     //  If there is a 5' overlap, the read isn't missing 5' coverage
      readCoverage3    |= is3prime;
      readContained    |= isContained;  //  Read is contained in something else
      readContainer    |= isContainer;  //  Read is a container of somethign else
      readPartial      |= isPartial;

      cov.add(overlaps[oo].a_bgn(), overlaps[oo].a_end() - overlaps[oo].a_bgn());
    }

    //  If we filtered all the overlaps, just get out of here.  Yeah, some code duplication,
    //  but cleaner than sticking an if block around the rest of the loop.

    if (cov.numberOfIntervals() == 0) {
      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    cov.merge();

    //  Analyze the intervals, save per-read information to the log.


    uint32  lastInt    = cov.numberOfIntervals() - 1;
    uint32  bgn        = cov.lo(0);
    uint32  end        = cov.hi(lastInt);
    bool    contiguous = (lastInt == 0) ? true : false;

    readFullCoverage  = (lastInt == 0) && (bgn == 0) && (end == readLen);
    readMissingMiddle = (lastInt != 0);

    //  Check our analysis.  Well, we were going to assert that readCoverage5==false implies
    //  that bgn>0, but it is quite possible to have overlaps covering the whole read with
    //  no dovetail overlaps extending off (which is what the readCoverage flags indicate).
    //
    //assert((bgn > 0)         == (readCoverage5 == false));
    //assert((end < readLen-1) == (readCoverage3 == false));

    //  Report the read

    //  Classify the read, possibly skip output.

    if (false == (((readFullCoverage  == true)  && (readSelect & READ_FULL_COVERAGE)) ||
                  ((readCoverage5     == false) && (readSelect & READ_MISSING_5)) ||
                  ((readCoverage3     == false) && (readSelect & READ_MISSING_3)) ||
                  ((readMissingMiddle == true)  && (readSelect & READ_MISSING_MIDDLE)) ||
                  ((readContained     == false) && (readSelect & READ_CONTAINED)) ||
                  ((readContainer     == false) && (readSelect & READ_CONTAINER)) ||
                  ((readPartial       == false) && (readSelect & READ_PARTIAL)))) {
      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    fprintf(LOG, "%u\tlen\t%u\trange\t%u\t%u\tregions\t%u",
            readID, readLen,
            bgn, end, cov.numberOfIntervals());

    if (cov.numberOfIntervals() > 1)
      for (uint32 ii=0; ii<cov.numberOfIntervals(); ii++)
        fprintf(LOG, "\t%u\t%u", cov.lo(ii), cov.hi(ii));

    fprintf(LOG, "\n");

    //  If we made it here, the read has overlaps we care about.  If we're plotting read length, add
    //  one point to the readLength data set.

    if (plotType | PLOT_READ_LENGTH)
      readLength->add(readLen);

    numReads++;

    //  For the other plot types, we need to run through all the overlaps again, classifying them
    //  again, then saving the overlap length.

    for (uint32 oo=0; oo<overlapsLen; oo++) {
      bool  is5prime    = (overlaps[oo].overlapAEndIs5prime()  == true) && (ovlSelect & OVL_5);
      bool  is3prime    = (overlaps[oo].overlapAEndIs3prime()  == true) && (ovlSelect & OVL_3);
      bool  isContained = (overlaps[oo].overlapAIsContained()  == true) && (ovlSelect & OVL_CONTAINED);
      bool  isContainer = (overlaps[oo].overlapAIsContainer()  == true) && (ovlSelect & OVL_CONTAINER);
      bool  isPartial   = (overlaps[oo].overlapIsPartial()     == true) && (ovlSelect & OVL_PARTIAL);

      if ((is5prime    == false) &&
          (is3prime    == false) &&
          (isContained == false) &&
          (isContainer == false) &&
          (isPartial   == false))
        continue;

      if (overlaps[oo].evalue() < ovlAtLeast)
        continue;

      if (overlaps[oo].evalue() > ovlAtMost)
        continue;

      //  No changes for 5' overlaps, plot from bgn==0 to end==whatever.
      if (is5prime) {
        coverage5.add(overlaps[oo].a_bgn(), overlaps[oo].a_end() - overlaps[oo].a_bgn());
      }

      //  For 3' overlaps, we need to reverse-complement, so the 3' end becomes position 0.
      if (is3prime) {
        uint32  bgn = readLen - overlaps[oo].a_end();
        uint32  end = readLen - overlaps[oo].a_bgn();

        coverage3.add(bgn, end - bgn);
      }

      //  For container reads (showing the location of contained reads
      //  in the container), normalize the A read length to 0..1.
      if (isContainer) {
        uint32  bgn = overlaps[oo].a_bgn();
        uint32  end = overlaps[oo].a_end();

        coverageC.add(bgn, end-bgn);
      }


      if (plotType | PLOT_OVERLAP_LENGTH) {
        overlapLength->add(overlaps[oo].span());
      }

      if (plotType | PLOT_OVERLAP_ERROR) {
        overlapError->add(overlaps[oo].evalue());
      }
    }

    overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
  }

  fclose(LOG);  //  Done with logging.


  readLength->finalizeData();
  overlapLength->finalizeData();
  overlapError->finalizeData();

  intervalList<uint32>  depth3(coverage3);  //  Convert the list of overlap intervals to depths.
  intervalList<uint32>  depth5(coverage5);
  intervalList<uint32>  depthC(coverageC);


  //char  N[FILENAME_MAX];
  FILE *F;

  fprintf(stderr, "\n\n");

  if (plotType | PLOT_READ_LENGTH) {
    fprintf(stderr, "readLength:     %9"F_U64P" objs  %8.2f +- %8.2f  median %8"F_U64P" mad %8"F_U64P"\n",
            readLength->numberOfObjects(),
            readLength->mean(),   readLength->stddev(),
            readLength->median(), readLength->mad());

    sprintf(N, "%s.read-length.histogram", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    readLength->writeHistogram(F, "read-length");
    fclose(F);

    sprintf(N, "%s.read-length.gp", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    fprintf(F, "set terminal 'png'\n");
    fprintf(F, "set output '%s.read-length.png'\n", outPrefix);
    fprintf(F, "plot '%s.read-length.histogram' using 1:2 with lines\n", outPrefix);
    fclose(F);

    sprintf(N, "gnuplot < %s.read-length.gp", outPrefix);
    system(N);
  }


  if (plotType | PLOT_OVERLAP_LENGTH) {
    fprintf(stderr, "overlapLength:  %9"F_U64P" objs  %8.2f +- %8.2f  median %8"F_U64P" mad %8"F_U64P"\n",
            overlapLength->numberOfObjects(),
            overlapLength->mean(),   overlapLength->stddev(),
            overlapLength->median(), overlapLength->mad());

    sprintf(N, "%s.overlap-length.histogram", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    overlapLength->writeHistogram(F, "overlap-length");
    fclose(F);

    sprintf(N, "%s.overlap-length.gp", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    fprintf(F, "set terminal 'png'\n");
    fprintf(F, "set output '%s.overlap-length.png'\n", outPrefix);
    fprintf(F, "plot '%s.overlap-length.histogram' using 1:2 with lines\n", outPrefix);
    fclose(F);

    sprintf(N, "gnuplot < %s.overlap-length.gp", outPrefix);
    system(N);
  }


  if (plotType | PLOT_OVERLAP_ERROR) {
    fprintf(stderr, "overlapError:   %9"F_U64P" objs  %8.4f +- %8.4f  median %8.4f mad %8.4f\n",
            overlapError->numberOfObjects(),
            AS_OVS_decodeEvalue(overlapError->mean()),   AS_OVS_decodeEvalue(overlapError->stddev()),
            AS_OVS_decodeEvalue(overlapError->median()), AS_OVS_decodeEvalue(overlapError->mad()));

    sprintf(N, "%s.overlap-error.histogram", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    overlapError->writeHistogram(F, "overlap-error");
    fclose(F);

    sprintf(N, "%s.overlap-error.gp", outPrefix);
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    fprintf(F, "set terminal 'png'\n");
    fprintf(F, "set output '%s.overlap-error.png'\n", outPrefix);
    fprintf(F, "plot '%s.overlap-error.histogram' using 1:2 with lines\n", outPrefix);
    fclose(F);

    sprintf(N, "gnuplot < %s.overlap-error.gp", outPrefix);
    system(N);
  }


  plotCoverageHistogram(depth5, numReads, outPrefix, "5prime");
  plotCoverageHistogram(depth3, numReads, outPrefix, "3prime");
  plotCoverageHistogram(depthC, numReads, outPrefix, "contained");


  delete overlapError;
  delete overlapLength;
  delete readLength;

  delete ovlStore;
  delete gkpStore;

  exit(0);
}
