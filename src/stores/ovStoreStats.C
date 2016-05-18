
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
 *    Brian P. Walenz beginning on 2015-OCT-27
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-31
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "stddev.H"
#include "intervalList.H"


#define OVL_5                 0x01
#define OVL_3                 0x02
#define OVL_CONTAINED         0x04
#define OVL_CONTAINER         0x08
#define OVL_PARTIAL           0x10


//  Should count unique-contained and repeat-contained separately from unique and repeat
//  uniq-anchor is also 'plausible chimera'

//  no-5-prime includes things that entirely cover the read, just no overhang

int
main(int argc, char **argv) {
  char           *gkpName        = NULL;
  char           *ovlName        = NULL;
  char           *outPrefix      = NULL;

  uint32          bgnID          = 0;
  uint32          endID          = UINT32_MAX;

  uint32          ovlSelect      = 0;
  double          ovlAtMost      = AS_OVS_encodeEvalue(1.0);
  double          ovlAtLeast     = AS_OVS_encodeEvalue(0.0);

  double          expectedMean   = 30.0;
  double          expectedStdDev =  7.0;

  bool            toFile         = true;

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


    else if (strcmp(argv[arg], "-C") == 0) {
      expectedMean   = atof(argv[++arg]);
      expectedStdDev = atof(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-c") == 0)
      toFile = false;


    else if (strcmp(argv[arg], "-b") == 0)
      bgnID = atoi(argv[++arg]);

    else if (strcmp(argv[arg], "-e") == 0)
      endID = atoi(argv[++arg]);


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
    fprintf(stderr, "  -C mean stddev           Expect coverage at mean +- stddev\n");
    fprintf(stderr, "  -c                       Write stats to stdout, not to a file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  outPrefix.per-read.log   One line per read, giving readID, read length and classification.\n");
    fprintf(stderr, "  outPrefix.summary        The primary statistical output.\n");
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
    fprintf(stderr, "  A contained read has at least one container overlap.  Container read    -> ---------------\n");
    fprintf(stderr, "  A container read has at least one contained overlap.  Contained overlap ->      -----\n");
    fprintf(stderr, "\n");

    exit(1);
  }

  //  Set the default to 'all' if nothing set.

  if (ovlSelect == 0)
    ovlSelect = 0xff;

  //  Open inputs, find limits.

  gkStore    *gkpStore = gkStore::gkStore_open(gkpName);
  ovStore    *ovlStore = new ovStore(ovlName, gkpStore);

  if (endID > gkpStore->gkStore_getNumReads())
    endID = gkpStore->gkStore_getNumReads();

  if (endID < bgnID)
    fprintf(stderr, "ERROR: invalid bgn/end range bgn=%u end=%u; only %u reads in the store\n", bgnID, endID, gkpStore->gkStore_getNumReads()), exit(1);

  ovlStore->setRange(bgnID, endID);

  //  Allocate output histograms.

  histogramStatistics   *readNoOlaps         = new histogramStatistics;  //  Bad reads!  (read length)
  histogramStatistics   *readHole            = new histogramStatistics;
  histogramStatistics   *readHump            = new histogramStatistics;
  histogramStatistics   *readNo5             = new histogramStatistics;
  histogramStatistics   *readNo3             = new histogramStatistics;

  histogramStatistics   *olapHole            = new histogramStatistics;  //  Hole size (sum of holes if more than one)
  histogramStatistics   *olapHump            = new histogramStatistics;  //  Hump size (sum of humps if more than one)
  histogramStatistics   *olapNo5             = new histogramStatistics;  //  5' uncovered size
  histogramStatistics   *olapNo3             = new histogramStatistics;  //  3' uncovered size

  histogramStatistics   *readLowCov          = new histogramStatistics;  //  Good reads!  (read length)
  histogramStatistics   *readUnique          = new histogramStatistics;
  histogramStatistics   *readRepeatCont      = new histogramStatistics;
  histogramStatistics   *readRepeatDove      = new histogramStatistics;
  histogramStatistics   *readSpanRepeat      = new histogramStatistics;
  histogramStatistics   *readUniqRepeatCont  = new histogramStatistics;
  histogramStatistics   *readUniqRepeatDove  = new histogramStatistics;
  histogramStatistics   *readUniqAnchor      = new histogramStatistics;

  histogramStatistics   *covrLowCov          = new histogramStatistics;  //  Good reads!  (overlap length)
  histogramStatistics   *covrUnique          = new histogramStatistics;
  histogramStatistics   *covrRepeatCont      = new histogramStatistics;
  histogramStatistics   *covrRepeatDove      = new histogramStatistics;
  histogramStatistics   *covrSpanRepeat      = new histogramStatistics;
  histogramStatistics   *covrUniqRepeatCont  = new histogramStatistics;
  histogramStatistics   *covrUniqRepeatDove  = new histogramStatistics;
  histogramStatistics   *covrUniqAnchor      = new histogramStatistics;

  histogramStatistics   *olapLowCov          = new histogramStatistics;  //  Good reads!  (overlap length)
  histogramStatistics   *olapUnique          = new histogramStatistics;
  histogramStatistics   *olapRepeatCont      = new histogramStatistics;
  histogramStatistics   *olapRepeatDove      = new histogramStatistics;
  histogramStatistics   *olapSpanRepeat      = new histogramStatistics;
  histogramStatistics   *olapUniqRepeatCont  = new histogramStatistics;
  histogramStatistics   *olapUniqRepeatDove  = new histogramStatistics;
  histogramStatistics   *olapUniqAnchor      = new histogramStatistics;

  //  Coverage interval lists, of all overlaps selected.

  //  Open outputs.

  char N[FILENAME_MAX];
  sprintf(N, "%s.per-read.log", outPrefix);

  FILE  *LOG = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

  //  Compute!

  uint32                 overlapsMax = 1024;

  uint32                 overlapsLen = 0;
  ovOverlap             *overlaps    = ovOverlap::allocateOverlaps(gkpStore, overlapsMax);

  overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);

  while (overlapsLen > 0) {
    uint32  readID  = overlaps[0].a_iid;
    uint32  readLen = gkpStore->gkStore_getRead(readID)->gkRead_sequenceLength();

    intervalList<uint32>   cov;
    uint32                 covID = 0;

    bool    readCoverage5     = false;
    bool    readCoverage3     = false;
    bool    readContained     = false;
    bool    readContainer     = false;
    bool    readPartial       = false;

    for (uint32 oo=0; oo<overlapsLen; oo++) {
      bool  is5prime    = (overlaps[oo].overlapAEndIs5prime()  == true) && (ovlSelect & OVL_5)         && (overlaps[oo].overlap5primeIsPartial() == false);
      bool  is3prime    = (overlaps[oo].overlapAEndIs3prime()  == true) && (ovlSelect & OVL_3)         && (overlaps[oo].overlap3primeIsPartial() == false);
      bool  isContained = (overlaps[oo].overlapAIsContained()  == true) && (ovlSelect & OVL_CONTAINED);
      bool  isContainer = (overlaps[oo].overlapAIsContainer()  == true) && (ovlSelect & OVL_CONTAINER);
      bool  isPartial   = (overlaps[oo].overlapIsPartial()     == true) && (ovlSelect & OVL_PARTIAL);

      //  Ignore the overlap?

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
      readNoOlaps->add(readLen);

      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    //  Generate a depth-of-coverage map, then merge intervals

    intervalList<uint32>  depth(cov);

    cov.merge();

    //  Analyze the intervals, save per-read information to the log.

    uint32  lastInt           = cov.numberOfIntervals() - 1;
    uint32  bgn               = cov.lo(0);
    uint32  end               = cov.hi(lastInt);
    bool    contiguous        = (lastInt == 0) ? true : false;

    bool    readFullCoverage  = (lastInt == 0) && (bgn == 0) && (end == readLen);
    bool    readMissingMiddle = (lastInt != 0);

    uint32  holeSize          = 0;
    uint32  no5Size           = bgn;
    uint32  no3Size           = readLen - end;

    for (uint32 ii=1; ii<cov.numberOfIntervals(); ii++)
      holeSize += cov.lo(ii) - cov.hi(ii-1);

    //  Handle bad cases.  If it's a partial overlap, ignore the is5prime and is3prime markings.


    if (readMissingMiddle == true) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "middle-missing");
      readHole->add(readLen);
      olapHole->add(holeSize);

      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    if ((readCoverage5 == false) && (readCoverage3 == false) && (readContained == false) && (readPartial == false)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "middle-only");
      readHump->add(readLen);
      olapHump->add(no5Size + no3Size);

      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    if ((readCoverage5 == false) && (readContained == false) && (readPartial == false)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "no-5-prime");
      readNo5->add(readLen);
      olapNo5->add(no5Size);

      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    if ((readCoverage3 == false) && (readContained == false) && (readPartial == false)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "no-3-prime");
      readNo3->add(readLen);
      olapNo3->add(no3Size);

      overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
      continue;
    }

    //  Handle good cases.  For partial overlaps, bgn and end are not the extent of the read.

    if (readPartial == false) {
      assert(bgn == 0);
      assert(end == readLen);
      assert(contiguous == true);
      assert(readFullCoverage == true);
    }

    //  Compute mean and std.dev of coverage.  From this, we decide if the read is 'unique',
    //  'repeat' or 'mixed'.  If 'mixed', we then need to decide if the read spans a repeat, or
    //  joins unique and repeat.

    double  covMean   = 0;
    double  covStdDev = 0;

    for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
      covMean += (depth.hi(ii) - depth.lo(ii)) * depth.depth(ii);

    covMean /= readLen;

    for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
      covStdDev += (depth.hi(ii) - depth.lo(ii)) * (depth.depth(ii) - covMean) * (depth.depth(ii) - covMean);

    covStdDev = sqrt(covStdDev / (readLen - 1));

    //  Classify each interval as either 'l'owcoverage, 'u'nique or 'r'epeat.

    char *classification = new char [depth.numberOfIntervals()];

    for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++) {
      if        (depth.depth(ii) < expectedMean - 3 * expectedStdDev) {
        classification[ii] = 'l';

      } else if (depth.depth(ii) < expectedMean + 3 * expectedStdDev) {
        classification[ii] = 'u';

      } else {
        classification[ii] = 'r';
      }
    }

    //  Try to detect if a read is part unique and part repeat.

    bool   isLowCov     = false;
    bool   isUnique     = false;
    bool   isRepeat     = false;
    bool   isSpanRepeat = false;
    bool   isUniqRepeat = false;
    bool   isUniqAnchor = false;

    int32  bgni = 0;
    int32  endi = depth.numberOfIntervals() - 1;

    char   type5 = classification[bgni];
    char   typem = 0;
    char   type3 = classification[endi];

    while ((bgni <= endi) && (type5 == classification[bgni]))
      bgni++;
    bgni--;

    while ((bgni <= endi) && (type3 == classification[endi]))
      endi--;
    endi++;

    delete[] classification;

    //  All the same classification?

    if (bgni == endi) {
      isLowCov = (type5 == 'l');
      isUnique = (type5 == 'u');
      isRepeat = (type5 == 'r');
    }

    //  Nope, if we aren't the same, assume it is uniqRepeat.

    else if (type5 != type3) {
      isUniqRepeat = true;
    }

    //  Nope, the same on both ends.  Assume we're just flipped.

    else {
      if (type5 == 'r')
        isUniqAnchor = true;
      else
        isSpanRepeat = true;
    }

    //  Now, do something with it.

    //  LOG - readID readLen classification

    if (isLowCov) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "low-cov");
      readLowCov->add(readLen);

      for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
        covrLowCov->add(depth.depth(ii), depth.hi(ii) - depth.lo(ii));
    }

    if (isUnique) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "unique");
      readUnique->add(readLen);

      for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
        covrUnique->add(depth.depth(ii), depth.hi(ii) - depth.lo(ii));
    }

    if ((isRepeat) && (readContained == true)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "contained-repeat");
      readRepeatCont->add(readLen);

      for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
        covrRepeatCont->add(depth.depth(ii), depth.hi(ii) - depth.lo(ii));
    }

    if ((isRepeat) && (readContained == false)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "dovetail-repeat");
      readRepeatDove->add(readLen);

      for (uint32 ii=0; ii<depth.numberOfIntervals(); ii++)
        covrRepeatDove->add(depth.depth(ii), depth.hi(ii) - depth.lo(ii));
    }

    if (isSpanRepeat) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "span-repeat");
      readSpanRepeat->add(readLen);
      olapSpanRepeat->add(depth.lo(endi) - depth.hi(bgni));
    }

    if ((isUniqRepeat) && (readContained == true)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "uniq-repeat-cont");
      readUniqRepeatCont->add(readLen);
    }

    if ((isUniqRepeat) && (readContained == false)) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "uniq-repeat-dove");
      readUniqRepeatDove->add(readLen);
    }

    if (isUniqAnchor) {
      fprintf(LOG, "%u\t%u\t%s\n", readID, readLen, "uniq-anchor");
      readUniqAnchor->add(readLen);
      olapUniqAnchor->add(depth.lo(endi) - depth.hi(bgni));
    }

    //  Done.  Read more data.

    overlapsLen = ovlStore->readOverlaps(overlaps, overlapsMax);
  }

  fclose(LOG);  //  Done with logging.

  readHole->finalizeData();
  olapHole->finalizeData();

  readHump->finalizeData();
  olapHump->finalizeData();

  readNo5->finalizeData();
  olapNo5->finalizeData();

  readNo3->finalizeData();
  olapNo3->finalizeData();


  readLowCov->finalizeData();
  olapLowCov->finalizeData();
  covrLowCov->finalizeData();

  readUnique->finalizeData();
  olapUnique->finalizeData();
  covrUnique->finalizeData();

  readRepeatCont->finalizeData();
  olapRepeatCont->finalizeData();
  covrRepeatCont->finalizeData();

  readRepeatDove->finalizeData();
  olapRepeatDove->finalizeData();
  covrRepeatDove->finalizeData();


  readSpanRepeat->finalizeData();
  olapSpanRepeat->finalizeData();

  readUniqRepeatCont->finalizeData();
  olapUniqRepeatCont->finalizeData();

  readUniqRepeatDove->finalizeData();
  olapUniqRepeatDove->finalizeData();

  readUniqAnchor->finalizeData();
  olapUniqAnchor->finalizeData();


  LOG = stdout;

  if (toFile == true) {
    sprintf(N, "%s.summary", outPrefix);

    LOG = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
  }

  fprintf(LOG, "category            reads     %%          read length        feature size or coverage  analysis\n");
  fprintf(LOG, "----------------  -------  -------  ----------------------  ------------------------  --------------------\n");
  fprintf(LOG, "middle-missing    %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (bad trimming)\n", readHole->numberOfObjects(), (float)readHole->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readHole->mean(), readHole->stddev(), olapHole->mean(), olapHole->stddev());
  fprintf(LOG, "middle-hump       %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (bad trimming)\n", readHump->numberOfObjects(), (float)readHump->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readHump->mean(), readHump->stddev(), olapHump->mean(), olapHump->stddev());
  fprintf(LOG, "no-5-prime        %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (bad trimming)\n", readNo5->numberOfObjects(), (float)readNo5->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readNo5->mean(), readNo5->stddev(), olapNo5->mean(), olapNo5->stddev());
  fprintf(LOG, "no-3-prime        %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (bad trimming)\n", readNo3->numberOfObjects(), (float)readNo3->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readNo3->mean(), readNo3->stddev(), olapNo3->mean(), olapNo3->stddev());
  fprintf(LOG, "\n");
  fprintf(LOG, "low-coverage      %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (easy to assemble, potential for lower quality consensus)\n", readLowCov->numberOfObjects(), (float)readLowCov->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readLowCov->mean(), readLowCov->stddev(), covrLowCov->mean(), covrLowCov->stddev());
  fprintf(LOG, "unique            %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (easy to assemble, perfect, yay)\n", readUnique->numberOfObjects(), (float)readUnique->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readUnique->mean(), readUnique->stddev(), covrUnique->mean(), covrUnique->stddev());
  fprintf(LOG, "repeat-cont       %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (potential for consensus errors, no impact on assembly)\n", readRepeatCont->numberOfObjects(), (float)readRepeatCont->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readRepeatCont->mean(), readRepeatCont->stddev(), covrRepeatCont->mean(), covrRepeatCont->stddev());
  fprintf(LOG, "repeat-dove       %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (hard to assemble, likely won't assemble correctly or even at all)\n", readRepeatDove->numberOfObjects(), (float)readRepeatDove->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readRepeatDove->mean(), readRepeatDove->stddev(), covrRepeatDove->mean(), covrRepeatDove->stddev());
  fprintf(LOG, "\n");
  fprintf(LOG, "span-repeat       %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (read spans a large repeat, usually easy to assemble)\n", readSpanRepeat->numberOfObjects(), (float)readSpanRepeat->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readSpanRepeat->mean(), readSpanRepeat->stddev(), olapSpanRepeat->mean(), olapSpanRepeat->stddev());
  fprintf(LOG, "uniq-repeat-cont  %7"F_U64P"  %6.2f  %10.2f +- %-8.2f                            (should be uniquely placed, low potential for consensus errors, no impact on assembly)\n", readUniqRepeatCont->numberOfObjects(), (float)readUniqRepeatCont->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readUniqRepeatCont->mean(), readUniqRepeatCont->stddev());
  fprintf(LOG, "uniq-repeat-dove  %7"F_U64P"  %6.2f  %10.2f +- %-8.2f                            (will end contigs, potential to misassemble)\n", readUniqRepeatDove->numberOfObjects(), (float)readUniqRepeatDove->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readUniqRepeatDove->mean(), readUniqRepeatDove->stddev());
  fprintf(LOG, "uniq-anchor       %7"F_U64P"  %6.2f  %10.2f +- %-8.2f   %10.2f +- %-8.2f   (repeat read, with unique section, probable bad read)\n", readUniqAnchor->numberOfObjects(), (float)readUniqAnchor->numberOfObjects()/gkpStore->gkStore_getNumReads()*100, readUniqAnchor->mean(), readUniqAnchor->stddev(), olapUniqAnchor->mean(), olapUniqAnchor->stddev());

  if (toFile == true)
    fclose(LOG);

  // Clean up the histograms
  delete readNoOlaps;
  delete readHole;
  delete readHump;
  delete readNo5;
  delete readNo3;

  delete olapHole;
  delete olapHump;
  delete olapNo5;
  delete olapNo3;

  delete readLowCov;
  delete readUnique;
  delete readRepeatCont;
  delete readRepeatDove;
  delete readSpanRepeat;
  delete readUniqRepeatCont;
  delete readUniqRepeatDove;
  delete readUniqAnchor;

  delete covrLowCov;
  delete covrUnique;
  delete covrRepeatCont;
  delete covrRepeatDove;
  delete covrSpanRepeat;
  delete covrUniqRepeatCont;
  delete covrUniqRepeatDove;
  delete covrUniqAnchor;

  delete olapLowCov;
  delete olapUnique;
  delete olapRepeatCont;
  delete olapRepeatDove;
  delete olapSpanRepeat;
  delete olapUniqRepeatCont;
  delete olapUniqRepeatDove;
  delete olapUniqAnchor;

  delete ovlStore;

  gkpStore->gkStore_close();

  exit(0);
}
