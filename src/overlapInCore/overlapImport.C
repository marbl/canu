
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

#include "system.H"
#include "strings.H"
#include "math.H"

#include "sqStore.H"
#include "ovStore.H"

#include <vector>

//  Life is so much easier if these are just global.

char const               *seqStoreName     = nullptr;

char const               *ovlFileName      = nullptr;
char const               *ovlStoreName     = nullptr;

enum class inType {
  asCoords,
  asHangs,
  asUnaligned,
  asPAF,
  asOVB,
  asRandom
};

inType                    inputType        = inType::asCoords;

double                    maxError         = 1.0;
uint32                    minReadLength    = 0;
uint32                    minOverlapLength = 0;

uint64                    rmin = 0, rmax = 0;
uint32                    abgn = 1, aend = 0;
uint32                    bbgn = 1, bend = 0;

std::vector<char const *> infiles;

uint64                    totalOverlaps   = 0;
uint64                    filteredOlapLen = 0;
uint64                    filteredReadLen = 0;
uint64                    filteredBoth    = 0;
uint64                    overlapOutput   = 0;

sqStore                  *ss = nullptr;

ovOverlap                 ov;

ovFile                   *of = nullptr;
ovStoreWriter            *os = nullptr;



void
outputOverlap(void) {
  bool  shortOvl  = false;
  bool  shortRead = false;

  //  Test if we want to emit this overlap.

  if (ov.length() < minOverlapLength)
    shortOvl = true;

  if ((ss->sqStore_getReadLength(ov.a_iid) < minReadLength) ||
      (ss->sqStore_getReadLength(ov.b_iid) < minReadLength))
    shortRead = true;

  //  Count the test results.

  if ((shortOvl) && (shortRead))   filteredBoth++;
  if  (shortOvl)                   filteredOlapLen++;
  if                (shortRead)    filteredReadLen++;

  //  Bail if we want to discard the overlap.

  totalOverlaps++;

  if ((shortOvl) || (shortRead))
    return;

  //  Count and output the good overlap.

  overlapOutput++;

  if (of)   of->writeOverlap(&ov);
  if (os)   os->writeOverlap(&ov);
}



void
makeRandomOverlaps(void) {
  mtRandom  mt;

  if (aend == 0)   aend = ss->sqStore_lastReadID();
  if (bend == 0)   bend = ss->sqStore_lastReadID();

  uint64  numRandom = rmin + floor(mt.mtRandomRealOpen() * (rmax - rmin));

  for (uint64 ii=0; ii<numRandom; ii++) {
    uint32   aID      = abgn + floor(mt.mtRandomRealOpen() * (aend - abgn));
    uint32   bID      = bbgn + floor(mt.mtRandomRealOpen() * (bend - bbgn));

    if (aID < bID) {
      uint32 t = aID;
      aID = bID;
      bID = aID;
    }

#if 0
    //  For testing when reads have no overlaps in store building.  Issue #302.
    aID = aID & 0xfffffff0;
    bID = bID & 0xfffffff0;

    if (aID == 0)   aID = 1;
    if (bID == 0)   bID = 1;
#endif

    uint32   aLen     = ss->sqStore_getReadLength(aID);
    uint32   bLen     = ss->sqStore_getReadLength(bID);

    bool     olapFlip = mt.mtRandom32() % 2;

    //  We could be fancy and make actual overlaps that make sense, or punt and make overlaps that
    //  are valid but nonsense.

    ov.a_iid = aID;
    ov.b_iid = bID;

    ov.flipped(olapFlip);

    ov.a_hang((int32)(mt.mtRandomRealOpen() * 2 * aLen - aLen));
    ov.b_hang((int32)(mt.mtRandomRealOpen() * 2 * bLen - bLen));

    ov.dat.ovl.forOBT = false;
    ov.dat.ovl.forDUP = false;
    ov.dat.ovl.forUTG = true;

    ov.erate(mt.mtRandomRealOpen() * 0.1);

    outputOverlap();
  }
}



void
readOVB(ovFile *in) {

  while (in->readOverlap(&ov)) {
    outputOverlap();
  }

  delete in;
}



void
readASCII(compressedFileReader *in) {
  char           S[1024];
  splitToWords   W;

  fgets(S, 1024, in->file());

  while (!feof(in->file())) {
    W.split(S);

    switch (inputType) {
      case inType::asCoords:
        ov.fromString(W, ovOverlapAsCoords);
        break;
      case inType::asHangs:
        ov.fromString(W, ovOverlapAsHangs);
        break;
      case inType::asUnaligned:
        ov.fromString(W, ovOverlapAsUnaligned);
        break;
      case inType::asPAF:
        ov.fromString(W, ovOverlapAsPaf);
        break;
      default:
        assert(0);
        break;
    }

    outputOverlap();

    fgets(S, 1024, in->file());
  }

  delete in;
}



int
main(int argc, char **argv) {

  argc = AS_configure(argc, argv, 1);

  std::vector<char const *> err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-S") == 0)
      seqStoreName = argv[++arg];
    else if (strcmp(argv[arg], "-o") == 0)
      ovlFileName = argv[++arg];
    else if (strcmp(argv[arg], "-O") == 0)
      ovlStoreName = argv[++arg];

    else if (strcmp(argv[arg], "-raw") == 0)
      sqRead_setDefaultVersion(sqRead_raw);
    else if (strcmp(argv[arg], "-obt") == 0)
      sqRead_setDefaultVersion(sqRead_corrected);
    else if (strcmp(argv[arg], "-utg") == 0)
      sqRead_setDefaultVersion(sqRead_trimmed);


    else if (strcmp(argv[arg], "-coords") == 0)
      inputType = inType::asCoords;
    else if (strcmp(argv[arg], "-hangs") == 0)
      inputType = inType::asHangs;
    else if (strcmp(argv[arg], "-unaligned") == 0)
      inputType = inType::asUnaligned;
    else if (strcmp(argv[arg], "-paf") == 0)
      inputType = inType::asPAF;
    else if (strcmp(argv[arg], "-ovb") == 0)
      inputType = inType::asOVB;

    else if (strcmp(argv[arg], "-maxerror") == 0) {
      if (*strtonumber(argv[++arg], maxError) == '%')
        maxError /= 100.0;
    }

    else if (strcmp(argv[arg], "-minreadlength") == 0)
      minReadLength = strtouint32(argv[++arg]);

    else if (strcmp(argv[arg], "-minoverlaplength") == 0)
      minOverlapLength = strtouint32(argv[++arg]);

    else if (strcmp(argv[arg], "-random") == 0) {
      inputType = inType::asRandom;
      decodeRange(argv[++arg], rmin, rmax);
    }

    else if (strcmp(argv[arg], "-a") == 0)
      decodeRange(argv[++arg], abgn, aend);
    else if (strcmp(argv[arg], "-b") == 0)
      decodeRange(argv[++arg], bbgn, bend);

    else if ((strcmp(argv[arg], "-") == 0) ||
             (fileExists(argv[arg])))
      infiles.push_back(argv[arg]);

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "%s: unknown option '%s'.\n", argv[0], argv[arg]);
      err.push_back(s);
    }
  }


  if (seqStoreName == nullptr)
    err.push_back("ERROR: no input seqStore (-S) supplied.\n");

  if ((infiles.size() == 0) && (inputType != inType::asRandom))
    err.push_back("ERROR: no input overlap files supplied.\n");


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [options] ascii-ovl-file-input.[.gz]\n", argv[0]);
    fprintf(stderr, "  -S name.seqStore    path to valid sequence store\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT FORMAT (pick exactly one):\n");
    fprintf(stderr, "  -o file.ovb         output file name\n");
    fprintf(stderr, "  -O name.ovlStore    output overlap store");
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUT FORMAT:\n");
    fprintf(stderr, "  -coords             as coordinates on each read (default)\n");
    fprintf(stderr, "  -hangs              as dovetail hangs\n");
    fprintf(stderr, "  -unaligned          as unaligned regions on each read\n");
    fprintf(stderr, "  -paf                as miniasm Pairwise mApping Format\n");
    fprintf(stderr, "  -ovb                as canu binary .ovb overlaps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "FILTER OPTIONS:\n");
    fprintf(stderr, "  -maxerror x         discard overlaps above x fraction error (e.g., '0.05').\n");
    fprintf(stderr, "  -maxerror x%%        discard overlaps above x percent error (e.g., '5.0%%').\n");
    fprintf(stderr, "  -minreadlength l    discard overlaps involving reads shorter than l bases.\n");
    fprintf(stderr, "  -minoverlaplength l discard overlaps shorter than l bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ VERSION:\n");
    fprintf(stderr, "  -raw                uncorrected raw reads\n");
    fprintf(stderr, "  -obt                corrected reads\n");
    fprintf(stderr, "  -utg                trimmed reads\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RANDOM OPTIONS:\n");
    fprintf(stderr, "  Doesn't read overlaps from an input file, instead generated\n");
    fprintf(stderr, "  random non-sense overlaps.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -random N           create N random overlaps, for store testing\n");
    fprintf(stderr, "  -a x-y              A read IDs will be between x and y\n");
    fprintf(stderr, "  -b x-y              B read IDs will be between x and y\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input file can be stdin ('-') or a gz/bz2/xz compressed file.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  ss = new sqStore(seqStoreName);

  of = (ovlFileName  == nullptr) ? nullptr : new ovFile(ss, ovlFileName, ovFileFullWrite);
  os = (ovlStoreName == nullptr) ? nullptr : new ovStoreWriter(ovlStoreName, ss);

  //  Make random inputs first.

  switch (inputType) {
    case inType::asCoords:
    case inType::asHangs:
    case inType::asUnaligned:
    case inType::asPAF:
      for (uint32 ff=0; ff<infiles.size(); ff++)
        readASCII(new compressedFileReader(infiles[ff]));
      break;

    case inType::asOVB:
      for (uint32 ff=0; ff<infiles.size(); ff++)
        readOVB(new ovFile(ss, infiles[ff], ovFileFull));
      break;

    case inType::asRandom:
      makeRandomOverlaps();
      break;

    default:
      assert(0);
      break;
  }

  //  Close outputs and inputs.

  delete os;
  delete of;
  delete ss;

  //  Bye.

  filteredOlapLen -= filteredBoth;
  filteredReadLen -= filteredBoth;

  fprintf(stderr, "Overlaps processed:\n");
  fprintf(stderr, "  %8lu\n", totalOverlaps);
  fprintf(stderr, "Overlaps discarded:\n");
  fprintf(stderr, "  %8lu %6.2f%% - overlap length  too short\n", filteredOlapLen, 100.0 * filteredOlapLen / totalOverlaps);
  fprintf(stderr, "  %8lu %6.2f%% - read    length  too short\n", filteredReadLen, 100.0 * filteredReadLen / totalOverlaps);
  fprintf(stderr, "  %8lu %6.2f%% - both    lengths too short\n", filteredBoth,    100.0 * filteredBoth    / totalOverlaps);
  fprintf(stderr, "Overlaps output:\n");
  fprintf(stderr, "  %8lu %6.2f%%\n", overlapOutput, 100.0 * overlapOutput / totalOverlaps);
  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}
