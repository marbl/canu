
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
 *    Brian P. Walenz from 2015-APR-09 to 2015-SEP-21
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-NOV-27
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "ovStore.H"
#include "tgStore.H"

#include "stashContains.H"

#include "splitToWords.H"
#include "intervalList.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_UTL_fasta.H"

#include <set>

using namespace std;

//  Debugging on which reads are filtered, which are used, and which are removed
//  to meet coverage thresholds.  Very big.
#undef DEBUG_LAYOUT


//  Duplicated in generateCorrectionLayouts.C
void
loadReadList(char *readListName, uint32 iidMin, uint32 iidMax, set<uint32> &readList) {
  char  L[1024];

  if (readListName == NULL)
    return;

  FILE *R = AS_UTL_openInputFile(readListName);

  for (fgets(L, 1024, R);
       feof(R) == false;
       fgets(L, 1024, R)) {
    splitToWords W(L);
    uint32       id = W(0);

    if ((iidMin <= id) &&
        (id     <= iidMax))
      readList.insert(W(0));
  }

  fclose(R);
}


void
generateLayout(
                        gkStore           *gkpStore,
                        tgTig             *tig,
                        bool               trimToAlign,
                        gkReadData        *readData,
                        uint32             minOutputLength, uint32 minOverlapLength) {

  //  Grab and save the raw read for the template.

  fprintf(stderr, "Processing read %u of length %u with %u evidence reads.\n",
          tig->tigID(), tig->length(), tig->numberOfChildren());

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  //  Now parse the layout and push all the sequences onto our seqs vector.
  if ( readData->gkReadData_getRead()->gkRead_sequenceLength() < minOutputLength) {
     return;
  }

  fprintf(stdout, "read%d %s\n", tig->tigID(), readData->gkReadData_getSequence());

  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getSequence(),
                                readData->gkReadData_getRead()->gkRead_sequenceLength());

    //  Trim the read to the aligned bit
    char   *seq    = readData->gkReadData_getSequence();
    uint32  seqLen = readData->gkReadData_getRead()->gkRead_sequenceLength();

    if (trimToAlign) {
      seq    += child->askip();
      seqLen -= child->askip() + child->bskip();

      seq[seqLen] = 0;
    }

    //  Used to skip if read length was less or equal to min_ovl_len
    if (seqLen < minOverlapLength) {
       continue;
    }

    fprintf(stdout, "%d %s\n", child->ident(), seq);
  }

  fprintf(stdout, "+ +\n");
}





int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *corName   = 0L;

  char             *scoreName = 0L;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;
  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  char              flgName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;
  FILE             *sumFile = 0L;

  uint32            expectedCoverage    = 40;    //  How many overlaps per read to save, global filter
  uint32            minEvidenceOverlap  = 40;
  uint32            minEvidenceCoverage = 4;

  uint32            iidMin       = 0;
  uint32            iidMax       = UINT32_MAX;

  uint32            minEvidenceLength   = 0;
  double            maxEvidenceErate    = 1.0;
  double            maxEvidenceCoverage = DBL_MAX;

  uint32            minCorLength        = 0;
  char             *readListName = NULL;
  set<uint32>       readList;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {   //  INPUTS
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {   //  OUTPUT FORMAT
      corName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      iidMin  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      iidMax  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-rl") == 0) {
      readListName = argv[++arg];

    } else if (strcmp(argv[arg], "-eL") == 0) {   //  EVIDENCE SELECTION
      minEvidenceLength  = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-eM") == 0) {
      minCorLength = atoi(argv[++arg]);


    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (corName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS\n");
    fprintf(stderr, "  -G gkpStore      mandatory path to gkpStore\n");
    fprintf(stderr, "  -S file          overlap score thresholds (from filterCorrectionOverlaps)\n");
    fprintf(stderr, "                     if not supplied, will be estimated from ovlStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS\n");
    fprintf(stderr, "  -C corStore      output layouts to store 'corStore'\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION\n");
    fprintf(stderr, "  -b bgnID         process reads starting at bgnID\n");
    fprintf(stderr, "  -e endID         process reads up to but not including endID\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EVIDENCE SELECTION\n");
    fprintf(stderr, "  -eL length       minimum length of evidence overlaps\n");
    fprintf(stderr, "  -eM length       minimum length of a corrected read\n");              //  not used in canu
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no input gkpStore (-G) supplied.\n");
    if (corName == NULL)
      fprintf(stderr, "ERROR: no output corStore (-C) supplied.\n");
    exit(1);
  }

  //  Open inputs and output tigStore.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore  *corStore = new tgStore(corName, 1);
  gkReadData        *rd = new gkReadData;

  uint32    numReads = gkpStore->gkStore_getNumReads();

  //  Threshold the range of reads to operate on.

  if (numReads < iidMin) {
    fprintf(stderr, "ERROR: only " F_U32 " reads in the store (IDs 0-" F_U32 " inclusive); can't process requested range -b " F_U32 " -e " F_U32 "\n",
            numReads,
            numReads-1,
            iidMin, iidMax);
    exit(1);
  }
  if (numReads < iidMax)
    iidMax = numReads;

  loadReadList(readListName, iidMin, iidMax, readList);

  for (uint32 ii=iidMin; ii<iidMax; ii++) {
    if ((readList.size() > 0) &&                     //  Skip reads not on the read list.  We need
        (readList.count(ii) == 0))
      continue;

    tgTig *layout = corStore->loadTig(ii);

    generateLayout(gkpStore, layout, true, rd, minCorLength, minEvidenceLength);

    corStore->unloadTig(ii);
  }
  fprintf(stdout, "- -\n");

  //  Close files and clean up.

  delete    corStore;

  gkpStore->gkStore_close();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
