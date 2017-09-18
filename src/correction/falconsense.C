
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

#include "AS_UTL_decodeRange.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_UTL_fasta.H"

#include "falconConsensus.H"

#include <set>

using namespace std;




//  Duplicated in generateCorrectionLayouts.C
void
loadReadList(char *readListName, uint32 iidMin, uint32 iidMax, set<uint32> &readList) {
  char  L[1024];

  if (readListName == NULL)
    return;

  errno = 0;
  FILE *R = fopen(readListName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", readListName, strerror(errno)), exit(1);

  fgets(L, 1024, R);
  while (!feof(R)) {
    splitToWords W(L);
    uint32       id = W(0);

    if ((iidMin <= id) &&
        (id     <= iidMax))
      readList.insert(W(0));

    fgets(L, 1024, R);
  }

  fclose(R);
}





//  A mash up of falcon_sense.C and outputFalcon.C

void
generateFalconConsensus(falconConsensus   *fc,
                        gkStore           *gkpStore,
                        tgTig             *tig,
                        bool               trimToAlign,
                        FILE              *F,
                        gkReadData        *readData,
                        uint32             minOutputLength) {

  //  Grab and save the raw read for the template.

  fprintf(stderr, "Processing read %u of length %u with %u evidence reads.\n",
          tig->tigID(), tig->length(), tig->numberOfChildren());

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  //  Now parse the layout and push all the sequences onto our seqs vector.

  falconInput   *evidence = new falconInput [tig->numberOfChildren() + 1];

  evidence[0].addInput(tig->tigID(),
                       readData->gkReadData_getSequence(),
                       readData->gkReadData_getRead()->gkRead_sequenceLength(),
                       0,
                       readData->gkReadData_getRead()->gkRead_sequenceLength());


  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getSequence(),
                                readData->gkReadData_getRead()->gkRead_sequenceLength());

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq    = readData->gkReadData_getSequence();
    uint32  seqLen = readData->gkReadData_getRead()->gkRead_sequenceLength();

    if (trimToAlign) {
      seq    += child->askip();
      seqLen -= child->askip() + child->bskip();

      seq[seqLen] = 0;
    }

    //  Used to skip if read length was less or equal to min_ovl_len

    evidence[cc+1].addInput(child->ident(), seq, seqLen, child->min(), child->max());
  }

  //  Loaded all reads, build consensus.

  uint32 splitSeqID = 0;

  //FConsensus::consensus_data *consensus_data_ptr = FConsensus::generate_consensus( seqs, min_cov, min_idt, min_ovl_len, max_read_len );

  falconData  *fd = fc->generateConsensus(evidence,
                                          tig->numberOfChildren() + 1);

#ifdef TRACK_POSITIONS
  //const std::string& sequenceToCorrect = seqs.at(0);
  char * originalStringPointer = consensus_data_ptr->sequence;
#endif

  char * split = strtok(fd->seq, "acgt");

  while (split != NULL) {
    if (strlen(split) > minOutputLength) {
      fprintf(stderr, "Generated read %u_%u of length %lu.\n",
              tig->tigID(), splitSeqID, strlen(split));

      AS_UTL_writeFastA(F, split, strlen(split), 60, ">read%u_%d\n", tig->tigID(), splitSeqID);

      splitSeqID++;

#ifdef TRACK_POSITIONS
      int distance_from_beginning = split - originalStringPointer;
      std::vector<int> relevantOriginalPositions(consensus_data_ptr->originalPos.begin() + distance_from_beginning, consensus_data_ptr->originalPos.begin() + distance_from_beginning + strlen(split));
      int firstRelevantPosition = relevantOriginalPositions.front();
      int lastRelevantPosition = relevantOriginalPositions.back();

      std::string relevantOriginalTemplate = seqs.at(0).substr(firstRelevantPosition, lastRelevantPosition - firstRelevantPosition + 1);

      // store relevantOriginalTemplate along with corrected read - not implemented
#endif
    }

    split = strtok(NULL, "acgt");
  }

  delete fd;  //FConsensus::free_consensus_data( consensus_data_ptr );
  delete [] evidence;
}




int
main(int argc, char **argv) {
  char             *gkpName   = 0L;
  char             *corName   = 0L;
  uint32            corVers   = 1;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;

  char              logName[FILENAME_MAX] = {0};
  char              sumName[FILENAME_MAX] = {0};
  char              flgName[FILENAME_MAX] = {0};
  FILE             *logFile = 0L;

  uint32            idMin = 0;
  uint32            idMax = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            numThreads         = 1;
  uint32            minAllowedCoverage = 4;
  double            minIdentity        = 0.5;
  uint32            minOutputLength    = 500;

  bool              trimToAlign        = true;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {   //  INPUTS
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-t") == 0) {   //  COMPUTE RESOURCES
      numThreads = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      idMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      idMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      readListName = argv[++arg];


    } else if (strcmp(argv[arg], "-cc") == 0) {   //  CONSENSUS
      minAllowedCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-ci") == 0) {
      minIdentity = atoi(argv[++arg]);


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
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "  -G gkpStore      mandatory path to gkpStore\n");
    fprintf(stderr, "  -C corStore      mandatory path to corStore\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESOURCE PARAMETERS\n");
    fprintf(stderr, "  -t numThreads    number of compute threads to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS\n");
    fprintf(stderr, "  -cc coverage     minimum consensus coverage to output corrected base\n");
    fprintf(stderr, "  -cl length       minimum length of corrected region\n");
    fprintf(stderr, "  -ci identity     minimum identity of an aligned evidence read\n");
    fprintf(stderr, "\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no gkpStore input (-G) supplied.\n");
    if (corName == NULL)
      fprintf(stderr, "ERROR: no corStore input (-C) supplied.\n");

    exit(1);
  }


  omp_set_num_threads(numThreads);



  //  Open inputs and output tigStore.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore  *corStore = new tgStore(corName, corVers);

  uint32    numReads = gkpStore->gkStore_getNumReads();

  //  Decide what reads to operate on.

  if (numReads < idMax)
    idMax = numReads;

  loadReadList(readListName, idMin, idMax, readList);

  //  Open logging and summary files

  logFile = AS_UTL_openOutputFile(outputPrefix, "log");

  //  Initialize processing.

  falconConsensus   *fc = new falconConsensus(minAllowedCoverage, minIdentity, minOutputLength);
  gkReadData        *rd = new gkReadData;

  //  And process.

  for (uint32 ii=idMin; ii<idMax; ii++) {
    if ((readList.size() > 0) &&                     //  Skip reads not on the read list.  We need
        (readList.count(ii) == 0))
      continue;

    tgTig *layout = corStore->loadTig(ii);

    generateFalconConsensus(fc, gkpStore, layout, trimToAlign, stdout, rd, minOutputLength);

    corStore->unloadTig(ii);
  }

  //  Close files and clean up.

  if (logFile != NULL)   fclose(logFile);

  delete    fc;
  delete    rd;
  delete    corStore;

  gkpStore->gkStore_close();

  return(0);
}
