
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
 *    Brian P. Walenz beginning on 2017-SEP-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "sqStore.H"
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

  //  To give the list some size - if the read list has no elements between iidMin and iidMax,
  //  nothing is inserted below, and then we _think_ no read list was supplied, and try to
  //  process every read, when in fact we should be processing no reads.

  readList.insert(0);

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

  AS_UTL_closeFile(R, readListName);
}





//  A mash up of falcon_sense.C and outputFalcon.C

void
generateFalconConsensus(falconConsensus   *fc,
                        sqStore           *seqStore,
                        tgTig             *tig,
                        bool               trimToAlign,
                        sqReadData        *readData,
                        uint32             minOlapLength) {

  //  Grab and save the raw read for the template.

  fprintf(stderr, "Processing read %u of length %u with %u evidence reads.\n",
          tig->tigID(), tig->length(), tig->numberOfChildren());

  seqStore->sqStore_loadReadData(tig->tigID(), readData);

  //  Now parse the layout and push all the sequences onto our seqs vector.

  falconInput   *evidence = new falconInput [tig->numberOfChildren() + 1];

  evidence[0].addInput(tig->tigID(),
                       readData->sqReadData_getRawSequence(),
                       readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw),
                       0,
                       readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw));


  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    seqStore->sqStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->sqReadData_getRawSequence(),
                                readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw));

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq    = readData->sqReadData_getRawSequence();
    uint32  seqLen = readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw);

    if (trimToAlign) {
      seq    += child->askip();
      seqLen -= child->askip() + child->bskip();

      seq[seqLen] = 0;
    }

    //  Ignore the read if it's less than the minimum overlap length - it'll have no chance to align correctly anyway!

    if (seqLen < minOlapLength) {
#ifdef BRI
      fprintf(stderr, "read %7u loaded %6u-%6u to template %6u-%6u -- TOO SHORT\n",
              child->ident(),
              child->askip(),
              child->askip() + seqLen,
              child->min(),
              child->max());
#endif
      continue;
    }

#ifdef BRI
    fprintf(stderr, "read %7u loaded %6u-%6u to template %6u-%6u\n",
            child->ident(),
            child->askip(),
            child->askip() + seqLen,
            child->min(),
            child->max());
#endif

    evidence[cc+1].addInput(child->ident(), seq, seqLen, child->min(), child->max());
  }

  //  Loaded all reads, build consensus.

  falconData  *fd = fc->generateConsensus(evidence, tig->numberOfChildren() + 1);

  //  Find the largest stretch of uppercase sequence.  Lowercase sequence denotes MSA coverage was below minOutputCoverage.

  uint32  bgn = 0;
  uint32  end = 0;

  for (uint32 in=0, bb=0, ee=0; ee<fd->len; ee++) {
    bool   isLower = (('a' <= fd->seq[ee]) && (fd->seq[ee] <= 'z'));
    bool   isLast  = (ee == fd->len - 1);

    if ((in == true) && (isLower || isLast))     //  Report the regions we could be saving.
      fprintf(stderr, "Read %u region %u-%u\n",
              tig->tigID(), bb, ee + isLast);

    if (isLower) {                                 //  If lowercase, declare that we're not in a
      in = 0;                                      //  good region any more.
    }

    else if (in == 0) {                            //  Otherwise, if not in a region (so the first
      bb = ee;                                     //  uppercase), remember the coordinate and
      in = 1;                                      //  switch to being 'in' a region.
    }

    if ((in == 1) && (ee + 1 - bb > end - bgn)) {  //  If 'in' a good region, remember the longest.
      bgn = bb;                                    //  'ee + 1': if the next letter is lower case
      end = ee + 1;                                //  our bgn,end interval will be set in this iteration
    }
  }

  //  Update the layout with consensus sequence, positions, et cetera.
  //  If the whole string is lowercase (grrrr!) then bgn == end == 0.

  tig->_sourceID    = tig->tigID();
  tig->_sourceBgn   = (end == 0) ? (0) : (fd->pos[bgn]);            //  Space based (probably).
  tig->_sourceEnd   = (end == 0) ? (0) : (fd->pos[end - 1] + 1);    //  Space based.

  resizeArrayPair(tig->_gappedBases, tig->_gappedQuals, tig->_gappedLen, tig->_gappedMax, end - bgn + 1, resizeArray_doNothing);

  for (uint32 ii=bgn; ii<end; ii++) {
    tig->_gappedBases[ii-bgn] = fd->seq[ii];
    tig->_gappedQuals[ii-bgn] = fd->eqv[ii];
  }

  tig->_gappedLen = end - bgn;

  tig->_gappedBases[tig->_gappedLen] = 0;
  tig->_gappedQuals[tig->_gappedLen] = 0;

#ifdef BRIOUT
  if (tig->_gappedBases[0] != 0)
    fprintf(stdout, "read%u %s\n", tig->tigID(), tig->_gappedBases);
#endif

  delete fd;
  delete [] evidence;
}




int
main(int argc, char **argv) {
  char             *seqName   = 0L;
  char             *corName   = 0L;
  uint32            corVers   = 1;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;

  uint32            idMin = 1;
  uint32            idMax = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            numThreads         = 1;

  uint32            minOutputCoverage  = 4;
  uint32            minOutputLength    = 1000;
  double            minOlapIdentity    = 0.5;
  double            minOlapLength      = 500;

  bool              trimToAlign        = true;
  bool              restrictToOverlap  = true;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUTS
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-t") == 0) {   //  COMPUTE RESOURCES
      numThreads = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-f") == 0) {   //  ALGORITHM OPTIONS
      restrictToOverlap = false;


    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      idMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      idMax = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      readListName = argv[++arg];


    } else if (strcmp(argv[arg], "-cc") == 0) {   //  CONSENSUS
      minOutputCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-oi") == 0) {
      minOlapIdentity = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minOlapLength = atof(argv[++arg]);


    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (seqName == NULL)
    err++;
  if (corName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S seqStore      mandatory path to seqStore\n");
    fprintf(stderr, "  -C corStore      mandatory path to corStore\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESOURCE PARAMETERS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t numThreads    number of compute threads to use\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ALGORITHM PARAMETERS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -f               align evidence to the full read, ignore overlap position\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Of outputs:\n");
    fprintf(stderr, "  -cc coverage     minimum consensus coverage to generate corrected base\n");
    fprintf(stderr, "  -cl length       minimum length of corrected region to output as a corrected read\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " Of inputs:\n");
    fprintf(stderr, "  -oi identity     minimum identity of an aligned evidence read overlap\n");
    fprintf(stderr, "  -ol length       minimum length of an aligned evidence read overlap\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR: no seqStore input (-S) supplied.\n");
    if (corName == NULL)
      fprintf(stderr, "ERROR: no corStore input (-C) supplied.\n");

    exit(1);
  }


  omp_set_num_threads(numThreads);

  //  Open inputs.

  sqRead_setDefaultVersion(sqRead_raw);

  sqStore  *seqStore = sqStore::sqStore_open(seqName);
  tgStore  *corStore = new tgStore(corName, corVers);

  uint32    numReads = seqStore->sqStore_getNumReads();

  //  Decide what reads to operate on.

  if (numReads < idMax)
    idMax = numReads;

  loadReadList(readListName, idMin, idMax, readList);

  //  Open logging and summary files

  FILE *logFile = AS_UTL_openOutputFile(outputPrefix, '.', "log",   false);   //  Not used.
  FILE *cnsFile = AS_UTL_openOutputFile(outputPrefix, '.', "cns",   true);
  FILE *seqFile = AS_UTL_openOutputFile(outputPrefix, '.', "fasta", false);   //  Not useful.

  //  Initialize processing.

  falconConsensus   *fc = new falconConsensus(minOutputCoverage, minOutputLength, minOlapIdentity, minOlapLength, restrictToOverlap);
  sqReadData        *rd = new sqReadData;

  //  And process.

  for (uint32 ii=idMin; ii<idMax; ii++) {
    if ((readList.size() > 0) &&                     //  Skip reads not on the read list.  We need
        (readList.count(ii) == 0))
      continue;

    tgTig *layout = corStore->loadTig(ii);

    generateFalconConsensus(fc, seqStore, layout, trimToAlign, rd, minOlapLength);

    if (cnsFile)
      layout->saveToStream(cnsFile);

    if (seqFile)
      layout->dumpFASTA(seqFile, false);

    corStore->unloadTig(ii);
  }

  //  Close files and clean up.

  AS_UTL_closeFile(logFile);
  AS_UTL_closeFile(cnsFile);
  AS_UTL_closeFile(seqFile);

  delete    fc;
  delete    rd;
  delete    corStore;

  seqStore->sqStore_close();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
