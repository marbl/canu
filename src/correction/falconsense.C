
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
                        gkStore           *gkpStore,
                        tgTig             *tig,
                        bool               trimToAlign,
                        gkReadData        *readData,
                        uint32             minOutputLength) {

  //  Grab and save the raw read for the template.

  fprintf(stderr, "Processing read %u of length %u with %u evidence reads.\n",
          tig->tigID(), tig->length(), tig->numberOfChildren());

  gkpStore->gkStore_loadReadData(tig->tigID(), readData);

  //  Now parse the layout and push all the sequences onto our seqs vector.

  falconInput   *evidence = new falconInput [tig->numberOfChildren() + 1];

  evidence[0].addInput(tig->tigID(),
                       readData->gkReadData_getRawSequence(),
                       readData->gkReadData_getRead()->gkRead_rawLength(),
                       0,
                       readData->gkReadData_getRead()->gkRead_rawLength());


  for (uint32 cc=0; cc<tig->numberOfChildren(); cc++) {
    tgPosition  *child = tig->getChild(cc);

    gkpStore->gkStore_loadReadData(child->ident(), readData);

    if (child->isReverse())
      reverseComplementSequence(readData->gkReadData_getRawSequence(),
                                readData->gkReadData_getRead()->gkRead_rawLength());

    //  For debugging/testing, skip one orientation of overlap.
    //
    //if (child->isReverse() == false)
    //  continue;
    //if (child->isReverse() == true)
    //  continue;

    //  Trim the read to the aligned bit
    char   *seq    = readData->gkReadData_getRawSequence();
    uint32  seqLen = readData->gkReadData_getRead()->gkRead_rawLength();

    if (trimToAlign) {
      seq    += child->askip();
      seqLen -= child->askip() + child->bskip();

      seq[seqLen] = 0;
    }

    //  Used to skip if read length was less or equal to min_ovl_len

    if (seqLen < 500) {
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

  falconData  *fd = fc->generateConsensus(evidence,
                                          tig->numberOfChildren() + 1);

  //  Find the largest stretch of uppercase sequence.  Lowercase sequence denotes MSA coverage was below minAllowedCoverage.

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
  char             *gkpName   = 0L;
  char             *corName   = 0L;
  uint32            corVers   = 1;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;

  uint32            idMin = 1;
  uint32            idMax = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            numThreads         = 1;
  uint32            minAllowedCoverage = 4;
  double            minIdentity        = 0.5;
  uint32            minOutputLength    = 500;

  bool              trimToAlign        = true;
  bool              restrictToOverlap  = true;

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


    } else if (strcmp(argv[arg], "-f") == 0) {   //  ALGORITHM OPTIONS
      restrictToOverlap = false;


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
      minIdentity = atof(argv[++arg]);


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
    fprintf(stderr, "ALGORITHM PARAMETERS\n");
    fprintf(stderr, "  -f               align evidence to the full read, ignore overlap position\n");
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

  //  Open inputs.

  gkStore  *gkpStore = gkStore::gkStore_open(gkpName);
  tgStore  *corStore = new tgStore(corName, corVers);

  uint32    numReads = gkpStore->gkStore_getNumReads();

  //  Decide what reads to operate on.

  if (numReads < idMax)
    idMax = numReads;

  loadReadList(readListName, idMin, idMax, readList);

  //  Open logging and summary files

  FILE *logFile = AS_UTL_openOutputFile(outputPrefix, '.', "log",   false);   //  Not used.
  FILE *cnsFile = AS_UTL_openOutputFile(outputPrefix, '.', "cns",   true);
  FILE *seqFile = AS_UTL_openOutputFile(outputPrefix, '.', "fasta", false);   //  Not useful.

  //  Initialize processing.

  falconConsensus   *fc = new falconConsensus(minAllowedCoverage, minIdentity, minOutputLength, restrictToOverlap);
  gkReadData        *rd = new gkReadData;

  //  And process.

  for (uint32 ii=idMin; ii<idMax; ii++) {
    if ((readList.size() > 0) &&                     //  Skip reads not on the read list.  We need
        (readList.count(ii) == 0))
      continue;

    tgTig *layout = corStore->loadTig(ii);

    generateFalconConsensus(fc, gkpStore, layout, trimToAlign, rd, minOutputLength);

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

  gkpStore->gkStore_close();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
