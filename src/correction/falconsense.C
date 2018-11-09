
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
 *    Sergey Koren beginning on 2018-SEP-10
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

#include "intervalList.H"

#include "sequence.H"

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
    uint32  id = strtouint32(L);

    if ((iidMin <= id) &&
        (id     <= iidMax))
      readList.insert(id);
  }

  AS_UTL_closeFile(R, readListName);
}




sqReadData *
loadReadData(uint32                     readID,
             sqStore                   *seqStore,
             map<uint32, sqRead *>     &reads,
             map<uint32, sqReadData *> &datas) {

  if (datas.count(readID) == 0) {
    datas[readID] = new sqReadData;
    seqStore->sqStore_loadReadData(readID, datas[readID]);
  }

  assert(datas[readID] != NULL);

  return(datas[readID]);
}



void
generateFalconConsensus(falconConsensus           *fc,
                        tgTig                     *layout,
                        sqStore                   *seqStore,
                        map<uint32, sqRead *>     &reads,
                        map<uint32, sqReadData *> &datas,
                        bool                       trimToAlign,
                        uint32                     minOlapLength) {

  //  What rolls down stairs
  //  alone or in pairs,
  //  rolls over your neighbor's dog?
  //  What's great for a snack,
  //  And fits on your back?
  //  It's log, log, log!

  fprintf(stdout, "%8u %7u %8u", layout->tigID(), layout->length(), layout->numberOfChildren());

  //  Parse the layout and push all the sequences onto our seqs vector.  The first 'evidence'
  //  sequence is the read we're trying to correct.

  falconInput   *evidence = new falconInput [layout->numberOfChildren() + 1];
  sqReadData    *readData = loadReadData(layout->tigID(), seqStore, reads, datas);

  evidence[0].addInput(layout->tigID(),
                       readData->sqReadData_getRawSequence(),
                       readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw),
                       0,
                       readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw));

  for (uint32 cc=0; cc<layout->numberOfChildren(); cc++) {
    tgPosition  *child = layout->getChild(cc);

    //  Grab the read data.

    readData = loadReadData(child->ident(), seqStore, reads, datas);

    //  Make a copy of the sequence.  Don't modify the original sequence data because it's potentially cached now.

    char    *seq    = duplicateString(readData->sqReadData_getRawSequence());
    uint32   seqLen = readData->sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw);

    //  Now screw up the sequence by reverse-complementing and trimming it.

    if (child->isReverse())
      reverseComplementSequence(seq, seqLen);

    uint32  b = 0;
    uint32  e = seqLen;

    if (trimToAlign) {
      b += child->askip();
      e -= child->bskip();
    }

    seq[e] = 0;

    //  Save the read if it is larger than the minimum overlap length.  Anything smaller than this will have zero chance of aligning.

    if (minOlapLength <= e - b)
      evidence[cc+1].addInput(child->ident(), seq + b, e - b, child->min(), child->max());

    delete [] seq;
  }

  //  Loaded all reads, build consensus.

  falconData  *fd = fc->generateConsensus(evidence, layout->numberOfChildren() + 1);

  //  Find the largest stretch of uppercase sequence.  Lowercase sequence denotes MSA coverage was below minOutputCoverage.

  uint32  bgn = 0;
  uint32  end = 0;

  for (uint32 in=0, bb=0, ee=0; ee<fd->len; ee++) {
    bool   isLower = (('a' <= fd->seq[ee]) && (fd->seq[ee] <= 'z'));
    bool   isLast  = (ee == fd->len - 1);

    if ((in == true) && (isLower || isLast))       //  Report the regions we could be saving.
      fprintf(stdout, " %6u-%-6u", bb, ee + isLast);

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

  fprintf(stdout, "\n");

  //  Update the layout with consensus sequence, positions, et cetera.
  //  If the whole string is lowercase (grrrr!) then bgn == end == 0.

  layout->_sourceID    = layout->tigID();
  layout->_sourceBgn   = (end == 0) ? (0) : (fd->pos[bgn]);            //  Space based (probably).
  layout->_sourceEnd   = (end == 0) ? (0) : (fd->pos[end - 1] + 1);    //  Space based.

  resizeArrayPair(layout->_gappedBases, layout->_gappedQuals, layout->_gappedLen, layout->_gappedMax, end - bgn + 1, resizeArray_doNothing);

  for (uint32 ii=bgn; ii<end; ii++) {
    layout->_gappedBases[ii-bgn] = fd->seq[ii];
    layout->_gappedQuals[ii-bgn] = fd->eqv[ii];
  }

  layout->_gappedLen = end - bgn;

  layout->_gappedBases[layout->_gappedLen] = 0;
  layout->_gappedQuals[layout->_gappedLen] = 0;

  //  One could dump bases and quals here, if so desired.

  ;

  //  Clean up.  Remvoe all the reads[] and datas[] we've loaded.

  for (map<uint32, sqRead     *>::iterator it=reads.begin(); it != reads.end(); ++it)
    delete it->second;

  for (map<uint32, sqReadData *>::iterator it=datas.begin(); it != datas.end(); ++it)
    delete it->second;

  reads.clear();
  datas.clear();

  delete    fd;
  delete [] evidence;
}




int
main(int argc, char **argv) {
  char             *seqName   = 0L;
  char             *corName   = 0L;
  uint32            corVers   = 1;

  char             *exportName = NULL;
  char             *importName = NULL;

  uint32            errorRate = AS_OVS_encodeEvalue(0.015);

  char             *outputPrefix = NULL;
  bool              outputCNS    = false;
  bool              outputFASTQ  = false;
  bool              outputLog    = false;

  uint32            idMin = 1;
  uint32            idMax = UINT32_MAX;
  char             *readListName = NULL;
  set<uint32>       readList;

  uint32            numThreads         = omp_get_max_threads();

  uint32            minOutputCoverage  = 4;
  uint32            minOutputLength    = 1000;
  double            minOlapIdentity    = 0.5;
  double            minOlapLength      = 500;

  bool              trimToAlign        = true;
  bool              restrictToOverlap  = true;

  argc = AS_configure(argc, argv);

  vector<char *>  err;
  int             arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUTS
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {
      corName = argv[++arg];


    } else if (strcmp(argv[arg], "-p") == 0) {   //  OUTPUTS
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-cns") == 0) {
      outputCNS = true;

    } else if (strcmp(argv[arg], "-fastq") == 0) {
      outputFASTQ = true;

    } else if (strcmp(argv[arg], "-log") == 0) {
      outputLog = true;


    } else if (strcmp(argv[arg], "-t") == 0) {   //  COMPUTE RESOURCES
      numThreads = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-f") == 0) {   //  ALGORITHM OPTIONS
      restrictToOverlap = false;


    } else if (strcmp(argv[arg], "-R") == 0) {   //  READ SELECTION
      readListName = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      decodeRange(argv[++arg], idMin, idMax);


    } else if (strcmp(argv[arg], "-cc") == 0) {   //  CONSENSUS
      minOutputCoverage = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-oi") == 0) {
      minOlapIdentity = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-ol") == 0) {
      minOlapLength = atof(argv[++arg]);


    } else if (strcmp(argv[arg], "-export") == 0) {   //  DEBUGGING
      exportName = argv[++arg];

    } else if (strcmp(argv[arg], "-import") == 0) {
      importName = argv[++arg];


    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if ((seqName == NULL) && (importName == NULL))
    err.push_back("ERROR: no seqStore input (-S) supplied.\n");

  if ((corName == NULL) && (importName == NULL))
    err.push_back("ERROR: no corStore input (-C) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "  -S seqStore        mandatory path to seqStore\n");
    fprintf(stderr, "  -C corStore        mandatory path to corStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUTS:\n");
    fprintf(stderr, "  -p prefix          output filename prefix\n");
    fprintf(stderr, "  -cns               enable primary output (to 'prefix.cns')\n");
    fprintf(stderr, "  -fastq             enable fastq output (to 'prefix.fastq')\n");
    fprintf(stderr, "  -log               enable (debug) logging output (to 'prefix.log')\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "RESOURCE PARAMETERS\n");
    fprintf(stderr, "  -t numThreads      number of compute threads to use (default: all)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ALGORITHM PARAMETERS\n");
    fprintf(stderr, "  -f                 align evidence to the full read, ignore overlap position\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "READ SELECTION\n");
    fprintf(stderr, "  -R readsToCorrect  only process reads listed in file 'readsToCorrect'\n");
    fprintf(stderr, "  -r bgn[-end]       only process reads from ID 'bgn' to 'end' (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS\n");
    fprintf(stderr, "  -cc coverage       output:   minimum consensus coverage needed call a corrected base\n");
    fprintf(stderr, "  -cl length         output:   minimum length of corrected region to output as a corrected read\n");
    fprintf(stderr, "  -oi identity       evidence: minimum identity of an aligned evidence read overlap\n");
    fprintf(stderr, "  -ol length         evidence: minimum length   of an aligned evidence read overlap\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DEBUGGING SUPPORT\n");
    fprintf(stderr, "  -export name       write the data used for the computation to file 'name'\n");
    fprintf(stderr, "  -import name       compute using the data in file 'name'\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }


  omp_set_num_threads(numThreads);

  //  Open inputs.

  sqRead_setDefaultVersion(sqRead_raw);

  sqStore *seqStore = NULL;
  tgStore *corStore = NULL;

  if (seqName) {
    fprintf(stderr, "-- Opening seqStore '%s'.\n", seqName);
    seqStore = sqStore::sqStore_open(seqName);
  }

  if (corName) {
    fprintf(stderr, "-- Opening corStore '%s' version %u.\n", corName, corVers);
    corStore = new tgStore(corName, corVers);
  }

  if ((seqStore) &&
      (seqStore->sqStore_getNumReads() < idMax))        //  Limit the range of processing to the
    idMax = seqStore->sqStore_getNumReads();            //  number of reads in the store.

  loadReadList(readListName, idMin, idMax, readList);   //  Further limit to a set of good reads.

  FILE *exportFile = NULL;
  FILE *importFile = NULL;

  if (exportName) {
    fprintf(stderr, "-- Opening export file '%s'.\n", exportName);
    exportFile = AS_UTL_openOutputFile(exportName);
  }

  if (importName) {
    fprintf(stderr, "-- Opening import file '%s'.\n", importName);
    importFile  = AS_UTL_openInputFile(importName);
  }

  //  Open logging and summary files

  FILE *logFile = NULL;
  FILE *cnsFile = NULL;
  FILE *seqFile = NULL;

  cnsFile = AS_UTL_openOutputFile(outputPrefix, '.', "cns",   outputCNS);
  seqFile = AS_UTL_openOutputFile(outputPrefix, '.', "fastq", outputFASTQ);
  logFile = AS_UTL_openOutputFile(outputPrefix, '.', "log",   outputLog);

  //  Initialize processing.

  falconConsensus           *fc = new falconConsensus(minOutputCoverage, minOutputLength, minOlapIdentity, minOlapLength, restrictToOverlap);
  map<uint32, sqRead *>      reads;
  map<uint32, sqReadData *>  datas;

  //  And process.

  fprintf(stdout, "    read    read evidence     corrected\n");
  fprintf(stdout, "      ID  length    reads       regions\n");
  fprintf(stdout, "-------- ------- -------- ------------- ...\n");

  //
  //  If input from a package file, load and process data until there isn't any more.
  //

  if (importFile) {
    tgTig                     *layout = new tgTig();

    FILE  *importedLayouts = AS_UTL_openOutputFile(importName, '.', "layout", (importName != NULL));
    FILE  *importedReads   = AS_UTL_openOutputFile(importName, '.', "fasta",  (importName != NULL));

    while (layout->importData(importFile, reads, datas, NULL, NULL) == true) {
      generateFalconConsensus(fc,
                              layout,
                              seqStore,
                              reads,
                              datas,
                              trimToAlign,
                              minOlapLength);

      if (cnsFile)
        layout->saveToStream(cnsFile);

      if (seqFile)
        layout->dumpFASTQ(seqFile, false);

      delete layout;
      layout = new tgTig();    //  Next loop needs an existing empty layout.
    }

    AS_UTL_closeFile(importedReads);
    AS_UTL_closeFile(importedLayouts);

    //  We properly should clean up the data loaded.

    delete layout;
  }

  //
  //  Otherwise, if we're just dumping data, just dump the data without processing.
  //

  else if (exportFile) {
    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout) {
        fprintf(stdout, "%8u %7u %8u", layout->tigID(), layout->length(), layout->numberOfChildren());

        layout->exportData(exportFile, seqStore, true);
        corStore->unloadTig(layout->tigID());

        fprintf(stdout, "        DUMPED\n");
      }
    }
  }

  //
  //  Otherwise, load and process from a store, the usual processing loop.
  //

  else {
    for (uint32 ii=idMin; ii<=idMax; ii++) {
      if ((readList.size() > 0) &&      //  Skip reads not on the read list,
          (readList.count(ii) == 0))    //  if there actually is a read list.
        continue;

      tgTig *layout = corStore->loadTig(ii);

      if (layout) {
        generateFalconConsensus(fc,
                                layout,
                                seqStore,
                                reads,
                                datas,
                                trimToAlign,
                                minOlapLength);

        if (cnsFile)
          layout->saveToStream(cnsFile);

        if (seqFile)
          layout->dumpFASTQ(seqFile, false);

        corStore->unloadTig(layout->tigID());
      }
    }
  }

  //  Close files and clean up.

  AS_UTL_closeFile(logFile);
  AS_UTL_closeFile(cnsFile);
  AS_UTL_closeFile(seqFile);

  AS_UTL_closeFile(exportFile);
  AS_UTL_closeFile(importFile);

  delete    fc;
  delete    corStore;

  seqStore->sqStore_close();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
