
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
 *    Brian P. Walenz beginning on 2018-FEB-08
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

#include "strings.H"
#include "files.H"

#include "AS_UTL_decodeRange.H"
#include "AS_UTL_reverseComplement.H"

#include <map>

using namespace std;




int
main(int argc, char **argv) {
  char             *seqName   = 0L;
  char             *corName   = 0L;
  uint32            corVers   = 1;

  char             *prefix = NULL;

  uint32            idMin = 0;
  uint32            idMax = UINT32_MAX;
  char             *haplotypeListPrefix = NULL;
  map<char*, FILE*> haplotypeList;

  uint32            minRatio           = 1;
  uint32            minOutputLength    = 500;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {   //  INPUTS
      seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];


    } else if (strcmp(argv[arg], "-cr") == 0) {
      minRatio = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-h") == 0) {
       ++arg; // skip the -h
       while (arg < argc && argv[arg][0] != '-') {
          haplotypeList[argv[arg++]] = NULL;
       }
       --arg;

    } else if (strcmp(argv[arg], "-cl") == 0) {
      minOutputLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-b") == 0) {   //  READ SELECTION
      idMin = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-e") == 0) {
      idMax = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (seqName == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -S seqStore ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS (all mandatory)\n");
    fprintf(stderr, "  -S seqStore      mandatory path to seqStore\n");
    fprintf(stderr, "  -p prefix        output prefix name, for logging and summary report\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONSENSUS PARAMETERS\n");
    fprintf(stderr, "  -cr ratio        minimum ratio between best and second best to classify\n");
    fprintf(stderr, "  -cl length       minimum length of output read\n");
    fprintf(stderr, "\n");

    if (seqName == NULL)
      fprintf(stderr, "ERROR: no sequence store input (-S) supplied.\n");
    exit(1);
  }


  //  Open inputs.

  sqStore  *seqStore = sqStore::sqStore_open(seqName);
  uint32    numReads = seqStore->sqStore_getNumReads();

  //  Decide what reads to operate on.

  if (numReads < idMax)
    idMax = numReads;



  // open all the haplotype read input and output files, assume we have few enough haplotypes that we won't hit max file limits
  map<char*, FILE*> outputFasta;
  char outputName[256];

  for (map<char*,FILE*>::iterator it=haplotypeList.begin(); it!=haplotypeList.end(); ++it) {
     it->second = AS_UTL_openInputFile(prefix, '.', it->first);
     sprintf(outputName, "%s.%s", prefix, it->first);
     outputFasta[it->first] = AS_UTL_openOutputFile(outputName, '.', "fasta");
  }
  outputFasta["unknown"] = AS_UTL_openOutputFile(prefix, '.', "unknown.fasta");

  // now loop reads and write
  char       *ovStr = new char [1024];

  fprintf(stderr, "Launched with range %d - %d\n", idMin, idMax);
  for (uint32 ii=idMin; ii<=idMax; ii++) {
     char *haplotype = NULL;
     double bestCount = 0;
     double secondBest = 0;
     double total = 0;
     for (map<char*,FILE*>::iterator it=haplotypeList.begin(); it!=haplotypeList.end(); ++it) {
        // read a line
        // make sure id matches, or die
        fgets(ovStr, 1024, it->second);
        if (ovStr == NULL) {
           fprintf(stderr, "Error: failed to read input line for haplotype %s\n", it->first);
           exit(1);
        }
        splitToWords  W(ovStr);
        // skip the read in the name (that is seqStore output readX clr= so W[0] needs to ignore the "read" text
        char *rid = W[0] + 4;
        fprintf(stderr, "For haplotype %s read %s with count %s\n", it->first, rid, W[4]);

        uint32 aid = strtouint32(rid);
        uint64 totalMers = strtouint64(W[3]);	// haplotype specific mers
        uint64 count = strtouint64(W[4]);		// haplotype mers in the read
        double scaledCount = (double) count / totalMers;
        fprintf(stderr, "After scaling for %s by " F_U64 " the count is %f\n", it->first, totalMers, scaledCount);

        if (aid != ii) {
           fprintf(stderr, "Error: expected %d but got read %d\n", ii, aid);
           exit(1);
        }
        total += scaledCount;
        if (scaledCount > 0) {
        	if (scaledCount <= bestCount && scaledCount > secondBest)
        		secondBest = scaledCount;
        	else if (scaledCount > bestCount) {
			   secondBest = bestCount;
			   bestCount = scaledCount;
			   haplotype = it->first;
			}
        }
     }
     sqReadData read;
     seqStore->sqStore_loadReadData(ii, &read);

     if (read.sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw) < minOutputLength) {
        continue;
     }

     // now classify read and write it
     if ((secondBest == 0 && bestCount != 0) || ((double)bestCount / secondBest > minRatio)) {
        // write to haplotype
     } else {
        // ambiguous
        haplotype="unknown";
     }
     fprintf(stderr, "Processing read %d classified as %s with counts %f and %f\n", ii, haplotype, bestCount, secondBest);

     AS_UTL_writeFastA(outputFasta[haplotype], read.sqReadData_getRawSequence(), read.sqReadData_getRead()->sqRead_sequenceLength(sqRead_raw), 0,
                       ">read" F_U32 "\n",
                       ii);
  }

  for (map<char*,FILE*>::iterator it=haplotypeList.begin(); it!=haplotypeList.end(); ++it) {
     fclose(it->second);
     fclose(outputFasta[it->first]);
  }
  fclose(outputFasta["unknown"]);

  seqStore->sqStore_close();

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
