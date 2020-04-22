
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
 *    Brian P. Walenz on 2015-JUN-16
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-MAR-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Nurk beginning on 2019-OCT-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "correctOverlaps.H"

#include "Binomial_Bound.H"




void
Read_Olaps(coParameters *G, sqStore *seqStore);

void
Correct_Frags(coParameters *G, sqStore *seqStore, FILE *correctedReads = NULL);

void
Redo_Olaps(coParameters *G, sqStore *seqStore);


int
main(int argc, char **argv) {
  coParameters  *G = new coParameters();

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      G->seqStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      G->bgnID = atoi(argv[++arg]);
      G->endID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-O") == 0) {  //  -F?  -S Olap_Path
      G->ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      G->errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      G->minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      G->checkTrivialDNA = true;
      fprintf(stderr, "Correcting HIFI reads\n");

    } else if (strcmp(argv[arg], "-c") == 0) {  //  For 'corrections' file input
      G->correctionsName = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {  //  Corrected reads file
      G->correctedName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {  //  For 'erates' output
      G->eratesName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {  //  But we're not threaded!
      G->numThreads = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if (G->seqStorePath == NULL)
    fprintf(stderr, "ERROR: no input sequence store (-S) supplied.\n"), err++;
  if (G->ovlStorePath == NULL)
    fprintf(stderr, "ERROR: no input overlap store (-O) supplied.\n"), err++;
  if (G->correctionsName == NULL)
    fprintf(stderr, "ERROR: no input read corrections file (-c) supplied.\n"), err++;
  if (G->eratesName == NULL)
    fprintf(stderr, "ERROR: no output erates file (-o) supplied.\n"), err++;


  if (err) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore -R bgn end ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S   seqStore           path to a sequence store\n");
    fprintf(stderr, "  -O   ovlStore           path to an overlap store\n");
    fprintf(stderr, "  -R   bgn end            only compute for reads bgn-end\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c   input-name         read corrections from 'input-name'\n");
    fprintf(stderr, "  -o   output-name        write updated error rates to 'output-name'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t   num-threads        not used; only one thread used\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -l   min-len            ignore overlaps shorter than this\n");
    fprintf(stderr, "  -e   max-erate s        ignore overlaps higher than this error\n");
    fprintf(stderr, "  -s   check trival dna   ignore alignment errors in simple sequence\n");
    exit(1);
  }

  //fprintf (stderr, "Quality Threshold = %.2f%%\n", 100.0 * Quality_Threshold);

  //
  //  Initialize Globals
  //

  fprintf(stderr, "Initializing.\n");

  double MAX_ERRORS = 1 + (uint32)(G->errorRate * AS_MAX_READLEN);

  Initialize_Match_Limit(G->Edit_Match_Limit, G->errorRate, MAX_ERRORS);

  for (int32 i=0;  i <= AS_MAX_READLEN;  i++)
    G->Error_Bound[i] = (int)ceil(i * G->errorRate);

  //
  //
  //

  fprintf(stderr, "Opening seqStore '%s'.\n", G->seqStorePath);

  sqStore *seqStore = new sqStore(G->seqStorePath);

  if (G->bgnID < 1)
    G->bgnID = 1;

  if (seqStore->sqStore_lastReadID() < G->endID)
    G->endID = seqStore->sqStore_lastReadID();

  //
  //  Load the reads for the overlaps we are going to be correcting, and apply corrections to them

  fprintf(stderr, "Correcting reads " F_U32 " to " F_U32 ".\n", G->bgnID, G->endID);

  compressedFileWriter *correctedReadsWriter = (G->correctedName == NULL) ? NULL : new compressedFileWriter(G->correctedName);

  //FILE *correctedReads = (G->correctedName == NULL) ? NULL : fopen(G->correctedName, "w");
  Correct_Frags(G, seqStore, (correctedReadsWriter == NULL) ? NULL : correctedReadsWriter->file());
  if (correctedReadsWriter != NULL) {
    //fclose(correctedReads);
    delete correctedReadsWriter;
    fprintf(stderr, "Exiting");
    exit(239);
  }

  //  Load overlaps we're going to correct

  fprintf(stderr, "Loading overlaps.\n");

  Read_Olaps(G, seqStore);

  //  Now sort them on the B iid.

  fprintf(stderr, "Sorting overlaps.\n");

  sort(G->olaps, G->olaps + G->olapsLen, Olap_Info_t_by_bID());

  //  Recompute overlaps

  fprintf(stderr, "Recomputing overlaps.\n");

  Redo_Olaps(G, seqStore);

  delete seqStore;
  seqStore = NULL;

  //  Sort the overlaps back into the original order

  fprintf(stderr, "Sorting overlaps.\n");

  sort(G->olaps, G->olaps + G->olapsLen, Olap_Info_t_by_Order());

  //  Dump the new erates

  fprintf (stderr, "Saving corrected error rates to file %s\n", G->eratesName);

  FILE *fp = AS_UTL_openOutputFile(G->eratesName);

  writeToFile(G->bgnID,    "loid", fp);
  writeToFile(G->endID,    "hiid", fp);
  writeToFile(G->olapsLen, "num",  fp);

  fprintf(stderr, "--Allocate " F_U64 " MB for output error rates.\n",
          (sizeof(uint16) * G->olapsLen) >> 20);

  uint16 *evalue = new uint16 [G->olapsLen];

  for (int32 i=0; i<G->olapsLen; i++)
    evalue[i] = G->olaps[i].evalue;

  writeToFile(evalue, "evalue", G->olapsLen, fp);

  delete [] evalue;

  AS_UTL_closeFile(fp, G->eratesName);

  //  Finished.

  //fprintf (stderr, "%d/%d failed/total alignments (%.1f%%)\n",
  //         Failed_Alignments_Ct, Total_Alignments_Ct,
  //         Total_Alignments_Ct == 0 ? 0.0 : (100.0 * Failed_Alignments_Ct) / Total_Alignments_Ct);

  delete G;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}


