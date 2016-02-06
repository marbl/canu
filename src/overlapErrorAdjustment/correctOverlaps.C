
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "correctOverlaps.H"

#include "Binomial_Bound.H"




void
Read_Olaps(coParameters *G, gkStore *gkpStore);

void
Correct_Frags(coParameters *G, gkStore *gkpStore);

void
Redo_Olaps(coParameters *G, gkStore *gkpStore);


int
main(int argc, char **argv) {
  coParameters  *G = new coParameters();

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      G->gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      G->bgnID = atoi(argv[++arg]);
      G->endID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-O") == 0) {  //  -F?  -S Olap_Path
      G->ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      G->errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      G->minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-c") == 0) {  //  For 'corrections' file input
      G->correctionsName = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {  //  For 'erates' output
      G->eratesName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {  //  But we're not threaded!
      G->numThreads = atoi(argv[++arg]);

    } else {
      err++;
    }

    arg++;
  }

  if (G->gkpStorePath == NULL)
    fprintf(stderr, "ERROR: no input gatekeeper store (-G) supplied.\n"), err++;
  if (G->ovlStorePath == NULL)
    fprintf(stderr, "ERROR: no input overlap store (-O) supplied.\n"), err++;
  if (G->correctionsName == NULL)
    fprintf(stderr, "ERROR: no input read corrections file (-c) supplied.\n"), err++;
  if (G->eratesName == NULL)
    fprintf(stderr, "ERROR: no output erates file (-o) supplied.\n"), err++;


  if (err) {
    fprintf(stderr, "USAGE:  %s [-d <dna-file>] [-o <ovl_file>] [-q <quality>]\n", argv[0]);
    fprintf(stderr, "            [-x <del_file>] [-F OlapFile] [-S OlapStore]\n");
    fprintf(stderr, "            [-c <cgb_file>] [-e <erate_file>\n");
    fprintf(stderr, "           <gkpStore> <CorrectFile> <lo> <hi>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Recalculates overlaps for frags  <lo> .. <hi>  in\n");
    fprintf(stderr, " <gkpStore>  using corrections in  <CorrectFile> \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-e <erate-file>  specifies binary file to dump corrected erates to\n");
    fprintf(stderr, "                 for later updating of olap store by  update-erates \n");
    fprintf(stderr, "-F             specify file of sorted overlaps to use (in the format\n");
    fprintf(stderr, "               produced by  get-olaps\n");
    fprintf(stderr, "-o <ovl_file>  specifies name of file to which OVL messages go\n");
    fprintf(stderr, "-q <quality>   overlaps less than this error rate are\n");
    fprintf(stderr, "               automatically output\n");
    fprintf(stderr, "-S             specify the binary overlap store containing overlaps to use\n");
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

  fprintf(stderr, "Opening gkpStore '%s'.\n", G->gkpStorePath);

  gkStore *gkpStore = gkStore::gkStore_open(G->gkpStorePath);

  if (G->bgnID < 1)
    G->bgnID = 1;

  if (gkpStore->gkStore_getNumReads() < G->endID)
    G->endID = gkpStore->gkStore_getNumReads();

  //  Load the reads for the overlaps we are going to be correcting, and apply corrections to them

  fprintf(stderr, "Correcting reads "F_U32" to "F_U32".\n", G->bgnID, G->endID);

  Correct_Frags(G, gkpStore);

  //  Load overlaps we're going to correct

  Read_Olaps(G, gkpStore);

  //  Now sort them on the B iid.

  sort(G->olaps, G->olaps + G->olapsLen, Olap_Info_t_by_bID());

  //  Recompute overlaps

  Redo_Olaps(G, gkpStore);

  //  Sort the overlaps back into the original order

  sort(G->olaps, G->olaps + G->olapsLen, Olap_Info_t_by_Order());

  //  Dump the new erates

  fprintf (stderr, "Saving corrected error rates to file %s\n", G->eratesName);

  {
    errno = 0;
    FILE *fp = fopen(G->eratesName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", G->eratesName, strerror(errno)), exit(1);

    AS_UTL_safeWrite(fp, &G->bgnID,    "loid", sizeof(int32),  1);
    AS_UTL_safeWrite(fp, &G->endID,    "hiid", sizeof(int32),  1);
    AS_UTL_safeWrite(fp, &G->olapsLen, "num",  sizeof(uint64), 1);

    uint16 *evalue = new uint16 [G->olapsLen];

    for (int32 i=0; i<G->olapsLen; i++)
      evalue[i] = G->olaps[i].evalue;

    AS_UTL_safeWrite(fp, evalue, "evalue", sizeof(uint16), G->olapsLen);

    delete [] evalue;

    fclose(fp);
  }

  //  Finished.

  //fprintf (stderr, "%d/%d failed/total alignments (%.1f%%)\n",
  //         Failed_Alignments_Ct, Total_Alignments_Ct,
  //         Total_Alignments_Ct == 0 ? 0.0 : (100.0 * Failed_Alignments_Ct) / Total_Alignments_Ct);

  exit(0);
}


