

const char *mainid = "$Id:  $";

#include "correctOverlaps.H"

#include "Binomial_Bound.H"




void
Read_Olaps(coParameters &G);

void
Correct_Frags(coParameters &G, gkStore *gkpStore);

void
Redo_Olaps(coParameters &G, gkStore *gkpStore);


int
main(int argc, char **argv) {
  coParameters  G;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      G.gkpStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {  //  -F?  -S Olap_Path
      G.ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-C") == 0) {  //  For 'corrections'
      G.correctionsName = argv[++arg];

    } else if (strcmp(argv[arg], "-E") == 0) {  //  For 'erates'
      G.eratesName = argv[++arg];


    } else if (strcmp(argv[arg], "-b") == 0) {
      G.bgnID = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-e") == 0) {
      G.endID = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-t") == 0) {
      G.numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-error") == 0) {
      AS_OVL_ERROR_RATE = atof(argv[++arg]);



    } else {
      err++;
    }

    arg++;
  }

  if (G.gkpStorePath == NULL)
    fprintf(stderr, "ERROR: no input gatekeeper store (-G) supplied.\n"), err++;
  if (G.ovlStorePath == NULL)
    fprintf(stderr, "ERROR: no input overlap store (-O) supplied.\n"), err++;
  if (G.correctionsName == NULL)
    fprintf(stderr, "ERROR: no input read corrections file (-C) supplied.\n"), err++;
  if (G.eratesName == NULL)
    fprintf(stderr, "ERROR: no output erates file (-E) supplied.\n"), err++;


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

  {
    for (int32 i=0;  i <= ERRORS_FOR_FREE;  i++)
      G.Edit_Match_Limit[i] = 0;

    int32 start = 1;
    for (int32 e = ERRORS_FOR_FREE + 1;  e < MAX_ERRORS;  e++) {
      start = Binomial_Bound(e - ERRORS_FOR_FREE, AS_OVL_ERROR_RATE, start);

      G.Edit_Match_Limit[e] = start - 1;

      assert(G.Edit_Match_Limit[e] >= G.Edit_Match_Limit[e - 1]);
    }

    for (int32 i=0;  i <= AS_MAX_READLEN;  i++)
      G.Error_Bound[i] = (int) (i * AS_OVL_ERROR_RATE);
  }

  //
  //
  //

  fprintf(stderr, "Opening gkpStore '%s'.\n", G.gkpStorePath);

  gkStore *gkpStore = new gkStore(G.gkpStorePath);

  if (G.bgnID < 1)
    G.bgnID = 1;

  if (gkpStore->gkStore_getNumReads() < G.endID)
    G.endID = gkpStore->gkStore_getNumReads();

  //  Load the reads for the overlaps we are going to be correcting, and apply corrections to them

  fprintf(stderr, "Correcting reads "F_U32" to "F_U32".\n", G.bgnID, G.endID);

  Correct_Frags(G, gkpStore);

  //  Load overlaps we're going to correct

  Read_Olaps(G);
 
  //  Now sort them on the B iid.

  sort(G.olaps, G.olaps + G.olapsLen, Olap_Info_t_by_bID());

  //  Recompute overlaps

  Redo_Olaps(G, gkpStore);

  //  Sort the overlaps back into the original order

  sort(G.olaps, G.olaps + G.olapsLen, Olap_Info_t_by_Order());

  //  Dump the new erates

  fprintf (stderr, "Saving corrected error rates to file %s\n", G.eratesName);

  {
    errno = 0;
    FILE *fp = fopen(G.eratesName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", G.eratesName, strerror(errno)), exit(1);

    AS_UTL_safeWrite(fp, &G.bgnID,    "loid", sizeof(int32),  1);
    AS_UTL_safeWrite(fp, &G.endID,    "hiid", sizeof(int32),  1);
    AS_UTL_safeWrite(fp, &G.olapsLen, "num",  sizeof(uint64), 1);

    uint16 *evalue = new uint16 [G.olapsLen];

    for (int32 i=0; i<G.olapsLen; i++)
      evalue[i] = G.olaps[i].evalue;

    AS_UTL_safeWrite(fp, evalue, "evalue", sizeof(uint16), G.olapsLen);

    delete [] evalue;

    fclose(fp);
  }

  //  Finished.

  //fprintf (stderr, "%d/%d failed/total alignments (%.1f%%)\n",
  //         Failed_Alignments_Ct, Total_Alignments_Ct,
  //         Total_Alignments_Ct == 0 ? 0.0 : (100.0 * Failed_Alignments_Ct) / Total_Alignments_Ct);

  exit(0);
}


