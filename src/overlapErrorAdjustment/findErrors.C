
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

#include "findErrors.H"

#include "Binomial_Bound.H"

void
Process_Olap(Olap_Info_t        *olap,
             char               *b_seq,
             bool                shredded,
             Thread_Work_Area_t *wa);

void
Read_Frags(feParameters   *G,
           sqStore        *seqStore);

void
Read_Olaps(feParameters   *G,
           sqStore        *seqStore);

void
Output_Corrections(feParameters *G);




//  Read fragments lo_frag..hi_frag (INCLUSIVE) from store and save the ids and sequences of those
//  with overlaps to fragments in global Frag .

static
void
extractReads(feParameters *G,
             sqStore      *seqStore,
             Frag_List_t  *fl,
             uint64       &nextOlap) {

  //  Clear the buffer.

  fl->readsLen = 0;
  fl->basesLen = 0;

  //  The original converted to lowercase, and made non-acgt be 'a'.

  char  filter[256];

  for (uint32 i=0; i<256; i++)
    filter[i] = 'a';

  filter['A'] = filter['a'] = 'a';
  filter['C'] = filter['c'] = 'c';
  filter['G'] = filter['g'] = 'g';
  filter['T'] = filter['t'] = 't';

  //  Return if we've exhausted the overlaps.

  if (nextOlap >= G->olapsLen)
    return;

  //  Count the amount of stuff we're loading.

  uint64 lastOlap = nextOlap;
  uint32 loID     = G->olaps[lastOlap].b_iid;  //  Actual ID we're extracting
  uint32 hiID     = loID;
  uint64 maxBases = 512 * 1024 * 1024;

  //  Find the highest read ID that we can load without exceeding maxBases.

  while ((fl->basesLen < maxBases) &&
         (lastOlap     < G->olapsLen)) {
    hiID = G->olaps[lastOlap].b_iid;                        //  Grab the ID of the overlap we're at.

    fl->readsLen += 1;                                      //  Add the read to our set.
    fl->basesLen += seqStore->sqStore_getReadLength(hiID) + 1;

    lastOlap++;                                             //  Advance to the next overlap
    while ((lastOlap < G->olapsLen) &&                      //
           (G->olaps[lastOlap].b_iid == hiID))              //  If we've exceeded the max size or hit the last overlap,
      lastOlap++;                                           //  the loop will stop on the next iteration.
  }

  //  If nothing to load, just return.

  if (fl->readsLen == 0)
    return;

  //  Report what we're going to do.

  fprintf(stderr, "extractReads()-- Loading reads " F_U32 " to " F_U32 " (" F_U32 " reads with " F_U64 " bases) overlaps " F_U64 " through " F_U64 ".\n",
          loID, hiID, fl->readsLen, fl->basesLen, nextOlap, lastOlap);

  //  Ensure there is space.

  if (fl->readsMax < fl->readsLen) {
    delete [] fl->readIDs;
    delete [] fl->readBases;

    fl->readsMax  = 12 * fl->readsLen / 10;
    fl->readIDs   = new uint32 [fl->readsMax];
    fl->readBases = new char * [fl->readsMax];
  }

  if (fl->basesMax < fl->basesLen) {
    delete [] fl->bases;

    fl->basesMax    = 12 * fl->basesLen / 10;
    fl->bases       = new char [fl->basesMax];
  }

  //  Load the sequence data for reads loID to hiID, as long as the read has an overlap.

  sqRead *read = new sqRead;

  fl->readsLen = 0;
  fl->basesLen = 0;

  while ((loID <= hiID) &&
         (nextOlap < G->olapsLen)) {
    seqStore->sqStore_getRead(loID, read);

    fl->readIDs[fl->readsLen]   = loID;                          //  Save the ID of _this_ read.
    fl->readBases[fl->readsLen] = fl->bases + fl->basesLen;      //  Set the data pointer to where this read should start.

    uint32  readLen    = read->sqRead_length();
    char   *readBases  = read->sqRead_sequence();

    for (uint32 bb=0; bb<readLen; bb++)
      fl->readBases[fl->readsLen][bb] = filter[readBases[bb]];

    fl->readBases[fl->readsLen][readLen] = 0;                    //  All good reads end.

    fl->basesLen += read->sqRead_length() + 1;                   //  Update basesLen to account for this read.
    fl->readsLen += 1;                                           //  And note that we loaded a read.

    nextOlap++;                                                  //  Advance past all the overlaps for this read.
    while ((nextOlap < G->olapsLen) &&
           (G->olaps[nextOlap].b_iid == loID))
      nextOlap++;

    if (nextOlap < G->olapsLen)                                  //  If we have valid overlap, grab the read ID.
      loID = G->olaps[nextOlap].b_iid;                           //  If we don't have a valid overlap, the loop will stop.
  }

  delete read;


  fprintf(stderr, "extractReads()-- Loaded.\n");
}



//  Process all old fragments in  Internal_seqStore. Only
//  do overlaps/corrections with fragments where
//    frag_iid % Num_PThreads == thread_id

void *
processThread(void *ptr) {
  Thread_Work_Area_t  *wa = (Thread_Work_Area_t *)ptr;

  for (int32 i=0; i<wa->frag_list->readsLen; i++) {
    int32  skip_id = -1;

    while (wa->frag_list->readIDs[i] > wa->G->olaps[wa->nextOlap].b_iid) {
      if (wa->G->olaps[wa->nextOlap].b_iid != skip_id) {
        fprintf(stderr, "SKIP:  b_iid = %d\n", wa->G->olaps[wa->nextOlap].b_iid);
        skip_id = wa->G->olaps[wa->nextOlap].b_iid;
      }
      wa->nextOlap++;
    }

    if (wa->frag_list->readIDs[i] != wa->G->olaps[wa->nextOlap].b_iid) {
      fprintf (stderr, "ERROR:  Lists don't match\n");
      fprintf (stderr, "frag_list iid = %d  nextOlap = %d  i = %d\n",
               wa->frag_list->readIDs[i],
               wa->G->olaps[wa->nextOlap].b_iid, i);
      exit (1);
    }

    wa->rev_id = UINT32_MAX;

    while ((wa->nextOlap < wa->G->olapsLen) && (wa->G->olaps[wa->nextOlap].b_iid == wa->frag_list->readIDs[i])) {
      if (wa->G->olaps[wa->nextOlap].a_iid % wa->G->numThreads == wa->thread_id) {
        Process_Olap(wa->G->olaps + wa->nextOlap,
                     wa->frag_list->readBases[i],
                     false,  //  shredded
                     wa);
      }

      wa->nextOlap++;
    }
  }

  pthread_exit(ptr);

  return(NULL);
}



//  Read old fragments in  seqStore  that have overlaps with
//  fragments in  Frag. Read a batch at a time and process them
//  with multiple pthreads.  Each thread processes all the old fragments
//  but only changes entries in  Frag  that correspond to its thread
//  ID.  Recomputes the overlaps and records the vote information about
//  changes to make (or not) to fragments in  Frag .


static
void
processReads(feParameters *G,
             sqStore      *seqStore,
             uint64       &passedOlaps,
             uint64       &failedOlaps) {

  pthread_attr_t  attr;

  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr, THREAD_STACKSIZE);

  pthread_t           *thread_id = new pthread_t          [G->numThreads];
  Thread_Work_Area_t  *thread_wa = new Thread_Work_Area_t [G->numThreads];

  for (uint32 i=0; i<G->numThreads; i++) {
    thread_wa[i].thread_id    = i;
    thread_wa[i].nextOlap     = 0;
    thread_wa[i].G            = G;
    thread_wa[i].frag_list    = NULL;
    thread_wa[i].rev_id       = UINT32_MAX;
    thread_wa[i].passedOlaps  = 0;
    thread_wa[i].failedOlaps  = 0;

    memset(thread_wa[i].rev_seq, 0, sizeof(char) * AS_MAX_READLEN);

    double MAX_ERRORS = 1 + (uint32)(G->errorRate * AS_MAX_READLEN);

    thread_wa[i].ped.initialize(G, G->errorRate);
  }

  uint64 frstOlap = 0;
  uint64 nextOlap = 0;

  Frag_List_t   frag_list_1;
  Frag_List_t   frag_list_2;

  Frag_List_t  *curr_frag_list = &frag_list_1;
  Frag_List_t  *next_frag_list = &frag_list_2;

  extractReads(G, seqStore, curr_frag_list, nextOlap);

  while (curr_frag_list->readsLen > 0) {

    // Process fragments in curr_frag_list in background

    fprintf(stderr, "processReads()-- Launching compute.\n");

    for (uint32 i=0; i<G->numThreads; i++) {
      thread_wa[i].nextOlap  = frstOlap;
      thread_wa[i].frag_list = curr_frag_list;

      int status = pthread_create(thread_id + i, &attr, processThread, thread_wa + i);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    // Read next batch of fragments

    frstOlap = nextOlap;

    extractReads(G, seqStore, next_frag_list, nextOlap);

    // Wait for background processing to finish

    fprintf(stderr, "processReads()-- Waiting for compute.\n");

    for (uint32 i=0; i<G->numThreads; i++) {
      void  *ptr;

      int status = pthread_join(thread_id[i], &ptr);

      if (status != 0)
        fprintf(stderr, "pthread_join error: %s\n", strerror(status)), exit(1);
    }

    //  Swap the lists and compute another block

    {
      Frag_List_t *s = curr_frag_list;
      curr_frag_list = next_frag_list;
      next_frag_list = s;
    }
  }

  //  Threads all done, sum up stats.

  passedOlaps = 0;
  failedOlaps = 0;

  for (uint32 i=0; i<G->numThreads; i++) {
    passedOlaps += thread_wa[i].passedOlaps;
    failedOlaps += thread_wa[i].failedOlaps;
  }

  delete [] thread_id;
  delete [] thread_wa;
}









int
main(int argc, char **argv) {
  feParameters  *G = new feParameters();

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-S") == 0) {
      G->seqStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-R") == 0) {
      G->bgnID = atoi(argv[++arg]);
      G->endID = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-O") == 0) {
      G->ovlStorePath = argv[++arg];

    } else if (strcmp(argv[arg], "-e") == 0) {
      G->errorRate = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-l") == 0) {
      G->minOverlap = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {  //  For 'corrections' file output
      G->outputFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      G->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-d") == 0) {
      G->Degree_Threshold = strtol(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-k") == 0) {
      G->Kmer_Len = strtol(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-p") == 0) {
      G->Use_Haplo_Ct = false;

    } else if (strcmp(argv[arg], "-V") == 0) {
      G->Vote_Qualify_Len = strtol(argv[++arg], NULL, 10);

    } else if (strcmp(argv[arg], "-x") == 0) {
      G->End_Exclude_Len = strtol(argv[++arg], NULL, 10);

    } else {
      fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (G->seqStorePath == NULL)
    err++;
  if (G->ovlStorePath == NULL)
    err++;
  if (G->numThreads == 0)
    err++;
  if (G->bgnID > G->endID)
    err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s -S seqStore -O ovlStore -R bgn-end ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -S   seqStore           path to a sequence store\n");
    fprintf(stderr, "  -O   ovlStore           path to an overlap store\n");
    fprintf(stderr, "  -R   bgn end            only compute for reads bgn-end\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o   output-name        write corrections to 'output-name'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e   error-rate         expected error rate in alignments\n");
    fprintf(stderr, "  -l   min-overlap        \n");
    fprintf(stderr, "  -t   num-threads        \n");
    fprintf(stderr, "  -d   degree-threshold   set keep flag if fewer than this many overlaps\n");
    fprintf(stderr, "  -k   kmer-size          minimum exact-match region to prevent change\n");
    fprintf(stderr, "  -p                      don't use the haplo_ct\n");
    fprintf(stderr, "  -V   vote-len           number of exact match bases around an error to vote to change\n");
    fprintf(stderr, "  -x   end-exclude-len    length of end of exact match to exclude in preventing change\n");

    if (G->seqStorePath == NULL)
      fprintf(stderr, "ERROR: no sequence store (-S) supplied.\n");
    if (G->ovlStorePath == NULL)
      fprintf(stderr, "ERROR: no overlap store (-O) supplied.\n");
    if (G->numThreads == 0)
      fprintf(stderr, "ERROR: number of compute threads (-t) must be larger than zero.\n");
    if (G->bgnID > G->endID)
      fprintf(stderr, "ERROR: read range (-R) %u-%u invalid.\n", G->bgnID, G->endID);

    exit(1);
  }

  //  Initialize Globals

  double MAX_ERRORS = 1 + (uint32)(G->errorRate * AS_MAX_READLEN);

  Initialize_Match_Limit(G->Edit_Match_Limit, G->errorRate, MAX_ERRORS);

  for  (uint32 i = 0;  i <= AS_MAX_READLEN;  i++)
    G->Error_Bound[i] = (int)ceil(i * G->errorRate);

  //  Load data.

  sqStore *seqStore = new sqStore(G->seqStorePath);

  if (G->bgnID < 1)
    G->bgnID = 1;

  if (seqStore->sqStore_lastReadID() < G->endID)
    G->endID = seqStore->sqStore_lastReadID();

  Read_Frags(G, seqStore);
  Read_Olaps(G, seqStore);

  //  Sort overlaps, process each.

  sort(G->olaps, G->olaps + G->olapsLen);

  uint64  passedOlaps = 0;
  uint64  failedOlaps = 0;

  processReads(G, seqStore, passedOlaps, failedOlaps);

  //  All done.  Sum up what we did.

  fprintf(stderr, "\n");
  fprintf(stderr, "Passed overlaps = %10" F_U64P " %8.4f%%\n", passedOlaps, 100.0 * passedOlaps / (failedOlaps + passedOlaps));
  fprintf(stderr, "Failed overlaps = %10" F_U64P " %8.4f%%\n", failedOlaps, 100.0 * failedOlaps / (failedOlaps + passedOlaps));

  //  Dump output.

  //Output_Details(G);
  Output_Corrections(G);

  //  Cleanup and exit!

  delete seqStore;

  delete G;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}

