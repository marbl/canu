
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
 *    Brian P. Walenz from 2015-MAY-29 to 2015-JUL-01
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
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
           gkStore        *gkpStore);

void
Read_Olaps(feParameters   *G,
           gkStore        *gkpStore);

void
Output_Corrections(feParameters *G);




//  From overlapInCore.C
int
Binomial_Bound (int e, double p, int Start, double Limit);



//  Read fragments lo_frag..hi_frag (INCLUSIVE) from store and save the ids and sequences of those
//  with overlaps to fragments in global Frag .

static
void
Extract_Needed_Frags(feParameters *G,
                     gkStore      *gkpStore,
                     uint32        loID,
                     uint32        hiID,
                     Frag_List_t  *fl,
                     uint64       &nextOlap) {

  //  The original converted to lowercase, and made non-acgt be 'a'.

  char  filter[256];

  for (uint32 i=0; i<256; i++)
    filter[i] = 'a';

  filter['A'] = filter['a'] = 'a';
  filter['C'] = filter['c'] = 'c';
  filter['G'] = filter['g'] = 'g';
  filter['T'] = filter['t'] = 't';

  //  Count the amount of stuff we're loading.

  fl->readsLen = 0;
  fl->basesLen = 0;

  uint64 lastOlap = nextOlap;
  uint32 ii       = 0;                        //  Index into reads arrays
  uint32 fi       = G->olaps[lastOlap].b_iid;  //  Actual ID we're extracting

  assert(loID <= fi);

  fprintf(stderr, "Extract_Needed_Frags()--  Loading used reads between "F_U32" and "F_U32".\n",
          fi, hiID);

  fprintf(stderr, "Extract_Needed_Frags()--  At overlap "F_U64"\n", lastOlap);

  while (fi <= hiID) {
    gkRead *read = gkpStore->gkStore_getRead(fi);

    fl->readsLen += 1;
    fl->basesLen += read->gkRead_sequenceLength() + 1;

    //  Advance to the next overlap

    lastOlap++;
    while ((lastOlap < G->olapsLen) && (G->olaps[lastOlap].b_iid == fi))
      lastOlap++;
    fi = (lastOlap < G->olapsLen) ? G->olaps[lastOlap].b_iid : hiID + 1;
  }

  fprintf(stderr, "Extract_Needed_Frags()--  Loading reads for overlaps "F_U64" to "F_U64"\n",
          nextOlap, lastOlap);
  fprintf(stderr, "Extract_Needed_Frags()--  reads "F_U32" bases "F_U64"\n", fl->readsLen, fl->basesLen);

  //  Ensure there is space.

  if (fl->readsMax < fl->readsLen) {
    delete [] fl->readIDs;
    delete [] fl->readBases;

    fprintf(stderr, "realloc reads from "F_U32" to "F_U32"\n", fl->readsMax, 12 * fl->readsLen / 10);

    fl->readIDs   = new uint32 [12 * fl->readsLen / 10];
    fl->readBases = new char * [12 * fl->readsLen / 10];

    fl->readsMax  = 12 * fl->readsLen / 10;
  }

  if (fl->basesMax < fl->basesLen) {
    delete [] fl->bases;

    fprintf(stderr, "realloc bases from "F_U64" to "F_U64"\n", fl->basesMax, 12 * fl->basesLen / 10);

    fl->bases       = new char [12 * fl->basesLen / 10];

    fl->basesMax    = 12 * fl->basesLen / 10;
  }

  //  Load.  This is complicated by loading only the reads that have overlaps we care about.

  fl->readsLen = 0;
  fl->basesLen = 0;

  gkReadData *readData = new gkReadData;

  ii = 0;
  fi = G->olaps[nextOlap].b_iid;

  assert(loID <= fi);

  while (fi <= hiID) {
    gkRead *read       = gkpStore->gkStore_getRead(fi);

    fl->readIDs[ii]     = fi;
    fl->readBases[ii]   = fl->bases + fl->basesLen;
    fl->basesLen       += read->gkRead_sequenceLength() + 1;

    gkpStore->gkStore_loadReadData(read, readData);

    uint32  readLen    = read->gkRead_sequenceLength();
    char   *readBases  = readData->gkReadData_getSequence();

    for (uint32 bb=0; bb<readLen; bb++)
      fl->readBases[ii][bb] = filter[readBases[bb]];

    fl->readBases[ii][readLen] = 0;  //  All good reads end.

    ii++;

    //  Advance to the next overlap.

    nextOlap++;
    while ((nextOlap < G->olapsLen) && (G->olaps[nextOlap].b_iid == fi))
      nextOlap++;
    fi = (nextOlap < G->olapsLen) ? G->olaps[nextOlap].b_iid : hiID + 1;
  }

  delete readData;

  fl->readsLen = ii;

  fprintf(stderr, "Extract_Needed_Frags()--  Loaded "F_U32" reads (%.4f%%).  Loaded IDs "F_U32" through "F_U32".\n",
          fl->readsLen, 100.0 * fl->readsLen / (hiID - 1 - loID),
          fl->readIDs[0], fl->readIDs[fl->readsLen-1]);
}



//  Process all old fragments in  Internal_gkpStore. Only
//  do overlaps/corrections with fragments where
//    frag_iid % Num_PThreads == thread_id

void *
Threaded_Process_Stream(void *ptr) {
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

  //pthread_mutex_lock(& Print_Mutex);
  //fprintf(stderr, "Thread %d processed %d olaps\n", wa->thread_id, olap_ct);
  //pthread_mutex_unlock(& Print_Mutex);

  pthread_exit(ptr);

  return(NULL);
}



//  Read old fragments in  gkpStore  that have overlaps with
//  fragments in  Frag. Read a batch at a time and process them
//  with multiple pthreads.  Each thread processes all the old fragments
//  but only changes entries in  Frag  that correspond to its thread
//  ID.  Recomputes the overlaps and records the vote information about
//  changes to make (or not) to fragments in  Frag .


static
void
Threaded_Stream_Old_Frags(feParameters *G, gkStore *gkpStore) {

  pthread_attr_t  attr;

  pthread_mutex_init(&G->Print_Mutex, NULL);

  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr, THREAD_STACKSIZE);

  pthread_t           *thread_id = new pthread_t         [G->numThreads];
  Thread_Work_Area_t  *thread_wa = new Thread_Work_Area_t[G->numThreads];

  for (uint32 i=0; i<G->numThreads; i++) {
    thread_wa[i].thread_id    = i;
    thread_wa[i].loID         = 0;
    thread_wa[i].hiID         = 0;
    thread_wa[i].nextOlap     = 0;
    thread_wa[i].G            = G;
    thread_wa[i].frag_list    = NULL;
    thread_wa[i].rev_id       = UINT32_MAX;
    thread_wa[i].failedOlaps  = 0;

    memset(thread_wa[i].rev_seq, 0, sizeof(char) * AS_MAX_READLEN);

    double MAX_ERRORS = 1 + (uint32)(G->errorRate * AS_MAX_READLEN);

    thread_wa[i].ped.initialize(G, G->errorRate);
  }

  uint32 loID  = G->olaps[0].b_iid;
  uint32 hiID  = loID + FRAGS_PER_BATCH - 1;

  uint32 endID = G->olaps[G->olapsLen - 1].b_iid;

  if (hiID > endID)
    hiID = endID;

  uint64 frstOlap = 0;
  uint64 nextOlap = 0;

  Frag_List_t   frag_list_1;
  Frag_List_t   frag_list_2;

  Frag_List_t  *curr_frag_list = &frag_list_1;
  Frag_List_t  *next_frag_list = &frag_list_2;

  Extract_Needed_Frags(G, gkpStore, loID, hiID, curr_frag_list, nextOlap);

  while (loID <= endID) {

    // Process fragments in curr_frag_list in background

    for (uint32 i=0; i<G->numThreads; i++) {
      thread_wa[i].loID      = loID;
      thread_wa[i].hiID      = hiID;
      thread_wa[i].nextOlap  = frstOlap;
      thread_wa[i].frag_list = curr_frag_list;

      int status = pthread_create(thread_id + i, &attr, Threaded_Process_Stream, thread_wa + i);

      if (status != 0)
        fprintf(stderr, "pthread_create error:  %s\n", strerror(status)), exit(1);
    }

    // Read next batch of fragments

    loID = hiID + 1;

    if (loID <= endID) {
      hiID = loID + FRAGS_PER_BATCH - 1;

      if (hiID > endID)
        hiID = endID;

      frstOlap = nextOlap;

      Extract_Needed_Frags(G, gkpStore, loID, hiID, next_frag_list, nextOlap);
    }

    // Wait for background processing to finish

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
}









int
main(int argc, char **argv) {
  feParameters  *G = new feParameters();

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      G->gkpStorePath = argv[++arg];

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
      G->Use_Haplo_Ct = FALSE;

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

  if (G->gkpStorePath == NULL)
    err++;
  if (G->ovlStorePath == NULL)
    err++;
  if (G->numThreads == 0)
    err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s[-ehp][-d DegrThresh][-k KmerLen][-x ExcludeLen]\n", argv[0]);
    fprintf(stderr, "        [-F OlapFile][-S OlapStore][-o CorrectFile]\n");
    fprintf(stderr, "        [-t NumPThreads][-v VerboseLevel]\n");
    fprintf(stderr, "        [-V Vote_Qualify_Len]\n");
    fprintf(stderr, "          <FragStore> <lo> <hi>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Makes corrections to fragment sequence based on overlaps\n");
    fprintf(stderr, "and recomputes overlaps on corrected fragments\n");
    fprintf(stderr, "Fragments come from <FragStore> <lo> and <hi> specify\n");
    fprintf(stderr, "the range of fragments to modify\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-d   set keep flag on end of frags with less than this many olaps\n");
    fprintf(stderr, "-F   specify file of sorted overlaps to use (in the format produced\n");
    fprintf(stderr, "     by  get-olaps\n");
    fprintf(stderr, "-h   print this message\n");
    fprintf(stderr, "-k   minimum exact-match region to prevent change\n");
    fprintf(stderr, "-o   specify output file to hold correction info\n");
    fprintf(stderr, "-p   don't use haplotype counts to correct\n");
    fprintf(stderr, "-S   specify the binary overlap store containing overlaps to use\n");
    fprintf(stderr, "-t   set number of p-threads to use\n");
    fprintf(stderr, "-v   specify level of verbose outputs, higher is more\n");
    fprintf(stderr, "-V   specify number of exact match bases around an error to vote to change\n");
    fprintf(stderr, "-x   length of end of exact match to exclude in preventing change\n");

    if (G->gkpStorePath == NULL)
      fprintf(stderr, "ERROR: no gatekeeper store (-G) supplied.\n");
    if (G->ovlStorePath == NULL)
      fprintf(stderr, "ERROR: no overlap store (-O) supplied.\n");
    if (G->numThreads == 0)
      fprintf(stderr, "ERROR: number of compute threads (-t) must be larger than zero.\n");

    exit(1);
  }


  //
  //  Initialize Globals
  //

  double MAX_ERRORS = 1 + (uint32)(G->errorRate * AS_MAX_READLEN);

  Initialize_Match_Limit(G->Edit_Match_Limit, G->errorRate, MAX_ERRORS);

  for  (uint32 i = 0;  i <= AS_MAX_READLEN;  i++)
    G->Error_Bound[i] = (int)ceil(i * G->errorRate);

  //
  //
  //

  gkStore *gkpStore = gkStore::gkStore_open(G->gkpStorePath);

  if (G->bgnID < 1)
    G->bgnID = 1;

  if (gkpStore->gkStore_getNumReads() < G->endID)
    G->endID = gkpStore->gkStore_getNumReads();

  Read_Frags(G, gkpStore);
  Read_Olaps(G, gkpStore);

  //  Now sort them!

  sort(G->olaps, G->olaps + G->olapsLen);

  //fprintf (stderr, "Before Stream_Old_Frags  Num_Olaps = "F_S64"\n", Num_Olaps);

  Threaded_Stream_Old_Frags(G, gkpStore);

  //fprintf (stderr, "                   Failed overlaps = %d\n", Failed_Olaps);

  gkpStore->gkStore_close();

  //Output_Details(G);
  Output_Corrections(G);

  delete G;

  exit(0);
}

