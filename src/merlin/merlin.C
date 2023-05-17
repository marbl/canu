
/******************************************************************************
 *
 *  This is a tool for linking haplotype specific k-mers in long-reads.
 *
 *  This software is based on:
 *    'Meryl'                  (https://github.com/marbl/meryl)
 *
 *  This is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "merlin-globals.H"

//  sweatShop sequence loading method.
//   - make a new merlinInput object.
//   - try to load data into it.  if there is no data left
//     in the iput, return nullptr to indicate all input
//     has been loaded.
//   - otherwise, do any initilization needed.
//   - return the new object; the sweatShop will then pass
//     the object to a compute function
//
void *
loadSequence(void *G) {
  merlinGlobal  *g = (merlinGlobal *)G;
  merlinInput   *s = new merlinInput;

  //  Try to load a new input sequence.  If it fails, destroy the merlinInput
  //  object we're trying to load and signal that we're done.

  if (g->seqFile->loadSequence(s->seq) == false) {
    delete s;
    return(nullptr);
  }

  //  Loaded something, initialize for processing.

  s->kiter.addSequence(s->seq.bases(), s->seq.length());

  //  We could, but do not, allocate space for undr and over here.  There's
  //  no gain to allocating here; we'd just be reserving lots of unused
  //  memory.

  fprintf(stderr, "Loaded sequence %s\n", s->seq.ident());
  return(s);
}


//  sweatShop compute and output functions.  All are in other files, just to
//  organize things a bit.
//

// void processBin      (void *G, void *T, void *S); // place holder
void processBuild    (void *G, void *T, void *S);
void updateGraph     (void *G, void *S);
void outputGraph     (void *G, void *S);



int
main(int32 argc, char **argv) {
  merlinGlobal  *G = new merlinGlobal(argv[0]);

  argc = AS_configure(argc, argv);

  std::vector<const char *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      G->seqName = argv[++arg];

    } else if (strcmp(argv[arg], "-marker") == 0) {
      G->markerDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-hapA") == 0) {
      G->hapADBname = argv[++arg];

    } else if (strcmp(argv[arg], "-hapB") == 0) {
      G->hapBDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      G->outName = argv[++arg];

    } else if (strcmp(argv[arg], "-peak") == 0) {
      G->peak = strtodouble(argv[++arg]); // we don't need these yet, but keeping as we will  later.

    } else if (strcmp(argv[arg], "-prob") == 0) {
      G->pLookupTable = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      G->minV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      G->maxV = strtouint64(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      G->threads = setNumThreads(argv[++arg]);   //  See comment below about setNumThreads().

    } else if (strcmp(argv[arg], "-memory") == 0) {
      G->maxMemory = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-build") == 0) { 
      G->reportType = OP_BUILD;

    } else if (strcmp(argv[arg], "-debug") == 0) {
      G->debug = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  //  Check inputs are present for the various modes.

  if (G->reportType != OP_BUILD) {
    if (G->seqName == nullptr)   err.push_back("No input sequences (-sequence) supplied.\n");
    if (G->outName == nullptr)   err.push_back("No output (-output) supplied.\n");
  }


  if  (G->reportType == OP_NONE) {
    err.push_back("No report type (-build) supplied.\n");
  }

  if (G->markerDBname == nullptr)  err.push_back("No marker meryl database (-marker) supplied.\n");


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type>              \\\n", argv[0]);
    fprintf(stderr, "         -sequence <seq.fast[a|q]>        \\\n");
    fprintf(stderr, "         -marker   <marker.meryl>         \\\n");
//  fprintf(stderr, "         -peak     <haploid_peak>         \\\n");
//  fprintf(stderr, "         -prob     <lookup_table>         \\\n");
    fprintf(stderr, "         -output   <output>               \n\n");
    fprintf(stderr, "  Lookup <markers> and build a graph given connections found in <reads>.\n");
    fprintf(stderr, "  Output files will be <output>.gfa and <output>.debug.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input -sequence can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory  m     Don't use more than m GB memory for loading mers\n");
    fprintf(stderr, "    -threads t     Multithreading for meryl lookup table construction.\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  -build\n");
    fprintf(stderr, "   Build a marker graph given connections found in <reads>.\n");
    fprintf(stderr, "   Required: -sequence, -markers and -output\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   Output: <output>.gfa : The graph built.\n");
    fprintf(stderr, "\n\n");
    fprintf(stderr, "  Optional output from -debug:\n");
    fprintf(stderr, "   <output>.THREAD_ID.debug.gz : some useful info for debugging.\n");
    fprintf(stderr, "      seqName <tab> <positions> <tab> <markers>\n");
    fprintf(stderr, "      seqName   - name of the sequence\n");
    fprintf(stderr, "      positions - first base position in the sequence where the markers are found\n");
    fprintf(stderr, "                  delemeted with \":\".\n");
    fprintf(stderr, "      markers   - marker k-mer found\n");
    fprintf(stderr, "                  delemeted with \":\".\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  //  Open read kmers, build a lookup table.

  // G->load_Kmetric();  // load prob. table
  G->load_Kmers(G->markerDBname);    // load lookup tables: markerDB
  G->load_Sequence();

  //  Configure the sweatShop.

  sweatShop       *ss = nullptr;
  merlinThrData   *td = new merlinThrData [G->threads];

  //  Check report type.
  //  the thread limit is set by setNumberOfWorkers() below.

  //  BUILD graph
  if (G->reportType == OP_BUILD) {
    fprintf(stderr, "-- Build marker graph to '%s'.gfa.\n", G->outName);
    ss = new sweatShop(loadSequence, processBuild, updateGraph);
    ss->setInOrderOutput(false);
  }

  //  Run the compute.

  if (ss) {
    ss->setNumberOfWorkers(G->threads);

    for (uint32 tt=0; tt<G->threads; tt++) {
      td[tt].threadID = tt;
      ss->setThreadData(tt, td + tt);
    }

    ss->setLoaderBatchSize(1);
    ss->setLoaderQueueSize(G->threads * 2);
    ss->setWorkerBatchSize(1);
    ss->setWriterQueueSize(16384);

    ss->run(G, false);
    
  }

  //  Output graph
  if (G->reportType == OP_BUILD) {
    G->outputGraph(); //  output entire .gfa

    if (G->hapADBname != nullptr) {
      G->load_Kmers(G->hapADBname);
      G->colorGraph("#FF8888");  //  peach
    }

    if (G->hapBDBname != nullptr) {
      G->load_Kmers(G->hapBDBname);
      G->colorGraph("#8888FF");  // purple
    }

    G->outputCsv();             //  writes the colored .csv file

  }

  //  Cleanup and quit.
  delete [] td;
  delete    ss;

  delete    G;

  fprintf(stderr, "Bye!\n");
  return(0);
}
