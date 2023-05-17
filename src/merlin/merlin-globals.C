
/******************************************************************************
 *
 *  This is a k-mer based variant evaluation tool for polishing assemblies.
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
#include "strings.H"
#include <libgen.h>


#if 1
//This came from /work/meryl/src/utility/src/files/reading-v1.H
bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

  if (F == nullptr)
    return(false);

  if ((L == nullptr) || (Lmax == 0))
    allocateArray(L, Lmax, 1024);

  Llen = 0;

  int32   ch     = getc(F);
  uint32  growth = 1024;

  if (feof(F))
    return(false);

  //  Keep reading characters until EOF or a line terminator is encountered.

  while ((feof(F) == false) && (ch != '\n')) {
    if (Llen + 1 >= Lmax)
      resizeArray(L, Llen, Lmax, Lmax + growth, _raAct::copyData | _raAct::clearNew);  //  Grow the array.

    L[Llen++] = ch;

    ch = getc(F);
  }

  L[Llen] = 0;

  //  Trim trailing whitespace.

  while ((Llen > 0) && (isWhiteSpace(L[Llen-1])))
    L[--Llen] = 0;

  return(true);
}
#endif



//  Look up kmers for its existence in markerLookup
kmvalu
merlinGlobal::lookupMarker(kmer     fmer,
                           kmer     rmer) {
  kmvalu val = markerLookup->value(fmer) + markerLookup->value(rmer);
  return val;
}


//  Read probabilities lookup table for 1-4 copy kmers.
void
merlinGlobal::load_Kmetric(void) {

  if (pLookupTable == nullptr)    //  No input supplied, no lookup table to make.
    return;

  fprintf(stderr, "-- Loading probability table '%s'.\n\n", pLookupTable);

  if (fileExists(pLookupTable) == false) {
    fprintf(stderr, "ERROR: Probability table (-prob) file '%s' doesn't exist!\n", pLookupTable);
    exit(1);
  }

  compressedFileReader  F(pLookupTable);

  uint32   lineMax = 0;
  uint32   lineLen = 0;
  uint32   lineNum = 0;
  char    *line    = nullptr;

  while (readLine(line, lineLen, lineMax, F.file())) {
    splitToWords  S(line, ',');

    if (S.numWords() == 2) {
      uint32  k = S.touint32(0);
      double  p = S.todouble(1);

      copyKmerK.push_back(k);
      copyKmerP.push_back(p);

      lineNum++;

      //  fprintf(stderr, "Copy-number: %u\t\tReadK: %u\tProbability: %f\n", lineNum, k, p);
    }

    else {
      fprintf(stderr, "Copy-number: invalid line %u:  '%s'\n", lineNum, line);
    }
  }

  delete [] line;
}

//  Load meryl DBs as lookups.
//  This will be re-used for hapA and hapB
void
merlinGlobal::load_Kmers(char* dbName) {
  double            reqMemory = 0.0;

  //  Make markerDB first so we know the k size
  merylFileReader*  markerDB     = new merylFileReader(dbName);
  
  //  Since estimateMemoryUsage() is now including space for temporary
  //  buffers that are used only when loading, this estimate is significantly
  //  too large for small datasets.  If table1 and table2 need only 5 GB
  //  memory (each), the estimate for each will also include several GB for
  //  buffers (based on the number of threads); 16 threads = 8 GB buffers.
  //  So while the data needs 10 GB memory, meryl claims it needs 2x 13 GB =
  //  26 GB memory.  Since the tables are loaded sequentially, it really only
  //  needs 13 - 8 + 13 - 8 = 18 GB peak, 10 GB final.
#warning estimate is too high

  fprintf(stderr, "-- Estimating required space for loading '%s'\n", markerDBname);
  markerLookup = new merylExactLookup();
  reqMemory += markerLookup->estimateMemoryUsage(markerDB, maxMemory, 0, minV, maxV);

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Memory needed: %.3f GB\n", reqMemory);
  fprintf(stderr, "-- Memory limit:  %.3f GB\n", maxMemory);
  fprintf(stderr, "--\n");

  if (reqMemory > maxMemory) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", markerDBname);
  markerLookup->load(markerDB, maxMemory, 0, minV, maxV);

  delete markerDB;    //  Not needed anymore.
}

void
merlinGlobal::load_Sequence(void) {

  //  Open input sequence.
  if (seqName != nullptr) {
    fprintf(stderr, "-- Opening sequences in '%s'.\n", seqName);
    seqFile = new dnaSeqFile(seqName);
  } 

}



void
merlinGlobal::open_Inputs(void) {

  if (reportType == OP_BUILD) {
    // if anything needs to be processed before running the core
    return;
  }

}

void
merlinGlobal::outputGraph(void) {

  //  output what will go to debug and graph

  //  Initialize ...
  if (graphGfaFile == nullptr) {
    char  name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s.gfa", outName);
    graphGfaFile = new compressedFileWriter(name);
    fprintf(graphGfaFile->file(), "H\tVN:Z:1.2\n");
  }

  //  Output gfa. First, iterate through all the nodes to print segments.
  //  Then output all the edges as links. Nothing sophisticated.
  //  Do this for nodes in FWD and REV.
  char kmerstringF[65];
  char kmerstringT[65];

  for ( auto itr = nodes.begin(); itr != nodes.end(); ++itr ) {
    itr->second->nodeName(kmerstringF);
    fprintf(graphGfaFile->file(), "S\t%s\n", kmerstringF);  //  itr->second is the node

    //  collect edges, Fwd then Rev
    for ( int ii = 0; ii < itr->second->edgesFwd.size(); ii++ ) {
      edges.push_back(itr->second->edgesFwd.at(ii));
    }
    for ( int ii = 0; ii < itr->second->edgesRev.size(); ii++ ) {
      edges.push_back(itr->second->edgesRev.at(ii));
    }
  }

  //  Output edges (L)
  fprintf(stderr, "Writing to gfa: %d S nodes and %d L edges\n\n", nodes.size(), edges.size());

  for ( int ii = 0; ii < edges.size(); ii++ ) {
    edges[ii]->fromNode->nodeName(kmerstringF);
    edges[ii]->toNode->nodeName(kmerstringT);
    fprintf(graphGfaFile->file(), "L\t%s\t%c\t%s\t%c\t*\tRC:i:%d\n",
        kmerstringF,
        edges[ii]->orientFrom,
        kmerstringT,
        edges[ii]->orientTo,
        edges[ii]->readCnt);
  }

}

//  Assign node group with colors

void
merlinGlobal::colorGraph(const char* col) {

  //  for DEBUG only
  /*
  char kmerstringF[65];
  char kmerstringR[65];
  */

  kmvalu val = 0;

  for ( auto itr = nodes.begin(); itr != nodes.end(); ++itr ) {
    //  itr->second is the node

    //  DEBUG : check reverse mer is correct
    /*
    itr->second->nodeName(kmerstringF);
    itr->second->reverseComplement().toString(kmerstringR);
    fprintf(stderr, "[ DEBUG ] :: checking both %s and %s\n", kmerstringF, kmerstringR);
    */
    val = markerLookup->value(itr->second->nodeId) + markerLookup->value(itr->second->reverseComplement());
    if (val > 0) {
      //  set col and count
      itr->second->col=col;
      itr->second->val=val;
    }
  }
}

void
merlinGlobal::outputCsv(void) {

  if (graphCsvFile == nullptr) {
    char  name[FILENAME_MAX+1];
    snprintf(name, FILENAME_MAX, "%s.csv", outName);
    graphCsvFile = new compressedFileWriter(name);
    fprintf(graphCsvFile->file(), "node\tcount\tcolor\n");
  }

  char kmerstring[65];
  for ( auto itr = nodes.begin(); itr != nodes.end(); ++itr ) {
    itr->second->nodeName(kmerstring);
    fprintf(graphCsvFile->file(), "%s\t%u\t%s\n", kmerstring, itr->second->val, itr->second->col);
  }

}
