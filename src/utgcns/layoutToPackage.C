
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

#include "runtime.H"
#include "system.H"
#include "strings.H"

#include "sqStore.H"
#include "sqCache.H"
#include "tgStore.H"

#include "unitigPartition.H"

#include <vector>
#include <map>
#include <string>
#include <fstream>



void
loadVerkkoLayouts(sqCache              *reads,
                  std::vector<tgTig *> &tigs,
                  compressedFileReader *layoutFile,
                  const char           *outputPrefix) {

  uint32        err = 0;
  splitToWords  W;
  uint32        lineLen = 0;
  uint32        lineMax = 0;
  char         *line    = nullptr;

  tgTig        *tig     = new tgTig;
  uint32        nReads  = 0;
  char          fname[FILENAME_MAX+1];

  snprintf(fname, FILENAME_MAX, "%s.tig_names.txt", outputPrefix);
  std::ofstream tigNames(fname, std::ofstream::out);
  snprintf(fname, FILENAME_MAX, "%s.readNames.txt", outputPrefix);
  std::ofstream readNames(fname, std::ofstream::out);

  //  A real tigStore doesn't have a zeroth tig, and we maintain that here by
  //  adding a nullptr to our list.
  tigs.push_back(nullptr);

  while (AS_UTL_readLine(line, lineLen, lineMax, layoutFile->file()) == true) {
    W.split(line);

    if      (nReads > 0) {
      uint32  id = reads->sqCache_mapNameToID(W[0]);

      if (id == 0) {
         fprintf(stderr, "ERROR: While processing tig %d, Read '%s' was not found.\n", tig->tigID(), W[0]);
         err++;
      }
      else {
         tig->addChild()->set(id, 0, 0, 0, strtouint32(W[1]), strtouint32(W[2]));
         readNames << W[0] << "\t" << id << "\t" << tig->tigID() << std::endl; 
      }

      nReads--;
    }

    else if (strcmp(W[0], "tig") == 0) {
#warning still need to save tigName somewhere
      tig->_tigID = tigs.size();
      tigNames << tig->tigID() << "\t" << W[1] << std::endl;
      //tig->_tigName = duplicateString(W[1]);
    }

    else if (strcmp(W[0], "len") == 0) {
      tig->_layoutLen = strtouint32(W[1]);
    }

    else if (strcmp(W[0], "rds") == 0) {
      nReads = strtouint32(W[1]);
    }

    else if (strcmp(W[0], "end") == 0) {
      if (nReads != 0) {
         fprintf(stderr, "ERROR: Tig '%d' reads doesn't match number expected\n", tig->tigID());
         err++;
      }
      fprintf(stderr, "-- Loading layouts - tig %6u of length %9u bp with %7u reads\n", tig->tigID(), tig->length(), tig->numberOfChildren());
      tigs.push_back(tig);
      tig = new tgTig;
    }

    else {
      fprintf(stderr, "ERROR: Unexpected input line %s\n", line);
      err++;
    }
  }

  delete [] line;
  delete    tig;
  tigNames.close();
  readNames.close();

  if (err > 0) {
    fprintf(stderr, "ERROR: loading layouts failed, check your input sequences and layout file!");
    fprintf(stderr, "\n");
    assert(0); 
  }
}



int
main(int argc, char **argv) {
  char const                *layoutFilename = nullptr;
  char const                *outputPrefix   = nullptr;

  std::vector<char const *>  readFilenames;

  double   partitionSize    = 1.00;   //  Size partitions to be 100% of the largest tig.
  double   partitionScaling = 1.00;   //  Estimated tig length is 100% of actual tig length.
  double   partitionReads   = 0.05;   //  5% of all reads can end up in a single partition.


  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-layout") == 0) {      //  Input ASCII layout file
      layoutFilename = argv[++arg];
    }

    else if (strcmp(argv[arg], "-reads") == 0) {       //  Input sequence files,
      while (fileExists(argv[arg+1]) == true)          //  multiple files per switch
        readFilenames.push_back(argv[++arg]);          //  supported.
    }

    else if (strcmp(argv[arg], "-output") == 0) {      //  Output package filename
      outputPrefix = argv[++arg];                      //  prefix; prefix.####.package
    }

    else if (strcmp(argv[arg], "-partition") == 0) {
      partitionSize    = strtodouble(argv[++arg]);     //  Partitioning parameters,
      partitionScaling = strtodouble(argv[++arg]);     //  the same as for utgcns.
      partitionReads   = strtodouble(argv[++arg]);
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if (layoutFilename == nullptr)
    err.push_back("ERROR:  No input layout (-layout) supplied!\n");

  if ((layoutFilename != nullptr) && (fileExists(layoutFilename) == false)) {
    char *es = new char [1024];
    snprintf(es, 1024, "ERROR:  Input layout (-layout) '%s' doesn't exist!\n", layoutFilename);
    err.push_back(es);
  }

  if (readFilenames.size() == 0)
    err.push_back("ERROR:  No input read sequence files (-reads) supplied!\n");

  for (uint32 rr=0; rr<readFilenames.size(); rr++)
    if (fileExists(readFilenames[rr]) == false) {
      char *es = new char [1024];
      snprintf(es, 1024, "ERROR:  Input read sequence file (-reads) '%s' doesn't exist!\n", readFilenames[rr]);
      err.push_back(es);
    }

  if (outputPrefix == nullptr)
    err.push_back("ERROR:  No output prefix (-output) supplied!\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  INPUT\n");
    fprintf(stderr, "    -layout l            Input layouts.\n");
    fprintf(stderr, "    -reads a [b ...]     Input reads, fasta/fasta, uncompressed/gz/bz2/xz.\n");
    fprintf(stderr, "    -output p            Output prefix.\n");
    fprintf(stderr, "    -partition s x n     Partitioning parameters.\n");
    fprintf(stderr, "                           s - max size of each partition\n");
    fprintf(stderr, "                           x - scaling\n");
    fprintf(stderr, "                           n - max number of reads per partition\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    return(1);
  }

  //  Load read sequences and build a map from read-name to read-id.

  sqCache *reads = new sqCache();

  for (char const *filename : readFilenames) {
    fprintf(stderr, "-- Loading reads from '%s'.\n", filename);
    reads->sqCache_loadReads(filename);
  }

  //
  //  Load the layout first so we can fail quickly if needed.
  //

  fprintf(stderr, "-- Loading layouts from '%s'.\n", layoutFilename);

  std::vector<tgTig *>   tigs;
  compressedFileReader  *layoutFile = new compressedFileReader(layoutFilename);

  loadVerkkoLayouts(reads, tigs, layoutFile, outputPrefix);

  delete layoutFile;

  fprintf(stderr, "-- Loading layouts from '%s': %lu layouts loaded.\n", layoutFilename, tigs.size() - 1);

  //
  //  Compute partitions.
  //

  tigPartitioning  tp;

  tp.loadTigInfo(tigs);

  tp.greedilyPartition(partitionSize,
                       partitionScaling,
                       partitionReads,
                       true);

  //
  //  Output packages for each partition.
  //

  writeBuffer  **package = new writeBuffer * [tp._nPartitions];

  fprintf(stderr, "-- Opening %u output packages.\n", tp._nPartitions);

  for (uint32 pi=0; pi<tp._nPartitions; pi++) {
    char packageName[FILENAME_MAX+1];
    snprintf(packageName, FILENAME_MAX, "%s.%04u.cnspack", outputPrefix, pi);

    package[pi] = new writeBuffer(packageName, "w");
  }

  //
  //  Output each tig to the correct package.
  //

  fprintf(stderr, "-- Creating packages with %lu tigs.\n", tigs.size() - 1);

  for (uint32 ti=0; ti<tigs.size(); ti++)
    if (tigs[ti])
      tigs[ti]->exportData(package[ tp._tigInfo[ti].partition ], reads, false);

  //
  //  Close outputs, cleanup and go home.
  //

  fprintf(stderr, "-- Closing output packages and cleaning up.\n");

  for (uint32 pi=0; pi<tp._nPartitions; pi++)
    delete package[pi];
  delete [] package;

  for (uint32 ti=0; ti<tigs.size(); ti++)
    delete tigs[ti];

  delete reads;

  fprintf(stderr, "Bye.\n");
  return(0);
}