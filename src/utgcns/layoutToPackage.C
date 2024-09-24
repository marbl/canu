
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
                  const char           *mapPrefix) {

  uint32        err = 0;
  splitToWords  W;
  uint32        lineLen = 0;
  uint32        lineMax = 0;
  char         *line    = nullptr;

  tgTig        *tig     = new tgTig;
  uint32        nReads  = 0;
  char          fname[FILENAME_MAX+1];

  FILE         *readMap = merylutil::openOutputFile(mapPrefix, '.', "readName_to_ID.map");
  FILE         *tigMap  = merylutil::openOutputFile(mapPrefix, '.', "tigName_to_ID.map");

  if (readMap) {
    fprintf(readMap, "    readID       tigID  read name\n");
    fprintf(readMap, "----------  ----------  --------------------\n");
  }

  if (tigMap) {
    fprintf(tigMap, "     tigID  tig name\n");
    fprintf(tigMap, "----------  --------------------\n");
  }

  //snprintf(fname, FILENAME_MAX, "%s.tig_names.txt", outputPrefix);
  //std::ofstream tigNames(fname, std::ofstream::out);
  //snprintf(fname, FILENAME_MAX, "%s.readNames.txt", outputPrefix);
  //std::ofstream readNames(fname, std::ofstream::out);

  //  A real tigStore doesn't have a zeroth tig, and we maintain that here by
  //  adding a nullptr to our list.
  tigs.push_back(nullptr);

  //  We're expecting a layout with format:
  //    tig  piece000001
  //    len  3434115
  //    rds  8738
  //    m54316_180808_005743/7668064/ccs  10831  0     200  900
  //    m54316_180808_005743/12583071/ccs  1510  11740   0  2000
  //    m54316_180808_005743/9240893/ccs   2002  12086
  //
  //  The whitespace is (probably) tabs, but we don't care.
  //
  //  The first pair of numbers is the position the trimmed read is expected
  //  to go.  The second pair is the amount to trim from the read in the
  //  orientation the read is placed in the tig (unitigConsensus.C:146).  So,
  //  the first read will be placed reversed, with 200 bp trimmed from the
  //  start (the original end) and 900 bp from the end (the original start):
  //
  //               0             10831
  //        [200bp]<------------------[900bp]
  //                  [0bp]-------------->[2000bp]
  //                          ------------->
  //
  while (merylutil::readLine(line, lineLen, lineMax, layoutFile->file()) == true) {
    W.split(line);

    //  nReads is more than zero if we've encountered a 'rds' line, and we're
    //  expecting to find that many read lines.  We first map the read name
    //  to the id of the read in our sqCache, then add the read to the tig,
    //  with or without skip information.
    //
    //  Errors are reported but the read is still processed.  We will crash
    //  asserting err == 0 after all the input is processed.

    if      (nReads > 0) {
      tgPosition *ch = tig->addChild();

      uint32      id     = reads->sqCache_mapNameToID(W[0]);                        //  First word: read ID.
      uint32      pa     = 0;                                                       //  Parent never used.
      int32       bgn    = (W.numWords() >= 3) ?  strtoint32(W[1])       : 0;       //  2nd: bgn position
      int32       end    = (W.numWords() >= 3) ?  strtoint32(W[2])       : 0;       //  3rd: end position
      bool        ignore = (W.numWords() == 4) ? (strtoint32(W[3]) != 0) : false;   //  4th: ignore this read?  (verkko)
      int32       askp   = (W.numWords() == 5) ?  strtoint32(W[3])       : 0;       //  4th: ignore bases?      (canu)
      int32       bskp   = (W.numWords() == 5) ?  strtoint32(W[4])       : 0;       //  5th: ignore bases?      (canu)

      if (id == 0)
        fprintf(stderr, "ERROR: While processing tig %d, Read '%s' was not found.\n", tig->tigID(), W[0]), err++;

      if (readMap)
        fprintf(readMap, "%10u  %10u  %s\n", id, tig->tigID(), W[0]);

      ch->set(id, pa, 0, 0, bgn, end, askp, bskp);
      ch->skipConsensus(ignore);

      if ((W.numWords() < 3) || (W.numWords() > 5))
        fprintf(stderr, "ERROR: While processing tig %d, expected 3, 4 or 5 words, got %u in line '%s'.\n", tig->tigID(), W.numWords(), line), err++;

      nReads--;
    }

    else if (strcmp(W[0], "tig") == 0) {
      tig->_tigID = tigs.size();

      if (tigMap)
        fprintf(tigMap, "%10u  %s\n", tig->tigID(), W[1]);
    }

    else if (strcmp(W[0], "len") == 0) {
      tig->_layoutLen = strtouint32(W[1]);
    }

    else if (strcmp(W[0], "rds") == 0) {
      nReads = strtouint32(W[1]);
    }

	else if (strcmp(W[0], "trm") == 0) {
	   tig->_suggestNoTrim = (strtouint32(W[1]) != 0);
    }

    else if (strcmp(W[0], "end") == 0) {
      if (nReads != 0)
        fprintf(stderr, "ERROR: Tig '%d' reads doesn't match number expected\n", tig->tigID()), err++;

      tig->cleanup();

      fprintf(stderr, "-- Loading layouts - tig %6u of length %9u bp with %7u reads\n", tig->tigID(), tig->length(), tig->numberOfChildren());
      tigs.push_back(tig);
      tig = new tgTig;
    }

    else {
      fprintf(stderr, "ERROR: Unexpected input line %s\n", line), err++;
    }
  }

  delete [] line;
  delete    tig;

  merylutil::closeFile(tigMap);
  merylutil::closeFile(readMap);

  if (err > 0) {
    fprintf(stderr, "ERROR: loading layouts failed, check your input sequences and layout file!");
    fprintf(stderr, "\n");
    assert(0); 
  }
}



int
main(int argc, char **argv) {
  char const                *layoutFilename = nullptr;
  char const                *outputPattern  = nullptr;
  char const                *mapPrefix      = nullptr;

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
      outputPattern = argv[++arg];                     //  pattern.
    }

    else if (strcmp(argv[arg], "-idmap") == 0) {       //  Output read/tig id map
      mapPrefix = argv[++arg];                         //  filename prefix.
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

  if (outputPattern == nullptr)
    err.push_back("ERROR:  No output pattern (-output) supplied!\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  INPUT\n");
    fprintf(stderr, "    -layout l                  Input layouts.\n");
    fprintf(stderr, "    -reads a [b ...]           Input reads, fasta/fasta, uncompressed/gz/bz2/xz.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -output name####.cnspack   Output file pattern.  Must include one or more consecutive '#'\n");
    fprintf(stderr, "                               symbols.  These will be replaced with the zero-leading number of\n");
    fprintf(stderr, "                               the package, expanded to fit the number of packages if needed.\n");
    fprintf(stderr, "                                 Examples:  (assume 173 total packages)\n");
    fprintf(stderr, "                                  'name.cnspack'     ->  'name.cnspack001'  through 'name.cnspack173'\n");
    fprintf(stderr, "                                  'name#.cnspack'    ->  'name001.cnspack'  through 'name173.cnspack'\n");
    fprintf(stderr, "                                  'name####.cnspack' ->  'name0001.cnspack' through 'name0173.cnspack'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -idmap prefix              Write a mapping of read/tig name <-> read/tig ID to files\n");
    fprintf(stderr, "                                 prefix.readName_to_ID.map\n");
    fprintf(stderr, "                                 prefix.tigName_to_ID.map\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -partition s x n           Partitioning parameters.\n");
    fprintf(stderr, "                                  s - max size of each partition\n");
    fprintf(stderr, "                                  x - scaling\n");
    fprintf(stderr, "                                  n - max number of reads per partition\n");
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

  loadVerkkoLayouts(reads, tigs, layoutFile, mapPrefix);

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
  //  Open packages for each partition.
  //

  writeBuffer  **package = new writeBuffer * [tp._nPartitions];

  {
    char    packageFormat[FILENAME_MAX+1];
    char    packageName[FILENAME_MAX+1];
    char    nDig = '0';
    char    nSym = '0';

    if      (tp._nPartitions < 10)      nDig = '1';   //  Count the number of
    else if (tp._nPartitions < 100)     nDig = '2';   //  digits needed.  Simple,
    else if (tp._nPartitions < 1000)    nDig = '3';   //  but ugly.
    else if (tp._nPartitions < 10000)   nDig = '4';
    else                                nDig = '5';

    char const   *pp = outputPattern;
    uint32        pi = 0;

    while ((*pp != 0) && (*pp != '#'))    //  Copy the output format up
      packageFormat[pi++] = *pp++;        //  to the first '#' symbol.

    while ((*pp != 0) && (*pp == '#'))    //  Count and skip '#' symbols.
      nSym++, pp++;

    if (nDig < nSym)                      //  If user requested more, use
      nDig = nSym;                        //  their value.

    if (nDig > '9')                       //  But never more than 9.
      nDig = '9';

    packageFormat[pi++] = '%';            //  Insert the proper format string.
    packageFormat[pi++] = '0';
    packageFormat[pi++] = nDig;
    packageFormat[pi++] = 'u';

    while (*pp != 0)                      //  Copy the rest of the string.
      packageFormat[pi++] = *pp++;
    packageFormat[pi] = 0;

    fprintf(stderr, "-- Opening %u output packages.\n", tp._nPartitions);

    for (uint32 pi=0; pi<tp._nPartitions; pi++) {
      snprintf(packageName, FILENAME_MAX, packageFormat, pi);

      fprintf(stderr, "--  '%s'\n", packageName);
      package[pi] = new writeBuffer(packageName, "w");
    }
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
