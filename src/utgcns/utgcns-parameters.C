
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

#include "utgcns.H"

//
//  Various functions related to command line parameters.
//


////////////////////
//
//  Loads reads from a sqStore partition into 'seqReads',
//
void
cnsParameters::loadPartitionedReads(void) {

  if (seqFile == nullptr)
    return;

  //  Allocate space for the reads, and buffers to load them.

  readBuffer  *rb = new readBuffer(seqFile);
  sqRead      *rd = new sqRead;

  uint64 magc;
  uint64 vers;
  uint64 defv;

  //  Read the header.

  if (rb->readIFFobject("MAGC", magc) == false)
    fprintf(stderr, "ERROR: File '%s' isn't a utgcns seqFile: no magic number found.\n", seqFile), exit(1);

  if (magc != 0x5f5f656c69467173llu)
    fprintf(stderr, "ERROR: File '%s' isn't a utgcns seqFile: found magic 0x%016lx.\n", seqFile, magc), exit(1);

  if (rb->readIFFobject("VERS", vers) == false)
    fprintf(stderr, "ERROR: File '%s' isn't a utgcns seqFile: no file version found.\n", seqFile), exit(1);

  if (vers != 0x0000000000000001llu)
    fprintf(stderr, "ERROR: File '%s' is a utgcns seqFile, but an unsupported version %lu.\n", seqFile, vers), exit(1);

  if (rb->readIFFobject("DEFV", defv) == false)
    fprintf(stderr, "ERROR: File '%s' isn't a utgcns seqFile: no default version found.\n", seqFile), exit(1);

  sqRead_defaultVersion = (sqRead_which)defv;

  fprintf(stderr, "-- Loading %s reads from seqFile '%s'\n", toString(sqRead_defaultVersion), seqFile);

  //  Read the reads.

  while (sqStore::sqStore_loadReadFromBuffer(rb, rd) == true) {
    seqReads[rd->sqRead_readID()] = rd;

    rd = new sqRead;
  }

  delete rd;
  delete rb;
}


////////////////////
//
//  Discovers which tigs are in the selected partition, loads into
//  std:set<uint32> processList, so we can later ignore tigs not
//  in the partition.
//
void
cnsParameters::loadProcessList(void) {
  uint32             Lmax = 1024;
  uint32             Llen = 0;
  char              *L    = new char [Lmax];
  char              *N    = new char [FILENAME_MAX + 1];

  snprintf(N, FILENAME_MAX, "%s/partitioning", tigName);

  if ((tigPart > 0) &&             //  Partitioning requested, and
      (fileExists(N) == true)) {   //  partitioning file exists, load it.
    FILE *F = merylutil::openInputFile(N);

    while (merylutil::readLine(L, Llen, Lmax, F)) {
      splitToWords S(L);

      if (S.touint32(5) == tigPart)
        processList.insert(S.touint32(0));
    }

    merylutil::closeFile(F, N);
  }

  delete [] N;
  delete [] L;
}


////////////////////
//
//  Loads a COPY of the specified tig, or returns nullptr if it isn't in the
//  partition we're operating on.  The COPY is important, as it lets us treat
//  tigs from the store and tigs from the package the same way.
//
tgTig *
cnsParameters::copyTig(uint32 ti) {

  if ((tigName != nullptr) && (processList.size() == 0))
    loadProcessList();

  if ((processList.size() > 0) &&       //  Ignore tigs not in our partition.
      (processList.count(ti) == 0))     //  (if a partition exists)
    return nullptr;

  return tigStore->copyTig(ti, new tgTig);
}


////////////////////
//
//  Returns true if we should NOT process this tig according to:
//   - lack of a tig or no reads
//   - length of the tig
//   - command line selection
//      - only unassem
//      - only contigs
//      - no singletons
//      - no repeats
//      - no bubbles
//
bool
cnsParameters::skipTig(tgTig *tig) {
  int32  skip = 0;

  if (tig == nullptr)
    return true;

  tig->_utgcns_verboseLevel = verbosity;   //  Propagate verbosity to low-level algorithms.

  if (tig->numberOfChildren() == 0)
    skip = fprintf(stdout, "%7u %9u %7u  skip: no children\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if ((tig->length() < minLen) ||
      (tig->length() > maxLen))
    skip = fprintf(stdout, "%7u %9u %7u  skip: length outside range %u-%u\n", tig->tigID(), tig->length(), tig->numberOfChildren(), minLen, maxLen);

  if ((tig->tigID() < tigBgn) ||
      (tig->tigID() > tigEnd))
    skip = fprintf(stdout, "%7u %9u %7u  skip: id outside range %u-%u\n", tig->tigID(), tig->length(), tig->numberOfChildren(), tigBgn, tigEnd);

  if ((onlyUnassem == true) && (tig->_class != tgTig_unassembled))
    skip = fprintf(stdout, "%7u %9u %7u  skip: not unassembled\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if ((onlyContig  == true) && (tig->_class != tgTig_contig))
    skip = fprintf(stdout, "%7u %9u %7u  skip: not contig\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if ((noSingleton == true) && (tig->numberOfChildren() == 1))
    skip = fprintf(stdout, "%7u %9u %7u  skip: not singleton\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if ((noRepeat    == true) && (tig->_suggestRepeat == true))
    skip = fprintf(stdout, "%7u %9u %7u  skip: not repeat\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if ((noBubble    == true) && (tig->_suggestBubble == true))
    skip = fprintf(stdout, "%7u %9u %7u  skip: not bubble\n", tig->tigID(), tig->length(), tig->numberOfChildren());

  if (skip)        //  If skipped, we'll immediately next get another
    return true;   //  tig from whatever input we're reading from.

  if (tig->numberOfChildren() > 1)
    fprintf(stdout, "%7u %9u %7u",         //  The end of this line is printed
            tig->tigID(), tig->length(),   //  by processTigs() after it
            tig->numberOfChildren());      //  filters contains.

  return false;
}


////////////////////
//
//  Delete any reads we've saved.
//
void
cnsParameters::unloadReads(void) {

  for (auto it=seqReads.begin(); it != seqReads.end(); ++it)
    delete it->second;

  seqReads.clear();
}


////////////////////
//
//  Finsihed processing; release memory, close inputs and outputs.
//
void
cnsParameters::closeAndCleanup(void) {
  delete seqStore;   seqStore = nullptr;
  delete tigStore;   tigStore = nullptr;

  delete importFile; importFile = nullptr;

  unloadReads();

  merylutil::closeFile(outResultsFile, outResultsName);
  merylutil::closeFile(outLayoutsFile, outLayoutsName);

  merylutil::closeFile(outSeqFileA, outSeqNameA);
  merylutil::closeFile(outSeqFileQ, outSeqNameQ);
}
