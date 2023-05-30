
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
std::map<uint32, sqRead *> *
loadPartitionedReads(char *seqFile) {

  if (seqFile == NULL)
    return(NULL);

  //  Allocate space for the reads, and buffers to load them.

  std::map<uint32, sqRead *>  *reads = new std::map<uint32, sqRead *>;
  readBuffer                  *rb    = new readBuffer(seqFile);
  sqRead                      *rd    = new sqRead;

  uint64 magc;
  uint64 vers;
  uint64 defv;

  //  Read the header.

  if (rb->readIFFobject("MAGC", magc) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no magic number found.\n", seqFile), exit(1);

  if (magc != 0x5f5f656c69467173llu)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: found magic 0x%016lx.\n", seqFile, magc), exit(1);

  if (rb->readIFFobject("VERS", vers) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no file version found.\n", seqFile), exit(1);

  if (vers != 0x0000000000000001llu)
    fprintf(stderr, "File '%s' is a utgcns seqFile, but an unsupported version %lu.\n", seqFile, vers), exit(1);

  if (rb->readIFFobject("DEFV", defv) == false)
    fprintf(stderr, "File '%s' isn't a utgcns seqFile: no default version found.\n", seqFile), exit(1);

  sqRead_defaultVersion = (sqRead_which)defv;

  fprintf(stderr, "Loading %s reads from seqFile '%s'\n", toString(sqRead_defaultVersion), seqFile);

  //  Read the reads.

  while (sqStore::sqStore_loadReadFromBuffer(rb, rd) == true) {
    (*reads)[rd->sqRead_readID()] = rd;

    rd = new sqRead;
  }

  delete rd;
  delete rb;

  //  Return the reads.

  return(reads);
}


////////////////////
//
//  Discovers which tigs are in the selected partition, loads into
//  std:set<uint32> processList, so we can later ignore tigs not
//  in the partition.
//
std::set<uint32>
loadProcessList(char *prefix, uint32 tigPart) {
  std::set<uint32>   processList;
  uint32             Lmax = 1024;
  uint32             Llen = 0;
  char              *L    = new char [Lmax];
  char              *N    = new char [FILENAME_MAX + 1];

  snprintf(N, FILENAME_MAX, "%s/partitioning", prefix);

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

  return(processList);
}
