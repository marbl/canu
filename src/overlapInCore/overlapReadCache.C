
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

#include "overlapReadCache.H"

#include <set>
#include <vector>
#include <algorithm>

using namespace std;


overlapReadCache::overlapReadCache(sqStore *seqStore_, uint64 memLimit) {
  seqStore    = seqStore_;
  nReads      = seqStore->sqStore_lastReadID();

  readAge     = new uint32 [nReads + 1];
  readLen     = new uint32 [nReads + 1];

  memset(readAge, 0, sizeof(uint32) * (nReads + 1));
  memset(readLen, 0, sizeof(uint32) * (nReads + 1));

  readSeqFwd  = new char * [nReads + 1];

  memset(readSeqFwd, 0, sizeof(char *) * (nReads + 1));

  memoryLimit = memLimit * 1024 * 1024 * 1024;
}



overlapReadCache::~overlapReadCache() {
  delete [] readAge;
  delete [] readLen;

  for (uint32 rr=0; rr<=nReads; rr++)
    delete [] readSeqFwd[rr];

  delete [] readSeqFwd;
}



void
overlapReadCache::loadRead(uint32 id) {
  seqStore->sqStore_getRead(id, &read);

  readLen[id] = read.sqRead_length();

  readSeqFwd[id] = new char [readLen[id] + 1];

  memcpy(readSeqFwd[id], read.sqRead_sequence(), sizeof(char) * readLen[id]);

  readSeqFwd[id][readLen[id]] = 0;
}



//  Make sure that the reads in 'reads' are in the cache.
//  Ideally, these are just the reads we need to load.
void
overlapReadCache::loadReads(set<uint32> reads) {
  uint32  nn = 0;
  uint32  nc = reads.size() / 25;

  //  For each read in the input set, load it.

  //if (reads.size() > 0)
  //  fprintf(stderr, "loadReads()--  Need to load %u reads.\n", reads.size());

  for (set<uint32>::iterator it=reads.begin(); it != reads.end(); ++it) {
    //if ((++nn % nc) == 0)
    //  fprintf(stderr, "loadReads()-- %6.2f%% finished.\n", 100.0 * nn / reads.size());

    if (readLen[*it] != 0)
      continue;

    loadRead(*it);
  }

  //fprintf(stderr, "loadReads()-- %6.2f%% finished.\n", 100.0);

  //  Age all the reads in the cache.

  uint32  nLoaded = 0;

  for (uint32 id=0; id<nReads; id++) {
    if (readLen[id] > 0)
      nLoaded++;
    readAge[id]++;
  }

  //fprintf(stderr, "loadReads()-- loaded %u -- %u in cache\n", reads.size(), nLoaded);
}


void
overlapReadCache::markForLoading(set<uint32> &reads, uint32 id) {

  //  Note that it was just used.
  readAge[id] = 0;

  //  Already loaded?  Done!
  if (readLen[id] != 0)
    return;

  //  Already pending?  Done!
  if (reads.count(id) != 0)
    return;

  //  Mark it for loading.
  reads.insert(id);
}



void
overlapReadCache::loadReads(ovOverlap *ovl, uint32 nOvl) {
  set<uint32>     reads;

  for (uint32 oo=0; oo<nOvl; oo++) {
    markForLoading(reads, ovl[oo].a_iid);
    markForLoading(reads, ovl[oo].b_iid);
  }

  loadReads(reads);
}



void
overlapReadCache::loadReads(tgTig *tig) {
  set<uint32>     reads;

  markForLoading(reads, tig->tigID());

  for (uint32 oo=0; oo<tig->numberOfChildren(); oo++)
    if (tig->getChild(oo)->isRead() == true)
      markForLoading(reads, tig->getChild(oo)->ident());

  loadReads(reads);
}



void
overlapReadCache::purgeReads(void) {
  uint32  maxAge     = 0;
  uint64  memoryUsed = 0;

  //  Find maxAge, and sum memory used

  for (uint32 rr=0; rr<=nReads; rr++) {
    if (maxAge < readAge[rr])
      maxAge = readAge[rr];

    memoryUsed += readLen[rr];
  }

  //  Purge oldest until memory is below watermark

  while ((memoryLimit < memoryUsed) &&
         (maxAge > 1)) {
    fprintf(stderr, "purgeReads()--  used " F_U64 "MB limit " F_U64 "MB -- purge age " F_U32 "\n", memoryUsed >> 20, memoryLimit >> 20, maxAge);

    for (uint32 rr=0; rr<=nReads; rr++) {
      if (maxAge == readAge[rr]) {
        memoryUsed -= readLen[rr];

        delete [] readSeqFwd[rr];  readSeqFwd[rr] = NULL;

        readLen[rr] = 0;
        readAge[rr] = 0;
      }
    }

    maxAge--;
  }
}

