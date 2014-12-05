#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "sweatShop.H"

#include "libmeryl.H"

#include <algorithm>

using namespace std;

//  var, old, new -- returns true if "(var == old) and var <- new"
//
//  CAS - #elif (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) > 40100

const uint32 pileMax     = 32768;

const uint32 kmerSize    = 22;
const uint32 kmerBits    = 2 * kmerSize;

const uint32 pilePreSize = 6;
const uint32 pilePreBits = 2 * pilePreSize;

const uint32 sortPreSize = 10;
const uint32 sortPreBits = 2 * sortPreSize;


class kmerPile {
public:
  kmerPile(uint32 prefix) {
    pileLen    = 0;
    pilePrefix = prefix;
  };
  ~kmerPile() {
  };

  void     initialize(uint32 prefix) {
    pileLen    = 0;
    pilePrefix = prefix;
  };

  void     addMer(uint64 mer) {
    pileDat[pileLen++] = mer;
  };

  void     sort(void) {
    ::sort(pileDat, pileDat + pileLen);
  };

  uint32   pileLen;
  uint32   pilePrefix;

  uint64   pileDat[pileMax];
};




class kmerSorter {
public:
  kmerSorter() {
    sorterLocked = 0;
    sorterLen    = 0;
    sorterMax    = 4;
    sorterMer    = new uint64 [sorterMax];
    sorterCnt    = new uint32 [sorterMax];
  };
  ~kmerSorter() {
    delete [] sorterMer;
    delete [] sorterCnt;
  };

  void  merge(uint64 *pileDat, uint32 pileLen) {
    uint32   nmax = MAX(16, sorterLen + pileLen / 4);
    uint64  *nmer = new uint64 [nmax];
    uint32  *ncnt = new uint32 [nmax];
    uint32   npos = 0;

    assert(nmax > 0);

    uint32   spos = 0;
    uint32   ppos = 0;

    bool     useSorterFirst = false;

    if ((sorterLen > 0) && (pileLen > 0)) {
      useSorterFirst = (sorterMer[0] < pileDat[0]);

    } else if (spos < sorterLen) {
      useSorterFirst = true;

    } else if (ppos < pileLen) {
      useSorterFirst = false;

    } else {
      assert(0);
    }

    if (useSorterFirst) {
      nmer[0] = sorterMer[spos];
      ncnt[0] = sorterCnt[spos];
      spos++;
    } else {
      nmer[0] = pileDat[ppos];
      ncnt[0] = 1;
      ppos++;
    }

    while ((spos < sorterLen) && (ppos < pileLen)) {

      if (nmax <= npos + 1) {
        nmax += (pileLen - ppos) + (sorterLen - spos) + 1;

        uint64 *nmermore = new uint64 [nmax];
        uint32 *ncntmore = new uint32 [nmax];

        memcpy(nmermore, nmer, sizeof(uint64) * (npos + 1));
        memcpy(ncntmore, ncnt, sizeof(uint32) * (npos + 1));

        delete [] nmer;  nmer = nmermore;
        delete [] ncnt;  ncnt = ncntmore;
      }

      if        (nmer[npos] == sorterMer[spos]) {
        ncnt[npos] += sorterCnt[spos];
        spos++;

      } else if (nmer[npos] == pileDat[ppos]) {
        ncnt[npos] += 1;
        ppos++;

      } else if (sorterMer[spos] < pileDat[ppos]) {
        npos++;
        nmer[npos] = sorterMer[spos];
        ncnt[npos] = sorterCnt[spos];
        spos++;

      } else {
        npos++;
        nmer[npos] = pileDat[ppos];
        ncnt[npos] = 1;
        ppos++;
      }
    }

    uint32 remain = (sorterLen - spos) + (pileLen - ppos);

    if (nmax < npos + 1 + remain) {
      nmax = npos + 1 + remain;

      uint64 *nmermore = new uint64 [nmax];
      uint32 *ncntmore = new uint32 [nmax];

      memcpy(nmermore, nmer, sizeof(uint64) * (npos + 1));
      memcpy(ncntmore, ncnt, sizeof(uint32) * (npos + 1));

      delete [] nmer;  nmer = nmermore;
      delete [] ncnt;  ncnt = ncntmore;
    }



    while (spos < sorterLen) {
      if        (nmer[npos] == sorterMer[spos]) {
        ncnt[npos] += sorterCnt[spos];
      } else {
        npos++;
        nmer[npos] = sorterMer[spos];
        ncnt[npos] = sorterCnt[spos];
      }

      spos++;
    }


    while (ppos < pileLen) {
      if (nmer[npos] == pileDat[ppos]) {
        ncnt[npos] += 1;
      } else {
        npos++;
        nmer[npos] = pileDat[ppos];
        ncnt[npos] = 1;
      }

      ppos++;
    }

    delete [] sorterMer;
    delete [] sorterCnt;

    sorterMer = nmer;
    sorterCnt = ncnt;
    sorterLen = npos + 1;
    sorterMax = nmax;

#if 1
    bool broken = false;

    for (uint32 i=1; i<sorterLen; i++) {
      assert(sorterMer[i-1] < sorterMer[i]);
      if (sorterMer[i-1] >= sorterMer[i])
        broken = true;
    }
#endif

  };

  void     write(uint32 prefix, FILE *F, merylStreamWriter *W) {
    char   km[64] = {0};
    uint32 kp = pilePreSize;
    uint32 np = 0;

    {
      uint32 pre = prefix;

      for (uint32 pp=0; pp<pilePreSize; pp++) {
        km[--kp] = bitsToLetter[pre & 0x03];
        pre >>= 2;
      }
    }

    np = kmerSize - pilePreSize;

    for (uint32 ii=0; ii<sorterLen; ii++) {
      uint64 mer = sorterMer[ii];

      kp = kmerSize;

      for (uint32 pp=0; pp<np; pp++) {
        km[--kp] = bitsToLetter[mer & 0x03];
        mer >>= 2;
      }

      fprintf(F, ">"uint32FMT"\n%s\n", sorterCnt[ii], km);

      if (W)
        W->addMer(prefix,        pilePreBits,
                  sorterMer[ii], kmerBits - pilePreBits,
                  sorterCnt[ii],
                  0L);
    }
  };

  volatile uint32   sorterLocked;
  uint32   sorterLen;
  uint32   sorterMax;

  uint64  *sorterMer;
  uint32  *sorterCnt;
};




class kmerGlobal {
public:
  kmerGlobal() {
    inName    = NULL;
    inFile    = NULL;

#if 0
    inputBufferMax  = 131072;
    inputBufferLen  = 0;
    inputBufferPos  = 0;
    inputBuffer     = new char [inputBufferMax];
#endif

    inputBufferMax  = 0;
    inputBufferLen  = 0;
    inputBufferPos  = 0;
    inputBuffer     = NULL;

    outPrefix = NULL;
    outFile   = NULL;

    fkPre     = 0;
    fkMer     = 0;

    rkPre     = 0;
    rkMer     = 0;

    kLen      = 0;

    pilesFreeLock = 0;
    pilesFreeLen  = 2048;
    pilesFreeMax  = 2 << pilePreBits;
    pilesFree     = new kmerPile * [pilesFreeMax];

    memset(pilesFree, 0, sizeof(kmerPile *) * pilesFreeMax);

    piles     = new kmerPile * [1 << pilePreBits];
    sorters   = new kmerSorter [1 << sortPreBits];

    memset(piles, 0, sizeof(kmerPile *) * (1 << pilePreBits));

    for (uint32 i=0; i<pilesFreeLen; i++)
      pilesFree[i] = new kmerPile(0);

    for (uint32 i=0; i< (1 << pilePreBits); i++)
      piles[i] = new kmerPile(i);

    pilesToSortLen  = 0;
    pilesToSortMax  = 2 * (1 << pilePreBits);
    pilesToSort     = new kmerPile * [pilesToSortMax];
  };
  ~kmerGlobal() {
    delete [] piles;
    delete [] sorters;
    delete [] pilesToSort;
    //delete [] inputBuffer;
  };


  void    initialize(void) {
    //inBuffer = new readBuffer(inName, 0);

#if 1
    inputBufferMax = 0;
    inputBufferLen = 0;
    inputBufferPos = 0;
    inputBuffer    = (char *)mapFile(inName, &inputBufferLen, 'r');
#endif

    naptime.tv_sec      = 0;
    naptime.tv_nsec     = 166666666ULL;  //  1/6 second
    naptime.tv_nsec     = 250000ULL;
  };

  kmerPile *getFreePile(uint32 prefix) {
    kmerPile *pp;

    while (__sync_bool_compare_and_swap(&pilesFreeLock, 0, 1) == false)
      nanosleep(&naptime, 0L);

    assert(pilesFreeLock == 1);

    if (pilesFreeLen == 0) {
      pilesFreeLock = 0;
      //fprintf(stderr, "ALLOCATE PILE!\n");
      pp = new kmerPile(prefix);

    } else {
      pp = pilesFree[--pilesFreeLen];
      pilesFreeLock = 0;
    }

    pp->initialize(prefix);

    return(pp);
  };

  void    releasePile(kmerPile *pile) {

    if (pilesFreeLen >= pilesFreeMax) {
      //fprintf(stderr, "DELETE PILE!\n");
      delete pile;

    } else {
      while (__sync_bool_compare_and_swap(&pilesFreeLock, 0, 1) == false)
        nanosleep(&naptime, 0L);

      assert(pilesFreeLock == 1);

      pilesFree[pilesFreeLen++] = pile;

      pilesFreeLock = 0;
    }

  };


  void    addToPile(uint64 pre, uint64 mer) {

    assert(piles[pre] != NULL);
    //if (piles[pre] == NULL)
    //  piles[pre] = getFreePile(pre);

    if (piles[pre]->pileLen < pileMax) {
      piles[pre]->addMer(mer);
      return;
    }

    if (pilesToSortMax <= pilesToSortLen) {
      fprintf(stderr, "realloc\n");
      exit(1);
    }

    pilesToSort[pilesToSortLen++] = piles[pre];

    piles[pre] = getFreePile(pre);
    piles[pre]->addMer(mer);
  };


  kmerPile *getFullPile(void) {
    if (pilesToSortLen == 0)
      return(NULL);

    //fprintf(stderr, "return pile "uint32FMT"\n", pilesToSort[pilesToSortLen-1]->pilePrefix);
    return(pilesToSort[--pilesToSortLen]);
  };


  kmerPile *allDataLoaded(void) {

    for (uint32 pp=0; pp < (1 << pilePreBits); pp++) {
      if ((piles[pp] != NULL) &&
          (piles[pp]->pileLen > 0)) {
        //fprintf(stderr, "Add pile "uint32FMT" to list.\n", pp);
        pilesToSort[pilesToSortLen++] = piles[pp];
      } else {
        delete piles[pp];
      }

      piles[pp] = NULL;
    }

    fprintf(stderr, "allDataLoaded()-- pilesToSortLen = "uint32FMT"\n", pilesToSortLen);

    return(getFullPile());
  };


  void    addBases(uint32 bgn, uint32 len) {
    uint32   kp2 = kmerBits - pilePreBits - 2;
    uint32   pp2 = pilePreBits - 2;

    uint64   mpp = uint64MASK(pilePreBits);
    uint64   mkp = uint64MASK(kmerBits - pilePreBits);

    for (uint32 pos=0; pos<len; pos++) {
      uint64  bt = letterToBits[ inputBuffer[bgn+pos] ];

      if (bt > 4) {
        kLen = 0;
        continue;
      }

      uint64 tm = 0;

      tm  = fkMer >> kp2;
      tm &= 0x00000003;

      fkPre <<= 2;
      fkPre  |= tm;

      fkMer <<= 2;
      fkMer  |= bt;

      tm  = rkMer & 0x00000003;

      rkPre >>= 2;
      rkPre  |= tm << pp2;

      rkMer >>= 2;
      rkMer  |= bt << kp2;

      kLen++;

      if (kLen < kmerSize)
        continue;

      kLen = kmerSize;

      fkPre  &= mpp;
      fkMer  &= mkp;

      rkPre  &= mpp;
      rkMer  &= mkp;

      addToPile(fkPre, fkMer);
      addToPile(rkPre, rkMer);
    }
  }

  bool    addBaseToKmer(char base) {
    uint64  bt = letterToBits[base];

    if (bt > 4) {
      kLen = 0;
      return(false);
    }

    uint64 tm = 0;

    tm  = fkMer >> (kmerBits - pilePreBits - 2);
    tm &= 0x00000003;

    fkPre <<= 2;
    fkPre  |= tm;

    fkMer <<= 2;
    fkMer  |= bt;

    tm  = rkMer & 0x00000003;

    rkPre >>= 2;
    rkPre  |= tm << (pilePreBits - 2);

    rkMer >>= 2;
    rkMer  |= bt << (kmerBits - pilePreBits - 2);

    kLen++;

    if (kLen < kmerSize) {
      return(false);
    }

    kLen = kmerSize;

    fkPre  &= uint64MASK(pilePreBits);
    fkMer  &= uint64MASK(kmerBits - pilePreBits);

    rkPre  &= uint64MASK(pilePreBits);
    rkMer  &= uint64MASK(kmerBits - pilePreBits);

    addToPile(fkPre, fkMer);
    addToPile(rkPre, rkMer);

    return(true);
  };


  void  write(void) {
    char outName[FILENAME_MAX];

    sprintf(outName, "%s.fasta", outPrefix);

    errno = 0;
    FILE                *F = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", outName, strerror(errno)), exit(1);

    //merylStreamWriter   *W = new merylStreamWriter(outPrefix, kmerSize, 0, sortPreBits, false);

    for (uint32 ss=0; ss < (1 << sortPreBits); ss++)
      sorters[ss].write(ss, F, NULL);

    fclose(F);
    //delete W;
  }

  char  *inName;
  FILE  *inFile;

  readBuffer  *inBuffer;

  uint64       inputBufferMax;
  uint64       inputBufferLen;
  uint64       inputBufferPos;
  char        *inputBuffer;

  char  *outPrefix;
  FILE  *outFile;

  uint64       fkPre;  //  Forward loaded kmer
  uint64       fkMer;

  uint64       rkPre;  //  Reverse loaded kmer
  uint64       rkMer;

  uint32       kLen;

  uint32       pilesFreeLock;
  uint32       pilesFreeLen;
  uint32       pilesFreeMax;
  kmerPile   **pilesFree;

  kmerPile   **piles;
  kmerSorter  *sorters;

  struct timespec   naptime;

  uint32       pilesToSortLen;
  uint32       pilesToSortMax;
  kmerPile   **pilesToSort;
};





uint64        bytesLoaded = 0;
uint64        basesLoaded = 0;
speedCounter  bytes(" bytes %8.0f Mbytes (%8.5f Mbytes/sec\r", 1048576, 1048576, true);

  //  Reads input, constructs kmers, adds kmers to piles of kmers.
void*
sifterThread(void *global) {
  kmerGlobal *glob = (kmerGlobal *)global;
  kmerPile   *pile = glob->getFullPile();

  if (pile)
    return(pile);

  //if ((glob->inFile == NULL) && (glob->inBuffer == NULL))
  //  return(NULL);

 anotherBase:
  //bytesLoaded++;
  //if ((bytesLoaded % (16 * 1048576)) == 0)
  //  fprintf(stderr, "sifterThread()-- loaded "uint64FMT" MB.\n", bytesLoaded >> 20);

#if 0
  //  Uses the readBuffer in char-by-char mode
  //
  char ch = glob->inBuffer->read();
  bytes.tick();

  if (glob->inBuffer->eof()) {
    delete glob->inBuffer;
    glob->inBuffer = NULL;
    return(glob->allDataLoaded());
  }

  if (glob->addBaseToKmer(ch) == false)
    goto anotherBase;

#endif


#if 0
  //  Uses the readBuffer in block-copy mode
  //
  uint32  len = glob->inBuffer->read(glob->inputBuffer, glob->inputBufferMax);

  if (len == 0) {
    delete glob->inBuffer;
    glob->inBuffer = NULL;
    return(glob->allDataLoaded());
  }

  glob->addBases(0, len);
  bytes.tick(len);

#endif

#if 1
  //  Uses a direct mmap'd file
  //
  uint64 len = glob->inputBufferLen - glob->inputBufferPos;

  if (len == 0)
    return(NULL);

  if (len > 16 * 1048576)
    len = 16 * 1048576;

  //fprintf(stderr, "Add "uint64FMT" bases.\n", len);

  glob->addBases(glob->inputBufferPos, len);
  bytes.tick(len);

  glob->inputBufferPos += len;
#endif

  pile = glob->getFullPile();

  if (pile == NULL)
    goto anotherBase;

  return(pile);
}



  //  Takes a pile of kmers, sorts it, and them merges into the appropriate kmerSorter objects.
void
sorterThread(void *global, void *thread, void *thing) {
  kmerGlobal *glob = (kmerGlobal *)global;
  kmerPile   *pile = (kmerPile   *)thing;

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 166666666ULL;  //  1/6 second
  naptime.tv_nsec     = 250000ULL;

  if (pile->pileLen == 0)
    //  Nothing to add.
    return;

  pile->sort();

  uint32  pileBgn       = 0;
  uint32  pileEnd       = 1;

  uint32  pileMaskShift = sortPreBits - pilePreBits;
  uint32  pileDataShift = kmerBits - sortPreBits;

  uint64  pileToSortPreMask = uint64MASK(sortPreBits - pilePreBits) << (kmerBits - sortPreBits);
  uint64  pileToSortMask    = uint64MASK(kmerBits - sortPreBits);

  uint32  sortPre           = 0;
  uint64  pileToSort        = 0;

  while (pileBgn < pile->pileLen) {
    sortPre    = (pile->pilePrefix << pileMaskShift) | (pile->pileDat[pileBgn] >> pileDataShift);
    pileToSort = pile->pileDat[pileBgn] & pileToSortPreMask;

    //fprintf(stderr, "0x"uint64HEX"\n", pileToSortPreMask);
    //fprintf(stderr, "0x"uint64HEX"\n", pileToSortMask);

    while ((pileEnd < pile->pileLen) &&
           ((pile->pileDat[pileEnd] & pileToSortPreMask) == pileToSort)) {
      //fprintf(stderr, "0x"uint64HEX" -> 0x"uint64HEX" "uint64FMT"\n",
      //        pile->pileDat[pileEnd],
      //        pile->pileDat[pileEnd] & pileToSortMask,
      //        pile->pileDat[pileEnd] & pileToSortMask);
      pile->pileDat[pileEnd] &= pileToSortMask;
      pileEnd++;
    }

    while (__sync_bool_compare_and_swap(&glob->sorters[sortPre].sorterLocked, 0, 1) == false)
      nanosleep(&naptime, 0L);

    assert(glob->sorters[sortPre].sorterLocked == 1);

    glob->sorters[sortPre].merge(pile->pileDat + pileBgn, pileEnd - pileBgn);

    glob->sorters[sortPre].sorterLocked = 0;

    pileBgn = pileEnd;
  }
}



//  Does nothing but delete the pile object.  We don't output till the end.
void
nullThread(void *global, void *thing) {
  kmerGlobal *glob = (kmerGlobal *)global;
  kmerPile   *pile = (kmerPile   *)thing;

  glob->releasePile(pile);
}



int
main(int argc, char **argv) {
  kmerGlobal  *kg = new kmerGlobal;

  int arg=1;
  int err=0;

  while (arg < argc) {
    if      (strcmp(argv[arg], "-i") == 0)
      kg->inName = argv[++arg];

    else if (strcmp(argv[arg], "-o") == 0)
      kg->outPrefix = argv[++arg];

    else
      err++;

    arg++;
  }
  if (kg->inName == NULL)
    err++;
  if (kg->outPrefix == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -i in.sequence -i prefix\n", argv[0]);
    exit(1);
  }

  kg->initialize();

  sweatShop  *ss = new sweatShop(sifterThread,
                                 sorterThread,
                                 nullThread);

  ss->setLoaderBatchSize(512);

  ss->setNumberOfWorkers(1);
  ss->setWriterQueueSize(16384);

  //for (uint32 i=0; i<config._numSearchThreads; i++)
  //  ss->setThreadData(i, new searcherState(i));

  ss->run(kg, true);

  delete ss;

  kg->write();

  delete kg;

  exit(0);
}
