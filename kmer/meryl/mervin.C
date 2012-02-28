#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "bio++.H"
#include "sweatShop.H"

#include "libmeryl.H"

#include <algorithm>

using namespace std;


const u32bit pileMax     = 4096;

const u32bit kmerSize    = 22;
const u32bit kmerBits    = 2 * kmerSize;

const u32bit pilePreSize = 8;
const u32bit pilePreBits = 2 * pilePreSize;

pthread_mutex_t        sorterLockMutex;


//const u32bit kmerPreSize = 12;
//const u32bit kmerPreBits = 2 * kmerPreSize;


//  There are 4^pilePrefix of these allocated at a cost of
//    16 + 8 * pileMax bytes each  (4gb with max=8192 and prefix=8)
class kmerPile {
public:
  kmerPile(u32bit prefix) {
    pileLen    = 0;
    pilePrefix = prefix;
  };
  ~kmerPile() {
  };

  void     addMer(u64bit mer) {
    pileDat[pileLen++] = mer;
    //fprintf(stderr, "Add mer 0x"u64bitHEX" to pile "u32bitFMT" now of length "u32bitFMT"\n",
    //        mer, pilePrefix, pileLen);
  };

  void     sort(void) {
    ::sort(pileDat, pileDat + pileLen);
  };

  u32bit   pileLen;
  u32bit   pilePrefix;

  u64bit   pileDat[pileMax];
};




class kmerSorter {
public:
  kmerSorter() {
    sorterLocked = 0;
    sorterLen    = 0;
    sorterMax    = 1024;
    sorterMer    = new u64bit [sorterMax];
    sorterCnt    = new u32bit [sorterMax];
  };
  ~kmerSorter() {
    delete [] sorterMer;
    delete [] sorterCnt;
  };

  void  merge(kmerPile *pile) {

    if (pile->pileLen == 0)
      //  Nothing to add.
      return;

    pile->sort();

    u32bit   nmax = MAX(16, sorterLen + pile->pileLen / 4);
    u64bit  *nmer = new u64bit [nmax];
    u32bit  *ncnt = new u32bit [nmax];
    u32bit   npos = 0;

    assert(nmax > 0);

    u32bit   spos = 0;
    u32bit   ppos = 0;

    bool     useSorterFirst = false;

    if ((sorterLen > 0) && (pile->pileLen > 0)) {
      useSorterFirst = (sorterMer[0] < pile->pileDat[0]);

    } else if (spos < sorterLen) {
      useSorterFirst = true;

    } else if (ppos < pile->pileLen) {
      useSorterFirst = false;

    } else {
      assert(0);
    }

    if (useSorterFirst) {
      nmer[0] = sorterMer[spos];
      ncnt[0] = sorterCnt[spos];
      spos++;
    } else {
      nmer[0] = pile->pileDat[ppos];
      ncnt[0] = 1;
      ppos++;
    }

    while ((spos < sorterLen) && (ppos < pile->pileLen)) {

      if (nmax <= npos + 1) {
        nmax += (pile->pileLen - ppos) + (sorterLen - spos) + 1;

        u64bit *nmermore = new u64bit [nmax];
        u32bit *ncntmore = new u32bit [nmax];

        memcpy(nmermore, nmer, sizeof(u64bit) * (npos + 1));
        memcpy(ncntmore, ncnt, sizeof(u32bit) * (npos + 1));

        delete [] nmer;  nmer = nmermore;
        delete [] ncnt;  ncnt = ncntmore;
      }

      if        (nmer[npos] == sorterMer[spos]) {
        ncnt[npos] += sorterCnt[spos];
        spos++;

      } else if (nmer[npos] == pile->pileDat[ppos]) {
        ncnt[npos] += 1;
        ppos++;

      } else if (sorterMer[spos] < pile->pileDat[ppos]) {
        npos++;
        nmer[npos] = sorterMer[spos];
        ncnt[npos] = sorterCnt[spos];
        spos++;

      } else {
        npos++;
        nmer[npos] = pile->pileDat[ppos];
        ncnt[npos] = 1;
        ppos++;
      }
    }

    u32bit remain = (sorterLen - spos) + (pile->pileLen - ppos);

    if (nmax < npos + 1 + remain) {
      nmax = npos + 1 + remain;

      u64bit *nmermore = new u64bit [nmax];
      u32bit *ncntmore = new u32bit [nmax];

      memcpy(nmermore, nmer, sizeof(u64bit) * (npos + 1));
      memcpy(ncntmore, ncnt, sizeof(u32bit) * (npos + 1));

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


    while (ppos < pile->pileLen) {
      if (nmer[npos] == pile->pileDat[ppos]) {
        ncnt[npos] += 1;
      } else {
        npos++;
        nmer[npos] = pile->pileDat[ppos];
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

    for (u32bit i=1; i<sorterLen; i++) {
      assert(sorterMer[i-1] < sorterMer[i]);
      if (sorterMer[i-1] >= sorterMer[i])
        broken = true;
    }
#endif

  };

  void     write(u32bit prefix, FILE *F, merylStreamWriter *W) {
    char   km[64] = {0};
    u32bit kp = pilePreSize;
    u32bit np = 0;

    {
      u32bit pre = prefix;

      for (u32bit pp=0; pp<pilePreSize; pp++) {
        km[--kp] = bitsToLetter[pre & 0x03];
        pre >>= 2;
      }
    }

    np = kmerSize - pilePreSize;

    for (u32bit ii=0; ii<sorterLen; ii++) {
      u64bit mer = sorterMer[ii];

      kp = kmerSize;

      for (u32bit pp=0; pp<np; pp++) {
        km[--kp] = bitsToLetter[mer & 0x03];
        mer >>= 2;
      }

      fprintf(F, ">"u32bitFMT"\n%s\n", sorterCnt[ii], km);

      W->addMer(prefix,        pilePreBits,
                sorterMer[ii], kmerBits - pilePreBits,
                sorterCnt[ii],
                0L);
    }
  };

  u32bit   sorterLocked;
  u32bit   sorterLen;
  u32bit   sorterMax;

  u64bit  *sorterMer;
  u32bit  *sorterCnt;
};




class kmerGlobal {
public:
  kmerGlobal() {
    inName    = NULL;
    inFile    = NULL;

    outPrefix = NULL;
    outFile   = NULL;

    fkPre     = 0;
    fkMer     = 0;

    rkPre     = 0;
    rkMer     = 0;

    kLen      = 0;

    piles     = new kmerPile * [1 << pilePreBits];
    sorters   = new kmerSorter [1 << pilePreBits];

    memset(piles, 0, sizeof(kmerPile *) * (1 << pilePreBits));

    pilesToSortLen  = 0;
    pilesToSortMax  = 2 * (1 << pilePreBits);
    pilesToSort     = new kmerPile * [pilesToSortMax];
  };
  ~kmerGlobal() {
    delete [] piles;
    delete [] sorters;
    delete [] pilesToSort;
  };


  void    initialize(void) {
#if 0
    if (strcmp(inName, "-") == 0)
      inFile = stdin;
    else
      inFile = fopen(inName, "r");
#endif
    inBuffer = new readBuffer(inName, 128 * 1024);


    int err = pthread_mutex_init(&sorterLockMutex, NULL);
    if (err)
      fprintf(stderr, "sweatShop::run()--  Failed to configure pthreads (state mutex): %s.\n", strerror(err)), exit(1);
  };


  void    addToPile(u64bit pre, u64bit mer) {

    if (piles[pre] == NULL)
      piles[pre] = new kmerPile(pre);

    if (piles[pre]->pileLen >= pileMax) {
      //fprintf(stderr, "PILE TO LIST\n");
      if (pilesToSortMax <= pilesToSortLen) {
        fprintf(stderr, "realloc\n");
        exit(1);
      }

      pilesToSort[pilesToSortLen++] = piles[pre];
      piles[pre] = new kmerPile(pre);
    }

    piles[pre]->addMer(mer);
  };


  kmerPile *getPile(void) {
    if (pilesToSortLen == 0)
      return(NULL);

    //fprintf(stderr, "return pile "u32bitFMT"\n", pilesToSort[pilesToSortLen-1]->pilePrefix);
    return(pilesToSort[--pilesToSortLen]);
  };


  kmerPile *allDataLoaded(void) {

    for (u32bit pp=0; pp < (1 << pilePreBits); pp++) {
      if ((piles[pp] != NULL) &&
          (piles[pp]->pileLen > 0)) {
        //fprintf(stderr, "Add pile "u32bitFMT" to list.\n", pp);
        pilesToSort[pilesToSortLen++] = piles[pp];
      } else {
        delete piles[pp];
      }

      piles[pp] = NULL;
    }

    fprintf(stderr, "allDataLoaded()-- pilesToSortLen = "u32bitFMT"\n", pilesToSortLen);

    return(getPile());
  };


  bool    addBaseToKmer(char base) {
    u64bit  bt = letterToBits[base];

    if (bt > 4) {
      kLen = 0;
      //fprintf(stderr, "addBaseToKmer()-- invalid base\n");
      return(false);
    }

    u64bit tm = 0;

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
      //fprintf(stderr, "addBaseToKmer()-- kLen "u32bitFMT"\n", kLen);
      return(false);
    }

    kLen = kmerSize;

    fkPre  &= u64bitMASK(pilePreBits);
    fkMer  &= u64bitMASK(kmerBits - pilePreBits);

    rkPre  &= u64bitMASK(pilePreBits);
    rkMer  &= u64bitMASK(kmerBits - pilePreBits);

    addToPile(fkPre, fkMer);
    //addToPile(rkPre, rkMer);

    //fprintf(stderr, "addBaseToKmer()-- added to pile.\n");
    return(true);
  };


  void  write(void) {
    char outName[FILENAME_MAX];

    sprintf(outName, "%s.fasta", outPrefix);

    errno = 0;
    FILE                *F = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", outName, strerror(errno)), exit(1);

    merylStreamWriter   *W = new merylStreamWriter(outPrefix, kmerSize, 0, pilePreBits, false);

    for (u32bit ss=0; ss < (1 << pilePreBits); ss++)
      sorters[ss].write(ss, F, W);

    fclose(F);
    delete W;
  }

  char  *inName;
  FILE  *inFile;

  readBuffer  *inBuffer;

  char  *outPrefix;
  FILE  *outFile;

  u64bit       fkPre;  //  Forward loaded kmer
  u64bit       fkMer;

  u64bit       rkPre;  //  Reverse loaded kmer
  u64bit       rkMer;

  u32bit       kLen;

  kmerPile   **piles;
  kmerSorter  *sorters;

  u32bit       pilesToSortLen;
  u32bit       pilesToSortMax;
  kmerPile   **pilesToSort;
};





u64bit   bytesLoaded = 0;
u64bit   basesLoaded = 0;

  //  Reads input, constructs kmers, adds kmers to piles of kmers.
void*
sifterThread(void *global) {
  kmerGlobal *glob = (kmerGlobal *)global;
  kmerPile   *pile = glob->getPile();

  if (pile)
    return(pile);

  if ((glob->inFile == NULL) && (glob->inBuffer == NULL))
    return(NULL);

 anotherBase:
  bytesLoaded++;
  if ((bytesLoaded % (16 * 1048576)) == 0)
    fprintf(stderr, "sifterThread()-- loaded "u64bitFMT" MB.\n", bytesLoaded >> 20);

#if 0
  char    ch = getc(glob->inFile);

  if (feof(glob->inFile)) {
    fclose(glob->inFile);
    glob->inFile = NULL;
    return(glob->allDataLoaded());
  }
#else
  char ch = glob->inBuffer->read();

  if (glob->inBuffer->eof()) {
    delete glob->inBuffer;
    glob->inBuffer = NULL;
    return(glob->allDataLoaded());
  }

#endif

  if (glob->addBaseToKmer(ch) == false)
    goto anotherBase;

  pile = glob->getPile();

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


  while (glob->sorters[pile->pilePrefix].sorterLocked == 1)
    nanosleep(&naptime, 0L);

  int err = pthread_mutex_lock(&sorterLockMutex);
  if (err != 0)
    fprintf(stderr, "Failed to lock mutex (%d).  Fail.\n", err), exit(1);

  assert(glob->sorters[pile->pilePrefix].sorterLocked == 0);
  glob->sorters[pile->pilePrefix].sorterLocked = 1;

  err = pthread_mutex_unlock(&sorterLockMutex);
  if (err != 0)
    fprintf(stderr, "Failed to unlock mutex (%d).  Fail.\n", err), exit(1);

  glob->sorters[pile->pilePrefix].merge(pile);

  assert(glob->sorters[pile->pilePrefix].sorterLocked == 1);
  glob->sorters[pile->pilePrefix].sorterLocked = 0;

  //fprintf(stderr, "sorterThread()-- pile 0x"u32bitHEX" with "u32bitFMT" items -- sorter has "u32bitFMT" things now.\n",
  //        pile->pilePrefix, pile->pileLen, glob->sorters[pile->pilePrefix].sorterLen);
}



//  Does nothing but delete the pile object.  We don't output till the end.
void
nullThread(void *global, void *thing) {
  //kmerGlobal *glob = (kmerGlobal *)global;
  kmerPile   *pile = (kmerPile   *)thing;

  delete pile;
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

  ss->setNumberOfWorkers(4);
  ss->setWriterQueueSize(16384);

  //for (u32bit i=0; i<config._numSearchThreads; i++)
  //  ss->setThreadData(i, new searcherState(i));

  ss->run(kg, true);

  delete ss;

  kg->write();

  delete kg;

  exit(0);
}
