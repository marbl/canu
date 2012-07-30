#include "AS_BAT_Datatypes.H"
#include "AS_BAT_OverlapCache.H"

#include "memoryMappedFile.H"

const char *mainid = "$Id: rewriteCache.C,v 1.2 2012-07-30 01:21:01 brianwalenz Exp $";

uint64  ovlCacheMagic = 0x65686361436c766fLLU;  //0102030405060708LLU;

//  BOGART GLOBALS

FragmentInfo     *FI  = 0L;
OverlapCache     *OC  = 0L;
BestOverlapGraph *OG  = 0L;
ChunkGraph       *CG  = 0L;
InsertSizes      *IS  = 0L;

//  BOGART GLOBALS

int
main(int argc, char **argv) {

  char     name[FILENAME_MAX];
  FILE    *file;


  gkStore        *gkpStore  = new gkStore("salmon.gkpStore", FALSE, FALSE);
  FragmentInfo   *FI        = new FragmentInfo(gkpStore, "salmon");


  sprintf(name, "salmon.ovlCache");
  if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
    return(false);

  fprintf(stderr, "OverlapCache()-- Loading graph from '%s'.\n", name);

  errno = 0;

  file = fopen(name, "r");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   baterrbits = AS_BAT_ERRBITS;
  uint32   ovserrbits = AS_OVS_ERRBITS;
  uint32   ovshngbits = AS_OVS_HNGBITS;

  AS_UTL_safeRead(file, &magic,      "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeRead(file, &baterrbits, "overlapCache_baterrbits", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &ovserrbits, "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &ovshngbits, "overlapCache_ovshngbits", sizeof(uint32), 1);

  if (magic != ovlCacheMagic)
    fprintf(stderr, "OverlapCache()-- ERROR:  File '%s' isn't a bogart ovlCache.\n", name), exit(1);

  uint64  _memLimit;
  uint64  _memUsed;
  uint32  _maxPer;
  uint32  _batMax;

  AS_UTL_safeRead(file, &_memLimit, "overlapCache_memLimit", sizeof(uint64), 1);
  AS_UTL_safeRead(file, &_memUsed, "overlapCache_memUsed", sizeof(uint64), 1);

  AS_UTL_safeRead(file, &_maxPer, "overlapCache_maxPer", sizeof(uint32), 1);
  AS_UTL_safeRead(file, &_batMax, "overlapCache_batMax", sizeof(uint32), 1);

  //_bat = new BAToverlap [_batMax];

  uint32 *_OVSerate = new uint32 [1 << AS_OVS_ERRBITS];
  double *_BATerate = new double [1 << AS_BAT_ERRBITS];

  AS_UTL_safeRead(file,  _OVSerate, "overlapCache_OVSerate", sizeof(uint32), 1 << AS_OVS_ERRBITS);
  AS_UTL_safeRead(file,  _BATerate, "overlapCache_BATerate", sizeof(double), 1 << AS_BAT_ERRBITS);

  uint32 *_cacheLen = new uint32          [FI->numFragments() + 1];

  size_t numRead = AS_UTL_safeRead(file,  _cacheLen, "overlapCache_cacheLen", sizeof(uint32), FI->numFragments() + 1);

  if (numRead != FI->numFragments() + 1)
    fprintf(stderr, "OverlapCache()-- Short read loading graph '%s'.  Fail.\n", name), exit(1);

  fclose(file);

  //  Memory map the overlaps

  sprintf(name, "salmon.ovlCacheDat");
  FILE *infile = fopen(name, "r");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  sprintf(name, "salmon.ovlCacheDat.reduced");
  FILE *otfile = fopen(name, "w");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  uint32          batsMax = 16 * 1024 * 1024;
  BAToverlapInt  *bats    = new BAToverlapInt [batsMax];

  //  Update pointers into the overlaps

  uint64   nDel = 0;
  uint64   nMod = 0;
  uint64   nOvl = 0;

  writeLog("OverlapCache()-- Freshly deleted fragments detected.  Cleaning overlaps.\n");

  char  N[FILENAME_MAX];

  sprintf(N, "salmon.overlapsRemoved.log");

  errno = 0;
  FILE *F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "OverlapCache()--  Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

  for (uint32 fi=1; fi<FI->numFragments() + 1; fi++) {
    if (_cacheLen[fi] == 0)
      continue;

    AS_UTL_safeRead(infile, bats, "bats", sizeof(BAToverlapInt), _cacheLen[fi]);

    if ((FI->fragmentLength(fi) == 0) &&
        (_cacheLen[fi] > 0)) {
      nDel++;
      fprintf(F, "Removing "F_U32" overlaps from deleted deleted fragment "F_U32"\n", _cacheLen[fi], fi);
      _cacheLen[fi] = 0;
      continue;
    }

    uint32  on = 0;

    for (uint32 oi=0; oi<_cacheLen[fi]; oi++) {
      uint32  iid = bats[oi].b_iid;
      bool    del = (FI->fragmentLength(iid) == 0);

      if ((del == false) &&
          (on < oi))
        bats[on] = bats[oi];

      if (del == false)
        on++;
    }

    if (_cacheLen[fi] != on) {
      nMod++;
      nOvl += _cacheLen[fi] - on;
      fprintf(F, "Removing "F_U32" overlaps from living fragment "F_U32"\n", _cacheLen[fi] - on, fi);
    }

    _cacheLen[fi] = on;

    AS_UTL_safeWrite(otfile, bats, "bats", sizeof(BAToverlapInt), _cacheLen[fi]);
  }

  fclose(F);

  fclose(infile);
  fclose(otfile);

  fprintf(stderr, "OverlapCache()-- Removed all overlaps from "F_U64" deleted fragments.  Removed "F_U64" overlaps from "F_U64" alive fragments.\n",
          nDel, nOvl, nMod);



  //  DUMP new data file

  {
  char  name[FILENAME_MAX];
  FILE *file;

  sprintf(name, "salmon.ovlCache.reduced");

  fprintf(stderr, "OverlapCache()-- Saving graph to '%s'.\n", name);

  errno = 0;

  file = fopen(name, "w");
  if (errno)
    fprintf(stderr, "OverlapCache()-- Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

  uint64   magic      = ovlCacheMagic;
  uint32   baterrbits = AS_BAT_ERRBITS;
  uint32   ovserrbits = AS_OVS_ERRBITS;
  uint32   ovshngbits = AS_OVS_HNGBITS;

  AS_UTL_safeWrite(file, &magic,      "overlapCache_magic",      sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &baterrbits, "overlapCache_baterrbits", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &ovserrbits, "overlapCache_ovserrbits", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &ovshngbits, "overlapCache_ovshngbits", sizeof(uint32), 1);

  AS_UTL_safeWrite(file, &_memLimit, "overlapCache_memLimit", sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &_memUsed, "overlapCache_memUsed", sizeof(uint64), 1);

  AS_UTL_safeWrite(file, &_maxPer, "overlapCache_maxPer", sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &_batMax, "overlapCache_batMax", sizeof(uint32), 1);

  AS_UTL_safeWrite(file,  _OVSerate, "overlapCache_OVSerate", sizeof(uint32), 1 << AS_OVS_ERRBITS);
  AS_UTL_safeWrite(file,  _BATerate, "overlapCache_BATerate", sizeof(double), 1 << AS_BAT_ERRBITS);

  AS_UTL_safeWrite(file,  _cacheLen, "overlapCache_cacheLen", sizeof(uint32), FI->numFragments() + 1);

  fclose(file);
  }




  exit(0);
}
