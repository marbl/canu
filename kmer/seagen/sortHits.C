#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  $Id$

//  Define this to get some extra console output about things that might
//  be of concern.
//
//#define DEVELOPMENT_CODE


typedef struct {
  u32bit dir;
  u32bit estID;
  u32bit scfID;
  u32bit scfLo;
  u32bit scfHi;
} hit_s;

char *tmpPath;


int
hitcmp(const void *a, const void *b) {
  hit_s  *A = (hit_s *)a;
  hit_s  *B = (hit_s *)b;

  if (A->scfID < B->scfID) return(-1);
  if (A->scfID > B->scfID) return(1);
  if (A->estID < B->estID) return(-1);
  if (A->estID > B->estID) return(1);
  if (A->scfLo < B->scfLo) return(-1);
  if (A->scfLo > B->scfLo) return(1);
  if (A->dir   < B->dir)   return(-1);
  if (A->dir   > B->dir)   return(1);
#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "WARN: Sort not stable!\n");
#endif
  return(0);
}


void
sortHits(hit_s *h, u32bit hp) {
#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "Sorting %u hits.\n", hp);
#endif
  qsort(h, hp, sizeof(hit_s), hitcmp);
}


void
writeHits(hit_s *h, u32bit& hp, u32bit& nt) {
  char  tmpname[1024];
  FILE *tmp;

#ifdef TRUE64BIT
  sprintf(tmpname, "%s.%03u", tmpPath, nt);
#else
  sprintf(tmpname, "%s.%03lu", tmpPath, nt);
#endif

  errno = 0;
  tmp = fopen(tmpname, "wb");
  if (tmp == 0L) {
    fprintf(stderr, "ERROR opening '%s' for writing.\n%s\n", tmpname, strerror(errno));
    exit(1);
  }

#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "Writing %u hits to '%s'\n", hp, tmpname);
#endif
  fwrite(h, sizeof(hit_s), hp, tmp);
  fclose(tmp);

  nt++;
  hp = 0;
}


void
writeHit(hit_s *h) {
  fprintf(stdout,
#ifdef TRUE64BIT
          "-%c -e %u -D %u %u %u\n",
#else
          "-%c -e %lu -D %lu %lu %lu\n",
#endif
          h->dir ? 'f' : 'r', h->estID, h->scfID, h->scfLo, h->scfHi);
}


int
main(int argc, char **argv) {

  if (argc < 4) {
    fprintf(stderr, "usage: %s [-v] [-m memorylimit] [-t temppath] hitfile1 hitfile2 ... > sorted-hits\n", argv[0]);
    fprintf(stderr, "       memory limit is MB\n");
    exit(1);
  }

  bool    beVerbose   = false;
  u64bit  memoryLimit = 128 * 1024 * 1024;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-m", 2) == 0) {
      arg++;
      memoryLimit   = atoi(argv[arg]);
      memoryLimit <<= 20;
    } else if (strncmp(argv[arg], "-t", 2) == 0) {
      arg++;
      tmpPath = argv[arg];
    } else {
      //  Must be at the first file name.  Break.
      break;
    }
    arg++;
  }


  //  Allocate a bunch of spaces to store hits.  This should be a
  //  command line option; we currently just allocate 2GB.
  //
  u64bit   hitsMax = memoryLimit / 20;
  u32bit   hitsPos = 0;
  hit_s   *hits    = new hit_s [hitsMax];

#ifdef DEVELOPMENT_CODE
  fprintf(stderr, "using memory limit of %u bytes --> %u hits.\n", memoryLimit, hitsMax);
  fprintf(stderr, "using tmp path of '%s'\n", tmpPath);
#endif

  //  Things for reading hits
  //
  FILE    *file;
  char     b[1024];
  aHit     a;
  bool     isBINARY;

  //  Things for merging
  //
  u32bit   numTemporary = 0;

  while (arg < argc) {

    //  Open the file, fatally failing if we cannot do it.
    //
    errno = 0;
    file = fopen(argv[arg], "r");
    if (file == 0L) {
      fprintf(stderr, "ESTmapper/filterEST-- ERROR opening '%s'\n%s\n", argv[arg], strerror(errno));
      exit(1);
    }

    //  Binary or ASCII input?
    //
    char x = (char)fgetc(file);
    ungetc(x, file);

    isBINARY = (x != '-');

#ifdef DEVELOPMENT_CODE
    if (isBINARY)
      fprintf(stderr, "reading BINARY hits from '%s'\n", argv[arg]);
    else
      fprintf(stderr, "reading ASCII hits from '%s'\n", argv[arg]);
#endif

    //  Read hits until we run out of space
    //
    while (!feof(file)) {
      if (isBINARY) {
        ahit_readBinary(&a, file);
      } else {
        fgets(b, 1024, file);
        ahit_parseString(&a, b);
      }

      if (beVerbose && ((hitsPos & 0xffff) == 0x0)) {
#ifdef TRUE64BIT
        fprintf(stderr, "  %u\r", hitsPos);
#else
        fprintf(stderr, "  %lu\r", hitsPos);
#endif
        fflush(stderr);
      }

      if (!feof(file)) {
        hits[hitsPos].dir   = a._direction;
        hits[hitsPos].estID = a._qsIdx;
        hits[hitsPos].scfID = a._dsIdx;
        hits[hitsPos].scfLo = a._dsLo;
        hits[hitsPos].scfHi = a._dsHi;

        hitsPos++;

        //  If we've run out of space, sort and write to a temporary file.
        //
        if (hitsPos == hitsMax) {
          sortHits(hits, hitsPos);
          writeHits(hits, hitsPos, numTemporary);
        }
      }
    }

    fclose(file);

    arg++;
  }


  //  All done reading.  If we have stuff to sort, sort it.
  //
  if (hitsPos > 0) {
    sortHits(hits, hitsPos);

    //  Just to make life easier, we write out the final piece.  This
    //  could be done with it in core, but that's harder.
    //
    if (numTemporary > 0)
      writeHits(hits, hitsPos, numTemporary);
  }


  //  If we've written temporary files, do a merge over all the files,
  //  and the stuff we have in core.  Otherwise, just write out the
  //  stuff we have in core.
  //
  if (numTemporary > 0) {
#ifdef DEVELOPMENT_CODE
    fprintf(stderr, "Merging.\n");
#endif

    FILE    **InF = new FILE * [numTemporary];
    hit_s    *InH = new hit_s [numTemporary];

    //  Open the input files
    //
    for (u32bit i=0; i<numTemporary; i++) {
      char tmpname[1024];
#ifdef TRUE64BIT
      sprintf(tmpname, "%s.%03u", tmpPath, i);
#else
      sprintf(tmpname, "%s.%03lu", tmpPath, i);
#endif
#ifdef DEVELOPMENT_CODE
      fprintf(stderr, "Opening '%s' for merging.\n", tmpname);
#endif
      errno = 0;
      InF[i] = fopen(tmpname, "rb");
      if (InF[i] == 0L) {
        fprintf(stderr, "ERROR opening '%s' for reading.\n%s\n", tmpname, strerror(errno));
        exit(1);
      }

      //  Read the first hit from each file.
      //
      if (fread(&InH[i], sizeof(hit_s), 1, InF[i]) == 0) {
#ifdef TRUE64BIT
        fprintf(stderr, "ESTmapper/sortHits-- WARNING: Temporary file %u has no hits?!\n", i);
#else
        fprintf(stderr, "ESTmapper/sortHits-- WARNING: Temporary file %lu has no hits?!\n", i);
#endif
        InH[i].dir   = ~u32bitZERO;
        InH[i].estID = ~u32bitZERO;
        InH[i].scfID = ~u32bitZERO;
        InH[i].scfLo = ~u32bitZERO;
        InH[i].scfHi = ~u32bitZERO;
      }
    }

    //  While there is still input, merge to the output
    //
    bool moreInput = true;
    hitsPos = 0;

    while (moreInput) {
      hitsPos++;
      if (beVerbose && ((hitsPos & 0xffff) == 0x0)) {
#ifdef TRUE64BIT
        fprintf(stderr, "Merging %u\r", hitsPos);
#else
        fprintf(stderr, "Merging %lu\r", hitsPos);
#endif
        fflush(stderr);
      }

      //  Pick the smallest hit -- if file [i] is finished, then hit[i]
      //  is bogus and all the values are set to maximal values.
      //
      u32bit smallestHit = 0;

      for (u32bit nextHit = smallestHit+1; nextHit < numTemporary; nextHit++) {
        if (hitcmp(InH+smallestHit, InH+nextHit) > 0)
          smallestHit = nextHit;
      }

      //  If the smallest hit is invalid, we're done.
      //
      if (InH[smallestHit].dir == ~u32bitZERO) {
        moreInput = false;
      } else {

        //  Otherwise, write the hit, and read a new one, invalidating
        //  it if the file is empty.
        //
        writeHit(InH + smallestHit);

        if (fread(&InH[smallestHit], sizeof(hit_s), 1, InF[smallestHit]) == 0) {
          InH[smallestHit].dir   = ~u32bitZERO;
          InH[smallestHit].estID = ~u32bitZERO;
          InH[smallestHit].scfID = ~u32bitZERO;
          InH[smallestHit].scfLo = ~u32bitZERO;
          InH[smallestHit].scfHi = ~u32bitZERO;
        }
      }
    }




    //  Should clean up, I know.


    //  Unlink the temporary files
    //
    for (u32bit i=0; i<numTemporary; i++) {
      char tmpname[1024];
#ifdef TRUE64BIT
      sprintf(tmpname, "%s.%03u", tmpPath, i);
#else
      sprintf(tmpname, "%s.%03lu", tmpPath, i);
#endif

      errno = 0;
      if (unlink(tmpname) == -1)
        fprintf(stderr, "ESTmapper/sortHits-- Can't remove '%s': %s\n", tmpname, strerror(errno));
    }

  } else {
#ifdef DEVELOPMENT_CODE
    fprintf(stderr, "No merge necessary.\n");
#endif
    for (u32bit i=0; i<hitsPos; i++)
      writeHit(hits+i);
  }

  delete [] hits;

  return(0);
}
