#include "posix.H"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "aHit.H"

//  $Id$
//  $Log$
//  Revision 1.2  2003/01/03 15:57:13  walenz
//  added cvs stuff
//

#define MAX_ESTS    (16 * 1024 * 1024)

typedef struct {
  bool     stillMore;
  bool     reload;
  FILE    *file;
  char     b[1024];
  aHit     a;
  bool     isBINARY;
} hitFile_s;

typedef struct {
  aHit     a;
  u32bit   estid;

  float    coverage;
  float    multiplicity;
} hit_s;

void
loadHit(hitFile_s *HF) {

  if (HF->stillMore && HF->reload) {
    HF->reload = false;

    if (HF->isBINARY) {
      ahit_readBinary(&HF->a, HF->file);
    } else {
      fgets(HF->b, 1024, HF->file);
      ahit_parseString(&HF->a, HF->b);
    }

    if (feof(HF->file))
      HF->stillMore = false;
  }
}


FILE *openOutputFile(char *n) {
  FILE *f = 0L;
  if (n) {
    errno = 0;
    while ((f = fopen(n, "w")) == 0L) {
      fprintf(stderr, "filterEST: ERROR: couldn't open '%s' for writing.\n%s\n", n, strerror(errno));
      sleep(1);
    }
  }
  return(f);
}


int
main(int argc, char **argv) {

  if (argc < 3) {
    fprintf(stderr, "ESTmapper utility function -- not for human use.\n");
    exit(1);
  }

  u32bit       uniqThresh    = 100;
  u32bit       reptThresh    = 200;
  float        qualityThresh = 0.2;

  u32bit     *hitCounts      = new u32bit [MAX_ESTS];
  u32bit       numESTs       = 0;

  FILE        *uniqueFile    = 0L;
  FILE        *filteredFile  = 0L;
  FILE        *filteredJunk  = 0L;
  FILE        *repeatsFile   = 0L;
  FILE        *repeatsJunk   = 0L;

  hitFile_s   *HF            = new hitFile_s [argc];
  u32bit       numFiles      = 0;

  //  count the numbers of stuff:
  //  ESTs in each category -- unique, filtered and repeat
  //  hits in each category -- unique, filteredSaved, filteredJunked, repeatSaved and repeatJunked
  //
  u32bit       estsUniq     = 0;
  u32bit       estsFilt     = 0;
  u32bit       estsRept     = 0;

  u32bit       hitsUniq     = 0;
  u32bit       hitsFiltSave = 0;
  u32bit       hitsFiltJunk = 0;
  u32bit       hitsReptSave = 0;
  u32bit       hitsReptJunk = 0;

  bool         beVerbose = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-uniquethreshold", 2) == 0) {
      uniqThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-repeatthreshold", 2) == 0) {
      reptThresh  = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-qualitythreshold", 2) == 0) {
      qualityThresh = atof(argv[++arg]);
    } else if (strncmp(argv[arg], "-counts", 2) == 0) {
      errno = 0;
      FILE *F = fopen(argv[++arg], "r");
      if (F == 0L) {
        fprintf(stderr, "filterEST: ERROR: couldn't open '%s' for reading.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }
      while (!feof(F)) {
#ifdef TRUE64BIT
        fscanf(F, " %u", hitCounts + numESTs);
#else
        fscanf(F, " %lu", hitCounts + numESTs);
#endif
        if (!feof(F))
          numESTs++;
      }
      fclose(F);
    } else if (strncmp(argv[arg], "-Funique", 3) == 0) {
      uniqueFile   = openOutputFile(argv[++arg]);
    } else if (strncmp(argv[arg], "-Ffiltered", 3) == 0) {
      filteredFile   = openOutputFile(argv[++arg]);
    } else if (strncmp(argv[arg], "-Jfiltered", 3) == 0) {
      filteredJunk   = openOutputFile(argv[++arg]);
    } else if (strncmp(argv[arg], "-Frepeats", 3) == 0) {
      repeatsFile   = openOutputFile(argv[++arg]);
    } else if (strncmp(argv[arg], "-Jrepeats", 3) == 0) {
      repeatsJunk   = openOutputFile(argv[++arg]);
    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;
    } else {
      errno = 0;

      HF[numFiles].stillMore = true;
      HF[numFiles].reload    = true;
      HF[numFiles].file      = fopen(argv[arg], "r");

      if (HF[numFiles].file == 0L) {
        fprintf(stderr, "filterEST: ERROR: couldn't open '%s' for reading.\n%s\n", argv[arg], strerror(errno));
        exit(1);
      }

      //  Binary or ASCII input?
      //
      char x = (char)fgetc(HF[numFiles].file);
      ungetc(x, HF[numFiles].file);

      HF[numFiles].isBINARY = (x != '-');

      if (HF[numFiles].file == 0L) {
        fprintf(stderr, "filterEST: ERROR opening '%s'\n", argv[arg]);
      } else {
        numFiles++;
      }
    }

    arg++;
  }

  //
  //  Check the args
  //
  if (numESTs == 0) {
    fprintf(stderr, "filterEST: ERROR:  Didn't get any hitCounts!\n");
    exit(1);
  }

  if (beVerbose) {
#ifdef TRUE64BIT
    fprintf(stderr, "ESTmapper/filterEST-- uniqThresh=%4u  reptThresh=%4u  qualityThresh=%4.2f\n",
            uniqThresh, reptThresh, qualityThresh);
#else
    fprintf(stderr, "ESTmapper/filterEST-- uniqThresh=%4lu  reptThresh=%4lu  qualityThresh=%4.2f\n",
            uniqThresh, reptThresh, qualityThresh);
#endif
  }


  //  While there is input, merge
  //
  bool  keepGoing = true;

  while (keepGoing) {
    keepGoing = false;

    //  Load hits, if needed.
    //
    for (u32bit i=0; i<numFiles; i++)
      loadHit(HF+i);

    //  See if we're done.
    //
    for (u32bit i=0; i<numFiles; i++)
      keepGoing |= HF[i].stillMore;

    if (keepGoing) {

      //  Find the lowest ESTid
      //
      u32bit estOfInterest = 1 << 30;
      for (u32bit i=0; i<numFiles; i++)
        if ((HF[i].stillMore) && (estOfInterest > HF[i].a._qsIdx))
          estOfInterest = HF[i].a._qsIdx;

      //  Create a list of hits.  We know the size it's supposed to be,
      //  so we can allocate exactly that space.
      //
      //  For each file, save all hits with the estOfInterest to a list.
      //
      u32bit       listLen = 0;
      u32bit       listMax = hitCounts[estOfInterest];
      hit_s       *list    = new hit_s [listMax];

      float        bestScore  = 0.0;
      float        worstScore = 1.0;

      for (u32bit i=0; i<numFiles; i++) {
        while ((HF[i].stillMore) && (HF[i].a._qsIdx == estOfInterest)) {
          if (listLen >= listMax) {
#ifdef TRUE64BIT
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  EST of interest = %u\n", estOfInterest);
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  Too many hits (%u instead of %u)\n", listLen, listMax);
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  hitCounts might be invalid!\n");
#else
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  EST of interest = %lu\n", estOfInterest);
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  Too many hits (%lu instead of %lu)\n", listLen, listMax);
            fprintf(stderr, "ESTmapper/filterEST-- ERROR:  hitCounts might be invalid!\n");
#endif
            exit(1);
          }

          memcpy(&list[listLen].a, &HF[i].a, sizeof(aHit));

          list[listLen].coverage     = (float)HF[i].a._covered / (float)HF[i].a._numMers;
          list[listLen].multiplicity = (float)HF[i].a._matched / (float)HF[i].a._covered;

          //  aHit->_covered is in bases, but aHit->_numMers is the
          //  number of mers.  Possible for coverage to be > 1.0.
          //
          if (list[listLen].coverage > 1.0)
            list[listLen].coverage = 1.0;

          if (list[listLen].coverage > bestScore)
            bestScore = list[listLen].coverage;

          if (list[listLen].coverage < worstScore)
            worstScore = list[listLen].coverage;

          listLen++;

          HF[i].reload = true;
          loadHit(HF+i);
        }
      }

      if (listLen != hitCounts[estOfInterest]) {
#ifdef TRUE64BIT
        fprintf(stderr, "ESTmapper/filterEST-- ERROR: listLen != hitCounts[%lu] (%lu != %lu)\n",
                estOfInterest, listLen, hitCounts[estOfInterest]);
#else
        fprintf(stderr, "ESTmapper/filterEST-- ERROR: listLen != hitCounts[%lu] (%lu != %lu)\n",
                estOfInterest, listLen, hitCounts[estOfInterest]);
#endif
      }

      float qualityThreshUse = qualityThresh;

      //  Pick a quality threshold based on the input scores.  This is
      //  only used if there are more than uniqThresh hits for this
      //  EST.
      //
      qualityThreshUse = 0.2 * (bestScore - worstScore) + worstScore;
      if (qualityThreshUse > qualityThresh)
        qualityThreshUse = qualityThresh;

      if (listLen <= uniqThresh) {
        estsUniq++;
        hitsUniq += listLen;

        //  If we've got fewer hits than the repeat threshold, save
        //  ALL the hits to the unique file.
        //
        for (u32bit i=0; i < listLen; i++) {
          ahit_printASCII(&list[i].a, uniqueFile);
        }
      } else {
        //  Otherwise, filter the hits, and save them to the filtered file (if
        //  the number of filtered hits is less than the repeat threshold),
        //  or save the EST idx to the repeats file.
        //
        u32bit  count = 0;
        for (u32bit i=0; i < listLen; i++)
          if (qualityThreshUse <= list[i].coverage)
            count++;

        if (count <= reptThresh) {
          estsFilt++;
          hitsFiltSave += count;
          hitsFiltJunk += listLen - count;

          for (u32bit i=0; i < listLen; i++) {
            if (qualityThreshUse <= list[i].coverage) {
              ahit_printASCII(&list[i].a, filteredFile);
            } else {
              if (filteredJunk)
                ahit_printASCII(&list[i].a, filteredJunk);
            }
          }
        } else {
          estsRept++;
          hitsReptSave += count;
          hitsReptJunk += listLen - count;

#ifdef TRUE64BIT
          fprintf(repeatsFile, "%u\n", estOfInterest);
#else
          fprintf(repeatsFile, "%lu\n", estOfInterest);
#endif
          fflush(repeatsFile);


          //  Save the best max(reptThresh, x * count)
          //    What's a good value for x?
          //
          if (repeatsJunk) {
#if 0
            u32bit save = (int)(0.2 * count);
            if (save < reptThresh)
              save = reptThresh;
#else
            u32bit save = reptThresh;
#endif

            for (u32bit i=0; i < save; i++)
              ahit_printASCII(&list[i].a, repeatsJunk);
            fflush(repeatsJunk);
          }
        }
      }

      if (beVerbose && ((estOfInterest % 1000) == 0)) {
#ifdef TRUE64BIT
        fprintf(stderr, "ESTmapper/filterEST-- UNIQ:%8u(%8u)  FILT:%8u(%8u/%8u)  REPT:%8u(%8u/%8u)\r",
                estsUniq, hitsUniq,
                estsFilt, hitsFiltSave, hitsFiltJunk,
                estsRept, hitsReptSave, hitsReptJunk);
#else
        fprintf(stderr, "ESTmapper/filterEST-- UNIQ:%8lu(%8lu)  FILT:%8lu(%8lu/%8lu)  REPT:%8lu(%8lu/%8lu)\r",
                estsUniq, hitsUniq,
                estsFilt, hitsFiltSave, hitsFiltJunk,
                estsRept, hitsReptSave, hitsReptJunk);
#endif
        fflush(stderr);
      }

      delete [] list;
    }
  }

  if (beVerbose) {
#ifdef TRUE64BIT
    fprintf(stderr, "ESTmapper/filterEST-- UNIQ:%8u(%8u)  FILT:%8u(%8u/%8u)  REPT:%8u(%8u/%8u)\n",
            estsUniq, hitsUniq,
            estsFilt, hitsFiltSave, hitsFiltJunk,
            estsRept, hitsReptSave, hitsReptJunk);
#else
    fprintf(stderr, "ESTmapper/filterEST-- UNIQ:%8lu(%8lu)  FILT:%8lu(%8lu/%8lu)  REPT:%8lu(%8lu/%8lu)\n",
            estsUniq, hitsUniq,
            estsFilt, hitsFiltSave, hitsFiltJunk,
            estsRept, hitsReptSave, hitsReptJunk);
#endif
    fflush(stderr);
  }

  fclose(uniqueFile);
  fclose(filteredFile);
  fclose(repeatsFile);

  return(0);
}
