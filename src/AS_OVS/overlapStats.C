
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id$";

//  install.packages(c("akima"))
//
//  XX = read.table("xxx.UniqueEnds.lengtherror")
//  XX = read.table("xxx.RepeatEnds.lengtherror")
//  library(akima)
//
//  x=XX[,1]
//  y=XX[,2]
//  z=XX[,3]
//
//  YY = interp(x, y, z, xo=seq(min(x), max(x), length = 100), yo=seq(min(y), max(y), length = 100))
//
//  contour(YY, nlevels=50)
//  filled.contour(YY, nlevels=50, color=rainbow)
//
//  EE  = read.table("xxx.UniqueEnds.enderror")
//  EL  = read.table("xxx.UniqueEnds.endlength")

//
//  plot "xxx.UniqueEnds.enderror" with lines, "xxx.RepeatEnds.enderror" with lines
//  plot "xxx.UniqueEnds.endlength" with lines, "xxx.RepeatEnds.endlength" with lines
//

    ////////////////////////////////////////
    //
    //  1) Repeat ends - for repeat-ends, examine error rates and
    //  lengths.  (e.g., error rate vs length)
    //

    ////////////////////////////////////////
    //
    //  2) Polymorphism - for non-repeat ends, examine error rates and
    //  lengths.
    //

    ////////////////////////////////////////
    //
    //  3) Short insert check - count number of times a read overlaps
    //  with its mate.  Is it an innie overlap?  Is it a repeat end?
    //

    ////////////////////////////////////////
    //
    //  4) Genome length
    //

    ////////////////////////////////////////
    //
    //  5) Library randomness
    //

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "AS_OVS_overlap.H"
#include "AS_OVS_overlapFile.H"
#include "AS_OVS_overlapStore.H"

#include "overlapStore.H"
#include "AS_UTL_histogram.H"


//  Any fragment end with more than this number of overlaps is labeled
//  a repeat.  One could probably find a better model using the length
//  of the overlap too.
//
typedef struct {
  AS_UTL_histogram    hist5;
  AS_UTL_histogram    hist3;
  uint64              repeatThreshold;
} RepeatModel;


//  [0] -- global
//  [1] -- 5'
//  [2] -- 3'
typedef struct {
  AS_UTL_histogram   overlapLength[3];
  AS_UTL_histogram   errorRates[3];
  AS_UTL_histogram3d lengthError[3];
} FragmentEndData;


AS_IID   *fragLibrary     = NULL;
AS_IID   *fragMateIID     = NULL;
uint16   *fragClearLength = NULL;

char      outputPrefix[FILENAME_MAX] = {0};

typedef struct {
   uint64          * readsPerLibrary;
   uint64         ** libraryVsLibraryOverlaps;
   HashTable_AS    * readsSeen;
   uint64            contained;
   uint64            totalOverlaps;
} LibraryOverlapData;

#define MIN_LIBRARY_SIZE 1000


//
//  Start of BoringStuff
//


void
loadClearLengths(gkStore *gkp) {

  if (fragClearLength != NULL)
    return;

  gkStream        *frgStream = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment       fr;
  uint64           maxIID    = gkp->gkStore_getNumFragments() + 1;

  fragLibrary     = (AS_IID *)safe_calloc(maxIID, sizeof(AS_IID));
  fragMateIID     = (AS_IID *)safe_calloc(maxIID, sizeof(AS_IID));
  fragClearLength = (uint16 *)safe_calloc(maxIID, sizeof(uint16));

  int  typ = AS_READ_CLEAR_LATEST;

#if 0
  //  In general, we need to load a different clear length based on
  //  the type of overlap.  Since we currently only compute stats on
  //  OVL overlaps, this is disabled.

  switch (ovl->dat.ovl.type) {
    case AS_OVS_TYPE_OVL:
      typ = AS_READ_CLEAR_OBT;
      break;
    case AS_OVS_TYPE_OBT:
      typ = AS_READ_CLEAR_OBTINI;
      break;
    case AS_OVS_TYPE_MER:
      typ = AS_READ_CLEAR_UNTRIM;
      break;
    default:
      fprintf(stderr, "Unknown type %d in overlap.\n", ovl->dat.ovl.type);
      exit(1);
      break;
  }
#endif

  while (frgStream->next(&fr)) {
    AS_IID  iid = fr.gkFragment_getReadIID();
    uint32  b, e;

    fragLibrary[iid]     = fr.gkFragment_getLibraryIID();
    fragMateIID[iid]     = fr.gkFragment_getMateIID();

    fr.gkFragment_getClearRegion(b, e, typ);

    fragClearLength[iid] = e - b;
  }

  delete frgStream;
}




//  a_hang   b_hang     type    label
//  (compared to zero)          (describes A)
//
//  any      any        0x00    all overlaps
//                      %0000
//
//  =        =          0x05    degenerate     ----------
//                      %0101                  ----------
//
//  =        <          0x06    5' containee   ----------
//                      %0110                  ---
//
//  =        >          0x07    5' contained   ---
//                      %0111                  ----------
//
//  <        =          0x09    3' contained          ---
//                      %1001                  ----------
//
//  <        <          0x0a    5' dovetal        ----------
//                      %1010                  ----------
//
//  <        >          0x0b    contained         ----
//                      %1011                  ----------
//
//  >        =          0x0d    3' containee   ----------
//                      %1101                         ---
//
//  >        <          0x0e    containee      ----------
//                      %1110                     ----
//
//  >        >          0x0f    3' dovetail    ----------
//                      %1111                     ----------
//
//  not an overlap -- 0x00
//  degenerate     -- 0x05
//  5' types       -- 0x06, 0x07, 0x0a
//  3' types       -- 0x09, 0x0d, 0x0f
//  C  types       -- 0x0b
//  C  types       -- 0x0e
//  unused types   -- 0x01, 0x02, 0x03, 0x04, 0x08, 0x0c
//
//
//
uint32
computeTypeOfOverlap(OVSoverlap ovl) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  tp = 0;

  if (ah == 0)
    tp |= 0x00000004;
  else if (ah < 0)
    tp |= 0x00000008;
  else
    tp |= 0x0000000c;

  if (bh == 0)
    tp |= 0x00000001;
  else if (bh < 0)
    tp |= 0x00000002;
  else
    tp |= 0x00000003;

  return(tp);
}

uint32
overlapTypeIs3prime(uint32 tp) {
  return((tp == 0x06) || (tp == 0x07) || (tp == 0x0a));
}

uint32
overlapTypeIs5prime(uint32 tp) {
  return((tp == 0x09) || (tp == 0x0d) || (tp == 0x0f));
}


//  Swiped from AS_BOG/AS_BOG_BestOverlapGraph.cc::olapLength()
uint32
computeLengthOfOverlap(OVSoverlap ovl) {
  int32   ah = ovl.dat.ovl.a_hang;
  int32   bh = ovl.dat.ovl.b_hang;
  uint32  le = 0;

  if (ah < 0) {
    if (bh < 0)
      le = fragClearLength[ovl.a_iid] + bh;
    else
      le = fragClearLength[ovl.b_iid] + ah - bh;
  } else {
    if (bh < 0)
      le = fragClearLength[ovl.a_iid] + bh - ah;
    else
      le = fragClearLength[ovl.a_iid] - ah;
  }

  return(le);
}



RepeatModel *
computeRepeatModels(OverlapStore *ovs, gkStore *gkp) {
  int          i;
  uint32       ovl5   = 0;
  uint32       ovl3   = 0;
  AS_IID       lastID = 0;
  OVSoverlap   ovl;

  //  The [0] repeat model is a global model, all the others are
  //  specific to a single library.

  //  Allocate N models, one for each library.
  RepeatModel *rm = (RepeatModel *)safe_calloc(gkp->gkStore_getNumLibraries() + 1, sizeof(RepeatModel));

  //  Then populate the histograms.
  for (i=0; i <= gkp->gkStore_getNumLibraries(); i++) {
    AS_UTL_histogramAllocate(&rm[i].hist5);
    AS_UTL_histogramAllocate(&rm[i].hist3);
  }

  AS_OVS_resetRangeOverlapStore(ovs);
  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OVL) == TRUE) {

    if (ovl.a_iid != lastID) {

      //  Update the stats for this fragment.
      if (lastID > 0) {
        AS_IID  libIID = fragLibrary[lastID];

        AS_UTL_histogramAdd(&rm[0].hist5, ovl5);
        AS_UTL_histogramAdd(&rm[0].hist3, ovl3);

        AS_UTL_histogramAdd(&rm[libIID].hist5, ovl5);
        AS_UTL_histogramAdd(&rm[libIID].hist3, ovl3);
      }

      ovl5   = 0;
      ovl3   = 0;
      lastID = ovl.a_iid;
    }

    uint32 tp = computeTypeOfOverlap(ovl);

    if (overlapTypeIs5prime(tp))  ovl5++;
    if (overlapTypeIs3prime(tp))  ovl3++;
  }

  //  Update the stats for the last fragment.
  if (lastID > 0) {
    AS_IID  libIID = fragLibrary[lastID];

    AS_UTL_histogramAdd(&rm[0].hist5, ovl5);
    AS_UTL_histogramAdd(&rm[0].hist3, ovl3);

    AS_UTL_histogramAdd(&rm[libIID].hist5, ovl5);
    AS_UTL_histogramAdd(&rm[libIID].hist3, ovl3);
  }

  //  Examine the histograms, build a model.

  for (i=0; i <= gkp->gkStore_getNumLibraries(); i++) {
    char  label[256] = {0};
    char  name[FILENAME_MAX] = {0};
    FILE *file = NULL;

    AS_UTL_histogramCompute(&rm[i].hist5);
    AS_UTL_histogramCompute(&rm[i].hist3);

#warning bogus compute of repeatThreshold
    rm[i].repeatThreshold = ((rm[i].hist5.mode + 2 * rm[i].hist5.mad) +
                             (rm[i].hist3.mode + 2 * rm[i].hist3.mad)) / 2;

    sprintf(name, "%s.repeatmodel.lib.%03d.stats", outputPrefix, i);

    errno = 0;
    file = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "Couldn't open '%s' for write: %s\n", name, strerror(errno));
      exit(1);
    }

    fprintf(file, "repeatThreshold = "F_U64"\n", rm[i].repeatThreshold);

    sprintf(label, "Lib "F_IID" 5'", i);
    AS_UTL_histogramShow(&rm[i].hist5, file, label);

    sprintf(label, "Lib "F_IID" 3'", i);
    AS_UTL_histogramShow(&rm[i].hist3, file, label);

    fclose(file);

    sprintf(name,  "%s.repeatmodel.lib.%03d.5prime.dat", outputPrefix, i);
    sprintf(label, "%s.repeatmodel.lib.%03d.5prime.dat", outputPrefix, i);
    AS_UTL_histogramDump(&rm[i].hist5, name, label);

    sprintf(name,  "%s.repeatmodel.lib.%03d.3prime.dat", outputPrefix, i);
    sprintf(label, "%s.repeatmodel.lib.%03d.3prime.dat", outputPrefix, i);
    AS_UTL_histogramDump(&rm[i].hist3, name, label);
  }

#if 1
  //  Dump the repeat thresholds
  fprintf(stderr, "== Repeat Model ==\n");
  fprintf(stderr, "\n");
  for (i=0; i <= gkp->gkStore_getNumLibraries(); i++)
    fprintf(stderr, "repeatThreshold[%2d] = "F_U64"\n", i, rm[i].repeatThreshold);
  AS_UTL_histogramShow(&rm[0].hist5, stderr, "Global 5'");
  AS_UTL_histogramShow(&rm[0].hist3, stderr, "Global 3'");
#endif

  return(rm);
}


//  Returns 0 if neither end is a repeat;
//          1 if both ends are repeats;
//          5 if the 5' end is;
//          3 if the 3' end is.
int
isRepeatEnd(OVSoverlap *ovls, uint64 ovlsLen, RepeatModel *rm) {
  uint64  i, n5=0, n3=0;

  //  If no chance of being a repeat, get out of here.
  if (ovlsLen < rm[0].repeatThreshold)
    return(0);

  for (i=0; i<ovlsLen; i++) {
    uint32  tp = computeTypeOfOverlap(ovls[i]);

    if (overlapTypeIs5prime(tp))
      n5++;

    if (overlapTypeIs3prime(tp))
      n3++;
  }

  if ((n5 >= rm[0].repeatThreshold) &&
      (n3 >= rm[0].repeatThreshold))
    return(1);

  if (n5 >= rm[0].repeatThreshold)
    return(5);

  if (n3 >= rm[0].repeatThreshold)
    return(3);

  return(0);
}


//
//  End of BoringStuff
//




FragmentEndData *
process_FragmentEnds(OVSoverlap *ovls, uint64 ovlsLen, gkStore *gkp,
                     RepeatModel *rm,
                     FragmentEndData *red,
                     int isRepeat) {
  int i;
  int repeatend = isRepeatEnd(ovls, ovlsLen, rm);
  int stats5    = 0;
  int stats3    = 0;

  //  Sometimes, I wish I'd just duplicate the code and make the minor
  //  tweaks needed.

  //  Depending on the isRepeat flag, decide if we want to accumulate
  //  stats for either end.
  //
  if (isRepeat) {
    //  Want to examine the repeats.
    if ((repeatend == 1) || (repeatend == 5))  stats5 = 1;
    if ((repeatend == 1) || (repeatend == 3))  stats3 = 1;
  } else {
    //  Want to examine the uniques.
    if ((repeatend == 0) || (repeatend == 3))  stats5 = 1;
    if ((repeatend == 0) || (repeatend == 5))  stats3 = 1;
  }

  //  Get outta here if we're not examining either end.
  //
  if ((stats5 == 0) && (stats3 == 0))
    return(red);


  if (red == NULL) {
    red = (FragmentEndData *)safe_calloc(1, sizeof(FragmentEndData));

    for (i=0; i<3; i++) {
      AS_UTL_histogramAllocate(&red->overlapLength[i]);
      AS_UTL_histogramAllocate(&red->errorRates[i]);
      AS_UTL_histogram3dAllocate(&red->lengthError[i], AS_READ_MAX_NORMAL_LEN+1, 300);
    }
  }


  for (i=0; i<ovlsLen; i++) {
    uint32  tp = computeTypeOfOverlap(ovls[i]);
    uint32  le = computeLengthOfOverlap(ovls[i]);

    double  origerr = AS_OVS_decodeQuality(ovls[i].dat.ovl.orig_erate) * 100.0;
    double  correrr = AS_OVS_decodeQuality(ovls[i].dat.ovl.corr_erate) * 100.0;

    uint64  e = (uint64)floor(origerr * 10);

    //  Update the length x errorrate histogram.

    AS_UTL_histogramAdd(&red->overlapLength[0], le);
    AS_UTL_histogramAdd(&red->errorRates[0], e);
    AS_UTL_histogram3dAdd(&red->lengthError[0], le, e);

    if (stats5) {
      AS_UTL_histogramAdd(&red->overlapLength[1], le);
      AS_UTL_histogramAdd(&red->errorRates[1], e);
      AS_UTL_histogram3dAdd(&red->lengthError[1], le, tp);
    }

    if (stats3) {
      AS_UTL_histogramAdd(&red->overlapLength[2], le);
      AS_UTL_histogramAdd(&red->errorRates[2], e);
      AS_UTL_histogram3dAdd(&red->lengthError[2], le, e);
    }
  }

  return(red);
}

void
finalize_FragmentEnds(gkStore *gkp, RepeatModel *rm, FragmentEndData *red, char *label) {
  char  name[FILENAME_MAX];
  int   i;

  if (red == NULL)
    return;

  for (i=0; i<3; i++) {
    AS_UTL_histogramCompute(red->overlapLength + i);
    AS_UTL_histogramCompute(red->errorRates + i);
    AS_UTL_histogram3dCompute(red->lengthError + i);
  }


  sprintf(name, "%s.%s.endlength", outputPrefix, label);
  AS_UTL_histogramDump(&red->overlapLength[0], name, "both end length");

  sprintf(name, "%s.%s.enderror", outputPrefix, label);
  AS_UTL_histogramDump(&red->errorRates[0], name, "both end error (% x10)");


  sprintf(name, "%s.%s.endlength.5", outputPrefix, label);
  AS_UTL_histogramDump(&red->overlapLength[1], name, "5' end length");

  sprintf(name, "%s.%s.enderror.5", outputPrefix, label);
  AS_UTL_histogramDump(&red->errorRates[1], name, "5' end error (% x10)");


  sprintf(name, "%s.%s.endlength.3", outputPrefix, label);
  AS_UTL_histogramDump(&red->overlapLength[2], name, "3' end length");

  sprintf(name, "%s.%s.enderror.3", outputPrefix, label);
  AS_UTL_histogramDump(&red->errorRates[2], name, "3' end error (% x10)");


  sprintf(name, "%s.%s.lengtherror", outputPrefix, label);
  AS_UTL_histogram3dDump(&red->lengthError[0], name, "both end length X error");

  sprintf(name, "%s.%s.lengtherror.5", outputPrefix, label);
  AS_UTL_histogram3dDump(&red->lengthError[1], name, "5' end length X error");

  sprintf(name, "%s.%s.lengtherror.3", outputPrefix, label);
  AS_UTL_histogram3dDump(&red->lengthError[2], name, "3' end length X error");
}





void
process_ShortInsert(OVSoverlap *ovls, uint64 ovlsLen, gkStore *gkp, RepeatModel *rm) {
}

void
finalize_ShortInsert(gkStore *gkp, RepeatModel *rm) {
}



void
process_GenomeLength(OVSoverlap *ovls, uint64 ovlsLen, gkStore *gkp, RepeatModel *rm) {
}

void
finalize_GenomeLength(gkStore *gkp, RepeatModel *rm) {
}



LibraryOverlapData *
process_LibraryRandomness(OVSoverlap *ovls, uint64 ovlsLen, gkStore *gkp,
                          RepeatModel *rm,
                          LibraryOverlapData *ovl) {
   int i = 0, j = 0;
   int repeatend = isRepeatEnd(ovls, ovlsLen, rm);
   gkFragment fr;

   // initialize data if necessary
   // library 0 represents reads with no library association
   if (ovl == NULL) {
      ovl = (LibraryOverlapData *)safe_calloc(1, sizeof(LibraryOverlapData));

      int32 numLibraries = gkp->gkStore_getNumLibraries()+1;
      ovl->readsPerLibrary = (uint64 *)safe_malloc(numLibraries * sizeof(uint64));
      for (i = 0; i < numLibraries; i++) {
         ovl->readsPerLibrary[i] = 0;
      }

      ovl->libraryVsLibraryOverlaps = (uint64 **)safe_malloc(numLibraries * sizeof(uint64 *));
      for (i = 0; i < numLibraries; i++) {
        ovl->libraryVsLibraryOverlaps[i] = (uint64 *)safe_malloc(numLibraries * sizeof(uint64));
         for (j = 0; j < numLibraries; j++) {
            ovl->libraryVsLibraryOverlaps[i][j] = 0;
         }
      }

      ovl->readsSeen = CreateScalarHashTable_AS();
      ovl->contained = 0;
      ovl->totalOverlaps = 0;
   }

   if (ovlsLen > 0) {
      // all the records have the same a_iid so we only need to get it once
      gkp->gkStore_getFragment(ovls[0].a_iid, &fr, GKFRAGMENT_INF);
      AS_IID libOne = fr.gkFragment_getLibraryIID();

      ovl->totalOverlaps += ovlsLen;

      // get the libraries of the reads from the gkp store and increment appropriate counters
      for (i=0; i<ovlsLen; i++) {
         uint32 contained = computeTypeOfOverlap(ovls[i]);

         // get the contained overlap and don't count ones that are contained
         if (contained != 0x0a && contained != 0x0f) {
            ovl->contained++;
         }
         else {
            gkp->gkStore_getFragment(ovls[i].b_iid, &fr, GKFRAGMENT_INF);
            AS_IID libTwo = fr.gkFragment_getLibraryIID();

            if (!ExistsInHashTable_AS(ovl->readsSeen, ovls[i].a_iid, 0)) {
               ovl->readsPerLibrary[libOne]++;
               InsertInHashTable_AS(ovl->readsSeen, ovls[i].a_iid, 0, 1, 0);
            }

            if (!ExistsInHashTable_AS(ovl->readsSeen, ovls[i].b_iid, 0)) {
               ovl->readsPerLibrary[libTwo]++;
               InsertInHashTable_AS(ovl->readsSeen, ovls[i].b_iid, 0, 1, 0);
            }

            ovl->libraryVsLibraryOverlaps[libOne][libTwo]++;
         }
      }
   }

   return (ovl);
}

void
finalize_LibraryRandomness(gkStore *gkp, RepeatModel *rm, LibraryOverlapData *ovl) {
   // output the results of our library analysis
   char  name[FILENAME_MAX];
   int i = 0, j = 0;

   if (ovl == NULL)
      return;

   sprintf(name, "%s.libVsLib", outputPrefix);
   FILE   *F;
   errno = 0;
   F = fopen(name, "w");
   if (errno) {
      fprintf(stderr, "Failed to open '%s' for write: %s.\n", name, strerror(errno));
      return;
   }

   int32 numLibraries = gkp->gkStore_getNumLibraries()+1;
   uint64 uncontained = ovl->totalOverlaps - ovl->contained;
   uint64 numReads = ovl->readsSeen->numNodes;

   fprintf(F, "Total overlaps: %ld\tcontained: %ld\tcontained percent: %.2f\n", ovl->totalOverlaps, ovl->contained, 100*(ovl->contained/(double)ovl->totalOverlaps));
   fprintf(F, "Overlaps per lib in [libID] [%% reads] format\n");
   for (i = 0; i < numLibraries; i++) {
      fprintf(F, "%d\t%.2f\n", i, 100*(ovl->readsPerLibrary[i]/(double)numReads));
   }

   for (i = 1; i < numLibraries; i++) {
      for (j = 1; j < numLibraries; j++) {
         if (ovl->readsPerLibrary[i] < MIN_LIBRARY_SIZE || ovl->readsPerLibrary[j] < MIN_LIBRARY_SIZE) {
            fprintf(F, "%d\t%d\t%s\n", i, j, "NA");
         }
         else {
            uint64 overlaps = ovl->libraryVsLibraryOverlaps[i][j];

            double expected = (ovl->readsPerLibrary[i]/(double)numReads)*(ovl->readsPerLibrary[j]/(double)(numReads));
            double expectedCount = uncontained*expected;

            double actual = (overlaps/(double)uncontained);
            double diffRatio = (overlaps - expectedCount) / expectedCount;

            fprintf(F, "%d\t%d\t%.2f\t"F_U64"\t%.2f\n", i, j, expectedCount, overlaps, diffRatio*100);
         }
      }
   }
}







int
main(int argc, char **argv) {
  char   *ovsName = NULL;
  char   *gkpName = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-O") == 0) {
      ovsName = argv[++arg];
    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];
    } else if (strcmp(argv[arg], "-o") == 0) {
      strcpy(outputPrefix, argv[++arg]);
    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (ovsName == NULL) || (gkpName == NULL)) {
    fprintf(stderr, "usage: %s -G gkpStore -O ovlStore -o outputPrefix\n", argv[0]);
    exit(1);
  }

  //  Open the stores.

  OverlapStore    *ovs = AS_OVS_openOverlapStore(ovsName);
  gkStore         *gkp = new gkStore(gkpName, FALSE, FALSE);

  //  Make sure the store contains OVLs, and not OBTs or MERs.

  {
    OVSoverlap ovl;

    AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_ANY);
    if (ovl.dat.ovl.type != AS_OVS_TYPE_OVL) {
      fprintf(stderr, "stats only supported for OVL type overlaps.\n");
      exit(1);
    }
  }

  loadClearLengths(gkp);

  RepeatModel *rm = computeRepeatModels(ovs, gkp);

  //  The goal is to make one pass through the overlap store,
  //  computing and/or remembering whatever we need to know for the
  //  various statistics.

  OVSoverlap   ovl;
  uint64       ovlsLen = 0;
  uint64       ovlsMax = 1048576;
  OVSoverlap  *ovls    = (OVSoverlap *)safe_malloc(ovlsMax * sizeof(OVSoverlap));;
  AS_IID       ovlsIID = 0;

  FragmentEndData *repeat   = NULL;
  FragmentEndData *unique   = NULL;

  LibraryOverlapData *overlaps = NULL;

  AS_OVS_resetRangeOverlapStore(ovs);
  while (AS_OVS_readOverlapFromStore(ovs, &ovl, AS_OVS_TYPE_OVL) == TRUE) {

    //  Suck in all the overlaps for A_IID, then pass those to
    //  something that will collect stats.  Why?  The repeat ends
    //  cannot be computed until we have all overlaps for a single
    //  read.
    //
    //  It'd be nice to separate them into three classes: 5', 3' and
    //  contained.

    if (ovl.a_iid != ovlsIID) {
      repeat = process_FragmentEnds(ovls, ovlsLen, gkp, rm, repeat, 1);
      unique = process_FragmentEnds(ovls, ovlsLen, gkp, rm, unique, 0);

      process_ShortInsert      (ovls, ovlsLen, gkp, rm);
      process_GenomeLength     (ovls, ovlsLen, gkp, rm);
      overlaps = process_LibraryRandomness(ovls, ovlsLen, gkp, rm, overlaps);

      ovlsLen = 0;
      ovlsIID = ovl.a_iid;
    }

    if (ovlsLen >= ovlsMax) {
      ovlsMax *= 2;
      ovls     = (OVSoverlap *)safe_realloc(ovls, ovlsMax * sizeof(OVSoverlap));
    }

    ovls[ovlsLen++] = ovl;
  }

  if (ovlsLen > 0) {
    repeat = process_FragmentEnds(ovls, ovlsLen, gkp, rm, repeat, 1);
    unique = process_FragmentEnds(ovls, ovlsLen, gkp, rm, unique, 0);

    process_ShortInsert      (ovls, ovlsLen, gkp, rm);
    process_GenomeLength     (ovls, ovlsLen, gkp, rm);
    overlaps = process_LibraryRandomness(ovls, ovlsLen, gkp, rm, overlaps);
  }

  //  After our one pass, we can finish up the statistics.

  finalize_FragmentEnds(gkp, rm, repeat, "RepeatEnds");
  finalize_FragmentEnds(gkp, rm, unique, "UniqueEnds");

  finalize_ShortInsert      (gkp, rm);
  finalize_GenomeLength     (gkp, rm);
  finalize_LibraryRandomness(gkp, rm, overlaps);

  //  Clean up.

  return(0);
}
