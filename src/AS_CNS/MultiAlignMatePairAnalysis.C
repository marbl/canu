
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: MultiAlignMatePairAnalysis.C,v 1.1 2011-12-16 03:21:45 brianwalenz Exp $";

#include "MultiAlignMatePairAnalysis.H"

#include <math.h>

#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;


#define MPA_ORIENT_INNIE      0x00
#define MPA_ORIENT_OUTTIE     0x01
#define MPA_ORIENT_NORMAL     0x02
#define MPA_ORIENT_ANTINORMAL 0x03

const char MPA_ORIENT_NAMES[4] = { 'I', 'O', 'N', 'A' };

class mpaLibrary {
public:
  mpaLibrary(gkLibrary *lib) {
    library    = lib;

    strcpy(libraryName, lib->libraryName);

    origMean   = lib->mean;
    origStdDev = lib->stddev;

    newSamples = 0;
    newMean    = 0;
    newStdDev  = 0;

    for (uint32 i=0; i<4; i++) {
      externHappy[i] = 0;
      externSad[i]   = 0;
    }
  };

  ~mpaLibrary() {
  };

  gkLibrary      *library;
  char            libraryName[LIBRARY_NAME_SIZE];

  double          origMean;
  double          origStdDev;

  uint64          newSamples;
  double          newMean;
  double          newStdDev;

  //  For externally mated fragments, happy?
  uint64          externHappy[4];
  uint64          externSad[4];

  //  Distances, one per pair orientation
  vector<int32>   dist[4];
};




matePairAnalysis::matePairAnalysis(char *gkpName) {

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  libdata  = new mpaLibrary * [gkpStore->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkpStore->gkStore_getNumLibraries() + 1; i++)
    libdata[i] = new mpaLibrary(gkpStore->gkStore_getLibrary(i));

  pairing = new AS_IID [gkpStore->gkStore_getNumFragments() + 1];
  library = new AS_IID [gkpStore->gkStore_getNumFragments() + 1];

  gkStream   *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment  fr;

  while(fs->next(&fr)) {
    uint32 iid = fr.gkFragment_getReadIID();

    pairing[iid] = fr.gkFragment_getMateIID();
    library[iid] = fr.gkFragment_getLibraryIID();

    if ((iid % 10000000) == 0)
      fprintf(stderr, "Loading fragment information %9d out of %9d\n", iid, gkpStore->gkStore_getNumFragments());
  }
  
  delete fs;
}



matePairAnalysis::~matePairAnalysis() {

  for (uint32 i=0; i<gkpStore->gkStore_getNumLibraries() + 1; i++)
    delete libdata[i];
  delete [] libdata;

  delete [] pairing;
  delete [] library;

  delete gkpStore;
}





void
matePairAnalysis::evaluateTig(MultiAlignT *ma) {
  map<AS_IID, uint32>  fragIdx;
  uint32               numFrag = GetNumIntMultiPoss(ma->f_list);
  uint32               maLen   = 0;

  for (uint32 i=0; i<numFrag; i++) {
    IntMultiPos *frg = GetIntMultiPos(ma->f_list, i);

    fragIdx[frg->ident] = i;

    maLen = MAX(maLen, frg->position.bgn);
    maLen = MAX(maLen, frg->position.end);
  }


  for (uint32 i=0; i<numFrag; i++) {
    IntMultiPos *frg = GetIntMultiPos(ma->f_list, i);
    AS_IID       mmi = pairing[frg->ident];

    if (mmi == 0)
      //  Not mated
      continue;

    AS_IID      libiid    = library[frg->ident];
    mpaLibrary *lib       = libdata[libiid];
    uint32      liborient = lib->library->orientation;
    int32       libmin    = lib->origMean - 3 * lib->origStdDev;
    int32       libmax    = lib->origMean + 3 * lib->origStdDev;

    bool        frgforward = (frg->position.bgn < frg->position.end);
    int32       frgmin     = MIN(frg->position.bgn, frg->position.end);
    int32       frgmax     = MAX(frg->position.bgn, frg->position.end);

    //  EXTERNALLY MATED
    if (fragIdx.find(mmi) == fragIdx.end()) {
      switch (liborient) {
        case AS_READ_ORIENT_UNKNOWN:
          assert(0);
          break;

        case AS_READ_ORIENT_INNIE:
          if (frgforward == true)   ((frgmin < maLen - libmax) ? lib->externSad[MPA_ORIENT_INNIE] : lib->externHappy[MPA_ORIENT_INNIE])++;
          if (frgforward == false)  ((frgmin >         libmax) ? lib->externSad[MPA_ORIENT_INNIE] : lib->externHappy[MPA_ORIENT_INNIE])++;
          break;

        case AS_READ_ORIENT_OUTTIE:
          if (frgforward == false)   ((frgmin < maLen - libmax) ? lib->externSad[MPA_ORIENT_OUTTIE] : lib->externHappy[MPA_ORIENT_OUTTIE])++;
          if (frgforward == true)    ((frgmin >         libmax) ? lib->externSad[MPA_ORIENT_OUTTIE] : lib->externHappy[MPA_ORIENT_OUTTIE])++;
          break;

        case AS_READ_ORIENT_NORMAL:
          assert(0);
          break;

        case AS_READ_ORIENT_ANTINORMAL:
          assert(0);
          break;

        default:
          break;
      }

      continue;
    }


    IntMultiPos  *mat        = GetIntMultiPos(ma->f_list, fragIdx[mmi]);
    bool          matforward = (mat->position.bgn < mat->position.end);
    int32         matmin     = MIN(mat->position.bgn, mat->position.end);
    int32         matmax     = MAX(mat->position.bgn, mat->position.end);

    if (frgmin > matmin)
      //  Skip it, wait for frg to be before mat
      continue;

    int32 distance = matmax - frgmin;
    
    if        ((frgforward == true)  && (matforward == false)) {
      lib->dist[MPA_ORIENT_INNIE].push_back(distance);

    } else if ((frgforward == false) && (matforward == true)) {
      lib->dist[MPA_ORIENT_OUTTIE].push_back(distance);

    } else if ((frgforward == true)  && (matforward == true)) {
      lib->dist[MPA_ORIENT_NORMAL].push_back(distance);

    } else {//((frgforward == false) && (matforward == false))
      lib->dist[MPA_ORIENT_ANTINORMAL].push_back(distance);
    }
  }
}



void
matePairAnalysis::summarize(char *prefix) {
  char   datName[FILENAME_MAX];

  if (prefix[strlen(prefix)-1] == '/')
    prefix[strlen(prefix)-1] = 0;

  sprintf(datName, "%s.gp", prefix);
  errno = 0;
  FILE *gp = fopen(datName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  Failed to open '%s': %s\n", datName, strerror(errno)), exit(1);

  for (uint32 li=1; li<gkpStore->gkStore_getNumLibraries() + 1; li++) {
    mpaLibrary *ld = libdata[li];

    fprintf(stderr, "%s\n", ld->library->libraryName);

    ld->newSamples = 0;
    ld->newMean    = 0;
    ld->newStdDev  = 0;

    for (uint32 oi=0; oi<4; oi++) {
      vector<int32>   &dist = ld->dist[oi];
      uint32           np   = dist.size();

      if (np == 0)
        continue;

      sort(dist.begin(), dist.end());

      int32  median    = dist[np * 1 / 2];
      int32  oneThird  = dist[np * 1 / 3];
      int32  twoThird  = dist[np * 2 / 3];

      int32  approxStd = MAX(median - oneThird, twoThird - median);

      int32  biggest   = median + approxStd * 5;
      int32  smallest  = median - approxStd * 5;

      int64  newSamples = 0;
      double newMean    = 0;
      double newStdDev  = 0;

      for (uint64 x=0; x<dist.size(); x++)
        if ((smallest  <= dist[x]) &&
            (dist[x]   <= biggest)) {
          newSamples += 1;
          newMean    += dist[x];
        }

      if (newSamples == 0)
        continue;

      newMean   = newMean / newSamples;
      newStdDev = 0.0;

      for (uint64 x=0; x<dist.size(); x++)
        if ((smallest  <= dist[x]) &&
            (dist[x]   <= biggest))
          newStdDev += (dist[x] - newMean) * (dist[x] - newMean);
          
      if (newSamples > 1)
        newStdDev = sqrt(newStdDev / (newSamples - 1));

      if (((ld->library->orientation == AS_READ_ORIENT_INNIE)      && (oi == MPA_ORIENT_INNIE)) ||
          ((ld->library->orientation == AS_READ_ORIENT_OUTTIE)     && (oi == MPA_ORIENT_OUTTIE)) ||
          ((ld->library->orientation == AS_READ_ORIENT_NORMAL)     && (oi == MPA_ORIENT_NORMAL)) ||
          ((ld->library->orientation == AS_READ_ORIENT_ANTINORMAL) && (oi == MPA_ORIENT_ANTINORMAL)))
        fprintf(stderr, "%9.2f +- %8.2f -> %9.2f +- %8.2f %c "F_U64"/"F_U64" samples external happy "F_U64" sad "F_U64"\n",
                ld->origMean, ld->origStdDev,
                newMean,  newStdDev,
                MPA_ORIENT_NAMES[oi],
                newSamples,
                dist.size(),
                ld->externHappy[oi],
                ld->externSad[oi]);
      else
        fprintf(stderr, "%9s +- %8s -> %9.2f +- %8.2f %c "F_U64"/"F_U64" samples\n",
                "N/A", "N/A",
                newMean,  newStdDev,
                MPA_ORIENT_NAMES[oi],
                newSamples,
                dist.size());

      sprintf(datName, "%s.%03d.%c.%s.dat",
              prefix, li, MPA_ORIENT_NAMES[oi], ld->library->libraryName);
      errno = 0;
      FILE *dat = fopen(datName, "w");
      if (errno)
        fprintf(stderr, "ERROR:  Failed to open '%s': %s\n", datName, strerror(errno)), exit(1);

      int32  step = (biggest - smallest) / 100;
      int32  min  = dist[0];
      int32  max  = dist[0] + step;
      int32  num  = 0;

      for (uint64 x=0; x<dist.size(); x++) {
        if (dist[x] < max) {
          num++;
        } else {
          fprintf(dat, "%d\t%.1f\t%d\t%d\n", min, (min + max) / 2.0, max, num);
          min = max;
          max = min + step;
          num = 0;
        }
      }

      fclose(dat);
    
      fprintf(gp, "\n");
      fprintf(gp, "set terminal png\n");
      fprintf(gp, "set output \"%s.png\"\n", datName);
      fprintf(gp, "set title \"%s bucketSize=%d\"\n", ld->library->libraryName, step);
      fprintf(gp, "plot \"%s\" using 2:4 with lines title \"%s %c\"\n", datName, ld->library->libraryName, MPA_ORIENT_NAMES[oi]);
      fprintf(gp, "\n");
      fprintf(gp, "\n");
      fprintf(gp, "\n");
      fprintf(gp, "\n");
      fprintf(gp, "\n");
      fprintf(gp, "\n");
      fprintf(gp, "\n");
    }
  }

  fclose(gp);

  sprintf(datName, "gnuplot < %s.gp", prefix);
  system(datName);
}
