
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

static const char *rcsid = "$Id: MultiAlignMatePairAnalysis.C,v 1.3 2011-12-19 00:51:47 brianwalenz Exp $";

#include "MultiAlignMatePairAnalysis.H"

#include <math.h>

#include <map>
#include <set>
#include <vector>
#include <algorithm>

using namespace std;


class mpaLibraryData {
public:
  mpaLibraryData() {
    median      = 0;
    oneThird    = 0;
    twoThird    = 0;

    approxStd   = 0;

    biggest     = 0;
    smallest    = 0;

    numSamples  = 0;
    mean        = 0;
    stddev      = 0;

    externHappy = 0;
    externSad   = 0;
  };

  ~mpaLibraryData() {
  };

  void    finalize(void);
  void    writeUpdate(FILE *output, int32 libOrient, gkLibrary *library);
  void    printSummary(FILE *output, int32 libOrient, gkLibrary *library);

  int32   median;
  int32   oneThird;
  int32   twoThird;

  double  approxStd;

  int32   biggest;
  int32   smallest;

  uint64  numSamples;
  double  mean;
  double  stddev;

  //  For externally mated fragments, happy?
  uint64  externHappy;
  uint64  externSad;

  //  Distances, one per pair orientation
  vector<int32>   dist;
};


class mpaLibrary {
public:
  gkLibrary      *library;
  mpaLibraryData  orient[5];
};



matePairAnalysis::matePairAnalysis(char *gkpName) {

  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  libdata  = new mpaLibrary [gkpStore->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkpStore->gkStore_getNumLibraries() + 1; i++)
    libdata[i].library = gkpStore->gkStore_getLibrary(i);

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
  delete [] libdata;
  delete [] pairing;
  delete [] library;
  delete    gkpStore;
}





void
matePairAnalysis::evaluateTig(MultiAlignT *ma) {
  map<AS_IID, uint32>  fragIdx;
  uint32               numFrag = GetNumIntMultiPoss(ma->f_list);
  int32                maLen   = 0;

  //  We could use GetMultiAlignLength(ma), but we need to go throught the list of
  //  fragments anyway (to build the fragIdx map) and we'll just find the length
  //  at the same time.
  //
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
    mpaLibrary *lib       = libdata + libiid;
    uint32      liborient = lib->library->orientation;
    int32       libmin    = lib->library->mean - 3 * lib->library->stddev;
    int32       libmax    = lib->library->mean + 3 * lib->library->stddev;

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
          if (frgforward == true)   ((frgmin < maLen - libmax) ? lib->orient[AS_READ_ORIENT_INNIE].externSad : lib->orient[AS_READ_ORIENT_INNIE].externHappy)++;
          if (frgforward == false)  ((frgmin >         libmax) ? lib->orient[AS_READ_ORIENT_INNIE].externSad : lib->orient[AS_READ_ORIENT_INNIE].externHappy)++;
          break;

        case AS_READ_ORIENT_OUTTIE:
          if (frgforward == false)   ((frgmin < maLen - libmax) ? lib->orient[AS_READ_ORIENT_OUTTIE].externSad : lib->orient[AS_READ_ORIENT_OUTTIE].externHappy)++;
          if (frgforward == true)    ((frgmin >         libmax) ? lib->orient[AS_READ_ORIENT_OUTTIE].externSad : lib->orient[AS_READ_ORIENT_OUTTIE].externHappy)++;
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
      lib->orient[AS_READ_ORIENT_INNIE].dist.push_back(distance);

    } else if ((frgforward == false) && (matforward == true)) {
      lib->orient[AS_READ_ORIENT_OUTTIE].dist.push_back(distance);

    } else if ((frgforward == true)  && (matforward == true)) {
      lib->orient[AS_READ_ORIENT_NORMAL].dist.push_back(distance);

    } else {//((frgforward == false) && (matforward == false))
      lib->orient[AS_READ_ORIENT_ANTINORMAL].dist.push_back(distance);
    }
  }
}




void
mpaLibraryData::finalize(void) {

  if (dist.size() == 0)
    return;

  sort(dist.begin(), dist.end());

  median    = dist[dist.size() * 1 / 2];
  oneThird  = dist[dist.size() * 1 / 3];
  twoThird  = dist[dist.size() * 2 / 3];

  approxStd = MAX(median - oneThird, twoThird - median);

  biggest   = median + approxStd * 5;
  smallest  = median - approxStd * 5;

  for (uint64 x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest)) {
      numSamples += 1;
      mean    += dist[x];
    }

  if (numSamples == 0)
    return;

  mean   = mean / numSamples;
  stddev = 0.0;

  for (uint64 x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest))
      stddev += (dist[x] - mean) * (dist[x] - mean);
          
  if (numSamples > 1)
    stddev = sqrt(stddev / (numSamples - 1));
}


void
matePairAnalysis::finalize(void) {
  for (uint32 li=1; li<gkpStore->gkStore_getNumLibraries() + 1; li++) {
    libdata[li].orient[0].finalize();  //  Unoriented
    libdata[li].orient[1].finalize();  //  Innie
    libdata[li].orient[2].finalize();  //  Outtie
    libdata[li].orient[3].finalize();  //  Normal
    libdata[li].orient[4].finalize();  //  Antinormal
  }
}



void
mpaLibraryData::printSummary(FILE *output, int32 libOrient, gkLibrary *library) {
  if (dist.size() == 0)
    return;

  if (library->orientation == libOrient)
    fprintf(stderr, "%9.2f +- %8.2f -> %9.2f +- %8.2f %s "F_U64"/"F_U64" samples external happy "F_U64" sad "F_U64"\n",
            library->mean, library->stddev,
            mean,  stddev,
            AS_READ_ORIENT_NAMES[libOrient],
            numSamples,
            dist.size(),
            externHappy,
            externSad);
    else
      fprintf(stderr, "%9s +- %8s -> %9.2f +- %8.2f %s "F_U64"/"F_U64" samples\n",
              "N/A", "N/A",
              mean,  stddev,
              AS_READ_ORIENT_NAMES[libOrient],
              numSamples,
              dist.size());
}



void
matePairAnalysis::printSummary(FILE *output) {
  for (uint32 li=1; li<gkpStore->gkStore_getNumLibraries() + 1; li++) {
    fprintf(output, "\n%s\n", libdata[li].library->libraryName);
    libdata[li].orient[0].printSummary(output, 0, libdata[li].library);  //  Unoriented
    libdata[li].orient[1].printSummary(output, 1, libdata[li].library);  //  Innie
    libdata[li].orient[2].printSummary(output, 2, libdata[li].library);  //  Outtie
    libdata[li].orient[3].printSummary(output, 3, libdata[li].library);  //  Normal
    libdata[li].orient[4].printSummary(output, 4, libdata[li].library);  //  Antinormal
  }
}


void
mpaLibraryData::writeUpdate(FILE *output, int32 libOrient, gkLibrary *library) {
  if (dist.size() == 0)
    return;

  if ((library->orientation == libOrient) &&
      (numSamples >= 100))
    fprintf(output, "lib uid %s distance %.2f %.2f\n", library->libraryName, mean, stddev);
}



void
matePairAnalysis::writeUpdate(char *prefix) {
  char   datName[FILENAME_MAX];

  if (prefix[strlen(prefix)-1] == '/')
    prefix[strlen(prefix)-1] = 0;

  sprintf(datName, "%s.distupdate", prefix);
  errno = 0;
  FILE *F = fopen(datName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  Failed to open '%s': %s\n", datName, strerror(errno)), exit(1);

  for (uint32 li=1; li<gkpStore->gkStore_getNumLibraries() + 1; li++) {
    libdata[li].orient[0].writeUpdate(F, 0, libdata[li].library);  //  Unoriented
    libdata[li].orient[1].writeUpdate(F, 1, libdata[li].library);  //  Innie
    libdata[li].orient[2].writeUpdate(F, 2, libdata[li].library);  //  Outtie
    libdata[li].orient[3].writeUpdate(F, 3, libdata[li].library);  //  Normal
    libdata[li].orient[4].writeUpdate(F, 4, libdata[li].library);  //  Antinormal
  }

  fclose(F);
}

void
matePairAnalysis::drawPlots(char *prefix) {
  char   datName[FILENAME_MAX];

  if (prefix[strlen(prefix)-1] == '/')
    prefix[strlen(prefix)-1] = 0;

  sprintf(datName, "%s.gp", prefix);
  errno = 0;
  FILE *gp = fopen(datName, "w");
  if (errno)
    fprintf(stderr, "ERROR:  Failed to open '%s': %s\n", datName, strerror(errno)), exit(1);


  for (uint32 li=1; li<gkpStore->gkStore_getNumLibraries() + 1; li++) {
    mpaLibrary *ld = libdata + li;

    for (uint32 oi=0; oi<4; oi++) {
      mpaLibraryData  *data = ld->orient + oi;

      if (data->dist.size() == 0)
        continue;

      sprintf(datName, "%s.%03d.%s.%s.dat",
              prefix, li, AS_READ_ORIENT_NAMES[oi], ld->library->libraryName);
      errno = 0;
      FILE *dat = fopen(datName, "w");
      if (errno)
        fprintf(stderr, "ERROR:  Failed to open '%s': %s\n", datName, strerror(errno)), exit(1);

      int32  step = (data->biggest - data->smallest) / 100;
      int32  min  = data->dist[0];
      int32  max  = data->dist[0] + step;
      int32  num  = 0;

      for (uint64 x=0; x<data->dist.size(); x++) {
        if (data->dist[x] < max) {
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
      fprintf(gp, "plot \"%s\" using 2:4 with lines title \"%s %s\"\n", datName, ld->library->libraryName, AS_READ_ORIENT_NAMES[oi]);
      fprintf(gp, "\n");
    }
  }

  fclose(gp);

  sprintf(datName, "gnuplot < %s.gp", prefix);
  system(datName);
}



double
matePairAnalysis::mean(uint32 libid) {
  mpaLibrary  *l = libdata + libid;
  return(l->orient[l->library->orientation].mean);
}



double
matePairAnalysis::stddev(uint32 libid) {
  mpaLibrary  *l = libdata + libid;
  return(l->orient[l->library->orientation].stddev);
}
