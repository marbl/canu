
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: analyzePosMap-libraryFate.C,v 1.1 2012-12-05 01:13:23 brianwalenz Exp $";

#include "analyzePosMap.H"

//  Describes read fate per library.  'placed-stone' means the read is in a stone, and was resolved to a contig.
//
//            --------------------CONTIGS----------------    ---------------UNITIGS--------------------   --READS--
//  LIBRARY   big-contig  small-contig  degenerate-contig    unqiue  rock  placed-stone  unplaced-stone   singleton
//

//  Describes mate fate per library.
//
//  LIBRARY   innie-short  innie-correct  innie-long
//            normal-short normal-correct normal-long
//            outtie-short outtie-correct outtie-long
//
//            extern-missing (should be in a contig in this scaffold)
//            extern-gap     (could be in a gap in this scaffold)
//            extern-end     (could be in a scaffold off the end)
//
//            typeRead1--typeRead2 -- unique--unique, unique--rock, unplaced-stone--singleton, etc.
//              No happiness here, just counts.
//



void
computeStdDev(vector<double> dist, double &mean, double &stddev) {
  mean   = 0.0;
  stddev = 0.0;

  if (dist.size() == 0)
    return;

  sort(dist.begin(), dist.end());

  double median     = dist[dist.size() * 1 / 2];
  double oneThird   = dist[dist.size() * 1 / 3];
  double twoThird   = dist[dist.size() * 2 / 3];

  double approxStd  = MAX(median - oneThird, twoThird - median);

  double biggest    = median + approxStd * 5;
  double smallest   = median - approxStd * 5;

  size_t numSamples = 0;

  for (size_t x=0; x<dist.size(); x++)
    if ((smallest  <= dist[x]) &&
        (dist[x]   <= biggest)) {
      numSamples += 1;
      mean       += dist[x];
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






class readFate {
public:
  readFate() {
    total               = 0;

    frgDeleted          = 0;
    frgSingleton        = 0;

    utgUnique           = 0;
    utgRock             = 0;
    utgStone            = 0;
    utgStoneUnresolved  = 0;
    utgPebble           = 0;
    utgSingleton        = 0;
    utgOther            = 0;

    ctgBig              = 0;
    ctgSmall            = 0;
    ctgDegenerate       = 0;
  };
  ~readFate() {};

  readFate &operator+=(readFate const &that) {
    total               += that.total;

    frgDeleted          += that.frgDeleted;
    frgSingleton        += that.frgSingleton;

    utgUnique           += that.utgUnique;
    utgRock             += that.utgRock;
    utgStone            += that.utgStone;
    utgStoneUnresolved  += that.utgStoneUnresolved;
    utgPebble           += that.utgPebble;
    utgSingleton        += that.utgSingleton;
    utgOther            += that.utgOther;

    ctgBig              += that.ctgBig;
    ctgSmall            += that.ctgSmall;
    ctgDegenerate       += that.ctgDegenerate;

    return(*this);
  }

  uint32   total;

  uint32   frgDeleted;
  uint32   frgSingleton;

  uint32   utgUnique;
  uint32   utgRock;
  uint32   utgStone;
  uint32   utgStoneUnresolved;
  uint32   utgPebble;
  uint32   utgSingleton;
  uint32   utgOther;

  uint32   ctgBig;
  uint32   ctgSmall;
  uint32   ctgDegenerate;
};


class mateFate {
public:
  mateFate() {
    total = 0;

    memset(sameScaffold, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);
    memset(diffScaffold, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);

    memset(scfDegenerate, 0, sizeof(uint32) * NUM_FRG_LABELS);
    memset(scfSingleton,  0, sizeof(uint32) * NUM_FRG_LABELS);
    memset(scfUnresolved, 0, sizeof(uint32) * NUM_FRG_LABELS);

    sameDegenerate = 0;
    diffDegenerate = 0;
  };
  ~mateFate() {};

  uint32    total;

  uint32    sameScaffold[NUM_FRG_LABELS][NUM_FRG_LABELS];
  uint32    diffScaffold[NUM_FRG_LABELS][NUM_FRG_LABELS];

  uint32    scfDegenerate[NUM_FRG_LABELS];
  uint32    scfSingleton[NUM_FRG_LABELS];
  uint32    scfUnresolved[NUM_FRG_LABELS];

  uint32    bothSingleton;
  uint32    bothUnresolved;

  uint32    sameDegenerate;
  uint32    diffDegenerate;
};



class insertSize {
public:
  insertSize() {
    mean   = 0;
    stddev = 0;
  };
  ~insertSize() {
  };

  void  addDistance(int32 dist) {
    distances.push_back(dist);
  };

  void  finalize(void) {
    computeStdDev(distances, mean, stddev);
  };

#if 0
  fprintf(O, "innie       "F_SIZE_T" %.2f +- %.2f bp\n", innie.distances.size(),  innie.mean,  innie.stddev);
  fprintf(O, "outtie      "F_SIZE_T" %.2f +- %.2f bp\n", outtie.distances.size(), outtie.mean, outtie.stddev);
  fprintf(O, "same        "F_SIZE_T" %.2f +- %.2f bp\n", same.distances.size(),   same.mean,   same.stddev);
#endif

  vector<double>   distances;
  double           mean;
  double           stddev;
};



class libraryFate {
public:
  libraryFate() {
  };
  ~libraryFate() {};

  readFate    read;
  mateFate    mate;

  insertSize  innie;
  insertSize  outtie;
  insertSize  same;
};




void
analyzeLibraryFate(char *outPrefix) {
  libraryFate  *libfate = new libraryFate [libDat.size()];
  libraryFate   allfate;

  //  Read Fate

  for (uint32 fi=0; fi<frgDat.size(); fi++) {
    uint32 ll = frgLibrary[fi];
 
    libfate[ll].read.total++;

    if (frgDat[fi].sta == 'D') {
      libfate[ll].read.frgDeleted++;
      continue;
    }

    if ((frgDat[fi].sta == 'C') ||
        (frgDat[fi].sta == 'S')) {
      libfate[ll].read.frgSingleton++;
      continue;
    }

    if (frgDat[fi].sta == 'd') {
      libfate[ll].read.ctgDegenerate++;
      continue;
    }

    if (frgDat[fi].sta == 'R') {
      libfate[ll].read.utgStoneUnresolved++;
      continue;
    }

    if (frgDat[fi].sta == 'o') {
      continue;
    }

    //  Otherwise, the read is in a untig/contig we care about.

    if      (frgDat[fi].sta == 'u')   libfate[ll].read.utgUnique++;
    else if (frgDat[fi].sta == 'r')   libfate[ll].read.utgRock++;
    else if (frgDat[fi].sta == 's')   libfate[ll].read.utgStone++;
    else if (frgDat[fi].sta == 'p')   libfate[ll].read.utgPebble++;
    else assert(0);

    if (ctgDat[frgDat[fi].ctg.con].len >= 10000)
      libfate[ll].read.ctgBig++;
    else
      libfate[ll].read.ctgSmall++;
  }

  for (uint32 li=0; li<libDat.size(); li++)
    allfate.read += libfate[li].read;

  for (uint32 li=0; li<libDat.size(); li++) {
    fprintf(stdout, "\n");
    fprintf(stdout, "libFate\t%s\n", libDat[li].name);
    fprintf(stdout, "reads\t"F_U32"\n", libfate[li].read.total);
    fprintf(stdout, "\n");
    fprintf(stdout, "category\tcount\tlibFrac\tcatFrac\n");
    fprintf(stdout, "UNPLACED\n");
    fprintf(stdout, "deleted\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.frgDeleted,
            libfate[li].read.frgDeleted * 100.0 / libfate[li].read.total,
            libfate[li].read.frgDeleted * 100.0 / allfate.read.frgDeleted);
    fprintf(stdout, "singleton\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.frgSingleton,
            libfate[li].read.frgSingleton * 100.0 / libfate[li].read.total,
            libfate[li].read.frgSingleton * 100.0 / allfate.read.frgSingleton);
    fprintf(stdout, "UNITIGS\n");
    fprintf(stdout, "degenerate\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.ctgDegenerate,
            libfate[li].read.ctgDegenerate * 100.0 / libfate[li].read.total,
            libfate[li].read.ctgDegenerate * 100.0 / allfate.read.ctgDegenerate);
    fprintf(stdout, "unique\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.utgUnique,
            libfate[li].read.utgUnique * 100.0 / libfate[li].read.total,
            libfate[li].read.utgUnique * 100.0 / allfate.read.utgUnique);
    fprintf(stdout, "rock\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.utgRock,
            libfate[li].read.utgRock * 100.0 / libfate[li].read.total,
            libfate[li].read.utgRock * 100.0 / allfate.read.utgRock);
    fprintf(stdout, "stone\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.utgStone,
            libfate[li].read.utgStone * 100.0 / libfate[li].read.total,
            libfate[li].read.utgStone * 100.0 / allfate.read.utgStone);
    fprintf(stdout, "stoneUnres\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.utgStoneUnresolved,
            libfate[li].read.utgStoneUnresolved * 100.0 / libfate[li].read.total,
            libfate[li].read.utgStoneUnresolved * 100.0 / allfate.read.utgStoneUnresolved);
    fprintf(stdout, "pebble\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.utgPebble,
            libfate[li].read.utgPebble * 100.0 / libfate[li].read.total,
            libfate[li].read.utgPebble * 100.0 / allfate.read.utgPebble);
    fprintf(stdout, "CONTIGS\n");
    fprintf(stdout, "degenerate\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.ctgDegenerate,
            libfate[li].read.ctgDegenerate * 100.0 / libfate[li].read.total,
            libfate[li].read.ctgDegenerate * 100.0 / allfate.read.ctgDegenerate);
    fprintf(stdout, "big\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.ctgBig,
            libfate[li].read.ctgBig * 100.0 / libfate[li].read.total,
            libfate[li].read.ctgBig * 100.0 / allfate.read.ctgBig);
    fprintf(stdout, "small\t"F_U32"\t%.2f\t%.2f\n",
            libfate[li].read.ctgSmall,
            libfate[li].read.ctgSmall * 100.0 / libfate[li].read.total,
            libfate[li].read.ctgSmall * 100.0 / allfate.read.ctgSmall);
  }


  //  Mate fate


  for (uint32 fi=0; fi<frgMate.size(); fi++) {
    uint32 mi = frgMate[fi];

    if (mi <= fi)
      //  Either not mated, or visited already.
      continue;

    assert(frgDat[fi].typ == 'f');
    assert(frgDat[mi].typ == 'f');

    uint32  fisi = frgDat[fi].scf.con;
    uint32  misi = frgDat[mi].scf.con;

    assert(fisi < scfDat.size());
    assert(misi < scfDat.size());

#if 0
    uint32    sameScaffold[NUM_FRG_LABELS][NUM_FRG_LABELS];
    uint32    diffScaffold[NUM_FRG_LABELS][NUM_FRG_LABELS];

    uint32    scfDegenerate[NUM_FRG_LABELS];
    uint32    scfSingleton[NUM_FRG_LABELS];
    uint32    scfUnresolved[NUM_FRG_LABELS];

    uint32    bothSingleton;
    uint32    bothUnresolved;

    uint32    sameDegenerate;
    uint32    diffDegenerate;
#endif

    //if (fisi == UINT32_MAX)


#if 0
      if ((frgDat[fi].sta == 'D') ||  //  deleted
          (frgDat[fi].sta == 'C') ||  //  chaff
          (frgDat[fi].sta == 'd') ||  //  degenerate
          (frgDat[fi].sta == 'R'))    //  unplaced surrogate read
        //  Not in a scaffold, fragment deleted or singleton or unplaced surrogate.
        continue;

    if ((mi <= fi) && (frgDat[mi].sta != 'C'))
      //  Catches both unmated (mi==0) and visited (mi < fi), in particular,
      //  unmated fi=0 is caught.  The test on status is to not skip mates where the
      //  lower id is a singleton.
      continue;


    //assert((scfDat[si].sta == 'p') || (scfDat[si].sta == 'd'));

    if (fisi == misi) {
    } else {
    }



    //  What type of pair is this?

    if (frgDat[fi].con == frgDat[mi].con) {
      //  Same Scaffold


      if (scfDat[si].sta == 'p')
        scfStat = scfStatSameScaffold;    //  Real scaffolds!
      else
        scfStat = scfStatSameDegenerate;  //  Degenerates.

      scfStat->number++;

      int32  dist = MAX(frgDat[fi].end, frgDat[mi].end) - MIN(frgDat[fi].bgn, frgDat[mi].bgn);

      if (dist < 0) {
        fprintf(stderr, "negative distance %d\n", dist);
        fprintf(stderr, "fi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c sta=%c %s\n",
                fi,
                frgDat[fi].con, frgDat[fi].bgn, frgDat[fi].end, frgDat[fi].len, frgDat[fi].ori, frgDat[fi].sta, frgNam[fi].c_str());
        fprintf(stderr, "mi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c sta=%c %s\n",
                mi,
                frgDat[mi].con, frgDat[mi].bgn, frgDat[mi].end, frgDat[mi].len, frgDat[mi].ori, frgDat[mi].sta, frgNam[mi].c_str());
        dist = 0;
      }
      assert(dist >= 0);

      if        (((frgDat[fi].ori == frgDat[mi].ori))) {
        scfStat->same.addDistance(dist);
      }

      else if (((frgDat[fi].bgn <=  frgDat[mi].bgn) && (frgDat[fi].ori == 'f')) ||
               ((frgDat[fi].bgn >=  frgDat[mi].bgn) && (frgDat[fi].ori == 'r'))) {
        scfStat->innie.addDistance(dist);
      }

      else if (((frgDat[fi].bgn <  frgDat[mi].bgn) && (frgDat[fi].ori == 'r')) ||
               ((frgDat[fi].bgn >  frgDat[mi].bgn) && (frgDat[fi].ori == 'f'))) {
        scfStat->outtie.addDistance(dist);
      }

      else {
        fprintf(stderr, "invalid position and orientation\n");
        fprintf(stderr, "fi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c sta=%c %s\n",
                fi,
                frgDat[fi].con, frgDat[fi].bgn, frgDat[fi].end, frgDat[fi].len, frgDat[fi].ori, frgDat[fi].sta, frgNam[fi].c_str());
        fprintf(stderr, "mi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c sta=%c %s\n",
                mi,
                frgDat[mi].con, frgDat[mi].bgn, frgDat[mi].end, frgDat[mi].len, frgDat[mi].ori, frgDat[mi].sta, frgNam[mi].c_str());
        //assert(0);
      }
    }
#endif
  }
}

