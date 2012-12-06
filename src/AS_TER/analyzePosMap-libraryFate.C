
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

static const char *rcsid = "$Id: analyzePosMap-libraryFate.C,v 1.2 2012-12-06 01:51:57 brianwalenz Exp $";

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

  if (dist.size() == 0) {
    fprintf(stderr, "NO SAMPLES\n");
    return;
  }

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

  if (numSamples == 0) {
    fprintf(stderr, "NO SAMPLES in range %f - %f\n", smallest, biggest);
    return;
  }

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
  };


  void       report(char *prefix, uint32 libid, char *libname, readFate &allfate) {
    char   N[FILENAME_MAX];
    FILE  *F;

    sprintf(N, "%s.readFate.lib%02d.%s", prefix, libid, libname);

    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Couldn't open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

    fprintf(F, "\n");
    fprintf(F, "libFate\t%s\n", libname);
    fprintf(F, "reads\t"F_U32"\n", total);
    fprintf(F, "\n");
    fprintf(F, "category\tcount\tlibFrac\tcatFrac\n");
    fprintf(F, "UNPLACED\n");
    fprintf(F, "deleted\t"F_U32"\t%.2f\t%.2f\n",
            frgDeleted,
            frgDeleted * 100.0 / total,
            frgDeleted * 100.0 / allfate.frgDeleted);
    fprintf(F, "singleton\t"F_U32"\t%.2f\t%.2f\n",
            frgSingleton,
            frgSingleton * 100.0 / total,
            frgSingleton * 100.0 / allfate.frgSingleton);
    fprintf(F, "UNITIGS\n");
    fprintf(F, "degenerate\t"F_U32"\t%.2f\t%.2f\n",
            ctgDegenerate,
            ctgDegenerate * 100.0 / total,
            ctgDegenerate * 100.0 / allfate.ctgDegenerate);
    fprintf(F, "unique\t"F_U32"\t%.2f\t%.2f\n",
            utgUnique,
            utgUnique * 100.0 / total,
            utgUnique * 100.0 / allfate.utgUnique);
    fprintf(F, "rock\t"F_U32"\t%.2f\t%.2f\n",
            utgRock,
            utgRock * 100.0 / total,
            utgRock * 100.0 / allfate.utgRock);
    fprintf(F, "stone\t"F_U32"\t%.2f\t%.2f\n",
            utgStone,
            utgStone * 100.0 / total,
            utgStone * 100.0 / allfate.utgStone);
    fprintf(F, "stoneUnres\t"F_U32"\t%.2f\t%.2f\n",
            utgStoneUnresolved,
            utgStoneUnresolved * 100.0 / total,
            utgStoneUnresolved * 100.0 / allfate.utgStoneUnresolved);
    fprintf(F, "pebble\t"F_U32"\t%.2f\t%.2f\n",
            utgPebble,
            utgPebble * 100.0 / total,
            utgPebble * 100.0 / allfate.utgPebble);
    fprintf(F, "CONTIGS\n");
    fprintf(F, "degenerate\t"F_U32"\t%.2f\t%.2f\n",
            ctgDegenerate,
            ctgDegenerate * 100.0 / total,
            ctgDegenerate * 100.0 / allfate.ctgDegenerate);
    fprintf(F, "big\t"F_U32"\t%.2f\t%.2f\n",
            ctgBig,
            ctgBig * 100.0 / total,
            ctgBig * 100.0 / allfate.ctgBig);
    fprintf(F, "small\t"F_U32"\t%.2f\t%.2f\n",
            ctgSmall,
            ctgSmall * 100.0 / total,
            ctgSmall * 100.0 / allfate.ctgSmall);

    fclose(F);
  };

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

    memset(sameScf, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);
    memset(diffScf, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);

    memset(sameCtg, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);
    memset(diffCtg, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);

    memset(sameUtg, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);
    memset(diffUtg, 0, sizeof(uint32) * NUM_FRG_LABELS * NUM_FRG_LABELS);
  };
  ~mateFate() {};

  mateFate &operator+=(mateFate const &that) {
    total += that.total;

    for (uint32 ii=0; ii<NUM_FRG_LABELS; ii++) {
      for (uint32 jj=0; jj<NUM_FRG_LABELS; jj++) {
        sameScf[ii][jj] += that.sameScf[ii][jj];
        diffScf[ii][jj] += that.diffScf[ii][jj];

        sameCtg[ii][jj] += that.sameCtg[ii][jj];
        diffCtg[ii][jj] += that.diffCtg[ii][jj];

        sameUtg[ii][jj] += that.sameUtg[ii][jj];
        diffUtg[ii][jj] += that.diffUtg[ii][jj];
      }
    }

    return(*this);
  };

  void       reportDetails(FILE    *F,
                           uint32   stat[NUM_FRG_LABELS][NUM_FRG_LABELS]) {
    fprintf(F, ".");
    for (uint32 jj=0; jj<NUM_FRG_LABELS; jj++)
      fprintf(F, "\t%s", frgLabelS[jj]);
    fprintf(F, "\n");

    for (uint32 ii=0; ii<NUM_FRG_LABELS; ii++) {
      fprintf(F, "%s", frgLabelS[ii]);
      for (uint32 jj=0; jj<NUM_FRG_LABELS; jj++)
        if (ii <= jj) {
          fprintf(F, "\t"F_U32, stat[ii][jj]);
        } else {
          fprintf(F, "\t");
          assert(stat[ii][jj] == 0);
        }
      fprintf(F, "\n");
    }
  };


  void       report(char *prefix, uint32 libid, char *libname, mateFate &allfate) {
    char   N[FILENAME_MAX];
    FILE  *F;

    sprintf(N, "%s.mateFate.lib%02d.%s", prefix, libid, libname);

    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Couldn't open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

    fprintf(F, "\n");
    fprintf(F, "libFate\t%s\n", libname);
    fprintf(F, "mates\t"F_U32"\n", total);
    fprintf(F, "\n");

    fprintf(F, "\n");
    fprintf(F, "SAME SCAFFOLD MATE COUNTS\n");
    reportDetails(F, sameScf);

    fprintf(F, "\n");
    fprintf(F, "DIFF SCAFFOLD MATE COUNTS\n");
    reportDetails(F, diffScf);

    fprintf(F, "\n");
    fprintf(F, "SAME CONTIG MATE COUNTS\n");
    reportDetails(F, sameCtg);

    fprintf(F, "\n");
    fprintf(F, "DIFF CONTIG MATE COUNTS\n");
    reportDetails(F, diffCtg);


    fprintf(F, "\n");
    fprintf(F, "SAME UNITIG MATE COUNTS\n");
    reportDetails(F, sameUtg);

    fprintf(F, "\n");
    fprintf(F, "DIFF UNITIG MATE COUNTS\n");
    reportDetails(F, diffUtg);
  };

  uint32    total;

  uint32    sameScf[NUM_FRG_LABELS][NUM_FRG_LABELS];
  uint32    diffScf[NUM_FRG_LABELS][NUM_FRG_LABELS];

  uint32    sameCtg[NUM_FRG_LABELS][NUM_FRG_LABELS];
  uint32    diffCtg[NUM_FRG_LABELS][NUM_FRG_LABELS];

  uint32    sameUtg[NUM_FRG_LABELS][NUM_FRG_LABELS];
  uint32    diffUtg[NUM_FRG_LABELS][NUM_FRG_LABELS];
};



class insertSize {
public:
  insertSize() {
    meanInnie  = sdevInnie  = 0.0;
    meanOuttie = sdevOuttie = 0.0;
    meanSame   = sdevSame   = 0.0;
  };
  ~insertSize() {
  };

  void  addDistance(uint32 A, int32 Abgn, int32 Aend, char Aori,
                    uint32 B, int32 Bbgn, int32 Bend, char Bori) {

    int32  dist = MAX(Aend, Bend) - MIN(Abgn, Bbgn);

    if (dist < 0) {
      fprintf(stderr, "negative distance %d\n", dist);
      fprintf(stderr, "A="F_U32" bgn="F_U32" end="F_U32" ori=%c %s\n",
              A, Abgn, Aend, Aori, frgNam[A].c_str());
      fprintf(stderr, "B="F_U32" bgn="F_U32" end="F_U32" ori=%c %s\n",
              B, Bbgn, Bend, Bori, frgNam[B].c_str());
      dist = 0;
    }
    assert(dist >= 0);

    if        (((Aori == Bori)))
      distSame.push_back(dist);

    else if (((Abgn <=  Bbgn) && (Aori == 'f')) ||
             ((Abgn >=  Bbgn) && (Aori == 'r')))
      distInnie.push_back(dist);

    else if (((Abgn <  Bbgn) && (Aori == 'r')) ||
             ((Abgn >  Bbgn) && (Aori == 'f')))
      distOuttie.push_back(dist);

    else {
      fprintf(stderr, "invalid position and orientation\n");
      fprintf(stderr, "A="F_U32" bgn="F_U32" end="F_U32" ori=%c %s\n",
              A, Abgn, Aend, Aori, frgNam[A].c_str());
      fprintf(stderr, "B="F_U32" bgn="F_U32" end="F_U32" ori=%c %s\n",
              B, Bbgn, Bend, Bori, frgNam[B].c_str());
      //assert(0);
    }
  };

  void  finalize(void) {

    fprintf(stderr, "FINALIZE %ld %ld %ld\n", distInnie.size(), distOuttie.size(), distSame.size());
    computeStdDev(distInnie,  meanInnie,  sdevInnie);
    computeStdDev(distOuttie, meanOuttie, sdevOuttie);
    computeStdDev(distSame,   meanSame,   sdevSame);
  };

  vector<double>   distInnie;
  double           meanInnie;
  double           sdevInnie;

  vector<double>   distOuttie;
  double           meanOuttie;
  double           sdevOuttie;

  vector<double>   distSame;
  double           meanSame;
  double           sdevSame;
};



class sizeFate {
public:
  sizeFate() {
  };
  ~sizeFate() {
  };

  void       report(char *prefix, uint32 libid, char *libname, sizeFate &allfate) {
    char   N[FILENAME_MAX];
    FILE  *F;

    sprintf(N, "%s.sizeFate.lib%02d.%s", prefix, libid, libname);

    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Couldn't open '%s' for writing: %s\n", N, strerror(errno)), exit(1);

    scf.finalize();
    ctg.finalize();
    utg.finalize();

    fprintf(F, "SCAFFOLD\n");
    fprintf(F, "innie       "F_SIZE_T" %.2f +- %.2f bp\n", scf.distInnie.size(),  scf.meanInnie,  scf.sdevInnie);
    fprintf(F, "outtie      "F_SIZE_T" %.2f +- %.2f bp\n", scf.distOuttie.size(), scf.meanOuttie, scf.sdevOuttie);
    fprintf(F, "same        "F_SIZE_T" %.2f +- %.2f bp\n", scf.distSame.size(),   scf.meanSame,   scf.sdevSame);

    fprintf(F, "CONTIG\n");
    fprintf(F, "innie       "F_SIZE_T" %.2f +- %.2f bp\n", ctg.distInnie.size(),  ctg.meanInnie,  ctg.sdevInnie);
    fprintf(F, "outtie      "F_SIZE_T" %.2f +- %.2f bp\n", ctg.distOuttie.size(), ctg.meanOuttie, ctg.sdevOuttie);
    fprintf(F, "same        "F_SIZE_T" %.2f +- %.2f bp\n", ctg.distSame.size(),   ctg.meanSame,   ctg.sdevSame);

    fprintf(F, "UNITIG\n");
    fprintf(F, "innie       "F_SIZE_T" %.2f +- %.2f bp\n", utg.distInnie.size(),  utg.meanInnie,  utg.sdevInnie);
    fprintf(F, "outtie      "F_SIZE_T" %.2f +- %.2f bp\n", utg.distOuttie.size(), utg.meanOuttie, utg.sdevOuttie);
    fprintf(F, "same        "F_SIZE_T" %.2f +- %.2f bp\n", utg.distSame.size(),   utg.meanSame,   utg.sdevSame);

    fclose(F);
  };

  insertSize  scf;
  insertSize  ctg;
  insertSize  utg;
};



class libraryFate {
public:
  libraryFate() {
  };
  ~libraryFate() {};

  readFate    read;
  mateFate    mate;
  sizeFate    size;
};




void
analyzeLibraryFate(char *outPrefix) {
  libraryFate  *libfate = new libraryFate [libDat.size()];
  libraryFate   allfate;

  //  Read Fate

  for (uint32 fi=0; fi<frgDat.size(); fi++) {
    uint32 li = frgLibrary[fi];
 
    libfate[li].read.total++;

    if (frgDat[fi].sta == 'D') {
      libfate[li].read.frgDeleted++;
      continue;
    }

    if ((frgDat[fi].sta == 'C') ||
        (frgDat[fi].sta == 'S')) {
      libfate[li].read.frgSingleton++;
      continue;
    }

    if (frgDat[fi].sta == 'd') {
      libfate[li].read.ctgDegenerate++;
      continue;
    }

    if (frgDat[fi].sta == 'R') {
      libfate[li].read.utgStoneUnresolved++;
      continue;
    }

    if (frgDat[fi].sta == 'o') {
      continue;
    }

    //  Otherwise, the read is in a untig/contig we care about.

    if      (frgDat[fi].sta == 'u')   libfate[li].read.utgUnique++;
    else if (frgDat[fi].sta == 'r')   libfate[li].read.utgRock++;
    else if (frgDat[fi].sta == 's')   libfate[li].read.utgStone++;
    else if (frgDat[fi].sta == 'p')   libfate[li].read.utgPebble++;
    else assert(0);

    if (ctgDat[frgDat[fi].ctg.con].len >= 10000)
      libfate[li].read.ctgBig++;
    else
      libfate[li].read.ctgSmall++;
  }

  for (uint32 li=0; li<libDat.size(); li++)
    allfate.read += libfate[li].read;

  allfate.read.report(outPrefix, 0, "ALLLIBRARIES", allfate.read);

  for (uint32 li=0; li<libDat.size(); li++)
    libfate[li].read.report(outPrefix, li+1, libDat[li].name, allfate.read);


  //  Mate fate


  for (uint32 fi=0; fi<frgMate.size(); fi++) {
    uint32 mi = frgMate[fi];
    uint32 li = frgLibrary[fi];

    if (mi <= fi)
      //  Either not mated, or visited already.
      continue;

    uint32 asta = frgLabelToInt[frgDat[fi].sta];
    uint32 bsta = frgLabelToInt[frgDat[mi].sta];

    if (asta > bsta) {
      asta = frgLabelToInt[frgDat[mi].sta];
      bsta = frgLabelToInt[frgDat[fi].sta];
    }

    libfate[li].mate.total++;

    if (frgDat[fi].scf.con == frgDat[mi].scf.con)
      libfate[li].mate.sameScf[asta][bsta]++;
    else
      libfate[li].mate.diffScf[asta][bsta]++;

    if (frgDat[fi].ctg.con == frgDat[mi].ctg.con)
      libfate[li].mate.sameCtg[asta][bsta]++;
    else
      libfate[li].mate.diffCtg[asta][bsta]++;

    if (frgDat[fi].utg.con == frgDat[mi].utg.con)
      libfate[li].mate.sameUtg[asta][bsta]++;
    else
      libfate[li].mate.diffUtg[asta][bsta]++;
  }


  for (uint32 li=0; li<libDat.size(); li++)
    allfate.mate += libfate[li].mate;

  allfate.mate.report(outPrefix, 0, "ALLLIBRARIES", allfate.mate);

  for (uint32 li=0; li<libDat.size(); li++)
    libfate[li].mate.report(outPrefix, li+1, libDat[li].name, allfate.mate);


  //  Insert sizes


  for (uint32 fi=0; fi<frgMate.size(); fi++) {
    uint32 mi = frgMate[fi];
    uint32 li = frgLibrary[fi];

    if (mi <= fi)
      //  Either not mated, or visited already.
      continue;

    uint32 fisi = frgDat[fi].scf.con;
    uint32 misi = frgDat[mi].scf.con;

    if ((fisi == UINT32_MAX) || (misi == UINT32_MAX))
      continue;

    if (fisi != misi)
      continue;

    if (scfDat[fisi].sta == 'd')
      continue;

    assert(scfDat[fisi].sta == 'p');

    libfate[li].size.scf.addDistance(fi, frgDat[fi].scf.bgn, frgDat[fi].scf.end, frgDat[fi].scf.ori,
                                     mi, frgDat[mi].scf.bgn, frgDat[mi].scf.end, frgDat[mi].scf.ori);
  }

  //for (uint32 li=0; li<libDat.size(); li++)
  //  allfate.size += libfate[li].size;

  allfate.size.report(outPrefix, 0, "ALLLIBRARIES", allfate.size);

  for (uint32 li=0; li<libDat.size(); li++)
    libfate[li].size.report(outPrefix, li+1, libDat[li].name, allfate.size);
}
