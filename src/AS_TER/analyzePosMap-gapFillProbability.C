
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

static const char *rcsid = "$Id: analyzePosMap-gapFillProbability.C,v 1.1 2012-12-05 01:13:23 brianwalenz Exp $";

#include "analyzePosMap.H"



class gapFill {
public:
  gapFill() {
    numPossible    = 0;
    expectedFilled = 0;

    for (uint32 i=0; i<1000; i++)
      probability[i] = 0.0;
  };
  ~gapFill() {
  };

  //  Compute the expected placement of the mate, compute probability that it
  //  will fall into a gap.

  void   addGapFill(uint32 fi) {
    uint32  si = frgDat[fi].scf.con;
    uint32  li = frgLibrary[fi];

    assert(scfDat[si].sta == 'p');

    if (li >= libDat.size())
      fprintf(stderr, "WARNING: fragment %s claims library "F_U32", invalid\n", frgNam[fi].c_str(), li);
    assert(li < libDat.size());

    int32 fbgn = frgDat[fi].scf.bgn;
    int32 fend = frgDat[fi].scf.end;

    int32 mbgn = 0;
    int32 mend = 0;

    if (frgDat[fi].scf.ori == 'f') {
      mbgn = fbgn + (libDat[li].mean - 3 * libDat[li].stddev);
      mend = fbgn + (libDat[li].mean + 3 * libDat[li].stddev);
    } else {
      mbgn = fend - (libDat[li].mean + 3 * libDat[li].stddev);
      mend = fend - (libDat[li].mean - 3 * libDat[li].stddev);
    }

    if (mbgn < 0)                mbgn = 0;
    if (mbgn > scfDat[si].len)   mbgn = scfDat[si].len;

    if (mend < 0)                mend = 0;
    if (mend > scfDat[si].len)   mend = scfDat[si].len;

    if (mend <= mbgn)
      return;

    double pgapfill = 0;

    for (int32 i=mbgn; i<mend; i++)
      if (scfDat[si].base[i] == 0)
        pgapfill += libDat[li].pdf[i-mbgn];

    if (pgapfill < 0.001)
      return;

    numPossible++;
    
    int32 pg = (int)floor(pgapfill * 1000);

    if (pg < 0)      pg = 0;
    if (pg >= 1000)  pg = 1000;

    assert(0  <= pg);
    assert(pg <  1000);

    probability[pg]++;
  };

  void   finalize(void) {
    expectedFilled = 0;

    for (uint32 i=0; i<1000; i++)
      expectedFilled += (i / 1000.0) * probability[i];
  };

  uint32     numPossible;
  double     expectedFilled;
  double     probability[1000];
};



class scaffoldStatistics {
public:
  scaffoldStatistics(char *label_) {
    strcpy(label, label_);

    number = 0;
  };
  ~scaffoldStatistics() {
  };

  void         finalize(void) {
    fill.finalize();
  };

  void         eraseReport(char *outPrefix) {
    char  C[FILENAME_MAX];
    FILE *O;
    FILE *F;

    sprintf(C, "%s.gapfill.summary", outPrefix);

    errno = 0;
    O = fopen(C, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", C, strerror(errno)), exit(1);

    fclose(O);
  };


  void         printReport(char *outPrefix) {
    char  C[FILENAME_MAX];
    FILE *O;
    FILE *F;

    sprintf(C, "%s.gapfill.summary", outPrefix);

    errno = 0;
    O = fopen(C, "a+");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", C, strerror(errno)), exit(1);

    fprintf(O, "----------------------------------------\n");
    fprintf(O, "%s\n", label);
    fprintf(O, "total       "F_U32" mate pairs involved\n", number);
    fprintf(O, "expect fill %.0f gaps from "F_U32" mates\n", floor(fill.expectedFilled), fill.numPossible);

    if (fill.numPossible > 0) {
      sprintf(C, "%s.gapfill.%s", outPrefix, label);

      errno = 0;
      F = fopen(C, "w");
      if (errno)
        fprintf(stderr, "Failed to open '%s' for writing: %s\n", C, strerror(errno)), exit(1);

      for (uint32 i=0; i<1000; i++)
        fprintf(F, "%.3f\t%f\n", (double)i / 1000.0, fill.probability[i]);
      fclose(F);
    }

    fclose(O);
  };

  char         label[1024];

  uint32       number;
  gapFill      fill;
};





//  For each mate pair, compute stats.  per library?
//
//  same scaffold (normal/degenerate)
//  same contig (normal/degenerate)
//  different scaffold (normal-normal, etc)
//  different contig (normal-normal, etc -- special case for same scaffold)
//
//  gap fill probability
//   - for each external mate, compute the probability the frag belongs in a gap, plot histogram
//
void
analyzeGapFillProbability(char *outPrefix) {
  scaffoldStatistics   *scfStatDiffScaffold    = new scaffoldStatistics("diffScaffold");
  scaffoldStatistics   *scfStatDiffScafDegen   = new scaffoldStatistics("diffScafDegen");
  scaffoldStatistics   *scfStatDiffScafSing    = new scaffoldStatistics("diffScafSing");
  scaffoldStatistics   *scfStatDiffScafSurr    = new scaffoldStatistics("diffScafSurr");

  scaffoldStatistics   *scfStat                = NULL;

  //  Fragment status is more-or-less the unitig status.  We care about:
  //    s - in a scaffold   -- u, r
  //    d - in a degenerate -- d
  //    S - in a surrogate  -- s R p
  //    C - in a singleton  -- S

  char  mapFrgStaToScfSta[256] = { 0 };

  mapFrgStaToScfSta['D'] = 0;    //  Deleted read
  mapFrgStaToScfSta['C'] = 'C';  //  Chaff read
  mapFrgStaToScfSta['d'] = 'd';  //  Degenerate contig
  mapFrgStaToScfSta['u'] = 's';  //  Unique unitig
  mapFrgStaToScfSta['r'] = 's';  //  Rock unitig
  mapFrgStaToScfSta['s'] = 'S';  //  Stone unitig
  mapFrgStaToScfSta['R'] = 'S';  //  Stone unitig, unresolved read
  mapFrgStaToScfSta['p'] = 'S';  //  Pebble unitig
  mapFrgStaToScfSta['S'] = 'C';  //  Singleton
  mapFrgStaToScfSta['o'] = 0;    //  Other unitig


  for (uint32 fi=0; fi<frgMate.size(); fi++) {
    uint32 mi = frgMate[fi];

    assert(frgDat[fi].typ == 'f');
    assert(frgDat[mi].typ == 'f');

    //  Skip this mate if fi is not in a scaffold that we care about.

    uint32  fisi = frgDat[fi].scf.con;
    uint32  misi = frgDat[mi].scf.con;

    if (fisi == misi)
      //  Same scaffold (or both no scaffold), won't close a gap.
      continue;

    if (fisi == UINT32_MAX)
      //  fi not in a scaffold.
      continue;

    if (scfDat[fisi].typ == 'd')
      //  Degenerate.
      continue;

    assert(fisi < scfDat.size());
    assert(scfDat[fisi].sta == 'p');

    //  What type of pair is this?  If the same scaffold, skip it, since it can't fill a gap.

    if (frgDat[fi].scf.con == frgDat[mi].scf.con)
      continue;

    char  fista = mapFrgStaToScfSta[frgDat[fi].sta];
    char  mista = mapFrgStaToScfSta[frgDat[mi].sta];

    //  For now, every read is correct in the scaffold.  Later we might want to isolate out
    //  rocks, stones, etc.
    fista = 's';

    //if (fista != 's')
    //  fprintf(stderr, "fista fi "F_U32" sta %c -> fista %c\n", fi, frgDat[fi].sta, fista);

    if ((fista == 0) || (mista == 0))
      continue;

    if      ((fista == 's') && (mista == 's'))
      scfStat = scfStatDiffScaffold;

    else if ((fista == 's') && (mista == 'C'))
      scfStat = scfStatDiffScafSing;

    else if ((fista == 's') && (mista == 'S'))
      scfStat = scfStatDiffScafSurr;
    
    else if ((fista == 's') && (mista == 'd'))
      scfStat = scfStatDiffScafDegen;

    else
      assert(0);


    scfStat->number++;
    scfStat->fill.addGapFill(fi);
  }


  scfStatDiffScaffold->finalize();
  scfStatDiffScafDegen->finalize();
  scfStatDiffScafSing->finalize();
  scfStatDiffScafSurr->finalize();

  //  Unlink whatever report is there; printReport() below appends to the file.

  scfStatDiffScaffold->eraseReport(outPrefix);

  scfStatDiffScaffold->printReport(outPrefix);
  scfStatDiffScafDegen->printReport(outPrefix);
  scfStatDiffScafSing->printReport(outPrefix);
  scfStatDiffScafSurr->printReport(outPrefix);

  delete scfStatDiffScaffold;
  delete scfStatDiffScafDegen;
  delete scfStatDiffScafSing;
  delete scfStatDiffScafSurr;
}
