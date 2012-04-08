
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

const char *mainid = "$Id: analyzePosMap.C,v 1.1 2012-04-08 13:36:06 brianwalenz Exp $";

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_splitToWords.H"

#include <map>
#include <vector>
#include <string>

using namespace std;


class pmEntry {
public:
  pmEntry() {
    con = UINT32_MAX;
    len = UINT32_MAX;
    bgn = UINT32_MAX;
    end = UINT32_MAX;
    ori = '?';
    typ = '?';
  };
  ~pmEntry() {
  };

  uint32   con;  //  The container parent
  uint32   len;  //  Length of object
  uint32   bgn;  //  Begin coord in parent
  uint32   end;  //  End coord in parent
  char     ori;  //  Orientation in parent
  char     typ;  //  Type of object; 's'caffold or 'd'egenerate
};



map<string,uint32>   libNames;   //  Map a name to mean/stddev
vector<uint32>       libMean;
vector<uint32>       libStdDev;

vector<uint32>       frgMate;
vector<uint32>       frgLibrary;



map<string,uint32>   scfNames;   //  Map a name to an index into vector<pmEntry>
map<string,uint32>   ctgNames;
map<string,uint32>   utgNames;
map<string,uint32>   frgNames;   //  (loaded from posmap.frags)

vector<pmEntry>      scfDat;
vector<pmEntry>      ctgDat;
vector<pmEntry>      utgDat;
vector<pmEntry>      frgDat;

vector<string>       scfNam;
vector<string>       ctgNam;
vector<string>       utgNam;
vector<string>       frgNam;



FILE *
openFile(char *prefix, char *type) {
  char  fileName[FILENAME_MAX];
  FILE *file;

  sprintf(fileName, "%s.posmap.%s", prefix, type);

  //fprintf(stderr, "open '%s'\n", fileName);

  errno = 0;
  file = fopen(fileName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", fileName, strerror(errno)), exit(1);

  return(file);
}



void
readPosMapLen(map<string,uint32> &names,
              vector<pmEntry>    &dat,
              vector<string>     &nam,
              char               *prefix,
              char               *posmap,
              char                scfdeg) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  pmEntry pme;

  fgets(line, lineMax, file);
  while (!feof(file)) {
    splitToWords W(line);

    string   name(W[0]);
    uint32   id = dat.size();

    assert(names.count(name) == 0);

    names[name] = id;

    pme.con = UINT32_MAX;
    pme.bgn = 0;
    pme.end = 0;
    pme.len = W(1);
    pme.ori = 0;
    pme.typ = scfdeg;

    dat.push_back(pme);
    nam.push_back(name);

    fgets(line, lineMax, file);
  }

  fclose(file);
  delete [] line;
}



void
readPosMapToScf(map<string,uint32> &scfNames,
                map<string,uint32> &names,
                vector<pmEntry>    &dat,
                char               *prefix,
                char               *posmap,
                char                degscf) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  fgets(line, lineMax, file);
 
  while (!feof(file)) {
    splitToWords W(line);

    string   objName(W[0]);
    string   scfName(W[1]);

    assert(names.count(objName) > 0);
    assert(scfNames.count(scfName) > 0);

    uint32   objId = names[objName];
    uint32   scfId = scfNames[scfName];

    dat[objId].con = scfId;
    dat[objId].bgn = W(2);
    dat[objId].end = W(3);
    dat[objId].len = W(3) - W(2);
    dat[objId].ori = W[4][0];
    dat[objId].typ = degscf;

    fgets(line, lineMax, file);
  }

  fclose(file);
  delete [] line;
}


//  After everything is loaded, read the frgutg posmap to determine which
//  fragment are unplaced in surrogates.
void
readPosMapToUtg(map<string,uint32> &scfNames,
                map<string,uint32> &names,
                vector<pmEntry>    &dat,
                char               *prefix,
                char               *posmap,
                char                degscf) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  fgets(line, lineMax, file);
 
  while (!feof(file)) {
    splitToWords W(line);

    string   objName(W[0]);

    assert(names.count(objName) > 0);

    uint32   objId = names[objName];

    if (dat[objId].typ == '?') {
      dat[objId].con = 0;  //  Could be a unitig id
      dat[objId].bgn = 0;
      dat[objId].end = 0;
      dat[objId].len = 0;
      dat[objId].ori = 0;
      dat[objId].typ = 'S';
    }

    fgets(line, lineMax, file);
  }

  fclose(file);
  delete [] line;
}








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

  vector<double>   distances;
  double           mean;
  double           stddev;
};


class gapFill {
public:
  gapFill() {
  };
  ~gapFill() {
  };

  //  Compute the expected placement of the mate, compute probability that it
  //  will fall into a gap.

  void   addGapFill(uint32 fi) {
  };

  void   finalize(void) {
  };

  vector<double>    probability;
};



class scaffoldStatistics {
public:
  scaffoldStatistics(char *label_) {
    strcpy(label, label_);

    number = 0;
  };
  ~scaffoldStatistics();

  void         finalize(void) {
    innie.finalize();
    outtie.finalize();
    same.finalize();
    fill.finalize();
  };

  void         printReport(void) {
    fprintf(stdout, "----------------------------------------\n");
    fprintf(stdout, "%s\n", label);
    fprintf(stdout, "total    "F_U32" mate pairs involved\n", number);
    fprintf(stdout, "innie    "F_SIZE_T" %.2f +- %.2f bp\n", innie.distances.size(),  innie.mean,  innie.stddev);
    fprintf(stdout, "outtie   "F_SIZE_T" %.2f +- %.2f bp\n", outtie.distances.size(), outtie.mean, outtie.stddev);
    fprintf(stdout, "same     "F_SIZE_T" %.2f +- %.2f bp\n", same.distances.size(),   same.mean,   same.stddev);
  };

  char         label[1024];

  uint32       number;

  insertSize   innie;
  insertSize   outtie;
  insertSize   same;

  gapFill      fill;
};




void
analyzeGapFillProbability() {

  scaffoldStatistics   *scfStatSameScaffold    = new scaffoldStatistics("sameScaffold");
  scaffoldStatistics   *scfStatSameDegenerate  = new scaffoldStatistics("sameDegenerate");
  scaffoldStatistics   *scfStatDiffScaffold    = new scaffoldStatistics("diffScaffold");
  scaffoldStatistics   *scfStatDiffDegenerate  = new scaffoldStatistics("diffDegenerate");
  scaffoldStatistics   *scfStatDiffScafDegen   = new scaffoldStatistics("diffScafDegen");

  scaffoldStatistics   *scfStatDiffScafSing    = new scaffoldStatistics("diffScafSin");
  scaffoldStatistics   *scfStatDiffDegenSing   = new scaffoldStatistics("diffDegenSing");

  scaffoldStatistics   *scfStat                = NULL;

  for (uint32 fi=0; fi<frgMate.size(); fi++) {
    uint32 mi = frgMate[fi];

    if ((frgDat[fi].typ == 'D') ||
        (frgDat[fi].typ == 'C') ||
        (frgDat[fi].typ == 'S'))
      //  Not in a scaffold, fragment deleted or singleton or unplaced surrogate.
      continue;

    if ((mi <= fi) && (frgDat[mi].typ != 'C'))
      //  Catches both unmated (mi==0) and visited (mi < fi), in particular,
      //  unmated fi=0 is caught.  The test on con is to not skip mates where the
      //  lower id is a singleton.
      continue;

    //  What type of pair is this?

    if (frgDat[fi].con == frgDat[mi].con) {
      //  Same Scaffold

      if (frgDat[fi].typ == 's')
        scfStat = scfStatSameScaffold;    //  Real scaffolds!
      else
        scfStat = scfStatSameDegenerate;  //  Degenerates.

      scfStat->number++;

      int32  dist = MAX(frgDat[fi].end, frgDat[mi].end) - MIN(frgDat[fi].bgn, frgDat[mi].bgn);

      if (dist <= 0) {
        fprintf(stderr, "fi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c typ=%c %s\n",
                fi,
                frgDat[fi].con, frgDat[fi].bgn, frgDat[fi].end, frgDat[fi].len, frgDat[fi].ori, frgDat[fi].typ, frgNam[fi].c_str());
        fprintf(stderr, "mi="F_U32" con="F_U32" bgn="F_U32" end="F_U32" len="F_U32" ori=%c typ=%c %s\n",
                mi,
                frgDat[mi].con, frgDat[mi].bgn, frgDat[mi].end, frgDat[mi].len, frgDat[mi].ori, frgDat[mi].typ, frgNam[mi].c_str());
      }
      assert(dist > 0);

      if        (((frgDat[fi].ori == frgDat[mi].ori))) {
        scfStat->same.addDistance(dist);

      } else if (((frgDat[fi].bgn <  frgDat[mi].bgn) && (frgDat[fi].ori == 'f')) ||
                 ((frgDat[fi].bgn >  frgDat[mi].bgn) && (frgDat[fi].ori == 'r'))) {
        scfStat->innie.addDistance(dist);
        
      } else if (((frgDat[fi].bgn <  frgDat[mi].bgn) && (frgDat[fi].ori == 'r')) ||
                 ((frgDat[fi].bgn >  frgDat[mi].bgn) && (frgDat[fi].ori == 'f'))) {
        scfStat->outtie.addDistance(dist);

      } else {
        assert(0);
      }

    } else {
      //  Different Scaffolds
      if      ((frgDat[fi].typ == 's') && (frgDat[mi].typ == 's'))
        scfStat = scfStatDiffScaffold;

      else if ((frgDat[fi].typ == 'd') && (frgDat[mi].typ == 'd'))
        scfStat = scfStatDiffDegenerate;

      else if ((frgDat[fi].typ == 's') && (frgDat[mi].typ == 'c'))
        scfStat = scfStatDiffScafSing;

      else if ((frgDat[fi].typ == 'd') && (frgDat[mi].typ == 'c'))
        scfStat = scfStatDiffDegenSing;

      else
        scfStat = scfStatDiffScafDegen;

      scfStat->number++;

      scfStat->fill.addGapFill(fi);
    }
  }


  scfStatSameScaffold->finalize();
  scfStatSameDegenerate->finalize();
  scfStatDiffScaffold->finalize();
  scfStatDiffDegenerate->finalize();
  scfStatDiffScafDegen->finalize();
  scfStatDiffScafSing->finalize();
  scfStatDiffDegenSing->finalize();

  scfStatSameScaffold->printReport();
  scfStatSameDegenerate->printReport();
  scfStatDiffScaffold->printReport();
  scfStatDiffDegenerate->printReport();
  scfStatDiffScafDegen->printReport();
  scfStatDiffScafSing->printReport();
  scfStatDiffDegenSing->printReport();
}





int
main(int argc, char **argv) {
  char    *prefix   = NULL;
  char    *gkpName  = NULL;

  uint32   lineMax  = 1048576;
  char    *line     = new char [1048576];

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if (prefix == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -o prefix [-h] [-i prefix.asm | < prefix.asm]\n", argv[0]);
    fprintf(stderr, "  -p prefix        prefix of posmap files; prefix of output files\n");

    fprintf(stderr, "\n");

    if  (prefix == NULL)
      fprintf(stderr, "ERROR: No assembly prefix (-p) supplied.\n");

    exit(1);
  }


  //  Read library information.

  fprintf(stderr, "Reading library information.\n");
  {
    FILE *file = openFile(prefix, "libraries");

    fgets(line, lineMax, file);
    while (!feof(file)) {
      splitToWords W(line);

      string   libName(W[0]);
      uint32   libId = libNames.size();

      libNames[libName] = libId;

      libMean.push_back(W(1));
      libMean.push_back(W(2));

      fgets(line, lineMax, file);
    }
  }

  //  Read fragment information.

  fprintf(stderr, "Reading fragment information.\n");
  {
    FILE *file = openFile(prefix, "frags");

    fgets(line, lineMax, file);
    while (!feof(file)) {
      splitToWords W(line);

      string   frgName(W[0]);
      string   libName(W[5]);

      uint32   frgId = frgNames.size();
      uint32   libId = libNames[libName];

      frgNames[frgName] = frgId;

      pmEntry pme;

      pme.con = UINT32_MAX;
      pme.bgn = 0;
      pme.end = 0;
      pme.len = 0;
      pme.ori = 0;
      pme.typ = 0;

      if      (W[3][0] == 'p')
        pme.typ = '?';
      else if (W[3][0] == 'd')
        pme.typ = 'D';  //  Deleted
      else if (W[3][0] == 'c')
        pme.typ = 'C';  //  Singleton (chaff)
      else
        assert(0);

      assert(frgDat.size() == frgId);
      assert(frgMate.size() == frgId);
      assert(frgNam.size() == frgId);

      frgDat.push_back(pme);
      frgMate.push_back(0);  //  All are initially unmated.
      frgNam.push_back(frgName);

      frgLibrary.push_back(libId);

      fgets(line, lineMax, file);
    }
  }

  //  Read mate information.

  fprintf(stderr, "Reading mate information.\n");
  {
    FILE *file = openFile(prefix, "mates");

    fgets(line, lineMax, file);
    while (!feof(file)) {
      splitToWords W(line);

      string   fr1Name(W[0]);
      string   fr2Name(W[1]);

      assert(frgNames.count(fr1Name) > 0);
      assert(frgNames.count(fr2Name) > 0);

      uint32   fr1Id = frgNames[fr1Name];
      uint32   fr2Id = frgNames[fr2Name];

      frgMate[fr1Id] = fr2Id;  //  Now they're mated!
      frgMate[fr2Id] = fr1Id;

      fgets(line, lineMax, file);
    }
  }

  //  Read scaffold information.

  fprintf(stderr, "Reading assembly information (scf, ctg, utg).\n");
  {
    readPosMapLen(scfNames, scfDat, scfNam, prefix, "scflen", 's');  //  Real scaffolds
    readPosMapLen(scfNames, scfDat, scfNam, prefix, "deglen", 'd');  //  Dreg scaffolds (promoted from dreg contigs)

    readPosMapLen(ctgNames, ctgDat, ctgNam, prefix, "ctglen", 's');  //  Real contigs
    readPosMapLen(ctgNames, ctgDat, ctgNam, prefix, "deglen", 'd');  //  Dreg contigs

    readPosMapLen(utgNames, utgDat, utgNam, prefix, "utglen", 'x');  //  Real and Dreg unitigs
  }

  fprintf(stderr, "Reading contig/unitig/fragment -> scaffold information.\n");
  {
    ctgDat.resize(ctgNames.size());
    utgDat.resize(utgNames.size());
    //frgDat.resize(frgNames.size());

    readPosMapToScf(scfNames, ctgNames, ctgDat, prefix, "ctgscf", 's');
    //  There is no degscf, we have to make it up, below.

    readPosMapToScf(scfNames, utgNames, utgDat, prefix, "utgscf", 's');
    readPosMapToScf(scfNames, utgNames, utgDat, prefix, "utgdeg", 'd');

    readPosMapToScf(scfNames, frgNames, frgDat, prefix, "frgscf", 's');
    readPosMapToScf(scfNames, frgNames, frgDat, prefix, "frgdeg", 'd');

    readPosMapToUtg(scfNames, frgNames, frgDat, prefix, "frgutg", 'S');  //  Surrogates
  }

  {
    for (uint32 si=0; si<scfDat.size(); si++) {
      if (scfDat[si].typ == 's')
        continue;

      uint32  ctgId = ctgDat.size();

      ctgDat[ctgId].con = si;
      ctgDat[ctgId].bgn = 0;
      ctgDat[ctgId].end = scfDat[si].len;
      ctgDat[ctgId].len = scfDat[si].len;
      ctgDat[ctgId].ori = 'f';
      ctgDat[ctgId].typ = 'd';
    }
  }


  fprintf(stderr, "Loaded "F_SIZE_T" scaffolds.\n",   scfNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" contigs.\n",     ctgNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" unitigs.\n",     utgNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" fragments.\n",   frgNames.size());

  //  Do something.



  //  From some single scaffold, count the number of mates that are external to it

#if 0
  for (uint32 si=0; si < scfNames.size(); si++) {
    uint32  internal = 0;
    uint32  external = 0;

    for (uint32 fi=0; fi < frgNames.size(); fi++) {
      if (frgDat[fi].con != si)
        continue;

      uint32  mi = frgMate[fi];

      if (mi == 0)
        continue;
      
      if (frgDat[fi].con == frgDat[mi].con)
        internal++;
      else
        external++;
    }

    fprintf(stderr, "scaffold "F_U32" %s %c internal "F_U32" external "F_U32"\n",
            si,
            scfNam[si].c_str(),
            scfDat[si].typ,
            internal,
            external);
  }
#endif



  //  For each mate pair, compute stats.  per library?
  //
  //  same scaffold (normal/degenerate)
  //  same contig (normal/degenerate)
  //  different scaffold (normal-normal, etc)
  //  different contig (normal-normal, etc -- special case for same scaffold)
  //
  //  gap fill probability
  //   - for each external mate, compute the probability the frag belongs in a gap, plot histogram


  analyzeGapFillProbability();

  exit(0);
}
