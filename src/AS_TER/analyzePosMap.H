
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

#ifndef ANALYZEPOSMAP_H
#define ANALYZEPOSMAP_H

static const char *rcsid_ANALYZEPOSMAP_H = "$Id: analyzePosMap.H,v 1.1 2012-12-05 01:13:23 brianwalenz Exp $";

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_splitToWords.H"
#include "AS_UTL_intervalList.H"

#include <map>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;



class libEntry {
public:
  libEntry() {
    memset(name, 0, sizeof(char) * 256);
    mean   = 0;
    stddev = 0;
  };
  ~libEntry() {
  };

  void  initialize(char *name_, int32 mean_, int32 stddev_) {
    strcpy(name, name_);
    mean   = mean_;
    stddev = stddev_;

    pdf.resize(6 * stddev + 1);

    double  c = 1 / (stddev * sqrt(2 * M_PI));
    double  d = 2 * stddev * stddev;

    for (int32 i=0; i<=6 * stddev; i++) {
      double x = i - 3 * stddev;

      pdf[i] = c * exp(-(x * x / d));
    }

    //FILE *F = fopen("libdat", "w");
    //for (uint32 i=0; i<=6*stddev; i++)
    //  fprintf(F, "%d\t%f\n", i, pdf[i]);
    //fclose(F);
  };

  char            name[256];

  int32           mean;
  int32           stddev;
  vector<double>  pdf;
};


class mapEntry {
public:
  mapEntry(uint32 con_, int32 bgn_, int32 end_, char ori_) {
    con = con_;
    bgn = bgn_;
    end = end_;
    ori = ori_;
  };
  mapEntry() {
    con = UINT32_MAX;
    bgn = 0;
    end = UINT32_MAX;
    ori = 'u';
  };
  ~mapEntry() {};

  uint32   con;  //  The container parent
  int32    bgn;  //  Begin coord in parent
  int32    end;  //  End coord in parent
  char     ori;  //  Orientation in parent
};


class frgEntry {
public:
  frgEntry() {
    len = UINT32_MAX;
    typ = '?';
    sta = '?';
  };
  ~frgEntry() {};

  uint32           len;  //  Length of object

  mapEntry         utg;
  mapEntry         ctg;
  mapEntry         scf;

  char             typ;  //  Type of object; fragment, unitig, degenerate, contig, scaffold
  char             sta;  //  Status of object; see below
};


class utgEntry {
public:
  utgEntry() {
    len = UINT32_MAX;
    typ = '?';
    sta = '?';
  };
  ~utgEntry() {};

  uint32           len;  //  Length of object

  vector<uint32>   frg;
  mapEntry         ctg;
  mapEntry         scf;

  char             typ;  //  Type of object; fragment, unitig, degenerate, contig, scaffold
  char             sta;  //  Status of object; see below
};


class ctgEntry {
public:
  ctgEntry() {
    len = UINT32_MAX;
    typ = '?';
    sta = '?';
  };
  ~ctgEntry() {};

  uint32           len;  //  Length of object

  vector<uint32>   frg;
  vector<uint32>   utg;
  mapEntry         scf;

  char             typ;  //  Type of object; fragment, unitig, degenerate, contig, scaffold
  char             sta;  //  Status of object; see below
};


class scfEntry {
public:
  scfEntry() {
    len     = 0;
    typ     = '?';
    sta     = '?';
    base    = NULL;
  };
  ~scfEntry() {
    delete [] base;
  }

  uint32           len;     //  Length of object

  vector<uint32>   frg;
  vector<uint32>   utg;
  vector<uint32>   ctg;

  char             typ;     //  Type of object; 's'caffold or 'd'egenerate
  char             sta;

  char            *base;    //  t/f if this is a gap or base
};


//  Objects
//
//    s  scaffold
//    c  contig
//    u  unitig
//    f  fragment
//
//  Status
//
//    fragment - initial label comes from 'posmap.frags'
//
//      The initial label comes from 'posmap.frags':
//        ?   from posmap 'placed'  (first letter only, not strcmp)
//        D   from posmap 'deleted' (first letter only)
//        C   from posmap 'chaff'   (first letter only)
//
//      The '?' category is refined by reading frgscf and frgdeg.  This changes ? into either
//        p   placed in a scaffold
//        d   placed in a degenerate
//
//      The 'p' category is further refined by labelFragments().  This changes 'p' into either the
//      unitig status or 'R' for unresolved repeat read.
//
//      The final list of fragment stauts codes:
//        D   deleted from assembly
//        C   singleton chaff
//        d   in degenerate contig
//        u   in unique unitig
//        r   in rock unitig
//        s   in stone unitig
//        R   in stone unitig, unresolved
//        p   in pebble (obsolete)
//        S   in singleton unitig
//        o   in other unitig
//
//    unitigs - from utgctg or utgscf
//      u   unique
//      r   rock
//      s   stone
//      p   pebble (obsolete)
//      S   singleton
//      o   other
//
//    contigs/scaffolds
//      p   in a scaffold ('placed')
//      d   in a degenerate


#define NUM_FRG_LABELS      10

extern char                 frgLabel[NUM_FRG_LABELS];
extern uint32               frgLabelToInt[256];

extern map<string,uint32>   libNames;    //  Index into libDat
extern vector<libEntry>     libDat;

extern vector<uint32>       frgMate;     //  Mapping of id -> mate id
extern vector<uint32>       frgLibrary;  //  Mapping of id -> library id

extern map<string,uint32>   frgNames;    //  Index into *Dat and *Nam
extern map<string,uint32>   utgNames;
extern map<string,uint32>   ctgNames;
extern map<string,uint32>   scfNames;

extern vector<frgEntry>     frgDat;
extern vector<utgEntry>     utgDat;
extern vector<ctgEntry>     ctgDat;
extern vector<scfEntry>     scfDat;

extern vector<string>       frgNam;
extern vector<string>       utgNam;
extern vector<string>       ctgNam;
extern vector<string>       scfNam;



void   loadPosMap(char *prefix, char *gkpName);

void   analyzeGapFillProbability(char *outPrefix);
void   analyzeLibraryFate(char *outPrefix);

void   analyzeSurrogates(char *outPrefix);
void   analyzeDegenerates(char *outPrefix);

#endif  //  ANALYZEPOSMAP_H
