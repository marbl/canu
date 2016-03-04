
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_bogusUtil.H
 *    src/bogart/AS_BAT_bogusUtil.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2014-DEC-23
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_BOGUSUTIL
#define INCLUDE_BOGUSUTIL

#include "AS_global.H"

#include "splitToWords.H"
#include "intervalList.H"

#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;


class genomeAlignment {
public:
  genomeAlignment() {
    frgIID = frgBgn = frgEnd = 0;
    genIID = 0;
    genBgn = genEnd = 0;
    identity = 0.0;
    isDeleted = isReverse = isSpanned = isRepeat = false;
  };

  int32   frgIID;
  int32   frgBgn;
  int32   frgEnd;

  int32   genIID;      //  Position in the actual sequence
  int32   genBgn;
  int32   genEnd;

  int32   chnBgn;      //  Position in the chained sequences
  int32   chnEnd;

  double  identity;    //  Percent identity of the alignment

  bool    isDeleted;   //  Used by bogusness

  bool    isReverse;
  bool    isSpanned;
  bool    isRepeat;
};


class referenceSequence {
public:
  referenceSequence(int32 cb, int32 ce, int32 rl, char *cn) {
    rschnBgn = cb;
    rschnEnd = ce;
    rsrefLen = rl;
    assert(strlen(cn) < 256);
    strcpy(rsrefName, cn);
  };

  int32  rschnBgn;
  int32  rschnEnd;
  int32  rsrefLen;
  char   rsrefName[256];
};


bool byFragmentID(const genomeAlignment &A, const genomeAlignment &B);
bool byGenomePosition(const genomeAlignment &A, const genomeAlignment &B);


void addAlignment(vector<genomeAlignment>   &genome,
                  int32  frgIID,
                  int32  frgBgn, int32  frgEnd, bool  isReverse,
                  int32  chnBgn, int32  chnEnd,
                  double identity,
                  int32  genIID,
                  int32  genBgn, int32  genEnd);

void loadNucmer(char                       *nucmerName,
                vector<genomeAlignment>    &genome,
                map<string, int32>         &IIDmap,
                vector<string>             &IIDname,
                vector<referenceSequence>  &refList,
                map<string,uint32>         &refMap,
                double                      minIdentity);

void loadSnapper(char                       *snapperName,
                 vector<genomeAlignment>    &genome,
                 map<string, int32>         &IIDmap,
                 vector<string>             &IIDname,
                 vector<referenceSequence>  &refList,
                 map<string,uint32>         &refMap,
                 double                      minIdentity);

void
loadReferenceSequence(char                       *refName,
                      vector<referenceSequence>  &refList,
                      map<string,uint32>         &refMap);

#endif  //  INCLUDE_BOGUSUTIL
