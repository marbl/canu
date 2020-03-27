
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef INCLUDE_BOGUSUTIL
#define INCLUDE_BOGUSUTIL

#include "runtime.H"

#include "files.H"

#include "strings.H"
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
    chnBgn = chnEnd = 0;
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
