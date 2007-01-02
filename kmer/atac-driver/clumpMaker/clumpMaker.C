#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "atac.H"

//  Aaron Halpern's clumpMaker algorithm.
//
//  To reproduce the original clumpMaker exactly, assuming that your
//  atac mapping is for a QUERYvsREFERENCE:
//
//  in=VISD6vsB35LC/VISD6vsB35LC
//  cut -d' ' -f 1-12 < $in.atac.ckpLast | grep "^M u" | sort -k5,5 -k6n > tmp.a.clumpMaker
//  $clumpMaker -S -c 50000 -2 -f tmp.a.clumpMaker > $in.50000.clumps
//
//  That is, use only the first 12 columns of info, only ungapped
//  matches (the original also allows gapped matches, but we don't
//  have any of those), then sort by the QUERY iid and position.  Yes,
//  the sort is supposed to be alphanumeric.
//
//  Then, run clumpMaker DISABLING it's sort (which sorts iids
//  numerically), using the second sequence as the reference.
//

class tClumpHit {
public:

  void set(atacMatch *m, bool seq1IsRef) {

    match    = *m;
    matchIID =  m->matchiid;

    if (seq1IsRef) {
      refIID       = m->iid1;
      refBeg       = m->pos1;
      refEnd       = m->pos1 + m->len1;
      qryIID       = m->iid2;
      qryBeg       = m->pos2;
      qryEnd       = m->pos2 + m->len2;
    } else {
      refIID       = m->iid2;
      refBeg       = m->pos2;
      refEnd       = m->pos2 + m->len2;
      qryIID       = m->iid1;
      qryBeg       = m->pos1;
      qryEnd       = m->pos1 + m->len1;
    }

    ori          = m->fwd2 ? 1 : -1;

    bestStart    = -1;
    bestExtend   = -1;
    scoreStart   = 0;
    scoreExtend  = 0;
    clump        = -1;
  };

  s64bit   get_bestScore() const {
    return(max(scoreStart, scoreExtend));
  };

  atacMatch  match;

  u32bit     matchIID;

  u32bit     refIID;
  s32bit     refBeg;
  s32bit     refEnd;

  u32bit     qryIID;
  s32bit     qryBeg;
  s32bit     qryEnd;

  s32bit     ori;

  s32bit     scoreStart;
  s32bit     bestStart;
  s32bit     scoreExtend;
  s32bit     bestExtend;
  s32bit     clump;
};



int
clumpHitCompareQry(const void *A, const void *B) {
  const tClumpHit *a = (const tClumpHit *)A;
  const tClumpHit *b = (const tClumpHit *)B;

  if (a->qryIID  > b->qryIID)   return(1);
  if (a->qryIID  < b->qryIID)   return(-1);
  if (a->qryBeg  > b->qryBeg)   return(1);
  if (a->qryBeg  < b->qryBeg)   return(-1);
  if (a->qryEnd  > b->qryEnd)   return(1);
  if (a->qryEnd  < b->qryEnd)   return(-1);

  if (a->refIID  > b->refIID)   return(1);
  if (a->refIID  < b->refIID)   return(-1);
  if (a->refBeg  > b->refBeg)   return(1);
  if (a->refBeg  < b->refBeg)   return(-1);
  if (a->refEnd  > b->refEnd)   return(1);
  if (a->refEnd  < b->refEnd)   return(-1);

  return(0);
}

int
clumpHitCompareIID(const void *A, const void *B) {
  const tClumpHit *a = (const tClumpHit *)A;
  const tClumpHit *b = (const tClumpHit *)B;

  if (a->matchIID > b->matchIID)    return(1);
  if (a->matchIID < b->matchIID)    return(-1);
  return(0);
}


bool
chainable(tClumpHit *a, tClumpHit *b, s32bit maxjump) {

  // return false if
  //    hits are to different chromosomes
  //    hits are not similarly oriented
  //    hits are too far apart on query axis
  //    hits are too far apart on reference axis
  //    hits are "out of order" (we're sorted by the qry)
  //
  return(!((a->refIID != b->refIID) || (a->qryIID != b->qryIID) ||
           (a->ori != b->ori) ||
           (b->qryBeg - a->qryEnd > maxjump) || 
           (a->ori * (b->refBeg - a->refEnd) > maxjump) ||
           (a->ori * (b->refBeg - a->refBeg) <  0)));
}


s32bit
score_all_hits(tClumpHit *hits,
               s32bit     clumpcost,
               s32bit     maxjump,
               u32bit     num_hits){

  // location of best score so far (to which we point whenever starting a new clump)
  s32bit bestEnd = -1;


  // best scores so far internal to a reference unit (scaffold, chromosome, etc)
  s32bit bestEndThis   = -1;
  s32bit bestScoreThis = -clumpcost;

  // furthest back still accessible ...
  u32bit furthest_back=0;
  
  for(u32bit i=0; i<num_hits; i++) {

    if ((i==0) || (hits[i].qryIID != hits[i-1].qryIID)) {
      bestEnd       = bestEndThis; // best of previous query unit
      bestEndThis   = i-1;
      bestScoreThis = -clumpcost;
    }

    // find best way of using this as start of a new clump
    if ((bestEndThis >= 0) &&
        (bestScoreThis >= 0)) {
      // start new clump that is not the first for this reference unit
      hits[i].scoreStart = hits[i].qryEnd - hits[i].qryBeg + bestScoreThis - clumpcost;
      hits[i].bestStart  = bestEndThis;
    } else {
      // clump would be first (to be used) for this reference unit
      hits[i].scoreStart = hits[i].qryEnd - hits[i].qryBeg - clumpcost;
      hits[i].bestStart  = bestEnd;
    }

    // find best way of extending a clump, if any
    if (furthest_back < i) {
      s32bit cutoff = hits[i].qryBeg - maxjump;

      while ((hits[furthest_back].qryIID != hits[i].qryIID) ||
             (hits[furthest_back].qryEnd < cutoff))
	furthest_back++;
    }

    s32bit extendScore = -clumpcost;
    s32bit extendprev  = -1;
    for (u32bit j=furthest_back; j<i; j++) {
      if (chainable(hits+j, hits+i, maxjump)) {
	s32bit tmpscore = hits[j].get_bestScore() + hits[i].qryEnd - hits[i].qryBeg;
	if(extendScore < tmpscore){
	  extendScore=tmpscore;
	  extendprev=j;
	}
      }
    }
    hits[i].scoreExtend = extendScore;
    hits[i].bestExtend  = extendprev;

    //figure out whether this is a new best ...
    s32bit tmpscore = hits[i].get_bestScore();
    if (tmpscore > bestScoreThis) {
      bestScoreThis = tmpscore;
      bestEndThis   = i;
    }
  }

  return(bestEndThis);
}








int
main(int argc, char **argv) {
  s32bit   clumpcost    = 50000;
  s32bit   maxjump      = 200000;
  bool     seq1IsRef    = false;
  char    *filename     = 0L;
  bool     isSorted     = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      clumpcost = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-j") == 0) {
      maxjump = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-1") == 0) {
      seq1IsRef = true;
    } else if (strcmp(argv[arg], "-2") == 0) {
      seq1IsRef = false;
    } else if (strcmp(argv[arg], "-f") == 0) {
      filename = argv[++arg];
    } else if (strcmp(argv[arg], "-S") == 0) {
      isSorted = true;
    } else {
      fprintf(stderr, "Unknown argument '%s'\n", argv[arg]);
    }
    arg++;
  }

  if (filename == 0L) {
    fprintf(stderr, "usage: %s [] -f filename\n", argv[0]);
    fprintf(stderr, "  -c x    penalty for clump start, default 50000\n");
    fprintf(stderr, "  -j x    max jump between consistent hits in a clump, default 200000\n");
    fprintf(stderr, "  -1      the reference assembly is the first one.\n");
    fprintf(stderr, "  -2      the reference assembly is the second one (default).\n");
    fprintf(stderr, "  -S      assume the input is already sorted by the query IID, position.\n");
    fprintf(stderr, "          this will also make the output sorted by queryIID, queryPosition\n");
    exit(1);
  }




  fprintf(stderr, "1 load the matches\n");

  atacFile    *AF      = new atacFile(filename);
  u32bit       hitsLen = AF->matches()->numberOfMatches();
  tClumpHit   *hits    = new tClumpHit [hitsLen];
  
  for (u32bit i=0; i<hitsLen; i++)
    hits[i].set(AF->matches()->getMatch(i), seq1IsRef);



  fprintf(stderr, "2 sort the matches\n");
  qsort(hits, hitsLen, sizeof(tClumpHit), clumpHitCompareQry);



  fprintf(stderr, "3 score the matches\n");
  s32bit bestEnd = score_all_hits(hits,
                                  clumpcost,
                                  maxjump,
                                  hitsLen);

  //  Mark the clumps
  //
  fprintf(stderr, "4 mark clumps\n");
  u32bit clump = 0;

  while(bestEnd >= 0) {
    hits[bestEnd].clump = clump;

    if (hits[bestEnd].scoreExtend > hits[bestEnd].scoreStart) {
      bestEnd = hits[bestEnd].bestExtend;
    } else {
      bestEnd = hits[bestEnd].bestStart;
      clump++;
    }
  }

  //  Sort the hits by iid, then merge into the output
  //
  fprintf(stderr, "5 sort the matches\n");
  qsort(hits, hitsLen, sizeof(tClumpHit), clumpHitCompareIID);


  //  For each clump, find the min/max extent in both sequences.  We
  //  use this to output the clump match record.
  //
  s32bit *clumpLoRef = new s32bit [clump];
  s32bit *clumpHiRef = new s32bit [clump];
  s32bit *clumpLoQry = new s32bit [clump];
  s32bit *clumpHiQry = new s32bit [clump];
  bool   *clumpOut   = new bool   [clump];

  for (u32bit xx=0; xx<clump; xx++) {
    clumpLoRef[xx] = 1000000000;
    clumpHiRef[xx] = 0;
    clumpLoQry[xx] = 1000000000;
    clumpHiQry[xx] = 0;
    clumpOut[xx] = false;
  }

  for (u32bit xx=0; xx<hitsLen; xx++) {
    s32bit  cc = hits[xx].clump;
    if (cc >= 0) {
      if (hits[xx].refBeg < clumpLoRef[cc])   clumpLoRef[cc] = hits[xx].refBeg;
      if (hits[xx].refEnd > clumpHiRef[cc])   clumpHiRef[cc] = hits[xx].refEnd;

      if (hits[xx].qryBeg < clumpLoQry[cc])   clumpLoQry[cc] = hits[xx].qryBeg;
      if (hits[xx].qryEnd > clumpHiQry[cc])   clumpHiQry[cc] = hits[xx].qryEnd;
    }
  }

  //  Dump the clumps
  //

  fprintf(stderr, "6 output matches with clumps\n");

  AF->writeHeader(stdout);

  for (u32bit mm=0; mm<hitsLen; mm++) {
    s32bit cc = hits[mm].clump;

    if ((cc >= 0) &&
        (clumpOut[cc] == false)) {
      atacMatch  C;
      sprintf(C.matchuid,  "clump"s32bitFMTW(06), cc);
      sprintf(C.parentuid, ".");
      C.matchiid = 0;
      C.type[0] = 'c';
      C.type[1] = 0;

      C.iid1 = hits[mm].match.iid1;
      C.iid2 = hits[mm].match.iid2;

      //  Set the position and length based on the correct reference
      //  -- in particular, since we get the IID and orientation from
      //  the copy of the match we don't need to listen to the
      //  seq1IsRef flag for those.

      if (seq1IsRef) {
        C.pos1 = clumpLoRef[cc];
        C.len1 = clumpHiRef[cc] - clumpLoRef[cc];
        C.pos2 = clumpLoQry[cc];
        C.len2 = clumpHiQry[cc] - clumpLoQry[cc];
      } else {
        C.pos1 = clumpLoQry[cc];
        C.len1 = clumpHiQry[cc] - clumpLoQry[cc];
        C.pos2 = clumpLoRef[cc];
        C.len2 = clumpHiRef[cc] - clumpLoRef[cc];
      }

      C.fwd1 = hits[mm].match.fwd1;
      C.fwd2 = hits[mm].match.fwd2;

      C.print(stdout, AF->labelA(), AF->labelB());

      clumpOut[cc] = true;
    }


    if (cc >= 0)
      sprintf(hits[mm].match.parentuid, "clump"s32bitFMTW(06), cc);
    else
      sprintf(hits[mm].match.parentuid, ".");

    hits[mm].match.print(stdout, AF->labelA(), AF->labelB());
  }

  return(0);
}
