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
  void set(u32bit matchiid,
           u32bit sid1, u32bit pos1, u32bit len1, u32bit fwd1,
           u32bit sid2, u32bit pos2, u32bit len2, u32bit fwd2, bool seq1IsRef) {
    iid          = matchiid;

    if (seq1IsRef) {
      refCh        = sid1;
      qryCh        = sid2;
      refBeg       = pos1;
      refEnd       = pos1 + len1;
      qryBeg       = pos2;
      qryEnd       = pos2 + len2;
    } else {
      refCh        = sid2;
      qryCh        = sid1;
      refBeg       = pos2;
      refEnd       = pos2 + len2;
      qryBeg       = pos1;
      qryEnd       = pos1 + len1;
    }

    ori          = fwd2 ? 1 : -1;

    bestStart    = -1;
    bestExtend   = -1;
    scoreStart   = 0;
    scoreExtend  = 0;
    clump        = -1;
  };

  s64bit   get_bestScore() const {
    return(max(scoreStart, scoreExtend));
  };

  u32bit   iid;
  s32bit   ori;
  u32bit   refCh;
  u32bit   qryCh;
  s32bit   refBeg;
  s32bit   refEnd;
  s32bit   qryBeg;
  s32bit   qryEnd;

  s32bit   scoreStart;
  s32bit   bestStart;
  s32bit   scoreExtend;
  s32bit   bestExtend;
  s32bit   clump;
};



int
clumpHitCompareQry(const void *A, const void *B) {
  const tClumpHit *a = (const tClumpHit *)A;
  const tClumpHit *b = (const tClumpHit *)B;

  if (a->qryCh  > b->qryCh)    return(1);
  if (a->qryCh  < b->qryCh)    return(-1);
  if (a->qryBeg > b->qryBeg)   return(1);
  if (a->qryBeg < b->qryBeg)   return(-1);
  if (a->qryEnd > b->qryEnd)   return(1);
  if (a->qryEnd < b->qryEnd)   return(-1);

  if (a->refCh  > b->refCh)    return(1);
  if (a->refCh  < b->refCh)    return(-1);
  if (a->refBeg > b->refBeg)   return(1);
  if (a->refBeg < b->refBeg)   return(-1);
  if (a->refEnd > b->refEnd)   return(1);
  if (a->refEnd < b->refEnd)   return(-1);

  return(0);
}

int
clumpHitCompareIID(const void *A, const void *B) {
  const tClumpHit *a = (const tClumpHit *)A;
  const tClumpHit *b = (const tClumpHit *)B;

  if (a->iid  > b->iid)    return(1);
  if (a->iid  < b->iid)    return(-1);
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
  return(!((a->refCh != b->refCh) || (a->qryCh != b->qryCh) ||
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

    if ((i==0) || (hits[i].qryCh != hits[i-1].qryCh)) {
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

      while ((hits[furthest_back].qryCh != hits[i].qryCh) ||
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





u32bit
matchesInInput(char *filename) {
  errno = 0;
  FILE *inFile = fopen(filename, "r");
  if (errno)
    fprintf(stderr, "Couldn't open '%s': %s\n", filename, strerror(errno)), exit(1);

  u32bit    iid = 0;
  char      inLine[1024];

  fgets(inLine, 1024, inFile);

  while (!feof(inFile)) {
    if ((inLine[0] == 'M') && (inLine[2] == 'u'))
      iid++;
    fgets(inLine, 1024, inFile);
  }

  fclose(inFile);

  return(iid);
}


tClumpHit*
readMatches(char *filename, u32bit hitsLen, bool seq1IsRef) {
  u32bit             iid  = 0;
  tClumpHit         *hits = new tClumpHit [hitsLen];

  atacFileStream  AF(filename);
  atacMatch      *m = AF.nextMatch('u');

  while (m) {
    hits[iid].set(iid,
                  m->iid1, m->pos1, m->len1, m->fwd1,
                  m->iid2, m->pos2, m->len2, m->fwd2,
                  seq1IsRef);
    iid++;

    m = AF.nextMatch('u');
  }

  return(hits);
}




int
main(int argc, char **argv) {
  s32bit   clumpcost    = 50000;
  s32bit   maxjump      = 200000;
  bool     refisone     = false;
  char    *filename     = 0L;
  bool     isSorted     = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-c") == 0) {
      clumpcost = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-j") == 0) {
      maxjump = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-1") == 0) {
      refisone = true;
    } else if (strcmp(argv[arg], "-2") == 0) {
      refisone = false;
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


  fprintf(stderr, "Pass 1: count the number of matches in the input.\n");
  u32bit             hitsLen = matchesInInput(filename);
  fprintf(stderr, "        found "u32bitFMT" hits\n", hitsLen);

  fprintf(stderr, "Pass 2: load the matches\n");
  tClumpHit         *hits    = readMatches(filename, hitsLen, refisone);

  if (isSorted == false) {
    fprintf(stderr, "        sort the matches\n");
    qsort(hits, hitsLen, sizeof(tClumpHit), clumpHitCompareQry);
  }

  fprintf(stderr, "        score the matches\n");
  s32bit bestEnd = score_all_hits(hits,
                                  clumpcost,
                                  maxjump,
                                  hitsLen);

  //  Mark the clumps
  //
  fprintf(stderr, "        mark clumps\n");
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
  fprintf(stderr, "        sort the matches\n");
  qsort(hits, hitsLen, sizeof(tClumpHit), clumpHitCompareIID);


  //  For each clump, find the min/max extent in both sequences.  We
  //  use this to output the clump match record.
  //
  s32bit *clumpLoA = new s32bit [clump];
  s32bit *clumpHiA = new s32bit [clump];
  s32bit *clumpLoB = new s32bit [clump];
  s32bit *clumpHiB = new s32bit [clump];
  bool   *clumpOut = new bool   [clump];
  for (u32bit xx=0; xx<clump; xx++) {
    clumpLoA[xx] = 1000000000;
    clumpHiA[xx] = 0;
    clumpLoB[xx] = 1000000000;
    clumpHiB[xx] = 0;
    clumpOut[xx] = false;
  }
  for (u32bit xx=0; xx<hitsLen; xx++) {
    if (hits[xx].clump >= 0) {
      s32bit  c = hits[xx].clump;
      if (hits[xx].refBeg < clumpLoA[c])   clumpLoA[c] = hits[xx].refBeg;
      if (hits[xx].refEnd > clumpHiA[c])   clumpHiA[c] = hits[xx].refEnd;

      if (hits[xx].qryBeg < clumpLoB[c])   clumpLoB[c] = hits[xx].qryBeg;
      if (hits[xx].qryEnd < clumpHiB[c])   clumpHiB[c] = hits[xx].qryEnd;
    }
  }

  //  Dump the clumps
  //

  fprintf(stderr, "Pass 3: output matches with clumps\n");


  atacFileStream  AF(filename);
  atacMatch      *m = AF.nextMatch('u');
  u32bit          xx = 0;

  while (m) {

    //  Make sure the iid agrees.
    if (xx != hits[xx].iid) {
      fprintf(stderr, "Augh!  Merge failure!  Not in sync!\n");
      fprintf(stderr, "xx = "u32bitFMT", hits[xx].iid = "u32bitFMT"\n", xx, hits[xx].iid);
      exit(1);
    }

    //  XXX This block was originally disabled!
    if (((hits[xx].refCh != m->iid1) && (hits[xx].qryCh != m->iid1)) ||
        ((hits[xx].refCh != m->iid1) && (hits[xx].qryCh != m->iid1))) {
      fprintf(stderr, "Augh!  Merge failure!  Not in sync!\n");
      fprintf(stderr, "hits["s32bitFMT"] = r: "u32bitFMT" "s32bitFMT" "s32bitFMT" q: "u32bitFMT" "s32bitFMT" "s32bitFMT"\n",
              hits[xx].iid,
              hits[xx].refCh, hits[xx].refBeg, hits[xx].refEnd,
              hits[xx].qryCh, hits[xx].qryBeg, hits[xx].qryEnd);
      exit(1);
    }

    if ((hits[xx].clump >= 0) &&
        (clumpOut[hits[xx].clump] == false)) {
      atacMatch  C;
      sprintf(C.matchuid,  "clump"s32bitFMTW(06), hits[xx].clump);
      sprintf(C.parentuid, ".");
      C.matchiid = 0;
      C.type[0] = 'c';
      C.type[1] = 0;
      C.iid1 = m->iid1;
      C.pos1 = clumpLoA[hits[xx].clump];
      C.len1 = clumpHiA[hits[xx].clump] - clumpLoA[hits[xx].clump];
      C.fwd1 = m->fwd1;
      C.iid2 = m->iid2;
      C.pos2 = clumpLoB[hits[xx].clump];
      C.len2 = clumpHiB[hits[xx].clump] - clumpLoB[hits[xx].clump];
      C.fwd2 = m->fwd2;

      C.print(stdout, AF.labelA(), AF.labelB());

      clumpOut[hits[xx].clump] = true;
    }

    if (hits[xx].clump >= 0)
      sprintf(m->parentuid, "clump"s32bitFMTW(06), hits[xx].clump);
    else
      sprintf(m->parentuid, ".");
    m->print(stdout, AF.labelA(), AF.labelB());

    xx++;
    m = AF.nextMatch('u');
  }

  return(0);
}
