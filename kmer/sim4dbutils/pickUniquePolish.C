#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "bio.h"
#include "sim4.H"

//  Derived from pickBestPolish.c.  We report only the single best
//  match, when it is obvious that there is EXACTLY one best match.
//
//  Example: we have ten matches, but one is 3%id better than everyone
//  else -- that is an obviously unique match.  The rest are noise.
//
//  Example: ten matches, but they're all about the same quality -- within
//  a few percent id, and about the same length.  We pick no match, and
//  silently discard all.
//

u32bit  statOneMatch      = 0;
u32bit  statConsistent    = 0;
u32bit  statInconsistent  = 0;
u32bit  statUnique        = 0;
u32bit  statLost          = 0;

u32bit  consistentTie         = 0;
u32bit  consistentMatches     = 0;
u32bit  consistentIdentity    = 0;
u32bit  consistentTooShort    = 0;
u32bit  consistentNot         = 0;

u32bit  totLQ = 0;
u32bit  totMQ = 0;
u32bit  totRQ = 0;

u32bit  qualityDifference = 5;

void
pickBestSlave(sim4polish **p, u32bit pNum) {
  u32bit        identitym = 0, nmatchesm = 0;  //  Best score for the mList
  u32bit        identityi = 0, nmatchesi = 0;  //  Best score the the iList
  u32bit        matchi = 0,    matchm = 0;

  //  Difficult choice here....
  //
  if (pNum == 1) {
    statOneMatch++;
    statUnique++;
    s4p_printPolish(stdout, p[0], S4P_PRINTPOLISH_FULL);
    return;
  }

  //  Find the best percentIdentity and best numberOfMatches.  
  //
  //  identityi is the best percent identity of all the matches for this EST, and
  //  nmatchesi is the number of matches for the longest best identity match(es).
  //  matchi    is the match index
  //
  //  nmatchesm is the best numMatches of all the matches for this EST, and 
  //  identitym is the highest percent identity for the best numMatches match(es).
  //  matchm    is the match index

  for (u32bit i=0; i<pNum; i++) {
    if ((p[i]->percentIdentity > identityi) || 
        (p[i]->percentIdentity == identityi && p[i]->numMatches > nmatchesi)) {
      identityi = p[i]->percentIdentity;
      nmatchesi = p[i]->numMatches;
      matchi    = i;
    }
   
    if ((p[i]->numMatches > nmatchesm) ||
        (p[i]->numMatches == nmatchesm && p[i]->percentIdentity > identitym)) {
      nmatchesm = p[i]->numMatches;
      identitym = p[i]->percentIdentity;
      matchm    = i;
    }
  }

  bool  matchIsOK = false;

  //  If we are in agreement on what the best quality match is,
  //  see if the best match is obviously unique.
  //
  if ((identityi == identitym) ||
      (nmatchesi == nmatchesm)) {
    statConsistent++;

    //  It's clear what the quality values of the best match is, but we
    //  don't know if those values are shared by more than one match.
    //  Count the number of matches with exactly those scores.  If
    //  there is more than one, then we cannot pick out a single best.
    //
    u32bit numBest = 0;
    for (u32bit i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) && (p[i]->numMatches == nmatchesi))
        numBest++;

    if (numBest > 1) {

      //  Dang, we mapped this guy more than once, exactly the same!
      //
      consistentTie++;

    } else {

      //  We claim to have a single best match.  See if any other
      //  matches are close to the quality of that one.
      //
      //  This says if (p[i]/ii >= 1.0 - Q), then we're close.

      u32bit  closeQuality = 0;
      for (u32bit i=0; i<pNum; i++)
        if (((p[i]->percentIdentity * 100) >= (identityi * (100 - qualityDifference))) ||
            ((p[i]->numMatches      * 100) >= (nmatchesi * (100 - qualityDifference))))
          closeQuality++;

      //  If only one match has close quality (the one we want to save!),
      //  save it.  Otherwise, label this query as multiple.

      u32bit  length = p[matchi]->exons[0].estFrom - p[matchi]->exons[0].estTo;

      if (closeQuality == 1) {
        matchIsOK = true;
        consistentMatches++;
      } else if ((length > 100) &&
                 (length / p[matchi]->estLen < 0.5)) {
        consistentTooShort++;
      } else {
        consistentNot++;
      }
    }

  } else {

    //  Otherwise, we disagree on what the best match is.
    //
    //  That is, the match with the highest identity is not the match
    //  with the highest number of matches -- a longer match exists, but
    //  at lower overall percent identity.

    statInconsistent++;

    //  Estimate the identity of the extended part, assuming the piece
    //  matched in common is matched at about the same identity.  Or
    //  just give up and say it's mapped to multiple places!

  }

  u32bit  best  = 0;
  u32bit  besti = 0;
  
  if (matchIsOK) {
    statUnique++;
    s4p_printPolish(stdout, p[matchi], S4P_PRINTPOLISH_FULL);

    assert(matchi == matchm);

    besti = matchi;
  } else {
    statLost++;

    //  Just pick the longest match, analyze that.

    for (u32bit i=0; i<pNum; i++) {
      u32bit  len = p[i]->exons[0].estFrom - p[i]->exons[0].estTo;

      if ((len  > best) ||
          ((len == best) && (p[i]->numMatches > p[besti]->numMatches))) {
        best  = len;
        besti = i;
      }
    }
  }

#if 1
  fprintf(stderr, "Uni:"u32bitFMTW(8)" Con:"u32bitFMTW(8)" (T:"u32bitFMTW(8)" M:"u32bitFMTW(8)" I:"u32bitFMTW(8)" N:"u32bitFMTW(8)") Inc:"u32bitFMTW(8)" -- Save:"u32bitFMTW(8)" Lost:"u32bitFMTW(8)"\r",
          statOneMatch,
          statConsistent, consistentTie, consistentMatches, consistentIdentity, consistentNot,
          statInconsistent,
          statUnique, statLost);
#endif
}






//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
pickBest(sim4polish **p, u32bit pNum) {

  pickBestSlave(p, pNum);

  for (u32bit i=0; i<pNum; i++)
    s4p_destroyPolish(p[i]);
}




int
main(int argc, char **argv) {
  u32bit       pNum         = 0;
  u32bit       pAlloc       = 8388608;
  sim4polish **p            = 0L;
  sim4polish  *q            = 0L;
  u32bit       estID        = ~0;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-n", 2) == 0) {
      pAlloc = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-q", 2) == 0) {
      qualityDifference = strtou32bit(argv[++arg], 0L);
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s < file > file\n", argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in the
  //  estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = s4p_readPolish(stdin)) != 0L) {
    if ((q->estID != estID) && (pNum > 0)) {
      pickBest(p, pNum);
      pNum  = 0;
    }

    //  Reallocate pointers?
    //
    if (pNum >= pAlloc) {
      p = (sim4polish **)realloc(p, sizeof(sim4polish *) * (pAlloc *= 2));
      if (p == 0L) {
        fprintf(stderr, "Out of memory: Couldn't allocate space for polish pointers.\n");
        exit(1);
      }
    }

    p[pNum++] = q;
    estID     = q->estID;
  }

  if (pNum > 0)
    pickBest(p, pNum);

  fprintf(stderr, "\n");

  return(0);
}

