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

u32bit statOneMatch      = 0;
u32bit statConsistent    = 0;
u32bit statInconsistent  = 0;
u32bit statUnique        = 0;
u32bit statLost          = 0;
u32bit uniqueMatches     = 0;
u32bit uniqueIdentity    = 0;
u32bit uniqueNot         = 0;





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

  if ((p[0]->estID % 1287) == 0) {
    fprintf(stderr, "Picking Best for estID="u32bitFMT" with %5d choices.\r", p[0]->estID, pNum);
    fflush(stderr);
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
  //
  //
  if ((identityi == identitym) ||
      (nmatchesi == nmatchesm)) {

    statConsistent++;

    //  look for the next best match.  If he is much shorter, or much lower identity,
    //  we can report our best match.

    //  Look for the next highest identity match
    //
    u32bit id = 0;
    u32bit nm = 0;

    //  Find the second largest percent identity and number of matches
    for (u32bit i=0; i<pNum; i++) {
      if ((p[i]->percentIdentity > id) && (p[i]->percentIdentity < identityi))
        id = p[i]->percentIdentity;
      if ((p[i]->numMatches > nm) && (p[i]->numMatches < nmatchesi))
        nm = p[i]->numMatches;
    }

    //  If the next longest is significantly shorter, print our best
    //  match.  We know this match has the highest identity, so the
    //  next longest match has lower identity, and is shorter.
    //
    //  Likewise, if the next highest identity match is significantly
    //  worse (we already know it's not longer), it's clear this is
    //  still the better match.
    //
    //  Otherwise, we declare it's not a unique mapping.
    //
    if (nmatchesi * 0.98 >= nm){
      uniqueMatches++;
      matchIsOK = true;
    } else if (identityi * 0.98 >= id) {
      uniqueIdentity++;
      matchIsOK = true;
    } else {
      uniqueNot++;
    }

  } else {

    //  Otherwise, we disagree on what the best match is.
    //
    //  That is, the match with the highest identity is not the match
    //  with the highest number of matches -- a longer match exists, but
    //  at lower overall percent identity.

    statInconsistent++;

    //  Estimate the identity of the extended part, assuming the piece
    //  matched in common is matched at about the same identity.

  }


  if (matchIsOK) {
    statUnique++;
    s4p_printPolish(stdout, p[matchi], S4P_PRINTPOLISH_FULL);
  } else {
    statLost++;
  }

  fprintf(stderr, "Uni:"u32bitFMTW(8)" Con:"u32bitFMTW(8)" (M:"u32bitFMTW(8)" I:"u32bitFMTW(8)" N:"u32bitFMTW(8)") Inc:"u32bitFMTW(8)" -- Save:"u32bitFMTW(8)" Lost:"u32bitFMTW(8)"\r",
          statOneMatch,
          statConsistent, uniqueMatches, uniqueIdentity, uniqueNot,
          statInconsistent,
          statUnique, statLost);
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
main(u32bit argc, char **argv) {
  u32bit       pNum   = 0;
  u32bit       pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  u32bit       estID  = ~0;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-n", 2) == 0) {
      pAlloc = strtou32bit(argv[++arg], 0L);
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

