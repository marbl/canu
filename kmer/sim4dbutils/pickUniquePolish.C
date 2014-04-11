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

uint32  statOneMatch      = 0;
uint32  statConsistent    = 0;
uint32  statInconsistent  = 0;
uint32  statUnique        = 0;
uint32  statLost          = 0;

uint32  consistentTie         = 0;
uint32  consistentMatches     = 0;
uint32  consistentIdentity    = 0;
uint32  consistentTooShort    = 0;
uint32  consistentNot         = 0;

uint32  totLQ = 0;
uint32  totMQ = 0;
uint32  totRQ = 0;

uint32  qualityDifference = 5;
uint32  minQuality        = 95;

sim4polishWriter *W = 0L;



void
pickUniqueSlave(sim4polish **p, uint32 pNum) {
  uint32        identitym = 0, nmatchesm = 0;  //  Best score for the mList
  uint32        identityi = 0, nmatchesi = 0;  //  Best score the the iList
  uint32        matchi = 0,    matchm = 0;

 //  Difficult choice here....
  //
  if (pNum == 1) {
    statOneMatch++;
    statUnique++;
    W->writeAlignment(p[0]);
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

  for (uint32 i=0; i<pNum; i++) {
    if ((p[i]->_percentIdentity > identityi) || 
        (p[i]->_percentIdentity == identityi && p[i]->_numMatches > nmatchesi)) {
      identityi = p[i]->_percentIdentity;
      nmatchesi = p[i]->_numMatches;
      matchi    = i;
    }
   
    if ((p[i]->_numMatches > nmatchesm) ||
        (p[i]->_numMatches == nmatchesm && p[i]->_percentIdentity > identitym)) {
      nmatchesm = p[i]->_numMatches;
      identitym = p[i]->_percentIdentity;
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
    uint32 numBest = 0;
    for (uint32 i=0; i<pNum; i++)
      if ((p[i]->_percentIdentity == identityi) && (p[i]->_numMatches == nmatchesi))
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

      uint32  closeQuality = 0;
      for (uint32 i=0; i<pNum; i++)
        if (((p[i]->_percentIdentity * 100) >= (identityi * (100 - qualityDifference))) ||
            ((p[i]->_numMatches      * 100) >= (nmatchesi * (100 - qualityDifference))))
          closeQuality++;

      //  If only one match has close quality (the one we want to save!),
      //  save it.  Otherwise, label this query as multiple.

      uint32  length = p[matchi]->_exons[0]._estFrom - p[matchi]->_exons[0]._estTo;

      if (closeQuality == 1) {
        matchIsOK = true;
        consistentMatches++;
      } else if ((length > 100) &&
                 (length / p[matchi]->_estLen < 0.5)) {
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

 
  if (matchIsOK) {
    statUnique++;
    assert(matchi == matchm);
    W->writeAlignment(p[matchi]);
  } else {
    statLost++;
  }
}





//  Delete all matches that are spanned, report everything else.
//  Matches that are close ties in span, but are clearly lower quality are deleted.
//
void
pickCoveringSlave(sim4polish **p, uint32 pNum, char doCovering) {
  uint32 *bgn = new uint32 [pNum];
  uint32 *end = new uint32 [pNum];

  for (uint32 i=0; i<pNum; i++) {
    if (doCovering == 'q') {
      if (p[i]->_matchOrientation == SIM4_MATCH_FORWARD) {
        bgn[i] = p[i]->_exons[0]._estFrom - 1;
        end[i] = p[i]->_exons[0]._estTo;
      } else {
        bgn[i] = p[i]->_estLen - p[i]->_exons[0]._estTo;
        end[i] = p[i]->_estLen - p[i]->_exons[0]._estFrom + 1;
      }
    }

    if (doCovering == 'g') {
      bgn[i] = p[i]->_exons[0]._genFrom - 1;
      end[i] = p[i]->_exons[0]._genTo;
    }
  }


  for (uint32 i=0; i<pNum; i++) {
    if (p[i] == NULL)
      continue;

    assert(p[i]->_numExons == 1);

    for (uint32 j=i+1; j<pNum; j++) {
      if (p[j] == NULL)
        continue;

      //  i contained in j
      //         ----
      //       ---------
      if ((bgn[j] <= bgn[i]) && (end[i] <= end[j])) {
        delete p[i];  p[i] = NULL;
        break;  //  This i is finished.
      }

      //  j contained in i
      //       ---------
      //         ----
      if ((bgn[i] <= bgn[j]) && (end[j] <= end[i])) {
        delete p[j];  p[j] = NULL;
        continue;  //  This j is finished.
      }

      // i almost contained in j
      //      ----                    ----
      //       ---------  OR   ----------
      if (((bgn[j] <= bgn[i] + 5) && (end[i] <= end[j])) ||
          ((bgn[j] <= bgn[i])     && (end[i] <= end[j] + 5))) {
        delete p[i];  p[i] = NULL;
        break;  //  This i is finished.
      }

      // j almost contained in i
      //       ---------  OR   ----------
      //      ----                    ----
      if (((bgn[i] <= bgn[j] + 5) && (end[j] <= end[i])) ||
          ((bgn[i] <= bgn[j])     && (end[j] <= end[i] + 5))) {
        delete p[j];  p[j] = NULL;
        continue;  //  This j is finished.
      }
    }
  }
  
  for (uint32 i=0; i<pNum; i++) {
    if (p[i] == NULL)
      continue;

    W->writeAlignment(p[i]);
  }

  delete [] bgn;
  delete [] end;
}






//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
pickUnique(sim4polish **p, uint32 pNum, char doCovering) {

  if (doCovering != 0)
    pickCoveringSlave(p, pNum, doCovering);
  else
    pickUniqueSlave(p, pNum);

  for (uint32 i=0; i<pNum; i++)
    delete p[i];
}




int
main(int argc, char **argv) {
  char         doCovering     = 0;
  uint32       pNum           = 0;
  uint32       pAlloc         = 1048576;
  uint32       lastID         = ~uint32ZERO;

  sim4polishStyle  style = sim4polishStyleDefault;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-cq") == 0) {
      doCovering = 'q';

    } else if (strcmp(argv[arg], "-cg") == 0) {
      doCovering = 'g';

    } else if (strcmp(argv[arg], "-q") == 0) {
      qualityDifference = strtouint32(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      style = sim4polishGFF3;

    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, "usage: %s [-q qualDiff] [-c] [-1] [-gff3] < file > file\n", argv[0]);
    fprintf(stderr, "  -q qualDiff    Only report alignments where the best is qualDiff better\n");
    fprintf(stderr, "                 in percent identity and coverage\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -cq            Only report alignments that are not contained in some\n");
    fprintf(stderr, "                 other alignment in the QUERY SEQUENCE.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -cg            Only report alignments that are not contained in some\n");
    fprintf(stderr, "                 other alignment in the GENOMIC SEQUENCE.\n");
    fprintf(stderr, "\n");

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in the
  //  estID.

  sim4polishReader  *R = new sim4polishReader("-");
  sim4polish       **p = new sim4polish * [pAlloc];
  sim4polish        *q = 0L;

  W = new sim4polishWriter("-", style);

  if (R->getsim4polishStyle() != style)
    fprintf(stderr, "warning: input format and output format differ.\n");

  while (R->nextAlignment(q)) {
    bool doPick = false;

    if ((doCovering == 'q') && (q->_estID != lastID))
      doPick = true;

    if ((doCovering == 'g') && (q->_genID != lastID))
      doPick = true;

    if ((doCovering == 0) && (q->_estID != lastID))
      doPick = true;

    if ((doPick == true) && (pNum > 0)) {
      pickUnique(p, pNum, doCovering);
      pNum  = 0;
    }

    if (pNum >= pAlloc) {
      sim4polish **P = new sim4polish * [pAlloc * 2];
      memcpy(p, P, sizeof(sim4polish *) * pAlloc);
      delete [] p;
      p = P;
      pAlloc *= 2;
    }

    p[pNum++] = q;
    lastID    = (doCovering == 'g') ? q->_genID : q->_estID;

    q = 0L;  //  Otherwise we delete the alignment we just saved!
  }

  if (pNum > 0)
    pickUnique(p, pNum, doCovering);

#if 0
  fprintf(stderr, "Uni:"uint32FMTW(8)" Con:"uint32FMTW(8)" (T:"uint32FMTW(8)" M:"uint32FMTW(8)" I:"uint32FMTW(8)" N:"uint32FMTW(8)") Inc:"uint32FMTW(8)" -- Save:"uint32FMTW(8)" Lost:"uint32FMTW(8)"\n",
          statOneMatch,
          statConsistent, consistentTie, consistentMatches, consistentIdentity, consistentNot,
          statInconsistent,
          statUnique, statLost);
#endif

  delete [] p;

  delete R;
  delete W;

  return(0);
}

