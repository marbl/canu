#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sim4reader.h"

//  Picks the best polish (or set of polishes that are all of the same
//  best quality) for each cDNA.


#define    EPS_X	1
#define    EPS_N_ESTS	10
#define    EPS_N_MRNA	15
#define    EPS_I	3


char const *usage =
"usage: %s [-mrna|-ests] < file > file\n"
"\n";


static int EPS_N;

//  Validate mode will print out ALL input matches, in the following
//  format
//
//  estid gaid nummatches percentid (genFr genTo %) () ()
//
//  With a * somewhere to denote the best ones.  Separate ESTs with
//  a dashed line.
//
//#define VALIDATE

#ifdef VALIDATE
void
printPolishValidate(FILE *O, sim4polish *p, int isBest) {
  int i;

  fprintf(O, "%8d %8d %4d %4d", p->estID, p->genID, p->percentIdentity, p->numMatches);
  for (i=0; i<p->numExons; i++)
    fprintf(O, " (%6d/%6d %6d/%6d %3d)", p->exons[i].estFrom, p->genLo + p->exons[i].genFrom, p->exons[i].estTo, p->genLo + p->exons[i].genTo, p->exons[i].percentIdentity);

  if (isBest)
    fprintf(O, " *");

  fprintf(O, "\n");
}
#endif

void
pickBestSlave(sim4polish **p, int pNum) {
  int           i = 0;
  int           identitym = 0, nmatchesm = 0;  //  Best score for the mList
  int           identityi = 0, nmatchesi = 0;  //  Best score the the iList
  int           numExons = 0, numExonsi = 0, numExonsm = 0;
  int		tmp_nmatches = 0;
  double        alpha;

  //  Difficult choice here....
  //
  if (pNum == 1) {
#ifndef VALIDATE
    printPolish(stdout, p[0]);
#endif
    return;
  }

  if ((p[0]->estID % 1287) == 0) {
    fprintf(stderr, "Picking Best for estID=%8d with %5d choices.\r", p[0]->estID, pNum);
    fflush(stderr);
  }

  //  Find the best percentIdentity and best numberOfMatches.  
  //
  //  identityi is the best percent identity of all the matches for this EST, and
  //  nmatchesi is the number of matches for the longest best identity match(es).
  //
  //  nmatchesm is the best numMatches of all the matches for this EST, and 
  //  identitym is the highest percent identity for the best numMatches match(es).

  for (i=0; i<pNum; i++) {

    if ((p[i]->percentIdentity > identityi) || 
        (p[i]->percentIdentity == identityi && p[i]->numMatches > nmatchesi)) {
      identityi = p[i]->percentIdentity;
      nmatchesi = p[i]->numMatches;
    }
   
    if ((p[i]->numMatches > nmatchesm) ||
        (p[i]->numMatches == nmatchesm && p[i]->percentIdentity > identitym)) {
      nmatchesm = p[i]->numMatches;
      identitym = p[i]->percentIdentity;
    }

  }


  //  Otherwise, if the best scores on both lists are the same, pick
  //  the matches with the largest number of exons
  //
  if ((identityi == identitym) &&
      (nmatchesi == nmatchesm)) {

    //  Find the largest number of exons, allowing some margin in numMatches
    //
    numExonsi = 0;
    for (i=0; i<pNum; i++) 
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches >= nmatchesi) &&
          (numExonsi < p[i]->numExons))
          numExonsi = p[i]->numExons;

    numExons = numExonsi;
    tmp_nmatches = nmatchesi;
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches >= nmatchesi - EPS_N) &&
          (numExons < p[i]->numExons - EPS_X)) {
          tmp_nmatches = p[i]->numMatches;
          numExons = p[i]->numExons;
      }

    //  Scan the entire list, printing the best stuff.  We cannot just
    //  scan both the mList and iList, as those probably contain
    //  duplicates.

#ifndef VALIDATE
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches == tmp_nmatches) &&
          (p[i]->numExons == numExons))
        printPolish(stdout, p[i]);
#else
    if (tmp_nmatches == nmatchesi) 
    fprintf(stdout, "--------------------1 (Clear Winner)\n");
    else
    fprintf(stdout, "--------------------2 (Exon Clear Winner)\n");
    for (i=0; i<pNum; i++)
      printPolishValidate(stdout, p[i], ((p[i]->percentIdentity == identityi) &&
                                         (p[i]->numMatches == tmp_nmatches) &&
                                         (p[i]->numExons == numExons)));
#endif

    return;
  }

  //  Start over.  Find the best two percentIdentities.  Break ties
  //  with numMatches.
  //
  //  i will be the best,
  //  m will be the second best
  //
  identityi = identitym = 0;
  nmatchesi = nmatchesm = 0;

  for (i=0; i<pNum; i++) {

    //  Pick the two matches with the highest (different) percent
    //  identities; for each, pick the highest number of matches.
    //
    //  First block:  Have we found a new best percent identity?
    //  If so, save it, and shift former best to second best.
    //
    //  Second and third blocks: make sure that we save the
    //  best numMatches for each.
    //
    if (p[i]->percentIdentity > identityi) {
      identitym = identityi;
      nmatchesm = nmatchesi;

      identityi = p[i]->percentIdentity;
      nmatchesi = p[i]->numMatches;
    } else if ((p[i]->percentIdentity == identityi) &&
               (p[i]->numMatches > nmatchesi)) {
      nmatchesi = p[i]->numMatches;
    } else if ((p[i]->percentIdentity < identityi) &&
               ((p[i]->percentIdentity > identitym) ||
                ((p[i]->percentIdentity == identitym) &&
                (p[i]->numMatches       > nmatchesm)))) {
      nmatchesm = p[i]->numMatches;
      identitym = p[i]->percentIdentity;
    }
  }

  //  Now, 'i' is the highest percent identity, 'm' is the second
  //  highest.  By definition, numMatches for 'i' is less than
  //  numMatches for 'm'.

  //  If the number of matches is different, output everything with the
  //  top score.
  //
  //  We are guaranteed that the identities are the same.  (I think)

  if (nmatchesi >= nmatchesm) {

    //  Find the match(es) with the largest number of exons
    
    numExonsi = 0;
    for (i=0; i<pNum; i++) 
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches >= nmatchesi) &&
          (numExonsi < p[i]->numExons)) 
          numExonsi = p[i]->numExons;


    numExons = numExonsi;
    tmp_nmatches = nmatchesi;
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches >= nmatchesi - EPS_N) &&
          (numExons < p[i]->numExons - EPS_X)) {
           numExons = p[i]->numExons;
           tmp_nmatches = p[i]->numMatches;
      }

#ifndef VALIDATE
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches      == tmp_nmatches) &&
          (p[i]->numExons	 == numExons))
        printPolish(stdout, p[i]);
#else
    if (tmp_nmatches == nmatchesi)
    fprintf(stdout, "--------------------3 (?)\n");
    else 
    fprintf(stdout, "--------------------4 (Exon ?)\n");
    for (i=0; i<pNum; i++)
      printPolishValidate(stdout, p[i], ((p[i]->percentIdentity == identityi) &&
                                         (p[i]->numMatches      == tmp_nmatches) &&
                                         (p[i]->numExons 	== numExons)));
#endif
    return;
  }

  //  Otherwise, compute alpha

  alpha = ((nmatchesm - nmatchesi) / 
           ((nmatchesm / (double)identitym) -
            (nmatchesi / (double)identityi)))/100;

  //  If alpha below a magic threshold, pick the shorter match.
  //
  if (alpha < 0.8) {

    //  Find the match(es) with the largest number of exons

    numExons = tmp_nmatches = 0;
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches >= nmatchesi) &&
          (numExons < p[i]->numExons))
           numExons = p[i]->numExons;

#ifndef VALIDATE
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identityi) &&
          (p[i]->numMatches      == nmatchesi) &&
          (p[i]->numExons 	 == numExons))
        printPolish(stdout, p[i]);
#else
    fprintf(stdout, "--------------------5 (alpha < 0.8)\n");
    for (i=0; i<pNum; i++)
      printPolishValidate(stdout, p[i], ((p[i]->percentIdentity == identityi) &&
                                         (p[i]->numMatches      == nmatchesi) &&
                                         (p[i]->numExons	== numExons)));
#endif
    return;
  }

  //  Otherwise, pick the longer one.

  //  XXX: We can still check:
  //    if an internal gap is in N's
  //    the number of exons
  //    etc, etc.



  //  See if the smaller one has an internal gap that corresponds to
  //  N's in the genome.  If so, assume that the exon mapped to the
  //  N's and pick the smaller.
  //
  //  Need code to process genome, finding N's larger than some threshold.
  //  Output as 'genID beg end'


  //  Find the largest number of exons for each of the contenders
   
  numExonsi = numExonsm = 0;
  for (i=0; i<pNum; i++) {
    if ((p[i]->percentIdentity == identitym) &&
        (p[i]->numMatches == nmatchesm) &&
        (numExonsm < p[i]->numExons))
         numExonsm = p[i]->numExons;
    else if ((p[i]->percentIdentity == identityi) &&
        (p[i]->numMatches == nmatchesi) &&
        (numExonsi < p[i]->numExons))
         numExonsi = p[i]->numExons;
  }

  if ((numExonsi > numExonsm + EPS_X) || (identityi > identitym + EPS_I)) {

#ifndef VALIDATE
  for (i=0; i<pNum; i++)
    if ((p[i]->percentIdentity == identityi) &&
        (p[i]->numMatches      == nmatchesi) &&
        (p[i]->numExons	       == numExonsi))
      printPolish(stdout, p[i]);
#else
  if (numExonsi > numExonsm + EPS_X)
  fprintf(stdout, "--------------------6 (Exon Plus alpha > 0.8)\n");
  else 
  fprintf(stdout, "--------------------7 (Pctid Plus alpha > 0.8)\n");

  for (i=0; i<pNum; i++)
    printPolishValidate(stdout, p[i], ((p[i]->percentIdentity == identityi) &&
                                       (p[i]->numMatches      == nmatchesi) &&
                                       (p[i]->numExons	      == numExonsi)));
#endif
  } else {
    numExons = numExonsm;
    tmp_nmatches = nmatchesm;
    for (i=0; i<pNum; i++)
      if ((p[i]->percentIdentity == identitym) &&
          (p[i]->numMatches >= nmatchesm - EPS_N) &&
          (numExons < p[i]->numExons - EPS_X)) {
          tmp_nmatches = p[i]->numMatches;
          numExons = p[i]->numExons;
      }
 
#ifndef VALIDATE
  for (i=0; i<pNum; i++)
    if ((p[i]->percentIdentity == identitym) &&
        (p[i]->numMatches      == tmp_nmatches) &&
        (p[i]->numExons        == numExons))
      printPolish(stdout, p[i]);
#else
  if (numExons == numExonsm)
  fprintf(stdout, "--------------------8 (alpha > 0.8)\n");
  else 
  fprintf(stdout, "--------------------9 (Exon alpha > 0.8)\n");
  for (i=0; i<pNum; i++)
    printPolishValidate(stdout, p[i], ((p[i]->percentIdentity == identitym) &&
                                       (p[i]->numMatches      == tmp_nmatches) &&
                                       (p[i]->numExons        == numExons)));
#endif
  }
}


//  Just a wrapper around the real best picker, so that we can easily
//  destroy polishes when we're done.
//
void
pickBest(sim4polish **p, int pNum) {
  int i;

  pickBestSlave(p, pNum);

  for (i=0; i<pNum; i++)
    destroyPolish(p[i]);
}


int
main(int argc, char **argv) {
  int          pNum   = 0;
  int          pAlloc = 8388608;
  sim4polish **p      = 0L;
  sim4polish  *q      = 0L;
  int          estID  = ~0;

  int          isEST  = -1;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-n", 2) == 0) {
      pAlloc = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-mrna", 2) == 0) {
      if (isEST == 1) {
        fprintf(stderr, "error: input type already set to EST!\n\n"); exit(1);
      }
      isEST = 0;
      EPS_N = EPS_N_MRNA;
    } else if (strncmp(argv[arg], "-ests", 2) == 0) {
      if (isEST == 0) {
        fprintf(stderr, "error: input type already set to mRNA!\n\n"); exit(1);
      }
      isEST = 1;
      EPS_N = EPS_N_ESTS;
    } else {
      fprintf(stderr, "unknown option: %s\n", argv[arg]);
    }
    arg++;
  }

  if (isEST == -1) {
    fprintf(stderr, "warning: input type unspecified; setting to EST.\n\n"); 
    isEST = 1;
    EPS_N = EPS_N_ESTS;
  }

  if (isatty(fileno(stdin)) || isatty(fileno(stdout))) {
    fprintf(stderr, usage, argv[0]);

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n\n");

    if (isatty(fileno(stdout)))
      fprintf(stderr, "error: Please redirect the output to a file.\n\n");

    exit(1);
  }

  //  Read polishes, picking the best when we see a change in
  //  the estID.

  p = (sim4polish **)malloc(sizeof(sim4polish *) * pAlloc);

  while ((q = readPolish(stdin)) != 0L) {
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

  return(0);
}

