
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

/* Module for finding micro-satellites at the end of sequence reads.
   Designed to find micros with models of not more than 6 bases at
   not more than 15% variation/error.  Speeds in thousands of sequence
   ends per second have been timed on my IBM 600E laptop as follows:

                Random Seq: 82.7K

   End Containing a: 6-mer:  2.2K
                     5-mer:  1.5K
                     4-mer:  1.0K
                     3-mer:   .7K
                     2-mer:   .55K
                     1-mer:   .40K

   These are the speeds just to check that there is a match to the first
   40bp's of a sequence end.  The time to run the branchpoint finder on
   the sequence, in the event that a micro-sat is found, must further be
   factored into the numbers above.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "AS_ALN_aligners.h"

#undef  TESTMICRO
//#define TESTMICRO
#undef TESTRANDOM

#undef   PRINT_MODELS
#undef   DEBUG_DP
#undef   PRINT_MATCHES

#define TESTLEN 40  /* Length of suffix/prefix to be tested for micro-sat. */
#define KMERLEN  3  /* Tuple length to base frequency count filter upon. */
#define ERRORS   6  /* Maximum # of errors permitted in a match. */

#define KMERNUM   64  /* == 4^KMERLEN */
#define KMERMSK   0xF /* == 4^(KMERLEN-1) - 1 */
#define TESTMER   38  /* == TESTLEN - (KMERLEN-1) */
#define CUTOFF    20  /* == TESTMER - KMERLEN*ERRORS */
#define MICRNUM 8192  /* == 4^(MICROMX+1) */

/* Compute the sum of the frequencies of the MICROMX(6) most
   frequent KMERLEN(3)-mers in the first (iff chk_prefix != 0) or
   last TESTLEN characters of "sequence".  Return this sum and also
   an array of the frequency counts in "bucket".  A string of length
   TESTLEN cannot possibly contain a micro-sat match with ERRORS errors
   if this sum is less than CUTOFF.
*/

static int Micro_Counter_Screen(char *sequence, int *bucket, int chk_prefix)
{ static int firstime = 1;
  static int code[128];

  int len, pop;

  if (firstime)
    { int i;

      firstime = 0;
      for (i = 0; i < 128; i++)   /* Setup DNA -> radix 4 map */
        code[i] = -1;
      code['a'] = code['A'] = 0;
      code['c'] = code['C'] = 1;
      code['g'] = code['G'] = 2;
      code['t'] = code['T'] = 3;
    }

  len = strlen(sequence);
  if (len < TESTLEN) return (0);

  /* Count # of occurences of each K-mer (3-mer) in the prefix
     or suffix of the sequence of length TESTLEN (40)          */

  { int i;

    for (i = 0; i < KMERNUM; i++)
      bucket[i] = 0;
  }

  { char *seq;
    int i, tuple;

    if (chk_prefix)
      seq = sequence;
    else
      seq = sequence+(len-TESTLEN);

    seq -= 1;
    tuple = 0;
    for (i = 1; i < KMERLEN; i++)
      tuple = (tuple << 2) | code[(int) seq[i]];
    for (i = KMERLEN; i <= TESTLEN; i++)
      { tuple = ((tuple & KMERMSK) << 2) | code[(int) seq[i]];
        bucket[tuple] += 1;
      }
  }

  /* Compute the sum of the top MICROMX(6) frequency counts in "pop" */

  { int i, x, cnt;
    int sort[TESTMER+1];

    for (i = 0; i <= TESTMER; i++) 
      sort[i] = 0;
    for (i = 0; i < KMERNUM; i++)
      sort[bucket[i]] += 1;

    cnt = pop = 0;
    for (i = TESTMER; i > 0; i--)
      for (x = sort[i]; x > 0; x--)
        { pop += i;
          cnt += 1;
          if (cnt >= MICROMX || pop >= TESTMER) goto done;
        }
  }
done:

  return (pop);
}

/* Compare the first (iff prefix != 0) or last TESTLEN characters of "seq"
   against the micro-sat "micro" of length "mlen".  NB: "micro" is not a
   '\0'-terminated string.  Methodology is wrap-around dynamic programming.
   Return # of differences in best match or a value > ERRORS if more than
   ERRORS differences are in a best match.
*/

static int Wrap_Around(char *seq, char *micro, int mlen, int ispref)
{ int x, y;
  int i, j;
  int *D, *C;
  int min;

  static int row1[MICROMX], row2[MICROMX];

  if (!ispref)
    seq = seq + (strlen(seq) - TESTLEN);

  C = row1;
  D = row2;
  for (i = 0; i < mlen; i++)
    D[i] = 0;

#ifdef DEBUG_DP
  printf("\n   ");
  for (i = 0; i < mlen; i++)
    printf("  %c",micro[i]);
  printf("\n");
#endif

  for (j = 0; j < TESTLEN; j++)
    { x = D[mlen-1] + (micro[0] != seq[j]);
      y = D[0]+1;
      if (y < x) x = y;
      C[0] = x;
      for (i = 1; i < mlen; i++)
        { x = D[i-1] + (micro[i] != seq[j]);
          y = D[i]+1;
          if (y < x) x = y;
          y = C[i-1]+1;
          if (y < x) x = y;
          C[i] = x;
        }
      min = x = C[mlen-1]+1;
      for (i = 0; i < mlen; i++)
        { if (x < C[i]) C[i] = x;
          if (C[i] < min) min = C[i];
          x = x+1;
        }

#ifdef DEBUG_DP
      printf("%c: ",seq[j]);
      for (i = 0; i < mlen; i++)
        printf("%3d",C[i]);
      printf(" --> %d\n",min);
#endif

      if (min > ERRORS) return (min);  /* Get on out early if already have
                                          more than ERRORS errors */
      { int *x; x = C; C = D; D = x; }
    }

  return (min);
}

/* Convert integer code "val" for a DNA string of length "len" back
   into a '\0'-terminated string.
*/

#define MAXDECODE 16

static char *DecodeTuple(int val, int len)
{ static int  decode[4];
  static int  first = 1;
  static char construct[MAXDECODE+1];

  if (first)
    { first = 0;
      decode[0] = 'a';
      decode[1] = 'c';
      decode[2] = 'g';
      decode[3] = 't';
    }
  if (len > MAXDECODE)
    { fprintf(stderr,"DecodeTuple buffer is too small -- enlarge\n");
      exit (1);
    }

  { int i;
    for (i = --len; i >= 0; i--)
      construct[len-i] = decode[(val >> 2*i) & 0x3];
    construct[len+1] = '\0';
  }

  return (construct);
}

/* Given sequence "seq", array of 3-mer frequency counts in "bucket",
   and the 3-mer "seed" that must be in the micro-sat model/prototype,
   find any matching models and return a pointer to a list of them.
   The list is -1-terminated and each preceeding non-negative integer N
   encodes a model as follows: (N>>12) is the length of the model, and
   (N&0xFFF) is the integer code of the model.  Thus an ASCII encoding
   of the encoded model is given by "DecodeTuple(N&0xFFF,N>>12)".
      If "newlist" is zero then Explore should add new models to the
   list of candidates it returns, otherwise it starts a new list.  If
   "newlist" is zero, then the "bucket" frequency count a previous
   choice of "seed" is made negative to indicate that models with the
   triple denoted by that seed have already been explored and that this
   invocation of call should not explore models containing the previously
   explored triple as a substring.
      Finally, "ispref" is non-zero if the model is to match a prefix of
   "seq" as opposed to a suffix.
      This routine is tailored specifically for KMERLEN = 3 and
   MICROMX = 6.  It does not work for other values.
*/

static int *Explore(char *seq, int *bucket, int seed, int newlist, int ispref)
{ int suf2, suf1;
  int pre2, pre1;
  int avd4 = 0, avd5 = 0, avd6 = 0;
  int tup1, tup2, tup3, tup4, tup5;
  int cnt0, cnt1, cnt2, cnt3, cnt4, cnt5;
  int i, period, dif;
  char *model;

  static int first = 1;
  static int used[KMERNUM];

  static int match[MICRNUM];
  static int matop, madif;

#define MATCH(len,code,dif)		\
{ if (dif <= ERRORS && dif <= madif)	\
    { if (dif < madif) matop = 0;	\
      match[matop] = (len<<12) + code;	\
      matop += 1;			\
      madif  = dif;			\
   }					\
}

  if (first)
    { first = 0;
      for (i = 0; i < KMERNUM; i++)
        used[i] = 0;
    }

  if (newlist)
    { matop = 0;
      madif = ERRORS+1;
    }

#define BUCKET(x) (bucket[x])

  /* Suppose seed == a1.a2.a3 as a DNA string */

  suf2 = (seed & 0xF) << 2;  /* a2.a3.__ */
  suf1 = (seed & 0x3) << 4;  /* a3.__.__ */
  pre2 = seed >> 2;          /* __.a1.a2 */
  pre1 = seed >> 4;          /* __.__.a1 */

  /* Periodicity of seed, if not set to 1 or 2, then its 3. */

  period = 3;

  /* Check implied 1-model (if any) */

  cnt0 = BUCKET(seed);
  if (seed % 21 == 0)       /* ==> (a1)* a model */

    { period = 1;           /* Seed has period 1 */
      avd4 = (seed & 0x3);  /* Don't let i = a3 for 4-models */
      avd5 = (seed & 0xF);  /* Don't let i = a2.a3 for 5-models */

      /* Sum of els above cutoff? */

      if (cnt0 >= CUTOFF)
        { model = DecodeTuple(pre1,1);
          dif   = Wrap_Around(seq,model,1,ispref);
	  MATCH(1,pre1,dif)
#ifdef PRINT_MODELS
          printf("Cycle of length 1: %s",model);
          printf(" -> %d\n",dif);
#endif
        }
    }

  /* Check implied 2-model (if any) */
  
  if (period != 1 && (seed & 0x3) == pre1)  /* ==> (a1.a2)* a model */

    { period = 2;            /* Seed has period 2 */
      avd4 = (pre2 & 0x3);   /* Don't let i = a2 for 4-models */
      avd6 = (suf2 | pre2);  /* Don't let i = a2.a1.a2 for 6-models */

      /* Doesn't involve a previous seed & sum of els above cutoff? */

      cnt1 = BUCKET(avd6);
      if (cnt1 >= 0)
        if (cnt0 + cnt1 >= CUTOFF)

        { model = DecodeTuple(pre2,2);
          dif   = Wrap_Around(seq,model,2,ispref);
          MATCH(2,pre2,dif)
#ifdef PRINT_MODELS
          printf("Cycle of length 2: %s",model);
          printf(" -> %d\n",dif);
#endif
        }
    }

  /* Check implied 3-model */

  if (period != 1)    /* ==> (a1.a2.a3)* is a model not yet explored */

    { cnt1 = BUCKET(suf2 | pre1);
      cnt2 = BUCKET(suf1 | pre2);

      /* Doesn't involve a previous seed & sum of els above cutoff? */

      if (cnt1 >=0 && cnt2 >=0)
        if (cnt0 + cnt1 + cnt2 >= CUTOFF)

        { model = DecodeTuple(seed,3);
          dif   = Wrap_Around(seq,model,3,ispref);
	  MATCH(3,seed,dif)
#ifdef PRINT_MODELS
          printf("Cycle of length 3: %s",model);
          printf(" -> %d\n",dif);
#endif
        }
    }

  /* Check 4-models a1.a2.a3.i for all i in [0,3] */

  for (i = 0; i < 4; i++)

    /* Avoid re-exploring smaller models */

    if (period == 3 || i != avd4)

      { cnt1 = BUCKET(suf2 | i);
        cnt2 = BUCKET(suf1 | (i << 2) | pre1);
        cnt3 = BUCKET((i << 4) | pre2);

        /* Doesn't involve a previous seed & sum of els above cutoff? */

        if (cnt1 >=0 && cnt2 >=0 && cnt3 >=0)
          if (cnt0 + cnt1 + cnt2 + cnt3 >= CUTOFF)

          { model = DecodeTuple(seed*4+i,4);
            dif   = Wrap_Around(seq,model,4,ispref);
	    MATCH(4,seed*4+i,dif)
#ifdef PRINT_MODELS
            printf("Cycle of length 4: %s",model);
            printf(" -> %d\n",dif);
#endif
          }
      } 

  /* Check 5-models a1.a2.a3.i for all i in [0,15] */

  for (i = 0; i < 16; i++)

    /* Avoid re-exploring smaller models */

    if (period != 1 || i != avd5)
      { if (period == 1 && pre1 == (i >> 2)) continue;
        if (period == 2 && pre2 == i) continue; 

        tup1 = (suf2 | (i >> 2));
        tup2 = (suf1 | i);
        tup3 = ((i << 2) | pre1);
        tup4 = (((i & 0x3) << 4) | pre2);
        cnt1 = BUCKET(tup1);
        cnt2 = BUCKET(tup2);
        cnt3 = BUCKET(tup3);
        cnt4 = BUCKET(tup4);

        /* Don't use a previous seed & count repeated els only once */

        if (cnt1 < 0 || cnt2 < 0 || cnt3 < 0 || cnt4 < 0) continue;
        used[seed] = 1;
        if (used[tup1]) cnt1 = 0; else used[tup1] = 1;
        if (used[tup2]) cnt2 = 0; else used[tup2] = 1;
        if (used[tup3]) cnt3 = 0; else used[tup3] = 1;
        if (used[tup4]) cnt4 = 0; else used[tup4] = 1;

        if (cnt0 + cnt1 + cnt2 + cnt3 + cnt4 >= CUTOFF)
          { model = DecodeTuple(seed*16+i,5);
            dif   = Wrap_Around(seq,model,5,ispref);
	    MATCH(5,seed*16+i,dif)
#ifdef PRINT_MODELS
            printf("Cycle of length 5: %s",model);
            printf(" -> %d\n",dif);
#endif
          }

        used[seed] = 0;
        used[tup1] = 0;
        used[tup2] = 0;
        used[tup3] = 0;
        used[tup4] = 0;
      } 

  /* Check 6-models a1.a2.a3.i for all i in [0,63] */

  for (i = 0; i < 64; i++)

    /* Avoid re-exploring smaller models */

    if (i != seed && (period != 2 || i != avd6))
      { if (period == 1 && pre1 == (i >> 4)) continue;
        if (period == 2 && pre2 == (i & 0xF)) continue;

        tup1 = (suf2 | (i >> 4));
        tup2 = (suf1 | (i >> 2));
        tup3 = i;
        tup4 = (((i & 0xF) << 2) | pre1);
        tup5 = (((i & 0x3) << 4) | pre2);
        cnt1 = BUCKET(tup1);
        cnt2 = BUCKET(tup2);
        cnt3 = BUCKET(tup3);
        cnt4 = BUCKET(tup4);
        cnt5 = BUCKET(tup5);

        /* Don't use a previous seed & count repeated els only once */

        if (cnt1 < 0 || cnt2 < 0 || cnt3 < 0 || cnt4 < 0 || cnt5 < 0) continue;
        used[seed] = 1;
        if (used[tup1]) cnt1 = 0; else used[tup1] = 1;
        if (used[tup2]) cnt2 = 0; else used[tup2] = 1;
        if (used[tup3]) cnt3 = 0; else used[tup3] = 1;
        if (used[tup4]) cnt4 = 0; else used[tup4] = 1;
        if (used[tup5]) cnt5 = 0; else used[tup5] = 1;

        if (cnt0 + cnt1 + cnt2 + cnt3 + cnt4 + cnt5 >= CUTOFF)
          { model = DecodeTuple(seed*64+i,6);
            dif   = Wrap_Around(seq,model,6,ispref);
	    MATCH(6,seed*64+i,dif)
#ifdef PRINT_MODELS
            printf("Cycle of length 6: %s",model);
            printf(" -> %d\n",dif);
#endif
          }

        used[seed] = 0;
        used[tup1] = 0;
        used[tup2] = 0;
        used[tup3] = 0;
        used[tup4] = 0;
        used[tup5] = 0;
      } 

  match[matop] = -1;
  return (match);
}

/* Given sequence "seq" and array of frequencies of its KMERLEN-mers in
   "bucket", find all models that match the first (iff "ispref" != 0)
   or last TESTLEN of the sequence with not more than ERRORS errors and
   call "handler" with each model.
*/

static void Micro_Essential_Els(char *seq, int *bucket,
                                void (*handler)(char *), int ispref)
{ int ptop, ess;
  int next[KMERNUM];
  int heap[TESTMER+1];
  int perm[KMERNUM];
  int *match = NULL;

  /* Determine the sort permutation of all non-zero entries of bucket, i.e.
     compute perm[0..ptop-1] such that bucket[perm[i+1]] <= bucket[perm[i]]
     for all i < ptop-1, bucket[perm[i]] > 0 for all i < ptop, and
     ptop = | { i : bucket[i] > 0 }.                                      */

  { int i, c, p;

    for (i = 0; i <= TESTMER; i++)
      heap[i] = -1;
    for (i = 0; i < KMERNUM; i++)
      { c = bucket[i];
        if (c > 0)
          { next[i] = heap[c];
            heap[c] = i;
          }
      }
    p = 0;
    for (c = TESTMER; c >= 0; c--)
      for (i = heap[c]; i >= 0; i = next[i])
        perm[p++] = i;
    ptop = p;
  }

  /* Determine which of the first few buckets in the sorted order are
     such that their frequency plus that of the following 5 buckets
     (or those that remain) is not less than CUTOFF.  Any matching
     model must contain at least one of these "essential" buckets. */

  { int i, s, clen;

    s = 0;
    for (i = 0; i < MICROMX; i++)
      { if (i >= ptop) break;
        s += bucket[perm[i]];
      }
    clen = i;
    for (i = 0; s >= CUTOFF; i++)
      if (i+clen >= ptop)
        { s -= bucket[perm[i]];
          clen -= 1;
        }
      else
        s += bucket[perm[i+MICROMX]] - bucket[perm[i]];
    ess = i;
  }

  /* Explore each essential seed to see if there are models, not
     containing a preceeding seed as a substring, that match.     */

#ifdef PRINT_MODELS
  printf("Exploring %d alternates\n",ess);
#endif

  { int i;
    for (i = 0; i < ess; i++)
      { match = Explore(seq,bucket,perm[i],i == 0,ispref);
        bucket[perm[i]] = -bucket[perm[i]];
      }
  }

#ifdef PRINT_MODELS
    { int i;
      printf("\nTop 3-tuple counts\n");
      for (i = 0; i < ptop; i++)
        printf(" %s",DecodeTuple(perm[i],3));
      printf("\n");
      for (i = 0; i < ptop; i++)
        printf(" %3d",bucket[perm[i]]);
      printf("\n");
    }
#endif

  /* Call the user's handler on each matching model (if any) */

  { int i, len, code;

    for (i = 0; match[i] >= 0; i++)
      { code = match[i] & 0xFFF;
        len  = match[i] >> 12;
        handler(DecodeTuple(code,len));
      }
  }

  return;
}

/* Find micro-sats of length up to MICROMX(6) that matches a prefix
   (iff "ispref" != 0) or suffix of length TESTLEN(40) of "seq" with
   not more than ERRORS(6) errors and call "handler" with each one
   (if any).  Return the first level filter cutoff value (good for
   testing only)
*/

int MicroFinder_AS(char *seq, int ispref, void (*handler)(char *))
{ int x, bucket[KMERNUM];

  x = Micro_Counter_Screen(seq,bucket,ispref);
  if (x >= CUTOFF)
    Micro_Essential_Els(seq,bucket,handler,ispref);
  return (x);
}

/* -------------------- TEST DRIVERS -------------------------------------*/

#define TRIALS 100000

static int ModelHit;

void TestHand(char *model)
{ ModelHit = 1; 
#ifdef PRINT_MATCHES
  printf("'%s'\n        ",model);
#endif
}

#ifdef TESTRANDOM

int main(int argc, char *argv[])
{ static int DNA[] = { 'a', 'c', 'g', 't' };

  char testseq[TESTLEN+1];
  int  histo[TESTLEN+1];
  int  histo2[TESTLEN+1];

  unsigned short ranbits[3];

  { unsigned int pid;

    pid = getpid();
    ranbits[0] = pid >> 16;
    ranbits[1] = pid & 0xFF;
    ranbits[2] = 0x330E;
  }

  { int i;

    for (i = 0; i <= TESTLEN; i++)
      histo[i] = histo2[i] = 0;
  }

  { int i, j, x;

    testseq[TESTLEN] = '\0';
    for (i = 0; i < TRIALS; i++)
      { for (j = 0; j < TESTLEN; j++)
          testseq[j] = DNA[(int) (4.*erand48(ranbits))];
        ModelHit = 0;
        x = MicroFinder_AS(testseq,1,TestHand);
        histo[x] += 1;
        if (x >= CUTOFF && ModelHit)
          histo2[x] += 1;
      }
  }

  { int i, cum, cum2;

    printf("\nHistogram of scores\n\n");
    printf("Score:  Count    Total  %%ofTotal  Count    Total  %%ofTotal\n");
    for (i = TESTLEN; i > 0; i--)
      if (histo[i] > 0 || histo2[i] > 0) break;
    cum = cum2 = 0;
    while (i > 0)
      { cum += histo[i];
        cum2 += histo2[i];
        printf(" %2d:  %7d  %7d  %6.2f%%",i,histo[i],cum,(100.*cum)/TRIALS);
        printf("  %7d  %7d  %6.2f%%",histo2[i],cum2,(100.*cum2)/TRIALS);
        if ((TESTMER-i) % KMERLEN == 0)
          printf(" <- %5.2f%% Error Lower Bound",
                 (100.*(TESTMER-i))/(KMERLEN*TESTLEN));
        printf("\n");
        i -= 1;
      }
  }

  exit (0);
}

#endif

#ifdef TESTMICRO

#define MICLEN 6

int main(int argc, char *argv[])
{ static int DNA[] = { 'a', 'c', 'g', 't' };

  char testseq[TESTLEN+1];
  int  histo[TESTLEN+1];
  int  histo2[TESTLEN+1];
  int  bucket[KMERNUM];

  unsigned short ranbits[3];

  { unsigned int pid;

    pid = getpid();
    ranbits[0] = pid >> 16;
    ranbits[1] = pid & 0xFF;
    ranbits[2] = 0x330E;
  }

  { int i;

    for (i = 0; i <= TESTLEN; i++)
      histo[i] = histo2[i] = 0;
  }

  { int i, j, errors;
    char micro[MICLEN+1];

    testseq[TESTLEN] = '\0';
    micro[MICLEN] = '\0';
    for (i = 0; i < TRIALS; i++)
      { int c, x;
     
        for (j = 0; j < MICLEN; j++)
          micro[j] = DNA[(int) (4.*erand48(ranbits))];
        j = c = 0;
        errors = ERRORS;
        while (j < TESTLEN)
          { if (erand48(ranbits)*(TESTLEN-j) < errors)
              { errors -= 1;
                switch ((int) (3.*erand48(ranbits)))
                { case 0:      /* Deletion */
                    c += 1;
                    break;
                  case 1:      /* Insertion */
                    testseq[j++] = DNA[(int) (4.*erand48(ranbits))];
                    break;
                  case 2:      /* Substitution */
                    do
                      x = DNA[(int) (4.*erand48(ranbits))];
                    while (x == micro[c]);
                    testseq[j++] = x;
                    c += 1;
                    break;
                }
              }
            else
              testseq[j++] = micro[c++];
            if (c >= MICLEN) c = 0;
          }
#ifdef PRINT_MATCHES
        printf("\nSeqnc = '%s'\n",testseq);
        printf("Micro = '%s'\n",micro);
        printf("Match = ");
#endif
        ModelHit = 0;
        x = MicroFinder_AS(testseq,1,TestHand);
        histo[x] += 1;
        if (x >= CUTOFF && ModelHit)
          histo2[x] += 1;
        else
          { printf("Missed on trial %d [ %hx %hx %hx ]\n",i,
                   ranbits[0],ranbits[1],ranbits[2]);
          }

      }
  }

  { int i, cum, cum2;

    printf("\nHistogram of scores\n\n");
    printf("Score:  Count    Total  %%ofTotal  Count    Total  %%ofTotal\n");
    for (i = TESTLEN; i > 0; i--)
      if (histo[i] > 0 || histo2[i] > 0) break;
    cum = cum2 = 0;
    while (i > 0)
      { cum += histo[i];
        cum2 += histo2[i];
        printf(" %2d:  %7d  %7d  %6.2f%%",i,histo[i],cum,(100.*cum)/TRIALS);
        printf("  %7d  %7d  %6.2f%%",histo2[i],cum2,(100.*cum2)/TRIALS);
        if ((TESTMER-i) % KMERLEN == 0)
          printf(" <- %5.2f%% Error Lower Bound",
                 (100.*(TESTMER-i))/(KMERLEN*TESTLEN));
        printf("\n");
        i -= 1;
      }
  }

  exit (0);
}

#endif
