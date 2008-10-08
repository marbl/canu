
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

static const char *rcsid = "$Id: CA_ALN_local.c,v 1.9 2008-10-08 22:02:54 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "AS_global.h"
#include "AS_UTL_reverseComplement.h"

#include "CA_ALN_local.h"
#include "AS_ALN_aligners.h"

/* Conditional control of diagnostic and debug print outs */

#undef  TEST_TABLE
#undef  TEST_DIAG
#undef TEST_TRAP
#undef PRINT_TRAP
#undef  TEST_DPREACH
#undef  TEST_TRAPTRIM
#undef  STAT_DPREACH
#undef REPORT_DPREACH
#undef  REPORT_SIZES
#undef  WARN_MISSED_SEGMENT

#undef REUSE_TABLE_UNSAFE   /* if the user gaurantees that this is acceptable,
			       reuse the "table" built on last call of Find_Local_Segments()
			       iff the address of the A sequence is the same between
			       last call and this call; in general, this is risky
			       business: bogus results ensue if in fact only the
			       address and not the string is the same! */

/* Index-based match parameters */


/* Note, KMERLEN 5, MINMATCH 20, MAXERROR 2, KTHRESH 6 is reasonable for performing fragment against fragment comparisons; it permits relatively small segments to be found; but it will not give acceptable run time for large comparisons such as a BAC against a BAC etc.  So, ... */


// default values controlling sensitivity and performance:
#define LOCAL_TUNED_FOR_FRAGMENTS
#ifdef LOCAL_TUNED_FOR_FRAGMENTS
#define KMERLEN   5   /* Must be >= 1 */
#define MINMATCH 20
#define MAXERROR  2
#define KTHRESH   6   /*  MINMATCH - (KMERLEN-1) - KMERLEN*MAXERROR */

#else

#define LOCAL_TUNED_FOR_MODERATE_ALIGNS
#ifdef LOCAL_TUNED_FOR_MODERATE_ALIGNS
#define KMERLEN   8   /* Must be >= 1 */
#define MINMATCH 36
#define MAXERROR  2
#define KTHRESH  13   /*  MINMATCH - (KMERLEN-1) - KMERLEN*MAXERROR */
#else  // near-identity long matches
#define KMERLEN   12   /* Must be >= 1 */
#define MINMATCH 100
#define MAXERROR  2
#define KTHRESH  65   /*  MINMATCH - (KMERLEN-1) - KMERLEN*MAXERROR */
#endif

#endif


int kmerlen=KMERLEN;
int minmatch=MINMATCH;
int maxerror=MAXERROR;
int kthresh=KTHRESH;



/* D.P. extension alignment scoring */

#ifdef LOCAL_TUNED_FOR_FRAGMENTS
/* N.B.: the larger MAXIGAP, the longer subpart of a trapezoid (potential
   segment) can be missed; this occurs if (a) Align_Recursion starts
   in the subpart, with a bad region on both sides of the subpart;
   the forward pass extends across one of the low-quality regions
   into a long high-quality region but then the reverse pass doesn't
   get back across the low-quality region--or rather, it does, but the
   best score doesn't.

   Of course, the smaller MAXIGAP, the shorter the segments and the
   less the chance of extending alignments across gaps between trapezoids.

   Setting MAXIGAP small allows alignment extension to end relatively
   early, so that a short high-quality region bounded by low-quality
   regions can end up being its own segment.  And, it also has the
   effect that if a low-quality region is small enough to get through,
   then it doesn't take that long a high-quality segment to increase
   the score back up to a new maximum value.
*/

#define MAXIGAP  3
#else
#define MAXIGAP  5
#endif

/*amount to subtract from score for mismatch */
//#define DIFFCOST 14
// we can define it to 14 in order not to extend the alignment
// only at a high level of stringency
#define DIFFCOST 3

/*amount to add to score for match*/
#define SAMECOST 1

static int diffcost=DIFFCOST;
static int samecost=SAMECOST;

#undef VARIABLE_SCORE_SCHEME  /* defining this causes scoring to
				 be set so that nearly
				 all extensions terminate within the
				 user-defined error rate; this means
			         we don't completely miss a high-quality
			         segment that is part of a larger segment
			         of lower quality--e.g., we could find
			         more-conserved portions of a repeat
			         even if the repeat as a whole was lower
			         fidelity, due to varying selectional
			         pressure at different positions;
				 the down-side to this is that we won't
				 chain together perfect matches across
				 a few bases of lower-fidelity sequence.
			            Leaving it undefined uses the
			         #define values of DIFFCOST and SAMECOST
*/
#ifdef VARIABLE_SCORE_SCHEME
  #define ALTERNATE_PCNT
#else
  #define ALTERNATE_PCNT
#endif


/* Trapezoid merging padding */

#define DPADDING   2
int bpadding;


static int BLOCKCOST = DIFFCOST*MAXIGAP;
static int MATCHCOST = DIFFCOST+SAMECOST;
#ifndef ALTERNATE_PCNT
static double RMATCHCOST = DIFFCOST+1.;
#endif

/* Major data types */

  /* Hit Record: Description of indexed based seed-match */

typedef struct {
  int diagonal;   /* Diagonal of hit */
  int bstart;     /* B position of start of hit */
  int bfinish;    /* B position of end of hit */
} HitRecord;

  /* Trapezoid Record: Description of trapezoidal match zone */

typedef struct _Trap_Tag {
  struct _Trap_Tag *next;   /* Organized in a list linked on this field */
  int top, bot;             /* B-coords of top and bottom of trapzoidal zone */
  int lft, rgt;             /* Left and right diagonals of trapzoidal zone */
} Trapezoid;

/*** INDEX CONSTRUCTION AND APPLICATION TO FILTERING ***/

/* Shared index and filter arrays used in this subsection */

typedef struct {
  int minim;
  int maxim;
  int count;
} DiagRecord;

static int  Kmask = -1;
static int *Table;          /* [0..Kmask+1] */
static int *Tuples = NULL;  /* [0..<Seqlen>-kmerlen] */
static int  Map[128];

static DiagRecord *DiagVec; /* [-(Alen-kmerlen)..(Blen-kmerlen) + maxerror] */

/* Reverse complement sequences -- so we do not recompute them over and over */
static char *ArevC,*BrevC;


/* Build index table for sequence S of length Slen. */

static void TableBuild(char *S, int Slen)
{ int   i, c;
  int   x, h;
  char *s;

  s = S+(kmerlen-1);

  for (c = 0; c <= Kmask; c++)
    Table[c] = 0;

  h = -kmerlen;
  c = 0;
  for (i = 0; i < kmerlen-1; i++)
    { x = Map[(int) (S[i])];
      if (x >= 0)
        c = (c << 2) | x;
      else
        { c <<= 2; h = i-(kmerlen-1); }
    }
  for (i = 0; i <= Slen-kmerlen; i++)
    { x = Map[(int) (s[i])];
      if (x >= 0)
        c = ((c << 2) | x) & Kmask;
      else
        { c = (c << 2) & Kmask; h = i; }
      if (i >= h+kmerlen)
        Table[c+1] += 1;
    }

  for (c = 2; c <= Kmask; c++)
    Table[c] += Table[c-1];

  h = -kmerlen;
  c = 0;
  for (i = 0; i < kmerlen-1; i++)
    { x = Map[(int) (S[i])];
      if (x >= 0)
        c = (c << 2) | x;
      else
        { c <<= 2; h = i-(kmerlen-1); }
    }
  for (i = 0; i <= Slen-kmerlen; i++)
    { x = Map[(int) (s[i])];
      if (x >= 0)
        c = ((c << 2) | x) & Kmask;
      else
        { c = (c << 2) & Kmask; h = i; }
      if (i >= h+kmerlen)
        Tuples[Table[c]++] = i;
    }

  for (c = Kmask; c >= 0; c--)
    Table[c+1] = Table[c];
  Table[0] = 0;

#ifdef TEST_TABLE
  { int i, c, j;

    static char Convert[] = { 'a', 'c', 'g', 't' };

    for (c = 0; c <= Kmask; c++)
      { fprintf(stderr, "Table[%d = ",c);
        for (i = kmerlen-1; i >= 0; i--)
          fprintf(stderr, "%c",Convert[c>>(i*2) & 0x3]);
        fprintf(stderr, "]\n");
        for (i = Table[c]; i < Table[c+1]; i++)
          fprintf(stderr, "  %d (%.*s)\n",Tuples[i],kmerlen,S+Tuples[i]);
      }
  }
#endif
}

/* Apply index to find filtered hits between sequences, returning pointer to
   array of HitRecords of length in the integer pointed at by Hitlen       */

static int HSORT(const void *l, const void *r)
{ HitRecord *x, *y;
  x = (HitRecord *) l;
  y = (HitRecord *) r;
  return (x->bstart - y->bstart);
}

static HitRecord *Find_Hits(char *A, int Alen, char *B, int Blen, int *Hitlen)
{ static int        HitMax = -1;
  static HitRecord *HitList;
  int hits, disconnect;
#ifdef REPORT_SIZES
  int sum;
#endif

  if (HitMax < 0)
    { HitMax = 10000;
      HitList = (HitRecord *) safe_malloc(sizeof(HitRecord)*HitMax);
    }

  { int i, j, c;
    int x, h;
    char *b;

    for (j = -Alen; j <= Blen+maxerror; j++)
      { DiagRecord *dp;
        dp = DiagVec + j;
        dp->count = dp->maxim = 0;
      }

    hits = 0;
    disconnect = minmatch - kmerlen;
#ifdef REPORT_SIZES
    sum = 0;
#endif
    h = -kmerlen;
    c = 0;
    for (i = 0; i < kmerlen-1; i++)
      { x = Map[(int) (B[i])];
        if (x >= 0)
          c = (c << 2) | x;
        else
          { c <<= 2; h = i-(kmerlen-1); }
      }
    b = B + (kmerlen-1);
    for (i = 0; i <= Blen-kmerlen; i++)
      { x = Map[(int) (b[i])];
        if (x >= 0)
          c = ((c << 2) | x) & Kmask;
        else
          { c = (c << 2) & Kmask; h = i; }
        if (i >= h+kmerlen)
          for (j = Table[c]; j < Table[c+1]; j++)
            { DiagRecord *dp;
              int e, k;
              k  = i-Tuples[j];
              dp = DiagVec + k;
              for (e = 0; e <= maxerror; e++)
                { if (dp->maxim < i-disconnect)
                    { if (dp->count >= kthresh)
                        { HitRecord *hp;
                          if (hits >= HitMax)
                            { HitMax = (int)(1.2*hits) + 5000;
                              HitList = (HitRecord *) safe_realloc(HitList, sizeof(HitRecord)*HitMax);
                            }
                          hp = HitList + hits;
                          hp->diagonal = k;
                          hp->bstart   = dp->minim;
                          hp->bfinish  = dp->maxim + kmerlen;
                          hits += 1;
                        }
                      dp->count = 0;
                    }
                  if (dp->count == 0)
                    dp->minim = i;
                  dp->count += 1;
                  dp->maxim = i;
                  dp += 1;
                }
            }
#ifdef REPORT_SIZES
        sum += Table[c+1] - Table[c];
#endif
      }

    for (j = -Alen; j <= Blen+maxerror; j++)
      { DiagRecord *dp;
        dp = DiagVec + j;
        if (dp->count >= kthresh)
          { HitRecord *hp;
            if (hits >= HitMax)
              { HitMax = (int)(1.2*hits) + 5000;
                HitList = (HitRecord *)safe_realloc(HitList,sizeof(HitRecord)*HitMax);
              }
            hp = HitList + hits;
            hp->diagonal = j;
            hp->bstart   = dp->minim;
            hp->bfinish  = dp->maxim + kmerlen;
            hits += 1;
          }
      }
  }

  qsort(HitList,hits,sizeof(HitRecord),HSORT);

#ifdef REPORT_SIZES
  fprintf(stderr, "\n  %9d %d-mers (%f%% of matrix)\n",sum,kmerlen,(100.*sum/Alen)/Blen);
  fprintf(stderr, "  %9d seed hits (%f%% of matrix)\n",hits,(100.*hits/Alen)/Blen);
#endif

#ifdef TEST_DIAG
  for (i = 0; i < hits; i++)
    fprintf(stderr, "Diag %d [%d,%d]\n",
           HitList[i].diagonal,HitList[i].bstart,HitList[i].bfinish);
#endif

  *Hitlen = hits;
  return (HitList);
}


/*** FORWARD AND REVERSE D.P. EXTENSION ROUTINES ***/
/*      Called at the mid-point of trapezoid -- mid X [lo,hi], the extension
        is computed to an end point and the lowest and highest diagonals
        are recorded.  These are returned in a partially filled Local_Segment
        record, that will be merged with that returned for extension in the
        opposite direction.
*/

static Local_Segment *TraceForwardPath(char *A, int Alen, char *B, int Blen,
                                       int mid, int lo, int hi)
{ static Local_Segment rez;
  int *V;
  int  mxv, mxl, mxr, mxi, mxj;
  int  i, j;
  int *Base1, *Base2;
#ifdef STAT_DPREACH
  long dparea;
#endif

  Base1 = ((int *) DiagVec);
  Base2 = Base1 + (Blen+1);

  /* Set basis from (mid,lo) .. (mid,hi) */

  V = Base1;
  if (lo < 0)    lo = 0;
  if (hi > Blen) hi = Blen;

  for (j = lo; j <= hi; j++)
    V[j] = 0;
  hi += MAXIGAP;
  if (hi > Blen) hi = Blen;
  for (; j <= hi; j++)
    V[j] = V[j-1] - diffcost;
  mxv = 0;
  mxr = mid - lo;
  mxl = mid - hi;
  mxi = mid;
  mxj = lo;

#ifdef STAT_DPREACH
  dparea = (hi-lo+1);
#endif

  /* Advance to next row */

  for (i = mid; lo <= hi && i < Alen; i++)
    { int  c, v;
      int *W;

      W = V;
      if (V == Base1)
        V = Base2;
      else
        V = Base1;

#ifdef TEST_DPREACH
      fprintf(stderr, "\n%d [%d,%d](%d,%d) x = %d\n        ",i,lo,hi,mxi,mxj,mxv);
      for (j = lo; j <= hi; j++)
        fprintf(stderr, "  %c  ",B[j]);
      fprintf(stderr, "\n ");
      for (j = lo; j <= hi; j++)
        fprintf(stderr, " %4d",W[j]);
      fprintf(stderr, "\n%c",A[i]);
#endif

      v = W[lo];
      c = V[lo] = v - diffcost;
#ifdef TEST_DPREACH
      fprintf(stderr, " %4d",c);
#endif
      for (j = lo+1; j <= hi; j++)
        { int r, t;

          t = c;
          c = v;
          v = W[j];
          if (Map[(int)A[i]] == Map[(int)B[j-1]] && Map[(int) (A[i])] >= 0)
	    c += MATCHCOST;

          r = c;
          if (v > r) r = v;
          if (t > r) r = t;

          V[j] = c = r - diffcost;
          if (c >= mxv)
            { mxv = c;
              mxi = i+1;
              mxj = j;
	      //fprintf(stderr, "reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
            }
#ifdef TEST_DPREACH
          fprintf(stderr, " %4d",c);
#endif
        }

      if (j <= Blen)
        { int r;

          if (Map[(int)A[i]] == Map[(int)B[j-1]] && Map[(int) (A[i])] >= 0)
	    v += MATCHCOST;

          r = v;
          if (c > r) r = c;

          V[j] = v = r - diffcost;
          if (v > mxv)
            { mxv = v;
              mxi = i+1;
              mxj = j;
	      //fprintf(stderr, "reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
            }
#ifdef TEST_DPREACH
          fprintf(stderr, " %4d",v);
#endif

          for (j++; j <= Blen; j++)
            { v -= diffcost;
              if (v < mxv - BLOCKCOST) break;
              V[j] = v;
#ifdef TEST_DPREACH
              fprintf(stderr, " %4d",v);
#endif
            }
        }
#ifdef TEST_DPREACH
      fprintf(stderr, "\n");
#endif

      hi = j-1;

      while (lo <= hi && V[lo] < mxv - BLOCKCOST)
        lo += 1;
      while (lo <= hi && V[hi] < mxv - BLOCKCOST)
        hi -= 1;

      if ((i+1) - lo > mxr)
        mxr = (i+1) - lo;
      if ((i+1) - hi < mxl)
        mxl = (i+1) - hi;

#ifdef STAT_DPREACH
      dparea += (hi-lo+1);
#endif
    }

#ifdef STAT_DPREACH
   fprintf(stderr, "  DP_Area = %ld  Peak is %d @ (%d,%d) in [%d,%d]\n",
          dparea,mxv,mxi,mxj,mxl,mxr);
#endif

  rez.aepos = mxj;
  rez.bepos = mxi;
  rez.ldiag = mxl;
  rez.hdiag = mxr;
  rez.score = mxv;
  return (&rez);
}

static Local_Segment *TraceReversePath(char *A, int Alen, char *B, int Blen,
                                       int top, int lo, int hi, int bot,
                                       int xfactor)
{ static Local_Segment rez;
  int *V;
  int  mxv, mxl, mxr, mxi, mxj;
  int  i, j;
  int *Base1, *Base2;
#ifdef STAT_DPREACH
  long dparea;
#endif

  Base1 = ((int *) DiagVec);
  Base2 = Base1 + (Blen+1);

  /* Set basis from (top,lo) .. (top,hi) */

  V = Base1;
  if (lo < 0)    lo = 0;
  if (hi > Blen) hi = Blen;

  for (j = hi; j >= lo; j--)
    V[j] = 0;
  lo -= MAXIGAP;
  if (lo < 0) lo = 0;
  for (; j >= lo; j--)
    V[j] = V[j+1] - diffcost;
  mxv = 0;
  mxr = top - lo;
  mxl = top - hi;
  mxi = top;
  mxj = lo;

#ifdef STAT_DPREACH
  dparea = (hi-lo+1);
#endif

  /* Advance to next row */

  if (top-1 <= bot) xfactor = BLOCKCOST;

  for (i = top-1; lo <= hi && i >= 0; i--)
    { int  c, v;
      int *W;

      W = V;
      if (V == Base1)
        V = Base2;
      else
        V = Base1;

#ifdef TEST_DPREACH
      fprintf(stderr, "\n%d [%d,%d](%d,%d) x = %d\n        ",i,lo,hi,mxi,mxj,mxv);
      for (j = hi; j >= lo; j--)
        fprintf(stderr, "  %c  ",B[j-1]);
      fprintf(stderr, "\n ");
      for (j = hi; j >= lo; j--)
        fprintf(stderr, " %4d",W[j]);
      fprintf(stderr, "\n%c",A[i]);
#endif

      v = W[hi];
      c = V[hi] = v - diffcost;
#ifdef TEST_DPREACH
      fprintf(stderr, " %4d",c);
#endif
      for (j = hi-1; j >= lo; j--)
        { int r, t;

          t = c;
          c = v;
          v = W[j];
          if (Map[(int)A[i]] == Map[(int)B[j]] && Map[(int) (A[i])] >= 0)
	    c += MATCHCOST;

          r = c;
          if (v > r) r = v;
          if (t > r) r = t;

          V[j] = c = r - diffcost;
          if (c >= mxv)
            { mxv = c;
              mxi = i;
              mxj = j;
	      //fprintf(stderr, "reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
            }
#ifdef TEST_DPREACH
          fprintf(stderr, " %4d",c);
#endif
        }

      if (j >= 0)
        { int r;

          if (Map[(int)A[i]] == Map[(int)B[j]] && Map[(int) (A[i])] >= 0)
	    v += MATCHCOST;

          r = v;
          if (c > r) r = c;

          V[j] = v = r - diffcost;
          if (v > mxv)
            { mxv = v;
              mxi = i;
              mxj = j;
	      //fprintf(stderr, "reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
            }
#ifdef TEST_DPREACH
          fprintf(stderr, " %4d",v);
#endif

          for (j--; j >= 0; j--)
            { v -= diffcost;
              if (v < mxv - xfactor) break;
              V[j] = v;
#ifdef TEST_DPREACH
              fprintf(stderr, " %4d",v);
#endif
            }
        }
#ifdef TEST_DPREACH
      fprintf(stderr, "\n");
#endif

      lo = j+1;

      while (lo <= hi && V[lo] < mxv - xfactor)
        lo += 1;
      while (lo <= hi && V[hi] < mxv - xfactor)
        hi -= 1;

      if (i == bot) xfactor = BLOCKCOST;

      if (i-lo > mxr)
        mxr = i-lo;
      if (i-hi < mxl)
        mxl = i-hi;

#ifdef STAT_DPREACH
      dparea += (hi-lo+1);
#endif
    }

#ifdef STAT_DPREACH
   fprintf(stderr, "  DP_Area = %ld  Peak is %d @ (%d,%d) in [%d,%d]\n",
          dparea,mxv,mxi,mxj,mxl,mxr);
#endif

  rez.abpos = mxj;
  rez.bbpos = mxi;
  rez.ldiag = mxl;
  rez.hdiag = mxr;
  rez.score = mxv;
  return (&rez);
}


/*** MERGING INDEX HITS INTO TRAPEZOIDAL ZONES ***/

static Trapezoid *Build_Trapezoids(char *A, int Alen, char *B, int Blen,
                                   HitRecord *list, int Hitlen, int *Traplen)
{ static Trapezoid  *free = NULL;

  Trapezoid *traporder, *traplist, *tailend;
  Trapezoid *b, *f, *t;
  int i, inserted;
  int trapcount, traparea;

  trapcount = 0;
  traparea  = 0;
  traporder = NULL;
  traplist  = NULL;
  for (i = 0; i < Hitlen; i++)
    { inserted = 0;
#ifdef TEST_TRAP
      fprintf(stderr, "  Diag %d [%d,%d]\n",
             list[i].diagonal,list[i].bstart,list[i].bfinish);
#endif
      f = NULL;
      for (b = traporder; b != NULL; b = t)
        { t = b->next;
          if (b->top < list[i].bstart - bpadding)
            { trapcount += 1;
              traparea  += (b->top - b->bot + 1) * (b->rgt - b->lft + 1);
              if (f == NULL)
                traporder = t;
              else
                f->next = t;
              b->next = traplist;
              traplist = b;
#ifdef PRINT_TRAP
	      fprintf(stderr, "pre-trim trapezoid: [%d,%d]X[%d,%d]\n",
		     b->bot,b->top,b->lft,b->rgt);
#endif
            }
          else if (list[i].diagonal > b->rgt + DPADDING)
            f = b;
          else if (list[i].diagonal >= b->lft - DPADDING)
            { if (list[i].diagonal < b->lft)
                b->lft = list[i].diagonal;
              if (list[i].diagonal > b->rgt)
                b->rgt = list[i].diagonal;
              if (list[i].bfinish > b->top)
                b->top = list[i].bfinish;

              if (f != NULL && f->rgt + DPADDING >= b->lft)
                { f->rgt = b->rgt;
                  if (f->bot > b->bot) f->bot = b->bot;
                  if (f->top < b->top) f->top = b->top;
                  f->next = t;
                  b->next = free;
                  free = b;
                }
              else if (t != NULL && t->lft - DPADDING <= b->rgt)
                { b->rgt = t->rgt;
                  if (b->bot > t->bot) b->bot = t->bot;
                  if (b->top < t->top) b->top = t->top;
                  b->next = t->next;
                  t->next = free;
                  free = t;
                  t = b->next;
                  f = b;
                }
              else
                f = b;
              inserted = 1;
            }
          else if (! inserted)
            { if (free == NULL)
                { free = (Trapezoid  *)safe_malloc(sizeof(Trapezoid));
                  free->next = NULL;
                }
              if (f == NULL)
                f = traporder = free;
              else
                f = f->next = free;
              free = f->next;
              f->next = b;
              f->top = list[i].bfinish;
              f->bot = list[i].bstart;
              f->lft = f->rgt = list[i].diagonal;
              f = b;
              inserted = 1;
            }
          else
            f = b;
        }
      if (! inserted)
        { if (free == NULL)
            { free = (Trapezoid  *)safe_malloc(sizeof(Trapezoid));
              free->next = NULL;
            }
          if (f == NULL)
            f = traporder = free;
          else
            f = f->next = free;
          free = f->next;
          f->next = b;
          f->top = list[i].bfinish;
          f->bot = list[i].bstart;
          f->lft = f->rgt = list[i].diagonal;
        }
#ifdef TEST_TRAP
      fprintf(stderr, "  Blist:");
      for (b = traporder; b != NULL; b = b->next)
        fprintf(stderr, " [%d,%d]x[%d,%d]",
               b->bot,b->top,b->lft,b->rgt);
      fprintf(stderr, "\n");
#endif
    }

  for (b = traporder; b != NULL; b = t)
    { t = b->next;
      trapcount += 1;
      traparea  += (b->top - b->bot + 1) * (b->rgt - b->lft + 1);
      b->next  = traplist;
      traplist     = b;
    }

#ifdef REPORT_SIZES
  fprintf(stderr, "\n  %9d trapezoids of area %d (%f%% of matrix)\n",
         trapcount,traparea,(100.*trapcount/Alen)/Blen);
#endif

  { int lag, lst, lclip;
    int abot, atop;

#ifdef TEST_TRAPTRIM
    fprintf(stderr, "B trimming:\n");
#endif
    for (b = traplist; b != NULL; b = b->next)
      { lag = (b->bot-MAXIGAP)+1;
        if (lag < 0) lag = 0;
        lst = b->top+MAXIGAP;
        if (lst > Blen) lst = Blen;

#ifdef TEST_TRAPTRIM
        fprintf(stderr, "   [%d,%d]x[%d,%d] = %d\n",
               b->bot,b->top,b->lft,b->rgt,b->top - b->bot + 1);
#endif

        for (i = lag; i < lst; i++)
          { if (Map[(int) (B[i])] >= 0)
              { if (i-lag >= MAXIGAP)
                  { if (lag - b->bot > 0)
                      { if (free == NULL)
                          { free = (Trapezoid  *)safe_malloc(sizeof(Trapezoid));
                            free->next = NULL;
                          }
                        t = free->next;
                        *free = *b;
                        b->next = free;
                        free = t;
                        b->top = lag;
                        b = b->next;
                        b->bot = i;
                        trapcount += 1;
                      }
                    else
                      b->bot = i;
#ifdef TEST_TRAPTRIM
                    fprintf(stderr, "  Cut trap B[%d,%d]\n",lag,i);
#endif
                  }
                lag = i+1;
              }
          }
        if (i-lag >= MAXIGAP)
          b->top = lag;
      }

#ifdef TEST_TRAPTRIM
    fprintf(stderr, "A trimming:\n");
#endif
    tailend = NULL;
    for (b = traplist; b != NULL; b = b->next)
      { if (b->top - b->bot < kmerlen) continue;

        abot = b->bot - b->rgt;
        atop = b->top - b->lft;

#ifdef TEST_TRAPTRIM
        fprintf(stderr, "   [%d,%d]x[%d,%d] = %d\n",
               b->bot,b->top,b->lft,b->rgt,b->top - b->bot + 1);
#endif

        lag = (abot - MAXIGAP) + 1;
        if (lag < 0) lag = 0;
        lst = atop + MAXIGAP;
        if (lst > Alen) lst = Alen;

        lclip = abot;
        for (i = lag; i < lst; i++)
          { if (Map[(int) (A[i])] >= 0)
              { if (i-lag >= MAXIGAP)
                  { if (lag > lclip)
                      { if (free == NULL)
                          { free = (Trapezoid  *)safe_malloc(sizeof(Trapezoid));
                            free->next = NULL;
                          }
                        t = free->next;
                        *free = *b;
                        b->next = free;
                        free = t;

#ifdef TEST_TRAPTRIM
                        fprintf(stderr, "     Clip to %d,%d\n",lclip,lag);
#endif
                        { int x, m;
                          x = lclip + b->lft;
                          if (b->bot < x)
                            b->bot = x;
                          x = lag + b->rgt;
                          if (b->top > x)
                            b->top = x;
                          m = (b->bot + b->top) / 2;
                          x = m - lag;
                          if (b->lft < x)
                            b->lft = x;
                          x = m - lclip;
                          if (b->rgt > x)
                            b->rgt = x;
#ifdef TEST_TRAPTRIM
                          fprintf(stderr, "        [%d,%d]x[%d,%d] = %d\n",
                                 b->bot,b->top,b->lft,b->rgt,b->top-b->bot+1);
#endif
                        }

                        b = b->next;
                        trapcount += 1;
                      }
                    lclip = i;
                  }
                lag = i+1;
              }
          }

        if (i-lag < MAXIGAP)
          lag = atop;

#ifdef TEST_TRAPTRIM
        fprintf(stderr, "     Clip to %d,%d\n",lclip,lag);
#endif
        { int x, m;
          x = lclip + b->lft;
          if (b->bot < x)
            b->bot = x;
          x = lag + b->rgt;
          if (b->top > x)
            b->top = x;
          m = (b->bot + b->top) / 2;
          x = m - lag;
          if (b->lft < x)
            b->lft = x;
          x = m - lclip;
          if (b->rgt > x)
            b->rgt = x;
#ifdef TEST_TRAPTRIM
          fprintf(stderr, "        [%d,%d]x[%d,%d] = %d\n",
                 b->bot,b->top,b->lft,b->rgt,b->top-b->bot+1);
#endif
        }

        tailend = b;
      }
  }

  if (tailend != NULL)
    { tailend->next = free;
      free = traplist;
    }

#ifdef REPORT_SIZES
  fprintf(stderr, "  %9d trimmed trap.s of area %d (%f%% of matrix)\n",
         trapcount,traparea,(100.*trapcount/Alen)/Blen);
#endif

  *Traplen = trapcount;
#ifdef PRINT_TRAP
  for(b=traplist;b!=NULL;b=b->next)
    fprintf(stderr, "Output trapezoid: [%d,%d]X[%d,%d]\n",b->bot,b->top,b->lft,b->rgt);
#endif
  return (traplist);
}


/*** FINDING ALIGNMENTS WITHIN A TRAPEZOIDAL ZONE ***/

static int TSORT(const void *l, const void *r)
{ Trapezoid *x, *y;
  x = *((Trapezoid **) l);
  y = *((Trapezoid **) r);
  return (x->bot - y->bot);
}

static int StSORT(const void *l, const void *r)
{ Local_Segment *x, *y;
  x = (Local_Segment *) l;
  y = (Local_Segment *) r;
  if (x->abpos < y->abpos)
    return (-1);
  else if (x->abpos > y->abpos)
    return (1);
  else
    return (x->bbpos - y->bbpos);
}

static int FnSORT(const void *l, const void *r)
{ Local_Segment *x, *y;
  x = (Local_Segment *) l;
  y = (Local_Segment *) r;
  if (x->aepos < y->aepos)
    return (-1);
  else if (x->aepos > y->aepos)
    return (1);
  else
    return (x->bepos - y->bepos);
}

static Trapezoid **Tarray = NULL;
static int        *Covered;
static Local_Segment *SegSols = NULL;
static int            SegMax = -1;
static int            NumSegs;

#ifdef REPORT_DPREACH
static int  Al_depth;
#endif

static void Align_Recursion(char *A, int Alen, char *B, int Blen,
                            Trapezoid *b, int current, int comp,
                            int MinLen, float MaxDiff, int Traplen)
{ int j, mid, indel;
  float pcnt;
  Local_Segment *hend, *lend;
  Trapezoid ltrp, htrp;

#undef START_AT_BEGINNING_OF_TRAP
#ifdef START_AT_BEGINNING_OF_TRAP
  mid = b->bot;
#else
  mid = (b->bot + b->top) / 2;
#endif

#ifdef REPORT_DPREACH
  fprintf(stderr, " [%d,%d]x[%d,%d] = %d (Depth = %d)\n",
         b->bot,b->top,b->lft,b->rgt,b->top - b->bot + 1,Al_depth);
#endif


  lend = TraceForwardPath(B,Blen,A,Alen,mid,mid-b->rgt,mid-b->lft);

  { int x;

    x = 0;
    do
      { x += 1;

      //fprintf(stderr, "Trying reverse pass\n");
        hend = TraceReversePath(B,Blen,A,Alen,
                                lend->bepos,lend->aepos,lend->aepos,
                                mid+MAXIGAP,BLOCKCOST+2*x*diffcost);
	//fprintf(stderr, "End reverse pass\n");
      }
    while (hend->bbpos > mid + x*MAXIGAP && hend->score < lend->score);

  hend->aepos = lend->aepos;
  hend->bepos = lend->bepos;


    /* We can miss a small segment here!

	the segment is at the beginning of a trapezoid;
	it is followed by a run of bad luck which is
		- not bad enough to terminate an extension but
		- long enough to have a negative score with abs. value
		  greater than the positive value of the segment that
	          will be missed
	after the run of bad luck is a larger good run which *does*
	   exceed the bad run, so that the best value for the
	  forward extension goes past the bad run

	  What happens is that when we trace backwards,
	  the best value occurs after the bad segment.

	  Thus, even if we start the search before the bad run,
	  the returned segment starts after the bad run.

	  I guess basically this means that we can miss a segment
	  if it has a positive value smaller than BLOCKCOST.

	  So, if we want small minimum segments,
	  could we lower BLOCKCOST?
	  This seems not to work; for instance,
	  with a scoring scheme of 1:10,
	  but a desire to find a segment consisting of, e.g.,
	  10 matches with one mismatch in the middle (i.e. score = 1),
	  we'd have to have BLOCKCOST = 0---not a good idea!

	  So, how about testing for the possible case and running the
	  search backwards when it occurs?


	  New case on which this same logic is attempted: if we are trying
	  a recursive alignment (based on ltrp or htrp) but the
	  TraceForwardPath step got nowhere, then try in reverse orientation.

    */

    if(hend->bbpos > mid+x*MAXIGAP || hend->bepos==mid )
      {
#ifdef WARN_MISSED_SEGMENT
	fprintf(stderr, "WARNING: might have missed a small local segment (possible segment with score < blockcost)!\n");
#endif

#define CHECK_FOR_MISSING_SEGMENT
#ifdef  CHECK_FOR_MISSING_SEGMENT

#ifdef WARN_MISSED_SEGMENT
	fprintf(stderr, "WARNING: will try to reverse direction of search!\n");
#endif

	/* Need to:
	   reverse both sequences
	   reverse mid, top bottom left and right
	   run forward
	   run backward until converged
	   reverse resulting segment
	   reverse both sequences
	*/

	mid=Blen-mid-1;
	{ int tmp;
  	  tmp=b->top;
	  b->top=Blen-b->bot-1;
	  b->bot=Blen-tmp-1;
	  tmp=b->rgt;
	  b->rgt=Blen-Alen-b->lft;
	  b->lft=Blen-Alen-tmp;
	}

	lend = TraceForwardPath(BrevC,Blen,ArevC,Alen,mid,mid-b->rgt,mid-b->lft);

	{ int x;

	  x = 0;
  	  do
	    { x += 1;
	      hend = TraceReversePath(BrevC,Blen,ArevC,Alen,
				      lend->bepos,lend->aepos,lend->aepos,
				      mid+MAXIGAP,BLOCKCOST+2*x*diffcost);
	    }
	  while (hend->bbpos > mid + x*MAXIGAP && hend->score < lend->score);

	  hend->aepos = lend->aepos;
	  hend->bepos = lend->bepos;

	}

	mid=Blen-mid-1;
	{ int tmp;

	  tmp=b->top;
	  b->top=Blen-b->bot-1;
	  b->bot=Blen-tmp-1;
	  tmp=b->rgt;
	  b->rgt=Blen-Alen-b->lft;
	  b->lft=Blen-Alen-tmp;

	  tmp=hend->ldiag;
	  hend->ldiag=Blen-Alen-hend->hdiag;
	  hend->hdiag=Blen-Alen-tmp;

	  // indices: start, end = positions in between bases,
	  //   so reversing is newpos=len-oldpos
	  tmp=hend->abpos;
	  hend->abpos=Alen-hend->aepos;
	  hend->aepos=Alen-tmp;
	  tmp=hend->bbpos;
	  hend->bbpos=Blen-hend->bepos;
	  hend->bepos=Blen-tmp;

	}


#endif  /* CHECK_FOR_MISSING_SEGMENT */

      }


  }





  ltrp = htrp = *b;
  ltrp.top = MIN(b->top,hend->bbpos) - MAXIGAP;


  htrp.bot = MAX(b->bot,hend->bepos) + MAXIGAP;

  if (hend->bepos - hend->bbpos >= MinLen &&
      hend->aepos - hend->abpos >= MinLen   )

    { indel = abs( (hend->abpos - hend->bbpos)
                 - (hend->aepos - hend->bepos) );

    /* original formula for pcnt doesn't seem to be robust to scoring scheme variation; use ALTERNATE_PCNT until Gene gets a fix in */

#ifndef ALTERNATE_PCNT
      pcnt = (1/RMATCHCOST)
           - (hend->score - indel)
           / (RMATCHCOST*(hend->bepos - hend->bbpos));
#else
      pcnt = (-hend->score+samecost*(hend->bepos-hend->bbpos))*1./
	(1.*(MATCHCOST)*(hend->bepos-hend->bbpos));
#endif

      if (pcnt <= MaxDiff) {
	hend->error = pcnt;

          for (j = current+1; j < Traplen; j++)
            { Trapezoid *t;
              int   ta, tb, ua, ub;

              t = Tarray[j];
              if (t->bot >= hend->bepos) break;

              tb = t->top - t->bot + 1;
              ta = t->rgt - t->lft + 1;
              if (t->lft < hend->ldiag)
                ua = hend->ldiag;
              else
                ua = t->lft;
              if (t->rgt > hend->hdiag)
                ub = hend->hdiag;
              else
                ub = t->rgt;

              if (ua > ub) continue;

              ua = ub - ua + 1;
              if (t->top > hend->bepos)
                ub = hend->bepos - t->bot + 1;
              else
                ub = tb;

              if (((1.*ua)/ta)*((1.*ub)/tb) > .99)
                Covered[j] = 1;
            }

          if (NumSegs >= SegMax)
            { SegMax = (int)(1.2*NumSegs) + 500;
              SegSols = (Local_Segment *) safe_realloc(SegSols, sizeof(Local_Segment)*SegMax);
            }

          { int d;

            d = hend->hdiag;  /*  Oops, diags to this point are b-a, not a-b. */
            hend->hdiag = - (hend->ldiag);
            hend->ldiag = - d;
            if (comp)
              { hend->bbpos = Blen - hend->bbpos;
                hend->bepos = Blen - hend->bepos;
                hend->ldiag = Blen + hend->ldiag;
                hend->hdiag = Blen + hend->hdiag;
              }
          }

          SegSols[NumSegs++] = *hend;

#ifdef REPORT_DPREACH
          fprintf(stderr, "  Hit from (%d,%d) to (%d,%d) within [%d,%d] score %d\n",
                 hend->abpos,hend->bbpos,hend->aepos,hend->bepos,
                 hend->ldiag,hend->hdiag,hend->score);
#endif
      }else{ // SAK
#ifdef REPORT_DPREACH
          fprintf(stderr, "  FAILED (%g > %g) Hit from (%d,%d) to (%d,%d) within [%d,%d] score %d\n",
		 pcnt, MaxDiff,
                 hend->abpos,hend->bbpos,hend->aepos,hend->bepos,
                 hend->ldiag,hend->hdiag,hend->score);
#endif

      }
    }

#ifdef REPORT_DPREACH
  Al_depth += 1;
#endif
  if (ltrp.top - ltrp.bot > MinLen && ltrp.top < b->top - MAXIGAP){
    Align_Recursion(A,Alen,B,Blen,&ltrp,current,comp,MinLen,MaxDiff,Traplen);
  }
  if (htrp.top - htrp.bot > MinLen){
    Align_Recursion(A,Alen,B,Blen,&htrp,current,comp,MinLen,MaxDiff,Traplen);
  }
#ifdef REPORT_DPREACH
  Al_depth -= 1;
#endif
}

static Local_Segment *Align_Trapezoids(char *A, int Alen, char *B, int Blen,
                                       Trapezoid *Traplist, int Traplen,
                                       int start, int comp,
                                       int MinLen, float MaxDiff, int *Seglen)
{ static int fseg;
  static int TarMax = -1;

  Trapezoid *b;
  int i;

  if (Traplen >= TarMax)
    { TarMax = (int)(1.2*Traplen) + 500;
      Tarray = (Trapezoid **) safe_realloc(Tarray,(sizeof(Trapezoid *) + sizeof(int))*TarMax);
      Covered = (int *) (Tarray + TarMax);
    }
  if (SegMax < 0)
    { SegMax = 1000;
      SegSols = (Local_Segment *) safe_malloc(sizeof(Local_Segment)*SegMax);
    }

  i = 0;
  b = Traplist;
  for (i = 0; i < Traplen; i++)
    { Tarray[i] = b;
      Covered[i] = 0;
      b = b->next;
    }

  qsort(Tarray,Traplen,sizeof(Trapezoid *),TSORT);

#ifdef REPORT_DPREACH
  Al_depth = 0;
#endif
  if (start) NumSegs = 0;
  fseg = NumSegs;
  for (i = 0; i < Traplen; i++)
    if (! Covered[i])
      { b = Tarray[i];
        if (b->top - b->bot < kmerlen) continue;
	//fprintf(stderr, "Trying hit %d\n",i);
        Align_Recursion(A,Alen,B,Blen,b,i,comp,MinLen,MaxDiff,Traplen);
      }

  if (NumSegs > fseg)
    { int i, j;

      qsort(SegSols+fseg,NumSegs-fseg,sizeof(Local_Segment),StSORT);
      for (i = fseg; i < NumSegs; i = j)
        { for (j = i+1; j < NumSegs; j++)

            { if (SegSols[j].abpos != SegSols[i].abpos) break;
              if (SegSols[j].bbpos != SegSols[i].bbpos) break;
	      if (/* segments in opposite orientations */
		  ((SegSols[j].bepos-SegSols[j].bbpos) > 0 &&
		  (SegSols[i].bepos-SegSols[i].bbpos) < 0 )
		  ||
		  ((SegSols[j].bepos-SegSols[j].bbpos) < 0 &&
		  (SegSols[i].bepos-SegSols[i].bbpos) > 0 ) )break;
	  #define KEEP_LONGEST_OVER_MAXDIFF
	  #ifdef KEEP_LONGEST_OVER_MAXDIFF
              if (SegSols[j].error <= MaxDiff &&
		  SegSols[i].error <= MaxDiff){
		if(

		   abs(SegSols[i].bepos-SegSols[i].bbpos)+abs(SegSols[i].aepos-SegSols[i].abpos)
		   <
		   abs(SegSols[j].bepos-SegSols[j].bbpos)+abs(SegSols[j].aepos-SegSols[j].abpos)

		   ){
		  SegSols[i].score=-1;i=j;
		} else {
		  SegSols[j].score=-1;
		}
	      } else {
		if(SegSols[j].error<=MaxDiff){
		  SegSols[i].score=-1;i=j;
		} else {
		  SegSols[j].score=-1;
		}
	      }
	  #else
              if (SegSols[j].score > SegSols[i].score)
                { SegSols[i].score = -1; i = j; }
              else
                SegSols[j].score = -1;
	  #endif
            }
        }

      qsort(SegSols+fseg,NumSegs-fseg,sizeof(Local_Segment),FnSORT);
      for (i = fseg; i < NumSegs; i = j)
        { for (j = i+1; j < NumSegs; j++)
            { if (SegSols[j].abpos != SegSols[i].abpos) break;
              if (SegSols[j].bbpos != SegSols[i].bbpos) break;
              if (SegSols[j].score > SegSols[i].score)
                { SegSols[i].score = -1; i = j; }
              else
                SegSols[j].score = -1;
            }
        }

      for (i = fseg; i < NumSegs; i++)
        if (SegSols[i].score >= 0)
          SegSols[fseg++] = SegSols[i];
      NumSegs = fseg;
    }

#ifdef REPORT_SIZES
  fprintf(stderr, "\n  %9d segments\n",NumSegs);
#endif

  *Seglen = NumSegs;
  return (SegSols);
}


/*** MASTER ROUTINE ***/

Local_Segment *Find_Local_Segments
                  (char *A, int Alen, char *B, int Blen, int Action,
                   int MinLen, float MaxDiff, int *Seglen)
{ static int   DagMax = -1;
  static int AseqLen = -1, BseqLen = -1;
  static char *Alast = NULL;
  int        numhit;
  HitRecord *hits;
  int        numtrap;
  Trapezoid *traps;
  int            numseg;
  Local_Segment *segs = NULL;


  bpadding=   kmerlen+2;

#ifdef VARIABLE_SCORE_SCHEME
  samecost=ceil(100.*MaxDiff);
  diffcost=100-samecost;
  MATCHCOST=samecost+diffcost;
  BLOCKCOST=diffcost*MAXIGAP;
#endif

  if(AseqLen<Alen){
    ArevC = (char*) safe_realloc(ArevC,sizeof(char)*(Alen+1));
    AseqLen=Alen;
  }
  strcpy(ArevC,A);
  reverseComplementSequence(ArevC, Alen);

  if(BseqLen<Blen){
    BrevC = (char*) safe_realloc(BrevC,sizeof(char)*(Blen+1));
    BseqLen=Blen;
  }
  strcpy(BrevC,B);
  reverseComplementSequence(BrevC, Blen);



  if (Alen >= DagMax || Blen >= DagMax)
    { if (Kmask < 0)
        { int i;
          for (i = 0; i < 128; i++)
            Map[i] = -1;
          Map['a'] = Map['A'] = 0;
          Map['c'] = Map['C'] = 1;
          Map['g'] = Map['G'] = 2;
          Map['t'] = Map['T'] = 3;
          Kmask = (1 << (2*kmerlen)) - 1;
          Table = (int *)safe_malloc(sizeof(int)*(Kmask+2));
        }
      if (Alen > Blen)
        DagMax = (int)(1.2*Alen) + 5000;
      else
        DagMax = (int)(1.2*Blen) + 5000;
      DagMax += sizeof(DiagRecord) - (DagMax % sizeof(DiagRecord));
      Tuples = (int *)safe_realloc(Tuples,sizeof(int)*DagMax + sizeof(DiagRecord)*(2*DagMax+maxerror+1));
      DiagVec = ((DiagRecord *) (Tuples + DagMax)) + (DagMax+1);
    }

#undef REUSE_TABLE_UNSAFE
#ifdef REUSE_TABLE_UNSAFE
  if (A != Alast)
#endif
    TableBuild(A,Alen);

#ifdef REPORT_SIZES
  fprintf(stderr, "\nFind_Local Stats for %d x %d comparison:\n",Alen,Blen);
#endif

  { int start;

    start = 1;
    if (Action != LOCAL_REVR)
      {
#ifdef REPORT_SIZES
        fprintf(stderr, "\n  Forward:\n");
#endif
        hits  = Find_Hits(A,Alen,B,Blen,&numhit);
        traps = Build_Trapezoids(A,Alen,B,Blen,hits,numhit,&numtrap);
        segs  = Align_Trapezoids(A,Alen,B,Blen,traps,numtrap,
                                 start,0,MinLen,MaxDiff,&numseg);
        start = 0;
      }

    if (Action != LOCAL_FORW)
      {
#ifdef REPORT_SIZES
        fprintf(stderr, "\n  Reverse:\n");
#endif
        hits  = Find_Hits(A,Alen,BrevC,Blen,&numhit);
        traps = Build_Trapezoids(A,Alen,BrevC,Blen,hits,numhit,&numtrap);
        segs  = Align_Trapezoids(A,Alen,BrevC,Blen,traps,numtrap,
                                 start,1,MinLen,MaxDiff,&numseg);
      }
  }

  Alast = A;
  *Seglen = numseg;
  return (segs);
}
