// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Clark Mobarry
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "GF_ALN_local.h"

// Note, KMERLEN 5, MINMATCH 20, MAXERROR 2, KTHRESH 6 is reasonable
// for performing fragment against fragment comparisons; it permits
// relatively small segments to be found; but it will not give
// acceptable run time for large comparisons such as a BAC against a
// BAC etc.  So, ...

#define KMERLEN   6   // Must be >= 1
#define MINMATCH 20   // (MINMATCH-KMERLEN) is the maximum jump distance.
#define MAXERROR  2   // maximum slop in diagnols for chaining KMERLEN hits.

// The minimum number kmer hits that constitutes an acceptable chain.
#define KTHRESH   (MINMATCH - (KMERLEN-1) - KMERLEN*MAXERROR)

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)




/* D.P. extension alignment scoring */

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

#define MAXIGAP_DEFAULT 3

static int MAXIGAP=MAXIGAP_DEFAULT;


/*amount to subtract from score for mismatch */
//#define DIFFCOST 14
// we can define it to 14 in order not to extend the alignment
// only at a high level of stringency
#define DIFFCOST 3

/*amount to add to score for match*/
#define SAMECOST 1

static int diffcost=DIFFCOST;
static int samecost=SAMECOST;

/* Trapezoid merging padding */

#define DPADDING   2
#define BPADDING   KMERLEN+2

static int BLOCKCOST = DIFFCOST*MAXIGAP_DEFAULT;
static int MATCHCOST = DIFFCOST+SAMECOST;


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


/*** UTILITY ROUTINES ***/

static
void
OutOfMemory(char const * const where) {
  fprintf(stderr,"COMPARE_LOCAL: Out of memory (%s)\n",where);
  exit (1);
}

static void Complement(char * const seq, int const len) {
  static char WCinvert[256];
  static int Firstime = 1;

  if (Firstime) {          /* Setup complementation array */
    int i;

    Firstime = 0;
    for(i = 0; i < 256;i++){
      WCinvert[i] = '?';
    }
    WCinvert[(int)'a'] = 't';
    WCinvert[(int)'c'] = 'g';
    WCinvert[(int)'g'] = 'c';
    WCinvert[(int)'t'] = 'a';
    WCinvert[(int)'n'] = 'n';
    WCinvert[(int)'A'] = 'T';
    WCinvert[(int)'C'] = 'G';
    WCinvert[(int)'G'] = 'C';
    WCinvert[(int)'T'] = 'A';
    WCinvert[(int)'N'] = 'N';
    WCinvert[(int)'-'] = '-'; // added this to enable alignment of gapped consensi
  }

  /* Complement and reverse sequence */

  {
    register char *s, *t;
    int c;

    s = seq;
    t = seq + (len-1);
    while (s < t)
      { c = *s;
      *s++ = WCinvert[(int) *t];
      *t-- = WCinvert[c];
      }
    if (s == t)
      *s = WCinvert[(int) *s];
  }
}


/*** INDEX CONSTRUCTION AND APPLICATION TO FILTERING ***/

/* Shared index and filter arrays used in this subsection */

typedef struct {
  int minim;
  int maxim;
  int count;
} DiagRecord;

static int  Kmask = -1;
static int *Table = NULL;   /* [0..Kmask+1] */
static int *Tuples = NULL;  /* [0..<Seqlen>-KMERLEN] */
static int  Map[128];

static DiagRecord *DiagVec; /* [-(Alen-KMERLEN)..(Blen-KMERLEN) + MAXERROR] */


/* Reverse complement sequences -- so we do not recompute them over and over */
static char *BrevC=NULL;

/* Build index table for sequence S of length Slen. */

static void TableBuild(char const * const S, int const Slen)
{ int   i, c;
 int   x, h;
 char const * const s = S+(KMERLEN-1);

 for (c = 0; c <= Kmask; c++)
   Table[c] = 0;

 h = -KMERLEN;
 c = 0;
 for (i = 0; i < KMERLEN-1; i++)
   { x = Map[(int) (S[i])];
   if (x >= 0)
     c = (c << 2) | x;
   else
     { c <<= 2; h = i-(KMERLEN-1); }
   }
 for (i = 0; i <= Slen-KMERLEN; i++)
   { x = Map[(int) (s[i])];
   if (x >= 0)
     c = ((c << 2) | x) & Kmask;
   else
     { c = (c << 2) & Kmask; h = i; }
   if (i >= h+KMERLEN)
     Table[c+1] += 1;
   }

 for (c = 2; c <= Kmask; c++)
   Table[c] += Table[c-1];

 h = -KMERLEN;
 c = 0;
 for (i = 0; i < KMERLEN-1; i++)
   { x = Map[(int) (S[i])];
   if (x >= 0)
     c = (c << 2) | x;
   else
     { c <<= 2; h = i-(KMERLEN-1); }
   }
 for (i = 0; i <= Slen-KMERLEN; i++)
   { x = Map[(int) (s[i])];
   if (x >= 0)
     c = ((c << 2) | x) & Kmask;
   else
     { c = (c << 2) & Kmask; h = i; }
   if (i >= h+KMERLEN)
     Tuples[Table[c]++] = i;
   }

 for (c = Kmask; c >= 0; c--)
   Table[c+1] = Table[c];
 Table[0] = 0;
}

/* Apply index to find filtered hits between sequences, returning pointer to
   array of HitRecords of length in the integer pointed at by Hitlen       */

static int HSORT(const void *l, const void *r)
{ HitRecord *x, *y;
 x = (HitRecord *) l;
 y = (HitRecord *) r;
 return (x->bstart - y->bstart);
}

static HitRecord *Find_Hits
(char const * const A, int const Alen,
 char const * const B, int const Blen, int * const Hitlen)
{ static int        HitMax = -1;
 static HitRecord *HitList;
 int hits, disconnect;

 if (HitMax < 0)
   { HitMax = 10000;
   HitList = (HitRecord *) malloc(sizeof(HitRecord)*HitMax);
   if (HitList == NULL)
     OutOfMemory("Hit list");
   }

 { int i, j, c;
 int x, h;
 char const * const b = B + (KMERLEN-1);

 for (j = -Alen; j <= Blen+MAXERROR; j++)
   { DiagRecord *dp;
   dp = DiagVec + j;
   dp->count = dp->maxim = 0;
   }

 hits = 0;
 disconnect = MINMATCH - KMERLEN;
 h = -KMERLEN;
 c = 0;
 for (i = 0; i < KMERLEN-1; i++)
   { x = Map[(int) (B[i])];
   if (x >= 0)
     c = (c << 2) | x;
   else
     { c <<= 2; h = i-(KMERLEN-1); }
   }

 for (i = 0; i <= Blen-KMERLEN; i++)
   { x = Map[(int) (b[i])];
   if (x >= 0)
     c = ((c << 2) | x) & Kmask;
   else
     { c = (c << 2) & Kmask; h = i; }
   if (i >= h+KMERLEN)
     for (j = Table[c]; j < Table[c+1]; j++)
       { DiagRecord *dp;
       int e, k;
       k  = i-Tuples[j];
       dp = DiagVec + k;
       for (e = 0; e <= MAXERROR; e++)
         { if (dp->maxim < i-disconnect)
           { if (dp->count >= KTHRESH)
             { HitRecord *hp;
             if (hits >= HitMax)
               { HitMax = (int)(1.2*hits) + 5000;
               HitList = (HitRecord *) realloc(HitList,
                                               sizeof(HitRecord)*HitMax);
               if (HitList == NULL)
                 OutOfMemory("Hit list");
               }
             hp = HitList + hits;
             hp->diagonal = k;
             hp->bstart   = dp->minim;
             hp->bfinish  = dp->maxim + KMERLEN;
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
   }

 for (j = -Alen; j <= Blen+MAXERROR; j++)
   { DiagRecord *dp;
   dp = DiagVec + j;
   if (dp->count >= KTHRESH)
     { HitRecord *hp;
     if (hits >= HitMax)
       { HitMax = (int)(1.2*hits) + 5000;
       HitList = (HitRecord *)realloc(HitList,sizeof(HitRecord)*HitMax);
       if (HitList == NULL)
         OutOfMemory("Hit list");
       }
     hp = HitList + hits;
     hp->diagonal = j;
     hp->bstart   = dp->minim;
     hp->bfinish  = dp->maxim + KMERLEN;
     hits += 1;
     }
   }
 }

 qsort(HitList,hits,sizeof(HitRecord),HSORT);

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

Local_Segment *TraceForwardPath
( char const * const A, int const Alen,
  char const * const B, int const Blen,
  int const mid, int lo, int hi)
{ static Local_Segment rez;
 int *V;
 int  mxv, mxl, mxr, mxi, mxj;
 int  i, j;
 int *Base1, *Base2;

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

 /* Advance to next row */

 for (i = mid; lo <= hi && i < Alen; i++)
   { int  c, v;
   int *W;

   W = V;
   if (V == Base1)
     V = Base2;
   else
     V = Base1;

   v = W[lo];
   c = V[lo] = v - diffcost;
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
       //printf("reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
       }
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
       //printf("reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
       }
     
     for (j++; j <= Blen; j++)
       { v -= diffcost;
       if (v < mxv - BLOCKCOST) break;
       V[j] = v;
       }
     }

   hi = j-1;

   while (lo <= hi && V[lo] < mxv - BLOCKCOST)
     lo += 1;
   while (lo <= hi && V[hi] < mxv - BLOCKCOST)
     hi -= 1;

   if ((i+1) - lo > mxr)
     mxr = (i+1) - lo;
   if ((i+1) - hi < mxl)
     mxl = (i+1) - hi;
   }

 rez.aepos = mxj;
 rez.bepos = mxi;
 rez.ldiag = mxl;
 rez.hdiag = mxr;
 rez.score = mxv;
 return (&rez);
}

Local_Segment *TraceReversePath(char const * const A, int const Alen,
                                char const * const B, int const Blen,
                                int const top, int lo, int hi,
                                int const bot, int xfactor) {
  static Local_Segment rez;
  int *V;
  int  mxv, mxl, mxr, mxi, mxj;
  int  i, j;
  int *Base1, *Base2;

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


    v = W[hi];
    c = V[hi] = v - diffcost;
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
        //printf("reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
        }
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
        //printf("reset mxv = %d at [%d,%d]\n",mxv,mxi,mxj);
        }
     
      for (j--; j >= 0; j--)
        { v -= diffcost;
        if (v < mxv - xfactor) break;
        V[j] = v;
        }
      }

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
    }

  rez.abpos = mxj;
  rez.bbpos = mxi;
  rez.ldiag = mxl;
  rez.hdiag = mxr;
  rez.score = mxv;
  return (&rez);
}


/*** MERGING INDEX HITS INTO TRAPEZOIDAL ZONES ***/

static Trapezoid *Build_Trapezoids(char const * const A, int const Alen,
                                   char const * const B, int const Blen,
                                   HitRecord const * const list, int const Hitlen, int * const Traplen) {
  static Trapezoid  *free = NULL;

  Trapezoid *traporder, *traplist, *tailend;
  Trapezoid *b, *f, *t;
  int i, inserted;
  int trapcount, traparea;

  trapcount = 0;
  traparea  = 0;
  traporder = NULL;
  traplist  = NULL;
  for (i = 0; i < Hitlen; i++) {
    inserted = 0;
    f = NULL;
    for (b = traporder; b != NULL; b = t) {
      t = b->next;
      if (b->top < list[i].bstart - BPADDING) {
        trapcount += 1;
        traparea  += (b->top - b->bot + 1) * (b->rgt - b->lft + 1);
        if (f == NULL)
          traporder = t;
        else
          f->next = t;
        b->next = traplist;
        traplist = b;
      }
      else if (list[i].diagonal > b->rgt + DPADDING)
        f = b;
      else if (list[i].diagonal >= b->lft - DPADDING) {
        if (list[i].diagonal < b->lft)
          b->lft = list[i].diagonal;
        if (list[i].diagonal > b->rgt)
          b->rgt = list[i].diagonal;
        if (list[i].bfinish > b->top)
          b->top = list[i].bfinish;

        if (f != NULL && f->rgt + DPADDING >= b->lft) {
          f->rgt = b->rgt;
          if (f->bot > b->bot) f->bot = b->bot;
          if (f->top < b->top) f->top = b->top;
          f->next = t;
          b->next = free;
          free = b;
        }
        else if (t != NULL && t->lft - DPADDING <= b->rgt) {
          b->rgt = t->rgt;
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
      } else if (! inserted) {
        if (free == NULL) {
          free = (Trapezoid  *)malloc(sizeof(Trapezoid));
          if (free == NULL)
            OutOfMemory("Trapezoid scan list");
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
      } else
        f = b;
    }
    if (! inserted) {
      if (free == NULL) {
        free = (Trapezoid  *)malloc(sizeof(Trapezoid));
        if (free == NULL)
          OutOfMemory("Trapezoid scan list");
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
  }

  for (b = traporder; b != NULL; b = t) {
    t = b->next;
    trapcount += 1;
    traparea  += (b->top - b->bot + 1) * (b->rgt - b->lft + 1);
    b->next  = traplist;
    traplist     = b;
  }

  {
    int lag, lst, lclip;
    int abot, atop;

    for (b = traplist; b != NULL; b = b->next) {
      lag = (b->bot-MAXIGAP)+1;
      if (lag < 0) lag = 0;
      lst = b->top+MAXIGAP;
      if (lst > Blen)
        lst = Blen;

      for (i = lag; i < lst; i++) {
        if (Map[(int) (B[i])] >= 0) {
          if (i-lag >= MAXIGAP) {
            if (lag - b->bot > 0) {
              if (free == NULL) {
                free = (Trapezoid  *)malloc(sizeof(Trapezoid));
                if (free == NULL)
                  OutOfMemory("Trapezoid cutter");
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
          }
          lag = i+1;
        }
      }
      if (i-lag >= MAXIGAP)
        b->top = lag;
    }

    tailend = NULL;
    for (b = traplist; b != NULL; b = b->next) {
      if (b->top - b->bot < KMERLEN) continue;

      abot = b->bot - b->rgt;
      atop = b->top - b->lft;


      lag = (abot - MAXIGAP) + 1;
      if (lag < 0) lag = 0;
      lst = atop + MAXIGAP;
      if (lst > Alen) lst = Alen;

      lclip = abot;
      for (i = lag; i < lst; i++) {
        if (Map[(int) (A[i])] >= 0) {
          if (i-lag >= MAXIGAP) {
            if (lag > lclip) {
              if (free == NULL) {
                free = (Trapezoid  *)malloc(sizeof(Trapezoid));
                if (free == NULL)
                  OutOfMemory("Trapezoid cutter");
                free->next = NULL;
              }
              t = free->next;
              *free = *b;
              b->next = free;
              free = t;

              {
                int x, m;
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

      {
        int x, m;
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
      }

      tailend = b;
    }
  }

  if (tailend != NULL) {
    tailend->next = free;
    free = traplist;
  }

  *Traplen = trapcount;
  return (traplist);
}


/*** FINDING ALIGNMENTS WITHIN A TRAPEZOIDAL ZONE ***/

static int TSORT(const void *l, const void *r) {
  Trapezoid *x, *y;
  x = *((Trapezoid **) l);
  y = *((Trapezoid **) r);
  return (x->bot - y->bot);
}

static int StSORT(const void *l, const void *r) {
  Local_Segment *x, *y;
  x = (Local_Segment *) l;
  y = (Local_Segment *) r;
  if (x->abpos < y->abpos)
    return (-1);
  else if (x->abpos > y->abpos)
    return (1);
  else
    return (x->bbpos - y->bbpos);
}

static int FnSORT(const void *l, const void *r) {
  Local_Segment *x, *y;
  x = (Local_Segment *) l;
  y = (Local_Segment *) r;
  if (x->aepos < y->aepos)
    return (-1);
  else if (x->aepos > y->aepos)
    return (1);
  else
    return (x->bepos - y->bepos);
}

static Trapezoid      **Tarray = NULL;
static int             *Covered;
static Local_Segment   *SegSols = NULL;
static int              SegMax = -1;
static int              NumSegs;

static void Align_Recursion(char const * const A, int const Alen,
                            char const * const B, int const Blen,
                            Trapezoid const * const b, int const current, int const comp,
                            int const MinLen, double const MaxDiff, int const Traplen) {
  int j, mid, indel;
  double pcnt;
  Local_Segment *hend, *lend;
  Trapezoid ltrp, htrp;

  mid = (b->bot + b->top) / 2;

  lend = TraceForwardPath(B,Blen,A,Alen,mid,mid-b->rgt,mid-b->lft);

  {
    int x = 0;

    do {
      x += 1;

      hend = TraceReversePath(B,Blen,A,Alen,
                              lend->bepos,lend->aepos,lend->aepos,
                              mid+MAXIGAP,BLOCKCOST+2*x*diffcost);
    } while (hend->bbpos > mid + x*MAXIGAP && hend->score < lend->score);

    hend->aepos = lend->aepos;
    hend->bepos = lend->bepos;
  }


  ltrp = htrp = *b;
  ltrp.top = min(b->top,hend->bbpos) - MAXIGAP;


  htrp.bot = max(b->bot,hend->bepos) + MAXIGAP;

  if (hend->bepos - hend->bbpos >= MinLen &&
      hend->aepos - hend->abpos >= MinLen) {

    indel = abs( (hend->abpos - hend->bbpos)
                 - (hend->aepos - hend->bepos) );

    pcnt = (-hend->score+samecost*(hend->bepos-hend->bbpos))*1./
      (1.*(MATCHCOST)*(hend->bepos-hend->bbpos));

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
        SegSols = (Local_Segment *) realloc(SegSols,
                                            sizeof(Local_Segment)*SegMax);
        if (SegSols == NULL)
          OutOfMemory("Segment Alignment array");
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
    }
  }

  if (ltrp.top - ltrp.bot > MinLen && ltrp.top < b->top - MAXIGAP)
    Align_Recursion(A,Alen,B,Blen,&ltrp,current,comp,MinLen,MaxDiff,Traplen);

  if (htrp.top - htrp.bot > MinLen)
    Align_Recursion(A,Alen,B,Blen,&htrp,current,comp,MinLen,MaxDiff,Traplen);
}



static Local_Segment *Align_Trapezoids(char const * const A, int const Alen,
                                       char const * const B, int const Blen,
                                       Trapezoid const * const Traplist, int const Traplen,
                                       int const start, int const comp,
                                       int const MinLen, double const MaxDiff, int * const Seglen) {
  static int fseg;
  static int TarMax = -1;

  if (Traplen >= TarMax) {
    TarMax = (int)(1.2*Traplen) + 500;
    Tarray = (Trapezoid **)
      realloc(Tarray,(sizeof(Trapezoid *) + sizeof(int))*TarMax);
    if (Tarray == NULL)
      OutOfMemory("Trapezoid array");
    Covered = (int *) (Tarray + TarMax);
  }
  if (SegMax < 0) {
    SegMax = 1000;
    SegSols = (Local_Segment *) malloc(sizeof(Local_Segment)*SegMax);
    if (SegSols == NULL)
      OutOfMemory("Segment Alignment array");
  }

  {
    Trapezoid * b = (Trapezoid *)Traplist;
    int i;
    for (i = 0; i < Traplen; i++) {
      Tarray[i] = b;
      Covered[i] = 0;
      b = b->next;
    }
  }

  qsort(Tarray,Traplen,sizeof(Trapezoid *),TSORT);

  if (start)
    NumSegs = 0;
  fseg = NumSegs;
  {
    int i;
    for (i = 0; i < Traplen; i++)
      if (! Covered[i]) {
        Trapezoid * b = Tarray[i];
        if (b->top - b->bot < KMERLEN) continue;
        //printf("Trying hit %d\n",i);
        Align_Recursion(A,Alen,B,Blen,b,i,comp,MinLen,MaxDiff,Traplen);
      }
  }

  if (NumSegs > fseg) {
    int i;
    int j=0;
    qsort(SegSols+fseg,NumSegs-fseg,sizeof(Local_Segment),StSORT);
    assert(j==0);
    for (i = fseg; i < NumSegs; i = j) {
      for (j = i+1; j < NumSegs; j++) {
	
        if (SegSols[j].abpos != SegSols[i].abpos) break;
        if (SegSols[j].bbpos != SegSols[i].bbpos) break;
        if (/* segments in opposite orientations */
            ((SegSols[j].bepos-SegSols[j].bbpos) > 0 &&
             (SegSols[i].bepos-SegSols[i].bbpos) < 0 ) 
            ||
            ((SegSols[j].bepos-SegSols[j].bbpos) < 0 &&
             (SegSols[i].bepos-SegSols[i].bbpos) > 0 ) )break;

        if (SegSols[j].error <= MaxDiff &&
            SegSols[i].error <= MaxDiff){

          if (abs(SegSols[i].bepos-SegSols[i].bbpos)+abs(SegSols[i].aepos-SegSols[i].abpos) <
              abs(SegSols[j].bepos-SegSols[j].bbpos)+abs(SegSols[j].aepos-SegSols[j].abpos)) {
            SegSols[i].score=-1;i=j;
          } else {
            SegSols[j].score=-1;
          }
        } else {
          if(SegSols[j].error<=MaxDiff){
            SegSols[i].score=-1;
            i=j;
          } else {
            SegSols[j].score=-1;
          }
        }
      }
    }
    
    qsort(SegSols+fseg,NumSegs-fseg,sizeof(Local_Segment),FnSORT);
    for ( i = fseg; i < NumSegs; i = j) {
      for (j = i+1; j < NumSegs; j++) {
        if (SegSols[j].abpos != SegSols[i].abpos) break;
        if (SegSols[j].bbpos != SegSols[i].bbpos) break;
        if (SegSols[j].score > SegSols[i].score) {
          SegSols[i].score = -1;
          i = j;
        } else
          SegSols[j].score = -1;
      }
    }
    
    for (i = fseg; i < NumSegs; i++)
      if (SegSols[i].score >= 0)
        SegSols[fseg++] = SegSols[i];
    NumSegs = fseg;
  }

  *Seglen = NumSegs;
  return (SegSols);
}












/*** MASTER ROUTINE ***/

Local_Segment *Find_Local_Segments(char const * const A, int const Alen,
                                   char const * const B, int const Blen,
                                   int const Action,
                                   int const MinLen, double const MaxDiff, int * const Seglen) {
  static int   DagMax = -1;
  static int   BseqLen = -1;

  int              numhit;
  HitRecord       *hits;
  int              numtrap;
  Trapezoid       *traps;
  int              numseg;
  Local_Segment   *segs;


  //  Defining this causes scoring to be set so that nearly all
  //  extensions terminate within the user-defined error rate; this
  //  means we don't completely miss a high-quality segment that is
  //  part of a larger segment of lower quality--e.g., we could find
  //  more-conserved portions of a repeat even if the repeat as a
  //  whole was lower fidelity, due to varying selectional pressure at
  //  different positions; the down-side to this is that we won't
  //  chain together perfect matches across a few bases of
  //  lower-fidelity sequence.  Leaving it undefined uses the #define
  //  values of DIFFCOST and SAMECOST
  //
  samecost  = (int)ceil(100.0 * MaxDiff);
  diffcost  = 100 - samecost;
  MATCHCOST = samecost + diffcost;
  BLOCKCOST = diffcost * MAXIGAP;

  if (Action != LOCAL_FORW) {
    if(BseqLen<Blen){
      BrevC = (char*) realloc(BrevC,sizeof(char)*(Blen+1));
      BseqLen=Blen;
      if(BrevC==NULL){OutOfMemory("B sequence reverse complement");}
    }
    strcpy(BrevC,B);
    Complement(BrevC,Blen);
  }

  if (Alen >= DagMax || Blen >= DagMax) {
    if (Kmask < 0) {
      int i;
      for (i = 0; i < 128; i++)
        Map[i] = -1;
      Map[(int)'a'] = Map[(int)'A'] = 0;
      Map[(int)'c'] = Map[(int)'C'] = 1;
      Map[(int)'g'] = Map[(int)'G'] = 2;
      Map[(int)'t'] = Map[(int)'T'] = 3;
      Kmask = (1 << (2*KMERLEN)) - 1;
      Table = (int *)malloc(sizeof(int)*(Kmask+2));
      if (Table == NULL)
        OutOfMemory("K-mer index");
    }
    if (Alen > Blen)
      DagMax = (int)(1.2*Alen) + 5000;
    else
      DagMax = (int)(1.2*Blen) + 5000;
    DagMax += sizeof(DiagRecord) - (DagMax % sizeof(DiagRecord));
    Tuples = (int *)realloc(Tuples,sizeof(int)*DagMax +
                            sizeof(DiagRecord)*(2*DagMax+MAXERROR+1));
    if (Tuples == NULL)
      OutOfMemory("K-mer index");
    DiagVec = ((DiagRecord *) (Tuples + DagMax)) + (DagMax+1);
  }

  TableBuild(A,Alen);

  int start = 1;

  if (Action != LOCAL_REVR) {
    hits  = Find_Hits(A,Alen,B,Blen,&numhit);
    traps = Build_Trapezoids(A,Alen,B,Blen,hits,numhit,&numtrap);
    segs  = Align_Trapezoids(A,Alen,B,Blen,traps,numtrap,
                             start,0,MinLen,MaxDiff,&numseg);
    start = 0;
  }

  if (Action != LOCAL_FORW) {
    hits  = Find_Hits(A,Alen,BrevC,Blen,&numhit);
    traps = Build_Trapezoids(A,Alen,BrevC,Blen,hits,numhit,&numtrap);
    segs  = Align_Trapezoids(A,Alen,BrevC,Blen,traps,numtrap,
                             start,1,MinLen,MaxDiff,&numseg);
  }
  
  *Seglen = numseg;
  return (segs);
} 
