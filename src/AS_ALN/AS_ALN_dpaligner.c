
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

/* Dynamic programming sequence comparison of two fragments.  General
   purpose utility that uses bit-vector d.p. for detection (see, "A Fast
   Bit-Vector Algorithm for Approximate String Matching on Dynamic
   Programming" J. ACM., to appear, by Gene Myers.) and the O(kn) greedy
   algorithm for alignment delivery (see "An O(ND) Difference Algorithm
   and Its Variations" Algorithmica 1 (1986), 251-266, by Gene Myers.)
   Both papers can be downloaded from
           "http://www.cs.arizona.edu/people/gene/vita.html"
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_ALN_aligners.h"

#undef     THRESH_DEBUG
#undef     DP_DEBUG
#undef     BOUND_DEBUG
#undef     WAVE_DEBUG
#undef     AFFINE_DEBUG

#undef   BP_DEBUG
#undef   BP_THRESH
#undef   BP_RATDP
#undef   BP_TAILS

#undef     WARN_SHORT

#define  BP_RATIO .25

#define  LD_RATIO .3 /* longer alignments favored over shorter alignments with error up to this;
			but note that this includes comparing no alignment against any alignment,
			thus this ratio limits the maximum error rate that can meaningfully be
			requested of the DP_Compare aligner */

static double ld_ratio = LD_RATIO;

#define  SUBCOST  2 /* Gap and Substitution costs for affine gap aligner */
#define  GAPCOST  1

#define WORD unsigned long  /* Bit vector unit */

static int WordSize;        /* Size in bits of vector size */

static int  WorkLimit = 0;  /* Current size of 2 arrays below */
static int *HorzDelta;      /* Holds horizontal deltas during d.p. */
static int *DistThresh;     /* Difference threshold values */
static float *DPMatrix;     /* Holds ratio values during branch point d.p. */

/* Probability that there are d or more errors in an alignment of
   length n (sum of substring lengths) over sequences at error rate e */

static double BinomialProb(int n, int d, double e)
{ static int     Nlast = -1, Dlast = -1; /* Last n- and d-values */
  static double  Slast, Elast = -1.;     /* Last answer and e-value */
  static double  LogE, LogC;          /* log e and log (1-e) of last e-value */
  static double *LogTable;            /* LogTable[i] = log(i!) */
  static int     LogMax = -1;         /* Max index for current LogTable */

  if (d == 0) return (1.);

  if (e < 1.e-50) return (0.);

  if (n > LogMax)    /* Need a bigger LogTable */
    { int     max;
      double *newp;
				/* Re-allocate */
      max = (int)((float)n*1.2 + 2048);
      //fprintf(stderr,"DP_COMPARE (BinomialProb): reallocing " F_SIZE_T " bytes\n",(max+1)*sizeof(double));
      newp = (double *) realloc(LogTable,(max+1)*sizeof(double));
      if (newp == NULL) return (-1.);

      { int k;			/* Fill in new part */

        if (LogMax < 0)
          newp[LogMax = 0] = 0.;
        for (k = LogMax+1; k <= max; k++)
          newp[k] = log((double) k) + newp[k-1];
#ifdef THRESH_DEBUG
        printf("\nLog Table Update\n");
        for (k = 0; k <= 9; k += 1)
          printf("    %4d: %g (%g)\n",k,newp[k],exp(newp[k]));
        for (k = 10; k <= max; k += 10)
          printf("    %4d: %g\n",k,newp[k]);
#endif
      }

      LogMax   = max;
      LogTable = newp;
    }

#define LOG_COMB(n,k)  ((LogTable[n] - LogTable[(n)-(k)]) - LogTable[k])

  if (e != Elast)       /* e-value is different than on last call */
    { LogE  = log(e);
      LogC  = log(1.-e);
      Elast = e;
      Nlast = n+1;  /* ==> next if-test below is true */
    }

  { double sum;
    int    k;

    if (n < Nlast || d < Dlast || d < (n-Nlast) + (d-Dlast))
      { sum = 0.;                /* compute from scratch */
        for (k = 0; k < d; k++)
          sum += exp(LOG_COMB(n,k) + (n-k)*LogC + k*LogE);
      }
    else                         /* compute incrementally from last value */
      { sum = Slast;
        for (k = Nlast; k < n; k++)
          sum -= exp(LOG_COMB(k,Dlast-1) + ((k+1)-Dlast)*LogC + Dlast*LogE);
        for (k = Dlast; k < d; k++)
          sum += exp(LOG_COMB(n,k) + (n-k)*LogC + k*LogE);
      }
    if (sum > 1.) sum = 1.;
    Slast = sum;
  }

#ifdef THRESH_DEBUG
  printf("BP(%d,%d,%g) = %g\n",n,d,e,1.-Slast);
#endif
  Nlast = n;
  Dlast = d;
  return (1. - Slast);
}

static int Space_n_Tables(int max, double erate, double thresh)
{ static double LastErate, LastThresh;
  static int    Firstime = 1;

  if (Firstime)  /* Setup bitvector parameters if first call. */
    WordSize = 8*sizeof(WORD);

  { int *newd;  /* If new maximum length, update working structures:
                       DistThresh, HorzDelta, BoundPos, BoundVal
                   and compute new DistThresh entries (if not about
                   to be computed below.                             */
    if (max > WorkLimit)
      { 
	max  = (int)(1.2*max);
        max  = ((max + 2048)/WordSize + 1)*WordSize;
	BinomialProb(2*max,1,.1);
        //fprintf(stderr,"DP_COMPARE (Space_n_Tables): reallocing " F_SIZE_T " bytes\n",(2*max+2)*sizeof(int)+(max+1)*sizeof(float));
        newd = (int *) realloc(DistThresh,
                            (2*max+2)*sizeof(int) + (max+1)*sizeof(float));
        if (newd == NULL) return (1);

        if (!Firstime && LastErate == erate && LastThresh == thresh)
          { int n, d;
            double p;

            d = newd[WorkLimit];
 #undef ERATE_ONE_SEQ
 #ifdef ERATE_ONE_SEQ
 	    /* this doesn't look quite right, even assuming erate is for
 	       error in one seq rather than mismatch; the problem is
 	       that two errors cancel each other only if the changes
 	       happen to be the same ... */
 	    p = 2.*erate - erate*erate; 
 #else
             p = erate;
 #endif
            for (n = WorkLimit + 1; n <= max; n += 1)
 #undef NplusDchooseD
 #ifdef NplusDchooseD
               { while (d <= n && BinomialProb(n+d,d,p) >= thresh) /* why n+d,d, rather than just n,d? */
 #else
               { while (d <= n && BinomialProb(n,d,p) >= thresh) 
 #endif
                 {
		   d += 1;
		   //		   fprintf(stderr,"BP of d=%d\n",d);
		 }

                newd[n] = d;
              }
#ifdef THRESH_DEBUG
            printf("\nThresh Table Extension\n");
            for (n = 1; n <= WorkLimit; n += 1)
              printf("    %5d: %5d\n",n,newd[n]);
#endif
          }

        WorkLimit  = max;
        DistThresh = newd;
        HorzDelta  = DistThresh + (max+1);
        DPMatrix   = (float *) (HorzDelta + (max+1));
      }
  }

  /* If error rate or threshold parameters have changed, or first time
     called, then recompute new probability threshold table.           */

  if (Firstime || LastErate != erate || LastThresh != thresh)
    { int n, d;
      double p;

      LastErate  = erate;
      LastThresh = thresh;
      Firstime   = 0;

      DistThresh[0] = d = 1;
#ifdef ERATE_ONE_SEQ
      /* this doesn't look quite right, even assuming erate is for
	 error in one seq rather than mismatch; the problem is
	 that two errors cancel eachother only if the changes
	 happen to be the same ... */
      p = 2.*erate - erate*erate;
#else
      p = erate;
#endif
      for (n = 1; n <= WorkLimit; n += 1)
#ifdef NplusDchooseD
      { while (d <= n && BinomialProb(n+d,d,p) >= thresh)
#else
        { while (d <= n && BinomialProb(n,d,p) >= thresh)
#endif
	  {
	    //	    fprintf(stdout,"BP of d=%d\n",d);
            d += 1;
	  }
          DistThresh[n] = d;
        }
#ifdef THRESH_DEBUG
      printf("\nNew Thresh Table\n");
      for (n = 1; n <= WorkLimit; n += 1)
        printf("    %5d: %5d\n",n,DistThresh[n]);
#endif
    }

  return (0);
}

/* O(kn) identity-based alignment algorithm.  Find alignment between
   a and b (of lengths alen and blen), that begins at finishing 
   boundary position *spnt.  Return at *spnt the diagonal at which the
   alignment starts.                                                   */

int *AS_ALN_OKNAlign(char *a, int alen, char *b, int blen, int *spnt, int diff)
{ int diag, wpos, level;
  int fcell, infinity;

  static int    Wtop = -1;
  static int   *Wave;
  static int   *TraceBuffer;

  if (diff >= Wtop)        /* Space for diff wave? */
    { int max, del, *newp;

      max = (int)(1.2*diff) + 50;
      del = (max+5)*(max+1);
      //fprintf(stderr,"DP_COMPARE (AS_ALN_OKNAlign): reallocing " F_SIZE_T " bytes\n",del*sizeof(int)+(max+1)*sizeof(int));
      newp = (int *) realloc(Wave,del*sizeof(int) + (max+1)*sizeof(int));
      if (newp == NULL) return (NULL);
      Wtop = max-1;
      Wave = newp;
      TraceBuffer = (int *) (Wave + del);
    }

  diag     = (alen-blen) + (*spnt); /* Finish diagonal. */
  infinity = blen+2;

  /* Process 0-wave. */

  { int i, j;

    if (diff == 0) goto zeroscript;

    if ((*spnt) < 0) /* (i,j) = initial boundary pt. */
      j = blen;
    else
      j = blen - (*spnt);
    i = diag + j;

    while (1)
      { if (i <= 0 || j <= 0) goto zeroscript;
        if (a[i] != b[j]) break;
        i -= 1;
        j -= 1;
      }

    Wave[0] = Wave[1] = infinity;
    Wave[2] = j;
    Wave[3] = Wave[4] = infinity;

#ifdef WAVE_DEBUG
    { int k;
      printf("\nLevel %2d:%*s",0,4*diff,"");
      for (k = 0; k < 5; k++)
        printf(" %3d",Wave[k]);
      printf("\n");
    }
#endif
  }

  /* Compute waves 1 through d-1 do, each wave has
     two boundary cells at each of its ends.       */

  { int m, n, k;

    m = 5;
    n = 0;
    for (level = 1; 1; level++)
      { Wave[m++] = infinity;
        Wave[m++] = infinity;
#ifdef WAVE_DEBUG
        printf("Level %2d:%*s",level,4*(diff-level),"");
        printf(" %3d %3d",infinity,infinity);
#endif
        n += 1;
        for (k = -level; k <= level; k++)
          { int i, j;

            j = Wave[n] - 1;
            if ((i = Wave[n-1]-1) < j)
              j = i;
            if ((i = Wave[n+1]) < j)
              j = i;
            i = (diag+k) + j;
            while (1)
              { if (i <= 0 || j <= 0)
                  { if (i <= 0)
                      *spnt = -j;
                    else
                      *spnt = i;
                    goto madeit;
                  }
                if (a[i] != b[j]) break;
                i -= 1;
                j -= 1;
              }
#ifdef WAVE_DEBUG
            printf(" %3d",j);
#endif
            Wave[m++] = j;
            n += 1;
          }
        Wave[m++] = infinity;
        Wave[m++] = infinity;
        n += 1;
#ifdef WAVE_DEBUG
        printf(" %3d %3d",infinity,infinity);
        printf(" .. %d..%d\n",m-(2*(level+2)+1),m-1);
#endif
      }

madeit:

#ifdef WAVE_DEBUG
    printf(" %3d\n\n",0);
#endif

    fcell = n;
    wpos  = k;
  }

  /* Trace back through wave structure and record
     trace of the alignment traced.                */

  { int d, n, k, t;

    t = 0;
    n = fcell;
    k = wpos;
    for (d = level-1; d >= 0; d--)
      { int i, j, m;

        j = Wave[m=n]-1;
        if ((i = Wave[n-1]-1) < j)
          { j = i; m = n-1; }
        if ((i = Wave[n+1]) < j)
          { j = i; m = n+1; }
        if (m < n)
          { TraceBuffer[t++] = - ((diag+k) + (j+1));
#ifdef WAVE_DEBUG
            printf("Delete b[%d] = %c\n",j+1,b[j+1]); 
#endif
            k -= 1;
          }
        else if (m > n)
          { TraceBuffer[t++] = j+1;
#ifdef WAVE_DEBUG
            printf("Insert a[%d] = %c\n",(diag+k) + (j+1),a[(diag+k) + (j+1)]); 
#endif
            k += 1;
          }
#ifdef WAVE_DEBUG
        else
          { printf("Substitute b[%d] = %c to %c = a[%d]\n",
                   j+1,b[j+1],a[(diag+k) + j+1],(diag+k) + j+1); 
          }
#endif
        n = m - (2*d+4);
      }
    TraceBuffer[t] = 0;
  }

  return (TraceBuffer);

  /* If perfect match, your done. */

zeroscript:
  TraceBuffer[0] = 0;
  *spnt = diag;
  return (TraceBuffer);
}

/* O(kn) affine gap cost alignment algorithm.  Find best alignment between
   a and b (of lengths alen and blen), within a band of width 2*diff centered
   on the diagonal containing finishing boundary position *epnt.  Return at
   *bpnt the diagonal at which the alignment starts.  A quick implementation
   that is space inefficient, takes O(kn) space as opposed to the O(n) that
   is possible.                                                             */

int *AS_ALN_OKNAffine(char *a, int alen, char *b, int blen,
                      int *bpnt, int *epnt, int diff)
{ int diag, jcrd;
  int infinity;
  int bwide, atop;
  int *C, *I, *TraceBuffer, *TraceTwo; 
  int best, bdag = 0;

  static int    Amax  = -1;
  static int   *Afarr = NULL;

  bwide = 2*diff + 1;
  if ((blen+1)*(2*bwide+2) >= Amax)
    {
      Amax = (blen+501)*(2*bwide+202);
      Afarr = (int *) ckrealloc(Afarr,Amax*sizeof(int));
#ifdef AFFINE_DEBUG
      fprintf(stderr,"Affine align allocating %d bytes\n",max*2*sizeof(int));
#endif
    }
  atop = (blen+1)*bwide;
  TraceBuffer = Afarr + 2*atop;
  TraceTwo    = TraceBuffer + (blen+1);

  diag     = (alen-blen) + (*epnt); /* Finish diagonal. */
  infinity = blen*(GAPCOST+SUBCOST+1)+2;

  { int i, k, m;

    if ((*epnt) <= 0) /* (i,j) = initial boundary pt. */
      { jcrd = blen;
        C = Afarr + jcrd*bwide + diff;
        I = C + atop;

        i = diag + jcrd;
        for (k = -diff; k <= diff; k++)
          { m = i+k;
            if (m < 0 || m > alen)
              C[k] = infinity;
            else
              C[k] = 0;
            I[k] = infinity;
          }
      }
    else
      { jcrd = blen - (*epnt);
        C = Afarr + jcrd*bwide + diff;
        I = C + atop;

        i = diag + jcrd;
        for (k = -diff; k <= diff; k++)
          { m = i+k;
            if (m < 0 || m > alen)
              C[k] = infinity;
            else if (m == alen)
              C[k] = 0;
            else
              C[k] = (alen-m) + GAPCOST;
            I[k] = infinity;
          }
      }
  }

  best = infinity;
  for (jcrd--; jcrd >= 0; jcrd--)
    { int i, m, k, x;
      int deljk = 0;

      i = diag + jcrd;
      if (i+diff < 0) break;

      C -= bwide;
      I -= bwide;

      for (k = diff; k >= -diff; k--)
        { m = i+k;
          if (k == -diff)
            I[k] = infinity;
          else
            { x = C[bwide + (k-1)] + GAPCOST;
              if (x > I[bwide + (k-1)])
                x = I[bwide + (k-1)];
              I[k] = x+1;
            }

          if (k == diff)
            deljk = infinity;
          else
            { x = C[k+1] + GAPCOST;
              if (x > deljk)
                x = deljk;
              deljk = x+1;
            }

          if (m < 0 || m > alen)
            C[k] = infinity;
          else if (m == alen)
            C[k] = 0;
          else
            { x = C[bwide + k];
              if (b[jcrd+1] != a[m+1])
                x += SUBCOST;
              if (x > I[k])
                x = I[k];
              if (x > deljk)
                x = deljk;
              C[k] = x;
            }
 
          if (m == 0 && C[k] < best)
            { best = C[k]; bdag = m - jcrd; }
        }

#ifdef AFFINE_DEBUG
      { printf("%3d:  ",jcrd+1);
        for (k = diff; k >= -diff; k--) 
          { if (C[k+bwide] >= infinity)
              printf("     *");
            else
              printf(" %5d",C[k+bwide]);
          }
        printf("\n      ");
        for (k = diff; k >= -diff; k--) 
          { if (I[k+bwide] >= infinity)
              printf("     *");
            else
              printf(" %5d",I[k+bwide]);
          }
        printf("\n  %c     ",b[jcrd+1]);
        for (k = diff+1; k >= -diff; k--) 
          { m = i+k;
            if (m > alen || m <= 0)
              printf(" %3d:*",m,a[m]);
            else
              printf(" %3d:%c",m,a[m]);
          }
        printf("\n%3d:        ",jcrd);
        for (k = diff; k >= -diff; k--) 
          { if (C[k] >= infinity)
              printf("     *");
            else
              printf(" %5d",C[k]);
          }
        printf("\n            ");
        for (k = diff; k >= -diff; k--) 
          { if (I[k] >= infinity)
              printf("     *");
            else
              printf(" %5d",I[k]);
          }
        printf("\n\n");
      }
#endif
    }

  if (jcrd < 0)
    { int m, k;

      C = Afarr + diff;
      for (k = diff; k >= -diff; k--)
        { m = diag+k;
          if (m >= 0 && m <= alen && C[k] < best)
            { best = C[k]; bdag = m; }
        }
    }

#ifdef AFFINE_DEBUG
  printf("Best = %d  bdag = %d\n\n",best,bdag);
#endif

#define C_STATE 0
#define D_STATE 1
#define I_STATE 2

  { int i, k, s, x;
    int deljk = 0, top;

    if (bdag >= 0)
      jcrd = 0;
    else
      jcrd = -bdag;
    i = bdag + jcrd;
    k = bdag - diag;
    s = C_STATE;
    C = Afarr + jcrd*bwide + diff;
    I = C + atop;
    top = 0;
    while (jcrd != blen && i != alen)
      { switch (s)
        { case C_STATE:
            if (b[jcrd+1] != a[i+1])
              x = SUBCOST;
            else
              x = 0;
            if (C[k] == C[bwide+k] + x)
              { i += 1; jcrd += 1; C += bwide; I += bwide; }
            else if (C[k] == I[k])
              s = I_STATE;
            else
              { s = D_STATE; deljk = C[k]; }
            break;
          case I_STATE:
            if (I[k] == C[bwide+(k-1)] + GAPCOST+1) s = C_STATE;
            jcrd += 1;
            C += bwide;
            I += bwide;
            k -= 1;
            TraceBuffer[top++] = -(i+1);
            break;
          case D_STATE:
            if (deljk == C[k+1] + GAPCOST+1) s = C_STATE;
            i += 1;
            k += 1;
            deljk -= 1;
            TraceBuffer[top++] = (jcrd+1);
        }
      }
    TraceBuffer[top] = 0;

#define ADJUST_EPNT
#ifdef ADJUST_EPNT
    *epnt=0;
    if(i!=alen){
      assert(i<alen);
      assert(jcrd==blen);
      *epnt=i-alen;
    } 
    if(jcrd!=blen){
      assert(jcrd<blen);
      assert(i==alen);
      *epnt=blen-jcrd;
    }
#endif
  }

  { int i, j;
    int score, mxlev, mnlev;
    int mxi, mxj;
    int mni = 0, mnj = 0, mnp = 0;
    int top = 0;
    int c, p;       /* Reparse alignment */

#define UPTICK  1
#define DWNTICK 1

    i = j = 1;
    if (bdag > 0)
      i = bdag+1;
    else if (bdag < 0)
      j = (-bdag)+1;
    mxlev = mnlev = score = 0;
    mxi = i-1;
    mxj = j-1;
  
    p = 0;
reloop:
    while ((c = TraceBuffer[p++]) != 0)
    { if (c < 0)
        { c = -c;
          while (i != c)
            if (a[i++] == b[j++])
              { score += UPTICK;
                if (score > mxlev)
                  { if (mnlev < mxlev)
                      { if (mni-mxi != 1 || mnj-mxj != 1)
                          { int k;
                            for (k = mxi+1; k <= mni; k++)
                              TraceTwo[top++] = mxj+1;
                            for (k = mxj+1; k <= mnj; k++)
                              TraceTwo[top++] = -(mni+1);
                          }
                        mxlev = mnlev = score = 0;
                        mxi = mni;
                        mxj = mnj;
                        i = mni+1;
                        j = mnj+1;
                        p = mnp;
                        goto resync;
                      }
                    mxlev = mnlev = score;
                    mxi = i-1;
                    mxj = j-1;
                  }
              }
            else
              { score -= DWNTICK;
                if (score <= mnlev)
                  { mnlev = score;
                    mni = i-1;
                    mnj = j-1;
                    mnp = p-1;
                  }
              }
          score -= DWNTICK;
          j += 1;
          if (score <= mnlev)
            { mnlev = score;
              mni = i-1;
              mnj = j-1;
              if (TraceBuffer[p] == 0)
                mnp = -1;
              else
                mnp = p;
            }
        }
      else
        { while (j != c)
            if (a[i++] == b[j++])
              { score += UPTICK;
                if (score > mxlev)
                  { if (mnlev < mxlev)
                      { if (mni-mxi != 1 || mnj-mxj != 1)
                          { int k;
                            for (k = mxi+1; k <= mni; k++)
                              TraceTwo[top++] = mxj+1;
                            for (k = mxj+1; k <= mnj; k++)
                              TraceTwo[top++] = -(mni+1);
                          }
                        mxlev = mnlev = score = 0;
                        mxi = mni;
                        mxj = mnj;
                        i = mni+1;
                        j = mnj+1;
                        p = mnp;
                        goto resync;
                      }
                    mxlev = mnlev = score;
                    mxi = i-1;
                    mxj = j-1;
                  }
              }
            else
              { score -= DWNTICK;
                if (score <= mnlev)
                  { mnlev = score;
                    mni = i-1;
                    mnj = j-1;
                    mnp = p-1;
                  }
              }
          score -= DWNTICK;
          i += 1;
          if (score <= mnlev)
            { mnlev = score;
              mni = i-1;
              mnj = j-1;
              if (TraceBuffer[p] == 0)
                mnp = -1;
              else
                mnp = p;
            }
        }
      resync:
        ;
      }
  
    while (a[i] != 0)
      { if (b[j] != 0)
          if (a[i++] == b[j++])
            { score += UPTICK;
              if (score > mxlev)
                { if (mnlev < mxlev)
                    { if (mni-mxi != 1 || mnj-mxj != 1)
                        { int k;
                          for (k = mxi+1; k <= mni; k++)
                            TraceTwo[top++] = mxj+1;
                          for (k = mxj+1; k <= mnj; k++)
                            TraceTwo[top++] = -(mni+1);
                        }
                      mxlev = mnlev = score = 0;
                      mxi = mni;
                      mxj = mnj;
                      i = mni+1;
                      j = mnj+1;
                      if (mnp < 0) continue;
                      p = mnp;
                      goto reloop;
                    }
                  mxlev = mnlev = score;
                  mxi = i-1;
                  mxj = j-1;
                }
            }
          else
            { score -= DWNTICK;
              if (score <= mnlev)
                { mnlev = score;
                  mni = i-1;
                  mnj = j-1;
                  mnp = -1;
                }
            }
	  else
            break;
      }

    if (score <= mxlev)
      if (i-mxi != 2 || j-mxj != 2)
      { int k;
        for (k = mxi+1; k <= i-1; k++)
          TraceTwo[top++] = mxj+1;
        for (k = mxj+1; k <= j-1; k++)
          TraceTwo[top++] = -(i);
      }
    TraceTwo[top] = 0;
  }

  *bpnt = bdag;
  return (TraceTwo);
}

/* Do bit-vector based d.p. on thresholded zone where a-boundary is 0's
   and zone starts in columns [beg,end] of problem.  If the zone reaches
   a boundary modeling a sufficiently large overlap, then the best overlap
   position is returned, and the difference in the overlap is placed in
   integer pointed at by diff.  Otherwise the impossible position blen
   is returned.                                                          */

static int Boundary(char *a, int alen, char *b, int blen,
                    int beg, int end, int minlen, int *diff)
{ int lft,  rgt;
  int lval, rval;
  int prob_thresh;
  int boundpos, boundval;
  int preminpos, preminval;
  int lastlocalminpos, lastlocalminscore,lastlft;

  static int Firstime = 1;
  static WORD bvect[256];	/* bvect[a] is equal-bit vector of symbol a */
  static int  slist[256], stop; /* slist[0..stop-1] == symbols in current
                                   segment of b being compared.           */
#ifdef DP_DEBUG
  printf("\nBoundary (%d,%d):\n",beg,end);
#endif

  if (end > alen) end = alen;

  if (Firstime)
    { int i;

      Firstime = 0;	/* Setup empty bvect array if first ever call */
      for (i = 0; i < 256; i++)
        bvect[i] = 0;
      stop = 0;
    }

  { int i;

    lft  = beg;		/* Setup horizontal delta on a-boundary */
    rgt  = end;
    lval = 0;
    rval = 0;
    for (i = lft; i <= rgt; i++)
      HorzDelta[i] = 0;
    prob_thresh = 1;
    boundpos = 0;
    boundval = 0;


    preminval=0;
    preminpos=0;

    lastlocalminpos=0;
    lastlocalminscore=alen-end;
    lastlft=0;
  }

  { int j, bmax;	/* For every WordSize'th row do */ 
			
    bmax = WordSize*((blen-1)/WordSize+1);
    for (j = 0; j < blen; j += WordSize)
      { WORD P, M, ebit = 0;
        int  row, bval;

        { WORD One; /* Compute list of chars in b[j+1..MAX(j+WordSize,blen)] */
          int i;    /* and non-zero bvect[a] for each one.                   */

          row = j + WordSize;
          if (row > blen) row = blen;
          One = 1;
          for (i = j+1; i <= row; i++)
            { if (bvect[(int)(b[i])] == 0)
                { slist[stop++] = b[i];
                  bvect[(int)(b[i])] = One;
                }
              else
                bvect[(int)(b[i])] |= One;
              ebit = One;
              One <<= 1;
            }
        }

#ifdef DP_DEBUG
        { int i, k;

          printf("\nSUBSEGMENT A(%d,%d) x B(%d,%d)\n",
                 lft,rgt,j+1,row);
          printf("\n            ");
          for (i = j+1; i <= row; i++)
            printf("  %c",b[i]);
          printf("\n %3d:     ",lft);
          k = lval;
          for (i = j+1; i <= row; i++)
            printf(" %2d",k++);
          printf(" %2d",k);
          printf("\n");
        }
#endif

        { int i, r, dj;          /* Must compute WordSize additional columns */
                                 /* beyond right border rgt.                 */
          dj = DistThresh[row];
          r = rgt + WordSize + (dj - prob_thresh);
          if (r > alen) r = alen;
          for (i = rgt+1; i <= r; i++)
            { HorzDelta[i] = 1;
              rval += 1;
            }
          rgt = r;
          prob_thresh = dj;
        }

        { int  i, h;            /* Do d.p. accross blocks to get deltas at  */
          WORD U, Y, X, mc, pc; /*   next level using bit-vector approach.  */
                                /*   of Myers (JACM, 1999)                  */
          M = 0;
#if 0 
          P = (WORD)(-1); // converted to unsigned via twos-complement.
#else
          P = 0; P = ~ P;  // for g++ to be happy.
#endif
          bval  = rval;
          lval += row-j;
	  rval  = lval;
          for (i = lft+1; i <= rgt; i++)
            { h  = HorzDelta[i];
#ifdef DP_DEBUG
              printf(" %3d: %2d %c",i,h,a[i]);
#endif
              pc = (h > 0);
              mc = (h < 0);

              U  = bvect[(int)(a[i])];
              Y  = U | mc;
              X  = (((Y & P) + P) ^ P) | Y;
              U |= M;

              Y = P;
              P = M | ~ (X | Y);
              M = Y & X;

              if (P & ebit)
                h = +1;
              else if (M & ebit)
                h = -1;
              else
                h = 0;

              Y = (P << 1) | pc;
              P = (M << 1) | mc | ~ (U | Y);
              M = Y & U;

              rval += (HorzDelta[i] = h);

#ifdef DP_DEBUG
              { WORD p, m;
                int  k, d;

                p = P;
                m = M;
                d = 0;
                for (k = j+1; k <= row; k++)
                  { if (p & 0x1)
                      d += 1;
                    else if (m & 0x1)
                      d -= 1;
                    m >>= 1;
                    p >>= 1;
                  }
                p = P;
                m = M;
                d = rval - d;
                printf(" %2d",d);
                for (k = j+1; k <= row; k++)
                  { if (p & 0x1)
                      d += 1;
                    else if (m & 0x1)
                      d -= 1;
                    printf(" %2d",d);
                    m >>= 1;
                    p >>= 1;
                  }
                printf("\n");
              }
#endif
            }
        }

        if (rgt == alen) /* Check values at right b-boundary */
          { int i, p;

            p = bval;
#ifdef BOUND_DEBUG
            printf("\n  Boundary Segment: P-thresh = %d row %d\n",prob_thresh,row);
#endif
            for (i = j+1; i <= row; i++) 
              { if (P & 0x1)
                  p += 1;
                else if (M & 0x1)
                  p -= 1;
                M >>= 1;
                P >>= 1;
#ifdef BOUND_DEBUG
                printf("  B%03d: %3d (%3d,%3d)",i,p,boundpos,boundval);
#endif




		/* changes to criteria for keeping overlaps : ALH 12/17/01 */

		// if p is small enough
		if (p < prob_thresh ){

		  // if overlap is shorter than minlen
		  if( i < minlen) {
		    // if better than best so far
		    if ((i-preminpos)*ld_ratio + preminval >= p) {


#ifdef BOUND_DEBUG
		      printf(".");
#endif

		      // store best overlap that is too short
		      preminpos = i;
		      preminval = p;
		    } 
		    //if consistent with adding gaps to pad prev. local best
		    if (i-lastlocalminpos==p-lastlocalminscore){
#ifdef BOUND_DEBUG		      
		      printf(",");
#endif
		    }else { 
		      // if new local best -- i.e. prev local best plus gaps is worse than here
		      if (i-lastlocalminpos>p-lastlocalminscore){
			lastlocalminpos=i;
			lastlocalminscore=p;
			lastlft=lft;
#ifdef BOUND_DEBUG		      
			printf("~");
#endif

		      } else { 

			//this should not happen???

			//... except that the
			//previous best might have been in a different
			//stripe (set of WordSize bases) and might
			//have involved an entry point that is now
			//lost due to restrictions on lft and rgt.
			//
			//... and perhaps other cases not yet forseen

			#if 0
			if( lft<=lastlft || i % WordSize != 1){
			  fprintf(stderr,"Possible logic problem: lft %d lastlft %d i %d WordSize %d i%%WordSize %d\n",
				  lft,lastlft,i,WordSize,i%WordSize);
			}
			assert( lft>lastlft && i % WordSize == 1);
			#endif

			//In this case, we really should consider this
			//location a new local best, so, even though
			//it is strictly worse,

			lastlocalminpos=i;
			lastlocalminscore=p;
			lastlft=lft;
#ifdef BOUND_DEBUG		      
			printf("^");
#endif

		      }
		    }
		  } else { // if overlap is long enough
		    // if better than best so far

#ifdef BOUND_DEBUG
		    printf(" [%d >?= %d (%f)] ",(int)((i-boundpos)*ld_ratio + boundval),p,p/(double)(i));
#endif

		    if ((i-boundpos)*ld_ratio + boundval >= p) {

		      // if equivalent to best_too_short plus external gaps
		      if( i-preminpos == p-preminval){
#ifdef BOUND_DEBUG
			printf("=");
#endif

#ifdef WARN_SHORT
			fprintf(stderr,
"WARNING: DP_Compare using overlap shorter than minlen\n"
"WARNING: because gaps external to the overlap can be placed to\n"
"WARNING: create an overlap that passes all tests and is long enough\n");
#endif
			// ... favor the short overlap
			// (the idea is ... the longer overlap was
			// good enough to be used, so we'll allow *some*
			// overlap, but we'd rather return the short one
			// than an overlap with external gaps
			boundpos=preminpos;
			boundval=preminval;

		      } else { 

			// if consistent with extension of previous local best
			if (i-lastlocalminpos == p-lastlocalminscore){

#ifdef WARN_SHORT
			  fprintf(stderr,
"WARNING: DP_Compare using overlap shorter than minlen\n"
"WARNING: because gaps external to the overlap can be placed to\n"
"WARNING: create an overlap that passes all tests and is long enough\n");
#endif

#ifdef BOUND_DEBUG
			  printf("^");
#endif

			  boundpos=lastlocalminpos;
			  boundval=lastlocalminscore;

			} else { // this overlap is the best so far
			  boundpos = i;
			  boundval = p;
#ifdef BOUND_DEBUG
			  printf("+\\");
#endif
			}
		      }
		    }
#ifdef BOUND_DEBUG
		  }
		  printf("*");
		}
		printf("\n");
#else
	      }
	  }
#endif



      }
          }

        while (stop > 0)         /* Reset bvect values for symbols in */
          bvect[slist[--stop]] = 0;  /*    current segment                */

        /* Trim back left and right boundaries using  */
        /*   threshold limits on d-values             */

#ifdef DP_DEBUG
        printf("  Thresh = %d\n",prob_thresh);
#endif
        while (rgt >= lft)
          { if (rval < prob_thresh)
              break;
	  //	  printf("rval %d still greater than prob_thresh %d for rgt %d\n",rval,prob_thresh,rgt);
            rval -= HorzDelta[rgt--];
          }

        if (rgt < lft) break;
    
        while (lval >= prob_thresh){
	  //	  printf("lval %d still greater than prob_thresh %d for lft %d\n",lval,prob_thresh,lft);
          lval += HorzDelta[++lft];
	}

#ifdef DP_DEBUG
        printf("  Range: (%d,%d)\n",lft,rgt);
#endif
      }

    boundpos = blen - boundpos;

    if (j >= blen && lft <= rgt)  /* Reached a-boundary: check values there */
      { int i, v, p;

        v = boundval + (int)(((double)boundpos)*ld_ratio);
#ifdef BOUND_DEBUG
	printf("v = %d = %d + [ %d = %d * %f\n",
	       v,boundval,(int)(((double)boundpos)*ld_ratio),
	       boundpos,ld_ratio);
#endif
        p = (int)rval;
        if (lft < minlen) lft = minlen;
        for (i = rgt; i >= lft; i--)
          {

#ifdef BOUND_DEBUG
            if (i == rgt)
              printf("\n  Boundary Segment: P-thresh = %d\n",prob_thresh);
            printf("  A%03d: %3d (%3d,%3d)",i,p,boundpos,boundval);
#endif

            if (p < prob_thresh)
              { if (v >= p)
                  { boundpos = i-alen;
                    boundval = v = p;
#ifdef BOUND_DEBUG
                    printf("+/");
                  }
                printf("*");
              }
            printf("\n");
#else
                  }
              }
#endif

            p -= HorzDelta[i];
          }
      }
  }

  *diff = boundval;

#ifdef BOUND_DEBUG
printf("Boundary returning (%d, %d)\n",boundval,boundpos);
#endif

  return (boundpos);
}

/* Given fragments a and b, find the best overlap between them subject
   to the following parameters/thresholds.  The overlap must start on
   one of the diagonals of the d.p. matrix in the interval [beg,end].
   For example if one gives the interval [-10,20], then the overlap
   either has less than the first 20bp of a unaligned or less than the
   first 10bp of b unaligned.  If the boolean variable opposite is nonzero
   then the fragments are to be aligned in the opposite orientation.  One
   is assuming an error rate of erate in the sequences, and is guaranteed
   to find only alignments for which the number of differences d in each
   prefix of length n is such that
           Sum_k=d^n (n choose k) erate^k (1-erate)^(n-k) < thresh.
   One should note carefully, that alignments not satisfying this property
   may be found, the point is that ones that don't may be missed.
   In addition, the alignment must involve at least minlen symbols of the
   prefix-sequence in the overlap.  The option what specifies what kind
   of comparison is to be performed as follows:

   AS_FIND_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.
   AS_FIND_ALIGN_NO_TRACE
      A hack by SAK...don't compute delta encoding.
   AS_FIND_AFFINE_ALIGN
      For this option, find the best alignment using an affine gap cost
      measure.  Substitutions cost SUBCOST and gaps of length len cost
      GAPCOST + len, where SUBCOST and GAPCOST are defined constants within
      AS_ALN_dpaligner.c
   AS_FIND_QVALIGN:
      NOT YET IMPLEMENTED.  Will ultimately use quality values to compute
      the best possible alignment.

   by the routine, and must be copied if it is to be retained beyond the
   given call.                                                            */

OverlapMesg *DP_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where)
{
  char *aseq, *bseq;
  int  ahang, bhang;
  int  *trace;
  Overlap *rawOverlap;
  static OverlapMesg QVBuffer;  //Note: return is static storage--do not free

  aseq = a->sequence;  /* Setup sequence access */
  bseq = b->sequence;

  rawOverlap=DP_Compare(aseq,bseq,beg,end,opposite,erate,thresh,minlen,what);

  if(rawOverlap==NULL)return NULL;

  ahang    = rawOverlap->begpos;
  bhang    = rawOverlap->endpos;
  trace    = rawOverlap->trace;

  if (ahang < 0 || (ahang == 0 && bhang > 0))
    { if (bhang >= 0)
        QVBuffer.overlap_type = AS_CONTAINMENT;
      else
        QVBuffer.overlap_type = AS_DOVETAIL;
      QVBuffer.ahg = -ahang;
      QVBuffer.bhg = -bhang;
      if (opposite)
        QVBuffer.orientation = AS_OUTTIE;
      else
        QVBuffer.orientation = AS_NORMAL;
      QVBuffer.aifrag = b->iaccession;
      QVBuffer.bifrag = a->iaccession;
      if (trace != NULL)
        { int i;
          for (i = 0; trace[i] != 0; i++)
            trace[i] = -trace[i];
        }
    }
  else
    { if (bhang <= 0)
        QVBuffer.overlap_type = AS_CONTAINMENT;
      else
        QVBuffer.overlap_type = AS_DOVETAIL;
      QVBuffer.ahg = ahang;
      QVBuffer.bhg = bhang;
      if (opposite)
        QVBuffer.orientation = AS_INNIE;
      else
        QVBuffer.orientation = AS_NORMAL;
      QVBuffer.aifrag = a->iaccession;
      QVBuffer.bifrag = b->iaccession;
    }
  if (what != AS_FIND_ALIGN_NO_TRACE && trace != NULL)
    QVBuffer.delta = Pack_Alignment_AS(trace,QVBuffer.ahg);
  else
    QVBuffer.delta = NULL;

  *where = ahang;

  return (&QVBuffer);
}




OverlapMesg *LD_DP_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where, double my_ld_ratio)
{
  OverlapMesg *O;
  double old_ld_ratio=ld_ratio;
  ld_ratio=my_ld_ratio;
  O=DP_Compare_AS(a,b,beg,end,opposite,erate,thresh,minlen,what,where);
  ld_ratio=old_ld_ratio;
  return(O);
}



/* WHAT WAS THE AS_CNS VERSION "DP_Compare()" IS NOW IN HERE */

/* Given fragments a and b, find the best overlap between them subject
   to the following parameters/thresholds.  The overlap must start on
   one of the diagonals of the d.p. matrix in the interval [beg,end].
   For example if one gives the interval [-10,20], then the overlap
   either has less than the first 20bp of a unaligned or less than the
   first 10bp of b unaligned.  If the boolean variable opposite is nonzero
   then the fragments are to be aligned in the opposite orientation.  One
   is assuming an error rate of erate in the sequences, and is guaranteed
   to find only alignments for which the number of differences d in each
   prefix of length n is such that
           Sum_k=d^n (n choose k) erate^k (1-erate)^(n-k) < thresh.
   One should note carefully, that alignments not satisfying this property
   may be found, the point is that ones that don't may be missed.
   In addition, the alignment must involve at least minlen symbols of the
   prefix-sequence in the overlap.  The option what specifies what kind
   of comparison is to be performed as follows:

   AS_FIND_OVERLAP:
      Just find a good alignment to the boundary of the d.p. matrix.  From
      this extrapolate a rough overlap relationship without an alignment
      and return the result (if there is one).
   AS_FIND_ALIGN:
      For this option, further go to the trouble of computing the alignment
      and store it in the overlap message.

   As for all other routines, the space for the overlap message is owned
   by the routine, and must be copied if it is to be retained beyond the
   given call.                                                            */

Overlap *DP_Compare(char *aseq, char *bseq,
                    int beg, int end, int opposite,
                    double erate, double thresh, int minlen,
                    CompareOptions what)
{ int   alen,  blen;
//  int   atop,  btop;
  int   pos1,  pos2;
  int   dif1,  dif2;

  static Overlap OVL;
  assert(erate>=0&&erate<1);
  ld_ratio=erate/(1.-erate); 
                             /* we want to balance ld_ratio with erate such that any alignment
				with mismatch rate <= erate will have mismatches <= ld_ratio * aligned length of one aligned substring;

				max mismatches for a given length is achieved by putting that number of gaps in the short side, i.e.
				max mismatches = d such that d/(n+d) = erate, i.e. d = n*erate/(1-erate)

				but we also want d <= ld_ratio * n, or n*erate/(1-erate) <= ld_ratio *n, so ld_ratio >= erate/(1-erate)
			     */

  alen = strlen(aseq);
  blen = strlen(bseq);
  aseq -= 1;
  bseq -= 1;

  if (end > alen) end = alen;
  if (beg < -blen) beg = -blen;

  { int max;

    if (alen > blen)
      max = alen;
    else
      max = blen;
    if (Space_n_Tables(max,erate,thresh)) return (NULL);
  }

  if (opposite)                 /* Compare in opposite orientation. */
    Complement_Seq(bseq+1);

  { int mid;

    /* Compute critical overlap finish pt. from a-boundary (if any) */

    if (end > 0 || (beg == 0 && end == 0))
      { mid = beg;
        if (beg < 0) mid = 0;
        pos1 = Boundary(aseq,alen,bseq,blen,mid,end,minlen,&dif1);
      }
    else
      pos1 = blen;

#ifdef BOUND_DEBUG
    printf("\nA vs. B Critical Point (%d,%d)\n",pos1,dif1);
#endif

    /* Compute critical overlap finish pts. from b-boundary (if any) */

    if (beg < 0)
      { mid = end;
        if (end > 0) mid = -1;
        pos2 = - Boundary(bseq,blen,aseq,alen,-mid,-beg,minlen,&dif2);
      }
    else
      pos2 = -alen;

#ifdef BOUND_DEBUG
    printf("\nB vs. A Critical Point (%d,%d)\n",pos2,dif2);
#endif
  }

  /* Need the better of the two best alignments from each boundary.
     Compute the better based on the slope rule and record the best
     in (pos1,dif1).                                                  */

  { int which, olen1, olen2;

    if (pos1 >= blen)
      { if (pos2 <= -alen)
          { if (opposite)
              Complement_Seq(bseq+1);
            return (NULL);
          }
        which = 1;
      }
    else if (pos2 <= -alen)
      which = 0;
    else
      { olen1 = pos1;
        if (olen1 < 0)
          olen1 = blen;
        else
          olen1 = blen-olen1;
        olen2 = pos2;
        if (olen2 > 0)
          olen2 = alen;
        else
          olen2 = alen+olen2;
        if (olen1 < olen2)
          which = (dif1 + ld_ratio*(olen2-olen1) >= dif2);
        else
          which = (dif2 + ld_ratio*(olen1-olen2) < dif1);
      }
    if (which)
      { pos1 = pos2; dif1 = dif2; }
  }

  /* Now that you have finish point, first compute complete alignment
     by calling OKNAlign or extrapolate start point as dictated by
     `what'.  Then build an overlap record modelling the overlap.   */

  OVL.endpos = pos1;
  OVL.diffs  = dif1;
  OVL.comp   = opposite;

  if (what == AS_FIND_ALIGN || what == AS_FIND_ALIGN_NO_TRACE)
    { OVL.trace  = AS_ALN_OKNAlign(aseq,alen,bseq,blen,&pos1,dif1);  
#define ELIM_TRAILING_SPACES
#ifdef ELIM_TRAILING_SPACES
      AS_ALN_clean_up_trace(OVL.trace,alen,blen,&pos1,&(OVL.endpos));
#endif
      OVL.begpos = pos1;
    }
  else if (what == AS_FIND_AFFINE_ALIGN)
    { OVL.trace  = AS_ALN_OKNAffine(aseq,alen,bseq,blen,&(OVL.begpos),&(OVL.endpos),dif1);
    } 
  else 
    { OVL.trace  = NULL;
      OVL.begpos = alen - (blen-OVL.endpos);
    }

  if(what == AS_FIND_ALIGN_NO_TRACE)
    { OVL.trace=NULL;
    }

  OVL.length = (alen + blen - (abs(OVL.begpos) + abs(OVL.endpos))) / 2;

  if (opposite)
    Complement_Seq(bseq+1);
  return (&OVL);
}

/* end of what was AS_CNS version of DP_Compare() */



/* Compute the longest tail path from point (imax,jmax) of the d.p. matrix
   assuming the score at that point is bpmax.  Only paths that are strictly
   positive at every vertex starting at (imax,jmax) are considered "alive".
   A tail path is one whose interior is alive and either reaches the border
   of the matrix or has a last vertex of score below 0.  The computation is
   sparse in that only the range of each row that contains a live path is
   retained in each row and the computation terminates if this range
   becomes empty.
*/

static double BPSuffixScore(char *a, int alen, char *b, int blen,
                            int imax, int jmax, double bpmax)
{ float *dc;
  int    lft, rgt;
  double bprat, bprem;
  int    longest;

  dc = DPMatrix;         /* Name the row space something short */

  bprat = BP_RATIO;      /* Setup the two possible edge scores */
  bprem = 1. - BP_RATIO;

  /* Initialize *entire* row to be 0.  For this routine, we maintain the
     invariant that any cells not in the [lft,rgt] active interval are
     zero or less.                                                     */

  { int i;
    for (i = 0; i <= alen; i++)
      dc[i] = 0.;
  }

  /* Compute the active part of the first row */

  { int i;

    dc[imax] = bpmax;
    lft = imax;
    for (i = imax+1; i <= alen; i++)
      { dc[i] = dc[i-1] - bprem;
        if (dc[i] <= 0.)          /* ==> reached end of active region */
          { longest = i-imax;
            rgt = i-1;
            break;
          }
      }
    if (i > alen)    /* ==> active region goes to end of matrix */
      { longest = alen-imax;
        rgt = alen;
      }
  }

#ifdef BP_TAILS
  { int i;

    printf("\nTAIL SEGMENTS A(%d,%d) x B(%d,%d)\n",imax,alen,jmax,blen);
    printf("           ");
    for (i = imax+1; i <= alen; i++)
      printf(" %5d",i);
    printf("\n");
    printf("           ");
    for (i = imax+1; i <= alen; i++)
      printf("     %c",a[i]);
    printf("\n");

    printf("     :");
    for (i = imax; i <= rgt; i++)
      { if (i < lft)
          printf("      ");
        if (dc[i] <= 0.)
          printf(" -----");
        else
          printf(" %5.2f",dc[i]);
      }
    printf("       >>> %d (%d,%d)\n",longest,lft,rgt);
  }
#endif

  /* Compute active parts of all subsequent rows till reach last row
     or active region becomes empty                                  */

  { int i, j, len;

    for (j = jmax+1; j <= blen; j++)
      { double c, e;
        int    posval = 0;

        /* Compute value in leftmost active column in new row */

        e = dc[lft];
        dc[lft] = c = e - bprem;
        if (c <= 0.)
          { if (lft - imax > j - jmax)
              len = lft-imax;
            else
              len = j-jmax;
            if (longest < len) longest = len;
          }

        /* Compute subsequent elements in new row until reach boundary
           of compute a value <= 0. beyond rgt boundary of previous row */

        for (i = lft+1; i <= alen; i++)
          { if (c > 0.)
              { posval = 1;        /* Note that a previous cell contributes  */
                c -= bprem;        /* to the current one only if it is still */
              }                    /* alive, i.e. >0.  If prior cell is      */
            if (e > 0.)            /* alive, then the cell is assigned a     */
              { posval = 1;        /* non-positive value and posval is zero. */
                if (a[i] == b[j])
                  e += bprat;
                else
                  e -= bprem;
              }
            if (c < e) c = e;
            e = dc[i];
            if (e > 0.)
              { posval = 1;
                if (e-bprem > c)
                  c = e-bprem;
              }
            dc[i] = c;
            if (c <= 0.)
              { if (posval)
                  { if (i - imax > j - jmax)
                      len = i-imax;
                    else
                      len = j-jmax;
                    if (longest < len) longest = len;
                  }
                if (i > rgt) break;  /* <= 0. and beyond rgt ==> quit */
              }
          }
        rgt = i-1;

        if (i > alen)   /* Ran off end of matrix ==> tail path */
          { if (alen - imax > j - jmax)
              len = alen-imax;
            else
              len = j-jmax;
            if (longest < len) longest = len;
          }

#ifdef BP_TAILS
        printf("%3d %c:",j,b[j]);
        for (i = imax; i <= rgt; i++)
          { if (i < lft)
              printf("      ");
            else if (dc[i] <= 0.)
              printf(" -----");
            else
              printf(" %5.2f",dc[i]);
          }
        printf("       >>> %d",longest);
#endif

        /* Adjust boundaries of active region by removing the <= 0
           prefix and suffix of the part of the row that was computed. */
  
        while (lft <= rgt && dc[lft] <= 0.)
          lft += 1;
        if (lft > rgt) break;   /* Active region is empty */
        while (dc[rgt] <= 0.)
          rgt -= 1;

#ifdef BP_TAILS
        printf(" (%d,%d)\n",lft,rgt);
#endif
      }

    /* Reached last row of matrix.  All non-positive entries in the
       active part of this row are tail paths.                        */

    if (j > jmax)
      { for (i = lft; i <= rgt; i++)
          if (dc[i] > 0.)
            { if (i - imax > blen - jmax)
                len = i-imax;
              else
                len = blen-jmax;
              if (longest < len) longest = len;
            }
      }
#ifdef BP_TAILS
    printf("\nFinal pass = %d\n",longest);
#endif
  }

  return (bpmax / longest);
}

/* Do a d.p. computation wherein matches score BP_RATIO and mismatches score
   1-BP_RATIO.  The boundary on the diagonals in [beg,end] are initialized
   to zero, all interior cells are computed by the basic recurrence.  Thus
   one has the score of the best path to the allowed part of the boundary
   in each cell.  The computation is further made sparse in that one only
   computes cells in the error-rate thresholded part of the matrix determined
   by a simultaneous unit-scoring computation.  The location of the peak
   value observed is taken as the branch point (if there is one).  Provided
   that this peak occurs in an alignment with at least prefix columns to
   the left and suffix columns to the right, a branchpoint is recorded and
   BPSuffixScore is called to determine the descent rate after the branch
   point.
*/

static BranchPointResult *BPDetector(char *a, int alen, char *b, int blen,
                                     int beg, int end, int prefix, int suffix)
{ int j, jfst;
  int lft, rgt;
  int ahg, bhg;
  int *dp;
  float *dc;
  double bprat, bprem, bpmax;
  int imax = 0, jmax = 0, dmax = 0; 

  static BranchPointResult Ans;

#ifdef BP_THRESH
  { int i;

    for (i = 0; i <= 1.5*alen; i++)
      printf("%3d: %3d\n",i,DistThresh[i]);
  }
#endif

  dp = HorzDelta;
  dc = DPMatrix;

  bprat = BP_RATIO;
  bprem = 1. - BP_RATIO;
  bpmax = 0.;

  if (end < 0)      /* Determine start row and # of 0-entries in that col. */
    { jfst = -end;
      ahg  = 0;
    }
  else
    { jfst = 0;
      ahg  = end;
    }

  /* Compute first row and its active region [lft,rgt] */
 
  if (jfst == 0)
    { int i;

      if (beg < 0)
        lft = 0;
      else
        lft = beg;
      rgt = end;
      for (i = lft; i <= rgt; i++)
        { dp[i] = 0;
          dc[i] = 0.;
        }
      for (i = 1; rgt+i <= alen && i < DistThresh[i]; i++)
        { dp[rgt+i] = i;
          dc[rgt+i] = -i*bprem;
        }
      rgt += (i-1);
    }
  else
    { int i;

      for (i = 0; i <= alen && i < DistThresh[i]; i++)
        { dp[i] = i;
          dc[i] = -i*bprem;
        }
      lft = 0;
      rgt = i-1;
    }

#ifdef BP_RATDP
   { int i;

#define BPFLEN  5
#define BPBOUND

     printf("\nSUBSEGMENT A(0,%d) x B(%d,%d)\n",alen,jfst,blen);
     printf("           ");
     for (i = 1; i <= alen; i++)
       printf(" %5d",i);
     printf("\n");
     printf("           ");
     for (i = 1; i <= alen; i++)
       printf("     %c",a[i]);
     printf("\n");
   }
#else
#if defined(BP_DEBUG) || defined(BP_THRESH)
   { int i;

#define BPFLEN 2
#define BPBOUND

     printf("\nSUBSEGMENT A(1,%d) x B(%d,%d)\n",alen,jfst,blen);
     printf("         ");
     for (i = 1; i <= alen; i++)
       printf(" %2d",i);
     printf("\n");
     printf("         ");
     for (i = 1; i <= alen; i++)
       printf("  %c",a[i]);
     printf("\n");
   }
#endif
#endif

#ifdef BP_DEBUG
   { int i;

     if (jfst > 0)
       printf("%3d %c:",jfst,b[jfst]);
     else
       printf("     :");
     for (i = 0; i <= rgt; i++)
       { if (i < lft)
           printf(" %*s",BPFLEN,"");
         else
           printf(" %*d",BPFLEN,dp[i]);
      }
    printf(" (%d,%d)\n",lft,rgt);
  }
#endif

#ifdef BP_RATDP
  { int i;

    if (jfst > 0)
      printf("%3d %c:",jfst,b[jfst]);
    else
      printf("     :");
    for (i = 0; i <= rgt; i++)
      { if (i < lft)
          printf(" %*s",BPFLEN,"");
        else
          printf(" %*.2f",BPFLEN,dc[i]);
      }
    printf("\n");
  }
#endif

  jfst += 1;
  if (beg < 0)
    bhg = jfst+beg;
  else
    bhg = jfst;

  for (j = jfst; j <= blen; j++)
    { { int e, c, i;
        float f, d;

        /* Compute lft-cell in current row */

        e = dp[lft];
        f = dc[lft];
        if (j <= -beg)
          { dp[lft] = c = 0;
            dc[lft] = d = 0.;
          }
        else
          { dp[lft] = c = e+1;
            d = f - bprem;
            if (d > bpmax)
              { bpmax = d;
                imax  = lft;
                jmax  = j;
                dmax  = c;
              } 
            dc[lft] = d;
          }

        /* Compute remaining cells corresponding to the columns of the
           active region of the previous row.                           */

        for (i = lft+1; i <= rgt; i++)
          { c += 1;
            d -= bprem;

            if (a[i] == b[j])
              f += bprat;
            else
              { f -= bprem;
                e += 1;
              }

            if (c > e) c = e;
            if (d < f) d = f;

            e = dp[i];
            if (e < c) c = e+1;
            f = dc[i];
            if (f-bprem > d) d = f-bprem;

            dp[i] = c;
            if (d > bpmax)
              { bpmax = d;
                imax  = i;
                jmax  = j;
                dmax  = c;
              } 
            dc[i] = d;
          }

        /* Do column rgt+1 */

        if (i <= alen)
          { c += 1;
            d -= bprem;

            if (a[i] == b[j])
              f += bprat;
            else
              { f -= bprem;
                e += 1;
              }

            if (c > e) c = e;
            if (d < f) d = f;

            dp[i] = c;
            if (d > bpmax)
              { bpmax = d;
                imax  = i;
                jmax  = j;
                dmax  = c;
              } 
            dc[i] = d;
            i += 1;
          }

        /* Do any remaining columns that are with error thresholded zone. */

        while (i <= alen)
          { c = c+1;
            if (i < j)
              if (i < bhg)
                { if (c >= DistThresh[bhg]) break; }
              else
                { if (c >= DistThresh[i]) break; }
            else
              if (i <= j+ahg)
                { if (c >= DistThresh[j]) break; }
              else
                { if (c >= DistThresh[i-ahg]) break; }

            dp[i] = c;
            d -= bprem;
            if (d > bpmax)
              { bpmax = d;
                imax  = i;
                jmax  = j;
                dmax  = c;
              } 
            dc[i] = d;
            i += 1;
          }

        rgt = i-1;
      }

#ifdef BPBOUND
    { int olft, orgt;

      olft = lft;
      orgt = rgt;
#endif

      /* Pull back left and right boundaries of active region
         based on error threshold.                             */

      while (1)
        { if (rgt < j)
            if (rgt < bhg)
              { if (dp[rgt] < DistThresh[bhg]) break; }
            else
              { if (dp[rgt] < DistThresh[rgt]) break; }
          else
            if (rgt <= j+ahg)
              { if (dp[rgt] < DistThresh[j]) break; }
            else
              { if (dp[rgt] < DistThresh[rgt-ahg]) break; }
          rgt -= 1;
        }

      if (rgt < lft) break;

      if (j > -beg)
        while (1)
          { if (lft < j)
              if (lft < bhg)
                { if (dp[lft] < DistThresh[bhg]) break; }
              else
                { if (dp[lft] < DistThresh[lft]) break; }
            else
              if (lft <= j+ahg)
                { if (dp[lft] < DistThresh[j]) break; }
              else
                { if (dp[lft] < DistThresh[lft-ahg]) break; }
            lft += 1;
          }

#ifdef BP_THRESH
      { int i;

        printf("%3d %c:",j,b[j]);
        for (i = 0; i <= orgt; i++)
          { if (i < olft)
              printf(" %*s",BPFLEN,"");
            else if (i < j)
              if (i < bhg)
                printf(" %*d",BPFLEN,DistThresh[bhg]);
              else
                printf(" %*d",BPFLEN,DistThresh[i]);
            else
              if (i <= j+ahg)
                printf(" %*d",BPFLEN,DistThresh[j]);
              else
                printf(" %*d",BPFLEN,DistThresh[i-ahg]);
          }
        printf("\n");
      }
#endif

#ifdef BP_RATDP
      { int i;

        printf("%3d %c:",j,b[j]);
        for (i = 0; i <= orgt; i++)
          { if (i < olft)
              printf(" %*s",BPFLEN,"");
            else
              printf(" %*.2f",BPFLEN,dc[i]);
          }
        printf("\n");
      }
#endif

#ifdef BP_DEBUG
      { int i;

        printf("%3d %c:",j,b[j]);
        for (i = 0; i <= orgt; i++)
          { if (i < olft)
              printf(" %*s",BPFLEN,"");
            else if (i == lft)
              printf(">%*d",BPFLEN,dp[i]);
            else if (i == rgt+1)
              printf("<%*d",BPFLEN,dp[i]);
            else
              printf(" %*d",BPFLEN,dp[i]);
          }
        if (i == rgt+1)
          printf("< (%d,%d)\n",lft,rgt);
        else
          printf(" (%d,%d)\n",lft,rgt);
      }
#endif

#ifdef BPBOUND
    }
#endif

      bhg += 1;
    }

#ifdef BP_RATDP
  printf("MAX(i,j) = %g(%d,%d)\n",bpmax,imax,jmax);
#endif

  { int plen, slen;

    if (imax < jmax)
      plen = imax + (int)(.5*dmax);
    else
      plen = jmax + (int)(.5*dmax);
    if (alen - imax < blen - jmax)
      slen = alen-imax;
    else
      slen = blen-jmax;
#ifdef BP_RATDP
  printf("  plen = %d slen = %d\n",plen,slen);
#endif
    if (bpmax <= 0. || plen < prefix || slen < suffix)
      Ans.apnt = Ans.bpnt = -1;
    else
      { Ans.apnt    = imax;
        Ans.bpnt    = jmax;
        Ans.ascent  = bpmax / plen;
        Ans.descent = BPSuffixScore(a,alen,b,blen,imax,jmax,bpmax);
      }
    return (&Ans);
  }
}

/* Given fragments a and b, determine if there is a branchpoint between
   them subject to the parameters beg, end, opposite, erate, and thresh as
   for DP_Compare_AS above.  Typically one should use twice the prevailing
   sequencing error rate for erate.  Moreover the branchpoint must occur
   at least minprefix symbols into the overlap and involve at least
   minsuffix symbols beyond the branch.  If there is insufficient memory
   then a NULL pointer is returned, otherwise a branchpoint record is
   returned.  The position of the branchpoint is returned as (-1,-1) if
   a branchpoint satisfying the given constraints was not found.  Otherwise
   a complete record as described with the typedef of the b.p. record
   is returned.  *This record is reused every time the routine is called.* */

BranchPointResult *BPnt_Compare_AS(InternalFragMesg *a, InternalFragMesg *b,
                                   int beg, int end, int opposite,
                                   double erate, double thresh,
                                   int minprefix, int minsuffix)
{ int   alen,  blen;
  char *aseq, *bseq;
  BranchPointResult *bpr;

  aseq = a->sequence;  /* Setup sequence access */
  bseq = b->sequence;
  alen = strlen(aseq);
  blen = strlen(bseq);
  aseq -= 1;
  bseq -= 1;

  if (Space_n_Tables(alen,erate,thresh)) return (NULL);

  if (opposite)                 /* Compare in opposite orientation. */
    Complement_Fragment_AS(b);

  bpr = BPDetector(aseq,alen,bseq,blen,beg,end,minprefix,minsuffix);

  if (opposite)                 /* Compare in opposite orientation. */
    Complement_Fragment_AS(b);

  return (bpr);
}



/* Identical to previous except that sequences  aseq  and  bseq  
*  have already been extracted, complemented if necessary,
*  and shifted so that the first characters
*  of each are  aseq [1]  and  bseq [1] , respectively.
*  Their lengths also have been determined to be  alen  and  blen . */

BranchPointResult *BPnt_Seq_Comp_AS(char * aseq, int alen,  char * bseq, int blen,
                                    int beg, int end, double erate, double thresh,
                                    int minprefix, int minsuffix)

{ if (Space_n_Tables(alen,erate,thresh)) return (NULL);
  return (BPDetector(aseq,alen,bseq,blen,beg,end,minprefix,minsuffix));
}





#define PRINT_WIDTH  50   /* Width of each line of a printed alignment */

/*** OVERLAP PRINT ROUTINE ***/

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

static void CNS_PrintAlign(FILE *file, int prefix, int suffix,
                       char *a, char *b, int *trace)
{ int i, j, o;
  static char Abuf[PRINT_WIDTH+1], Bbuf[PRINT_WIDTH+1];
  static int  Firstime = 1;

  if (Firstime)
    { Firstime = 0;
      Abuf[PRINT_WIDTH] = Bbuf[PRINT_WIDTH] = '\0';
    }
                                           /* buffer/output next column */
#define COLUMN(x,y)			\
{ if (o >= PRINT_WIDTH)			\
    { fprintf(file,"\n\t%s\n",Abuf);	\
      fprintf(file,"\t%s\n",Bbuf);	\
      o = 0;				\
    }					\
  Abuf[o] = (x);			\
  Bbuf[o] = (y);			\
  o += 1;				\
}

  a -= 1;
  b -= 1;
  o  = 0;
  i = j = 1;

  if (prefix > 25)
    { i = prefix-24;
      prefix = 25;
    }
  else if (prefix < -25)
    { j = (-prefix) - 24;
      prefix = -25;
    }

  if (prefix > 0)
    while (prefix-- > 0)     /* Output unaligned prefix */
      COLUMN(a[i++],' ')
  else
    while (prefix++ < 0)
      COLUMN(' ',b[j++])

  { int p, c;      /* Output columns of alignment til reach trace end */

    p = 0;
    while ((c = trace[p++]) != 0)
      if (c < 0)
        { c = -c;
          while (i != c)
            COLUMN(a[i++],b[j++])
          COLUMN('-',b[j++])
        }
      else
        { while (j != c)
            COLUMN(a[i++],b[j++])
          COLUMN(a[i++],'-')
        }
  }

  if (suffix < 0) suffix = -suffix;
  if (suffix > 25)
    suffix = 25;

  { int x, y, s;     /* Output remaining column including unaligned suffix */

    s = 0;
    y = 1;
    while ((x = a[i++]) != 0)
      { if ((y = b[j++]) != 0)
          COLUMN(x,y)
        else
          { do 
              { COLUMN(x,' ')
                s += 1;
              }
            while ((x = a[i++]) != 0 && s < suffix);
            break;
          }
      }
    if (y)
      while ((y = b[j++]) != 0 && s < suffix)
        { COLUMN(' ',y)
          s += 1;
        }
  }

  fprintf(file,"\n\t%.*s\n",o,Abuf);   /* Print remainder of buffered col.s */
  fprintf(file,"\t%.*s\n",o,Bbuf);
}

/* Print an ASCII representation of the overlap in align between fragments
   a and b to given file.                                                  */

void Print_Overlap(FILE *file, char *aseq, char *bseq, Overlap *align)
{ int sym1, sym2;

  if (align->comp)
    { sym1 = '<'; sym2 = '-'; }
  else
    { sym1 = '-'; sym2 = '>'; }
  if (align->begpos >= 0)
    if (align->endpos <= 0)
      { fprintf(file,"  A -----+------+---->");
        fprintf(file,"    dif/len = %d/%d = %5.2f%%\n",
                     align->diffs,align->length,
                     (100.*align->diffs)/align->length);
        fprintf(file,"  B %4d %c------%c %-4d\n",
                     align->begpos,sym1,sym2,-align->endpos);
      }
    else
      { fprintf(file,"  A -----+------> %-4d",align->endpos);
        fprintf(file,"       dif/len = %d/%d = %5.2f%%\n",
                     align->diffs,align->length,
                     (100.*align->diffs)/align->length);
        fprintf(file,"  B %4d %c------+----%c\n",align->begpos,sym1,sym2);
      }
  else
    if (align->endpos >= 0)
      { fprintf(file,"  A %4d -------> %-4d",
                     -align->begpos,align->endpos);
        fprintf(file,"       dif/len = %d/%d = %5.2f%%\n",
                     align->diffs,align->length,
                     (100.*align->diffs)/align->length);
        fprintf(file,"  B %c----+------+----%c\n",sym1,sym2);
      }
    else
      { fprintf(file,"  A %4d -------+----> A",-align->begpos);
        fprintf(file,"    dif/len = %d/%d = %5.2f%%\n",
                     align->diffs,align->length,
                     (100.*align->diffs)/align->length);
        fprintf(file,"  B %c----+------%c %-4d\n",sym1,sym2,-align->endpos);
      }
  fprintf(file,"\n");

  if (align->trace != NULL)
    { int *trace;

      trace = align->trace;

#ifdef DEBUG
      { int i;

        fprintf(file,"\nUncompressed trace:\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file,"  %3d\n",trace[i]);
      }
#endif

      if (align->comp)
        Complement_Seq(bseq);

      CNS_PrintAlign(file,align->begpos,align->endpos,aseq,bseq,trace);

      if (align->comp)
        Complement_Seq(bseq);
    }
} 

Overlap *Copy_Overlap(Overlap *ovl)
{ int      i, len;
  Overlap *new_ovl;

  new_ovl = (Overlap *) malloc(sizeof(Overlap));
  new_ovl->begpos = ovl->begpos;
  new_ovl->endpos = ovl->endpos;
  new_ovl->diffs  = ovl->diffs;
  new_ovl->length = ovl->length;
  new_ovl->comp   = ovl->comp;

  len = 0;
  while (ovl->trace[len] != 0)
    len += 1;

  new_ovl->trace = (int *) malloc(sizeof(int)*(len+1));

  for (i = 0; i <= len; i++)
    new_ovl->trace[i] = ovl->trace[i];

  return (new_ovl);
}

