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

#include "GF_ALN_local.h"

/* O(kn) identity-based alignment algorithm.  Find alignment between
   a and b (of lengths alen and blen), that begins at finishing
   boundary position *spnt.  Return at *spnt the diagonal at which the
   alignment starts.                                                   */

int *AS_ALN_OKNAlign(const char *a, int alen, const char *b, int blen, int *spnt, int diff) {
  int diag, wpos, level;
  int fcell, infinity;

  static int    Wtop = -1;
  static int   *Wave;
  static int   *TraceBuffer;

  if (diff >= Wtop)        /* Space for diff wave? */
    { int max, del, *newp;

      max = (int)(1.2*diff) + 50;
      del = (max+5)*(max+1);
      //fprintf(stderr,"DP_COMPARE (AS_ALN_OKNAlign): reallocing %ld bytes\n",del*sizeof(int)+(max+1)*sizeof(int));
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
  }

  /* Compute waves 1 through d-1 do, each wave has
     two boundary cells at each of its ends.       */

  { int m, n, k;

    m = 5;
    n = 0;
    for (level = 1; 1; level++)
      { Wave[m++] = infinity;
        Wave[m++] = infinity;
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
            Wave[m++] = j;
            n += 1;
          }
        Wave[m++] = infinity;
        Wave[m++] = infinity;
        n += 1;
      }

madeit:


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
            k -= 1;
          }
        else if (m > n)
          { TraceBuffer[t++] = j+1;
            k += 1;
          }
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
