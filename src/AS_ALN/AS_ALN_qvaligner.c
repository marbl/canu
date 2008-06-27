
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

/* Utility routines to complement, unpack and pack alignments, and print
   overlaps.  Also a routine for re-aligning an overlap using quality
   values based on a banding approach as described in "REAligner: A Program
   for Refining DNA Sequence Multi-Alignments," J. Comp. Bio. 4 (1997),
   369-383, by Eric Anson and Gene Myers. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_ALN_aligners.h"

#undef DEBUG

#define PRINT_WIDTH  50   /* Width of each line of a printed alignment */
#define BAND_WIDTH    3   /* Width of band about original for realignment */

/*** UTILITY ROUTINES ***/


/* Complement the sequence in fragment aseq.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Complement_Seq(char *aseq)
{ static char WCinvert[256];
  static int Firstime = 1;

  if (Firstime)          /* Setup complementation array */
    { int i;

      Firstime = 0;
      for(i = 0; i < 256;i++){
	WCinvert[i] = '?';
      }
      WCinvert['a'] = 't';
      WCinvert['c'] = 'g';
      WCinvert['g'] = 'c';
      WCinvert['t'] = 'a';
      WCinvert['n'] = 'n';
      WCinvert['A'] = 'T';
      WCinvert['C'] = 'G';
      WCinvert['G'] = 'C';
      WCinvert['T'] = 'A';
      WCinvert['N'] = 'N';
      WCinvert['-'] = '-'; // added this to enable alignment of gapped consensi
    }

  { int len;                    /* Complement and reverse sequence */
    len = strlen(aseq);

    { register char *s, *t;
      int c;

      s = aseq;
      t = aseq + (len-1);
      while (s < t)
	{ c = *s;
          *s++ = WCinvert[(int)*t];
          *t-- = WCinvert[(int)c];
        }
      if (s == t)
        *s = WCinvert[(int)*s];
    }
  }
}


/* Complement the sequence in fragment message a.  This include also
   revsersing the order of the quality values.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

void Complement_Fragment_AS(InternalFragMesg *a)
{ static char WCinvert[256];
  static int Firstime = 1;

  if (Firstime)          /* Setup complementation array */
    {
      int i;
      Firstime = 0;
      for(i = 0; i < 256;i++){
	WCinvert[i] = '?';
      }
      WCinvert['a'] = 't';
      WCinvert['c'] = 'g';
      WCinvert['g'] = 'c';
      WCinvert['t'] = 'a';
      WCinvert['n'] = 'n';
      WCinvert['A'] = 'T';
      WCinvert['C'] = 'G';
      WCinvert['G'] = 'C';
      WCinvert['T'] = 'A';
      WCinvert['N'] = 'N';
      WCinvert['-'] = '-'; // added this to enable alignment of gapped consensi
    }

  { int len;                    /* Complement and reverse sequence */
    len = strlen(a->sequence);

    { register char *s, *t;
      int c;

      s = a->sequence;
      t = a->sequence + (len-1);
      while (s < t)
        { // Sanity Check!
	  assert(WCinvert[(int) *t] != '?' &&
		 WCinvert[(int) *s] != '?');

	  c = *s;
          *s++ = WCinvert[(int) *t];
          *t-- = WCinvert[c];
        }
      if (s == t)
        *s = WCinvert[(int) *s];
    }

    if (a->quality != NULL)
      { register char *s, *t;   /* Reverse quality value array */
        int c;

        s = a->quality;
        t = a->quality + (len-1);
        while (s < t)
          { c = *s;
            *s++ = *t;
            *t-- = c;
          }
      }
  }
}

/* Convert an overlap message delta into an array of int that more
   directly encode the alignment.  For an unpacked trace < i1, i2, ... in, 0>
   a negative number j indicates that a dash should be placed before A[-j]
   and a positive number k indicates that a dash should be placed before
   B[k], where A and B are the two sequences of the overlap.  These indels
   occur in order along the alignment.

   A pointer to an array containing the unpacked trace is returned.  This
   array is owned by the routine and is reused by it with each subsequent
   call.  If the unpacked trace is needed beyond a subsequent call, the
   caller must copy its contents to a memory area they have allocated.   */

int *Unpack_Alignment_AS(OverlapMesg *align)
{ static int    *UnpackBuffer=NULL;
  static int buffalloc=0;
  int buffused=0;
  signed char  *calign;
  int           apos, bpos;

  if(UnpackBuffer==NULL){
    UnpackBuffer = (int*) safe_malloc (sizeof(int)*(2*AS_READ_MAX_LEN+1));
    buffalloc=(2*AS_READ_MAX_LEN+1);
  }

  apos   = align->ahg + 1;  /* ahg >= 0 for all overlaps */
  bpos   = 1;
  calign = align->delta;

  { int i, c;
    int *spt;

    spt = UnpackBuffer;
    for (i = 0; (c = calign[i]) != 0; i++)
      if (c == AS_LONG_DELTA_CODE)
        { apos += AS_LONGEST_DELTA;      /* Uninterrupted match of 126 bases */
          bpos += AS_LONGEST_DELTA;
        }
      else if (c == AS_POLY_DELTA_CODE)
        { c = calign[++i];  /* Multi-base gap */
          if (c < 0)
            while (c < 0)
              { c    += 1;
                bpos += 1;
		if(buffused>=buffalloc){
		  buffalloc=buffalloc*2+10;
		  UnpackBuffer = (int*)safe_realloc(UnpackBuffer,sizeof(int)*buffalloc);
		  spt=UnpackBuffer+buffused;
		}
		buffused++;
                *spt++ = -apos;
              }
          else
            while (c > 0)
              { c    -= 1;
                apos += 1;
		if(buffused>=buffalloc){
		  buffalloc=buffalloc*2+10;
		  UnpackBuffer = (int*)safe_realloc(UnpackBuffer,sizeof(int)*buffalloc);
		  spt=UnpackBuffer+buffused;
		}
		buffused++;
                *spt++ = bpos;
              }
        }
      else
        { if (c < 0)        /* Single gap */
            { bpos -= c;
              apos -= (c+1);
	      if(buffused>=buffalloc){
		buffalloc=buffalloc*2+10;
		UnpackBuffer = (int*)safe_realloc(UnpackBuffer,sizeof(int)*buffalloc);
		spt=UnpackBuffer+buffused;
	      }
	      buffused++;
              *spt++ = -apos;
            }
          else
            { bpos += (c-1);
              apos += c;
	      if(buffused>=buffalloc){
		buffalloc=buffalloc*2+10;
		UnpackBuffer = (int*)safe_realloc(UnpackBuffer,sizeof(int)*buffalloc);
		spt=UnpackBuffer+buffused;
	      }
	      buffused++;
              *spt++ = bpos;
            }
        }
    *spt = 0;
  }

  return (UnpackBuffer);
}

/*  Produce an overlap delta for an unpacked trace between two sequences
    A and B where the first prefix symbols of A are unaligned with B
    (prefix must be > 0).                                             */

signed char *Pack_Alignment_AS(int *trace, int prefix)
{ static signed char *PackBuffer = NULL;
  static int          PackSize   = -1;
  signed char  *spt;
  int           apos, bpos, i, c, size;

  size = 1;
  apos = prefix;
  bpos = 0;
  for (i = 0; (c = trace[i]) != 0; i++)
    { if (c < 0)
        { c = -c;
          while (c-apos > AS_LONGEST_DELTA)
            { size += 1;
              apos += AS_LONGEST_DELTA;
              bpos += AS_LONGEST_DELTA;
            }
          size += 1;
          bpos += c-apos;
          apos  = c-1;
        }
      else
        { while (c-bpos > AS_LONGEST_DELTA)
            { size += 1;
              apos += AS_LONGEST_DELTA;
              bpos += AS_LONGEST_DELTA;
            }
          size += 1;
          apos += c-bpos;
          bpos  = c-1;
        }
    }

  if (size > PackSize)
    { PackSize = (int)(1.4*size) + 1000;
      PackBuffer = (signed char *) safe_realloc(PackBuffer,PackSize);
    }

  spt  = PackBuffer;
  apos = prefix;
  bpos = 0;
  for (i = 0; (c = trace[i]) != 0; i++)
    { if (c < 0)
        { c = -c;
          while (c-apos > AS_LONGEST_DELTA)
            { *spt++ = AS_LONG_DELTA_CODE;
              apos += AS_LONGEST_DELTA;
              bpos += AS_LONGEST_DELTA;
            }
          *spt++ = -(c-apos);
          bpos += c-apos;
          apos  = c-1;
        }
      else
        { while (c-bpos > AS_LONGEST_DELTA)
            { *spt++ = AS_LONG_DELTA_CODE;
              apos += AS_LONGEST_DELTA;
              bpos += AS_LONGEST_DELTA;
            }
          *spt++ = (c-bpos);
          apos += c-bpos;
          bpos  = c-1;
        }
    }
  *spt = 0;

  return (PackBuffer);
}

/*** OVERLAP PRINT ROUTINE ***/

/* Print an alignment to file between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.  */

void PrintAlign(FILE *file, int prefix, int suffix,
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

  if (prefix > AS_READ_MAX_LEN)
    { i = prefix-24;
      prefix = 25;
    }

  while (prefix-- > 0)     /* Output unaligned prefix */
    COLUMN(a[i++],' ')

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
  if (suffix > AS_READ_MAX_LEN)
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

void Print_Overlap_AS(FILE *file, InternalFragMesg *a,
                                  InternalFragMesg *b, OverlapMesg *align)
{
  if (a->iaccession == align->bifrag)
    {
      InternalFragMesg *c;
      c = a;
      a = b;
      b = c;
    }

  fprintf(file,"\nOVERLAP BETWEEN");
  fprintf(file," A = (%s," F_IID ")",AS_UID_toString(a->eaccession),a->iaccession);
  fprintf(file," and");
  fprintf(file," B = (%s," F_IID ")",AS_UID_toString(b->eaccession),b->iaccession);
  fprintf(file,"\n\n");

  switch (align->orientation)
  { case AS_NORMAL:
      if (align->bhg <= 0)
        { fprintf(file,"  A -----+------+----> %-4d\n",align->bhg);
          fprintf(file,"    %4d -------> B\n",align->ahg);
        }
      else
        { fprintf(file,"  A -----+------> %-4d\n",align->bhg);
          fprintf(file,"    %4d -------+----> B\n",align->ahg);
        }
      break;
    case AS_INNIE:
      if (align->bhg <= 0)
        { fprintf(file,"  A -----+------+----> %-4d\n",align->bhg);
          fprintf(file,"    %4d <------- B\n",align->ahg);
        }
      else
        { fprintf(file,"  A -----+------> %-4d\n",align->bhg);
          fprintf(file,"    %4d <------+----- B\n",align->ahg);
        }
      break;
    case AS_OUTTIE:
      if (align->bhg <= 0)
        { fprintf(file,"  A <----+------+----- %-4d\n",align->bhg);
          fprintf(file,"    %4d -------> B\n",align->ahg);
        }
      else
        { fprintf(file,"  A <----+------- %-4d\n",align->bhg);
          fprintf(file,"    %4d -------+----> B\n",align->ahg);
        }
      break;
    case AS_ANTI:
    case AS_UNKNOWN:
    default:
      assert(0);
      break;
  }
  fprintf(file,"\n");

  if (align->delta != NULL)
    { int *trace;

      trace = Unpack_Alignment_AS(align);

#ifdef DEBUG
      { int i;

        fprintf(stderr,"Print_Overlap_AS()-- uncompressed trace\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file," %3d",trace[i]);
        fprintf(stderr, "\n");
      }
#endif

      if (align->orientation == AS_INNIE)
        Complement_Fragment_AS(b);
      else if (align->orientation == AS_OUTTIE)
        Complement_Fragment_AS(a);

      PrintAlign(file,align->ahg,align->bhg,a->sequence,b->sequence,trace);

      if (align->orientation == AS_INNIE)
        Complement_Fragment_AS(b);
      else if (align->orientation == AS_OUTTIE)
        Complement_Fragment_AS(a);
    }
}



// Analyze an alignment between a and b given in trace (unpacked).
// Prefix gives the length of the initial prefix of a that is unaligned.
// Suffix gives the length of the bhang: positive if there is a suffix of
//    b that is unaligned, negative if there is a suffix of a instead.
// Returns a -1 terminated list of the positions in the a sequences at
// which errors of type amode occur, as well as:
//   alen - # of a symbols in overlap,
//   blen - # of b symbols in overlap,
//   del  - # of unaligned symbols in a,
//   sub  - # of substitutions,
//   ins  - # of unaligned symbols in b,
//   affdel - # of runs of unaligned symbols in a,
//   affins - # of runs of unaligned symbols in b.
//   blockdel - # of runs of size > blocksize of unaligned symbols in a,
//   blockins - # of runs of size > blocksize of unaligned symbols in b,
//   const blocksize - min length of an indel to count as a block.
//   optional biggestBlock - the length of the largest mismatch
//
void
AnalyzeAffineAlign(int prefix, int suffix,
                   char *a, char *b, int *trace, int amode,
                   int *alen, int *blen, int *del, int *sub, int *ins,
                   int *affdel, int *affins,
                   int *blockdel, int *blockins, int blocksize,
                   int *biggestBlock) {
  int i, j, imax, jmax;
  int inserts, deletes, subtit;
  int affinserts, affdeletes, lastindel,blockcount;
  int alength, blength;

  if (biggestBlock)
    *biggestBlock =0;

  imax = strlen(a);
  jmax = strlen(b);

  a -= 1;
  b -= 1;
  i  = prefix+1;
  j  = 1;

  alength = 0;
  blength = 0;
  deletes = 0;
  subtit  = 0;
  inserts = 0;
  affinserts = 0;
  affdeletes = 0;
  blockcount=0;
  lastindel = 0;
  *blockins=0;
  *blockdel=0;

  {
    int p, c;      /* Output columns of alignment til reach trace end */

    p = 0;
    while ((c = trace[p++]) != 0) {
      if (c < 0) {
        if(c!=lastindel) {
          affinserts++;
          blockcount=1;
        } else {
          (*blockins)+=(++blockcount == blocksize);
          if ((biggestBlock) && (blockcount > *biggestBlock))
              *biggestBlock = blockcount;
        }
        lastindel=c;
        c = -c;
        while (i < c) {
          if (a[i++] != b[j++])
            subtit += 1;
          alength += 1;
          blength += 1;
        }
        inserts += 1;
        blength += 1;
        j += 1;
      } else {
        if(c!=lastindel) {
          affdeletes++;
          blockcount=1;
        } else {
          (*blockdel)+=(++blockcount == blocksize);
          if ((biggestBlock) && (blockcount > *biggestBlock))
              *biggestBlock = blockcount;
        }
        lastindel=c;
        while (j < c) {
          if (a[i++] != b[j++])
            subtit += 1;
          alength += 1;
          blength += 1;
        }
        deletes += 1;
        alength += 1;
        i += 1;
      }
    }
  }

  if ((i < imax) && (j < jmax)) {
    int x, y;     /* Output remaining column including unaligned suffix */

    while ((x = a[i++]) != 0) {
      if ((y = b[j++]) != 0){
        if (x != y)
          subtit += 1;
        alength += 1;
        blength += 1;
      } else
        break;
    }
  }

  *alen = alength;
  *blen = blength;
  *del  = deletes;
  *sub  = subtit;
  *ins  = inserts;
  *affdel = affdeletes;
  *affins = affinserts;
}






// Analyze the overlap between fragments a and b.
//     alen - # of a symbols in overlap,
//     blen - # of b symbols in overlap,
//     del  - # of unaligned symbols in a,
//     sub  - # of substitutions,
//     ins  - # of unaligned symbols in b,
//     affdel - # of runs of unaligned symbols in a,
//     affins - # of runs of unaligned symbols in b.
//
void Analyze_Affine_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                               OverlapMesg *align, int amode,
                               int *alen, int *blen, int *del, int *sub, int *ins,
                               int *affdel, int *affins,
                               int *blockdel, int *blockins, int blocksize,
                               int *biggestBlock) {
  int swap = 0;

  if (a->iaccession == align->bifrag) {
    InternalFragMesg *c;
    c = a;
    a = b;
    b = c;
    swap = 1;
  }

  if (align->delta != NULL) {
    int *trace = Unpack_Alignment_AS(align);

    if (align->orientation == AS_INNIE)
      Complement_Fragment_AS(b);
    else if (align->orientation == AS_OUTTIE)
      Complement_Fragment_AS(a);

    AnalyzeAffineAlign(align->ahg,align->bhg,
                       a->sequence,b->sequence,trace,amode,
                       alen,blen,del,sub,ins,
                       affdel,affins,
                       blockdel,blockins,blocksize,
                       biggestBlock);

    if (align->orientation == AS_INNIE)
      Complement_Fragment_AS(b);
    else if (align->orientation == AS_OUTTIE)
      Complement_Fragment_AS(a);
  }

  if (swap) {
    swap = *alen;
    *alen = *blen;
    *blen = swap;

    swap = *del;
    *del = *ins;
    *ins = swap;

    swap = *affdel;
    *affdel = *affins;
    *affins = swap;
  }
}






/***** Wrapper for bubble smoothing overlap detector on top of affine dp_compare *****/

/*
   Usage/arguments as for DP_Compare_AS, except that some of the parameters
   are hijacked:
   - "erate" is used as a filter after the fact, but not passed to DP_C*_AS
   - "what" is overridden
   - "minlen" is passed to DP_Compare_AS, but it is also used as a filter
     to require that minlen matches were found

   Procedure:
   - call DP_Compare_AS, using maximal erate and what=AS_FIND_ALIGN;
            (should use affine alignment option when Gene checks it in).
   - use Analyze_Affine_Overlap_AS() to evaluate the resulting alignment.
   - if the error rate in the alignment, adjusted for an affine scoring scheme,
     is better than user's erate, return the OverlapMesg*; else return NULL;
     optionally, also test number of large indels.

   Assumptions/caveats:
   - as with DP_Compare_AS, the returned message must be copied if it is to
   be retained.
   - this version does not test for sanity with respect to placement of gaps;
     in principle, a gap of several hundred bases at the end (most likely
     indicating a true branchpoint) would be accepted if proposed; this
     seems to be safe enough since DP_Compare_AS doesn't seem to find overlaps
     above about 12% simple (non-affine) error rate.  However, more minor
     versions of this could cause (false) overlaps of shallow branchpoints
     (true branchpoints occurring near the ends of fragments)

       Branchpoint:

               .........+++++++          ("." matches; +:# mismatch)
	       .........#######

          could be treated as overlapping with an affine gap:

               .........+++++++
	       .........-------######    ("-" a gap)

     Equally, a (short) bad fragment end has a better change of being
     overlapped if affine gaps allow it to find the best match within a
     modest window of uncertainty:

         Insufficient trimming of the second fragment:

               .....ACAGTAGACGAGATAGGATAGATAGAGTAGACAGATAGTTGACTAAC
	            ||||||||||||||||||||||
               .....ACAGTAGACGAGATAGGATAGACAGTTA

	 could be interpreted as:

               .....ACAGTAGACGAGATAGGATAGATAGAGTAGACAGATAGTTGACTAAC
	            ||||||||||||||||||||||         ||| ||
               .....ACAGTAGACGAGATAGGATAGA---------CAGTTA



*/


/* The following are given as external to allow for "expert" consumers to
   modify the values without changing the function's argument list
   [in order to maintain the parallelism with the DP_Compare_AS
   argument list].  Thus, for instance, if an initial attempt fails to
   find an overlap required by consensus, the criteria could be relaxed
   in a second call.

   Users who take advantage of this should probably be sure to restore
   the default values each time they modify them; see, e.g.,
   Local_Overlap_AS_forCNS() */

// larger erate than normal to encourage finding overlaps with large indels
double MAXDPERATE=.20;
// boolean test
int AS_ALN_TEST_NUM_INDELS=1;
// size of indel to count as "large" when testing number of large indels
int AFFINEBLOCKSIZE= 4;
// number of large indels allowed
int AFFINE_MAX_BLOCKS=3;

OverlapMesg *AS_ALN_affine_overlap(InternalFragMesg *a, InternalFragMesg *b,
                           int beg, int end, int opposite,
                           double erate, double thresh, int minlen,
                           CompareOptions what, int *where)
{

  OverlapMesg *O;

  { double my_ld_ratio=.25;
    O = LD_DP_Compare_AS(a,b,beg,end,opposite,MAX(MAXDPERATE,erate),thresh,minlen,AS_FIND_AFFINE_ALIGN,where,my_ld_ratio);
  }

  if (O != NULL){

    int del, sub, ins, affdel, affins, alen, blen, blockdel, blockins;
    float errRate, errRateAffine;

    Analyze_Affine_Overlap_AS(a,b,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
			      &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, NULL);

    errRate = (sub+ins+del)/(double)(alen+ins);

    errRateAffine = (sub+affins+affdel)/
                    (double)(alen+ins-(del-affdel+ins-affins));

    O->quality=errRateAffine;
#define AFFINE_OVERLAP_DEBUG 0
    #if AFFINE_OVERLAP_DEBUG > 0
      printf("Alen %d, Blen %d, del %d, sub %d, ins %d\n"
	     " affdel %d, affins %d, blockdel %d, blockins %d\n",
	     alen,blen,del,sub,ins,
	     affdel,affins,blockdel,blockins);
      printf("Simple mismatch rate %f\n",errRate);
      printf("Affine mismatch rate %f\n",errRateAffine);
      #if AFFINE_OVERLAP_DEBUG > 1
      Print_Overlap_AS(stdout,a,b,O);
      #endif
    #endif

    if(errRateAffine<=erate&&
       alen>=minlen&&blen>=minlen&&
       ( (!AS_ALN_TEST_NUM_INDELS) || blockins+blockdel<=AFFINE_MAX_BLOCKS )  )
      {
	return(O);
      }

    #if AFFINE_OVERLAP_DEBUG > 1
    printf("Affine overlap found but failed quality tests\n");
    #endif

  }

#if AFFINE_OVERLAP_DEBUG > 3
  printf("No affine overlap found\n");
#endif
  return(NULL);

}




void Compute_Olap_Version(InternalFragMesg* a,InternalFragMesg *b,OverlapMesg *O,int *ahang,int *bhang, char *ori){

  if (a->iaccession == O->bifrag)
    {
      InternalFragMesg *c;
      c = a;
      a = b;
      b = c;
      *ahang = -O->ahg;
      *bhang = -O->bhg;
    } else {
      *ahang = O->ahg;
      *bhang = O->bhg;
    }

  *ori = (( O->orientation == AS_INNIE || O->orientation == AS_OUTTIE) ?  'I' : 'N');

  return;
}



/*  AS_ALN_clean_up_trace removes leading and trailing gaps */

void AS_ALN_clean_up_trace(int *trace,int alen, int blen,int *spos,int *epos){
  { int i=0;
    int j=0;
    int changeahang=0;
    int changebhang=0;
    char c;
    //printf("Trace (lens %d %d):",alen,blen);
    while(trace[i]!=0){
      c='*';
      if(trace[i]<-alen){
	changebhang++;
      } else if (trace[i]>blen){
	changebhang--;
      } else if (trace[i]==-1){
	changeahang--;
      } else if (trace[i]==1){
	changeahang++;
      } else {
	c=' ';
	trace[j++]=trace[i];
      }
      //printf(" %c%d",c,trace[i]);
      i++;
    }
    //printf("\n");
    trace[j]=0;
    *spos+=changeahang;
    *epos+=changebhang;
  }
}

