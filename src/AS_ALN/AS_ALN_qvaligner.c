
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
{ static int    UnpackBuffer[2*AS_READ_MAX_LEN+1];
  signed char  *calign;
  int           apos, bpos;

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
                *spt++ = -apos;
              }
          else
            while (c > 0)
              { c    -= 1;
                apos += 1;
                *spt++ = -bpos;
              }
        }
      else
        { if (c < 0)        /* Single gap */
            { bpos -= c;
              apos -= (c+1);
              *spt++ = -apos;
            }
          else
            { bpos += (c-1);
              apos += c;
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
      PackBuffer = (signed char *) realloc(PackBuffer,PackSize);
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
  fprintf(file," A = (" F_UID "," F_IID ")",a->eaccession,a->iaccession);
  fprintf(file," and");
  fprintf(file," B = (" F_UID "," F_IID ")",b->eaccession,b->iaccession);
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

        fprintf(file,"\nUncompressed trace:\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file,"  %3d\n",trace[i]);
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

/*** OVERLAP ANALYSIS ROUTINE ***/

/* Analyze an alignment between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.
   Returns a -1 terminated list of the positions in the a sequences at
   which errors of type amode occur, as well as:
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b.
*/

static int *AnalyzeAlign(int prefix, int suffix,
			 char *a, char *b, int *trace, int amode,
                         int *alen, int *blen, int *del, int *sub, int *ins)
{ int i, j, oa, ob;
  int inserts, deletes, subtit;
  int alength, blength;
  int *amistake, *bmistake;
  int dodels, doins, dosubs;
  static int *mistake = NULL;
  static int mislen  = -1;

  alength = strlen(a);
  blength = strlen(b);
  alength += 1;
  blength += 1;
  if (mislen < alength+blength)
    { mislen = (int)(2*(alength+blength)) + 4;
      mistake = (int *)realloc(mistake,mislen*sizeof(int));
    } 
  amistake = mistake;
  bmistake = mistake + alength;

  dodels = doins = dosubs = 0;
  if (amode == AS_ANALYZE_ALL)
    dodels = doins = dosubs = 1;
  else if (amode == AS_ANALYZE_DELETES)
    dodels = 1;
  else if (amode == AS_ANALYZE_INSERTS)
    doins = 1;
  else
    dosubs = 1;
  
  a -= 1;
  b -= 1;
  oa = 0;
  ob = 0;
  i  = prefix+1;
  j  = 1;

  alength = 0;
  blength = 0;
  deletes = 0;
  subtit  = 0;
  inserts = 0;

  { int p, c;      /* Output columns of alignment til reach trace end */

    p = 0;
    while ((c = trace[p++]) != 0)
      if (c < 0)
        { c = -c;
          while (i != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          inserts += 1;
          if (dodels)
            amistake[oa++] = i;
          if (doins)
            bmistake[ob++] = j;
          blength += 1;
          j += 1;
        }
      else
        { while (j != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          deletes += 1;
          if (doins)
            amistake[oa++] = i;
          if (dodels)
            bmistake[ob++] = j;
          alength += 1;
          i += 1;
        }

  }

  { int x, y;     /* Output remaining column including unaligned suffix */

    while ((x = a[i++]) != 0)
      { if ((y = b[j++]) != 0)
	{ if (x != y)
            { subtit += 1;
              if (dosubs)
	      { amistake[oa++] = i-1;
	        bmistake[ob++] = j-1;
	      }
	    }
	  alength += 1;
	  blength += 1;
        }else{
          break;
        }
      }
  }

  amistake[oa] = -1;
  bmistake[ob] = -1;

  *alen = alength;
  *blen = blength;
  *del  = deletes;
  *sub  = subtit;
  *ins  = inserts;
  
  return (mistake);
}


/* Analyze an alignment between a and b given in trace (unpacked).
   Prefix gives the length of the initial prefix of a that is unaligned.
   Suffix gives the length of the bhang: positive if there is a suffix of
      b that is unaligned, negative if there is a suffix of a instead.
   Returns a -1 terminated list of the positions in the a sequences at
   which errors of type amode occur, as well as:
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b,
     affdel - # of runs of unaligned symbols in a,
     affins - # of runs of unaligned symbols in b.
     blockdel - # of runs of size > blocksize of unaligned symbols in a,
     blockins - # of runs of size > blocksize of unaligned symbols in b,
     const blocksize - min length of an indel to count as a block.

*/

static int *AnalyzeAffineAlign(int prefix, int suffix, 
			 char *a, char *b, int *trace, int amode,
                         int *alen, int *blen, int *del, int *sub, int *ins,
			 int *affdel, int *affins,
			 int *blockdel, int *blockins, int blocksize)
{ int i, j, oa, ob;
  int inserts, deletes, subtit;
  int affinserts, affdeletes, lastindel,blockcount;
  int alength, blength;
  int *amistake, *bmistake;
  int dodels, doins, dosubs;
  static int *mistake = NULL;
  static int mislen  = -1;

  alength = strlen(a);
  blength = strlen(b);
  alength += 1;
  blength += 1;
  if (mislen < alength+blength)
    { mislen = (int)(1.4*(alength+blength)) + 500;
      mistake = (int *) realloc(mistake,mislen*sizeof(int));
    } 
  amistake = mistake;
  bmistake = mistake + alength;

  dodels = doins = dosubs = 0;
  if (amode == AS_ANALYZE_ALL)
    dodels = doins = dosubs = 1;
  else if (amode == AS_ANALYZE_DELETES)
    dodels = 1;
  else if (amode == AS_ANALYZE_INSERTS)
    doins = 1;
  else
    dosubs = 1;
  
  a -= 1;
  b -= 1;
  oa = 0;
  ob = 0;
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

  { int p, c;      /* Output columns of alignment til reach trace end */

    p = 0;
    while ((c = trace[p++]) != 0)
      if (c < 0)
        { if(c!=lastindel)
	    {
	      affinserts++;
	      blockcount=1;
	    } else {
	      (*blockins)+=(++blockcount == blocksize);
	    }
	  lastindel=c;
	  c = -c;
          while (i != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          inserts += 1;
          if (dodels)
            amistake[oa++] = i;
          if (doins)
            bmistake[ob++] = j;
          blength += 1;
          j += 1;
        }
      else
        { if(c!=lastindel)
	    { affdeletes++;
	      blockcount=1;
  	    } else {
	      (*blockdel)+=(++blockcount == blocksize);
	    }
	  lastindel=c;
	  while (j != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          deletes += 1;
          if (doins)
            amistake[oa++] = i;
          if (dodels)
            bmistake[ob++] = j;
          alength += 1;
          i += 1;
        }

  }

  { int x, y;     /* Output remaining column including unaligned suffix */

    while ((x = a[i++]) != 0)
      { if ((y = b[j++]) != 0){
          if (x != y)
            { subtit += 1;
              if (dosubs)
                { amistake[oa++] = i-1;
                  bmistake[ob++] = j-1;
                }
            }
	  alength += 1;
	  blength += 1;
      } else
          break;
      }
  }

  amistake[oa] = -1;
  bmistake[ob] = -1;

  *alen = alength;
  *blen = blength;
  *del  = deletes;
  *sub  = subtit;
  *ins  = inserts;
  *affdel = affdeletes;
  *affins = affinserts;

  return (mistake);
}


static int *AnalyzeAffineAlign_MaxIndelSize(int prefix, int suffix, 
			 char *a, char *b, int *trace, int amode,
                         int *alen, int *blen, int *del, int *sub, int *ins,
			 int *affdel, int *affins,
			 int *blockdel, int *blockins, int blocksize,
			 int *biggestBlock)
{ int i, j, oa, ob;
  int inserts, deletes, subtit;
  int affinserts, affdeletes, lastindel,blockcount;
  int alength, blength;
  int *amistake, *bmistake;
  int dodels, doins, dosubs;
  static int *mistake = NULL;
  static int mislen  = -1;

  alength = strlen(a);
  blength = strlen(b);
  alength += 1;
  blength += 1;
  if (mislen < alength+blength)
    { mislen = (int)(1.4*(alength+blength)) + 500;
      mistake = (int *) realloc(mistake,mislen*sizeof(int));
    } 
  amistake = mistake;
  bmistake = mistake + alength;


  *biggestBlock =0;
  dodels = doins = dosubs = 0;
  if (amode == AS_ANALYZE_ALL)
    dodels = doins = dosubs = 1;
  else if (amode == AS_ANALYZE_DELETES)
    dodels = 1;
  else if (amode == AS_ANALYZE_INSERTS)
    doins = 1;
  else
    dosubs = 1;
  
  a -= 1;
  b -= 1;
  oa = 0;
  ob = 0;
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

  { int p, c;      /* Output columns of alignment til reach trace end */

    p = 0;
    while ((c = trace[p++]) != 0)
      if (c < 0)
        { if(c!=lastindel)
	    {
	      affinserts++;
	      blockcount=1;
	    } else {
	      (*blockins)+=(++blockcount == blocksize);
              if(blockcount>*biggestBlock)*biggestBlock=blockcount;
	    }
	  lastindel=c;
	  c = -c;
          while (i != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          inserts += 1;
          if (dodels)
            amistake[oa++] = i;
          if (doins)
            bmistake[ob++] = j;
          blength += 1;
          j += 1;
        }
      else
        { if(c!=lastindel)
	    { affdeletes++;
	      blockcount=1;
  	    } else {
	      (*blockdel)+=(++blockcount == blocksize);
              if(blockcount>*biggestBlock)*biggestBlock=blockcount;
	    }
	  lastindel=c;
	  while (j != c)
            { if (a[i++] != b[j++])
                { subtit += 1;
                  if (dosubs)
                    { amistake[oa++] = i-1;
                      bmistake[ob++] = j-1;
                    }
                }
              alength += 1;
              blength += 1;
            }
          deletes += 1;
          if (doins)
            amistake[oa++] = i;
          if (dodels)
            bmistake[ob++] = j;
          alength += 1;
          i += 1;
        }

  }

  { int x, y;     /* Output remaining column including unaligned suffix */

    while ((x = a[i++]) != 0)
      { if ((y = b[j++]) != 0){
          if (x != y)
            { subtit += 1;
              if (dosubs)
                { amistake[oa++] = i-1;
                  bmistake[ob++] = j-1;
                }
            }
	  alength += 1;
	  blength += 1;
      } else
          break;
      }
  }

  amistake[oa] = -1;
  bmistake[ob] = -1;

  *alen = alength;
  *blen = blength;
  *del  = deletes;
  *sub  = subtit;
  *ins  = inserts;
  *affdel = affdeletes;
  *affins = affinserts;

  return (mistake);
}


/* Analyze the overlap between fragments a and b.
   Returns a -1 terminated list of the positions in the a sequences at
   which errors of type amode occur, as well as:
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b.
*/

int *Analyze_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                        OverlapMesg *align, int amode,
                        int *alen, int *blen, int *del, int *sub, int *ins)
{ int *mistakes = NULL;
  int swap;

  swap = 0;
  if (a->iaccession == align->bifrag)
    { InternalFragMesg *c;
      c = a;
      a = b;
      b = c;
      swap = 1;
    }
 
  if (align->delta != NULL)
    { int *trace;

      trace = Unpack_Alignment_AS(align);

#ifdef DEBUG
      { int i;

        fprintf(file,"\nUncompressed trace:\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file,"  %3d\n",trace[i]);
      }
#endif

      if (align->orientation == AS_INNIE)
        Complement_Fragment_AS(b);
      else if (align->orientation == AS_OUTTIE)
        Complement_Fragment_AS(a);

      mistakes = AnalyzeAlign(align->ahg,align->bhg,
			      a->sequence,b->sequence,trace,amode,
                              alen,blen,del,sub,ins);

      if (align->orientation == AS_INNIE)
        { int alength, blength, i, j, c;
          Complement_Fragment_AS(b);
          alength = strlen(a->sequence)+1;
          blength = strlen(b->sequence)+1;
          for (i = alength; mistakes[i] >= 0; i++)
            mistakes[i] = blength - mistakes[i];
          j = alength;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
      else if (align->orientation == AS_OUTTIE)
        { int alength, i, j, c;
          Complement_Fragment_AS(a);
          alength = strlen(a->sequence)+1;
          for (i = 0; mistakes[i] >= 0; i++)
            mistakes[i] = alength - mistakes[i];
          j = 0;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
    }

  if (swap)
    { int x;
      x = *alen;
      *alen = *blen;
      *blen = x;
      x = *del;
      *del = *ins;
      *ins = x;
      return (mistakes + strlen(a->sequence) + 1);
    }
  return (mistakes);
} 




/* Analyze the overlap between fragments a and b.
   Returns a -1 terminated list of the positions in the a sequences at
   which errors of type amode occur, as well as:
     alen - # of a symbols in overlap,
     blen - # of b symbols in overlap,
     del  - # of unaligned symbols in a,
     sub  - # of substitutions,
     ins  - # of unaligned symbols in b,
     affdel - # of runs of unaligned symbols in a,
     affins - # of runs of unaligned symbols in b.
     
*/

int *Analyze_Affine_Overlap_AS(InternalFragMesg *a, InternalFragMesg *b,
                        OverlapMesg *align, int amode,
                        int *alen, int *blen, int *del, int *sub, int *ins,
                        int *affdel, int *affins, 
			int *blockdel, int *blockins, int blocksize)
{ int *mistakes = NULL;
  int swap;

  swap = 0;
  if (a->iaccession == align->bifrag)
    { InternalFragMesg *c;
      c = a;
      a = b;
      b = c;
      swap = 1;
    }
 
  if (align->delta != NULL)
    { int *trace;

      trace = Unpack_Alignment_AS(align);

#ifdef DEBUG
      { int i;

        fprintf(file,"\nUncompressed trace:\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file,"  %3d\n",trace[i]);
      }
#endif

      if (align->orientation == AS_INNIE)
        Complement_Fragment_AS(b);
      else if (align->orientation == AS_OUTTIE)
        Complement_Fragment_AS(a);

      mistakes = AnalyzeAffineAlign(align->ahg,align->bhg,
				    a->sequence,b->sequence,trace,amode,
				    alen,blen,del,sub,ins,
				    affdel,affins,
				    blockdel,blockins,blocksize);

      if (align->orientation == AS_INNIE)
        { int alength, blength, i, j, c;
          Complement_Fragment_AS(b);
          alength = strlen(a->sequence)+1;
          blength = strlen(b->sequence)+1;
          for (i = alength; mistakes[i] >= 0; i++)
            mistakes[i] = blength - mistakes[i];
          j = alength;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
      else if (align->orientation == AS_OUTTIE)
        { int alength, i, j, c;
          Complement_Fragment_AS(a);
          alength = strlen(a->sequence)+1;
          for (i = 0; mistakes[i] >= 0; i++)
            mistakes[i] = alength - mistakes[i];
          j = 0;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
    }

  if (swap)
    { int x;
      x = *alen;
      *alen = *blen;
      *blen = x;
      x = *del;
      *del = *ins;
      *ins = x;
      x = *affdel;
      *affdel = *affins;
      *affins = x;
      return (mistakes + strlen(a->sequence) + 1);
    }
  return (mistakes);
} 


/* Analyze_Affine_Overlap_IndelSize_AS: 
   Just like Analyze_Affin_Overlap_AS with addition of evaluation of biggest block mismatch */

int *Analyze_Affine_Overlap_IndelSize_AS(InternalFragMesg *a, InternalFragMesg *b,
                        OverlapMesg *align, int amode,
                        int *alen, int *blen, int *del, int *sub, int *ins,
			int *affdel, int *affins, 
			int *blockdel, int *blockins, int blocksize,int *biggestBlock)
{ int *mistakes = NULL;
  int swap;

  swap = 0;
  if (a->iaccession == align->bifrag)
    { InternalFragMesg *c;
      c = a;
      a = b;
      b = c;
      swap = 1;
    }
 
  if (align->delta != NULL)
    { int *trace;

      trace = Unpack_Alignment_AS(align);

#ifdef DEBUG
      { int i;

        fprintf(file,"\nUncompressed trace:\n");
        for (i = 0; trace[i] != 0; i++)
          fprintf(file,"  %3d\n",trace[i]);
      }
#endif

      if (align->orientation == AS_INNIE)
        Complement_Fragment_AS(b);
      else if (align->orientation == AS_OUTTIE)
        Complement_Fragment_AS(a);

      mistakes = AnalyzeAffineAlign_MaxIndelSize(align->ahg,align->bhg,
				    a->sequence,b->sequence,trace,amode,
				    alen,blen,del,sub,ins,
				    affdel,affins,
				    blockdel,blockins,blocksize,biggestBlock);

      if (align->orientation == AS_INNIE)
        { int alength, blength, i, j, c;
          Complement_Fragment_AS(b);
          alength = strlen(a->sequence)+1;
          blength = strlen(b->sequence)+1;
          for (i = alength; mistakes[i] >= 0; i++)
            mistakes[i] = blength - mistakes[i];
          j = alength;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
      else if (align->orientation == AS_OUTTIE)
        { int alength, i, j, c;
          Complement_Fragment_AS(a);
          alength = strlen(a->sequence)+1;
          for (i = 0; mistakes[i] >= 0; i++)
            mistakes[i] = alength - mistakes[i];
          j = 0;
          i -= 1;
          while (j < i)
            { c = mistakes[i];
              mistakes[i--] = mistakes[j];
              mistakes[j++] = c;
            }
        }
    }

  if (swap)
    { int x;
      x = *alen;
      *alen = *blen;
      *blen = x;
      x = *del;
      *del = *ins;
      *ins = x;
      x = *affdel;
      *affdel = *affins;
      *affins = x;
      return (mistakes + strlen(a->sequence) + 1);
    }
  return (mistakes);
} 





/*** QUALITY VALUE REALIGNMENT ROUTINE ***/

/* Path abstraction that models cells that a trace
     alignment passes through in a given row.       */

typedef struct {
    int   row;    /* current row */
    int   low;    /* trace passes thru (row,low)  */
    int   hgh;    /*                to (row,hgh) of row */
    int   tpt;    /* next indel of trace to go through (or e.o.t.) */
} PathSlice;

static int *RAtrace;   /* Trace path slices are following */

static void Extend_Row(PathSlice *f)   /* Traverse indels on current row */
{ while (RAtrace[f->tpt] == (f->row+1))
    { f->hgh += 1;
      f->tpt += 1;
    }
}

static void Advance_Path(PathSlice *f)  /* Move to next row and traverse */
{ register int i, j, c;                 /*   indels on that row as well. */

  /* N.B.: when outside the boundary of the d.p. matrix the effect of
           Advance_Path is to move along a diagonal edge and stop.    */

  i = f->hgh + 1;
  if (RAtrace[c = f->tpt] == -i)
    { i -= 1;
      c += 1;
    }
  j = (f->row += 1) + 1;
  f->low = i;
  while (RAtrace[c] == j)
    { i += 1;
      c += 1;
    }
  f->hgh = i;
  f->tpt = c;
}

/* Realign aseq and bseq in a band of width BAND_WIDTH about alignment
   encoded in itrace.  On input *ahang gives the existing a-overhang
   (see Overlap doc) and *bhang gives the exisiting b-overhang.  Both
   are non-negative integers.  A pointer to the unpacked trace representing
   the best alignment in the band is returned, *ahang is set to specify
   the a-overhang for the new alignment, and *bhang is set to specify its
   b-overhang.  Both may be *negative* integers.                          */

static int *ReAlign(InternalFragMesg *a, InternalFragMesg *b,
                      int *ahang, int *bhang, int *itrace)

  /* The portion of the d.p. matrix computed is stored in the array dpspace
     which is large enough for all possible alignments.  Upon completion of
     the d.p. computation, the rows of the d.p. matrix in the band will be
     in the interval [dpbot,dptop].  The smallest column in row j that is in
     the band is at dpidx[j].  The width of the banded region in row j is
     dploc[j+1] - dploc[j].  The value of a particular cell (i,j) in the
     band is dploc[j][i-dpidx[j]].                                        */

{ static int dpidx[AS_READ_MAX_LEN+1], *dploc[AS_READ_MAX_LEN+2];
  static int dpspace[(2*AS_READ_MAX_LEN + 2)*(2*BAND_WIDTH+1)];
  int        dptop, dpbot;

  /* Buffer for the output trace */

  static int RealignBuffer[2*AS_READ_MAX_LEN+1];
  int       *otrace;

  char *aseq, *bseq;  /* A and B sequences, qv-arrays, and lengths */
  char *aqlv, *bqlv;
  int   alen,  blen;

  aseq    = a->sequence - 1;
  bseq    = b->sequence - 1;
  aqlv    = a->quality;
  bqlv    = b->quality;
  alen    = strlen(aseq+1);
  blen    = strlen(bseq+1);
  RAtrace = itrace;          /*  Set global trace used by path operators */

  //#ifndef UNDEFINE_UNIT_INDEL
  //  #define DEL(a,qa)      1
  //  #define SUB(a,qa,b,qb) ((a)!=(b))
  //#endif

  { PathSlice  mid, lag, adv;  /* Path slices for current row,
                                  current row-BAND, & current row + BAND */
    int        plow, phgh;     /* cells in band of previous row have column
                                  index in the interval [plow,phgh]  */
    int       *dpc, *dpl;      /* ptrs to first d.p. cells in current and
                                  previous rows */
  
    /* Setup initial path slices */

    mid.tpt = 0;      /* Setup current row slice so it is the leftmost cell */
    mid.row = 0;      /*  in row 0 (always starts in row 0 as *ahang >= 0.  */
    mid.low = mid.hgh = *ahang;
  
    lag = mid;             /* Setup lag so it start BAND_WIDTH rows earlier */
    lag.row -= BAND_WIDTH;
    lag.low -= BAND_WIDTH;
    lag.hgh -= BAND_WIDTH;
  
    Extend_Row(&mid);  /* Catch any other cells on the trace in the 1st row */
  
    { int j;           /* Setup adv so it starts BAND_WIDTH rows later */
  
      adv = mid;
      for (j = 1; j <= BAND_WIDTH; j++)
        Advance_Path(&adv);
    }
  
#ifdef DEBUG
    printf("\n");
    printf("  %3d: %3d-%-3d",lag.row,lag.low,lag.hgh);
    printf("  %3d: %3d-%-3d",mid.row,mid.low,mid.hgh);
    printf("  %3d: %3d-%-3d",adv.row,adv.low,adv.hgh);
#endif

    /* Compute the first row of the band */
  
    { register int low, hgh; /* cells in band of current row will have column
                                          indices in the interval [low,hgh] */
      low = mid.low-BAND_WIDTH;
      if (low > lag.low) low = lag.low;
      if (low < 0) low = 0;
  
      hgh = mid.hgh+BAND_WIDTH;
      if (hgh < adv.hgh) hgh = adv.hgh;
      if (hgh > alen) hgh = alen;
  
      dpl = dpc = dpspace;  /* setup d.p. row pointers */
      dpbot = mid.row;
  
      dpidx[mid.row] = low; /* establish band record for current row */
      dploc[mid.row] = dpc;
  
      { int i, e;
  
        *dpc++ = e = 0;
        if (mid.row == 0)                /* 1st row on A-boundary */
          for (i = low+1; i <= hgh; i++)
            *dpc++ = 0;
        else                             /* 1st row not on A-boundary */
          for (i = low+1; i <= hgh; i++)
            *dpc++ = e = e + DEL(aseq[i],aqlv[i]);
      }
    
      plow = low;
      phgh = hgh;
    }
  
#ifdef DEBUG
    printf("  %3d-%-3d",plow,phgh);
    printf("\n");
#endif

    /* Compute each subsequent row of the band til done. */
  
    while (1)
  
      { Advance_Path(&adv);   /* Advance path slices */
        Advance_Path(&mid);
        Advance_Path(&lag);
  
        /* Done if mid is thru B-boundary or BAND cells beyond A-boundary */

        if (mid.low - BAND_WIDTH > alen || bseq[mid.row] == 0) break;
  
#ifdef DEBUG
        printf("  %3d: %3d-%-3d",lag.row,lag.low,lag.hgh);
        printf("  %3d: %3d-%-3d",mid.row,mid.low,mid.hgh);
        printf("  %3d: %3d-%-3d",adv.row,adv.low,adv.hgh);
#endif
  
        { register int low, hgh; /* cells in band of current row will have
                                    column indices in the interval [low,hgh] */
          low = mid.low-BAND_WIDTH;
          if (low > lag.low) low = lag.low;
          if (low < 0) low = 0;
  
          hgh = mid.hgh+BAND_WIDTH;
          if (hgh < adv.hgh) hgh = adv.hgh;
          if (hgh > alen) hgh = alen;
  
          dpidx[mid.row] = low; /* establish band record for current row */
          dploc[mid.row] = dpc;
  
          /* Do the d.p.!  Invariant at the start of computing each cell (i,j)
             is: i = index of current column, e = dp(i-1,j), c = dp(i-1,j-1). */

          { int i, e, c, d;
  
            i = low;
            if (i == 0)          /* 1st cell on B-boundary */
              { *dpc++ = e = c = 0;
                dpl += 1;
              }
            else if (i == plow)  /* 1st cell directly above 1st of last row */
              *dpc++ = e = (c = *dpl++) + DEL(bseq[mid.row],bqlv[mid.row]);
            else                 /* 1st cell right of 1st cell of last row */
              { dpl += (i-plow)-1;
                e = (*dpl++) + SUB(aseq[i],aqlv[i],bseq[mid.row],bqlv[mid.row]);
                d = (c = *dpl++) + DEL(bseq[mid.row],bqlv[mid.row]);
                if (e > d) e = d;
                *dpc++ = e;
              }
            for (i++; i <= phgh; i++)  /* 2nd thru last column of last row */
              { e += DEL(aseq[i],aqlv[i]);
                c += SUB(aseq[i],aqlv[i],bseq[mid.row],bqlv[mid.row]);
                if (e > c) e = c;
                d = (c = *dpl++) + DEL(bseq[mid.row],bqlv[mid.row]);
                if (e > d) e = d;
                *dpc++ = e;
              } 
            if (i <= hgh)         /* Cells right of last column of last row */
              { e += DEL(aseq[i],aqlv[i]);                /* 1st one */
                c += SUB(aseq[i],aqlv[i],bseq[mid.row],bqlv[mid.row]);
                if (e > c) e = c;
                *dpc++ = e;
                for (i++; i <= hgh; i++)          /* The remainder */
                  *dpc++ = e = e + DEL(aseq[i],aqlv[i]);
              }
          }
  
          plow = low;
          phgh = hgh;
        }
  
#ifdef DEBUG
        printf("  %3d-%-3d",plow,phgh);
        printf("\n");
#endif
      }

    dptop = mid.row-1;      /* Finish up band description */
    dploc[mid.row] = dpc;
  }

#ifdef DEBUG
  { int i, *d;   /* Print out the d.p. band */

    printf("\n          ");
    for (i = dpidx[dpbot]+1; i <= alen; i++)
      printf("  %c",aseq[i]);
    printf("\n");

    for (i = dpbot; i <= dptop; i++)
      { printf("%3d: %c ",i,((i>0)?bseq[i]:' '));
        printf("%*s",3*(dpidx[i]-dpidx[dpbot]),"");
        for (d = dploc[i]; d < dploc[i+1]; d++)
          printf("%3d",*d);
        printf("\n");
      }
  }
#endif

  /* Traverse the boundary of the band to find the minimum cell (row,col).
     If there are several candidates, pick the one that is closest to the
     original finish.                                                    */

  { int row, col;

    { int min, trg, *d, b, i;

      d = dploc[dptop+1];
      b = dpidx[dptop];
      min = *--d;
      row = dptop;
      col = i = b + (d - dploc[dptop]);
#ifdef DEBUG
      printf("\nBoundary Minimum Computation:\n");
      printf("  (%d,%d) -> %d\n",row,col,min);
#endif
      if (dptop == blen)    /* Traverse part on B-boundary ... */
        { trg = alen + *bhang;  /* column of original finish */
          for (i--; i >= b; i--)
            { d -= 1;
              if (*d < min || (*d == min && abs(i-trg) < abs(col-trg)))
                { col = i;
                  min = *d;
                }
#ifdef DEBUG
              printf("  (%d,%d) -> %d\n",row,i,*d);
#endif
            }
        } 
      d = dploc[dptop];          /* ... and part on A-boundary */
      trg = blen - *bhang;              /* row of original finish */
      for (i = dptop-1; dpidx[i] + (d - dploc[i]) > alen; i--)
        { d -= 1;
          if (*d < min || (*d == min && abs(i-trg) < abs(row-trg)))
            { row = i;
              col = alen;
              min = *d;
            }
#ifdef DEBUG
          printf("  (%d,%d) -> %d\n",i,alen,*d);
#endif
          d = dploc[i];
        }
#ifdef DEBUG
      printf("-------------\n");
      printf("  (%d,%d) -> %d\n",row,col,min);
#endif
    }

    /* Compute new trace via backtrack through the d.p. matrix.  Accumulate
       unpacked trace in *reverse* order.                                   */

    if (col == alen)         /* Set B overhang */
      *bhang = blen-row;
    else
      *bhang = -(alen-col);

    { int *d, *e;   /* d = @dp(row,col), e == @dp(row-1,col)  */
      int *dc, *dl; /* dc = dploc[row], dl == dploc[row-1] */

      otrace = RealignBuffer + 2*AS_READ_MAX_LEN;  /* start at end of buffer */
      *otrace = 0;
      d = (dc = dploc[row]) + (col-dpidx[row]);
      e = (dl = dploc[row-1]) + (col-dpidx[row-1]);
      while (row > 0 && col > 0)
        { if (e > dl && e <= dc && e[-1] + SUB(aseq[col],aqvl[col],bseq[row],bqvl[row]) == *d)
            { row -= 1;
              col -= 1;    /* substitution */
              d    = e-1;
              dc   = dl;
              e    = (dl = dploc[row-1]) + (col-dpidx[row-1]);
            }
          else if (e >= dl && e < dc && *e + DEL(bseq[row],bqvl[row]) == *d)
            { row -= 1;
              d    = e;    /* delete in B */
              dc   = dl;
              e    = (dl = dploc[row-1]) + (col-dpidx[row-1]);
              *--otrace = -(col+1);
            }
          else if (d > dc && d[-1] + DEL(aseq[col],aqvl[col]) == *d)
            { col -= 1;
              d   -= 1;    /* delete in A*/
              e   -= 1;
              *--otrace = (row+1);
            }
        }
    }

    if (row == 0)         /* Set A overhang */
      *ahang = col;
    else
      *ahang = -row;
  }

  return (otrace);
}

/* Given fragments a and b and overlap align between them (computed by an
   identity-based method), refine the overlap using quality values for the
   fragments and return a pointer to a buffered overlap message containing
   the new result.                                                        */

OverlapMesg *QV_ReAligner_AS(InternalFragMesg *a, InternalFragMesg *b,
                             OverlapMesg *align)
{ int   *itrace, *otrace;
  int    ahang, bhang, orient;
  static OverlapMesg QVBuffer;

  itrace = Unpack_Alignment_AS(align);

  ahang  = align->ahg;
  bhang  = align->bhg;
  orient = align->orientation;

  if( a->iaccession == align->bifrag)
    {
      InternalFragMesg* t;
      t = a;
      a = b;
      b = t;
    }

  if (orient == AS_INNIE)        /* complement sequences as necessary */
    Complement_Fragment_AS(b);
  else if (orient == AS_OUTTIE)
    Complement_Fragment_AS(a);
 
  otrace = ReAlign(a,b,&ahang,&bhang,itrace);   /* compute new alignment */

  if (orient == AS_INNIE)       /* restore complemented sequence */
    Complement_Fragment_AS(b);
  else if (orient == AS_OUTTIE)
    Complement_Fragment_AS(a);

  //#define DEBUG
#ifdef DEBUG
  { int i;
    printf("\nOutput trace:  ahang = %d bhang = %d\n",ahang,bhang);
    for (i = 0; otrace[i] != 0; i++)
      printf("   %3d\n",otrace[i]);
  }
#endif

  /* Unfortunately, the new alignment may put A and B into a different
     orientation and overlap type with respect to each other and the new
     overlap record has to reflect that.  Indeed, if the A overhang has
     become negative then the orientation of A and B changes and the roles
     of A and B need to be reversed.                                       */

  QVBuffer = *align;
  if (ahang < 0)       /* Argh: Switch A and B and orientations */
    { 
      if (bhang >= 0)
        QVBuffer.overlap_type = AS_CONTAINMENT;
      else
        QVBuffer.overlap_type = AS_DOVETAIL;
      QVBuffer.ahg = -ahang;
      QVBuffer.bhg = -bhang;
      if (orient == AS_INNIE)
        QVBuffer.orientation = AS_OUTTIE;
      else if (orient == AS_OUTTIE)
        QVBuffer.orientation = AS_INNIE;

      QVBuffer.aifrag = b->iaccession;
      QVBuffer.bifrag = a->iaccession;

      { int i;
        for (i = 0; otrace[i] != 0; i++)
          otrace[i] = -otrace[i];
      }
    }
  else  /* The easy case: ahang is still positive so orientation is the same */
    { if (bhang <= 0)
        QVBuffer.overlap_type = AS_CONTAINMENT;
      else
        QVBuffer.overlap_type = AS_DOVETAIL;

      QVBuffer.ahg = ahang;
      QVBuffer.bhg = bhang;
    }

#ifdef DEBUG
  { int i;
    printf("\nOutput trace:  ahang = %d bhang = %d\n",ahang,bhang);
    for (i = 0; otrace[i] != 0; i++)
      printf("   %3d\n",otrace[i]);
  }
  printf("\nQV:  ahang = %d bhang = %d\n",QVBuffer.ahg,QVBuffer.bhg);
#endif

  /* Compress new trace to a delta */

  QVBuffer.delta = Pack_Alignment_AS(otrace,QVBuffer.ahg);

#ifdef DEBUG
  Print_Overlap_AS(stdout,a,b,&QVBuffer);
#endif

  return (&QVBuffer);
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
    O = LD_DP_Compare_AS(a,b,beg,end,opposite,max(MAXDPERATE,erate),thresh,minlen,AS_FIND_AFFINE_ALIGN,where,my_ld_ratio);
  }

  if (O != NULL){

    int del, sub, ins, affdel, affins, alen, blen, blockdel, blockins;
    float errRate, errRateAffine;

    Analyze_Affine_Overlap_AS(a,b,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
			      &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE);

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


#undef DEBUG_MEMALLOC

void *ckalloc(size_t size)	/* Guarded malloc utility */
{ void *newp;
  assert(size>0); 
  newp = (void *) malloc(size);
  if (newp == NULL)
    { fprintf(stderr,"Out of memory %s line %d\n",__FILE__,__LINE__);
      exit (1);
    }

#ifdef DEBUG_MEMALLOC
  fprintf(stderr,"\nAllocated %d to %X\n",size,newp);
#endif

  return (newp);
}

void *ckrealloc(void* ptr, size_t size)	/* Guarded realloc utility */
{ 
  void* newp;
  assert(ptr!=NULL);
  assert(size>0);
#ifdef DEBUG_MEMALLOC
  fprintf(stderr,"\nReallocating %X to size %d\n",ptr,size);
#endif
  newp = (void *) realloc(ptr,size);
  if (newp == NULL)
    { fprintf(stderr,"Out of memory %s line %d\n",__FILE__,__LINE__);
      exit (1);
    }

#ifdef DEBUG_MEMALLOC
  fprintf(stderr,"\n\tto %X\n",newp);
#endif

  return(newp);
}


void *ckreallocNullOK(void* ptr, size_t size)	/* Guarded realloc utility */
{ 
  void* newp;
  assert(size>0);
#ifdef DEBUG_MEMALLOC
  fprintf(stderr,"\nReallocating (null?) %X to size %d\n",ptr,size);
#endif
  newp = (void *) realloc(ptr,size);
  if (newp == NULL)
    { fprintf(stderr,"Out of memory %s line %d\n",__FILE__,__LINE__);
      exit (1);
    }
#ifdef DEBUG_MEMALLOC
  fprintf(stderr,"\n\tto %X\n",newp);
#endif

  return(newp);
}

