
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
/* MODULE FOR READING FASTA SEQUENCES, COMPUTING OVERLAPS, AND THE COLUMN
   SETS FOR THE CORRELATED-DIFFERENCES DETECTOR.
*/

#undef INPUT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_global.h"
#include "AS_ALN_aligners.h"

//void *safe_malloc(size_t size)	/* Guarded malloc utility */
//{ void *new;
//
//  new = (void *) malloc(size);
//  if (new == NULL)
//    { fprintf(stderr,"Out of memory\n");
//      exit (1);
//    }
//  return (new);
//}

/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */

#define LBUFLEN 512

char *get_sequence(FILE *input)
{ static char *seqbuf, linebuf[LBUFLEN];
  static int   first = 1;
  static int   top, nei;

  register char *newbuf;
  register size_t l;
  register int e, bol, beg;

  if (first)
    { first  = 0;
      top    = 2048;
      seqbuf = (char *) safe_malloc(sizeof(char)*top);
      if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
      if (*linebuf != '>')
        { fprintf(stderr,"First line must start with an >-sign\n");
          exit (1);
        }
    }
  else
    { if (!nei) return (NULL); }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  beg = 1;
  e   = 0;
  while((nei = (fgets(linebuf,LBUFLEN,input) != NULL)) != 0)
    { if (bol && *linebuf == '>')
        if (beg)
          { do
              { l = strlen(linebuf);
                if (linebuf[l-1] == '\n') break;
              }
            while (fgets(linebuf,LBUFLEN,input) != NULL);
          }
        else
          break;
      else
        { l = strlen(linebuf);
          if (e + l >= top)
            { top = (int) (1.5*(e+l) + 200);
              newbuf = (char *) safe_malloc(sizeof(char)*top);
              seqbuf[e] = '\0';
              strcpy(newbuf,seqbuf);
              safe_free(seqbuf);
              seqbuf = newbuf;
            }
          strcpy(seqbuf+e,linebuf);
          bol = (linebuf[l-1] == '\n');
          e = (e+l) - bol;
          beg = 0;
        }
    }
  seqbuf[e] = '\0';

  newbuf = (char *) safe_malloc(sizeof(char)*(e+1));
  strcpy(newbuf,seqbuf);
  
  return (newbuf);
}

/* Get_sequences gets all the FASTA formatted sequences from input, where
   it is assuming the input is a series of such sequences.  It sets *nseq
   to the number of sequences read less one, allocates space for each
   sequence, and returns an allocated array, seqa[0..k], of pointers to
   the sequences.
*/

char **get_sequences(FILE *input, int *nseq)
{ int    max, k;
  char **seqa, **seqn;

  max  = 32;
  seqa = (char **) safe_malloc(max*sizeof(char *));

  k = 0;
  while (1)
    { for (; k < max; k++)
        { seqa[k] = get_sequence(input);
          if (seqa[k] == NULL) break;
        }
      if (k < max) break;
      seqn = (char **) safe_malloc(2*max*sizeof(char *));
      for (k = 0; k < max; k++)
        seqn[k] = seqa[k];
      safe_free(seqa);
      seqa = seqn;
      max *= 2;
    }

  *nseq = k-1;
  return (seqa);
}

/* Write seq on the standard output, 50 symbols to a line. */

void show_sequence(char *seq)
{ size_t len, i;

  len = strlen(seq);
  for (i = 0; i < len; i += 50)
    if (i+50 < len)
      printf("%.50s\n",seq+i);
    else
      printf("%s\n",seq+i);
}

int main(int argc, char *argv[])
{ int    K;
  char **Seqs;
  InternalFragMesg  A, B;
  char qlty1[AS_READ_MAX_LEN+1], qlty2[AS_READ_MAX_LEN+1];

  Seqs = get_sequences(stdin,&K);

#ifdef INPUT
  { int i;

    printf("\nThe Sequences %d:\n\n",K+1);
    for (i = 0; i <= K; i++)
      { printf("> %d\n",i+1);
        show_sequence(Seqs[i]);
      }
  }
#endif

  { int i, j;

    A.quality = qlty1;
    B.quality = qlty2;
    for (j = 0; j < K; j++)
      for (i = j+1; i <= K; i++)
        { A.sequence = Seqs[j];
          B.sequence = Seqs[i];
          A.iaccession = A.eaccession = j+1;
          B.iaccession = B.eaccession = i+1;
          BPnt_Compare_AS(&A,&B,0,0,0,.02,1e-6,5,5);
/*
          BPnt_Compare_AS(&A,&B,-strlen(B.sequence),strlen(A.sequence),
                          0,.02,1e-6,40);
*/
        }
  }

  exit (0);
  return 0;
}
