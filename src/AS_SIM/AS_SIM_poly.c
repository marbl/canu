
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

static const char *rcsid = "$Id: AS_SIM_poly.c,v 1.4 2008-10-08 22:03:00 brianwalenz Exp $";

/*
   Module:      Polymorphism simulator, a UNIX filter.
   Description: See CDS/Assember/poly.doc
   Author:      Gene Myers
   Written:     October 10, 1998
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#undef  DEBUG           /* Causes extra diagnostic stuff to get printed */
#define INPUTMAX  2048  /* Maximum poly spec. line length (checked) */

static char *ProgName;  /* Program name 'poly' */

static int   Seed;      /* Seed for random number generation */

static int   Comments;  /* Don't print out any comments */

static FILE *SpecFile;  /* File to read poly spec. from */

static FILE *SeqFile;   /* File containing string to be poly'd */

                        /* Spec. record type scalar: */
#define SNP         0   /*   SNP type */
#define DELETE      1   /*   Deletion type */
#define TRANSLOCATE 2   /*   Translocation type */

typedef struct { int    type;   /* Type of spec record */
                 int    min;    /* Segment size range is [min,max] for   */
                 int    max;    /*   deletion and translocation records. */
                 double rate;   /* % of genome to be operated on. */
                 int    number; /* Number of ops (calculated) */
               } mutator;

static mutator *Operator; /* Array of poly spec. records */
static int      OpCnt;    /* Number of poly spec. records #'d 0,1,...OpCnt-1 */
                          /*   Operator[0] is always a SNP record and the
                               only SNP record, all others are Del. or Tran. */

static int *Sizes;
static int *Locs;
static int *Sort;
static int *SNPs;


/* >>>> BASIC UTILITIES <<<< */

#define SEGLEN    50    /* Output sequences, SEGLEN chars to a line */

static char SubChar[8][3] = { {'c', 'g', 't'},  /* What chars to substitute */
                              {'a', 'g', 't'},
                              {'a', 'c', 't'},
                              {'a', 'c', 'g'},
                              {'C', 'G', 'T'},
                              {'A', 'G', 'T'},
                              {'A', 'C', 'T'},
                              {'A', 'C', 'G'} };

static int Decode[128];     /* Converts chars a,c,g,t to integers 0,1,2,3;
                               and A,C,G,T to integers 4,5,6,7.             */

  /* (Re)Allocate space and exit if none */

static void *ckalloc(int size)
{ void *m;
  if ((m = malloc(size)) == NULL)
    { fprintf(stdout,"\n\t*** Out of Memory (%s requesting %d bytes)\n",
              ProgName,size);
      exit (1);
    }
  return (m);
}

static void *ckrealloc(void *m,int size)
{ if ((m = realloc(m,size)) == NULL)
    { fprintf(stdout,"\n\t*** Out of Memory (%s requesting %d bytes)\n",
              ProgName,size);
      exit (1);
    }
  return (m);
}

  /* Randomly mutate an input character. */

static char substitute(int orig)
{ register int i;
  i = 3.0*drand48();
  return (SubChar[Decode[orig]][i]);
}


/* >>>> POLYMORPHISM ENGINE <<<< */

/* 'qsort' comparator for sorting Indel array */

static int DEL_ORDER(const void *x, const void *y)
{ int a, b;
  a = *((int *) x);
  b = *((int *) y);
  return (Locs[a] - Locs[b]);
}

/* 'qsort' comparator for sorting SNP array */

static int SNP_ORDER(const void *x, const void *y)
{ int a, b;
  a = *((int *) x);
  b = *((int *) y);
  return (a - b);
}

/* Output next character of poly'd sequence */

void output(int c)
{ static int row = 0;

  if (row++ >= SEGLEN)
    { putchar('\n');
      row = 1;
    }
  putchar(c);
}

void generate_poly(char *dna, int len)
{ int      N, S, A, B;
  mutator *o;

  /* For each del. and trans. poly spec. record, determine number of ops to
     perform for that mutation and set in 'number' field of record.  On
     completion and for the remainder of the routine N = # of deletions
     and insertions to induce (a translocation is modeled as an indel pair.) */

{ register int i, n;
  register double a, x;

  N = 0;
  for (i = 1; i < OpCnt; i++)
    { o = Operator + i;
      a = o->rate * len;
      x = (2.*a) / (o->min + o->max);
      n = x;
      if (drand48() < x - n)
        n += 1;
      if (o->type == TRANSLOCATE) n *= 2;
      o->number = n;
      N += n;
    }
}

  /* Allocate (Sizes,Locs,Sort) records for indels, and SNPs array of SNPS */

  Sizes = (int *) ckalloc(sizeof(int)*N);
  Locs  = (int *) ckalloc(sizeof(int)*(N+1));
  Sort  = (int *) ckalloc(sizeof(int)*(N+1));

  /* Compute a list of indel operations as per the deletion and translocation
     spec. records.  At the finish A = size+1 of string after all deletions
     and B = size of string after all ops (deletion and insertion) */

{ register int i, j, k;

  j = 0;
  A = 0;
  B = 0;
  for (i = 1; i < OpCnt; i++)
    { o = Operator + i;
      for (k = 0; k < o->number; k++)
        { Sizes[j] = drand48()*(o->max - o->min + 1) + o->min;
          A += Sizes[j];
          if (o->type == TRANSLOCATE)
            { Sizes[++j] = -1; k += 1; }
          else
            B += Sizes[j];
          j += 1;
        }
    }
  A = len - A + 1;
  B = len - B;
  if (A <= 0)
    { fprintf(stdout,"\n\t*** You've asked to delete the entire sequence!\n");
      exit (1);
    }
}

  /* Compute the number, S, of SNPs to perform (had to wait as the number is
     a percentage of the *final* string length.  Allocate array for
     SNP locations. */

{ register int n;
  register double a;

  o = Operator;
  a = o->rate * B;
  n = a;
  if (drand48() < a - n)
    n += 1;
  S = o->number = n;

  SNPs  = (int *) ckalloc(sizeof(int)*(S+1));
}

  /* Compute a list of the positions (in the deleted image of the string) for
     each indel op. and set the Sort permuation to the identity.  Also set
     Locs[N] to a guardian value. */

{ register int j, b;

  for (j = 0; j <= N; j++)
    { Locs[j] = drand48() * A;  /* Uniform over [0,A-1]  */
      Sort[j] = j;              /*                  ^^ ! */
    }
  Locs[N] = len+1;

  /* Compute a list of the positions for each SNP op. */

  b = B - (S-1);
  for (j = 0; j < S; j++)
    SNPs[j] = drand48() * b;   /* Uniform over [0,b-1] */
  SNPs[S] = len+1;
}

  /* Determine the Sort permutation for the indel records, and sort the
     SNP positions directly. */

  qsort(Sort,N+1,sizeof(int),DEL_ORDER);
  qsort(SNPs,S+1,sizeof(int),SNP_ORDER);

#ifdef DEBUG
{ register int i, j, k;

  fprintf(stdout,"Reduced range = [0,%d]\n\n",A);
  for (j = 0; j < N; j++)
    { k = Sort[j];
      fprintf(stdout,"  %3d: %6d",j,Locs[k]);
      if (Sizes[k] > 0)
        fprintf(stdout," %4d-",Sizes[k]);
      else
        fprintf(stdout," %4d+",Sizes[k-1]);
      for (i = 1; i < OpCnt; i++)
        { o = Operator + i;
          if (k < o->number)
            break;
          k -= o->number;
        }
      fprintf(stdout," Member %d\n",i);
    }
}
#endif

  /* Now offset indel locations to positions in originating string,
     and offset SNPs by 1 per SNP (selection without replacement). */

{ register int j, k, offset;

  offset = 0;
  for (j = 0; j < N; j++)
    { k = Sort[j];
      Locs[k] += offset;
      if (Sizes[k] > 0)
        offset += Sizes[k];
    }

  for (j = 0; j < S; j++)
    SNPs[j] += j;
}

  /*  If requested print header giving all polymorphisms enacted on the
      string.  Note carefully that the SNP information is described in
      the coordinate of the string that results from first executing all
      indel operations. */

  if (Comments)
    { register int j, k, h;

      printf("# Did Structural Polymorphisms:\n");
      for (j = 0; j < N; j++)
        { k = Sort[j];
          if (Sizes[k] > 0)
            printf("#   Delete [%d,%d]\n",Locs[k],Locs[k]+Sizes[k]);
          else
            printf("#   Insert [%d,%d] at %d\n",
                   Locs[k-1],Locs[k-1]+Sizes[k-1],Locs[k]);
        }
      printf("#\n");

      printf("# Followed by SNPs:\n");
      for (h = 0; h < S; h++)
        printf("#   SNP at %d\n",SNPs[h]+1);
      printf("#\n");
    }

#ifdef DEBUG
{ register int i, j, k;

  fprintf(stdout,"Reduced range = [0,%d]\n\n",A);
  for (j = 0; j < N; j++)
    { k = Sort[j];
      fprintf(stdout,"  %3d: %6d",j,Locs[k]);
      if (Sizes[k] > 0)
        fprintf(stdout," %4d-",Sizes[k]);
      else
        fprintf(stdout," %4d+",Sizes[k-1]);
      for (i = 1; i < OpCnt; i++)
        { o = Operator + i;
          if (k < o->number)
            break;
          k -= o->number;
        }
      fprintf(stdout," Member %d\n",i);
    }
}
#endif

  /* Finally, enact the polymorphism in a single pass over the originating
     string outputting the result directly to stdout.  Basically merge sort
     traversal of string (i = 0,1,...len-1) with positions of indel ops
     (Locs[Sort[0]],Locs[Sort[1]],...) performing these ops as go and then
     further merging with SNPs list thereafter. */

{ register int i, j, k, h;
  register int a, b, n;

  printf(">\n");
  k = Sort[j = 0];
  b = h = 0;
  i = 0;
  while (i < len || k < N)
    if (i < Locs[k])
      { if (h == SNPs[b])
          { output(substitute(dna[i++])); b += 1; }
        else
          output(dna[i++]);
        h += 1;
      }
    else if (Sizes[k] > 0)
      { i += Sizes[k];
        k = Sort[++j];
      }
    else
      { a = Locs[k-1];
        n = a + Sizes[k-1];
        while (a < n)
          { if (h == SNPs[b])
              { output(substitute(dna[a++])); b += 1; }
            else
              output(dna[a++]);
            h += 1;
          }
        k = Sort[++j];
      }
  putchar('\n');
}

}


/* >>>> DNA SEQUENCE INPUTTER <<<< */

/* Read_seq gets the first FASTA formated sequence from the given input file.
   It is designed to handle lines and strings of arbitrary length.  It
   allocates memory for the sequence and returns a pointer to it.  */

#define LBUFLEN 512  /* Buffer length for sequence input (not a limit!) */

char *read_seq(FILE *input, int *len)
{ static   char  linebuf[LBUFLEN];
  register char *seqbuf;
  register int l, e, bol = 0, top;

  top = 2048;
  seqbuf = (char *) ckalloc(sizeof(char)*top);
  if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
  if (*linebuf == '#')
  {
    do
    {
      if (bol && *linebuf == '>')
        break;
      else if (bol)
      {
        if (*linebuf != '#')
        {
          fprintf(stderr,"Poly read_seq: No FASTA header after a comment\n\t%s", linebuf);
          exit (1);
        }
        else
          /* Do your # line processing here */;
      }
      l = strlen(linebuf);
      bol = (linebuf[l-1] == '\n');
    }
    while (fgets(linebuf,LBUFLEN,input) != NULL);
  }
  if (*linebuf != '>')
    { fprintf(stdout,"First line must start with an >-sign\n");
      exit (1);
    }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  e   = 0;
  while (fgets(linebuf,LBUFLEN,input) != NULL)
    { if (bol && *linebuf == '>')
        break;
      l = strlen(linebuf);
      if (e + l >= top)
        { top = 1.5*(e+l) + 200;
          seqbuf = (char *) ckrealloc((void *) seqbuf,sizeof(char)*top);
        }
      strcpy(seqbuf+e,linebuf);
      bol = (linebuf[l-1] == '\n');
      e = (e+l) - bol;
    }
  seqbuf[e] = '\0';

  if (e == 0)
    { fprintf(stdout,"Input sequence is empty\n");
      exit (1);
    }

  *len = e;
  return (seqbuf);
}


/* >>>> POLYMORPHISM SPECIFICATION INPUTTER <<<< */

/* Read and "parse" the specification file building the array
   Operator[0..OpCnt-1] of specification records for which record
   0 is always *the* SNP record. */

void getinput(void)
{ static char buffer[INPUTMAX];

  /* Get the name of file to read DNA from and open it. */

  if (fgets(buffer,INPUTMAX,SpecFile) == NULL)
    { fprintf(stdout,"\n\t*** Expecting input\n");
      exit (1);
    }
  if (buffer[strlen(buffer)-1] != '\n')
    { fprintf(stdout,"\n\t*** Input line is longer than %d chars\n",INPUTMAX);
      exit (1);
    }
  buffer[strlen(buffer)-1] = '\0';
  SeqFile = fopen(buffer,"r");
  if (SeqFile == NULL)
    { fprintf(stdout,"\n\t*** Cannot open file '%s'\n",buffer);
      exit (1);
    }
  if (Comments)
    { printf("#\n# Sequence from file: %s\n#\n",buffer);
      fflush(stdout);
    }

  /*  Read in each spec line and parse into an Operator record */

{ int otop;
  int minlen, maxlen;
  double mrate, srate;

  otop = 128;
  Operator = (mutator *) ckalloc(sizeof(mutator)*otop);
  OpCnt = 1;

  srate = 0.;
  while (fgets(buffer,INPUTMAX,SpecFile) != NULL)
    switch (*buffer)

    { case 'S':  /* SNP spec. */
        if (sscanf(buffer+1," %lf",&mrate) != 1)
          { fprintf(stdout,"Expect '%%lf' after S\n");
            exit (1);
          }
        srate += mrate;
        break;

      case 'D':  /* Delete spec. */
      case 'X':  /* Translocate sepc. */
        if (OpCnt >= otop)
          { otop *= 2;     /* Enlarge array if necessary */
            Operator = (mutator *)
                          ckrealloc((void *) Operator,sizeof(mutator)*otop);
          }
        if (sscanf(buffer+1," %d - %d %lf",&minlen,&maxlen,&mrate) != 3)
          { fprintf(stdout,"Expect '%%d - %%d %%lf' after D or X\n");
            exit (1);
          }
        Operator[OpCnt].min  = minlen;
        Operator[OpCnt].max  = maxlen;
        Operator[OpCnt].rate = mrate;
        if (*buffer == 'D')
          { if (Comments)
              printf("# Delete %g%% of the genome in blocks of %d-%d bp\n",
                     100.*mrate,minlen,maxlen);
            Operator[OpCnt].type = DELETE;
          }
        else
          { if (Comments)
              printf("# Translocate %g%% of the genome in blocks of %d-%d bp\n",
                     100.*mrate,minlen,maxlen);
            Operator[OpCnt].type = TRANSLOCATE;
          }
        OpCnt += 1;
        break;

      default:
        fprintf(stdout,"\n\t*** Unrecognized command %c\n",*buffer);
        exit (1);
    }

  Operator[0].type = SNP;  /* Entry 0 is always the SNP entry */
  Operator[0].min  = 1;
  Operator[0].max  = 1;
  Operator[0].rate = srate;
  if (Comments)
    fprintf(stdout,"# SNP rate = %g%%\n#\n",100.*srate);
}

#ifdef DEBUG
{ int i;
  for (i = 0; i < OpCnt; i++)
    { fprintf(stdout,"Op %2d: ",i);
      if (Operator[i].type == SNP)
        fprintf(stdout,"SNP ");
      else if (Operator[i].type == DELETE)
        fprintf(stdout,"Del ");
      else
        fprintf(stdout,"Xfr ");
      fprintf(stdout," %g [%d,%d]\n",Operator[i].rate,
                     Operator[i].min,Operator[i].max);
    }
}
#endif

}


/* >>>> MAIN / TOP <<<< */

int main(int argc, char *argv[])
{
  ProgName = argv[0];

/* Process program flags and setup random number generator seed as requested. */

{ register int i;
  int illegal, no_seed;

  no_seed  = 1;
  Comments = 1;
  illegal  = 0;
  while (argc > 1 && argv[1][0] == '-')
    { if (argv[1][1] == 's')
        { if (argc > 2)
            { Seed    = atoi(argv[2]);
              no_seed = 0;
              argc   -= 1;
              argv   += 1;
            }
        }
      else
        { for (i = 1; argv[1][i] != '\0'; i++)
            if (argv[1][i] == 'F')
              Comments = 0;
            else
              illegal = 1;
          if (i == 1) illegal = 1;
        }
      argv += 1;
      argc -= 1;
    }

  if (argc != 2 || illegal)
    { fprintf(stdout,"\n\t*** Usage is: %s [-s {seed:%%d}] [-F] specfile\n",
              ProgName);
      exit (1);
    }

  if (no_seed)
    Seed = getpid();
  srand48(Seed);

  if (Comments)
    fprintf(stdout,"#\n# Seed = %d\n",Seed);
}

/* Open spec. file */

  if ((SpecFile = fopen(argv[1],"r")) == NULL)
    { fprintf(stdout,"\n\t*** Cannot open spec file %s\n",argv[1]);
      exit (1);
    }

/* Setup mutation tables */

  Decode['a'] = 0;
  Decode['c'] = 1;
  Decode['g'] = 2;
  Decode['t'] = 3;
  Decode['A'] = 4;
  Decode['C'] = 5;
  Decode['G'] = 6;
  Decode['T'] = 7;

/* Read in the specification from SpecFile. */

  getinput();

/* Read and then permute DNA sequence.  */

{ char *dna;
  int   len;

  dna = read_seq(SeqFile,&len);

  generate_poly(dna,len);
}

  exit(0);
}
