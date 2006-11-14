
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
/* MODULE FOR READING FASTA SEQUENCES, AND COMPUTING OVERLAPS USING LOCAL
   ALIGNMENT METHODS
*/

#undef INPUT
#define DEBUG 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <AS_ALN_aligners.h>

/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */

#define LBUFLEN 512


void *ckalloc(size_t size)	/* Guarded malloc utility */
{ void *newp;
  assert(size>0); 
  newp = (void *) malloc(size);
  if (newp == NULL)
    { fprintf(stderr,"Out of memory\n");
      exit (1);
    }
  return (newp);
}

void *ckrealloc(void* ptr, size_t size)	/* Guarded realloc utility */
{ 
  void* newp;
  assert(ptr!=NULL);
  assert(size>0);
  newp = (void *) realloc(ptr,size);
  if (newp == NULL)
    { fprintf(stderr,"Out of memory\n");
      exit (1);
    }
  return(newp);
}


char *get_sequence(FILE *input, char **seq, char **name )
{ static char *seqbuf, *namebuf,nextname[LBUFLEN],linebuf[LBUFLEN];
  static int   first = 1;
  static int   top, nei;

  register char *newbuf,*newbuf2;
  register size_t l;
  register int e, bol, beg;

  if (first)
    { first  = 0;
      top    = 2048;
      seqbuf = (char *) ckalloc(sizeof(char)*top);
      namebuf = (char *) ckalloc(sizeof(char)*top);
      if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
      if (*linebuf != '>')
        { fprintf(stderr,"First line must start with an >-sign\n");
          exit (1);
        }
      else
	{
	  char *newname;
	  newname = (char*) ckalloc(2048*sizeof(char));
	  if(sscanf(linebuf,">%s",newname)!=1){
	    if(sscanf(linebuf,"> %s",newname)!=1){
	      fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	    }
	  } 
	  if(strlen(newname)>2047){
	    fprintf(stderr,"identifier %s too long -- abort!\n",
		    newname);
	  }
	  newname = (char *)realloc(newname,strlen(newname)+1);
	  assert(newname!=NULL);
	  *name = newname;
	}   

    }
  else
    { 
      if (!nei) return (NULL); 
      if(*nextname == '>'){
	char *newname;
	newname = (char*) ckalloc(2048*sizeof(char));
	if(sscanf(nextname,">%s",newname)!=1){
	  if(sscanf(nextname,"> %s",newname)!=1){
	    fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	  }
	} 
	if(strlen(newname)>2047){
	  fprintf(stderr,"identifier %s too long -- abort!\n",
		  newname);
	}
	newname = (char *)realloc(newname,strlen(newname)+1);
	assert(newname!=NULL);
	*name = newname;
      }   
    }

  do
    { l = strlen(linebuf);
      if (linebuf[l-1] == '\n') break;
    }
  while (fgets(linebuf,LBUFLEN,input) != NULL);

  bol = 1;
  beg = 1;
  e   = 0;
  while((nei = (fgets(linebuf,LBUFLEN,input) != NULL)) != 0)
    {
      if (bol && *linebuf == '>')
      {
        if (beg)
          { 
	    do
              { l = strlen(linebuf);
                if (linebuf[l-1] == '\n') break;
              }
            while (fgets(linebuf,LBUFLEN,input) != NULL);
          }
        else{
	  strcpy(nextname,linebuf);
          break;
	}
      }else
        { l = strlen(linebuf);
          if (e + l >= top)
            { top = (int) (1.5*(e+l) + 200);
              newbuf = (char *) ckalloc(sizeof(char)*top);
              seqbuf[e] = '\0';
              strcpy(newbuf,seqbuf);
              free(seqbuf);
              seqbuf = newbuf;
            }
          strcpy(seqbuf+e,linebuf);
          bol = (linebuf[l-1] == '\n');
          e = (e+l) - bol;
          beg = 0;
        }
    }
  seqbuf[e] = '\0';

  newbuf = (char *) ckalloc(sizeof(char)*(e+1));
  strcpy(newbuf,seqbuf);

  {
    int i;
    for(i=0;newbuf[i]!='\0';i++){
      newbuf[i]=toupper(newbuf[i]);
    }
  }
  
  *seq = newbuf;
  return newbuf;
}

/* Get_sequences gets all the FASTA formatted sequences from input, where
   it is assuming the input is a series of such sequences.  It sets *nseq
   to the number of sequences read less one, allocates space for each
   sequence, and returns an allocated array, seqa[0..k], of pointers to
   the sequences.
*/

void get_sequences(FILE *input, int *nseq,char ***seqs,char ***names)
{ int    max, k;
  char **seqa, **seqn;
  char **namea, **namen;

  max  = 32;
  seqa = (char **) ckalloc(max*sizeof(char *));
  namea = (char **) ckalloc(max*sizeof(char *));

  k = 0;
  while (1)
    { for (; k < max; k++)
        { if (get_sequence(input,&(seqa[k]),&(namea[k]))== NULL) break;
        }
      if (k < max) break;
      seqn = (char **) ckalloc(2*max*sizeof(char *));
      namen = (char **) ckalloc(2*max*sizeof(char *));
      for (k = 0; k < max; k++){
        seqn[k] = seqa[k];
        namen[k] = namea[k];
      }
      free(seqa);
      free(namea);
      seqa = seqn;
      namea = namen;
      max *= 2;
    }

  *nseq = k-1;
  *seqs=seqa;
  *names=namea;
}

static void Complement(char *seq, int len)
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

  /* Complement and reverse sequence */

  { register char *s, *t;
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

int main(int argc, char *argv[])
{ int    K;
  char **Seqs;
  char **Names;
  int i, j;
  int len;

  get_sequences(stdin,&K,&Seqs,&Names);
  for (j = 0; j <= K; j++){
    len=strlen(Seqs[j]);
    Complement(Seqs[j],len);
    printf(">%s revcompl\n",Names[j]);
    for(i=0;i<len;i+=60){
      //      printf("%.60s\n",Seqs[j]+i);
      int left;
      left=MIN(60,len-i);
      fwrite(Seqs[j]+i,sizeof(char),left,stdout);
      putc('\n',stdout);
    }
  }
  return(0);
}
