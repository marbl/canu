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
#include <math.h>
#include <ctype.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_CGB_all.h"
#include "AS_MSG_pmesg.h"
#include "AS_ALN_aligners.h"

/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */


#define LBUFLEN 512

extern int max_indel_AS_ALN_LOCOLAP_GLOBAL;
extern int aimforlongest;
extern int kmerlen;
extern int minmatch;
extern int maxerror;
extern int kthresh;

void usage(void)
{
	fprintf( stderr, "usage: AS_ALN_boxfiller <-a abnd> <-b bbnd> <-e error_rate> <-m minlen> <-P proto_output_file> [-f.orwardOnly] [-n.oAlign] [-I]\n");
	fprintf( stderr, "\t-n suppresses printing of alignment\n");
	fprintf( stderr, "\t-f causes only forward orientations to be searched\n");
	fprintf( stderr, "\t-I uses near-identity settings, appropriate for finding 2% errors in 100-bp matches (much faster)\n");
	exit(1);
}

char *get_sequence(FILE *input)
{ static char *seqbuf, linebuf[LBUFLEN];
  static int   first = 1;
  static int   top, nei;

  register char *newbuf;
  register size_t l;
  register int e, bol, beg;

  if (first)
    { first  = 0;
      top    = 1024;
      seqbuf = (char *) ckalloc(sizeof(char)*top);
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
  while (nei = (fgets(linebuf,LBUFLEN,input) != NULL))
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
  seqa = (char **) ckalloc(max*sizeof(char *));

  k = 0;
  while (1)
    { for (; k < max; k++)
        { seqa[k] = get_sequence(input);
          if (seqa[k] == NULL) break;
        }
      if (k < max) break;
      seqn = (char **) ckalloc(2*max*sizeof(char *));
      for (k = 0; k < max; k++)
        seqn[k] = seqa[k];
      free(seqa);
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
  OverlapMesg  *O;
  char qlty1[AS_READ_MAX_LEN+1], qlty2[AS_READ_MAX_LEN+1];
  FILE *OVLFile=NULL;
  int abnd, abndFromUser = FALSE;
  int bbnd, bbndFromUser = FALSE;
  int minlen=40;
  int first=1;
  int doRevToo=1;
  int noAlign=0;

  double err=.06;

  aimforlongest=1;

#ifdef OLDWAY	  
  if(argc>2&&strcmp(argv[1],"-P")==0){
    fprintf(stderr,"Printing OVLs to %s\n",argv[2]);
    OVLFile=fopen(argv[2],"w");
    assert(OVLFile!=NULL);
  }
  if(argc>2)
  {
	  if(strcmp(argv[1],"-P")!=0|| argc>=6) {
		  assert(argc>=4);
		  err=atof(argv[argc-3]);
		  abnd=atoi(argv[argc-1]);
		  bbnd=atoi(argv[argc-2]);
		  if((strcmp(argv[1],"-P")!=0&&argc>4)||argc>6){
			  minlen=atoi(argv[argc-4]);
		  }
	  }
  }
#else
	  { /* Parse the argument list using "man 3 getopt". */ 
		  int ch,errflg=0;
		  optarg = NULL;
		  while (!errflg && ((ch = getopt(argc, argv, "a:b:e:hm:P:fnI")) != EOF))
		  {
#if 0
			  fprintf(GlobalData->stderrc,"* ch = %c optopt= %c optarg = %s\n", ch, optopt, (optarg?optarg:"(NULL)"));
			  fflush(stderr);
#endif
			  switch(ch) {
				  case 'a':
					  abnd = atoi(optarg);
					  abndFromUser = TRUE;
					  break;
				  case 'b':
					  bbnd = atoi(optarg);
					  bbndFromUser = TRUE;
					  break;
				  case 'e':
					  err = atof(optarg);
					  break;
           			  case 'f':
				    doRevToo=0;
				    break;
				  case 'h':
					  usage();
					  break;
			          case 'I':  //tunings for near-id searches
				          kmerlen=12;
					  minmatch=100;
					  maxerror=2;
					  kthresh=65;
				          break;
				  case 'm':
					  minlen = atoi(optarg);
					  break;
			          case 'n':
				    noAlign=1;
				    break;
				  case 'P':
					  // outputPath = strdup(optarg);
					  fprintf(stderr,"Printing OVLs to %s\n", optarg);
					  OVLFile=fopen( optarg, "w");
					  assert(OVLFile!=NULL);
					  break;
				  case '?':
					  fprintf( stderr, "Unrecognized option -%c\n", optopt);
					  usage();
				  default :
					  errflg++;
					  usage();
			  }
		  }
	  }
#endif

  Seqs = get_sequences(stdin,&K);
  fprintf(stderr,"Read in %d sequences\n",K+1);

#ifdef INPUT
  { int i;

    printf("\nThe Sequences %d:\n\n",K+1);
    for (i = 0; i <= K; i++)
      { printf("> %d\n",i+1);
        show_sequence(Seqs[i]);
      }
  }
#endif

  { int i, j, where;
    int olaps, tlaps;
    int ori;

    A.quality = NULL;
    B.quality = NULL;
    tlaps = olaps = 0;


    for (j = 0; j < K; j++)
      for (i = j+1; i <= K; i++)
        { 
	  //	  printf("\n\nComparing sequences %d and %d\n\n",j,i);
	  tlaps += 1;
          A.sequence = Seqs[j];
          B.sequence = Seqs[i];
          A.iaccession = A.eaccession = j+1;
          B.iaccession = B.eaccession = i+1;
	  //	  printf("\n\nForward comparison results:\n\n");

		  if ( !abndFromUser )
			  abnd = strlen(A.sequence);
		  if ( !bbndFromUser )
			  bbnd = -strlen(B.sequence);

	  if(first){
		  first=0;
		  fprintf(stderr,"Calling BoxFill_AS with band %d to %d,erate %f, minlen %d\n",
				  bbnd,abnd,err,minlen);
	  }
	  

	  for(ori=0;ori<1+doRevToo;ori++){

#ifdef DEBUG_BUBBLE_SMOOTHING
            fprintf(stderr,"BoxFill_AS: err=%f ori=%d minlen=%d\n",
                    err, ori, minlen);
#endif            
	    O = BoxFill_AS(&A,&B,bbnd,abnd,
			       ori,
			       err,
			       1e-6,
			       minlen,
			       AS_FIND_LOCAL_ALIGN,
			       &where);
#ifdef DEBUG_BUBBLE_SMOOTHING
            if( NULL != O) {
              fprintf(stderr,"O->aifrag=%ld O->bifrag=%ld\n", O->aifrag, O->bifrag);
              fprintf(stderr,"O->ahg=%d     O->bhg=%d\n", O->ahg, O->bhg);
              fprintf(stderr,"O->orientation=%c\n", O->orientation);
              fprintf(stderr,"O->overlap_type=%c\n", O->overlap_type);
              fprintf(stderr,"O->quality=%f\n", O->quality);
            }
#endif
            
	    if (O != NULL){
	      olaps += 1;
              if(!noAlign)Print_Overlap_AS(stdout,&A,&B,O);
	      { 
		int alen,blen,del,sub,ins,affdel,affins,blockdel,blockins;
		double errRate,errRateAffine;
		int AFFINEBLOCKSIZE=4;
		Analyze_Affine_Overlap_AS(&A,&B,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
					  &affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, NULL);
		
		errRate = (sub+ins+del)/(double)(alen+ins);
		
		errRateAffine = (sub+affins+affdel)/
		  (double)(alen-del+affins+affdel);
		
		printf("\n\nAlen %d, Blen %d, del %d, sub %d, ins %d\n"
		       " affdel %d, affins %d, blockdel %d, blockins %d\n",
		       alen,blen,del,sub,ins,
		       affdel,affins,blockdel,blockins);
		printf("Simple mismatch rate %f\n",errRate);
		printf("Affine mismatch rate %f\n",errRateAffine);

		printf("Largest block mismatch %d\n",max_indel_AS_ALN_LOCOLAP_GLOBAL);

		O->min_offset=O->max_offset=O->ahg;
		if(OVLFile!=NULL){
		  GenericMesg pmesg;
		  pmesg.m=O;
		  pmesg.t=MESG_OVL;
		  WriteProtoMesg_AS(OVLFile,&pmesg);
		}
	      }

            }
	  }
        }
    fprintf(stderr,"Performed %d x 2 compares, found %d overlaps\n",
            tlaps,olaps);
  }

  return(0);
}
