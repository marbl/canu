
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
#include <ctype.h> // for toupper

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "AS_ALN_aligners.h"
#include "AS_CGB_all.h"


/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */

#define LBUFLEN 512

void usage(void)
{
	fprintf( stderr, "usage: AS_ALN_testdriver1 <-a abnd> <-b bbnd> <-e error_rate> <-m minlen> <-o|-O> <-P proto_output_file> [-1 <Aseq.fasta> -2 <Bseq.fasta>]\n");
	exit(1);
}

static int firstget =1;

char *get_sequence(FILE *input, char **seq, char **name )
{ static char *seqbuf, *namebuf,nextname[LBUFLEN],linebuf[LBUFLEN];
  static int   top, nei;

  register char *newbuf,*newbuf2;
  register size_t l;
  register int e, bol, beg;

  if (firstget)
    { firstget  = 0;
      top    = 2048;
      seqbuf = (char *) safe_malloc(sizeof(char)*top);
      namebuf = (char *) safe_malloc(sizeof(char)*top);
      if (fgets(linebuf,LBUFLEN,input) == NULL) return (NULL);
      if (*linebuf != '>')
        { fprintf(stderr,"First line must start with an >-sign\n");
          exit (1);
        }
      else
	{
	  char *newname;
	  newname = (char*) safe_malloc(2048*sizeof(char));
	  if(sscanf(linebuf,">%s",newname)!=1){
	    if(sscanf(linebuf,"> %s",newname)!=1){
	      fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	    }
	  }
	  if(strlen(newname)>2047){
	    fprintf(stderr,"identifier %s too long -- abort!\n",
		    newname);
	  }
	  newname = (char *)safe_realloc(newname,strlen(newname)+1);
	  assert(newname!=NULL);
	  *name = newname;
	}

    }
  else
    {
      if (!nei) return (NULL);
      if(*nextname == '>'){
	char *newname;
	newname = (char*) safe_malloc(2048*sizeof(char));
	if(sscanf(nextname,">%s",newname)!=1){
	  if(sscanf(nextname,"> %s",newname)!=1){
	    fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	  }
	}
	if(strlen(newname)>2047){
	  fprintf(stderr,"identifier %s too long -- abort!\n",
		  newname);
	}
	newname = (char *)safe_realloc(newname,strlen(newname)+1);
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

  {
    int i;
    for(i=0;newbuf[i]!='\0';i++){
      newbuf[i]=toupper(newbuf[i]);
      switch(newbuf[i]){
      case 'A':
      case 'C':
      case 'G':
      case 'T':
	break;
      default:
	newbuf[i]='N';
      }
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
  seqa = (char **) safe_malloc(max*sizeof(char *));
  namea = (char **) safe_malloc(max*sizeof(char *));

  k = 0;
  while (1)
    { for (; k < max; k++)
        { if (get_sequence(input,&(seqa[k]),&(namea[k]))== NULL) break;
        }
      if (k < max) break;
      seqn = (char **) safe_malloc(2*max*sizeof(char *));
      namen = (char **) safe_malloc(2*max*sizeof(char *));
      for (k = 0; k < max; k++){
        seqn[k] = seqa[k];
        namen[k] = namea[k];
      }
      safe_free(seqa);
      safe_free(namea);
      seqa = seqn;
      namea = namen;
      max *= 2;
    }

  *nseq = k-1;
  *seqs=seqa;
  *names=namea;
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
{ int    K,KB;
  char **Seqs;
  char **Names;
  char **SeqsB;
  char **NamesB;
  InternalFragMesg  A, B;
  OverlapMesg  *O;
  FILE *OVLFile=NULL;
  int ori;
  int abnd=0;
  int bbnd=0;
  int abndFromUser=0;
  int bbndFromUser=0;
  int minlen=40;
  int first=1;
  int printOlaps=0;
  int printOlapsOnly=0;
  int printAligns=1;
  double err=.06;
  FILE *file1=stdin,*file2=NULL;


  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "a:b:e:hm:oOP:1:2:")) != EOF))
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
	case 'h':
	  usage();
	  break;
	case 'm':
	  minlen = atoi(optarg);
	  break;
	case 'O':
	  printOlapsOnly=1;
	  printOlaps=1;
	  break;
	case 'o':
	  printOlaps=1;
	  break;
	case 'P':
	  // outputPath = strdup(optarg);
	  fprintf(stderr,"Printing OVLs to %s\n", optarg);
	  OVLFile=fopen( optarg, "w");
	  assert(OVLFile!=NULL);
	  break;
	case '1':
	  file1=fopen(optarg,"r");
	  assert(file1!=NULL);
	  break;
	case '2':
	  file2=fopen(optarg,"r");
	  assert(file2!=NULL);
	  break;
	case '?':
	  errflg++;
	  usage();
	  break;
	default :
	  fprintf( stderr, "Unrecognized option -%c\n", optopt);
	  errflg++;
	  usage();
	}
      }
  }


  get_sequences(file1,&K,&Seqs,&Names);
  fprintf(stderr,"Read in %d sequences\n",K+1);
  if(file2!=NULL){
    firstget=1;
    get_sequences(file2,&KB,&SeqsB,&NamesB);
    fprintf(stderr,"Read in %d sequences\n",KB+1);
  }


#ifdef INPUT
  { int i;

    printf("\nThe Sequences %d:\n\n",K+1);
    for (i = 0; i <= K; i++)
      { printf("> %d\n",i+1);
        show_sequence(Seqs[i]);
      }


    if(file2!=NULL){
      printf("\nThe OTHER Sequences %d:\n\n",K+1);
      for (i = 0; i <= KB; i++)
	{ printf("> %d\n",i+1);
        show_sequence(Seqs[i]);
	}
    }

  }
#endif


  {
    int i, j, where=1;
    int olaps, tlaps;

    A.quality = NULL;
    B.quality = NULL;
    tlaps = olaps = 0;
    for (j = 0; j <= K; j++){
      for (i = (file2==NULL ? j+1 : 0); i <= (file2==NULL ? K : KB); i++){
	tlaps += 1;
	A.sequence = Seqs[j];
	A.iaccession = j+1;
        A.eaccession = AS_UID_fromInteger(A.iaccession);

	B.sequence = (file2==NULL ? Seqs[i] : SeqsB[i]);
	B.iaccession = i+1 + (file2!=NULL ? (K+1) : 0);
        B.eaccession = AS_UID_fromInteger(B.iaccession);

	for(ori=0;ori<=1;ori++){

	  if(strlen(B.sequence)<-bbnd)
	  bbnd=-strlen(B.sequence);

	  O = DP_Compare_AS(&A,&B,
			    bbndFromUser?bbnd : -strlen(B.sequence),
			    abndFromUser?abnd : strlen(A.sequence),
			    ori,err,1e-6,minlen,AS_FIND_ALIGN,&where);
          if (O != NULL){
            olaps += 1;
	    if(!printOlapsOnly)Print_Overlap_AS(stdout,&A,&B,O);
	    if(!printOlapsOnly)printf("Overlap quality: %f\n",O->quality);

	    { int del, sub, ins, affdel, affins, alen, blen, blockdel, blockins;
	      float errRate, errRateAffine;

#define AFFINEBLOCKSIZE 4
	      Analyze_Affine_Overlap_AS(&A,&B,O,AS_ANALYZE_ALL,&alen,&blen,&del,&sub,&ins,
					&affdel,&affins,&blockdel,&blockins,AFFINEBLOCKSIZE, NULL);

	      errRate = (sub+ins+del)/(double)(alen+ins);

	      errRateAffine = (sub+affins+affdel)/
		(double)(alen+ins-(del-affdel+ins-affins));

	      if(!printOlapsOnly)printf("Alen %d, Blen %d, del %d, sub %d, ins %d\n"
		     " affdel %d, affins %d, blockdel %d, blockins %d\n",
		     alen,blen,del,sub,ins,
		     affdel,affins,blockdel,blockins);
	      if(!printOlapsOnly)printf("Simple mismatch rate %f\n",errRate);
	      if(!printOlapsOnly)printf("Affine mismatch rate %f\n",errRateAffine);
	      if(!printOlapsOnly)printf("dp_olap: %s %s %e\n",
		       Names[j],(file2==NULL ? Names[i] : NamesB[i]),errRate);


	      if(printOlaps){
		char ori;
		int ahang,bhang;

		Compute_Olap_Version(&A,&B,O,&ahang,&bhang,&ori);
		printf("OLAP: %s %s %c %d %d %f %f Len= %d\n",
		       Names[j],(file2==NULL ? Names[i] : NamesB[i]),ori,ahang,bhang,errRate,errRateAffine, (alen < blen ) ? alen : blen);
	      }

	    }
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
  return (0);
}
