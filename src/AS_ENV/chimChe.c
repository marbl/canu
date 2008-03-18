#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "kmer_utils.h"
#include "AS_ALN_aligners.h"
#include "AS_UTL_alloc.h"

#undef DEBUG_SEGMENTATION

static int  Map[128];
int ksize=10;
int kmax;
int fullnames=0;
int dbfullnames=0;

int GRANULARITY=10;

/* Get_sequence gets the next FASTA formated sequence from input.  It
   assumes that it and only it is doing input and is designed to handle
   lines and strings of arbitrary length.  It allocates memory for the
   sequence and returns a pointer to it.  */

#define LBUFLEN 512
 
void usage(char *pgm)
{
	fprintf( stderr, 
		 "usage: %s -q <query.fasta> [-d <db.fasta>] -k <kmer size> [-f] [-g <skip>] [-m minlen] [-t]\n"
		 "\t-q specifies multi-fasta sequences to be analyzed\n"
		 "\t-d specifies multi-fasta sequences to be indexed (if missing, -q file used)\n"
		 "\t-k specifies the size of the substrings; <= 14, default 10\n"
		 "\t-g skip gives the number of bp to skip between segmentation tests\n"
		 "\t-m minlen specifies smallest overlap to be used\n"
		 "\t-f causes full deflines to be printed\n"
		 "\t-t causes best simple, best left and best right to be trimmed to the same window, to reduce partial-sequence artifacts\n",
		 pgm
		 );
	exit(1);
}


static int firstget;

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
	  if(fullnames){
	    int l;
	    if(linebuf[0] != '>'){
		fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
		exit(1);
	    }
	    l=strlen(linebuf);
	    if(l>2048){
	      fprintf(stderr,"identifier %s too long -- abort!\n",
		      linebuf);
	      exit(1);
	    }
	    strcpy(newname,linebuf+1);
	    newname[l-2]='\0';
	  } else {
	    if(sscanf(linebuf,">%s",newname)!=1){
	      if(sscanf(linebuf,"> %s",newname)!=1){
		fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
		exit(1);
	      }
	    } 
	  }
	  if(strlen(newname)>2047){
	    fprintf(stderr,"identifier %s too long -- abort!\n",
		    newname);
	    exit(1);
	  }
	  newname = (char *)safe_realloc(newname,strlen(newname)+1);
	  *name = newname;
	}   

    }
  else
    { 
      if (!nei) return (NULL); 
      if(*nextname == '>'){
	char *newname;
	newname = (char*) safe_malloc(2048*sizeof(char));
	if(fullnames){
	  int l; 
	  if(linebuf[0] != '>'){
	    fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	    exit(1);
	  }
	  l=strlen(linebuf);
	  if(l>2048){
	    fprintf(stderr,"identifier %s too long -- abort!\n",
		    linebuf);
	    exit(1);
	  }
	  strcpy(newname,linebuf+1);
	  newname[l-2]='\0';
	} else {
	  if(sscanf(linebuf,">%s",newname)!=1){
	    if(sscanf(linebuf,"> %s",newname)!=1){
	      fprintf(stderr,"Abort: Couldn't resolve defline %s\n",linebuf);
	      exit(1);
	    }
	  } 
	}
	if(strlen(newname)>2047){
	  fprintf(stderr,"identifier %s too long -- abort!\n",
		  newname);
	}
	newname = (char *)safe_realloc(newname,strlen(newname)+1);
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
{ int    maxl, k;
  char **seqa, **seqn;
  char **namea, **namen;

  maxl  = 32;
  seqa = (char **) safe_malloc(maxl*sizeof(char *));
  namea = (char **) safe_malloc(maxl*sizeof(char *));
  firstget=1;
  k = 0;
  while (1)
    { for (; k < maxl; k++)
        { if (get_sequence(input,&(seqa[k]),&(namea[k]))== NULL) break;
        }
      if (k < maxl) break;
      seqn = (char **) safe_malloc(2*maxl*sizeof(char *));
      namen = (char **) safe_malloc(2*maxl*sizeof(char *));
      for (k = 0; k < maxl; k++){
        seqn[k] = seqa[k];
        namen[k] = namea[k];
      }
      safe_free(seqa);
      safe_free(namea);
      seqa = seqn;
      namea = namen;
      maxl *= 2;
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


void calc_kmer_members(const char **Seqs,const int ksize,int **kmerCounts,int **kmerIndex,const int numseqs){
  int i,j;
  int kword;
  int *Counts, *Index;

  assert(ksize < 16);

  kmax=1;
  kmax = kmax << (ksize*2);
  kmax--;

  Counts = (int *)safe_malloc((kmax+1)*sizeof(int));

  for(i=0;i<=kmax;i++){
    Counts[i]=0;
  }

  // count all kmers
  for(i=0;i<numseqs;i++){
    int h = -ksize;
    kword=0;
    for(j=0;j<ksize-1;j++){
      int x = Map[(int) (Seqs[i][j])];
      if (x >= 0){
        kword = (kword << 2) | x;
      }else{
	kword <<= 2; 
	h = j-(ksize-1);
      } 
    }
    while(Seqs[i][j]!='\0'){
      int x = Map[(int) (Seqs[i][j])];
      if (x >= 0){
        kword = ((kword << 2) | x) & kmax;
      }else{
	kword <<= 2;
	h = j-(ksize-1);
      } 
      if (j >= h+ksize){
        Counts[kword+1] ++;
      }
      j++;
    }
  }

  // cumulative offset;
  for (kword = 2; kword <= kmax+1; kword++){
    Counts[kword] += Counts[kword-1];
  }

  {// fill in the index
    int numkmers = Counts[kmax+1];

    Index = (int *) safe_malloc(numkmers * sizeof(int));

#if 0
    // this should be unnecessary since we fill the whole table!
    for (i = 0; i <= numkmers; i++){
      Index[i] = -1;
    }
#endif

    for(i=0;i<numseqs;i++){
      int h = -ksize;
      j=0;
      kword=0;
      for(j=0;j<ksize-1;j++){
	int x = Map[(int) (Seqs[i][j])];
	if (x >= 0){
	  kword = (kword << 2) | x;
	}else{
	  kword <<= 2; 
	  h = j-(ksize-1);
	} 
      }
      while(Seqs[i][j]!='\0'){
	int x = Map[(int) (Seqs[i][j])];
	if (x >= 0){
	  kword = ((kword << 2) | x) & kmax;
	}else{
	  kword <<= 2; 
	  h = j-(ksize-1);
	} 
	if (j >= h+ksize){
	  Index[Counts[kword]++] = i;
	}
	j++;
      }
    }

    // since we advanced to end, pointers into index are actually
    // for the next kmer; shift 'em.
    for (kword = kmax; kword >= 0; kword--){
      Counts[kword+1] = Counts[kword];
    }
    Counts[0]=0;

  }

  *kmerCounts = Counts;
  *kmerIndex = Index;
  return;

}

typedef struct hittag {
  int count;
  int index;
} hits;

int hit_compare(const void *a, const void *b){
  hits *A=(hits*)a;
  hits *B=(hits*)b;
  return(B->count-A->count);
}




int main(int argc, char *argv[])
{ int    K=-1,KB=-1;
  char **Seqs;
  char **Names;
  char **SeqsB;
  char **NamesB;
  char *seqfilename=NULL,*dbfilename=NULL;
  int internalCompare=0; /* whether query and database sequences are the same */
  int *Profiles;
  int *ProfilesB=NULL;
  int ori;
  int first=1;
  FILE *seqfile=NULL, *dbfile=NULL;
  int *kmerCounts,*kmerIndex;
  int *len,*lenB;
  int *frontDiscount;
  int i,j;
  int maxlen=0;
  int minlen=40;
  int doTrimming=0;

  argc = AS_configure(argc, argv);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv, "fFg:hk:q:d:m:t")) != EOF))
      {
	switch(ch) {
	case 'f':
	  fullnames=1;
	  break;
	case 'F':
	  dbfullnames=1;
	  break;
	case 'g':
	  GRANULARITY = atoi(optarg);
	  assert(GRANULARITY>0&&GRANULARITY<500);
	  break;
	case 'h':
	  usage(argv[0]);
	  break;
	case 'q':
	  if(seqfile!=NULL)fclose(seqfile);
	  seqfile=fopen(optarg,"r");
	  assert(seqfile!=NULL);
	  seqfilename=optarg;
	  if(dbfile==NULL){
	    dbfile=fopen(optarg,"r");
	    assert(dbfile!=NULL);
	    dbfilename=optarg;
	  }
	  break;
	case 'd':
	  if(dbfile!=NULL){
	    fclose(dbfile);
	  }
	  dbfile=fopen(optarg,"r");
	  assert(dbfile!=NULL);
	  dbfilename=optarg;
	  break;
	case 'k':
	  ksize=atoi(optarg);
	  assert(ksize>0&&ksize<14);
	  break;
	case 'm': 
	  minlen=atoi(optarg);
	  assert(minlen>=40&&minlen<1400);
	  break;
	case 't':
	  doTrimming=1;
	  fprintf(stderr,"Will trim to common range\n");
	  break;
	case '?':
	  errflg++;
	  usage(argv[0]);
	  break;
	default :
	  fprintf( stderr, "Unrecognized option -%c\n", optopt);
	  errflg++;
	  usage(argv[0]);
	}
      }
  }

  assert(seqfile!=NULL);
  assert(dbfile!=NULL);
  if(seqfilename==dbfilename || strcmp(seqfilename,dbfilename)==0){
    internalCompare=1;
  } else {
    internalCompare=0;
  }

  get_sequences(seqfile,&K,&Seqs,&Names);
  assert(K>=0);
  fprintf(stderr,"Read in %d query sequences\n",K+1);
  len = (int *)safe_malloc(sizeof(int)*(K+1));
  for(i=0;i<=K;i++){
    len[i]=strlen(Seqs[i]);
    if(maxlen<len[i]){
      maxlen=len[i];
    }
  }

  if(dbfullnames)fullnames=1;

  get_sequences(dbfile,&KB,&SeqsB,&NamesB);
  assert(KB>=0);
  fprintf(stderr,"Read in %d database sequences\n",KB+1);
  lenB = (int *)safe_malloc(sizeof(int)*(KB+1));
  frontDiscount = (int *)safe_malloc(sizeof(int)*(KB+1));
  for(i=0;i<=KB;i++){
    lenB[i]=strlen(SeqsB[i]);
    frontDiscount[i]=0;
  }

  for (i = 0; i < 128; i++){
    Map[i] = -1;
  }
  Map['a'] = Map['A'] = 0;
  Map['c'] = Map['C'] = 1;
  Map['g'] = Map['G'] = 2;
  Map['t'] = Map['T'] = 3;

  calc_kmer_members((const char **)SeqsB,ksize,&kmerCounts,&kmerIndex,KB+1);

  fprintf(stderr,"Built index\n");


  { 
    int **hitCounts=NULL;

    hitCounts = (int **) safe_malloc(maxlen*sizeof(int*));
    
    for(j= 0; j<maxlen; j++){
      hitCounts[j] = (int *) safe_malloc((KB+1)*sizeof(int));
    }

    for (i = 0; i <= K; i++){
      int k;
      int kword = 0;
      int h = -ksize;
      int bestfront=-1, bestback=-1;
      int bestscore=-1,bestloc=-1;
      int bestsimple=-1,simplescore=-1;
      int ilen,ilenlessone;
      int startbestmatch=0;
      int endbestmatch = len[i]-1;

      ilen=len[i];
      ilenlessone=ilen-1;

      for(k=0;k<ksize;k++){
	for(j= 0; j<=KB; j++){
	  hitCounts[k][j] = 0;
	}
      }

      for(j=0;j<ksize-1;j++){
	int x = Map[(int) (Seqs[i][j])];
	if (x >= 0){
	  kword = (kword << 2) | x;
	}else{
	  kword <<= 2; 
	  h = j-(ksize-1);
	} 
      }


      while(Seqs[i][j]!='\0'){
	int x = Map[(int) (Seqs[i][j])];
	if (x >= 0){
	  kword = ((kword << 2) | x) & kmax;
	}else{
	  kword <<= 2;
	  h = j-(ksize-1);
	} 
	for(k=0;k<=KB;k++){
	  hitCounts[j][k]=hitCounts[j-1][k];
	}
	if (j >= h+ksize){
	  for(k=kmerCounts[kword];k<kmerCounts[kword+1];k++){
	    if(internalCompare && kmerIndex[k]==i)continue;
	    if(k==kmerCounts[kword]||(kmerIndex[k]!=kmerIndex[k-1])){
	      hitCounts[j][kmerIndex[k]]++;
	    }
	  }
	}
	j++;
      }

      for( k=0;k<=KB;k++){
	if(hitCounts[ilenlessone][k]>simplescore){
	  simplescore = hitCounts[ilenlessone][k];
	  bestsimple=k;
	}
      }  
      
      {
	Overlap *ovl=NULL;
	double erate=0.02;
	// below, .9 fudge factor may be necessary to handle cases where some matching kmers are random out of order matches
	int minovl=simplescore *.9;
	while(erate<.4){
	  ovl = DP_Compare(Seqs[i],SeqsB[bestsimple],-lenB[bestsimple]+minovl,len[i]-minovl,0,erate,1e-6,maxim(minlen,minovl),AS_FIND_LOCAL_ALIGN_NO_TRACE);
	  if(ovl!=NULL)break;
	  erate*=2.;
	}
	if(ovl!=NULL){
	  if(ovl->begpos>0){
	    startbestmatch=ovl->begpos;
	  }
	  if(ovl->endpos<0){
	    endbestmatch=len[i]+ovl->endpos;
	    assert(endbestmatch>0);
	  }
	}
      }

      //      printf("startbestmatch init at %d\n",startbestmatch);
      //      printf("endbestmatch init at %d\n",endbestmatch);

      for(k=0;k<=KB;k++){
	frontDiscount[k]=hitCounts[startbestmatch+ksize-2][k];
      }


      for(j=startbestmatch+ksize-1;j<=endbestmatch;j+=GRANULARITY){
	int maxfront=-1;
	int maxback=-1;
	int whichfront=-1,whichback=-1; 
	for(k=0;k<=KB;k++){
	  if(hitCounts[j][k]-frontDiscount[k]>maxfront){
	    maxfront=hitCounts[j][k]-frontDiscount[k];
	    whichfront=k;
	  }
	  if(hitCounts[endbestmatch][k]-hitCounts[j][k]>maxback){
	    maxback=hitCounts[endbestmatch][k]-hitCounts[j][k];
	    whichback=k;
	  }
	}
	if(maxfront+maxback>bestscore){
	  bestscore=maxfront+maxback;
	  bestfront=whichfront;
	  bestback=whichback;
	  bestloc=j;
#ifdef DEBUG_SEGMENTATION
	  fprintf(stderr,"New best %d: loc %d seqs %s / %s bestfront %d bestback %d\n",
		  bestscore,
		  bestloc,NamesB[bestfront],NamesB[bestback],bestfront,bestback);
#endif
	}
      }

      if(doTrimming){
	int frontstart=startbestmatch;
	int frontend=endbestmatch;

	int backstart=startbestmatch;
	int backend=endbestmatch;

	Overlap *ovl=NULL;
	double erate=0.02;
	// below, .9 fudge factor may be necessary to handle cases where some matching kmers are random out of order matches
	int minovl=hitCounts[bestloc][bestfront]*.9;
	if(bestfront!=bestsimple){
	  while(erate<.4){
	    ovl = DP_Compare(Seqs[i],SeqsB[bestfront],-lenB[bestfront]+minovl,len[i]-minovl,0,erate,1e-6,maxim(minlen,minovl),AS_FIND_LOCAL_ALIGN_NO_TRACE);
	    if(ovl!=NULL)break;
	    erate*=2.;
	  }
	  if(ovl!=NULL){
	    if( minim(frontend,len[i]+ovl->endpos) <= maxim(frontstart,ovl->begpos) ){
	      // complain
	      fprintf(stderr,
		      "trouble with overlap found at erate %f: ahang %d bhang %d => frontend %d <= frontstart %d\n"
		      ">Qseq\n%s\n>Dseq\n%s\n",
		      erate,
		      ovl->begpos, ovl->endpos, frontend, frontstart,Seqs[i],SeqsB[bestfront]);
	      // and do nothing!
	    } else {
	      // update 
	      frontstart = maxim(frontstart,ovl->begpos);
	      frontend = minim(frontend,len[i]+ovl->endpos);
	    }
	  }
	}
	erate=0.02;
	// below, .9 fudge factor may be necessary to handle cases where some matching kmers are random out of order matches
	minovl=(hitCounts[ilenlessone][bestback]-hitCounts[bestloc][bestback])*.9; 
	//	fprintf(stderr,"DEBUG: minovl = %d\n",hitCounts[ilenlessone][bestback]);
	ovl=NULL;
	if(bestback!=bestsimple){
	  //	  fprintf(stderr,"initial settings: backstart, backend to %d %d\n",backstart,backend);
	  while(erate<.4){
	    ovl = DP_Compare(Seqs[i],SeqsB[bestback],-lenB[bestback],len[i],0,erate,1e-6,maxim(minlen,minovl),AS_FIND_LOCAL_ALIGN_NO_TRACE);
	    //ovl = DP_Compare(Seqs[i],SeqsB[bestback],-lenB[bestback],len[i],0,erate,1e-6,40,AS_FIND_LOCAL_ALIGN_NO_TRACE);
	    if(ovl!=NULL)break;
	    erate*=2.;
	  }
	  if(ovl!=NULL){

	    if( minim(backend,len[i]+ovl->endpos) <= maxim(backstart,ovl->begpos) ){
	      // complain
	      fprintf(stderr,
		      "trouble with overlap found at erate %f: ahang %d bhang %d => backend %d <= backstart %d\n"
		      ">Qseq\n%s\n>Dseq\n%s\n",
		      erate,
		      ovl->begpos, ovl->endpos, backend, backstart,Seqs[i],SeqsB[bestback]);
	      // and do nothing
	    } else {
	      // update 
	      backstart=maxim(backstart,ovl->begpos);
	      backend=minim(backend,len[i]+ovl->endpos);
	      //	      fprintf(stderr,"Updating backstart, backend to %d %d\n",backstart,backend);
	    }
	  }
	}
	
	// things are problematic if the overlap is to a region that doesn't come close to the implied breakpoint (bestloc);
	// however, partial sequence issues can make an absolute test fail, hence the constant 100 below:
	if(frontstart<bestloc+100&&frontend>bestloc-100){
	  startbestmatch= (startbestmatch > frontstart ? startbestmatch : frontstart);
	  endbestmatch= (endbestmatch < frontend ? endbestmatch : frontend);
	}
	if(backstart<bestloc+100&&backend>bestloc-100){
	  startbestmatch= (startbestmatch > backstart ? startbestmatch : backstart);
	  endbestmatch= (endbestmatch < backend ? endbestmatch : backend);
	}

	//	fprintf(stderr,"Final bestmatch start %d end %d\n",startbestmatch,endbestmatch);
	
	simplescore=hitCounts[endbestmatch][bestsimple]-hitCounts[startbestmatch+ksize-2][bestsimple];
	bestscore=
	  hitCounts[endbestmatch][bestback]-hitCounts[bestloc][bestback]+
	  hitCounts[bestloc][bestfront]-hitCounts[startbestmatch+ksize-2][bestfront];


      }

      //      printf("startbestmatch final at %d\n",startbestmatch);
      //      printf("endbestmatch final at %d\n",endbestmatch);

      assert(startbestmatch<=bestloc+100&&endbestmatch>bestloc-100);


      if(simplescore <= hitCounts[endbestmatch][bestfront]-hitCounts[startbestmatch+ksize-2][bestfront]){
	bestsimple=bestfront;
	simplescore=hitCounts[endbestmatch][bestfront]-hitCounts[startbestmatch+ksize-2][bestfront];
	//	fprintf(stderr,"Resetting best simple, presumably due to trimming (now bestfront)\n");
      }

      if( simplescore <= hitCounts[endbestmatch][bestback]-hitCounts[startbestmatch+ksize-2][bestback]){
	bestsimple=bestback;
	simplescore=hitCounts[endbestmatch][bestback]-hitCounts[startbestmatch+ksize-2][bestback];
	//	fprintf(stderr,"Resetting best simple, presumably due to trimming (now bestback)\n");
      }

      if(bestscore>simplescore){
	printf("%s (len=%d): best split = %d %s : %s (score = %d from %d to %d ; separately, scores = %d and %d; best single %s scores %d )\n",
	       Names[i],ilen,bestloc,NamesB[bestfront],NamesB[bestback],
	       bestscore,startbestmatch,endbestmatch,
	       hitCounts[endbestmatch][bestfront]-hitCounts[startbestmatch+ksize-2][bestfront],
	       hitCounts[endbestmatch][bestback]-hitCounts[startbestmatch+ksize-2][bestback],
	       NamesB[bestsimple],simplescore);
      } else {
	printf("%s (len=%d): best match = %s (from %d to %d, score = %d )\n",
	       Names[i],ilen,NamesB[bestsimple],startbestmatch,endbestmatch,simplescore);
      }

    }

  }

  return (0);
}
