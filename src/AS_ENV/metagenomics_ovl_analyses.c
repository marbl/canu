#include <stdlib.h>
#include <unistd.h>
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"
#include  "OlapStoreOVL.h"

//static int cutoffs[10]={300,250,200,150,100,60,50,20,10,5};
//static int ncuts=10;
#define ncuts 4
static int cutoffs[ncuts]={250,100,20,0};


static int prev;
static int t;
static int **cnt,**cln;
static int maxsmp;
static Long_Olap_Data_t  olap;
static double ***self,***other,***clnself,***clnother;
static int *smp;
static int debug=0;
static int psmp;
static int **allbeg,**allend,**otherbeg,**otherend,**selfbeg,**selfend;
static int *alln,*othern,*selfn;
static int maxovls=0;
static int flen;
static FragStoreHandle  my_frg_store;
static GateKeeperStore my_gkp_store ; // See AS_PER_gkpStore.h
static OVL_Store_t  * my_ovl_store;
static OVL_Stream_t  * my_stream;
static ReadStructp fsread;
static int currmate;
static int ***samesmpantiovlhistogram;
static int ***alldepthhistogram;
static int ***otherdepthhistogram;
static int ***selfdepthhistogram;
static int maxdepth=500;
static int *smpnfrgs;
static int *smpbases;


void usage(char *pgm){
  fprintf(stderr,"Usage: %s [-A <diversity overlap histogram file>] [-B <binary coverages file>] [-C <depths of coverage file>] [-S <similarity scores file>] [-b startfrg] [-e endfrg] [-d] -f <frgStore> -g <gkpStore> -o <ovlStore> -s <samples>\n",pgm);
  exit(-1);
}

void setup_stores(char *OVL_Store_Path, char *Frg_Store_Path, char *Gkp_Store_Path){

  assert(OVL_Store_Path!=NULL);
  assert(Frg_Store_Path!=NULL);
  assert(Gkp_Store_Path!=NULL);

  assert(my_ovl_store==NULL);
  my_ovl_store = New_OVL_Store ();
  my_stream = New_OVL_Stream ();
  Open_OVL_Store (my_ovl_store, OVL_Store_Path);

  InitGateKeeperStore(&my_gkp_store, Gkp_Store_Path);
  OpenReadOnlyGateKeeperStore(&my_gkp_store);

  my_frg_store = openFragStore( Frg_Store_Path, "r");
}

void init_frg(){

  int i,k;
  prev=olap.a_iid;
  psmp=smp[olap.a_iid];
  t=0;
  for(i=1;i<=maxsmp;i++){
    for(k=0;k<ncuts;k++){
      cnt[i][k]=0;
      cln[i][k]=0;
    }
  }

  for(k=0;k<ncuts;k++){
    othern[k]=selfn[k]=alln[k]=0;
  }


  {
    /*****************/
    // get fragment length
    /*****************/
    uint clr_bgn,clr_end;
    getFragStore(my_frg_store,prev,FRAG_S_ALL,fsread);
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    flen = clr_end-clr_bgn;
  }

  currmate=-1;
  {
    /*************************/
    // check for an appropriate mate
    /*************************/
    GateKeeperFragmentRecord gkpFrag;
    int rv1 = getGateKeeperFragmentStore(my_gkp_store.frgStore,prev,&gkpFrag);
    assert(rv1==0);
    if(gkpFrag.numLinks==1){
      GateKeeperLinkRecordIterator iterator;
      GateKeeperLinkRecord link;
      CreateGateKeeperLinkRecordIterator(my_gkp_store.lnkStore, gkpFrag.linkHead,prev, &iterator);
      while(NextGateKeeperLinkRecordIterator(&iterator, &link))
	currmate = (link.frag1 == prev) ? link.frag2 : link.frag1;
    }
  }

  if(debug>0){
    fprintf(stderr,"%d mate %d len %d smp %d\n",prev,currmate,flen,psmp);
  }
}


int intcmp(const void *a, const void *b){
  int *A=(int*)a;
  int *B=(int*)b;
  return (*A-*B);
}


void depth_of_coverage( int *begs, int *ends, int n, int *depthhistogram , char type, int k){

  int p=0;
  int d=0;
  int i=0;
  int j=0;

  if(debug>0){
    fprintf(stderr,"depth frg %d cutoff %d class %c",prev,cutoffs[k],type);
  }

  if(n==0){ 
    depthhistogram[0]+=flen;
    if(debug>0){
      fprintf(stderr,"\t0-%d:0\n",flen);
    }
    return;
  }

  qsort(begs,n,sizeof(int),intcmp);
  qsort(ends,n,sizeof(int),intcmp);

  assert(begs[0]<ends[0]);
  assert(begs[0]>=0);

  while(p<flen){
    //    ## main cases:
    //    # next is one or more starts
    //    # next is one or more stops
    //    # next is a combination of starts and stops
    //    # no more 

    if(i<n&&(j==n||begs[i]<ends[j])){

      // print "\t$p-$begs[$i]:$d" if ($begs[$i]>0);
      if(begs[i]>0){
	if(debug>0){
	  fprintf(stderr, "\t%d-%d:%d",p,begs[i],d);
	}
	depthhistogram[( d > maxdepth ? 500 : d)]+= begs[i]-p;
      }

      p=begs[i++];
      d++;

      while(i<n&&begs[i]==p){
	i++;
	d++;
      }

    } else if (j<n &&(i==n||begs[i]>ends[j])){  

      //      print "\t$p-$ends[$j]:$d";
      if(debug>0){
	fprintf(stderr,"\t%d-%d:%d",p,ends[j],d);
      }

      depthhistogram[( d > maxdepth ? 500 : d)] += ends[j]-p;

      p=ends[j++];
      d--;

      while(j<n&&ends[j]==p){
	j++;
	d--;
      }

    } else if (i<n && j<n && begs[i]==ends[j]){
      int prevp=p;
      int prevd=d;

      p=begs[i];

      while(i<n && begs[i]==p){
	i++;
	d++;
      }

      while(j<n && ends[j]==p){
	j++;
	d--;
      }

      if(prevd==d){
	p=prevp;
	d=prevd;
      } else {

	//	print "\t$prevp-$p:$prevd";
	if(debug>0){
	  fprintf(stderr,"\t%d-%d:%d",prevp,p,prevd);
	}

	depthhistogram[( prevd > maxdepth ? 500 : prevd)]+=p-prevp;

      }

    } else {

      assert(i>=n && j>=n);

      //      print "\t$p-$l:0";
      if(debug>0){
	fprintf(stderr,"\t%d-%d:0",p,flen);
      }

      depthhistogram[0]+=flen-p;

      p=flen;

    }
  }
  if(debug>0){
    fprintf(stderr,"\n");
  }


}

void process_frg(){
  int i,j,k;

  smpnfrgs[psmp]++;
  smpbases[psmp]+=flen;

  if(debug>0){
    fprintf( stderr,"%d",prev);
    for(i=1;i<=maxsmp;i++){
      fprintf(stderr," %d/%d",cnt[i][0],cln[i][0]);
    }
    fprintf(stderr, " t=%d\n",t);
  }


  for(k=0;k<ncuts;k++){
    /*************************/
    // same-sample clean overlaps histogram
    /*************************/
    if(debug>0){
      fprintf(stderr,"frg %d cutoff %d clean same-sample antis = %d\n",
	      prev,cutoffs[k],cln[psmp][k]);
    }
    if(cln[psmp][k]>maxdepth){
      samesmpantiovlhistogram[psmp][k][maxdepth]++;
    } else {
      samesmpantiovlhistogram[psmp][k][cln[psmp][k]]++;
    }

    /*************************/
    // depth of coverage histograms
    /*************************/
    depth_of_coverage(allbeg[k],allend[k],alln[k],alldepthhistogram[psmp][k],'A',k);
    depth_of_coverage(otherbeg[k],otherend[k],othern[k],otherdepthhistogram[psmp][k],'O',k);
    depth_of_coverage(selfbeg[k],selfend[k],selfn[k],selfdepthhistogram[psmp][k],'S',k);
  }
  
  /*************************/
  // sample similarity fractions
  /*************************/
  for(i=1;i<=maxsmp;i++){
    if(i!=psmp){
      for(k=0;k<ncuts;k++){
	double x,y, denom;
	x=cnt[psmp][k];
	y=cnt[i][k];
	denom=x+y;
	if(denom>0){
	  self[psmp][i][k]+=x/denom;
	  other[psmp][i][k]+=y/denom;
	}
	x=cln[psmp][k];
	y=cln[i][k];
	denom=x+y;
	if(denom>0){
	    clnself[psmp][i][k]+=x/denom;
	    clnother[psmp][i][k]+=y/denom;
	}
      }
    }
  }

}

void process_null_frg(int iid){
  int s = smp[iid];
  int k;
  int flen;
  {
    uint clr_bgn,clr_end;
    getFragStore(my_frg_store,prev,FRAG_S_ALL,fsread);
    getClearRegion_ReadStruct(fsread, &clr_bgn,&clr_end, READSTRUCT_LATEST);
    flen=clr_end-clr_bgn;
  }

  smpbases[s]+=flen;
  smpnfrgs[s]++;

  for(k=0;k<ncuts;k++){
    alldepthhistogram[s][k][0]+=flen;
    otherdepthhistogram[s][k][0]+=flen;
    selfdepthhistogram[s][k][0]+=flen;
    samesmpantiovlhistogram[s][k][0]++;
  }

}

void output_sim_matrix_info(FILE* simFile){
  int i,j,k;
  double d;

  if(simFile==NULL)return;

  // output
  for(k=0;k<ncuts;k++){
    fprintf(simFile, "All %d:\n",cutoffs[k]);
    for(i=1;i<=maxsmp;i++){
      for(j=1;j<=maxsmp;j++){
	if(i==j){
	  fprintf(simFile, " self");
	} else {
	  if(other[i][j][k]>0){
	    if(self[i][j][k]*self[j][i][k]>0){
	      d=other[i][j][k]/sqrt(self[i][j][k]*self[j][i][k]);
	      if(debug>0){
		fprintf(simFile," %f=%f/sqrt(%f * %f)",d,other[i][j][k],self[i][j][k],self[j][i][k]);
	      }else{
		fprintf(simFile," %.3f",d);
	      }
	    } else {
	      if(debug>0){
		fprintf(simFile," NaN=%f/sqrt(%f * %f)",other[i][j][k],self[i][j][k],self[j][i][k]);
	      } else {
		fprintf(simFile, " NaN");
	      }
	    }
	  } else {
	    fprintf(simFile, " 0");
	  }
	}
      }
      fprintf(simFile, "\n");
    }
    fprintf(simFile, "\n");
    fprintf(simFile, "\n");
    fprintf(simFile, "\n");

    fprintf(simFile, "Good %d:\n",cutoffs[k]);
    for(i=1;i<=maxsmp;i++){
      for(j=1;j<=maxsmp;j++){
	if(i==j){
	  fprintf(simFile, " self");
	} else {
	  if(clnother[i][j][k]>0){
	    if(clnself[i][j][k]*clnself[j][i][k]>0){
	      d=clnother[i][j][k]/sqrt(clnself[i][j][k]*clnself[j][i][k]);
	      if(debug>0){
		fprintf(simFile, " %f=%f/sqrt(%f * %f)",d,clnother[i][j][k],clnself[i][j][k],clnself[j][i][k]);
	      }else{
		fprintf(simFile," %.3f",d);
	      }
	    } else {
	      if(debug>0){
		printf (" NaN=%f/sqrt(%f * %f)",clnother[i][j][k],clnself[i][j][k],clnself[j][i][k]);
	      } else {
		fprintf(simFile, " NaN");
	      }
	    }
	  } else {
	    fprintf(simFile, " 0");
	  }
	}
      }
      fprintf(simFile, "\n");
    }
    fprintf(simFile, "\n");
    fprintf(simFile, "\n");
    fprintf(simFile, "\n");

  }
}

void output_binary_cvg_info(FILE *binaryCvgFile){
  int i,j,k;
  for(i=1;i<=maxsmp;i++){
    fprintf(binaryCvgFile,"Sample %d binary coverage:\n",i);
    for(k=0;k<ncuts;k++){
      int nonzeroA=0;
      int nonzeroO=0;
      int nonzeroS=0;
      for(j=1;j<=maxdepth;j++){
	nonzeroA += alldepthhistogram[i][k][j];
	nonzeroO += otherdepthhistogram[i][k][j];
	nonzeroS += selfdepthhistogram[i][k][j];
      }
      assert(smpbases[i]-nonzeroA==alldepthhistogram[i][k][0]);
      assert(smpbases[i]-nonzeroO==otherdepthhistogram[i][k][0]);
      assert(smpbases[i]-nonzeroS==selfdepthhistogram[i][k][0]);
      if(debug>0){
	fprintf(binaryCvgFile,"%d\t%d:\t%d/%d=%f\t%d/%d=%f\t%d/%d=%f\n",
		i,
		cutoffs[k],
		nonzeroS,smpbases[i],
		((double)nonzeroS)/((double)smpbases[i]),
		nonzeroO,smpbases[i],
		((double)nonzeroO)/((double)smpbases[i]),
		nonzeroA,smpbases[i],
		((double)nonzeroA)/((double)smpbases[i]));
      } else {
	fprintf(binaryCvgFile,"%d\t%d:\t%f\t%f\t%f\n",
		i,
		cutoffs[k],
		((double)nonzeroS)/((double)smpbases[i]),
		((double)nonzeroO)/((double)smpbases[i]),
		((double)nonzeroA)/((double)smpbases[i]));
      }
    }
  }
}

void output_depth_of_cvg_info(FILE *depthCvgFile){
  int i,j,k;
  for(i=1;i<=maxsmp;i++){
    fprintf(depthCvgFile,"Sample %d depth of coverage:\n",i);
    for(k=0;k<ncuts;k++){
      fprintf(depthCvgFile,"\t%d:\n",cutoffs[k]);
      int ttlA=0;
      int ttlO=0;
      int ttlS=0;
      for(j=0;j<=maxdepth;j++){
	ttlA += alldepthhistogram[i][k][j];
	ttlO += otherdepthhistogram[i][k][j];
	ttlS += selfdepthhistogram[i][k][j];

	if(debug>0){
	  fprintf(depthCvgFile,"%d\t%d\t%d\t%d/%d=%f\t%d/%d=%f\t%d/%d=%f\n",
		  i,cutoffs[k],
		  j,
		  selfdepthhistogram[i][k][j],smpbases[i],
		  ((double)selfdepthhistogram[i][k][j])/((double)smpbases[i]),
		  otherdepthhistogram[i][k][j],smpbases[i],
		  ((double)otherdepthhistogram[i][k][j])/((double)smpbases[i]),
		  alldepthhistogram[i][k][j],smpbases[i],
		  ((double)alldepthhistogram[i][k][j])/((double)smpbases[i]));
	} else {
	  fprintf(depthCvgFile,"%d\t%d\t%d\t%f\t%f\t%f\n",
		  i,cutoffs[k],
		  j,
		  ((double)selfdepthhistogram[i][k][j])/((double)smpbases[i]),
		  ((double)otherdepthhistogram[i][k][j])/((double)smpbases[i]),
		  ((double)alldepthhistogram[i][k][j])/((double)smpbases[i]));
	}

      }
      assert(smpbases[i]==ttlA);
      assert(smpbases[i]==ttlO);
      assert(smpbases[i]==ttlS);

    }
  }
}


void   output_antiOvl_histos(FILE *antiOvlHistoFile){
  int i,j,k;
  for(i=1;i<=maxsmp;i++){
    fprintf(antiOvlHistoFile,"Sample %d anti-ovl depths on %d frgs:\n",i,smpnfrgs[i]);
    for(k=0;k<ncuts;k++){
      fprintf(antiOvlHistoFile,"\t%d:\n",cutoffs[k]);
      int ttl=0;
      for(j=0;j<=maxdepth;j++){
	ttl += samesmpantiovlhistogram[i][k][j];

	fprintf(antiOvlHistoFile,
		"%d\t%d\t%d\t%d\n",
		i,
		cutoffs[k],
		j,
		samesmpantiovlhistogram[i][k][j]);

      }
      assert(smpnfrgs[i]==ttl);
    }
  }
}


void alloc_sim_arrays(){
  int i,j,k;
  self = (double ***) safe_malloc(sizeof(double**)*(maxsmp+1));
  other = (double ***) safe_malloc(sizeof(double**)*(maxsmp+1));
  clnself = (double ***) safe_malloc(sizeof(double**)*(maxsmp+1));
  clnother = (double ***) safe_malloc(sizeof(double**)*(maxsmp+1));
  cnt = (int**) safe_malloc(sizeof(int*)*(maxsmp+1));
  cln = (int**) safe_malloc(sizeof(int*)*(maxsmp+1));

  smpnfrgs = (int *) safe_malloc(sizeof(int)*(maxsmp+1));
  smpbases = (int *) safe_malloc(sizeof(int)*(maxsmp+1));

  for(i=0;i<=maxsmp;i++){
    self[i]=(double**)safe_malloc(sizeof(double*)*(maxsmp+1));
    other[i]=(double**)safe_malloc(sizeof(double*)*(maxsmp+1));
    clnself[i]=(double**)safe_malloc(sizeof(double*)*(maxsmp+1));
    clnother[i]=(double**)safe_malloc(sizeof(double*)*(maxsmp+1));
    for(j=0;j<=maxsmp;j++){
      self[i][j]=(double*)safe_malloc(sizeof(double)*(ncuts));
      other[i][j]=(double*)safe_malloc(sizeof(double)*(ncuts));
      clnself[i][j]=(double*)safe_malloc(sizeof(double)*(ncuts));
      clnother[i][j]=(double*)safe_malloc(sizeof(double)*(ncuts));
      for(k=0;k<ncuts;k++){
	self[i][j][k]=0;
	other[i][j][k]=0;
	clnself[i][j][k]=0;
	clnother[i][j][k]=0;
      }
    }
    cnt[i]=(int*)safe_malloc(sizeof(int)*ncuts);
    cln[i]=(int*)safe_malloc(sizeof(int)*ncuts);
    smpnfrgs[i] = 0;
    smpbases[i] = 0;
  }

}

void alloc_depths(){
  int i,k,j;
  samesmpantiovlhistogram = (int ***) safe_malloc((maxsmp+1)*sizeof(int**));
  alldepthhistogram = (int ***) safe_malloc((maxsmp+1)*sizeof(int**));
  otherdepthhistogram = (int ***) safe_malloc((maxsmp+1)*sizeof(int**));
  selfdepthhistogram = (int ***) safe_malloc((maxsmp+1)*sizeof(int**));
  for(i=0;i<=maxsmp;i++){
    samesmpantiovlhistogram[i] = (int **) safe_malloc(ncuts*sizeof(int*));
    alldepthhistogram[i] = (int **) safe_malloc(ncuts*sizeof(int*));
    otherdepthhistogram[i] = (int **) safe_malloc(ncuts*sizeof(int*));
    selfdepthhistogram[i] = (int **) safe_malloc(ncuts*sizeof(int*));
    for(k=0;k<ncuts;k++){
      samesmpantiovlhistogram[i][k] = (int*) safe_malloc((maxdepth+1)*sizeof(int));
      alldepthhistogram[i][k] = (int*) safe_malloc((maxdepth+1)*sizeof(int));
      otherdepthhistogram[i][k] = (int*) safe_malloc((maxdepth+1)*sizeof(int));
      selfdepthhistogram[i][k] = (int*) safe_malloc((maxdepth+1)*sizeof(int));
      for(j=0;j<=maxdepth;j++){
	samesmpantiovlhistogram[i][k][j]=0;
	alldepthhistogram[i][k][j]=0;
	otherdepthhistogram[i][k][j]=0;
	selfdepthhistogram[i][k][j]=0;
      }
    }
  }
}

void alloc_ends(int initmax){

  int i,k;
  assert(maxovls==0);
  allbeg = (int **) safe_malloc( sizeof(int *)*ncuts);
  allend = (int **) safe_malloc( sizeof(int *)*ncuts);
  alln = (int *) safe_malloc( sizeof(int)*ncuts );
  selfbeg = (int **) safe_malloc( sizeof(int *)*ncuts);
  selfend = (int **) safe_malloc( sizeof(int *)*ncuts);
  selfn = (int *) safe_malloc( sizeof(int)*ncuts );
  otherbeg = (int **) safe_malloc( sizeof(int *)*ncuts);
  otherend = (int **) safe_malloc( sizeof(int *)*ncuts);
  othern = (int *) safe_malloc( sizeof(int)*ncuts );

  for(k=0;k<ncuts;k++){
    allbeg[k] = (int *) safe_malloc (  sizeof(int)*initmax );
    allend[k] = (int *) safe_malloc (  sizeof(int)*initmax );
    selfbeg[k] = (int *) safe_malloc ( sizeof(int)*initmax );
    selfend[k] = (int *) safe_malloc ( sizeof(int)*initmax );
    otherbeg[k] = (int *) safe_malloc (sizeof(int)*initmax );
    otherend[k] = (int *) safe_malloc (sizeof(int)*initmax );
  }
  maxovls=initmax;

}

void realloc_ends(int newmax){
  int k;
  assert(maxovls>0);
  for(k=0;k<ncuts;k++){
    allbeg[k] = (int *) safe_realloc ( allbeg[k], sizeof(int)*newmax );
    allend[k] = (int *) safe_realloc ( allend[k], sizeof(int)*newmax );
    selfbeg[k] = (int *) safe_realloc ( selfbeg[k], sizeof(int)*newmax );
    selfend[k] = (int *) safe_realloc ( selfend[k], sizeof(int)*newmax );
    otherbeg[k] = (int *) safe_realloc ( otherbeg[k], sizeof(int)*newmax );
    otherend[k] = (int *) safe_realloc ( otherend[k], sizeof(int)*newmax );
  }
  maxovls=newmax;
}

int main (int argc, char *argv[]){

  char full_ovlPath[1000];
  char full_frgPath[1000];
  char full_gkpPath[1000];
  char sampleFileName[1000];
  int setFullFrg=0;
  int setFullGkp=0;
  int setFullOvl=0;
  uint32  lastfrag,lastovlfrg;
  char *smpsFile=NULL;
  FILE *smps;
  CDS_IID_t iid1,iid2;
  int errflg=0;
  char ch;
  int i,j,k;
  int smpIdx;
  int clean;
  int sb;
  double d;
  int begidx=-1,endidx=-1;
  FILE *simFile=NULL;
  FILE *binaryCvgFile=NULL;
  FILE *depthCvgFile=NULL;
  FILE *antiOvlHistoFile=NULL;
  fsread=new_ReadStruct();

  while  (! errflg
	  && ((ch = getopt (argc, argv, "b:dDe:f:g:m:o:s:A:B:C:S:"))!=EOF)){
    switch  (ch)
      {
      case 'b':
	begidx=atoi(optarg);
	break;
      case 'd':
	debug=1;
	break;
      case 'D':
	debug=2;
	break;
      case 'e':
	endidx=atoi(optarg);
	break;
      case 'f':
	strcpy( full_frgPath, argv[optind - 1]);
	setFullFrg = TRUE;
	break;
      case 'g':
	strcpy( full_gkpPath, argv[optind - 1]);
	setFullGkp = TRUE;
	break;
      case 'o':
	strcpy(full_ovlPath,argv[(optind-1)]);
	setFullOvl = TRUE;
	break;
      case 's':
	smpsFile=argv[(optind-1)];
	break;
      case 'A':
	antiOvlHistoFile = fopen(argv[(optind-1)],"w");
	break;
      case 'B':
	binaryCvgFile = fopen(argv[(optind-1)],"w");
	break;
      case 'C':
	depthCvgFile = fopen(argv[(optind-1)],"w");
	break;
      case 'S':
	simFile = fopen(argv[(optind-1)],"w");
	break;
      default:
	errflg++;
      }
  }

  if(errflg ||
     setFullOvl!=1 ||
     setFullGkp!=1 ||
     setFullFrg!=1
     ){
    usage(argv[0]);
  }

  setup_stores(full_ovlPath,full_frgPath,full_gkpPath);
  lastfrag = getLastElemFragStore (my_frg_store);

  assert(smpsFile!=NULL);

  if(begidx==-1)begidx=1;
  lastovlfrg = Last_Frag_In_OVL_Store (my_ovl_store);
  if(endidx==-1)endidx=lastovlfrg;
  assert(begidx<=endidx);
  assert(endidx<=lastovlfrg);
  smp=(int*)safe_malloc(sizeof(int)*(lastfrag+1));

  maxsmp=0;
  smps=fopen(smpsFile,"r");
  while(fscanf(smps,F_IID " %d",&iid1,&smpIdx)==2){
    smp[iid1]=smpIdx;
    if(smpIdx>maxsmp){
      maxsmp=smpIdx;
    }
  }
  fclose(smps);
  fprintf( stderr, "read IID2SMP -- maxsmp = %d\n",maxsmp);


  alloc_sim_arrays();

  alloc_depths();

  alloc_ends(5000);

  
  prev=0;

  Init_OVL_Stream (my_stream, begidx, endidx, my_ovl_store);

  while  (Next_From_OVL_Stream (& olap, my_stream)){
    if(debug>0){
      fprintf (stderr,"    %8d %8d %c %5d %5d %4.1f %4.1f\n",
              olap . a_iid, olap . b_iid,
              olap . flipped ? 'I' : 'N',
              olap . a_hang, olap . b_hang,
              Expand_Quality (olap . orig_erate) * 100.0,
              Expand_Quality (olap . corr_erate) * 100.0);
    }

    if(olap.a_iid!=prev){
      if(prev!=0&&t>0){
	process_frg();
	prev++;
      }
      if(prev==0){prev=1;}
      for(i=prev;i<olap.a_iid;i++){
	if(debug>0){
	  fprintf(stderr,"No used overlaps for %d\n",i);
	}
	process_null_frg(i);
      }
      init_frg();
    }
    if (currmate == olap.b_iid) continue;

    clean=0;

    if( olap.flipped && olap.a_hang < 0 && olap.b_hang < 0 && currmate != olap.b_iid) {
      clean=1;
    }

    sb=smp[olap.b_iid];

    for(k=0;k<ncuts;k++){

      if(olap.orig_erate>cutoffs[k]){break;}

      cnt[sb][k]++;

      if(clean){
	cln[sb][k]++;
      }

      { 
	int b,e;
	b = ( olap.a_hang>0 ? olap.a_hang : 0 );
	e = flen + ( olap.b_hang < 0 ? olap.b_hang : 0 );

	assert(b<e);
	assert(b>=0&&b<flen);
	assert(e>0&&e<=flen);

	if(debug>0){
	  fprintf(stderr,"=> [%d,%d] asmp=%d bsmp=%d\n",b,e,psmp,sb);
	}

        allbeg[k][alln[k]]=b;
	allend[k][alln[k]]=e;
	alln[k]++;

	if(sb==psmp){
	  selfbeg[k][selfn[k]]=b;
	  selfend[k][selfn[k]]=e;
	  selfn[k]++;
	} else {
	  otherbeg[k][othern[k]]=b;
	  otherend[k][othern[k]]=e;
	  othern[k]++;
	}

	if(alln[k]>maxovls){
	  realloc_ends(maxovls*2);
	}

      }

    }
    if(k>t)t=k;

  }
  Free_OVL_Stream (my_stream);
  Free_OVL_Store (my_ovl_store);


  // handle last overlapped frag plus those following...
  if(prev!=0&&t>0){
    process_frg();
    prev++;
  }
  for(i=prev;i<=lastfrag;i++){
    if(debug>0){
      fprintf(stderr,"No used overlaps for %d\n",i);
    }
    process_null_frg(i);
  }

  output_sim_matrix_info(simFile);
  output_binary_cvg_info(binaryCvgFile);
  output_depth_of_cvg_info(depthCvgFile);
  output_antiOvl_histos(antiOvlHistoFile);

  return(0);

}

