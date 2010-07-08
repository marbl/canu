//Copyright (c) 2003 by Mihaela Pertea


#include "sim4.H"
#include "sites_score.h"
#include "sites_donor.h"
#include "sites_acceptor.h"

char DONOR_TREE[] = "( 0 2 4 10000 l( 1 2 9 7841 l( 3 0 8 5666 l( 5 0 3 3977 l( 7 0 7 2186 l( 9 -1 -1 995 l r ) r( 10 -1 -1 1191 l r ) ) r( 8 3 10 1791 l( 15 -1 -1 931 l r ) r( 16 -1 -1 860 l r ) ) ) r( 6 -1 -1 1689 l r ) ) r( 4 -1 -1 2175 l r ) ) r( 2 -1 -1 2159 l r ) )"; // \n5 20\n";

char ACCEPTOR_TREE[] = "( 0 1 23 10000 l( 1 3 21 6544 l( 3 3 20 3573 l( 5 3 16 2146 l( 7 -1 -1 1295 l r ) r( 8 -1 -1 851 l r ) ) r( 6 -1 -1 1427 l r ) ) r( 4 1 21 2971 l( 15 1 20 1914 l( 17 -1 -1 1009 l r ) r( 18 -1 -1 905 l r ) ) r( 16 -1 -1 1057 l r ) ) ) r( 2 -1 -1 3456 l r ) )"; // \n44 72\n";

#define  TRUE  1
#define  FALSE  0
#define  ACCEPTOR_LEN  29                    /* Positions +44,72 in a80 */
#define  ACCEPTOR_SIGNAL_OFFSET  24          /* Start of  AG  */

#define  DONOR_LEN  16                        /* Positions +5,20 in d80 */
#define  DONOR_SIGNAL_OFFSET  5               /* Start of  GT  */

#define  MARKOV_DEGREE  3
#define  MARKOV_LEN  64                     /* ALPHABET_SIZE ^ MARKOV_DEGREE */
#define  LOW_SCORE  -99.0  /* Score if pattern does not have GT or AG signal */


#define SITE_LEN 162

#define CODING_LEN 80

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

typedef struct tree {
	int val;
   	int consens;
   	int poz;
   	int no;
   	struct tree *left;
   	struct tree *right;
   } tree;

void postorder(tree *root)
{
  if(root)
  {
    postorder(root->left);
    postorder(root->right);
    printf("[%d %d %d %d] ", root->val, root->consens, root->poz, root->no);
  }
}

typedef unsigned int word;

int  Acc  (const int *, double *,tree *t,int ind);
int  Don  (const int *, double *, tree *t,int ind);
int  comp(const void *a, const void *b);
int  findfile(const int * S, tree *t);
int  readtree(char *line, tree *t, int start);
int  find(char *line, int start);
int  Is_Cod_NonCod  (const int * , double *, int ind);
float ****Load4dim(int d1, int d2, int d3, int d4);
void free4dim(float ****ptr,int d1, int d2, int d3);

#define  Start_PosEx 56
#define  Stop_PosEx 84

#define  Start_PosIn 75
#define  Stop_PosIn 90

#define  Start_Cod 0
#define  Stop_Cod 79

#define Start_NoCod 82
#define Stop_NoCod 161


int markov_degree;
int markov_len;
tree *tacc = NULL;
tree *tdon = NULL;
int readtacc=FALSE;
int readtdon=FALSE;
int accmax = 0;
int donmax = 0;
float  ****Acc_Positive_Table = NULL;
float  ****Acc_Negative_Table = NULL;
int *Acc_Tables_Loaded = NULL;
float  ****Don_Positive_Table = NULL;
float  ****Don_Negative_Table = NULL;
int *Don_Tables_Loaded = NULL;
float  Cod_Positive_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
float  Cod_Negative_Table [4][CODING_LEN] [ALPHABET_SIZE] [MARKOV_LEN];
int  Cod_Tables_Loaded[4] = {FALSE,FALSE,FALSE,FALSE};

void
Sim4::loadGeneSplicerModel()
{
  int i;

  markov_degree=1;
  markov_len=(int)pow(ALPHABET_SIZE,1);

  if(!readtdon) {
   
    tdon = (tree *) ckalloc(sizeof(tree));
    if (tdon == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
   
    donmax=readtree(DONOR_TREE, tdon, 0);
    readtdon=TRUE;

    // alloc memory for the tables
    Don_Positive_Table=Load4dim(donmax,DONOR_LEN,ALPHABET_SIZE,markov_len);
    Don_Negative_Table=Load4dim(donmax,DONOR_LEN,ALPHABET_SIZE,markov_len);
    Don_Tables_Loaded=(int *) ckalloc(donmax*sizeof(int));
    if(Don_Tables_Loaded == NULL) {
      fprintf(stderr,"Memory allocation for donor site tables failed.\n");
      abort();
    }
    for(i=0;i<donmax;i++) Don_Tables_Loaded[i]=FALSE;
  }

  if(!readtacc) {
   
    // read the structure of the acceptor tree
    tacc = (tree *) ckalloc(sizeof(tree));
    if (tacc == NULL) {fprintf(stderr," Memory allocation for tree failure.\n"); abort();}
    accmax=readtree(ACCEPTOR_TREE, tacc, 0);

#ifdef DEBUG
    printf("readtacc = %d when readtacc should be 0\n", readtacc);
    printf("accmax = %d\n", accmax);
    postorder(tacc);
    printf("\n");
#endif

    readtacc=TRUE;
   
    // alloc memory for the tables
    Acc_Positive_Table=Load4dim(accmax,ACCEPTOR_LEN,ALPHABET_SIZE,markov_len);
    Acc_Negative_Table=Load4dim(accmax,ACCEPTOR_LEN,ALPHABET_SIZE,markov_len);
    Acc_Tables_Loaded=(int *) ckalloc(accmax*sizeof(int));
    if(Acc_Tables_Loaded == NULL) {
      fprintf(stderr,"Memory allocation for acceptor site tables failed.\n");
      abort();
    }
    for(i=0;i<accmax;i++) Acc_Tables_Loaded[i]=FALSE;
  }
}

void free4dim(float ****ptr,int d1, int d2, int d3)
{
  int i,j,k;

  for(i=0;i<d1;i++) {
    for(j=0;j<d2;j++) {
      for(k=0;k<d3;k++) {
	if(ptr[i][j][k] != NULL ) 
          ckfree(ptr[i][j][k]);
      }
      if(ptr[i][j] != NULL ) 
        ckfree(ptr[i][j]);
    }
    if(ptr[i] != NULL ) 
      ckfree(ptr[i]);
  }
  ckfree(ptr);
}


void freetree(tree *t)
{
  if(t==NULL) return;
  freetree(t->left);
  freetree(t->right);
  ckfree(t);
  t=NULL;
}

void
Sim4::UnLoadSites_GeneSplicer()
{
  int i;

  if(readtacc) {
    free4dim(Acc_Positive_Table,accmax,ACCEPTOR_LEN,ALPHABET_SIZE);
    free4dim(Acc_Negative_Table,accmax,ACCEPTOR_LEN,ALPHABET_SIZE);
    if(Acc_Tables_Loaded  != NULL ) ckfree(Acc_Tables_Loaded);
  }

  if(readtdon) {
    free4dim(Don_Positive_Table,donmax,DONOR_LEN,ALPHABET_SIZE);
    free4dim(Don_Negative_Table,donmax,DONOR_LEN,ALPHABET_SIZE);
    if(Don_Tables_Loaded != NULL ) ckfree(Don_Tables_Loaded);
  }

#ifdef DEBUG
  printf("tacc:\n");
  postorder(tacc);
  printf("\n");
#endif

  if(readtacc) 
    freetree(tacc);

#ifdef DEBUG
  printf("tdon:\n");
  postorder(tdon);
  printf("\n");
#endif

  if(readtdon)
    freetree(tdon);

  readtacc=FALSE;
  readtdon=FALSE;

  for(i=0;i<4;i++) 
    Cod_Tables_Loaded[i]=FALSE;

}

float ****Load4dim(int d1, int d2, int d3, int d4)
{
  int i,j,k;
  float ****ptr;
  
  ptr = (float ****) ckalloc(d1 * sizeof(float ***));
  if(ptr==NULL) {
    fprintf(stderr,"Memory allocation for splice site tables failed.\n");
    abort();
  }
  for(i=0;i<d1;i++) {
    ptr[i] = (float ***) ckalloc(d2 * sizeof(float **));
    if(ptr[i]==NULL) {
      fprintf(stderr,"Memory allocation for splice site tables failed.\n");
      abort();
    }
    for(j=0;j<d2;j++) {
      ptr[i][j] = (float **) ckalloc(d3*sizeof(float *));
      if(ptr[i][j]==NULL) {
	fprintf(stderr,"Memory allocation for splice site tables failed.\n");
	abort();
      }
      for(k=0;k<d3;k++) {
	ptr[i][j][k] = (float *) ckalloc(d4*sizeof(float));
	if(ptr[i][j][k]==NULL) {
	  fprintf(stderr,"Memory allocation for splice site tables failed.\n");
	  abort();
	}
      }
    }
  }

  return(ptr);
}


double
Sim4::ScoreAcceptor_GeneSplicer(char *Data)
{
  double Score,S1,S2;
  int i,ind;
  int T[100];
  double score1,score2,score3;
  char *B = Data;  

#if 0
  assert( strlen(Data) >= SITE_LEN);
  
  for(i=0;i<SITE_LEN;i++) {
    switch (Data[i]){
    case 'A':
    case 'a': B[i]=0;break;
    case 'C':
    case 'c': B[i]=1;break;
    case 'G':
    case 'g': B[i]=2;break;
    case 'T': 
    case 't': B[i]=3;break;
    default: B[i]=0;
    }
  }
#endif

#if 0
  /* moved to loadGeneSplicerModel */
  markov_degree=1;
  markov_len=(int)pow(ALPHABET_SIZE,1);

  if(!readtacc) {
    
    // read the structure of the acceptor tree 
    tacc = (tree *) ckalloc(sizeof(tree));
    if (tacc == NULL) {fprintf(stderr," Memory allocation for tree failure.\n"); abort();}
    accmax=readtree(ACCEPTOR_TREE, tacc, 0);

#ifdef DEBUG
    printf("readtacc = %d when readtacc should be 0\n", readtacc);
    printf("accmax = %d\n", accmax);
    postorder(tacc);
    printf("\n");
#endif    

    readtacc=TRUE;
    
    // alloc memory for the tables
    Acc_Positive_Table=Load4dim(accmax,ACCEPTOR_LEN,ALPHABET_SIZE,markov_len);
    Acc_Negative_Table=Load4dim(accmax,ACCEPTOR_LEN,ALPHABET_SIZE,markov_len);
    Acc_Tables_Loaded=(int *) ckalloc(accmax*sizeof(int));
    if(Acc_Tables_Loaded == NULL) {
      fprintf(stderr,"Memory allocation for acceptor site tables failed.\n");
      abort();
    }
    for(i=0;i<accmax;i++) Acc_Tables_Loaded[i]=FALSE;      
  }
#endif


  for(i=0;i<=Stop_PosEx-Start_PosEx;i++)
    T[i]=B[i+Start_PosEx];

  ind=Acc(T, &S1, tacc,0);
  if(ind==0) return(0);

  if(accmax>1) Acc(T, &S2, tacc,1);
  else S2=S1;
  score1=(S1+S2)/2;

  //  if(score1<=THR_ACC) score1=-99;
  
  score2=0;
  score3=0;
    
  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];

  Is_Cod_NonCod(T,&score2,0);

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];
  

  Is_Cod_NonCod(T,&score3,1);
    
//  printf("score1 = %.5f, score2 = %.5f, score3 = %.5f\n", score1, score2, score3);
  Score=score1+score2+score3;

  return(Score);
	  
      
}  

double
Sim4::ScoreDonor_GeneSplicer(char *Data)
{
  double Score,S1,S2;
  int ind,i;
  int T[100];
  double score1,score2,score3;
  char *B = Data;  

#if 0
  assert( strlen(Data) >= SITE_LEN);
  
  for(i=0;i<SITE_LEN;i++) {
    switch (Data[i]){
    case 'A': 
    case 'a': B[i]=0;break;
    case 'C':
    case 'c': B[i]=1;break;
    case 'G':
    case 'g': B[i]=2;break;
    case 'T':
    case 't': B[i]=3;break;
    default: B[i]=0;
    }
  }
#endif

#if 1
  /* LLL moved to loadGeneSplicerModel */ 
  markov_degree=1;
  markov_len=(int)pow(ALPHABET_SIZE,1);

  if(!readtdon) {
    
    tdon = (tree *) ckalloc(sizeof(tree));
    if (tdon == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
    
    donmax=readtree(DONOR_TREE, tdon, 0);
    readtdon=TRUE;

    // alloc memory for the tables
    Don_Positive_Table=Load4dim(donmax,DONOR_LEN,ALPHABET_SIZE,markov_len);
    Don_Negative_Table=Load4dim(donmax,DONOR_LEN,ALPHABET_SIZE,markov_len);
    Don_Tables_Loaded=(int *) ckalloc(donmax*sizeof(int));
    if(Don_Tables_Loaded == NULL) {
      fprintf(stderr,"Memory allocation for donor site tables failed.\n");
      abort();
    }
    for(i=0;i<donmax;i++) Don_Tables_Loaded[i]=FALSE;   
  }
#endif

  for(i=0;i<=Stop_PosIn-Start_PosIn;i++)
    T[i]=B[i+Start_PosIn];

  ind=Don(T, &S1, tdon,0);
  if(ind==0) return(0);
  if(donmax>1) Don(T, &S2, tdon,1);
  else S2=S1;
  score1=(S1+S2)/2;


  score2=0;
  score3=0;

  for(i=0;i<=Stop_Cod-Start_Cod;i++)
    T[i]=B[i+Start_Cod];

  Is_Cod_NonCod(T,&score2,2);

  
  for(i=0;i<=Stop_NoCod-Start_NoCod;i++)
    T[i]=B[i+Start_NoCod];
  
  Is_Cod_NonCod(T,&score3,3);


  Score=score1+score2+score3;

  return Score;
	  
      
}  


    
int readtree(char *line, tree *t, int start)
{
 int len;
 int i,n;
 int val,valmax;
 char part[10];
 len=strlen(line);

 i=start;
 while((line[i]=='(')||(line[i]==' ')) i++;
 n=i;
 while(line[i]!=' ')
 {
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->val=atoi(part);
 valmax=t->val;

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->consens=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->poz=atoi(part);

 i++;
 n=i;
 while(line[i]!=' ')
 { 
	part[i-n]=line[i];
	i++;
 }
 part[i-n]='\0';
 t->no=atoi(part);

 t->left=NULL;
 t->right=NULL;

 i+=2;n=i;
 if(line[i]=='(') 
   {
     i=find(line,i+1);
     t->left = (tree *) ckalloc(sizeof(tree));
     if (t->left == NULL) {fprintf(stderr,"Memory allocation for tree failure.\n"); abort();}
     val=readtree(line,t->left,n);
     if(val>valmax) valmax=val;
   }
	
 i+=2;n=i;
 if(line[i]=='(') 
   {
     i=find(line,i+1);
     t->right = (tree *) ckalloc(sizeof(tree));
     if (t->right == NULL) {
       fprintf(stderr,"Memory allocation for tree failure.\n"); 
       abort();
     }
     val=readtree(line,t->right,n);
     if(val>valmax) valmax=val;
   }
 valmax++;
 return(valmax);
}

int find(char *line, int start)
{
 int stop,i;

 i=start;

 while(line[i]!=')')
 	if(line[i]=='(') i=find(line,i+1);
 	else i++;
 stop=i+1;
 return(stop);
}
 	

int comp(const void *a, const void *b)
{ 
  if(*(double *)a > *(double *)b) return(1);
  else if (*(double *)a==*(double *)b) return(0);
  else return(-1);

}  
  

int findfile(const int * S, tree *t)
{
	int val, cons, poz;
	val=t->val;

	cons=t->consens;
	if( cons !=-1)
	{ 
		poz=t->poz;
	    if(S[poz]==cons)
	    	val=findfile(S,t->left);
	    else val=findfile(S, t->right);
	}

	return(val);
}

int findleaf(tree *t, int n, int leaf, int *found) 
{
  int ret=n;

  if(t==NULL) { fprintf(stderr,"tree NULL\n");exit(0);}

  if(t->val == leaf) {*found=1; return(n+1);}

  if(t->left == NULL && t->right == NULL)  return(n+1);
  if(t->left != NULL) ret=findleaf(t->left,n,leaf,found);
  if(!(*found) && t->right != NULL) ret=findleaf(t->right,ret,leaf,found);
  
  return(ret);
}
   
  
  



int  Acc  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (ACCEPTOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely acceptor
*  site.  Also set  Return_Score  to the probability that it is an acceptor
*  site. */

{
  double  Positive_Sum, Negative_Sum, Score;
#if  RETURN_TRUE_PROB
  double  X, Y;
#endif
  int  i, j, k, Sub, no, idx;

/* see which acceptor you should use */

  if(ind) {
    no=findfile(S,t);
    k=0;
  }
  else
    no=0;
  
  idx = 0;
  if  (! Acc_Tables_Loaded[no])
    {
      for  (i = markov_degree - 1;  i < ACCEPTOR_LEN;  i ++)
	for  (k = 0;  k < markov_len;  k ++)
	  for  (j = 0;  j < ALPHABET_SIZE;  j ++)
	    {
              Acc_Positive_Table[no][i][j][k] = acc[no][idx++];
	    }
      
      for  (i = markov_degree - 1;  i < ACCEPTOR_LEN;  i ++)
	for  (k = 0;  k < markov_len;  k ++)
	  for  (j = 0;  j < ALPHABET_SIZE;  j ++)
	    {
              Acc_Negative_Table[no][i][j][k] = acc[no][idx++];
	    }
      
      Acc_Tables_Loaded[no]  = TRUE;
    }
  

  /*
  if  (S [ACCEPTOR_SIGNAL_OFFSET] != 0
  || S [ACCEPTOR_SIGNAL_OFFSET + 1] != 2)    // AG
  {
    * Return_Score = LOW_SCORE;
    return  FALSE;
  }
  */
    
  Sub = 0;
  for  (i = 0;  i < markov_degree;  i ++)
    Sub = ALPHABET_SIZE * Sub + S [i];
  
  Positive_Sum = Acc_Positive_Table [no][markov_degree - 1] [0] [Sub];
  Negative_Sum = Acc_Negative_Table [no][markov_degree - 1] [0] [Sub];
  
  for  (i = markov_degree;  i < ACCEPTOR_LEN;  i ++)
    {
      j = S [i];
      Positive_Sum += Acc_Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Acc_Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (markov_len / ALPHABET_SIZE)) + j;
    }
  


  Score = Positive_Sum - Negative_Sum;

   * Return_Score = Score;

   return(1);
  }



int  Don  (const int * S, double * Return_Score, tree *t,int ind)

/* Evaluate string  S [0 .. (DONOR_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */
{
   double  Positive_Sum, Negative_Sum, Score;
   int no;

#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Sub, idx;

   /* see which donor file you should use */
   if(ind) {
     no=findfile(S,t);
     k=0;
   }
   else 
     no=0;

   idx = 0;
   if  (! Don_Tables_Loaded[no] )
       {
        for  (i = markov_degree - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < markov_len;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Don_Positive_Table[no][i][j][k] = don[no][idx++]; 
              }

        for  (i = markov_degree - 1;  i < DONOR_LEN;  i ++)
          for  (k = 0;  k < markov_len;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Don_Negative_Table[no][i][j][k] = don[no][idx++];
              }
        Don_Tables_Loaded [no] = TRUE;
       }

   /*
   if  (S [DONOR_SIGNAL_OFFSET] != 2
           || S [DONOR_SIGNAL_OFFSET + 1] != 3)    // GT 
       {
        * Return_Score = LOW_SCORE;
        return  FALSE;
       }
   */

   Sub = 0;
   for  (i = 0;  i < markov_degree;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Don_Positive_Table [no] [markov_degree - 1] [0] [Sub];
   Negative_Sum = Don_Negative_Table [no] [markov_degree - 1] [0] [Sub];

   for  (i = markov_degree;  i < DONOR_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Don_Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Don_Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (markov_len / ALPHABET_SIZE)) + j;
     }
 
   Score = Positive_Sum - Negative_Sum;

   * Return_Score = Score;

   return(1);
  }


int  Is_Cod_NonCod  (const int * S, double * Return_Score, int ind)

/* Evaluate string  S [0 .. (CODING_LEN -1)] and
*  return  TRUE  or  FALSE  as to whether it is a likely donor
*  site.  Also set  Return_Score  to the probability that it is an donor
*  site. */

  {
   double  Positive_Sum, Negative_Sum, Score;
   double *scores;
   int no;


#if  RETURN_TRUE_PROB
   double  X, Y;
#endif
   int  i, j, k, Sub, idx;

   no=ind;

   switch (no) {
   case 0: // case of exon in acceptor
     scores = score_ex_acc;
     break;
   case 1: // case of intron in acceptor
     scores = score_in_acc;
     break;
   case 2: // case of exon in donor
     scores = score_ex_don;
     break;
   case 3: // case of intron in donor
     scores = score_in_don;
     break;
   }

   idx = 0;
   if  (! Cod_Tables_Loaded[no] )
       {
        for  (i = markov_degree - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < markov_len;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Cod_Positive_Table[no][i][j][k] = scores[idx++];
              }

        for  (i = markov_degree - 1;  i < CODING_LEN;  i ++)
          for  (k = 0;  k < markov_len;  k ++)
            for  (j = 0;  j < ALPHABET_SIZE;  j ++)
              {
               Cod_Negative_Table[no][i][j][k] = scores[idx++];
              }

        Cod_Tables_Loaded [no] = TRUE;
       }

   Sub = 0;
   for  (i = 0;  i < markov_degree;  i ++)
     Sub = ALPHABET_SIZE * Sub + S [i];

   Positive_Sum = Cod_Positive_Table [no] [markov_degree - 1] [0] [Sub];
   Negative_Sum = Cod_Negative_Table [no] [markov_degree - 1] [0] [Sub];

   for  (i = markov_degree;  i < CODING_LEN;  i ++)
     {
      j = S [i];
      Positive_Sum += Cod_Positive_Table [no] [i] [j] [Sub];
      Negative_Sum += Cod_Negative_Table [no] [i] [j] [Sub];
      Sub = ALPHABET_SIZE * (Sub % (markov_len / ALPHABET_SIZE)) + j;
     }
 


   Score = Positive_Sum - Negative_Sum;

   * Return_Score = Score;

   return (1);
  }

