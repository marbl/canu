#include <pthread.h>
#include "sim4.H"

#define GENESPLICER_SPAN    80
#define GLIMMER_XSPAN       30
#define GLIMMER_ISPAN       20
#define S4_SPAN              0
#define MAX_SPAN            80

static int   spl_encode[256] = { 0 };
static int   rev_compl[256] = { 0 };
static int   spliceInit = 0;

int const gt[5][5] = {{0, 0, 0, 2, 0},
                      {0, 0, 0, 2, 0},
                      {2, 3, 2, 5, 2},
                      {0, 0, 0, 2, 0},
                      {0, 0, 0, 2, 0}};
int const ct[5][5] = {{0, 0, 0, 2, 0},
                      {2, 2, 2, 5, 2},
                      {0, 0, 0, 2, 0},
                      {0, 0, 0, 2, 0},
                      {0, 0, 0, 2, 0}};
int const ag[5][5] = {{2, 2, 5, 2, 2},
                      {0, 0, 2, 0, 0},
                      {0, 0, 2, 0, 0},
                      {0, 0, 2, 0, 0},
                      {0, 0, 2, 0, 0}};
int const ac[5][5] = {{2, 5, 2, 2, 2},
                      {0, 2, 0, 0, 0},
                      {0, 3, 0, 0, 0},
                      {0, 2, 0, 0, 0},
                      {0, 2, 0, 0, 0}};


#if 0
int const gt[4][4] = {{0, 0, 0, 2},{0, 0, 0, 2},{2, 2, 2, 5},{0, 0, 0, 2}};
int const ct[4][4] = {{0, 0, 0, 2},{2, 2, 2, 5},{0, 0, 0, 2},{0, 0, 0, 2}};
int const ag[4][4] = {{2, 2, 5, 2},{0, 0, 2, 0},{0, 0, 2, 0},{0, 0, 2, 0}};
int const ac[4][4] = {{2, 5, 2, 2},{0, 2, 0, 0},{0, 2, 0, 0},{0, 2, 0, 0}};
#endif


/* GLIMMER functions - move to glimmer.h? */

static char  Glimmer_TRAIN_DIR[] = "./GlimmerModels/";
static char  Glimmer_posDonModelPath[] = "donors.162.pos.icm";
static char  Glimmer_negDonModelPath[] = "donors.162.neg.icm";
static char  Glimmer_posAccModelPath[] = "acceptors.162.pos.icm";
static char  Glimmer_negAccModelPath[] = "acceptors.162.neg.icm";

struct  Fixed_Length_ICM_t  donor_pos_model, donor_neg_model;
struct  Fixed_Length_ICM_t  acceptor_pos_model, acceptor_neg_model;
int     donor_pos_model_len, donor_neg_model_len;
int     acceptor_pos_model_len, acceptor_neg_model_len;
int     initGlimmerModel = 0;

void
Sim4::loadGlimmerModel (char *train_dir)
{
   char filename[1000];

   if (!initGlimmerModel) { /* LLL is this still needed? Yes, since it is initialized in the class Sim4*/
      sprintf(filename, "%s/%s", train_dir, Glimmer_posDonModelPath);
      readModel (&donor_pos_model, filename);

      sprintf(filename, "%s/%s", train_dir, Glimmer_negDonModelPath);
      readModel (&donor_neg_model, filename);

      sprintf(filename, "%s/%s", train_dir, Glimmer_posAccModelPath);
      readModel (&acceptor_pos_model, filename);

      sprintf(filename, "%s/%s", train_dir, Glimmer_negAccModelPath);
      readModel (&acceptor_neg_model, filename);

      donor_pos_model_len = getModelLength (donor_pos_model);
      donor_neg_model_len = getModelLength (donor_neg_model);
      acceptor_pos_model_len = getModelLength (acceptor_pos_model);
      acceptor_neg_model_len = getModelLength (acceptor_neg_model);

      if (donor_pos_model_len!=donor_neg_model_len)
         fatal ("ERROR:  Positive and negative donor model lengths differ\n");
      if (acceptor_pos_model_len!=acceptor_neg_model_len)
         fatal ("ERROR:  Positive and negative acceptor model lengths differ\n");

      initGlimmerModel = 1;
   }

   return;

}

double
Sim4::ScoreDonor_Glimmer (char *asegment, char *train_dir)
{
   double  pos_score, neg_score, diff;

#if 0
   /* LLL moved to loadGlimmerModel */
   if (!initGlimmerModel) {
      loadGlimmerModel(train_dir);
      initGlimmerModel = 1;
   }
#endif

   pos_score = Score_Window (donor_pos_model, asegment, GLIMMER_XSPAN);
   neg_score = Score_Window (donor_neg_model, asegment, GLIMMER_XSPAN);
   diff = pos_score - neg_score;

// printf ("%s  %9.5f  %9.5f  %9.5f\n", string, pos_score, neg_score, diff);

   return diff;
}

double
Sim4::ScoreAcceptor_Glimmer (char *asegment, char *train_dir)
{
   double  pos_score, neg_score, diff;

#if 0
   // LLL moved to loadGlimmerModel */
   if (!initGlimmerModel) {
      loadGlimmerModel(train_dir);
      initGlimmerModel = 1;
   }
#endif

   pos_score = Score_Window (acceptor_pos_model, asegment, GLIMMER_ISPAN);
   neg_score = Score_Window (acceptor_neg_model, asegment, GLIMMER_ISPAN);
   diff = pos_score - neg_score;

// printf ("%s  %9.5f  %9.5f  %9.5f\n", string, pos_score, neg_score, diff);

   return diff;
}


/* Generic splice scoring functions: new_splice(), splice_donor(), splice_donor_uni(),
   splice_acceptor(), splice_acceptor_uni(), splice_init()  */


Sim4::splice_t *
Sim4::new_splice(char c, int xs, int xe, int ys, int ye, double score, splice_t *next)
{
  splice_t *sp = (splice_t *)ckalloc(sizeof(splice_t));

  sp->type = c; sp->xs = xs; sp->xe = xe; 
  sp->ys = ys; sp->ye = ye; sp->score = score;
  sp->next = next;

  return sp; 
}

void
Sim4::splice_init(int spl_model)
{
   int i;

   pthread_mutex_lock(&(globalParams->_splice_mutex));

   if (!spliceInit) {

       for (i=0; i<256; spl_encode[i]=0, rev_compl[i]='T', i++);

       spl_encode[(int)'A'] = spl_encode[(int)'a'] = 0;
       spl_encode[(int)'C'] = spl_encode[(int)'c'] = 1;
       spl_encode[(int)'G'] = spl_encode[(int)'g'] = 2;
       spl_encode[(int)'T'] = spl_encode[(int)'t'] = 3;

       rev_compl[(int)'A'] = rev_compl[(int)'a'] = 'T';
       rev_compl[(int)'C'] = rev_compl[(int)'c'] = 'G';
       rev_compl[(int)'G'] = rev_compl[(int)'g'] = 'C';
       rev_compl[(int)'T'] = rev_compl[(int)'t'] = 'A';

       if (spl_model == SPLICE_GENESPLICER) 
          loadGeneSplicerModel();
       else if (spl_model == SPLICE_GLIMMER)
          loadGlimmerModel(Glimmer_TRAIN_DIR);

       spliceInit = 1;
   }
   pthread_mutex_unlock(&(globalParams->_splice_mutex));
}

void
Sim4::splice_donor(char *xseq, char *yseq, int M, int N, double *gt_score,
                   double *ct_score, double **max_Gf, double **max_Cf, 
                   int **start_Gi, int **start_Ci)
{
  int    *CCf, *Xt;
  double *mG, *mC, tmpf;
  int    *sC, *sG; 
  int     i, j, tmpi, ss, ssx, cx, c;
  char   *s, *t;

  CCf = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mG = *max_Gf = (double *)ckalloc((N+1)*sizeof(double));
  sG = *start_Gi = (int *)ckalloc((N+1)*sizeof(int));
  mC = *max_Cf = (double *)ckalloc((N+1)*sizeof(double));
  sC = *start_Ci = (int *)ckalloc((N+1)*sizeof(int));

  t = yseq; Xt[0] = CCf[0] = 0;
  for (j=1; j<=M; j++) { CCf[j] = j; Xt[j] = 0; } 

  mG[0] = mC[0] = -999999;
  for (j=0; j<=M; j++) {
    if ((100*gt_score[j])>mG[0]) { mG[0] = 100*gt_score[j]; sG[0] = j; }
    if ((100*ct_score[j])>mC[0]) { mC[0] = 100*ct_score[j]; sC[0] = j; } 
  }
  
  for (i=1; i<=N; i++, t++) {
    s = xseq;
    ss = CCf[0]; ssx = Xt[0];
    c = ++CCf[0]; cx = Xt[0];
    for (j=1; j<=M; j++, s++) {
      tmpi=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
      if (tmpi==c+1);
      else if (tmpi==CCf[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmpi; ss = CCf[j]; CCf[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_Gf and max_Cf */
    mG[i] = mC[i] = -999999; 
    for (j=0; j<=M; j++) {
      assert(Xt[j]+CCf[j]!=0);
      tmpf = (int)(stepct(j)*Xt[j]/(double)(Xt[j]+CCf[j])*100); 
      if ((tmpf+100*gt_score[j])>mG[i]) {
        mG[i] = tmpf+100*gt_score[j]; sG[i] = j;
#if 0
        fprintf(stderr, "%2d: mG[i]=%1.6f tmpf=%1.6f gt_score[%2d]=%1.6f\n",
                i, mG[i], tmpf, j, gt_score[j]);
#endif
      }
      if ((tmpf+100*ct_score[j])>mC[i]) {
        mC[i] = tmpf+100*ct_score[j]; sC[i] = j;
      }
    } 
  } 
  ckfree(CCf);
  ckfree(Xt); 
}  

void
Sim4::splice_donor_uni(char *xseq, char *yseq, int M, int N,
                       double *It_score, double **max_If, int **start_Ii)
{
  int *CCf, *Xt, tmpi;
  double *mI, tmpf;
  int *sI;
  int i, j, ss, ssx, cx, c;
  char *s, *t;

  CCf = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mI = *max_If = (double *)ckalloc((N+1)*sizeof(double));
  sI = *start_Ii = (int *)ckalloc((N+1)*sizeof(int));

  t = yseq; Xt[0] = CCf[0] = 0;
  for (j=1; j<=M; j++) { CCf[j] = j; Xt[j] = 0; }

  mI[0] = -999999;
  for (j=0; j<=M; j++)
    if ((100*It_score[j])>mI[0]) { mI[0] = 100*It_score[j]; sI[0] = j; }

  for (i=1; i<=N; i++, t++) {
    s = xseq;
    ss = CCf[0]; ssx = Xt[0];
    c = ++CCf[0]; cx = Xt[0];
    for (j=1; j<=M; j++, s++) {
      tmpi=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
      if (tmpi==c+1);
      else if (tmpi==CCf[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmpi; ss = CCf[j]; CCf[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }

    /* compute max_If */
    mI[i] = -999999;
    for (j=0; j<=M; j++) {
      assert(Xt[j]+CCf[j]!=0);
      tmpf = (int)(stepct(j)*Xt[j]/(double)(Xt[j]+CCf[j])*100)+100*It_score[j];
      if (tmpf>mI[i]) {
        mI[i] = tmpf; sI[i] = j;
      }
    }
  }
  ckfree(CCf); ckfree(Xt);
}


void
Sim4::splice_acceptor(char *xseq, char *yseq, int M, int N, 
                      double *ag_score, double *ac_score, double **max_Gb, 
                      double **max_Cb, int **end_Gi, int **end_Ci)
{
  int *CCb, *Xt;
  double *mC, *mG, tmpf;
  int *eC, *eG;
  int tmpi, i, j, ss, ssx, cx, c;
  char *t, *s;

  CCb = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mG = *max_Gb = (double *)ckalloc((N+1)*sizeof(double));
  eG = *end_Gi = (int *)ckalloc((N+1)*sizeof(int));
  mC = *max_Cb = (double *)ckalloc((N+1)*sizeof(double));
  eC = *end_Ci = (int *)ckalloc((N+1)*sizeof(int));

  t = yseq+N-1; CCb[M] = Xt[M] = 0;
  for (j=M-1; j>=0; j--) { CCb[j] = M-j; Xt[j] = 0; }
  
  mG[N] = mC[N] = -999999;
  for (j=M; j>=0; j--) {
    if ((100*ag_score[j])>mG[N]) { mG[N] = 100*ag_score[j]; eG[N] = j+1; }
    if ((100*ac_score[j])>mC[N]) { mC[N] = 100*ac_score[j]; eC[N] = j+1; }
  }
  
  for (i=N-1; i>=0; i--, t--) {
    s = xseq+M-1; 
    ss = CCb[M]; ssx = Xt[M];
    c = ++CCb[M]; cx = Xt[M];
    for (j=M-1; j>=0; j--, s--) {
      tmpi=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
      if (tmpi==c+1) ; 
      else if (tmpi==CCb[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmpi; ss = CCb[j]; CCb[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_Gb and max_Cb */
    mG[i] = -999999; mC[i] = -999999;
    for (j=M; j>=0; j--) {
      assert(CCb[j]+Xt[j]!=0);
      tmpf = (int)(stepct(M-j)*Xt[j]/(double)(CCb[j]+Xt[j])*100);
      if ((tmpf+100*ag_score[j])>mG[i]) {
        mG[i] = tmpf+100*ag_score[j]; eG[i] = j+1;
      }
      if ((tmpf+100*ac_score[j])>mC[i]) {
        mC[i] = tmpf+100*ac_score[j]; eC[i] = j+1;
      }
    } 
  } 
  ckfree(CCb);
  ckfree(Xt); 
}


void
Sim4::splice_acceptor_uni(char *xseq, char *yseq, int M, int N,
                          double *aI_score, double **max_Ib, int **end_Ii)
{
  int *CCb, *Xt;
  double *mI, tmpf;
  int *eI;
  int tmpi, i, j, ss, ssx, cx, c;
  char *t, *s;


  CCb = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mI = *max_Ib = (double *)ckalloc((N+1)*sizeof(double));
  eI = *end_Ii = (int *)ckalloc((N+1)*sizeof(int));

  t = yseq+N-1; CCb[M] = Xt[M] = 0;
  for (j=M-1; j>=0; j--) { CCb[j] = M-j; Xt[j] = 0; }

  mI[N] = -999999;
  for (j=M; j>=0; j--)
    if ((100*aI_score[j])>mI[N]) { mI[N] = 100*aI_score[j]; eI[N] = j+1; }

  for (i=N-1; i>=0; i--, t--) {
    s = xseq+M-1;
    ss = CCb[M]; ssx = Xt[M];
    c = ++CCb[M]; cx = Xt[M];
    for (j=M-1; j>=0; j--, s--) {
      tmpi=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
      if (tmpi==c+1) ;
      else if (tmpi==CCb[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);

      c = tmpi; ss = CCb[j]; CCb[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }

    /* compute max_Ib */
    mI[i] = -999999;
    for (j=M; j>=0; j--) {
      assert(CCb[j]+Xt[j]!=0);
      tmpf = (int)(stepct(M-j)*Xt[j]/(double)(CCb[j]+Xt[j])*100)+100*aI_score[j];
      if (tmpf>mI[i]) {
        mI[i] = tmpf; eI[i] = j+1;
      }
    }
  }
  ckfree(CCb); ckfree(Xt);
}



void
Sim4::splice(char *in_seqx, int ls, int us, int le, int ue,
             char *in_seqy, int ys, int ye, 
             splice_t **gcell, splice_t **ccell, int ori, int spl_model)
{
  double *gtscore=NULL, *ctscore=NULL, *agscore=NULL, *acscore=NULL;
  int i;
  double tmpf, maxCscore, maxGscore;
  int Gxs, Gxe, Gy, Cxs, Cxe, Cy;
  double *max_Cf=NULL, *max_Gf=NULL, *max_Cb=NULL, *max_Gb=NULL;
  int *start_Gi=NULL, *start_Ci=NULL, *end_Gi=NULL, *end_Ci=NULL;
  char *nsegmentL=NULL, *nsegmentR=NULL, *asegmentL=NULL, *asegmentR=NULL;

  //  Initialize the encoding.  This isn't quite as wonderful as
  //  it should be, as there is a chance that two different threads
  //  could initialize the encoding twice, but then again,
  //  it doesn't matter.
  //

//splice_init(spl_model);   LLL

  nsegmentL = (char *) ckalloc(2*MAX_SPAN + 2*MAX_SLIDE + 3);
  nsegmentR = (char *) ckalloc(2*MAX_SPAN + 2*MAX_SLIDE + 3);

  /* Obs: for Glimmer scoring, need only remember the reverse complemented
     segments; but for now we allocate two arrays */
  if (spl_model==SPLICE_GLIMMER) {
     asegmentL = (char *) ckalloc(2*MAX_SPAN + 2*MAX_SLIDE + 3);
     asegmentR = (char *) ckalloc(2*MAX_SPAN + 2*MAX_SLIDE + 3);
  }

  if (ori==FWD || ori==BOTH) {
    gtscore = (double *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(double));
    agscore = gtscore+(us-ls+2);
  }  
  if (ori==BWD || ori==BOTH) {
    ctscore = (double *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(double));
    acscore = ctscore+(us-ls+2);
  }
  
  switch (spl_model) {
    case SPLICE_ORIGINAL:
         splice_original(in_seqx,ls,us,le,ue,in_seqy,ys,ye,gtscore,agscore,ctscore,acscore,ori,nsegmentL,nsegmentR);
         break;

    case SPLICE_GENESPLICER:
         splice_GeneSplicer(in_seqx,ls,us,le,ue,in_seqy,ys,ye,gtscore,agscore,ctscore,acscore,ori,nsegmentL,nsegmentR);
         break;

    case SPLICE_GLIMMER:
         splice_Glimmer(in_seqx,ls,us,le,ue,in_seqy,ys,ye,gtscore,agscore,ctscore,acscore,ori,nsegmentL,nsegmentR,asegmentL,asegmentR);
         break;

    default:
         fprintf(stderr, "Unrecognized splice model (%d). Using original.\n", spl_model);
         splice_original(in_seqx,ls,us,le,ue,in_seqy,ys,ye,gtscore,agscore,ctscore,acscore,ori,nsegmentL,nsegmentR);
         break;
  }


  if (ori==FWD) {
    splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                     gtscore, &max_Gf, &start_Gi);
    splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                        agscore, &max_Gb, &end_Gi);
    ckfree(gtscore);  /* ckfree(agscore) */

  } else if (ori==BWD) {
    splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                     ctscore, &max_Cf, &start_Ci);
    splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                        acscore, &max_Cb, &end_Ci);
    ckfree(ctscore);  /* ckfree(acscore) */

  } else {
    splice_donor(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                 gtscore, ctscore, &max_Gf, &max_Cf, &start_Gi, &start_Ci);
    splice_acceptor(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                    agscore, acscore, &max_Gb, &max_Cb, &end_Gi, &end_Ci);
    ckfree(gtscore);  /* ckfree(agscore); */
    ckfree(ctscore);  /* ckfree(acscore); */
  }

#if 0
  for (i=0; i<=ye-ys+1; i++) {
    fprintf(stderr, "%3d: max_Gf=%8d  max_Cf=%8d  max_Gb=%8d  max_Cb=%8d\n",
            i,
            max_Gf[i], max_Cf[i], max_Gb[i], max_Cb[i]);
  }
#endif

  maxCscore = -999999; maxGscore = -999999;
  Gxs = Gxe = Gy = Cxs = Cxe = Cy = -1;
  if (ori==FWD || ori==BOTH) {
    for (i=0; i<=ye-ys+1; i++) {
      if ((tmpf=max_Gf[i]+max_Gb[i])>maxGscore) {
        maxGscore = tmpf;
        /* save (i, start_Gi[i], end_Gi[i]); */
        Gxs = ls+start_Gi[i]-1; Gxe = le+end_Gi[i]-1; Gy = ys+i-1;
      }
    }
    ckfree(max_Gf); ckfree(max_Gb);
    ckfree(start_Gi); ckfree(end_Gi);
  }
  if (ori==BWD || ori==BOTH) {
    for (i=0; i<=ye-ys+1; i++) {
      if ((tmpf=max_Cf[i]+max_Cb[i])>maxCscore) {
        maxCscore = tmpf;
        /* save (i, start_Ci[i], end_Ci[i]); */
        Cxs = ls+start_Ci[i]-1; Cxe = le+end_Ci[i]-1; Cy = ys+i-1;
      }
    }
    ckfree(max_Cf); ckfree(max_Cb);
    ckfree(start_Ci); ckfree(end_Ci);
  }

#if 0
  fprintf(stderr, "%8d %8d %8d %8d %8f\n%8d %8d %8d %8d %f\n",
          Gxs, Gxe, Gy, Gy+1, maxGscore,
          Cxs, Cxe, Cy, Cy+1, maxCscore);
#endif

  *gcell = new_splice('G', Gxs, Gxe, Gy, Gy+1, maxGscore, NULL);
  *ccell = new_splice('C', Cxs, Cxe, Cy, Cy+1, maxCscore, NULL);

#ifdef DEBUG
  printf("Type: %c  sx: %d  se: %d  ys: %d  score: %d\n",
         gcell.type, gcell.xs, gcell.xe, gcell.ys, gcell.score);

  printf("Type: %c  sx: %d  se: %d  ys: %d  score: %d\n",
         ccell.type, ccell.xs, ccell.xe, ccell.ys, ccell.score);
#endif

  ckfree(nsegmentL); ckfree(nsegmentR);

  if (spl_model==SPLICE_GLIMMER) {
     ckfree(asegmentL); ckfree(asegmentR);
  }

  return;
}


/* Customized splice signal scoring functions:
   splice_original(), splice_GeneSplicer(), splice_Glimmer() */

void
Sim4::splice_original(char *in_seqx, int ls, int us, int le, int ue,
                      char *in_seqy, int ys, int ye,
                      double *gtscore, double *agscore,
                      double *ctscore, double *acscore, int ori,
                      char *nsegmentL, char *nsegmentR)
{
   int p, q, i;
   char  *s,*t, ch;

  
   for (i=0, s=in_seqx+ls-MAX_SPAN-1; i<2*MAX_SPAN+us-ls+3; nsegmentL[i++] = spl_encode[(int)(*s++)]);
   for (i=0, s=in_seqx+le-2-MAX_SPAN-1; i<2*MAX_SPAN+ue-le+3; nsegmentR[i++] = spl_encode[(int)(*s++)]);
                                                                                                                         
   if (ori==FWD || ori==BOTH) { 
  
       if (globalParams->_dontForceCanonicalSplicing) {
          for (p=0, s=nsegmentL+MAX_SPAN; p<=us-ls+1; p++, s++)
            gtscore[p] = 0; 
          for (q=ue-le+1, s=nsegmentR+MAX_SPAN+ue-le+2; q>=0; q--, s--)
            agscore[q] = 0;
       } else { 
          for (p=0, s=nsegmentL+MAX_SPAN; p<=us-ls+1; p++, s++)
            gtscore[p] = gt[(int)*s][(int)*(s+1)];
          for (q=ue-le+1, s=nsegmentR+MAX_SPAN+ue-le+2; q>=0; q--, s--)
            agscore[q] = ag[(int)*(s-1)][(int)*s];
       }
   }

         
   if (ori==BWD || ori==BOTH) {

       /* reverse complement the nsegments, 0-3 alphabet */
       for (s=nsegmentL, t=nsegmentL+2*MAX_SPAN+us-ls+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }
       for (s=nsegmentR, t=nsegmentR+2*MAX_SPAN+ue-le+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }

       if (globalParams->_dontForceCanonicalSplicing) {
         for (p=0, s=nsegmentL+MAX_SPAN+us-ls+2; p<=us-ls+1; p++, s++)
           ctscore[p] = 0;
         for (q=ue-le+1, s=nsegmentR+MAX_SPAN; q>=0; q--, s--)
           acscore[q] = 0;
       } else {
         for (p=0, s=nsegmentL+MAX_SPAN+us-ls+2; p<=us-ls+1; p++, s--)
           ctscore[p] = ag[(int)*(s-1)][(int)*s];
         for (q=ue-le+1, s=nsegmentR+MAX_SPAN; q>=0; q--, s++)
           acscore[q] = gt[(int)*s][(int)*(s+1)];
       }
   }

   return;
}


void
Sim4::splice_GeneSplicer(char *in_seqx, int ls, int us, int le, int ue,
                         char *in_seqy, int ys, int ye,
                         double *gtscore, double *agscore,
                         double *ctscore, double *acscore, int ori,
                         char *nsegmentL, char *nsegmentR)
{
   int p, q, i;
   char  *s,*t, ch;


   for (i=0, s=in_seqx+ls-MAX_SPAN-1; i<2*MAX_SPAN+us-ls+3; nsegmentL[i++] = spl_encode[(int)(*s++)]);
   for (i=0, s=in_seqx+le-2-MAX_SPAN-1; i<2*MAX_SPAN+ue-le+3; nsegmentR[i++] = spl_encode[(int)(*s++)]);


   if (ori==FWD || ori==BOTH) {

       for (p=0, s=nsegmentL+MAX_SPAN; p<=us-ls+1; p++, s++) {
           gtscore[p] = ScoreDonor_GeneSplicer(s-GENESPLICER_SPAN);

           if (gtscore[p] < -14) gtscore[p] = -14.0;
           if (gtscore[p] > 19) gtscore[p] = 19;
           gtscore[p] = 5.0*(gtscore[p]+14.0)/33.0;
           gtscore[p] = 0.4*gtscore[p] + 0.6*gt[(int)*s][(int)*(s+1)];
       }
       for (q=ue-le+1, s=nsegmentR+MAX_SPAN+ue-le+2; q>=0; q--, s--) {
           agscore[q] = ScoreAcceptor_GeneSplicer(s-GENESPLICER_SPAN-1);

           if (agscore[q] < -23) agscore[q] = -23.0;
           if (agscore[q] > 20) agscore[q] = 20.0;
           agscore[q] = 5.0*(agscore[q]+23.0)/43.0;
           agscore[q] = 0.4*agscore[q] + 0.6*ag[(int)*(s-1)][(int)*s];
       }
   }


   if (ori==BWD || ori==BOTH) {

       /* reverse complement the nsegments, 0-3 alphabet */
       for (s=nsegmentL, t=nsegmentL+2*MAX_SPAN+us-ls+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }
       for (s=nsegmentR, t=nsegmentR+2*MAX_SPAN+ue-le+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }


       for (p=0, s=nsegmentL+MAX_SPAN+us-ls+2; p<=us-ls+1; p++, s--) {
           ctscore[p] = ScoreAcceptor_GeneSplicer(s-GENESPLICER_SPAN-1);


           if (ctscore[p] < -23) ctscore[p] = -23.0;
           if (ctscore[p] > 20) ctscore[p] = 20.0;
           ctscore[p] = 5.0*(ctscore[p]+23.0)/43.0;
           ctscore[p] = 0.4*ctscore[p] + 0.6*ag[(int)*(s-1)][(int)*s];
       }
       for (q=ue-le+1, s=nsegmentR+MAX_SPAN; q>=0; q--, s++) {
           acscore[q] = ScoreDonor_GeneSplicer(s-GENESPLICER_SPAN);


           if (acscore[q] < -14) acscore[q] = -14.0;
           if (acscore[q] > 19) acscore[q] = 19.0;
           acscore[q] = 5.0*(acscore[q]+14.0)/33.0;
           acscore[q] = 0.4*acscore[q] + 0.6*gt[(int)*s][(int)*(s+1)];
       }
   }

   return;
}


void
Sim4::splice_Glimmer(char *in_seqx, int ls, int us, int le, int ue,
                     char *in_seqy, int ys, int ye,
                     double *gtscore, double *agscore,
                     double *ctscore, double *acscore, int ori,
                     char *nsegmentL, char *nsegmentR, char *asegmentL, char *asegmentR)
{
   int p, q, i;
   char  *s,*t, ch;


   for (i=0, s=in_seqx+ls-MAX_SPAN-1; i<2*MAX_SPAN+us-ls+3; nsegmentL[i++] = spl_encode[(int)(*s++)]);
   for (i=0, s=in_seqx+le-2-MAX_SPAN-1; i<2*MAX_SPAN+ue-le+3; nsegmentR[i++] = spl_encode[(int)(*s++)]);


   /* Glimmer specific matrices */
   for (i=0, s=in_seqx+ls-MAX_SPAN-1; i<2*MAX_SPAN+us-ls+3; asegmentL[i++] = *s++);
   for (i=0, s=in_seqx+le-2-MAX_SPAN-1; i<2*MAX_SPAN+ue-le+3; asegmentR[i++] = *s++);

   if (ori==FWD || ori==BOTH) {

      for (p=0, s=nsegmentL+MAX_SPAN, t=asegmentL+MAX_SPAN; p<=us-ls+1; p++, s++, t++) {
          gtscore[p] = ScoreDonor_Glimmer(t-GLIMMER_XSPAN, Glimmer_TRAIN_DIR);

          if (gtscore[p] < 0) gtscore[p] = 0.0;
          if (gtscore[p] > 0.31) gtscore[p] = 0.31;
          gtscore[p] = 5.0*(gtscore[p]+0.0)/0.31;
          gtscore[p] = 0.2*gtscore[p] + 0.8*gt[(int)*s][(int)*(s+1)];
      }
      for (q=ue-le+1, s=nsegmentR+MAX_SPAN+ue-le+2, t=asegmentR+MAX_SPAN+ue-le+2; q>=0; q--, s--, t--) {
          agscore[q] = ScoreAcceptor_Glimmer(t-GLIMMER_ISPAN-1, Glimmer_TRAIN_DIR);


          if (agscore[q] < -0.16) agscore[q] = -0.16;
          if (agscore[q] > 0.23) agscore[q] = 0.23;
          agscore[q] = 5.0*(agscore[q]+0.16)/0.39;
          agscore[q] = 0.2*agscore[q] + 0.8*ag[(int)*(s-1)][(int)*s];
      }
   }


   if (ori==BWD || ori==BOTH) {

       /* reverse complement the nsegments, 0-3 alphabet */
       for (s=nsegmentL, t=nsegmentL+2*MAX_SPAN+us-ls+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }
       for (s=nsegmentR, t=nsegmentR+2*MAX_SPAN+ue-le+3-1; s<t; s++, t--)
           { ch = 3-(*s); *s = 3-(*t); *t = ch; }




       /* reverse complement the asegments, ACTG alphabet */
       for (s=asegmentL, t=asegmentL+2*MAX_SPAN+us-ls+3-1; s<t; s++, t--)
            { ch = rev_compl[(int)*s]; *s = rev_compl[(int)*t]; *t = ch; }
       for (s=asegmentR, t=asegmentR+2*MAX_SPAN+ue-le+3-1; s<t; s++, t--)
            { ch = rev_compl[(int)*s]; *s = rev_compl[(int)*t]; *t = ch; }
   

       for (p=0, s=nsegmentL+MAX_SPAN+us-ls+2, t=asegmentL+MAX_SPAN+us-ls+2; p<=us-ls+1; p++, s--, t--)  {
           ctscore[p] = ScoreAcceptor_Glimmer(t-GLIMMER_ISPAN-1, Glimmer_TRAIN_DIR);


           if (ctscore[p] < -0.16) ctscore[p] = -0.16;
           if (ctscore[p] > 0.23) ctscore[p] = 0.23;
           ctscore[p] = 5.0*(ctscore[p]+0.16)/0.39;
           ctscore[p] = 0.2*ctscore[p] + 0.8*ag[(int)*(s-1)][(int)*s];
       }
       for (q=ue-le+1, s=nsegmentR+MAX_SPAN, t=asegmentR+MAX_SPAN; q>=0; q--, s++, t++) {
           acscore[q] = ScoreDonor_Glimmer(t-GLIMMER_XSPAN, Glimmer_TRAIN_DIR);
   

           if (acscore[q] < 0) acscore[q] = 0.0;
           if (acscore[q] > 0.31) acscore[q] = 0.31;
           acscore[q] = 5.0*(acscore[q]+0.0)/0.31;
           acscore[q] = 0.2*acscore[q] + 0.8*gt[(int)*s][(int)*(s+1)];
       }
   }
   
   return;
}

void
Sim4::splice_close ()
{
   UnLoadSites_GeneSplicer();

   spliceInit = 0;
}

