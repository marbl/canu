#include "sim4.H"
#include "sim4db.H"

#ifndef __lint
/*@unused@*/
static const char rcsid[] =
"$Id$";
#endif

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


int encodeInitialized = 0;
int encode[256]       = { 0 };

#if 0
int const gt[4][4] = {{0, 0, 0, 2},{0, 0, 0, 2},{2, 2, 2, 5},{0, 0, 0, 2}};
int const ct[4][4] = {{0, 0, 0, 2},{2, 2, 2, 5},{0, 0, 0, 2},{0, 0, 0, 2}};
int const ag[4][4] = {{2, 2, 5, 2},{0, 0, 2, 0},{0, 0, 2, 0},{0, 0, 2, 0}};
int const ac[4][4] = {{2, 5, 2, 2},{0, 2, 0, 0},{0, 2, 0, 0},{0, 2, 0, 0}};
#endif



splice_t *
Sim4::new_splice(char c, int xs, int xe, int ys, int ye, int score, splice_t *next)
{
  splice_t *sp = (splice_t *)ckalloc(sizeof(splice_t));

  sp->type = c; sp->xs = xs; sp->xe = xe; 
  sp->ys = ys; sp->ye = ye; sp->score = score;
  sp->next = next;

  return sp; 
}

void
Sim4::splice_donor(uchar *xseq, uchar *yseq, int M, int N, int *gt_score,
                   int *ct_score, int **max_Gf, int **max_Cf, 
                   int **start_Gi, int **start_Ci)
{
  int *CCf, *mG, *mC, *sC, *sG, *Xt;
  int i, j, tmp, ss, ssx, cx, c;
  uchar *s, *t;

  CCf = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mG = *max_Gf = (int *)ckalloc((2*N+2)*sizeof(int));
  sG = *start_Gi = mG+(N+1);
  mC = *max_Cf = (int *)ckalloc((2*N+2)*sizeof(int));
  sC = *start_Ci = mC+(N+1);

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
      tmp=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
      if (tmp==c+1);
      else if (tmp==CCf[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmp; ss = CCf[j]; CCf[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_Gf and max_Cf */
    mG[i] = mC[i] = -999999; 
    for (j=0; j<=M; j++) {
      assert(Xt[j]+CCf[j]!=0);
      tmp = (int)(stepct(j)*Xt[j]/(double)(Xt[j]+CCf[j])*100); 
      if ((tmp+100*gt_score[j])>mG[i]) {
        mG[i] = tmp+100*gt_score[j]; sG[i] = j;
#if 0
        fprintf(stderr, "%2d: mG[i]=%8d tmp=%8d gt_score[%2d]=%8d\n",
                i, mG[i], tmp, j, gt_score[j]);
#endif
      }
      if ((tmp+100*ct_score[j])>mC[i]) {
        mC[i] = tmp+100*ct_score[j]; sC[i] = j;
      }
    } 
  } 
  free(CCf); free(Xt); 
}  

void
Sim4::splice_donor_uni(uchar *xseq, uchar *yseq, int M, int N,
                       int *It_score, int **max_If, int **start_Ii)
{
  int *CCf, *mI, *sI, *Xt;
  int i, j, tmp, ss, ssx, cx, c;
  uchar *s, *t;

  CCf = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mI = *max_If = (int *)ckalloc((2*N+2)*sizeof(int));
  sI = *start_Ii = mI+(N+1);
  
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
      tmp=min(min(CCf[j]+1, ss+(*t!=*s)),c+1);
      if (tmp==c+1);
      else if (tmp==CCf[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmp; ss = CCf[j]; CCf[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_If */           
    mI[i] = -999999;        
    for (j=0; j<=M; j++) {
      assert(Xt[j]+CCf[j]!=0);
      tmp = (int)(stepct(j)*Xt[j]/(double)(Xt[j]+CCf[j])*100)+100*It_score[j];
      if (tmp>mI[i]) {
        mI[i] = tmp; sI[i] = j;
      }
    }
  }
  free(CCf); free(Xt);
}


void
Sim4::splice_acceptor(uchar *xseq, uchar *yseq, int M, int N, 
                      int *ag_score, int *ac_score, int **max_Gb, 
                      int **max_Cb, int **end_Gi, int **end_Ci)
{
  int *CCb, *Xt, *mC, *mG, *eC, *eG;
  int tmp, i, j, ss, ssx, cx, c;
  uchar *t, *s;

  CCb = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mG = *max_Gb = (int *)ckalloc((2*N+2)*sizeof(int));
  eG = *end_Gi = mG+(N+1);
  mC = *max_Cb = (int *)ckalloc((2*N+2)*sizeof(int));
  eC = *end_Ci = mC+(N+1);
  
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
      tmp=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
      if (tmp==c+1) ; 
      else if (tmp==CCb[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      c = tmp; ss = CCb[j]; CCb[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_Gb and max_Cb */
    mG[i] = -999999; mC[i] = -999999;
    for (j=M; j>=0; j--) {
      assert(CCb[j]+Xt[j]!=0);
      tmp = (int)(stepct(M-j)*Xt[j]/(double)(CCb[j]+Xt[j])*100);
      if ((tmp+100*ag_score[j])>mG[i]) {
        mG[i] = tmp+100*ag_score[j]; eG[i] = j+1;
      }
      if ((tmp+100*ac_score[j])>mC[i]) {
        mC[i] = tmp+100*ac_score[j]; eC[i] = j+1;
      }
    } 
  } 
  free(CCb); free(Xt); 
}


void
Sim4::splice_acceptor_uni(uchar *xseq, uchar *yseq, int M, int N,
                          int *aI_score, int **max_Ib, int **end_Ii)
{
  int *CCb, *Xt, *mI, *eI;
  int tmp, i, j, ss, ssx, cx, c;
  uchar *t, *s;

  CCb = (int *)ckalloc((M+1)*sizeof(int));
  Xt = (int *)ckalloc((M+1)*sizeof(int));
  mI = *max_Ib = (int *)ckalloc((2*N+2)*sizeof(int));
  eI = *end_Ii = mI+(N+1);

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
      tmp=min(min(CCb[j]+1, ss+(*t!=*s)),c+1);
      if (tmp==c+1) ;
      else if (tmp==CCb[j]+1) cx = Xt[j];
      else cx = ssx + (*t==*s);
      
      c = tmp; ss = CCb[j]; CCb[j] = c; ssx = Xt[j]; Xt[j] = cx;
    }
    
    /* compute max_Ib */           
    mI[i] = -999999; 
    for (j=M; j>=0; j--) {
      assert(CCb[j]+Xt[j]!=0);
      tmp = (int)(stepct(M-j)*Xt[j]/(double)(CCb[j]+Xt[j])*100)+100*aI_score[j];
      if (tmp>mI[i]) {
        mI[i] = tmp; eI[i] = j+1;
      }
    }
  }
  free(CCb); free(Xt);
}



void
Sim4::splice(uchar *in_seqx, int ls, int us, int le, int ue,
             uchar *in_seqy, int ys, int ye, 
             splice_t **gcell, splice_t **ccell, int ori)
{
  int p, q, *gtscore=NULL, *ctscore=NULL, *agscore=NULL, *acscore=NULL;
  int i, tmp;
  int maxCscore, maxGscore, Gxs, Gxe, Gy, Cxs, Cxe, Cy, keep_Ci, keep_Gi;
  int *max_Cf=NULL, *max_Gf=NULL, *max_Cb=NULL, *max_Gb=NULL;
  int *start_Gi=NULL, *start_Ci=NULL, *end_Gi=NULL, *end_Ci=NULL;
  uchar *s;

#if 0
  fprintf(stderr, "Hello from splice()!\n");
  fprintf(stderr, "ls=%8d us=%8d le=%8d ue=%8d ys=%8d ye=%8d\n", ls, us, le, ue, ys, ye);
  fprintf(stderr, "seqX = %80.80s\n", in_seqx + ls - 1);
  fprintf(stderr, "seqY = %80.80s\n", in_seqy);
#endif

  //  Initialize the encoding.  This isn't quite as wonderful as
  //  it should be, as there is a chance that two different threads
  //  could initialize the encoding twice, but then again,
  //  it doesn't matter.
  //
  if (encodeInitialized == 0) {
    encodeInitialized = 1;

    for (unsigned int i=256; i;)
      encode[--i] = 4;

    encode['A'] = encode['a'] = 0;
    encode['C'] = encode['c'] = 1;
    encode['G'] = encode['g'] = 2;
    encode['T'] = encode['t'] = 3;
  }

  if (ori==FWD || ori==BOTH) {
    gtscore = (int *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(int));
    agscore = gtscore+(us-ls+2);

    if (dbParams._dontForceCanonicalSplicing) {
      for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++) 
        gtscore[p] = 0;
      for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
        agscore[q] = 0;
    } else {
      for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++) 
        gtscore[p] = gt[encode[*s]][encode[*(s+1)]];
      for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
        agscore[q] = ag[encode[*(s-1)]][encode[*s]]; 
    }
#if 0
    for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++) 
      fprintf(stderr, "gtscore[%2d] = %d\n", p, gtscore[p]);
#endif
  } 
  if (ori==BWD || ori==BOTH) {
    ctscore = (int *)ckalloc(((us-ls+2)+(ue-le+2))*sizeof(int));
    acscore = ctscore+(us-ls+2);
    if (dbParams._dontForceCanonicalSplicing) {
      for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++)  
        ctscore[p] = 0;
      for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
        acscore[q] = 0;
    } else {
      for (p=0, s=in_seqx+ls-1; p<=us-ls+1; p++, s++)  
        ctscore[p] = ct[encode[*s]][encode[*(s+1)]];  
      for (q=ue-le+1, s=in_seqx+ue-1; q>=0; q--, s--)
        acscore[q] = ac[encode[*(s-1)]][encode[*s]];
    }
  } 

  if (ori==FWD) {
    splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                     gtscore, &max_Gf, &start_Gi);
    splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                        agscore, &max_Gb, &end_Gi);
    free(gtscore);              /* free(agscore); */

  } else if (ori==BWD) {
    splice_donor_uni(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1,
                     ctscore, &max_Cf, &start_Ci);              
    splice_acceptor_uni(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1,
                        acscore, &max_Cb, &end_Ci);
    free(ctscore);              /* free(acscore); */

  } else {
    splice_donor(in_seqx+ls-1, in_seqy+ys-1, us-ls+1, ye-ys+1, 
                 gtscore, ctscore, &max_Gf, &max_Cf, &start_Gi, &start_Ci);

    splice_acceptor(in_seqx+le-1, in_seqy+ys-1, ue-le+1, ye-ys+1, 
                    agscore, acscore, &max_Gb, &max_Cb, &end_Gi, &end_Ci);

    free(gtscore);              /* free(agscore);     */
    free(ctscore);              /* free(acscore); */
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
      if ((tmp=max_Gf[i]+max_Gb[i])>maxGscore) {
        maxGscore = tmp;
        /* save (i, start_Gi[i], end_Gi[i]); */
        Gxs = ls+start_Gi[i]-1; Gxe = le+end_Gi[i]-1; Gy = ys+i-1;
        keep_Gi = i;
      }  
    }
    free(max_Gf); free(max_Gb); /* free(start_Gi); free(end_Gi); */
  }
  if (ori==BWD || ori==BOTH) {
    for (i=0; i<=ye-ys+1; i++) {
      if ((tmp=max_Cf[i]+max_Cb[i])>maxCscore) {
        maxCscore = tmp;
        /* save (i, start_Ci[i], end_Ci[i]); */
        Cxs = ls+start_Ci[i]-1; Cxe = le+end_Ci[i]-1; Cy = ys+i-1;
        keep_Ci = i;
      }
    }
    free(max_Cf); free(max_Cb); /* free(start_Ci); free(end_Ci); */
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

  return;
}
