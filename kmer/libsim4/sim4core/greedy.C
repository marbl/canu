#include "sim4.H"

//#define ANNOUNCEEXIT(S) fprintf(stdout, S);
#define ANNOUNCEEXIT(S)

int
Sim4::greedy(uchar *s1, uchar *s2, int m, int n0, int OFFSET1, int OFFSET2, Exon **lblock, Exon **rblock)
{
  int     col,                    /* column number */
  d,                      /* current distance */
  k,                      /* current diagonal */
  Cost,
  blower,flower,          /* boundaries for searching diagonals */
  bupper,fupper,
  row;                    /* row number */
  int     max_d;                  /* bound on size of edit script */
  int     back, forth;            /* backward and forward limits at exit */

  int     *blast_d, *flast_d;     /* rows containing the last d (at crt step, d-1) */
  int     *btemp_d, *ftemp_d;     /* rows containing tmp values for the last d */
  int     *min_row, *min_diag;    /* min (b)/ max (f) row (and diagonal) */
  int     *max_row, *max_diag;    /* reached for cost d=0, ... m.  */
  
  const int MAX_D = max_d = max(wordSize,(int)(P*m+1));

  if (n0 < m) {
    if (m < (int)min(wordSize, (1+P)*n0)) {
      *lblock = *rblock = new_exon(OFFSET2+1,OFFSET1+1,OFFSET2+n0,OFFSET1+m,
                                   m,n0-m+(int)(P*m+1),0,NULL);
      ANNOUNCEEXIT("greedy-1\n");
      return(m-n0+(int)(P*n0+1));
    } else if (m > (int)min(wordSize, (1+P)*n0)) {
      *lblock = *rblock = 0L;
      ANNOUNCEEXIT("greedy-2\n");
      return(MAX_D+1);
    }
  }       

  const int n1    = min(m+max_d+1, n0);
  const int n2    = n1;
  const int DELTA = n2-m;

  const int l_offset1 = OFFSET1;
  const int r_offset1 = OFFSET1;
  const int l_offset2 = OFFSET2;
  const int r_offset2 = OFFSET2 + n0 - n2;

  const int L_ORIGIN = MAX_D;
  const int R_ORIGIN = MAX_D - DELTA;

  const unsigned char *l_s1 = s1;
  const unsigned char *r_s1 = s1;
  const unsigned char *l_s2 = s2;
  const unsigned char *r_s2 = s2 + n0 - n2;


  for (row=m, col=n2; row>0 && col>0 && (r_s1[row-1]==r_s2[col-1]); row--,col--)
    /*LINTED empty loop body*/;
  
  if (row == 0) {
    /* hit last row; stop search */
    *lblock = *rblock = new_exon(r_offset2-m+n2+1,r_offset1+1,r_offset2+n2,
                                 r_offset1+m,m,0,0,NULL);
    ANNOUNCEEXIT("greedy-3\n");
    return 0;
  }
  

  //  Instead of doing eight calls to ckalloc, we do one, and dish out
  //  that in pieces.
  //

  int *allocdSpace = (int *)ckalloc((4*(MAX_D+n2+1) + 4*(MAX_D+1)) * sizeof(int));

  blast_d  = allocdSpace;              //  MAX_D+n2+1
  btemp_d  = blast_d  + (MAX_D+n2+1);  //  MAX_D+n2+1
  flast_d  = btemp_d  + (MAX_D+n2+1);  //  MAX_D+n2+1
  ftemp_d  = flast_d  + (MAX_D+n2+1);  //  MAX_D+n2+1
  max_row  = ftemp_d  + (MAX_D+n2+1);  //  MAX_D+1
  min_row  = max_row  + (MAX_D+1);     //  MAX_D+1
  max_diag = min_row  + (MAX_D+1);     //  MAX_D+1
  min_diag = max_diag + (MAX_D+1);     //  MAX_D+1


  for (k=0; k<=MAX_D+n2; ++k) {
    blast_d[k] = m+1;
    btemp_d[k] = m+1;
  }

  blast_d[R_ORIGIN+DELTA] = row;
  
  blower = R_ORIGIN + DELTA - 1;
  bupper = R_ORIGIN + DELTA + 1;


  for (row=0; row<n1 && row<m && (l_s1[row]==l_s2[row]); row++)
    /*LINTED empty loop body*/;

  if (row == m) {
    /* hit last row; stop search */
    *lblock = *rblock = new_exon(l_offset2+1,l_offset1+1,l_offset2+m,
                                 l_offset1+m,m,0,0,NULL);
    free(allocdSpace);

    ANNOUNCEEXIT("greedy-4\n");
    return 0;
  }

  for (k=0; k<=MAX_D+n1; ++k) {
    flast_d[k]=-1;
    ftemp_d[k]=-1;
  }
  flast_d[L_ORIGIN] = row;           
  
  flower = L_ORIGIN - 1;
  fupper = L_ORIGIN + 1;
  
  for (d=1; d<=MAX_D; d++) {
    min_row[d] = m+1;
    max_row[d] = -1;
  }
  min_row[0]  = blast_d[R_ORIGIN+DELTA];
  min_diag[0] = R_ORIGIN+DELTA;
  max_row[0]  = flast_d[L_ORIGIN];
  max_diag[0] = L_ORIGIN;
  
  back = forth = -1;
  
  d = 1;
  while (d <= max_d) {             
    
    /* for each relevant diagonal ... */
    for (k = blower; k <= bupper; k++) {
      /* get space for the next edit instruction */

      /* find a d on diagonal k */
      if (k==-d+DELTA+R_ORIGIN) {
        
        /* move left from the last d-1 on diagonal k+1 */
        row = blast_d[k+1];
      }
      else if (k==d+DELTA+R_ORIGIN) {
        
        /* move up from the last d-1 on diagonal k-1 */
        row = blast_d[k-1]-1;
      } else if ((blast_d[k]<=blast_d[k+1]) &&
                 (blast_d[k]-1<=blast_d[k-1])) {
        
        /* substitution */
        row = blast_d[k]-1;
        
      } else if ((blast_d[k-1]<=blast_d[k+1]-1) &&
                 (blast_d[k-1]<=blast_d[k]-1)) {
        /* move right from the last d-1 on diagonal k-1 */
        row = blast_d[k-1]-1;
      } else  {
        /* move left from the last d-1 on diagonal k+1 */
        row = blast_d[k+1];
      }
      /* code common to the three cases */
      col = row + k - R_ORIGIN;
      
      /* slide up the diagonal */
      while (row > 0 && col > 0 && (r_s1[row-1]==r_s2[col-1])) {
        --row;
        --col;
      }        
      btemp_d[k] = row;

#if 0
      if (row == 0 || col == 0)
        max_d = d;
#endif
    }     /* for k */
    
    min_row[d] = btemp_d[DELTA+R_ORIGIN];
    min_diag[d] = DELTA+R_ORIGIN;
    for (k=blower; k<=bupper; ++k) {
      blast_d[k] = btemp_d[k]; btemp_d[k] = m+1;
      if (blast_d[k]<min_row[d]) {  
        min_row[d] = blast_d[k];
        min_diag[d] = k;
      }
    }                                         
    
    /* record cell, if paths overlap with minimum combined cost */
    /* obs: it suffices to search up to Cost=min(d-1,(max_d-d)) */
    for (Cost=0; Cost<d; Cost++) {
      if ((min_row[d]<=max_row[Cost]) &&
          ((max_d > d+Cost) || (max_d==d+Cost && (forth<0)))) {
        max_d = d+Cost;
        back  = d;
        forth = Cost;
        break;
      }
    }
    
    --blower; ++bupper;
    
    /* for each relevant diagonal ... */
    for (k = flower; k <= fupper; k++) {
      /* get space for the next edit instruction */
      
      /* find a d on diagonal k */       
      if (k==-d+L_ORIGIN) {
        /* move down from the last d-1 on diagonal k+1 */
        row = flast_d[k+1]+1;
        
      } else if (k==d+L_ORIGIN) {
        /* move right from the last d-1 on diagonal k-1 */
        row = flast_d[k-1];
        
      } else if ((flast_d[k]>=flast_d[k+1]) &&
                 (flast_d[k]+1>=flast_d[k-1])) {
        
        /* substitution */
        row = flast_d[k]+1; 
        
      } else if ((flast_d[k+1]+1>=flast_d[k-1]) &&
                 (flast_d[k+1]>=flast_d[k])) {
        
        /* move left from the last d-1 on diagonal k+1 */
        row = flast_d[k+1]+1;
      } else { 
        /* move right from the last d-1 on diagonal k-1 */
        row = flast_d[k-1];
      } 
      /* code common to the three cases */
      col = row + k - L_ORIGIN;
      /* slide down the diagonal */
      if (row>=0)
        while (row < m && col < n1 && (l_s1[row]==l_s2[col])) {
          ++row;
          ++col;
        }        
      ftemp_d[k] = row;

#if 0
      if (row == m || col == n1)
        max_d = d;
#endif
    }     /* for k */    
    
    max_row[d] = ftemp_d[L_ORIGIN];
    max_diag[d] = L_ORIGIN;
    for (k=flower; k<=fupper; ++k) {
      flast_d[k] = ftemp_d[k]; ftemp_d[k] = -1;
      if (flast_d[k]>max_row[d]) {
        max_row[d] = flast_d[k];  
        max_diag[d] = k;
      }       
    }            
    
    /* record backward and forward limits, if minimum combined
     * cost in overlapping. Note: it suffices to search up to
     * Cost=min(d,(max_d-d)).
     */          
    for (Cost=0; Cost<=d; Cost++) {
      if ((min_row[Cost]<=max_row[d]) &&
          ((max_d>d+Cost) || (max_d==d+Cost && (forth<0)))) {
        max_d = d+Cost;
        back  = Cost;
        forth = d;
        break;
      }
    }
    --flower;
    ++fupper;
    
    ++d;  /* for d */
  }
  
  if (d>MAX_D) {
    *lblock = *rblock = NULL;
    free(allocdSpace);
    ANNOUNCEEXIT("greedy-5\n");
    return d;  
  }


  //  XXX:  Quick fix!
  //
  if ((back < 0) || (forth < 0)) {
    *rblock = *lblock = 0L;
    fprintf(stdout, "Choke!\n");
    return(MAX_D+1);
  }


  if (m-min_row[back]>=max_row[forth]) {

      if ((r_offset2+1+min_diag[back]-R_ORIGIN) < 
          (l_offset2+max_diag[forth]-L_ORIGIN)) {
        *rblock = *lblock = new_exon(l_offset2+1,l_offset1+1,
                                     l_offset2+n0,l_offset1+m,
                                     m,back+forth,0,NULL);
      } else {
        *rblock = new_exon(r_offset2+1+min_row[back]+min_diag[back]-R_ORIGIN,
                           r_offset1+1+min_row[back],
                           r_offset2+n2,r_offset1+m,
                           m-min_row[back],back,0,NULL);
        *lblock = new_exon(l_offset2+1,l_offset1+1,
                           l_offset2+min_row[back]+max_diag[forth]-L_ORIGIN,
                           l_offset1+min_row[back],
                           min_row[back],forth,0,*rblock);
      }
  } else {
      if ((r_offset2+1+min_diag[back]-R_ORIGIN) < 
          (l_offset2+max_diag[forth]-L_ORIGIN)) {
        *rblock = *lblock = new_exon(l_offset2+1,l_offset1+1,
                                       l_offset2+n0,l_offset1+m,
                                       m,back+forth,0,NULL);
      } else {
        *rblock = new_exon(r_offset2+1+max_row[forth]+min_diag[back]-R_ORIGIN,
                           r_offset1+1+max_row[forth],
                           r_offset2+n2,r_offset1+m,m-max_row[forth],back,0,NULL);
        *lblock = new_exon(l_offset2+1,l_offset1+1,
                           l_offset2+max_row[forth]+max_diag[forth]-L_ORIGIN,
                           l_offset1+max_row[forth],max_row[forth],forth,0,*rblock);
      }
  }
  
  free(allocdSpace);

  ANNOUNCEEXIT("greedy-6\n");
  return back+forth;
}
