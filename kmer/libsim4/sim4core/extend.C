#include "sim4.H"




int
Sim4::extend_bw(uchar *s1, uchar *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)
{
  int     col,                    /* column number */
  row,                    /* row number */
  max_d,                  /* bound on the length of the edit script
                           */
  d,                      /* current compressed distance */
  k,                      /* current diagonal */
  DELTA,                  /* n-m  */
  ORIGIN,
  lower,
  upper;
  int     *last_d, *temp_d;       /* column containing the last p */
  int     *min_row, *min_diag;    /* min (b)/ max (f) row (and diagonal) */
  /* reached for cost d=0, ... m.  */
  DELTA = n-m;
  max_d = m+1;

  ORIGIN = m;
  for (row=m, col=n; row>0 && col>0 && (s1[row-1]==s2[col-1]); row--,col--)
    /*LINTED empty loop body*/; 

  if ((row == 0) || (col == 0)) {
    *line1 = row+offset1;
    *line2 = col+offset2;

    return 0;
  }

  int *allocdSpace = (int *)ckalloc((m+n+1+m+n+1+m+1+m+1) * sizeof(int));

  last_d   = allocdSpace;     //  m+n+1
  temp_d   = last_d + m+n+1;  //  m+n+1
  min_row  = temp_d + m+n+1;  //  m+1
  min_diag = min_row + m+1;   //  m+1

  for (k=0; k<=m+n; ++k)
    last_d[k]=m+1;
  last_d[ORIGIN+DELTA] = row;

  lower = ORIGIN + DELTA - 1;
  upper = ORIGIN + DELTA + 1;

  for (d=1; d<=m; d++)
    min_row[d] = m+1;
  
  min_row[0] = last_d[ORIGIN+DELTA];
  min_diag[0] = ORIGIN + DELTA;

  d = 0;
  while ((++d<=max_d) && 
         ((d-1<=good_ratio(m-min_row[d-1])) ||
          ((d>=2) && (d-2<=good_ratio(m-min_row[d-2]))))) {

    /* for each relevant diagonal ... */
    for (k = lower; k <= upper; k++) {

      /* find a d on diagonal k */
      if (k==-d+DELTA+ORIGIN) {
        /* move down from the last d-1 on diagonal k+1 */
        row = last_d[k+1];
        /* op = INSERT; */

      } else if (k==d+DELTA+ORIGIN) {
        /* move right from the last d-1 on diagonal k-1 */
        row = last_d[k-1]-1;
        /* op = DELETE; */

      } else if ((last_d[k]-1<=last_d[k+1]) &&
                 (last_d[k]-1<=last_d[k-1]-1)) {
        /* substitution */
        row = last_d[k]-1;
        /* op = SUBSTITUTE; */

      } else if ((last_d[k-1]-1<=last_d[k+1]) &&
                 (last_d[k-1]-1<=last_d[k]-1)) {
        /* move right from the last d-1 on diagonal k-1 */
        row = last_d[k-1]-1;
        /* op = DELETE; */

      } else  {
        /* move left from the last d-1 on diagonal k+1 */
        row = last_d[k+1];
        /* op = INSERT; */

      }

      /* code common to the three cases */
      /* slide down the diagonal */

      col = row+k-ORIGIN;

      while ((row > 0) && (col > 0) && (s1[row-1]==s2[col-1]))
        { row--; col--; }

      temp_d[k] = row;

      if ((row == 0) && (col == 0)) {
        /* hit southeast corner; have the answer */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }
      if (row == 0) {
        /* hit first row; don't look further */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }

      if (col == 0) {
        /* hit last column; don't look further */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }
    }

    min_row[d] = last_d[ORIGIN+DELTA];
    min_diag[d] = ORIGIN+DELTA;
    for (k=lower; k<=upper; ++k)
      if (temp_d[k]<min_row[d]) {
        min_row[d] = temp_d[k];
        min_diag[d] = k;
      }

    for (k=lower; k<=upper; k++) {
      last_d[k] = temp_d[k];
    }

    --lower;
    ++upper;
  }

  /* report here the previous maximal match, stored in min_diag and min_row */
  while ((d>0) && (min_row[d-1]-min_row[d]<3))
    d--;

  *line1 = min_row[d]+offset1;
  *line2 = min_row[d]+min_diag[d]-ORIGIN+offset2;

        free(allocdSpace);

  return d;
}


int
Sim4::extend_fw(uchar *s1, uchar *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)
{
  int     col,                    /* column number */
  row,                    /* row number */
  max_d,                  /* bound on the length of the edit script
                           */
  d,                      /* current compressed distance */
  k,                      /* current diagonal */
  ORIGIN,
  lower,
  upper;
  int     *last_d, *temp_d;       /* column containing the last p */
  int     *max_row, *max_diag;    /* min (b)/ max (f) row (and diagonal) */
  /* reached for cost d=0, ... m.  */
  max_d = m+1;

  ORIGIN = m;
  for (row=0, col=0; col<n && row<m && (s1[row]==s2[col]); row++, col++)
    /*LINTED empty loop body*/; 

  if (row == m) {
    *line1 = row+offset1;
    *line2 = col+offset2;

    return 0;
  }
  if (col == n) {
    *line1 = row+offset1;
    *line2 = col+offset2;

    return 0;
  }

  int *allocdSpace = (int *)ckalloc((m+n+1+m+n+1+m+1+m+1) * sizeof(int));

  last_d   = allocdSpace;     //  m+n+1
  temp_d   = last_d + m+n+1;  //  m+n+1
  max_row  = temp_d + m+n+1;  //  m+1
  max_diag = max_row + m+1;   //  m+1

  for (k=0; k<=m+n; ++k) last_d[k]=-1;
  last_d[ORIGIN] = row;

  lower = ORIGIN - 1;
  upper = ORIGIN + 1;

  for (d=1; d<=m; d++)
    max_row[d] = -1;
  
  max_row[0] = last_d[ORIGIN];
  max_diag[0] = ORIGIN;

  d = 0;
  while ((++d<=max_d) && 
         ((d-1<=good_ratio(max_row[d-1])) || 
          ((d>=2) && (d-2<=good_ratio(max_row[d-2]))))) {

    /* for each relevant diagonal ... */
    for (k = lower; k <= upper; k++) {

      /* find a d on diagonal k */
      if (k==-d+ORIGIN) {

        /* move down from the last d-1 on diagonal k+1 */
        row = last_d[k+1]+1;
        /* op = DELETE; */
      } else if (k==d+ORIGIN) {

        /* move right from the last d-1 on diagonal k-1 */
        row = last_d[k-1];
        /* op = INSERT; */
      } else if ((last_d[k]>=last_d[k+1]) &&
                 (last_d[k]+1>=last_d[k-1])) {

        /* substitution */
        row = last_d[k]+1;
        /* op = SUBSTITUTE; */
      } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                 (last_d[k+1]>=last_d[k])) {

        /* move down from the last d-1 on diagonal k+1 */
        row = last_d[k+1]+1;
        /* op = DELETE; */
      } else {

        /* move right from the last d-1 on diagonal k-1 */
        row = last_d[k-1];
        /* op = INSERT; */
      }

      /* code common to the three cases */
      /* slide down the diagonal */

      col = row+k-ORIGIN;

      if (row>=0)
        while ((row < m) && (col < n) && (s1[row]==s2[col]))
          { row++; col++; }

      temp_d[k] = row;

      if ((row == m) && (col == n)) {
        /* hit southeast corner; have the answer */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }
      if (row == m) {
        /* hit last row; don't look further */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }

      if (col == n) {
        /* hit last column; don't look further */

        *line1 = row+offset1;
        *line2 = col+offset2;

        free(allocdSpace);

        return d;
      }
    }
    max_row[d] = last_d[ORIGIN];
    max_diag[d] = ORIGIN;
    for (k=lower; k<=upper; ++k)
      if (temp_d[k]>max_row[d]) {
        max_row[d] = temp_d[k];
        max_diag[d] = k;
      }

    for (k=lower; k<=upper; k++) {
      last_d[k] = temp_d[k];
    }

    --lower;
    ++upper;
  }

  /* report here the previous maximal match, stored in max_diag and max_row */

  while ((d>0) && (max_row[d]-max_row[d-1]<3))
    d--;

  *line1 = max_row[d]+offset1;
  *line2 = max_row[d]+max_diag[d]-ORIGIN+offset2;

  free(allocdSpace);

  return d;

  /*
     if ((d>2) && (max_row[d-1]-max_row[d-2]<3)) {
     *line1 = max_row[d-2]+offset1;
     *line2 = max_row[d-2]+max_diag[d-2]-ORIGIN+offset2;

        free(allocdSpace);
     
     return d-2;
     }
     
     *line1 = max_row[d-1]+offset1;
     *line2 = max_row[d-1]+max_diag[d-1]-ORIGIN+offset2;

        free(allocdSpace);

     return d-1;
     */
}
