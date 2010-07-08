#include "sim4.H"

//  This is used if _accurateSequences is enabled....and it's never
//  enabled.  The memory allocations here are NOT optimized.



#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "bio.h"
#include "ckalloc.H"


typedef  struct ValNode {
  void           *data;
  struct ValNode *next;
} *ValNodePtr;


void
link_to_data_list(void *data, ValNodePtr *head, ValNodePtr *prev) {
  ValNodePtr curr;

  curr = (ValNodePtr)ckalloc(sizeof(struct ValNode));
  curr->data = data;
  curr->next = NULL;

  if(*prev == NULL)  
    *head = curr;
  else
    (*prev)->next = curr;
  *prev = curr;
}


void
ValNodeFreeData(ValNodePtr data_list) {
  ValNodePtr   tmp_node;

  while ((tmp_node=data_list)!=NULL) {
    ckfree(tmp_node->data);
    data_list = data_list->next;
    ckfree(tmp_node); 
  }
}



int
Sim4::Xextend_bw(char *s1, char *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)  
{ 
  int     col,                    /* column number */
          row,                    /* row number */
          max_d,                  /* bound on the length of the edit script */
          d,                      /* current compressed distance */
          k,                      /* current diagonal */
          DELTA,                  /* n-m  */
          ORIGIN, 
          lower,
          upper;
  int     *last_d, *temp_d;       /* column containing the last p */
  int     *min_row, *min_diag;    /* min (b)/ max (f) row (and diagonal) */
                                  /* reached for cost d=0, ... m.  */
  coords ***trace_AG, ***trace_AC;
  coords  *AG_cell, *AC_cell, *newcoords;

  ValNodePtr  data_list = NULL, prev = NULL;

  DELTA = n-m;
  max_d = m+1;
          
  trace_AG = (coords ***)ckalloc((max_d+1)*sizeof(coords **)); 
  trace_AC = (coords ***)ckalloc((max_d+1)*sizeof(coords **)); 
  for (d=0; d<=max_d; d++) {
       trace_AG[d] = (coords **)ckalloc((m+n+1)*sizeof(coords *));
       trace_AC[d] = (coords **)ckalloc((m+n+1)*sizeof(coords *));
  }

  ORIGIN = m;

  trace_AG[0][ORIGIN+DELTA] = &last_AG;
  trace_AC[0][ORIGIN+DELTA] = &last_AC;

  for (row=m, col=n; row>0 && col>0 && (s1[row-1]==s2[col-1]); row--,col--)
        /*LINTED empty loop body*/; 
  for (k=n; (k>=2) && (k>=col); k--)
       if (!strncmp((char *)(s2+k-2),"AG",2)) {
           newcoords = (coords *)ckalloc(sizeof(coords)); 
           link_to_data_list((void *)newcoords, &data_list, &prev);

           newcoords->pos2 = k-DELTA+offset1 +1;    /* to compensate for -1 */
           newcoords->pos1 = k+offset2 +1;          /* refer to sim4b1.c */
           trace_AG[0][ORIGIN+DELTA] = newcoords;
       } else if (!strncmp((char *)(s2+k-2),"AC",2)) {
           newcoords = (coords *)ckalloc(sizeof(coords));
           link_to_data_list((void *)newcoords, &data_list, &prev);
           
           newcoords->pos2 = k-DELTA+offset1 +1;
           newcoords->pos1 = k+offset2 +1;
           trace_AC[0][ORIGIN+DELTA] = newcoords;
       } 

  if ((row == 0) || (col == 0)) {
       *line1 = row+offset1;
       *line2 = col+offset2; 
        (void)memcpy(&last_AG,trace_AG[0][ORIGIN+DELTA],sizeof(coords));
        (void)memcpy(&last_AC,trace_AC[0][ORIGIN+DELTA],sizeof(coords));
        ValNodeFreeData(data_list);
        free_coords(trace_AG,max_d+1);
        free_coords(trace_AC,max_d+1);

        return 0;
  }    
       
  last_d = (int *)ckalloc((m+n+1)*sizeof(int));
  temp_d = (int *)ckalloc((m+n+1)*sizeof(int));
  
  for (k=0; k<=m+n; ++k) last_d[k]=m+1;
  last_d[ORIGIN+DELTA] = row;
       
  lower = ORIGIN + DELTA - 1;
  upper = ORIGIN + DELTA + 1;
  
  min_row = (int *)ckalloc((m+1)*sizeof(int));
  min_diag = (int *)ckalloc((m+1)*sizeof(int));

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
                        AG_cell = trace_AG[d-1][k+1];
                        AC_cell = trace_AC[d-1][k+1];
               } else if (k==d+DELTA+ORIGIN) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1]-1;
                        /* op = DELETE; */
                        AG_cell = trace_AG[d-1][k-1];
                        AC_cell = trace_AC[d-1][k-1]; 
               } else if ((last_d[k]-1<=last_d[k+1]) &&
                          (last_d[k]-1<=last_d[k-1]-1)) {
                        /* substitution */
                        row = last_d[k]-1;
                        /* op = SUBSTITUTE; */
                        AG_cell = trace_AG[d-1][k];
                        AC_cell = trace_AC[d-1][k];
               } else if ((last_d[k-1]-1<=last_d[k+1]) &&
                          (last_d[k-1]-1<=last_d[k]-1)) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1]-1;
                        /* op = DELETE; */
                        AG_cell = trace_AG[d-1][k-1];
                        AC_cell = trace_AC[d-1][k-1];
               } else  { 
                        /* move left from the last d-1 on diagonal k+1 */
                        row = last_d[k+1];
                        /* op = INSERT; */
                        AG_cell = trace_AG[d-1][k+1];
                        AC_cell = trace_AC[d-1][k+1];
               }
  
               /* code common to the three cases */
               /* slide down the diagonal */
  
               col = row+k-ORIGIN;
  
               trace_AG[d][k] = AG_cell;
               trace_AC[d][k] = AC_cell;

               while ((row > 0) && (col > 0) && (s1[row-1]==s2[col-1])) {
                 if ((col>1) && !strncmp((char *)(s2+col-2),"AG",2)) {
                    newcoords = (coords *)ckalloc(sizeof(coords));
                    link_to_data_list((void *)newcoords, &data_list, &prev);

                    newcoords->pos1 = row + k - ORIGIN + offset2 +1;
                    newcoords->pos2 = row + offset1 +1; 
                    trace_AG[d][k] = newcoords;
                 } else if ((col>1) && !strncmp((char *)(s2+col-2),"AC",2)) {
                    newcoords = (coords *)ckalloc(sizeof(coords));
                    link_to_data_list((void *)newcoords, &data_list, &prev);
                    
                    newcoords->pos1 = row + k - ORIGIN + offset2 +1;
                    newcoords->pos2 = row + offset1 +1;
                    trace_AC[d][k] = newcoords;
                 } 
                 row--; col--; 
               }

               if ((col>1) && !strncmp((char *)(s2+col-2),"AG",2)) {
                    newcoords = (coords *)ckalloc(sizeof(coords));
                    link_to_data_list((void *)newcoords, &data_list, &prev);
                    
                    newcoords->pos1 = row + k - ORIGIN + offset2 +1;
                    newcoords->pos2 = row + offset1 +1;
                    trace_AG[d][k] = newcoords;
               } else if ((col>1) && !strncmp((char *)(s2+col-2),"AC",2)) {
                    newcoords = (coords *)ckalloc(sizeof(coords));
                    link_to_data_list((void *)newcoords, &data_list, &prev);
                    
                    newcoords->pos1 = row + k - ORIGIN + offset2 +1;
                    newcoords->pos2 = row + offset1 +1;
                    trace_AC[d][k] = newcoords;
               }

               temp_d[k] = row;
          
               if ((row == 0) && (col == 0)) {
                   /* hit southeast corner; have the answer */
          
                   (void)memcpy(&last_AG,trace_AG[d][k],sizeof(coords));
                   (void)memcpy(&last_AC,trace_AC[d][k],sizeof(coords));

                   ckfree(last_d);
                   ckfree(temp_d);
                   ckfree(min_row);
                   ckfree(min_diag);
                   ValNodeFreeData(data_list);
                   free_coords(trace_AG,max_d+1);
                   free_coords(trace_AC,max_d+1);
                        
                   *line1 = row+offset1;  
                   *line2 = col+offset2;  
       
                   return d;
               }
               if (row == 0) {
                   /* hit first row; don't look further */
  
                   (void)memcpy(&last_AG,trace_AG[d][k],sizeof(coords));
                   (void)memcpy(&last_AC,trace_AC[d][k],sizeof(coords));

                   ckfree(last_d);
                   ckfree(temp_d);
                   ckfree(min_row);
                   ckfree(min_diag);
                   ValNodeFreeData(data_list);
                   free_coords(trace_AG,max_d+1);
                   free_coords(trace_AC,max_d+1);

                   *line1 = row+offset1;  
                   *line2 = col+offset2;
       
                   return d;
               }          
  
               if (col == 0) {
                   /* hit last column; don't look further */
                   (void)memcpy(&last_AG,trace_AG[d][k],sizeof(coords));
                   (void)memcpy(&last_AC,trace_AC[d][k],sizeof(coords));

                   ckfree(last_d);
                   ckfree(temp_d);
                   ckfree(min_row);
                   ckfree(min_diag);
                   ValNodeFreeData(data_list);
                   free_coords(trace_AG,max_d+1);
                   free_coords(trace_AC,max_d+1);
        
                   *line1 = row+offset1;  
                   *line2 = col+offset2;  
  
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

     (void)memcpy(&last_AG,trace_AG[d][min_diag[d]],sizeof(coords));
     (void)memcpy(&last_AC,trace_AC[d][min_diag[d]],sizeof(coords));
                   
     ckfree(min_row);       
     ckfree(min_diag);
     ckfree(last_d);
     ckfree(temp_d); 
     ValNodeFreeData(data_list);
     free_coords(trace_AG,max_d+1);
     free_coords(trace_AC,max_d+1);
       
     return d;     
}                       

int
Sim4::Xextend_fw(char *s1, char *s2, int m, int n, int offset1, int offset2, int *line1, int *line2)              
{ 
  int     col,                    /* column number */
          row,                    /* row number */
          max_d,                  /* bound on the length of the edit script */       
          d,                      /* current compressed distance */
          k,                      /* current diagonal */
          ORIGIN,    
          lower,     
          upper;  
  int    *last_d, *temp_d;        /* column containing the last p */
  int    *max_row, *max_diag;     /* min (b)/ max (f) row (and diagonal) */
                                  /* reached for cost d=0, ... m.  */
  coords ***trace_GT, ***trace_CT;  
  coords  *GT_cell, *CT_cell, *newcoords;

  ValNodePtr  data_list = NULL, prev = NULL;

  max_d = m+1;

  trace_GT = (coords ***)ckalloc((max_d+1)*sizeof(coords **));
  trace_CT = (coords ***)ckalloc((max_d+1)*sizeof(coords **));
  for (d=0; d<=max_d; d++) {
       trace_GT[d] = (coords **)ckalloc((m+n+1)*sizeof(coords *));
       trace_CT[d] = (coords **)ckalloc((m+n+1)*sizeof(coords *));
  }    

  ORIGIN = m;
  trace_GT[0][ORIGIN] = &last_GT;
  trace_CT[0][ORIGIN] = &last_CT;

  for (row=0, col=0; col<n && row<m && (s1[row]==s2[col]); row++, col++)
        /*LINTED empty loop body*/; 
  for (k=0; (k<=n-2) && (k<=row); k++)
       if (!strncmp((char *)(s2+k),"GT",2)) {
              newcoords = (coords *)ckalloc(sizeof(coords));
              link_to_data_list((void *)newcoords, &data_list, &prev);
              newcoords->pos2 = k+offset1;
              newcoords->pos1 = k+offset2;
              trace_GT[0][ORIGIN] = newcoords;
       } else if (!strncmp((char *)(s2+k),"CT",2)) {
              newcoords = (coords *)ckalloc(sizeof(coords));
              link_to_data_list((void *)newcoords, &data_list, &prev);
              newcoords->pos2 = k+offset1;
              newcoords->pos1 = k+offset2;
              trace_CT[0][ORIGIN] = newcoords;
       }
       
  if ((row == m) || (col == n)){
       *line1 = row+offset1;
       *line2 = col+offset2;
        (void)memcpy(&last_GT,trace_GT[0][ORIGIN],sizeof(coords));
        (void)memcpy(&last_CT,trace_CT[0][ORIGIN],sizeof(coords));
        ValNodeFreeData(data_list);
        free_coords(trace_GT,max_d+1);
        free_coords(trace_CT,max_d+1);

        return 0;   
  }  
     
  last_d = (int *)ckalloc((m+n+1)*sizeof(int));
  temp_d = (int *)ckalloc((m+n+1)*sizeof(int));
     
  for (k=0; k<=m+n; ++k) last_d[k]=-1;
  last_d[ORIGIN] = row;
                   
  lower = ORIGIN - 1;
  upper = ORIGIN + 1;
  
  max_row = (int *)ckalloc((m+1)*sizeof(int));
  max_diag = (int *)ckalloc((m+1)*sizeof(int)); 
  
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
                        GT_cell = trace_GT[d-1][k+1];
                        CT_cell = trace_CT[d-1][k+1];
               } else if (k==d+ORIGIN) {
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1];
                        /* op = INSERT; */
                        GT_cell = trace_GT[d-1][k-1];
                        CT_cell = trace_CT[d-1][k-1];
               } else if ((last_d[k]>=last_d[k+1]) &&
                          (last_d[k]+1>=last_d[k-1])) {
                        /* substitution */
                        row = last_d[k]+1;
                        /* op = SUBSTITUTE; */ 
                        GT_cell = trace_GT[d-1][k];
                        CT_cell = trace_CT[d-1][k];
               } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                          (last_d[k+1]>=last_d[k])) {
                        /* move down from the last d-1 on diagonal k+1 */
                        row = last_d[k+1]+1;
                        /* op = DELETE; */
                        GT_cell = trace_GT[d-1][k+1];
                        CT_cell = trace_CT[d-1][k+1];
               } else  {
  
                        /* move right from the last d-1 on diagonal k-1 */
                        row = last_d[k-1];
                        /* op = INSERT; */
                        GT_cell = trace_GT[d-1][k-1];
                        CT_cell = trace_CT[d-1][k-1];
               }
  
               /* code common to the three cases */
               /* slide down the diagonal */
  
               col = row+k-ORIGIN;
  
               trace_GT[d][k] = GT_cell;
               trace_CT[d][k] = CT_cell;

               if (row>=0)
               while ((row < m) && (col < n) && (s1[row]==s2[col])) {
                  if ((col<n-1) && !strncmp((char *)(s2+col),"GT",2)) {
                     newcoords = (coords *)ckalloc(sizeof(coords));
                     link_to_data_list((void *)newcoords, &data_list, &prev);

                     newcoords->pos1 = row + k - ORIGIN + offset2;
                     newcoords->pos2 = row + offset1;
                     trace_GT[d][k] = newcoords;
                  } else if ((col<n-1) && !strncmp((char *)(s2+col),"CT",2)) { 
                     newcoords = (coords *)ckalloc(sizeof(coords));
                     link_to_data_list((void *)newcoords, &data_list, &prev);
              
                     newcoords->pos1 = row + k - ORIGIN + offset2;
                     newcoords->pos2 = row + offset1;
                     trace_CT[d][k] = newcoords;
                  } 
                     
                  row++; col++;
               }

               if ((col<n-1) && !strncmp((char *)(s2+col),"GT",2)) {
                     newcoords = (coords *)ckalloc(sizeof(coords));
                     link_to_data_list((void *)newcoords, &data_list, &prev);
                     
                     newcoords->pos1 = row + k - ORIGIN + offset2;
                     newcoords->pos2 = row + offset1;
                     trace_GT[d][k] = newcoords;
               } else if ((col<n-1) && !strncmp((char *)(s2+col),"CT",2)) {
                     newcoords = (coords *)ckalloc(sizeof(coords));
                     link_to_data_list((void *)newcoords, &data_list, &prev);
                     
                     newcoords->pos1 = row + k - ORIGIN + offset2;
                     newcoords->pos2 = row + offset1;
                     trace_CT[d][k] = newcoords;
               } 

               temp_d[k] = row;
          
               if ((row == m) && (col == n)) {
                   /* hit southeast corner; have the answer */
                   (void)memcpy(&last_GT,trace_GT[d][k],sizeof(coords));
                   (void)memcpy(&last_CT,trace_CT[d][k],sizeof(coords));

                   ValNodeFreeData(data_list);
                   free_coords(trace_GT,max_d+1);
                   free_coords(trace_CT,max_d+1);
                   ckfree(last_d);
                   ckfree(temp_d);
                   ckfree(max_row);
                   ckfree(max_diag);
                   *line1 = row+offset1;
                   *line2 = col+offset2;
                        
                   return d;
               } 
               if (row == m) {
                   /* hit last row; don't look further */
                   (void)memcpy(&last_GT,trace_GT[d][k],sizeof(coords));
                   (void)memcpy(&last_CT,trace_CT[d][k],sizeof(coords));
                   
                   ValNodeFreeData(data_list);
                   free_coords(trace_GT,max_d+1);
                   free_coords(trace_CT,max_d+1);
                   ckfree(temp_d);
                   ckfree(last_d);
                   ckfree(max_row);
                   ckfree(max_diag);
               
                   *line1 = row+offset1;
                   *line2 = col+offset2;
                        
                   return d;
               }        
               
               if (col == n) {
                   /* hit last column; don't look further */
                   (void)memcpy(&last_GT,trace_GT[d][k],sizeof(coords));
                   (void)memcpy(&last_CT,trace_CT[d][k],sizeof(coords));
                   
                   ValNodeFreeData(data_list);
                   free_coords(trace_GT,max_d+1);
                   free_coords(trace_CT,max_d+1);

                   ckfree(temp_d);
                   ckfree(last_d);
                   ckfree(max_row);
                   ckfree(max_diag);
                        
                   *line1 = row+offset1;
                   *line2 = col+offset2;
                        
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
     
     (void)memcpy(&last_GT,trace_GT[d][max_diag[d]],sizeof(coords));
     (void)memcpy(&last_CT,trace_CT[d][max_diag[d]],sizeof(coords));

     ckfree(max_row);     
     ckfree(max_diag);
     ckfree(last_d);
     ckfree(temp_d); 
     ValNodeFreeData(data_list);
     free_coords(trace_GT,max_d+1);
     free_coords(trace_CT,max_d+1);

                   
     return d;     
}                  
