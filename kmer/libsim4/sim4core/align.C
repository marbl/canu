#include "sim4.H"


#ifndef __lint
/*@unused@*/
static const char rcsid[] =
"$Id$";
#endif


void
Sim4::align_path(int i1, int j1, int i2, int j2, int dist, edit_script **head, edit_script **tail)
{
  int     *last_d, *temp_d,       /* forward vectors */
  *rlast_d, *rtemp_d;     /* backward vectors */

  edit_script *head1, *tail1, *head2, *tail2;
  int midc, rmidc;
  int start, lower, upper;
  int rstart, rlower, rupper;
  int c, k, row;
  int mi, mj, tmp, ll, uu;
  char flag;

  *head = *tail = NULL;

  /* Boundary cases */
  if (i1 == i2) {
    if (j1 == j2) *head = NULL;
    else {
      head1 = (edit_script *) ckalloc(sizeof(edit_script));
      head1->op_type = INSERT;
      head1->num = j2-j1;
      head1->next = NULL;
      *head = *tail = head1;
    }
    return;
  }

  if (j1 == j2) {
    head1 = (edit_script *) ckalloc(sizeof(edit_script));
    head1->op_type = DELETE;
    head1->num = i2-i1;
    head1->next = NULL;
    *head = *tail = head1;
    return;
  }

  if (dist <= 1) {
    start = j1-i1; 
    if (j2-i2 == j1-i1) {
      head1 = (edit_script *) ckalloc(sizeof(edit_script));
      head1->op_type = SUBSTITUTE;
      head1->num = i2-i1;
      head1->next = NULL;
      *head = *tail = head1;
    } else if (j2-j1 == i2-i1+1) {

      tmp = snake(start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *) ckalloc(sizeof(edit_script));
        head1->op_type = SUBSTITUTE;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *) ckalloc(sizeof(edit_script));
      head2->op_type = INSERT;
      head2->num = 1;

      if (*head) head1->next = head2; 
      else *head = head2; 
      *tail = head2;
      head2->next = NULL;

      if (i2-tmp) {
        head1 = head2;
        *tail = head2 = (edit_script *)ckalloc(sizeof(edit_script));
        head2->op_type = SUBSTITUTE;
        head2->num = i2-tmp;
        head2->next = NULL;
        head1->next = head2;
      } 
    } else if (j2-j1+1 == i2-i1) {

      tmp = snake(start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *) ckalloc(sizeof(edit_script));
        head1->op_type = SUBSTITUTE;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *) ckalloc(sizeof(edit_script));
      head2->op_type = DELETE;
      head2->num = 1;

      if (*head) head1->next = head2;
      else *head = head2;
      *tail = head2;
      head2->next = NULL;

      if (i2>tmp+1) {
        head1 = head2;
        *tail = head2 = (edit_script *)ckalloc(sizeof(edit_script));
        head2->op_type = SUBSTITUTE;
        head2->num = i2-tmp-1;
        head2->next = NULL;
        head1->next = head2;
      }
    } else {
      (void)fprintf(stderr, 
                    "align.c: warning: something wrong when aligning."); 
    }                
    return;
  }
  
  /* Divide the problem at the middle cost */
  midc = dist/2;
  rmidc = dist - midc;
  
  /* Compute the boundary diagonals */
  start = j1 - i1;
  lower = max(j1-i2, start-midc);
  upper = min(j2-i1, start+midc);
  rstart = j2-i2; 
  rlower = max(j1-i2, rstart-rmidc);
  rupper = min(j2-i1, rstart+rmidc);

  /* Allocate space for forward vectors */
  last_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
  temp_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;

  for (k=lower; k<=upper; k++)
    last_d[k] = -1;

  last_d[start] = snake(start,i1,i2,j2);
  
  /* Forward computation */  
  for (c=1; c<=midc; ++c) {
    ll = max(lower,start-c); 
    uu = min(upper,start+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll) {
        /* DELETE : down from (k+1,c-1) */
        row = last_d[k+1]+1;
      } else if (k == uu) {
        /* INSERT : right from (k-1,c-1) */
        row = last_d[k-1];
      } else if ((last_d[k]>=last_d[k+1]) &&
                 (last_d[k]+1>=last_d[k-1])) {
        /* SUBSTITUTE */
        row = last_d[k]+1;
      } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                 (last_d[k+1]>=last_d[k])) {
        /* DELETE */
        row = last_d[k+1]+1;
      } else {
        /* INSERT */
        row = last_d[k-1];
      }

      temp_d[k] = snake(k,row,i2,j2);
    } 
    for (k=ll; k<=uu; ++k)
      last_d[k] = temp_d[k];
  }

  /* Allocate space for backward vectors */
  rlast_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;
  rtemp_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;

  for (k=rlower; k<=rupper; k++)
    rlast_d[k] = i2+1;

  rlast_d[rstart] = rsnake(rstart,i2,i1,j1,len1);

  /* Backward computation */
  for (c=1; c<=rmidc; ++c) {
    ll = max(rlower,rstart-c);
    uu = min(rupper,rstart+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll) {
        /* INSERT : left from (k+1,c-1) */
        row = rlast_d[k+1];
      } else if (k == uu) {
        /* DELETE : up from (k-1,c-1) */  
        row = rlast_d[k-1]-1;
      } else if ((rlast_d[k]-1<=rlast_d[k+1]) &&
                 (rlast_d[k]-1<=rlast_d[k-1]-1)) {
        /* SUBSTITUTE */
        row = rlast_d[k]-1;
      } else if ((rlast_d[k-1]-1<=rlast_d[k+1]) &&
                 (rlast_d[k-1]-1<=rlast_d[k]-1)) {
        /* DELETE */
        row = rlast_d[k-1]-1;
      } else {
        /* INSERT */
        row = rlast_d[k+1];
      }
      
      rtemp_d[k] = rsnake(k,row,i1,j1,len1);
    }    
    for (k=ll; k<=uu; ++k)
      rlast_d[k] = rtemp_d[k];
  }

  /* Find (mi, mj) such that the distance from (i1, j1) to (mi, mj) is
     midc and the distance from (mi, mj) to (i2, j2) is rmidc.
     */

  flag = FALSE;
  mi = i1; mj = j1;
  ll = max(lower,rlower); uu = min(upper,rupper);
  for (k=ll; k<=uu; ++k) {
    if (last_d[k]>=rlast_d[k]) {
      if (last_d[k]-i1>=i2-rlast_d[k]) {
        mi = last_d[k]; mj = k+mi;
      } else {
        mi = rlast_d[k]; mj = k+mi; 
      } 
      flag = TRUE;

      break;
    }         
  }              
  free(last_d+lower); free(rlast_d+rlower);
  free(temp_d+lower); free(rtemp_d+rlower);
  
  if (flag) {    
    /* Find a path from (i1,j1) to (mi,mj) */
    align_path(i1,j1,mi,mj,midc,&head1,&tail1);
    
    /* Find a path from (mi,mj) to (i2,j2) */
    align_path(mi,mj,i2,j2,rmidc,&head2,&tail2);
    
    /* Join these two paths together */
    if (head1) tail1->next = head2;
    else head1 = head2;
  } else {  
    (void)fprintf(stderr, 
                  "align.c: warning: something wrong when dividing\n");
    head1 = NULL;
  }         
  *head = head1;
  if (head2) *tail = tail2;
  else *tail = tail1;
}



/* Condense_script - merge contiguous operations of the same type together */
void
Sim4::Condense_script(edit_script *head)
{
  edit_script *tp, *tp1;

  tp = head;
  while (tp != NULL) {
    while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
      tp->num = tp->num + tp1->num;
      tp->next = tp1->next;
      free(tp1);
    }
    tp = tp->next;
  }
}

/* Flip_script - reverse the script list */
void
Sim4::Flip_script(struct edit_script **script)
{
  struct edit_script *ep, *ahead, *behind;

  ahead = *script;
  ep = NULL;
  while (ahead!=NULL) {
    behind = ep;
    ep = ahead;
    ahead = ahead->next;
    ep->next = behind;
  }
  *script = ep;
}



void
Sim4::align_reverse(int *S)
{
  int auxi, *begi, *endi;
  
  begi = S; endi = S + *(S-1);
  while (begi < endi) {
    auxi = *begi;
    *begi = *--endi;
    *endi = auxi;
    begi++;      
  }               
  return;         
}                

/* Alignment display routine */



void
Sim4::Free_script(edit_script *head)
{
  edit_script *tp, *tp1;

  tp = head;
  while (tp != NULL) {
    tp1 = tp->next;
    free(tp);
    tp = tp1;
  }
}

