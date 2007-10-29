#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#include "bio.h"

//  Liliana Florea's halign (a sim4-derivitive).

#define DEL 0
#define INS 1
#define SUB 2

#ifdef min
#undef min
#endif
#define min(x,y)      ((x)<=(y) ? (x):(y))

#ifdef max
#undef max
#endif
#define max(x,y)      ((x)>=(y) ? (x):(y))


typedef struct edit_script {
  int  op_type;   /* SUB, INS or DEL */
  int  num;        /* Number of operations */
  struct edit_script *next;
} edit_script;

typedef struct edit_script_list {
  int          offset1, offset2;
  int          len1, len2;
  int          score;
  int          first;
  edit_script *script;
} edit_script_list;



static
int
snake(const char *seq1, const char *seq2, int k, int x, int endx, int endy) {
  int y;

  if (x<0) return x;
  y = x+k;
  while ((x < endx) &&
         (y < endy) &&
         (toUpper[seq1[x]] == toUpper[seq2[y]])) {
    ++x;
    ++y;
  }
  return x;
}


static
int
rsnake(const char *seq1, const char *seq2, int k, int x, int startx, int starty, int M, int N) {
  int y;

  if (x>M) return x;
  if ((startx<0) || (starty<0))
    fprintf(stderr, "halign::rsnake()-- TROUBLE!!! startx:  %5d,  starty:  %5d\n",startx, starty);
  if ((x>M) || (x+k>N))
    fprintf(stderr, "halign::rsnake()-- TROUBLE!!! x:  %5d,  y:  %5d\n",x,x+k);

  y = x+k;
  while ((x>startx) &&
         (y>starty) &&
         (toUpper[seq1[x-1]] == toUpper[seq2[y-1]])) {
    --x;
    --y;
  }
  return x;
}


static
int
align_get_dist(const char *seq1,
               const char *seq2,
               int i1, int j1,
               int i2, int j2,
               int limit,
               void *ph) {
  int *last_d, *temp_d;
  int goal_diag, ll, uu;
  int c, k, row;
  int start;
  int lower, upper;

  /* Compute the boundary diagonals */
  start = j1 - i1;
  lower = max(j1-i2, start-limit);
  upper = min(j2-i1, start+limit);
  goal_diag = j2-i2;

  if (goal_diag > upper || goal_diag < lower) {
    fprintf(stderr, "The two sequences are not really similar.\n");
    fprintf(stderr, "Please try an exact aligning method.\n");
    exit(1);
  }

  /* Allocate space for forward vectors */
  last_d = (int *)palloc2((upper-lower+1)*sizeof(int), ph) - lower;
  temp_d = (int *)palloc2((upper-lower+1)*sizeof(int), ph) - lower;

  /* Initialization */
  for (k=lower; k<=upper; ++k) last_d[k] = INT_MIN;
  last_d[start] = snake(seq1, seq2, start, i1, i2, j2);

  if (last_d[goal_diag] >= i2)
    return 0;

  for (c=1; c<=limit; ++c) {
    ll = max(lower,start-c); uu = min(upper, start+c);
    for (k=ll; k<=uu; ++k) {
      if (k == ll)
        row = last_d[k+1]+1;    /* DELETE */
      else if (k == uu)
        row = last_d[k-1];      /* INSERT */
      else if ((last_d[k]>=last_d[k+1]) &&
               (last_d[k]+1>=last_d[k-1]))
        row = last_d[k]+1;      /*SUBSTITUTE */
      else if ((last_d[k+1]+1>=last_d[k-1]) &&
               (last_d[k+1]>=last_d[k]))
        row = last_d[k+1]+1;    /* DELETE */
      else
        row = last_d[k-1];      /* INSERT */

      temp_d[k] = snake(seq1,seq2,k,row,i2,j2);
    }

    for (k=ll; k<=uu; ++k) last_d[k] = temp_d[k];

    if (last_d[goal_diag] >= i2)
      return c;
  }

  /* Ran out of distance limit */
  return -1;
}


static
void
align_path(const char *seq1,
           const char *seq2,
           int i1, int j1,
           int i2, int j2,
           int dist,
           edit_script **head,
           edit_script **tail,
           int M,
           int N,
           void *ph) {

  int     *last_d, *temp_d;       /* forward vectors */
  int    *rlast_d, *rtemp_d;     /* backward vectors */

  edit_script *head1, *tail1, *head2, *tail2;
  int midc, rmidc;
  int start;
  int lower, upper;
  int rstart, rlower, rupper;
  int c, k, row;
  int mi, mj, tmp, ll, uu;
  char flag;

  *head = *tail = NULL;

  /* Boundary cases */
  if (i1 == i2) {
    if (j1 == j2) *head = NULL;
    else {
      head1 = (edit_script *)palloc2(sizeof(edit_script), ph);
      head1->op_type = INS;
      head1->num = j2-j1;
      head1->next = NULL;
      *head = *tail = head1;
    }
    return;
  }

  if (j1 == j2) {
    head1 = (edit_script *)palloc2(sizeof(edit_script), ph);
    head1->op_type = DEL;
    head1->num = i2-i1;
    head1->next = NULL;
    *head = *tail = head1;
    return;
  }

  if (dist <= 1) {
    start = j1-i1;
    if (j2-i2 == j1-i1) {
      head1 = (edit_script *)palloc2(sizeof(edit_script), ph);
      head1->op_type = SUB;
      head1->num = i2-i1;
      head1->next = NULL;
      *head = *tail = head1;
    } else if (j2-j1 == i2-i1+1) {

      tmp = snake(seq1,seq2,start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *)palloc2(sizeof(edit_script), ph);
        head1->op_type = SUB;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *)palloc2(sizeof(edit_script), ph);
      head2->op_type = INS;
      head2->num = 1;

      if (*head) head1->next = head2;
      else *head = head2;
      *tail = head2;
      head2->next = NULL;

      if (i2-tmp) {
        head1 = head2;
        *tail = head2 = (edit_script *)palloc2(sizeof(edit_script), ph);
        head2->op_type = SUB;
        head2->num = i2-tmp;
        head2->next = NULL;
        head1->next = head2;
      }
    } else if (j2-j1+1 == i2-i1) {

      tmp = snake(seq1,seq2,start,i1,i2,j2);
      if (tmp>i1) {
        head1 = (edit_script *)palloc2(sizeof(edit_script), ph);
        head1->op_type = SUB;
        head1->num = tmp-i1;
        *head = head1;
      }
      head2 = (edit_script *)palloc2(sizeof(edit_script), ph);
      head2->op_type = DEL;
      head2->num = 1;

      if (*head) head1->next = head2;
      else *head = head2;
      *tail = head2;
      head2->next = NULL;

      if (i2>tmp+1) {
        head1 = head2;
        *tail = head2 = (edit_script *)palloc2(sizeof(edit_script), ph);
        head2->op_type = SUB;
        head2->num = i2-tmp-1;
        head2->next = NULL;
        head1->next = head2;
      }
    } else {
      fprintf(stderr, "halign::align_path()-- warning: something wrong when aligning.");
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
  last_d = (int *)palloc2((upper-lower+1)*sizeof(int), ph) - lower;
  temp_d = (int *)palloc2((upper-lower+1)*sizeof(int), ph) - lower;

  for (k=lower; k<=upper; k++) last_d[k] = -1;
  last_d[start] = snake(seq1,seq2,start,i1,i2,j2);

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

      temp_d[k] = snake(seq1,seq2,k,row,i2,j2);
    }
    for (k=ll; k<=uu; ++k)
      last_d[k] = temp_d[k];
  }

  /* Allocate space for backward vectors */
  rlast_d = (int *)palloc2((rupper-rlower+1)*sizeof(int), ph) - rlower;
  rtemp_d = (int *)palloc2((rupper-rlower+1)*sizeof(int), ph) - rlower;

  for (k=rlower; k<=rupper; k++) rlast_d[k] = i2+1;
  rlast_d[rstart] = rsnake(seq1,seq2,rstart,i2,i1,j1,M,N);

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

      rtemp_d[k] = rsnake(seq1,seq2,k,row,i1,j1,M,N);
    }
    for (k=ll; k<=uu; ++k)
      rlast_d[k] = rtemp_d[k];
  }

  /* Find (mi, mj) such that the distance from (i1, j1) to (mi, mj) is
     midc and the distance from (mi, mj) to (i2, j2) is rmidc.
  */

  flag = 0;
  mi = i1; mj = j1;
  ll = max(lower,rlower);
  uu = min(upper,rupper);

  for (k=ll; k<=uu; ++k) {
    if (last_d[k]>=rlast_d[k]) {
      if (last_d[k]-i1>=i2-rlast_d[k]) {
        mi = last_d[k]; mj = k+mi;
      } else {
        mi = rlast_d[k]; mj = k+mi;
      }
      flag = 1;

      break;
    }
  }

  if (flag) {
    /* Find a path from (i1,j1) to (mi,mj) */
    align_path(seq1,seq2,i1,j1,mi,mj,midc,&head1,&tail1,M,N,ph);

    /* Find a path from (mi,mj) to (i2,j2) */
    align_path(seq1,seq2,mi,mj,i2,j2,rmidc,&head2,&tail2,M,N,ph);

    /* Join these two paths together */
    if (head1) tail1->next = head2;
    else head1 = head2;
  } else {
    fprintf(stderr, "halign::align_path()-- warning: something wrong when dividing\n");
    head1 = NULL;
  }
  *head = head1;
  if (head2) *tail = tail2;
  else *tail = tail1;
}





void
halign(const char *seq1,
       const char *seq2,
       const int   len1,
       const int   len2,
       char *alnline1,
       char *alnline2) {
  edit_script      *head, *tail, *tp;
  int               i;
  void             *ph;

  ph = pallochandle(0);

  align_path(seq1, seq2,
             0, 0,
             len1, len2,
             align_get_dist(seq1, seq2, 0, 0, len1, len2, len1+len2, ph),
             &head, &tail,
             len1, len2,
             ph);
        
  /* generate the alignment(s) */

  *alnline1 = 0;
  *alnline2 = 0;

  for (tp=head; tp; tp=tp->next) {
    switch (tp->op_type) {
      case SUB:
        for (i=0; i<tp->num; i++) {
          if (toUpper[*seq1] == toUpper[*seq2]) {
            *alnline1 = toLower[*seq1];
            *alnline2 = toLower[*seq2];
          } else {
            *alnline1 = toUpper[*seq1];
            *alnline2 = toUpper[*seq2];
          }
          seq1++;
          seq2++;
          alnline1++;
          alnline2++;
        }
        break; 

      case INS:
        for (i=0; i<tp->num; i++) {
          *alnline1 = '-';
          *alnline2 = toUpper[*seq2];
          seq2++;
          alnline1++;
          alnline2++;
        }
        break;  

      case DEL:
        for (i=0; i<tp->num; i++) {
          *alnline2 = '-';
          *alnline1 = toUpper[*seq1];
          seq1++;
          alnline1++;
          alnline2++;
        }
        break;

      default:
        fprintf(stderr, "halign::halign()-- unrecognized op_type in script. %d\n", tp->op_type);
        exit(0);
    }
  }
  *alnline1 = 0;
  *alnline2 = 0;

  pfree2(ph);
  pfreehandle(ph);
}
