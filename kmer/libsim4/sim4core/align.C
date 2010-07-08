#include "sim4.H"

//  only for debugging
#include <sys/types.h>
#include <signal.h>

//  Define this to do bounds checking on the arrays used here
//#define CHECK_BOUNDS


#ifdef CHECK_BOUNDS
class boundedIntArray {
public:
  boundedIntArray(int offset, int size) {
    //fprintf(stderr, "boundedIntArray: offset=%d  size=%d\n", offset, size);
    _o  = offset;
    _m  = size;
    _a  = new int [_m];

    bzero(_a, sizeof(int) * _m);

    _a -= _o;
  };
  ~boundedIntArray() {
    _a += _o;
    delete [] _a;
  };

  int &operator[](int i) {
    if (i < _o) {
      fprintf(stderr, "********** i=%d o=%d\n", i, _o);
      exit(1);
    }
    if (i >= _o + _m) {
      fprintf(stderr, "********** i=%d o=%d m=%d\n", i, _o, _m);
      exit(1);
    }
    return(_a[i]);
  };

  int   _o;
  int   _m;
  int  *_a;
  int   _crud;
};
#endif



int
Sim4::align_get_dist(int i1, int j1, int i2, int j2, int limit) {
  
  //  Compute the boundary diagonals
  int start     = j1 - i1;
  int lower     = max(j1-i2, start-limit);
  int upper     = min(j2-i1, start+limit);
  int goal_diag = j2-i2;

  if (goal_diag > upper || goal_diag < lower)
    return(-1);

  //  Allocate space for forward vectors
#ifdef CHECK_BOUNDS
  boundedIntArray last_d(lower, upper-lower+1);
  boundedIntArray temp_d(lower, upper-lower+1);
#else
  int *last_d = (int *)ckalloc((upper-lower+1) * sizeof(int)) - lower;
  int *temp_d = (int *)ckalloc((upper-lower+1) * sizeof(int)) - lower;
#endif

  //  Initialization -- it's set to an easy to recognize value for
  //  debugging.
  for (int k=lower; k<=upper; ++k)
    last_d[k] = -2109876543;

  last_d[start] = snake(start, i1, i2, j2);

  if (last_d[goal_diag] >= i2) {
#ifndef CHECK_BOUNDS
    ckfree(last_d+lower);
    ckfree(temp_d+lower);
#endif
    return(0);
  }

  for (int c=1; c<=limit; ++c) {
    int ll = max(lower,start-c);
    int uu = min(upper, start+c);

    for (int k=ll; k<=uu; ++k) {
      int row;

      if (k == ll)
        row = last_d[k+1]+1;    //  DELETE
      else if (k == uu)
        row = last_d[k-1];      //  INSERT
      else if ((last_d[k]>=last_d[k+1]) &&
               (last_d[k]+1>=last_d[k-1]))
        row = last_d[k]+1;      //  SUBSTITUTE
      else if ((last_d[k+1]+1>=last_d[k-1]) &&
               (last_d[k+1]>=last_d[k]))
        row = last_d[k+1]+1;    //  DELETE
      else
        row = last_d[k-1];      //  INSERT
      
      temp_d[k] = snake(k,row,i2,j2);
    }     
    
    for (int k=ll; k<=uu; ++k)
      last_d[k] = temp_d[k];

    if (last_d[goal_diag] >= i2) {
      //  Free working vectors
#ifndef CHECK_BOUNDS
      ckfree(last_d+lower);
      ckfree(temp_d+lower);
#endif
      return(c);
    }
  }

#ifndef CHECK_BOUNDS
  ckfree(last_d+lower);
  ckfree(temp_d+lower);
#endif

  //  Ran out of distance limit
  return(-1);
}



int
Sim4::get_dist(int i1, int j1, int i2, int j2, int limit)
{
        int *SS, *DD, *II;
        int goal_diag;
        int c, k, t1, t2, t;
        int start, lower, upper;

        /* Compute the boundary diagonals */
        start = j1 - i1;
        lower = max(j1-i2, start-limit);
        upper = min(j2-i1, start+limit);
        goal_diag = j2-i2;

        if (goal_diag > upper || goal_diag < lower) {
                printf("The two sequences are not really similar.(1 %d; %d %d %d %d)\n", limit, i1, j1, i2, j2);
                printf("Please try exact phase 1 method\n.");
                exit(1);
        }

        /* Allocate space for forward vectors */
        SS = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
        DD = (int *)ckalloc((upper-lower+2)*sizeof(int)) - lower;
        II = (int *)ckalloc((upper-lower+2)*sizeof(int)) - lower + 1;

        /* Initialization */
        for (k=lower; k<=upper; ++k) SS[k] = -99999;
        for (k=lower; k<=upper+1; ++k) DD[k] = -99999;
        for (k=lower-1; k<=upper; ++k) II[k] = -99999;
        SS[start] = snake(start, i1, i2, j2);

        if (SS[goal_diag] >= i2) {
#ifdef STATS
                printf("get_dist = %d\n", 0);
#endif

                /* Free working vectors */
                free(SS+lower);
                free(DD+lower);
                free(II+lower-1);
                return 0;
        }

        for (c=1; c<=limit; ++c) {
                t = max(lower, start-c);
                t1 = II[t-1];
                for (k=t; k<=min(upper, start+c); ++k) {
                        t2 = II[k];
                        II[k] = max(t1, SS[k]);
                        t1 = t2;
                        DD[k] = max(DD[k+1]+1, SS[k]);
                        SS[k] = snake(k, min(j2-k,max(max(SS[k]+1, II[k]), DD[k]
)), i2, j2);
                }

                if (SS[goal_diag] >= i2) {
#ifdef STATS
                        printf("get_dist = %d\n", c);
#endif

                        /* Free working vectors */
                        free(SS+lower);
                        free(DD+lower);
                        free(II+lower-1);
                        return c;
                }
        }

        /* Ran out of distance limit */
        printf("Two sequences are not really similar.\n");
        printf("Please try exact phase 1.\n");
        exit(1);
}




void
Sim4::align_path(int i1, int j1,
                 int i2, int j2,
                 int dist,
                 edit_script **head,
                 edit_script **tail) {
#ifndef CHECK_BOUNDS
  int     *last_d = 0L;
  int     *temp_d = 0L;
  int     *rlast_d = 0L;
  int     *rtemp_d = 0L;
#endif
  edit_script *head1 = 0L;
  edit_script *tail1 = 0L;
  edit_script *head2 = 0L;
  edit_script *tail2 = 0L;

  //fprintf(stderr, "align_path()-- i1=%d j1=%d i2=%d j2=%d dist=%d\n", i1, j1, i2, j2, dist);

  int ll=0;
  int uu=0;

  *head = *tail = NULL;

  //  Boundary cases
  if (i1 == i2) {
    if (j1 == j2) {
      *head = NULL;
    } else {
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
    int start = j1-i1; 
    if (j2-i2 == j1-i1) {
      head1 = (edit_script *) ckalloc(sizeof(edit_script));
      head1->op_type = SUBSTITUTE;
      head1->num = i2-i1;
      head1->next = NULL;
      *head = *tail = head1;
    } else if (j2-j1 == i2-i1+1) {

      int tmp = snake(start,i1,i2,j2);
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

      int tmp = snake(start,i1,i2,j2);
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
      fprintf(stderr, "Sim4::align_path()-- warning: something wrong when aligning."); 
      *head = 0L;
      *tail = 0L;
    }
    return;
  }
  
  //  Divide the problem at the middle cost
  int midc  = dist/2;
  int rmidc = dist - midc;
  
  //  Compute the boundary diagonals
  int start = j1 - i1;
  int lower = max(j1-i2, start-midc);
  int upper = min(j2-i1, start+midc);
  int rstart = j2-i2; 
  int rlower = max(j1-i2, rstart-rmidc);
  int rupper = min(j2-i1, rstart+rmidc);


#if 0
  fprintf(stderr, "dist = %d\n", dist);
  fprintf(stderr, "midc = %d  rmidc = %d\n", midc, rmidc);
  fprintf(stderr, "j1 = %d\n", j1);
  fprintf(stderr, "i1 = %d\n", i1);
  fprintf(stderr, "j2 = %d\n", j2);
  fprintf(stderr, "i2 = %d\n", i2);
  fprintf(stderr, "start  = %d  lower  = %d  upper  = %d\n", start, lower, upper);
  fprintf(stderr, "rstart = %d  rlower = %d  rupper = %d\n", rstart, rlower, rupper);
#endif


  //  Allocate space for forward vectors
#ifdef CHECK_BOUNDS
  boundedIntArray last_d(lower, upper-lower+1);
  boundedIntArray temp_d(lower, upper-lower+1);
#else
  last_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
  temp_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
#endif

  for (int k=lower; k<=upper; k++)
    last_d[k] = -1;

  last_d[start] = snake(start,i1,i2,j2);
  
  //  Forward computation  
  for (int c=1; c<=midc; ++c) {
    ll = max(lower,start-c); 
    uu = min(upper,start+c);
    //fprintf(stderr, "c=%d ll=%d uu=%d\n", c, ll, uu);
    for (int k=ll; k<=uu; ++k) {
      int row;

      if (k == ll) {
        //  DELETE : down from (k+1,c-1)
        row = last_d[k+1]+1;
      } else if (k == uu) {
        //  INSERT : right from (k-1,c-1)
        row = last_d[k-1];
      } else if ((last_d[k]>=last_d[k+1]) &&
                 (last_d[k]+1>=last_d[k-1])) {
        //  SUBSTITUTE
        row = last_d[k]+1;
      } else if ((last_d[k+1]+1>=last_d[k-1]) &&
                 (last_d[k+1]>=last_d[k])) {
        //  DELETE
        row = last_d[k+1]+1;
      } else {
        //  INSERT
        row = last_d[k-1];
      }

      temp_d[k] = snake(k,row,i2,j2);

      //fprintf(stderr, "k=%d  row=%d  temp_d[k]=%d\n", k, row, temp_d[k]);
    } 
    for (int k=ll; k<=uu; ++k)
      last_d[k] = temp_d[k];
  }

  //  Allocate space for backward vectors
#ifdef CHECK_BOUNDS
  boundedIntArray rlast_d(rlower, rupper-rlower+1);
  boundedIntArray rtemp_d(rlower, rupper-rlower+1);
#else
  rlast_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;
  rtemp_d = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;
#endif

  for (int k=rlower; k<=rupper; k++)
    rlast_d[k] = i2+1;

  rlast_d[rstart] = rsnake(rstart,i2,i1,j1,i2+1);

  //  Backward computation
  for (int c=1; c<=rmidc; ++c) {
    ll = max(rlower,rstart-c);
    uu = min(rupper,rstart+c);
    for (int k=ll; k<=uu; ++k) {
      int row;

      if (k == ll) {
        //  INSERT : left from (k+1,c-1)
        row = rlast_d[k+1];
      } else if (k == uu) {
        //  DELETE : up from (k-1,c-1)  
        row = rlast_d[k-1]-1;
      } else if ((rlast_d[k]-1<=rlast_d[k+1]) &&
                 (rlast_d[k]-1<=rlast_d[k-1]-1)) {
        //  SUBSTITUTE
        row = rlast_d[k]-1;
      } else if ((rlast_d[k-1]-1<=rlast_d[k+1]) &&
                 (rlast_d[k-1]-1<=rlast_d[k]-1)) {
        //  DELETE
        row = rlast_d[k-1]-1;
      } else {
        //  INSERT
        row = rlast_d[k+1];
      }
      
      rtemp_d[k] = rsnake(k,row,i1,j1,i2+1);
    }    
    for (int k=ll; k<=uu; ++k)
      rlast_d[k] = rtemp_d[k];
  }

  //  Find (mi, mj) such that the distance from (i1, j1) to (mi, mj)
  //  is midc and the distance from (mi, mj) to (i2, j2) is rmidc.

  int flag = 0;
  int mi   = 0;
  int mj   = 0;

  ll = max(lower,rlower);
  uu = min(upper,rupper);

  //fprintf(stderr, "ll=%d uu=%d\n", ll, uu);

  for (int k=ll; k<=uu; ++k) {
    //fprintf(stderr, "last_d[%d] = %d   rlast_d[%d] = %d\n", k, last_d[k], k, rlast_d[k]);

    if (last_d[k] >= rlast_d[k]) {
      if (last_d[k] - i1 >= i2 - rlast_d[k]) {
        mi = last_d[k];
        mj = k+mi;
      } else {
        mi = rlast_d[k];
        mj = k+mi; 
      } 

      flag = 1;

      break;
    }
  }

#ifndef CHECK_BOUNDS
  ckfree(last_d  + lower);
  ckfree(rlast_d + rlower);
  ckfree(temp_d  + lower);
  ckfree(rtemp_d + rlower);
#endif

  //fprintf(stderr, "flag=%d mi=%d mj=%d\n", flag, mi, mj);

  if (flag == 0) {
    //fprintf(stderr, "Sim4::align_path()-- warning: something wrong when dividing\n");

#if 0
    //  Pick the middle k, keep going.

    int k= ll + (uu-ll) / 2;

    if (last_d[k] - i1 >= i2 - rlast_d[k]) {
      mi = last_d[k];
      mj = k+mi;
    } else {
      mi = rlast_d[k];
      mj = k+mi; 
    } 

#else
    //kill(getpid(), SIGSEGV);
    *head = 0L;
    *tail = 0L;
    return;
#endif
  }

  
  //  Find a path from (i1,j1) to (mi,mj)
  align_path(i1,j1,mi,mj,midc,&head1,&tail1);
    
  //  Find a path from (mi,mj) to (i2,j2)
  align_path(mi,mj,i2,j2,rmidc,&head2,&tail2);
    
  //  Join these two paths together
  if (head1)
    tail1->next = head2;
  else
    head1 = head2;

  *head = head1;

  if (head2)
    *tail = tail2;
  else
    *tail = tail1;
}


void
Sim4::path(int i1, int j1, char type1, int i2, int j2, char type2, int dist, edit_script **head, edit_script **tail)
{
        int *SS, *DD, *II;       /* Forward vectors */
        int *RS, *RD, *RI;       /* Backward vectors */

        edit_script *head1, *tail1, *head2, *tail2;
        int midc, rmidc;
        int start, lower, upper;
        int rstart, rlower, rupper;
        int c, k, t1, t2, t;
        int maxint;
        int mi, mj, mtype;
        char flag;

/*
printf("i1=%d,j1=%d,type1=%d,i2=%d,j2=%d,type2=%d,dist=%d\n",i1,j1,type1,i2,j2,type2,dist);
*/

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
           if (j2-i2 == j1-i1) {
                head1 = (edit_script *) ckalloc(sizeof(edit_script));
                head1->op_type = SUBSTITUTE;
                head1->num = i2-i1;
                head1->next = NULL;
                *head = *tail = head1;
           } else if (j2-i2 > j1-i1) {
                if (type1 == INSERT) {
                        head1 = (edit_script *) ckalloc(sizeof(edit_script));
                        head1->op_type = INSERT;
                        head1->num = 1;
                        head2 = (edit_script *) ckalloc(sizeof(edit_script));
                        head2->op_type = SUBSTITUTE;
                        head2->num = i2-i1;
                } else {
                        head1 = (edit_script *) ckalloc(sizeof(edit_script));
                        head1->op_type = SUBSTITUTE;
                        head1->num = i2-i1;
                        head2 = (edit_script *) ckalloc(sizeof(edit_script));
                        head2->op_type = INSERT;
                        head2->num = 1;
                }
                head1->next = head2;
                head2->next = NULL;
                *head = head1;
                *tail = head2;
           } else if (j2-i2 < j1-i1) {
                if (type1 == DELETE) {
                        head1 = (edit_script *) ckalloc(sizeof(edit_script));
                        head1->op_type = DELETE;
                        head1->num = 1;
                        head2 = (edit_script *) ckalloc(sizeof(edit_script));
                        head2->op_type = SUBSTITUTE;
                        head2->num = j2-j1;
                } else {
                        head1 = (edit_script *) ckalloc(sizeof(edit_script));
                        head1->op_type = SUBSTITUTE;
                        head1->num = j2-j1;
                        head2 = (edit_script *) ckalloc(sizeof(edit_script));
                        head2->op_type = DELETE;
                        head2->num = 1;
                }
                head1->next = head2;
                head2->next = NULL;
                *head = head1;
                *tail = head2;
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
        SS = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
        DD = (int *)ckalloc((upper-lower+2)*sizeof(int)) - lower;
        II = (int *)ckalloc((upper-lower+2)*sizeof(int)) - lower + 1;

        /* Forward computation */
        for (k=lower; k<=upper; ++k) SS[k] = -99999;
        for (k=lower; k<=upper+1; ++k) DD[k] = -99999;
        for (k=lower-1; k<=upper; ++k) II[k] = -99999;
        if (type1 == SUBSTITUTE) SS[start] = snake(start, i1, i2, j2);
        else if (type1 == DELETE) {
                DD[start] = i1;
                SS[start] = snake(start,i1,i2,j2);
        } else {
                II[start] = i1; 
                SS[start] = snake(start,i1,i2,j2);
        }

        for (c=1; c<=midc; ++c) {
                t = max(lower, start-c);
                t1 = II[t-1];
                for (k=t; k<=min(upper, start+c); ++k) {
                        t2 = II[k];
                        II[k] = max(t1, SS[k]);
                        t1 = t2;
                        DD[k] = max(DD[k+1]+1, SS[k]);
                        SS[k] = snake(k, min(j2-k,max(max(SS[k]+1, II[k]), DD[k])), i2, j2);
                }
        }

        /* Allocate space for backward vectors */
        RS = (int *)ckalloc((rupper-rlower+1)*sizeof(int)) - rlower;
        RD = (int *)ckalloc((rupper-rlower+2)*sizeof(int)) - rlower + 1;
        RI = (int *)ckalloc((rupper-rlower+2)*sizeof(int)) - rlower;

        /* Backward computation */
        maxint = i2 + dist + _estLen;
        for (k=rlower; k<=rupper; ++k) RS[k] = maxint;
        for (k=rlower-1; k<=rupper; ++k) RD[k] = maxint;
        for (k=rlower; k<=rupper+1; ++k) RI[k] = maxint;
        if (type2 == SUBSTITUTE)
           RI[rstart] = RD[rstart] = RS[rstart] = rsnake(rstart, i2, i1, j1, i2+1);
        else if (type2 == DELETE) RD[rstart] = i2;
        else RI[rstart] = i2;

        for (c=1; c<=rmidc; ++c) {
                t = max(rlower, rstart-c);
                t1 = RD[t-1];
                for (k=t; k<=min(rupper, rstart+c); ++k) {
#if 0
                        int x = min(min(RS[k]-1,RD[k]),RI[k]);
                        printf("<<<%d>>>", x);
                        assert(0<=x);
                        assert (x<=_estLen);
                        printf("%d", x);
#endif
                        RS[k] = rsnake(k, max(j1-k, min(min(RS[k]-1,RD[k]),RI[k])),i1,j1,i2+1);
                        t2 = RD[k];
                        RD[k] = min(t1-1, RS[k]);
                        t1 = t2;
                        RI[k] = min(RI[k+1], RS[k]);
                }
        }

        /* Find (mi, mj, mtype) such that
               the distance from (i1, j1, type1) to (mi, mj, mtype) is midc
           and the distance from (mi, mj, mtype) to (i2, j2, type2) is rmidc.
        */ 

        flag = 0;
        for (k=max(lower,rlower); k<=min(upper,rupper);++k) {

/*
printf("k=%d, SS=%d, RS=%d, DD=%d, RD=%d, II=%d, RI=%d\n",k,SS[k],RS[k],DD[k],RD[k],II[k],RI[k]);
*/

                if (SS[k]>=RS[k] || DD[k]>=RD[k] || II[k]>=RI[k]) {
                        if (DD[k]>=RD[k]) {
                                mi = DD[k];
                                mj = k+mi;
                                mtype = DELETE;
                        } else if (II[k] >= RI[k]) {
                                mi = II[k];
                                mj = k+mi;
                                mtype = INSERT;
                        } else {
                                mi = SS[k];
                                mj = k+mi;
                                mtype = SUBSTITUTE;
                        }

/*
                        printf("mi=%d, mj=%d, mtype=%d\n", mi, mj, mtype);
*/
                        flag = 1;
                        break;
                }
        }

        /* Free working vectors */
        free(SS+lower);
        free(DD+lower);
        free(II+lower-1);
        free(RS+rlower);
        free(RD+rlower-1);
        free(RI+rlower);

        if (flag) {
                /* Find a path from (i1,j1,type1) to (mi,mj,mtype) */
                path(i1,j1,type1,mi,mj,mtype,midc,&head1,&tail1);

                /* Find a path from (mi,mj,mtype) to (i2,j2,type2) */
                path(mi,mj,mtype,i2,j2,type2,rmidc,&head2,&tail2);

                /* Join these two paths together */
                if (head1) tail1->next = head2;
                else head1 = head2;
        } else {
                printf("Something wrong when dividing\n");
                head1 = NULL;
        }
        *head = head1;
        if (head2) *tail = tail2;
        else *tail = tail1;
}



//  Condense_script - merge contiguous operations of the same type together
void
Sim4::Condense_script(edit_script *head)
{
  edit_script *tp, *tp1;

  tp = head;
  while (tp != NULL) {
    while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
      tp->num = tp->num + tp1->num;
      tp->next = tp1->next;
      ckfree(tp1);
    }
    tp = tp->next;
  }
}

//  Flip_script - reverse the script list
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



void
Sim4::Free_script(edit_script *head)
{
  edit_script *tp, *tp1;

  tp = head;
  while (tp != NULL) {
    tp1 = tp->next;
    ckfree(tp);
    tp = tp1;
  }
}

