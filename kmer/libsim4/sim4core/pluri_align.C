//#include <limits.h>
#include "sim4.H"


int
Sim4::align_get_dist(int i1, int j1, int i2, int j2, int limit)
{
  int *last_d, *temp_d;
  int goal_diag, ll, uu;
  int c, k, row;
  int start, lower, upper;
  
  /* Compute the boundary diagonals */
  start = j1 - i1;
  lower = max(j1-i2, start-limit);
  upper = min(j2-i1, start+limit);
  goal_diag = j2-i2;
  
  if (goal_diag > upper || goal_diag < lower)
    return(-1);
  
  /* Allocate space for forward vectors */ 
  last_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;
  temp_d = (int *)ckalloc((upper-lower+1)*sizeof(int)) - lower;

  /* Initialization */
  for (k=lower; k<=upper; ++k)
    last_d[k] = -2147483647;  //  INT_MIN;

  last_d[start] = snake(start, i1, i2, j2);
  
  if (last_d[goal_diag] >= i2) {
    /* Free working vectors */
    free(last_d+lower);
    free(temp_d+lower);
    return 0;
  }

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
      
      temp_d[k] = snake(k,row,i2,j2);
    }     
    
    for (k=ll; k<=uu; ++k) last_d[k] = temp_d[k];

    if (last_d[goal_diag] >= i2) {
      /* Free working vectors */
      free(last_d+lower);
      free(temp_d+lower);
      return c;
    }
  }

  free(last_d+lower);
  free(temp_d+lower);

  /* Ran out of distance limit */
  return -1;
}








/* Condense_both_Ends  --  merge contiguous operations of the same type    */
/* together; return both new ends of the chain.                            */
void
Sim4::Condense_both_Ends(edit_script **head,
                         edit_script **tail,
                         edit_script **prev) {
  edit_script *tp, *tp1;

  tp = *head; *prev = NULL;
  while (tp != NULL) {
    while (((tp1 = tp->next) != NULL) && (tp->op_type == tp1->op_type)) {
      tp->num = tp->num + tp1->num;
      tp->next = tp1->next;
      free(tp1);
    }
    if (tp->next) *prev = tp;
    else *tail = tp;
    tp = tp->next;
  }
}




void
Sim4::pluri_align(int *dist_ptr,
                  Exon *theExons,
                  struct edit_script_list **Aligns,
                  sim4_stats_t *st) {
  int    i, end1, end2, diff, ali_dist;
  uchar *a, *b;

  Exon  *thisExon = theExons;
  Exon  *nextExon;

  int    EditDistance = 0;     //  Sum of all tmpi, previously known as TMPI
  int    AlignmentLength = 0;

  struct edit_script_list *enew;
  struct edit_script      *head;
  struct edit_script      *tmp_script;
  struct edit_script      *left;
  struct edit_script      *right;
  struct edit_script      *prev;

  st->numberOfMatches = 0;
  st->numberOfNs      = 0;

  head = NULL;
  *Aligns = NULL;
  *dist_ptr = ali_dist = 0;

  end1 = len1; 
  end2 = len2;

  nextExon = thisExon->next_exon;

  while (nextExon && nextExon->toGEN) {
#if 0
    if ((diff=thisExon->frEST-nextExon->toEST-1)!=0) {
      if (thisExon->toGEN) {
        enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
        enew->next_script = *Aligns;
        *Aligns = enew;
        (*Aligns)->script = head;
        (*Aligns)->offset1 = thisExon->frGEN;
        (*Aligns)->offset2 = thisExon->frEST;
        (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
        (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
        (*Aligns)->score = ali_dist;
        ali_dist = 0;
        head = NULL;
      }
      end1 = nextExon->toGEN;
      end2 = nextExon->toEST;  

    } else if (((diff=thisExon->frGEN-nextExon->toGEN-1)!=0) &&
               thisExon->toGEN) {
          struct edit_script  *newthing;

      newthing = (edit_script *) ckalloc(sizeof(edit_script));
      newthing->op_type = DELETE;
      newthing->num = diff;
      newthing->next = head;
      head = newthing;
    } else if (diff) 
      end1 = nextExon->toGEN;
#else
    diff = thisExon->frEST - nextExon->toEST - 1;
      
    if (diff != 0) {
      if (thisExon->toGEN) {
        enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
        enew->next_script = *Aligns;
        *Aligns = enew;
        (*Aligns)->script = head;
        (*Aligns)->offset1 = thisExon->frGEN;
        (*Aligns)->offset2 = thisExon->frEST;
        (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
        (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
        (*Aligns)->score = ali_dist;
        ali_dist = 0;
        head = NULL;
      }
      end1 = nextExon->toGEN;
      end2 = nextExon->toEST;  
    } else {
      diff = thisExon->frGEN - nextExon->toGEN - 1;
      if (diff != 0) {
        if (thisExon->toGEN) {
          struct edit_script  *newthing;

          newthing = (edit_script *) ckalloc(sizeof(edit_script));
          newthing->op_type = DELETE;
          newthing->num     = diff;
          newthing->next    = head;
          head = newthing;
        } else {
          end1 = nextExon->toGEN;
        }
      }
    }
#endif

    diff = align_get_dist(nextExon->frGEN-1, nextExon->frEST-1,
                          nextExon->toGEN, nextExon->toEST,
                          max(1000, .2*(nextExon->toEST - nextExon->frEST + 1)));


    //  Return if the alignment fails.
    //
    if (diff < 0) {
      st->numberOfMatches = 0;
      st->numberOfNs      = 0;
      st->percentID       = -1;

      //fprintf(stderr, "The two sequences are not really similar.\nPlease try an exact method.\n");
      
      //  XXX:  MEMORY LEAK!  Should free the Aligns!
      //
      *Aligns             = 0L;

      return;
    }

#ifdef STATS
    if (diff>P*(nextExon->toEST-nextExon->frEST+1))
      (void)printf("Warning: Distance threshold on segment exceeded.\n");
#endif

    align_path(nextExon->frGEN-1, nextExon->frEST-1, 
               nextExon->toGEN, nextExon->toEST, diff, &left, &right);

    Condense_both_Ends(&left, &right, &prev);
    
    if (!thisExon->toGEN && right->op_type == DELETE) {
      /* remove gaps at end of alignment */
      diff -= 0+right->num;         /* subtract GAP_OPEN = 0 */
      nextExon->toGEN -= right->num;
      end1 -= right->num; 
      if (head && (head->op_type == DELETE)) 
        head->num += right->num;
      free(right); prev->next = NULL; 
      right = prev;
    } 

    if ((!nextExon->next_exon || !nextExon->next_exon->toGEN) &&
        left && (left->op_type == DELETE)) {
      diff -= 0+left->num;          /* subtract GAP_OPEN = 0 */
      nextExon->frGEN += left->num;

      tmp_script = left->next; 
      if (right == left)
        right = tmp_script;
      free(left);
      left = tmp_script; 
    }
    
    *dist_ptr += diff;
    ali_dist += diff;

    a = seq1 + nextExon->frGEN - 1;
    b = seq2 + nextExon->frEST - 1;

    nextExon->numMatches = 0;
    nextExon->numNs      = 0;
    nextExon->numInDel   = 0;
    nextExon->numEdits   = 0;

    tmp_script = left;  

#ifdef DEBUG
    fprintf(stderr, "a=%s\nb=%s\n", a, b);
#endif

    while (tmp_script) {
      switch (tmp_script->op_type) {
      case  DELETE:
#ifdef DEBUG
        fprintf(stderr, "delete %d\n", tmp_script->num);
#endif
        nextExon->numInDel += tmp_script->num;
        nextExon->numEdits += tmp_script->num;
        a                  += tmp_script->num;
        break;
      case  INSERT:
#ifdef DEBUG
        fprintf(stderr, "insert %d\n", tmp_script->num);
#endif
        nextExon->numInDel += tmp_script->num;
        nextExon->numEdits += tmp_script->num;
        b                  += tmp_script->num;
        break;
      case  SUBSTITUTE:

        //  Count the number of matches and edits.
        //
        //  An edit is a true substitute -- a base for a different base,
        //  not a base for an 'n'.
        //
        for (i=0; i<tmp_script->num; ++i, ++a, ++b) {

#ifdef DEBUG
          fprintf(stderr, "pair=%c%c\n", *a, *b);
#endif

          if (((*a == 'N') || (*a == 'n')) && ((*b == 'N') || (*b == 'n'))) {
            //  Got an 'n'.  It isn't a match and it isn't an edit.
            //
            nextExon->numNs++;
          } else {
            if (*a == *b) {
              //  Got a match.
              nextExon->numMatches++;
            } else {
              //  Got a substitution
              nextExon->numEdits++;
            }
          }
        }
        break;
      }         
#ifdef DEBUG
      fprintf(stderr, "next script.\n");
#endif
      tmp_script = tmp_script->next;
    }

#ifdef DEBUG
    fprintf(stderr, "next exon.\n");
#endif

    nextExon->alignmentLength = (nextExon->toGEN - nextExon->frGEN + 1 +
                                 nextExon->toEST - nextExon->frEST + 1 +
                                 nextExon->numInDel);
    nextExon->percentID = (int)floor(100*(1 - 2 * nextExon->numEdits /
                                          (float)(nextExon->alignmentLength)));


    st->numberOfMatches += nextExon->numMatches;
    st->numberOfNs      += nextExon->numNs;

    EditDistance    += nextExon->numEdits;
    AlignmentLength += (nextExon->toGEN - nextExon->frGEN + 1 +
                        nextExon->toEST - nextExon->frEST + 1 +
                        nextExon->numInDel);

    right->next = head;
    head = left;

    thisExon = nextExon;
    nextExon = thisExon->next_exon;
  }


  /* at the beginning of the sequences */
  if (nextExon!=NULL) {

    if ((diff=thisExon->frEST-nextExon->toEST-1)!=0 && (diff!=len2)) {
      enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
      enew->next_script = *Aligns;
      *Aligns = enew;
      (*Aligns)->offset1 = thisExon->frGEN;
      (*Aligns)->offset2 = thisExon->frEST;
      (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
      (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
      (*Aligns)->script = head;
      (*Aligns)->score = ali_dist;
      
    } else if (diff!=len2) {

      /* modified to cut introns at the beginning of the sequence */
      enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
      enew->next_script = *Aligns;
      *Aligns = enew;
      (*Aligns)->offset1 = thisExon->frGEN;
      (*Aligns)->offset2 = 1;
      (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
      (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
      (*Aligns)->script = head;
      (*Aligns)->score = ali_dist;
    }
  }

  st->percentID = 0;
  if (AlignmentLength > 0)
    st->percentID = (int)floor(100 * (1 - 2 * (double)EditDistance / (double)AlignmentLength));
}




void
Sim4::updateStatistics(Exon *theExon,
                       sim4_stats_t *st) {

  theExon = theExon->next_exon;

  st->numberOfMatches = 0;
  st->numberOfNs      = 0;

  int    EditDistance    = 0;
  int    AlignmentLength = 0;

  while (theExon && theExon->toGEN) {
    st->numberOfMatches += theExon->numMatches;
    st->numberOfNs      += theExon->numNs;

    EditDistance    += theExon->numEdits;
    AlignmentLength += (theExon->toGEN - theExon->frGEN + 1 +
                        theExon->toEST - theExon->frEST + 1 +
                        theExon->numInDel);

    theExon = theExon->next_exon;
  }

  st->percentID = 0;
  if (AlignmentLength > 0)
    st->percentID = (int)floor(100 * (1 - 2 * (double)EditDistance / (double)AlignmentLength));
}
