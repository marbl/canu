#include <math.h>
#include "sim4.H"


//  Condense_both_Ends -- merge contiguous operations of the same type
//  together; return both new ends of the chain.
//
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
      ckfree(tp1);
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
  char *a, *b;

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

  head      = 0L;
  *Aligns   = 0L;
  *dist_ptr = ali_dist = 0;

  end1 = _genLen; 
  end2 = _estLen;

  nextExon = thisExon->next_exon;

  while (nextExon && nextExon->toGEN) {
    diff = thisExon->frEST - nextExon->toEST - 1;
      
    if (diff != 0) {
      if (thisExon->toGEN) {
        enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
        enew->next_script = *Aligns;
        *Aligns = enew;
        (*Aligns)->script  = head;
        (*Aligns)->offset1 = thisExon->frGEN;
        (*Aligns)->offset2 = thisExon->frEST;
        (*Aligns)->len1    = end1-(*Aligns)->offset1+1;
        (*Aligns)->len2    = end2-(*Aligns)->offset2+1;
        (*Aligns)->score   = ali_dist;
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

    diff = align_get_dist(nextExon->frGEN-1, nextExon->frEST-1,
                          nextExon->toGEN, nextExon->toEST,
                          max(1000, (int)(.2*(nextExon->toEST - nextExon->frEST + 1))));


    //  Return if the alignment fails.
    //
    if (diff < 0) {
      st->numberOfMatches = 0;
      st->numberOfNs      = 0;
      st->percentID       = -1;
      *Aligns             = 0L;
      return;
    }

#ifdef STATS
    if (diff > P * (nextExon->toEST - nextExon->frEST + 1))
      (void)printf("Warning: Distance threshold on segment exceeded.\n");
#endif

    align_path(nextExon->frGEN-1, nextExon->frEST-1, 
               nextExon->toGEN, nextExon->toEST, diff, &left, &right);

    //  Return if the alignment fails -- this occurred once aligning
    //  dros frags to dros using snapper.  Snapper was giving the wrong
    //  sequence for the seeds it also supplied.
    //
    if ((left == 0L) || (right == 0L)) {
      st->numberOfMatches = 0;
      st->numberOfNs      = 0;
      st->percentID       = -1;
      *Aligns             = 0L;
      return;
    }

    Condense_both_Ends(&left, &right, &prev);
    
    if (!thisExon->toGEN && right->op_type == DELETE) {
      /* remove gaps at end of alignment */
      diff -= 0+right->num;         /* subtract GAP_OPEN = 0 */
      nextExon->toGEN -= right->num;
      end1 -= right->num; 
      if (head && (head->op_type == DELETE)) 
        head->num += right->num;
      ckfree(right);
      prev->next = NULL; 
      right = prev;
    } 

    if ((!nextExon->next_exon || !nextExon->next_exon->toGEN) &&
        left && (left->op_type == DELETE)) {
      diff -= 0+left->num;          /* subtract GAP_OPEN = 0 */
      nextExon->frGEN += left->num;

      tmp_script = left->next; 
      if (right == left)
        right = tmp_script;
      ckfree(left);
      left = tmp_script; 
    }
    
    *dist_ptr += diff;
    ali_dist += diff;

    a = _genSeq + nextExon->frGEN - 1;
    b = _estSeq + nextExon->frEST - 1;

    nextExon->numMatches = 0;
    nextExon->numNs      = 0;
    nextExon->numInDel   = 0;
    nextExon->numEdits   = 0;

    tmp_script = left;  

    //  These are used during SUBSTITUTE below to tell if the base at
    //  a (b) is N (upper or lower case).
    //
    bool  an = false;
    bool  bn = false;

    while (tmp_script) {
      switch (tmp_script->op_type) {
      case  DELETE:
        nextExon->numInDel += tmp_script->num;
        nextExon->numEdits += tmp_script->num;
        a                  += tmp_script->num;
        break;
      case  INSERT:
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

          an = (*a == 'N') || (*a == 'n');
          bn = (*b == 'N') || (*b == 'n');
            

          if (an && bn) {
            //  Both are N.  It isn't a match and it isn't an edit.
            //
            nextExon->numNs++;
          } else if (an || bn) {
            //  One is an N.  Someone has low quality sequence, and we
            //  should penalize.  We need to special case this because
            //  IUPACidentity[][] claims N matches all.
            //
            nextExon->numEdits++;
          } else if (IUPACidentity[(int)*a][(int)*b]) {
            //  Got a match.
            nextExon->numMatches++;
          } else {
            //  Got a substitution
            nextExon->numEdits++;
          }
        }
        break;
      }
      tmp_script = tmp_script->next;
    }

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

    if ((diff=thisExon->frEST-nextExon->toEST-1)!=0 && (diff != _estLen)) {
      enew = (edit_script_list *)ckalloc(sizeof(edit_script_list));
      enew->next_script = *Aligns;
      *Aligns = enew;
      (*Aligns)->offset1 = thisExon->frGEN;
      (*Aligns)->offset2 = thisExon->frEST;
      (*Aligns)->len1 = end1-(*Aligns)->offset1+1;
      (*Aligns)->len2 = end2-(*Aligns)->offset2+1;
      (*Aligns)->script = head;
      (*Aligns)->score = ali_dist;
      
    } else if (diff != _estLen) {

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
