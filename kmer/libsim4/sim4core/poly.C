#include <math.h>
#include "sim4.H"

#define MIN_EXON 12

void
Sim4::get_polyAT(char *seq, int len, int *pT, int *pA, int flag)
{
  register int i, sum10, sum20;
  register char *s, *t, *v;
  int last10;

  int MAX10 = 2;
  int MAX20 = 5;

  char encodingA[128];
  char encodingT[128];


  if (flag!=T_ONLY) {
    memset(encodingA, (char)1, 128);
    encodingA[(int)'A'] = encodingA[(int)'X'] = encodingA[(int)'N'] = 0;

    for (i=0, s=seq+len, sum10=0, last10=len+1; i<10 && s>seq && sum10<=MAX20; i++)  {
      sum10 += encodingA[(int)*(--s)];
      /*          if (!encodingA[*s] && sum10<=MAX10) last10 = s-seq+1; */
    }

    t = v = seq+len;
    sum20 = sum10;
    for ( ; s>=seq && (sum10<=MAX10 || sum20<=MAX20); ) {
      if (!encodingA[(int)*s] && sum10<=MAX10 && (seq+len>=s+20 || sum20<MAX10))
        last10 = (int)(s-seq+1);
      if (--s>seq) {
        sum10 += encodingA[(int)*s] - encodingA[(int)*(--t)];
        sum20 += encodingA[(int)*s] -(((seq+len)-s>20) ? encodingA[(int)*(--v)] : 0);
      }
    }

    if (last10>len-10) *pA = len+1;
    else {
      s = seq+last10+8;
      while (s >= seq && !encodingA[(int)*s]) s--;
      if ((s-seq+1)-last10+1<=5)
        *pA = (int)(s-seq+2);
      else
        *pA = last10;
    }
  } else *pA = len+1;
  *pA = len-(*pA)+1;

  if (flag!=A_ONLY) {

    memset(encodingT, (char)1, 128);
    encodingT[(int)'T'] = encodingT[(int)'X'] = encodingT[(int)'N'] = 0;

    for (i=0, s=seq-1, sum10=0, last10=0; i<10 && i<len-1 && sum10<=MAX20; i++) {
      sum10 += encodingT[(int)*(++s)];
      /*          if (!encodingT[*s] && sum10<=MAX10) last10 = s-seq+1; */
    }

    t = v = seq-1;
    sum20 = sum10;
    for ( ; s<seq+len && (sum10<=MAX10 || sum20<=MAX20); ) {
      if (!encodingT[(int)*s] && sum10<=MAX10 && (s-seq>=19 || sum20<MAX10))
        last10 = (int)(s-seq+1);
      if (++s<seq+len) {
        sum10 += encodingT[(int)*s] - encodingT[(int)*(++t)];
        sum20 += encodingT[(int)*s] - ((s-seq>=20) ? encodingT[(int)*(++v)] : 0);
      }
    }

    if (last10<=10) *pT = 0;
    else {
      s = seq+last10-10;
      while (s < seq+len && !encodingT[(int)*s]) s++;
      if (last10-(s-seq)+1<=5)
        *pT = (int)(s-seq);
      else
        *pT = last10;
    }
  } else *pT = 0;
}

void
Sim4::trim_polyA_align(struct edit_script_list **Sptr, Exon *lblock, Exon **exons, const int bc, int *pA, char *s1,char *s2) 
{
  edit_script_list *head = *Sptr;
  edit_script *tp;
  int tmpi = 0, num, idents = 0, identsN = 0;
  char *a, *b;
  Exon *prev;

  int i, j;  /* i index in the cDNA */

  if (bc>head->offset2+head->len2-1) {
    *pA = bc;
    return;
  }

  if (bc==head->offset2) {
    /* cDNA gap: remove the entire script; is this properly sorted? LLL */
    *Sptr = head->next_script;
    Free_script(head->script);
    ckfree(head);
    while ((*exons)->frEST>=bc) {
      prev = find_previous(lblock,*exons);

      if (prev == 0L) {
        fprintf(stderr, "trim_polyA_align(): Corrupted exon list, cannot find the previous exon (remove entire script).\n");
        for (; lblock; lblock = lblock->next_exon)
          fprintf(stderr, "  GEN f=%8d t=%8d  EST f=%8d t=%8d   flag=%d\n",
                  lblock->frGEN, lblock->toGEN, lblock->frEST, lblock->toEST, lblock->flag);
        kill(getpid(), SIGKILL);
      }

      prev->next_exon = (*exons)->next_exon;
      //freeExon(*exons);  garbage collected
      *exons = prev;
    }
    *pA = bc;
    return;
  }

  Flip_script(&(head->script));
  i = head->offset2 + head->len2 -1;
  j = head->offset1 + head->len1 -1;
  tp = head->script;

  while (i>=bc && tp) {
    num = tp->num;
    switch (tp->op_type) {
      case INSERT:
        if (i>=bc && bc>i-num+1) {
          (*exons)->numInDel -= i - bc + 1;
          (*exons)->numEdits -= i - bc + 1;
          tmpi    += i-bc+1;
          tp->num -= i-bc+1;
          i        = bc-1;
        } else {
          (*exons)->numInDel -= num;
          (*exons)->numEdits -= num;
          tmpi += num;
          i    -= num;
          head->script = tp->next;
          ckfree(tp);
          tp = head->script;
        }
        break;
      case DELETE:
        (*exons)->numInDel -= num;
        (*exons)->numEdits -= num;
        j    -= num;
        tmpi += num;
        head->script = tp->next;
        ckfree(tp);
        tp = head->script;
        break;
      case SUBSTITUTE:
        if (i>=bc && bc>i-num+1) {
          a = s2+i-1; b = s1+j-1;
          while (a>=s2+bc-1) {
            if (*a != *b) {
              (*exons)->numEdits--;
              tmpi++;
            } else {
              if (*a == 'N') {
                (*exons)->numNs--;
                identsN++;
              } else {
                (*exons)->numMatches--;
                idents++;
              }
            }
            a--;
            b--;
          }
          j -= i-bc+1; tp->num -= i-bc+1; i = bc-1;
        } else {
          /* at most 1 nt remaining */
          a = s2+i-1; b = s1+j-1;
          while (a>=s2+i-num) {
            if (*a != *b) {
              (*exons)->numEdits--;
              tmpi++;
            } else {
              if (*a == 'N') {
                (*exons)->numNs--;
                identsN++;
              } else {
                (*exons)->numMatches--;
                idents++;
              }
            }
            a--;
            b--;
          }

          i -= num; j -= num;
          head->script = tp->next;
          ckfree(tp);
          tp = head->script;
        }
        break;
#if 0
      default:
        fatalf("Unrecognized opcode %d.\n",tp->op_type);
#endif
    }
    /* indel walk */
  }
  assert(i==bc-1);

  while ((tp != 0L) &&
         (tp->op_type != SUBSTITUTE) && (j+1 >= (*exons)->frGEN)) {
    if (tp->op_type==INSERT) {
      i -= tp->num;
      tmpi += tp->num;
      (*exons)->numInDel -= tp->num;
      (*exons)->numEdits -= tp->num;
    } else if (j<(*exons)->frGEN && i<(*exons)->frEST) {
      j -= tp->num;
    } else {
      j -= tp->num;
      tmpi += tp->num;
      (*exons)->numInDel -= tp->num;
      (*exons)->numEdits -= tp->num;
    }
    head->script = tp->next;
    ckfree(tp);
    tp = head->script;
  }

  if (head->script==NULL) {
    *Sptr = head->next_script;
    ckfree(head);
  } else {
    head->len1 = j-head->offset1+1;
    head->len2 = i-head->offset2+1;
    head->score -= tmpi;
    Flip_script(&(head->script));
  }

  if ((*exons)->frEST>i) {
    prev = find_previous(lblock,*exons);

    if (prev == 0L) {
      fprintf(stderr, "trim_polyA_align(): Corrupted exon list, cannot find the previous exon (frEST).\n");
      for (; lblock; lblock = lblock->next_exon)
        fprintf(stderr, "  GEN f=%8d t=%8d  EST f=%8d t=%8d   flag=%d\n",
                lblock->frGEN, lblock->toGEN, lblock->frEST, lblock->toEST, lblock->flag);
      kill(getpid(), SIGKILL);
    }

    prev->next_exon = (*exons)->next_exon;
    //freeExon(*exons);  garbage collected
    *exons = prev;
  } else {
    (*exons)->toEST = i;
    (*exons)->toGEN = j;
    (*exons)->length = (*exons)->toEST-(*exons)->frEST+1;

    (*exons)->alignmentLength = ((*exons)->toGEN - (*exons)->frGEN + 1 +
                                 (*exons)->toEST - (*exons)->frEST + 1 +
                                 (*exons)->numInDel);
    (*exons)->percentID   = computePercentIdentity((*exons)->numEdits,
                                                   (*exons)->alignmentLength);
  }
  *pA = i+1;

  return;
}



void
Sim4::remove_polyA_back(struct edit_script_list **Sptr, Exon *Exons,
                        char *s1, char *s2,
                        int l2, int *lastA) {
  Exon *t;
  Exon *exons_tail;
  char *b, *end;
  int numA, pA, dummy, trim_p, reverse_script=0;
  int startPos=0, cutAmount=0;

  *lastA = l2+1;  pA = 0;
  if (!Exons || ! Exons->next_exon || ! Exons->next_exon->toGEN) return;

  if ((*Sptr)->next_script &&
      (*Sptr)->offset1<(*Sptr)->next_script->offset1) {
    reverse_script = 1;
    script_flip_list(Sptr);
  }

  exons_tail = Exons->next_exon;
  while (exons_tail->next_exon && exons_tail->next_exon->toGEN)
    exons_tail=exons_tail->next_exon;

  trim_p = 1;

  if (exons_tail) {
    startPos = exons_tail->toEST;

    while ((t=exons_tail)!=NULL && t->toGEN && trim_p) {
      /* compute the 'A' contents of the exon */
      b = s2 + t->toEST-1; end = s2+t->frEST-1; numA = 0;
      while (b>=end && numA+(b-end)>=globalParams->_polyTailPercent*t->length) { 
        numA += (*b--=='A'); 
      }

      //  Determine how much of the cut stuff was actually
      //  poly-containing.  The first method below returns the number of
      //  bases cut from the end of the est, while the second return the
      //  number of bases cut from the end of the alignment.
      //
      //cutAmount = l2 - *lastA + 1;

      if (numA>=globalParams->_polyTailPercent*t->length) {
        /* remove the entire exon */
        trim_polyA_align(Sptr,Exons,&exons_tail,t->frEST,lastA,s1,s2);
        cutAmount = startPos - *lastA + 1;
      } else {
        get_polyAT(s2+(*Sptr)->offset2-1,(*Sptr)->len2,&dummy,&pA,A_ONLY);
        if (pA) {
          int ct_pA;
          /* first position to be removed */ 
          ct_pA = t->toEST-pA+1; 
          ct_pA = (ct_pA-t->frEST>MIN_EXON) ? ct_pA : t->frEST;
          /* note: pA is the last (innermost) position in the tail */
          trim_polyA_align(Sptr,Exons,&exons_tail,ct_pA,lastA,s1,s2);
          cutAmount = startPos - *lastA + 1;
        }
        if (t==exons_tail) trim_p = 0;
      }
    }
  }

  *lastA = cutAmount;

  if (reverse_script) script_flip_list(Sptr);
}



/* s2 is the cdna */
void
Sim4::trim_polyT_align(struct edit_script_list **Sptr, Exon **exons, const int ec, int *pT, char *s1, char *s2)
{
  edit_script_list *head = *Sptr;
  edit_script *tp;
  int tmpi = 0, num, idents = 0, identsN = 0;
  char *a, *b;
  Exon *t;

  int i, j;  /* i index in the cDNA */

  if (ec<head->offset2) {
    *pT = ec; 
    return;
  }

  if (ec==head->offset2+head->len2-1) {
    /* cDNA gap: remove the entire script */
    *Sptr = head->next_script;
    Free_script(head->script);
    ckfree(head);
    while ((*exons)->frEST<ec) {
      t = *exons;
      *exons = t->next_exon;
      //freeExon(t);  garbage collected
    }
    *pT = ec;
    return;
  }
 
  i = head->offset2;
  j = head->offset1;
  tp = head->script;

  while (i<=ec && tp) {
    num = tp->num;
    switch (tp->op_type) {
      case INSERT:
        if (i<=ec && ec<i+num-1) {
          (*exons)->numInDel -= ec - i + 1;
          (*exons)->numEdits -= ec - i + 1;
          tmpi    += ec-i+1;
          tp->num -= ec-i+1;
          i        = ec+1; 
        } else {
          (*exons)->numInDel -= num;
          (*exons)->numEdits -= num;
          tmpi += num;
          i    += num;
          head->script = tp->next; 
          ckfree(tp);
          tp = head->script; 
        }
        break;
      case DELETE:     
        (*exons)->numInDel -= num;
        (*exons)->numEdits -= num;
        j    += num;
        tmpi += num;
        head->script = tp->next;
        ckfree(tp);
        tp = head->script; 
        break;
      case SUBSTITUTE:
        if (i<=ec && ec<i+num-1) {
          a = s2+i-1; b = s1+j-1;
          while (a<s2+ec) {
            if (*a != *b) {
              (*exons)->numEdits--;
              tmpi++;
            } else {
              if (*a == 'N') {
                (*exons)->numNs--;
                identsN++;
              } else {
                (*exons)->numMatches--;
                idents++;
              }
            }
            a++;
            b++;
          }
          j += ec-i+1; tp->num -= ec-i+1; i = ec+1;
        } else {
          /* at most 1 nt remaining */
          a = s2+i-1; b = s1+j-1;
          while (a<s2+i+tp->num-1) {
            if (*a != *b) {
              (*exons)->numEdits--;
              tmpi++;
            } else {
              if (*a == 'N') {
                (*exons)->numNs--;
                identsN++;
              } else {
                (*exons)->numMatches--;
                idents++;
              }
            }
            a++;
            b++;
          }

          i +=num; j += num;
          head->script = tp->next;
          ckfree(tp);
          tp = head->script; 
        }
        break;
    }
    /* indel walk */
  }
  assert(i==ec+1);

  while ((tp != 0L) &&
         (tp->op_type!=SUBSTITUTE) && (j-1<=(*exons)->toGEN)) {
    if (tp->op_type==INSERT) {
      i += tp->num;
      tmpi += tp->num;
      (*exons)->numInDel -= tp->num;
      (*exons)->numEdits -= tp->num;
    } else if (j>=(*exons)->toGEN && i>=(*exons)->toEST) {
      j += tp->num;
    } else {
      j += tp->num;
      tmpi += tp->num;
      (*exons)->numInDel -= tp->num;
      (*exons)->numEdits -= tp->num;
    }
    head->script = tp->next;
    ckfree(tp);
    tp = head->script;
  }
    
  if (head->script==NULL) {
    *Sptr = head->next_script;
    ckfree(head);
  } else {
    head->len1 -= j-head->offset1;
    head->len2 -= i-head->offset2;
    head->offset2 = i;
    head->offset1 = j;          
    head->score -= tmpi;
  }

  if ((*exons)->toEST<i) {
    t = *exons;
    *exons = t->next_exon;
    //freeExon(t);  garbage collected
  } else {
    (*exons)->frEST = i; 
    (*exons)->frGEN = j;
    (*exons)->length = (*exons)->toEST-(*exons)->frEST+1;

    (*exons)->alignmentLength = ((*exons)->toGEN - (*exons)->frGEN + 1 +
                                 (*exons)->toEST - (*exons)->frEST + 1 +
                                 (*exons)->numInDel);
    (*exons)->percentID   = computePercentIdentity((*exons)->numEdits,
                                                   (*exons)->alignmentLength);
  }
  *pT = i-1;
  return;
}



void
Sim4::remove_polyT_front(struct edit_script_list **Sptr, Exon *Exons, char *s1, char *s2, int *lastT)
{
  Exon *t, *exons_head; /* start from Lblock */
  char *b, *end;
  int numT, dummy, trim_p, reverse_script=0, pT;
  int startPos=0, cutAmount=0;

  *lastT = pT = 0;
  if (!Exons || !Exons->next_exon || !Exons->next_exon->toGEN) return;
 
  if ((*Sptr)->next_script && 
      (*Sptr)->offset1>(*Sptr)->next_script->offset1) {
    script_flip_list(Sptr);
    reverse_script = 1;
  }

  exons_head = Exons->next_exon; trim_p = 1;

  if (exons_head) {
    startPos = exons_head->frEST;

    while ((t=exons_head)!=NULL && t->toGEN && trim_p) {
      /* compute the 'T' contents of the exon */
      b = s2 + t->frEST-1; end = s2+t->toEST; numT = 0;
      while (b<end && (numT+t->toEST-(b-s2)>=globalParams->_polyTailPercent*t->length)) {
        numT += (*b++=='T');
      }

      //  Determine how much of the cut stuff was actually
      //  poly-containing.  The first method below returns the number of
      //  bases cut from the end of the est, while the second return the
      //  number of bases cut from the end of the alignment.
      //
      //cutAmount = l2 - *lastT + 1;

      if (numT>=globalParams->_polyTailPercent*t->length) {
        /* remove the entire exon */
        trim_polyT_align(Sptr,&exons_head,t->toEST,lastT,s1,s2);
        cutAmount = *lastT - startPos + 1;
      } else {
        get_polyAT(s2+(*Sptr)->offset2-1,(*Sptr)->len2,&pT,&dummy,T_ONLY);
        if (pT) {
          int ct_pT;
          ct_pT = pT + (*Sptr)->offset2-1;
          ct_pT = (t->toEST-ct_pT>MIN_EXON) ? ct_pT : t->toEST;
          trim_polyT_align(Sptr,&exons_head,ct_pT,lastT,s1,s2);
          cutAmount = *lastT - startPos + 1;
        }
        if (t==exons_head) trim_p = 0;
      }
    }
  }

  Exons->next_exon = exons_head;

  *lastT = cutAmount;

  if (reverse_script) script_flip_list(Sptr);
}
