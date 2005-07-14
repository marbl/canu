#include "sim4.H"


//  Original call was if (!strncmp(S, "GT", 2)) {}
//  which is if (S == "GT")
//
#define DAcmp(S, A, B) (((S)[0] == A) && ((S)[1] == B))


void
Sim4::complement_exons(Exon **left, int M, int N) { 
  Exon  *tmp_block, *right;
  char   prev, ch;

#ifdef SPLSCORE
  double spl=0, prevspl=0;
#endif

  prev = 'U'; /* unknown, should trigger error */
  tmp_block = *left;
  while (tmp_block) {
    if (tmp_block->toGEN) {
      register int aux;
      
      if (tmp_block->next_exon && tmp_block->next_exon->toGEN) {
        ch = tmp_block->ori;
        tmp_block->ori = prev;

#ifdef SPLSCORE
        spl = tmp_block->splScore;
        tmp_block->splScore = prevspl;
        prevspl = spl;
#endif

        switch (ch) {
          case 'C':
            prev = 'G';
            break;
          case 'G':
            prev = 'C';
            break;
          case 'N':
            prev = 'N';
            break;
          case 'E':
            prev = 'E';
            break;
          default:
            fatal("sim4b1.c: Inconsistency. Check exon orientation at complementation.");
        }
      } else {
        tmp_block->ori = prev;
#ifdef SPLSCORE
        tmp_block->splScore = prevspl;
#endif
      }
      aux = tmp_block->frGEN;
      tmp_block->frGEN = M+1-tmp_block->toGEN;
      tmp_block->toGEN = M+1-aux;
      aux = tmp_block->frEST;
      tmp_block->frEST = N+1-tmp_block->toEST;
      tmp_block->toEST = N+1-aux;
    }
    tmp_block = tmp_block->next_exon;
    if (tmp_block && tmp_block->toGEN) 
      right = tmp_block;
  }
  flip_list(left,&right);
}


#if 0

//  Dead code, 05 apr 2004, bpw

void
Sim4::check_consistency_intron_ori(Exon *exons, int match_ori, char *gene)
{
  Exon *t=exons;
  int numG, numC, numE, numN;
  
  numG = numC = numE = numN = 0;

  if (!t->toGEN) t = t->next_exon;
  while (t && t->toGEN) {
    if (t->next_exon  && t->next_exon->toGEN) {
      switch (t->ori) {
        case 'G':
          numG++;
          break;
        case 'C':
          numC++;
          break;
        case 'N':
          numN++;
          break;
        case 'E':
          numE++;
          break;
        default:
          fatal("sim4b1.c: Unrecognized intron orientation.");
      } 
    }
    t = t->next_exon;
  }
  if (numG && numC) 
    fprintf(stderr, "Warning: Introns reported on both strands (%s).\n", gene);
  else if ((numG && (match_ori==BWD)) || (numC && (match_ori==FWD))) 
    fprintf(stderr, "Warning: Introns orientations inconsistent with the reported match strand (%s).\n", gene);
  
  if (numN) 
    fprintf(stderr, "Warning: Ambiguous intron orientation (%s).\n", gene);
  if (numE) 
    fprintf(stderr, "Warning: Internal gap in the mRNA (%s).\n", gene);

  return;
}

#endif


void
ABORT(int X, char *Y) {
  if (X) {
    fprintf(stderr, "ERROR: %s\n", Y);
    kill(getpid(), 10);
  }
}


void
Sim4::get_stats(Exon *lblock, sim4_stats_t *st) {
  Exon *t, *t1;

  st->icoverage = 0;
  st->internal  = 1;

  if ((lblock->next_exon == NULL) || !lblock->next_exon->toGEN)                            
    st->internal = 0;

  for (t=lblock->next_exon; t; t = t->next_exon)
    st->icoverage += t->length;

  t = lblock;
  while (t) {
    t1 = t->next_exon;

    if ((t->toGEN) &&
        (t1) && 
        (t1->frEST - t->toEST - 1 > 0) &&
        t1->toGEN)
      st->internal = 0;

    t = t1;
  }

  if (!globalParams->_forceStrandPrediction) {
    if (((st->orientation != BOTH) && (st->percentID < 90)) ||
        (st->internal == 0)) {
      st->orientation = BOTH;
    }
  }
}





void
Sim4::flip_list(Exon **left, Exon **right)
{
  Exon   *ep, *ahead, *behind;

  *right = *left;
  ahead = *left;
  ep = NULL;
  while (ahead!=NULL) {
    behind = ep;
    ep = ahead;
    ahead = ahead->next_exon;
    ep->next_exon = behind;
  }
  *left = ep;
}



/* operates on a list sorted in increasing order of exon coordinates */
void
Sim4::compact_list(Exon **Lblock, Exon **Rblock)
{  
  Exon *tmp_block=*Lblock, *tmp_block1;
  int diff;

  while ((tmp_block!=NULL) && 
         ((tmp_block1=tmp_block->next_exon)!=NULL) &&
         tmp_block1->toGEN) {
    if ((abs((tmp_block1->frEST-tmp_block1->frGEN) -
             (tmp_block->toEST-tmp_block->toGEN))<=wordSize) &&
        ((diff=tmp_block1->frEST-tmp_block->toEST-1)<=MAX_INTERNAL_GAP)) {
      /* merge blocks */
      tmp_block->toGEN = tmp_block1->toGEN;
      tmp_block->toEST = tmp_block1->toEST;
      tmp_block->length = tmp_block->toEST-tmp_block->frEST+1;
      tmp_block->edist += tmp_block1->edist;
      tmp_block->edist -= (int)(globalParams->_percentError * diff);
      tmp_block->next_exon = tmp_block1->next_exon;
      
      //freeExon(tmp_block1);  garbage collected
    } else
      tmp_block = tmp_block1;
  }
  /* reset right end of the list */
  *Rblock = tmp_block;
}

/* ------------------  memory management routines --------------- */


int
Sim4::good_ratio(int length)
{
  if (length<=wordSize/2) return 2;
  else if (length<2*wordSize) return DIST_CUTOFF;
  else return (int)(.75 * globalParams->_percentError * length + 1);
}


void
Sim4::merge(Exon **t0, Exon **t1) 
{   
  Exon  *tmp0, *tmp1;
  int    diff;

  if ((*t0) && !(*t0)->toGEN)
    tmp0 = (*t0)->next_exon;
  else    
    tmp0 = *t0;

  while (tmp0 && (tmp0!=*t1)) {
    tmp1 = tmp0->next_exon;
    assert(tmp1!=NULL);

    if (tmp1 && tmp1->toGEN && tmp0->toGEN && 
        (abs((tmp1->frEST-tmp1->frGEN)-(tmp0->toEST-tmp0->toGEN))<=wordSize) &&
        ((tmp1->frEST - tmp0->toEST - 1 <= wordSize))) {

      diff = tmp1->frEST - tmp0->toEST - 1;
      
      /* merge blocks tmp0 and tmp1 */
      tmp0->frGEN = min(tmp0->frGEN, tmp1->frGEN);
      tmp0->frEST = min(tmp0->frEST, tmp1->frEST);
      tmp0->toGEN = max(tmp1->toGEN, tmp0->toGEN);
      tmp0->toEST = max(tmp1->toEST, tmp0->toEST);
      tmp0->length = tmp0->toEST-tmp0->frEST+1;
      tmp0->flag = tmp1->flag;
      tmp0->edist += tmp1->edist;
      tmp0->edist -= (int)(globalParams->_percentError * diff);
      if (tmp1==*t1) {
        /*  tmp0->flag = (*t1)->flag; */
        *t1 = tmp0;
      }
      tmp0->next_exon = tmp1->next_exon;

      //freeExon(tmp1);  garbage collected
    } else {
      tmp0 = tmp0->next_exon;
    }
  }
}

void
Sim4::free_align(edit_script_list *aligns) {
  edit_script_list *head;

  head = aligns;

  while ((head=aligns)!=NULL) {
    aligns = aligns->next_script;
    Free_script(head->script);
    ckfree(head);
  }
}



Exon *
Sim4::bmatch (char *s1, char *s2, int l1, int l2, int offset1, int offset2)
{
  int  i, j, i1, score;
  Exon *newthing=NULL;
  
  for (i1=i=l1-3; i>=l2-3; i--, i1=i) {
    for (j=l2-3; j>=2; j--, i1--)
      if (*(s1+i1)!=*(s2+j))
        break;
    
    if (j<2) {
      /* exact match for CDS found; check signals */
      score = 0;
      if (*(s1+(i1--))==*(s2+(j--))) score++;
      if (*(s1+(i1--))==*(s2+(j--))) score++;
      if (*(s1+i1+l2-1)==*(s2+j+l2-1)) score++;
      if (*(s1+i1+l2)==*(s2+j+l2)) score++;
      if (score>=3) {
        newthing = _exonManager.newExon(i1+3+offset1, offset2, i1+3+offset1+l2-5,
                                        offset2+l2-5, l2-4, 0, 0, NULL);
        newthing->ori = (G_score >= abs(C_score)) ? 'G' : 'C';
        
        return newthing;
      }
    }
  }
  return NULL;
}

Exon *
Sim4::fmatch (char *s1, char *s2, int l1, int l2, int offset1, int offset2)
{
  int  i, j, i1, score;
  Exon *newthing=NULL;
  
  for (i1=i=2; i<l1-l2+3; i++, i1=i) {
    for (j=2; j<l2-2; j++, i1++)
      if (*(s1+i1)!=*(s2+j))
        break;

    if (j>=l2-2) {
      /* exact match found for internal part, look for signals */
      score = 0;
      if (*(s1+(i1++))==*(s2+(j++))) score++;
      if (*(s1+(i1++))==*(s2+(j++))) score++;
      if (*(s1+i1-l2)==*s2) score++;
      if (*(s1+i1-l2+1)==*(s2+1)) score++;
      if (score>=3) {
        newthing = _exonManager.newExon(i+offset1,offset2,i1+offset1-2,offset2+l2-5,
                                        l2-4,0,0,NULL);
        newthing->ori = (G_score >= abs(C_score)) ? 'G' : 'C';

        return newthing;
      }
    }  
  }    
  return NULL;
}


/* -------------------- to be added to psublast ---------------------- */

bool 
Sim4::get_sync_flag(Exon *lblock, Exon *rblock, int w)
{
  int numx=0, e2;
  Exon *t;

  if (((t=lblock->next_exon)==NULL) || !t->toGEN)
    return 0;
  numx++;
  e2 = t->toEST;

  while (((t=t->next_exon)!=NULL) && t->toGEN) {
    ++numx;
    if ((t->frEST-e2>1) || 
        (t!=rblock && ((t->toEST-t->frEST+1<2*w+2) || (t->toGEN-t->frGEN+1<2*w+2))))
      return 0;
    e2 = t->toEST;
  }
  return ((numx<3) ? 0:1);
}



void 
Sim4::sync_slide_intron(int in_w, Exon **lblock, sim4_stats_t *st) {
  Exon *t0, *t1, *head = *lblock;
  splice_t *g, *c, *cell;
  int w1, w2, i;

  //  Count the exons to allocate space for Glist, Clist and oris
  //
  int ni   = 1;
  t0 = head->next_exon;
  while (t0 && (t1=t0->next_exon) && t1->toGEN) {
    ni++;
    t0 = t1;
  }        

  //  These leak, but I don't know where.
  splice_t **Glist = (splice_t **)ckcalloc(ni * sizeof(splice_t *));
  splice_t **Clist = (splice_t **)ckcalloc(ni * sizeof(splice_t *));
  char      *oris  = (char      *)ckcalloc(ni * sizeof(char));

  if ((Glist == 0L) || (Clist == 0L) || (oris == 0L))
    fprintf(stderr, "Can't allocate memory for sync_slide_intron() with %d exons.\n", ni), exit(1);

  int numG = 0, Gscore = 0;
  int numC = 0, Cscore = 0;

  ni = 0;

  /* assume forward orientation */
  t0 = head->next_exon;
  while (t0 && (t1=t0->next_exon) && t1->toGEN) {
    g = c = NULL;
    if (t1->frEST-t0->toEST-1==0) {

      if (DAcmp(_genSeq+t0->toGEN, 'G', 'T') &&
          DAcmp(_genSeq+t1->frGEN-3, 'A', 'G')) {
        g = new_splice('G',t0->toGEN,t1->frGEN,t0->toEST,t1->frEST,-1,NULL);
        t0->ori = 'G';
        oris[ni] = 'G';
        numG++;
#ifdef SPLSCORE
        t0->splScore = 999999;
#endif
      } else if (DAcmp(_genSeq+t0->toGEN, 'C', 'T') &&
                 DAcmp(_genSeq+t1->frGEN-3, 'A', 'C')) {
        c = new_splice('C',t0->toGEN,t1->frGEN,t0->toEST,t1->frEST,-1,NULL);
        t0->ori = 'C';
        oris[ni] = 'C';
        numC++;
#ifdef SPLSCORE
        t0->splScore = 888888;
#endif
      } else {
        w1 = min(in_w, min(t0->length-1, t0->toGEN-t0->frGEN));
        w2 = min(in_w, min(t1->length-1, t1->toGEN-t1->frGEN));
        splice(_genSeq, t0->toGEN-w1, t0->toGEN+w1, t1->frGEN-w2, t1->frGEN+w2,
               _estSeq, t0->toEST-w1, t1->frEST+w2, &g, &c, BOTH);

        Gscore += g->score;
        Cscore += c->score;
        cell = NULL;
        oris[ni] = '*';
        if (g->score>c->score) {
          numG++;
          cell = g;
          oris[ni] = 'G';
        } else if (c->score>g->score) {
          numC++;
          cell = c;
          oris[ni] = 'C';
        } else if (c->score==g->score) {
          numG++;
          numC++;
          cell = g;
          oris[ni] = 'G';
        }
#ifdef SPLSCORE
        t0->splScore = cell->score;
#endif
        t0->ori = oris[ni];
        t0->toGEN = cell->xs;
        t0->toEST = cell->ys;
        t1->frGEN = cell->xe;
        t1->frEST = cell->ye;
        t0->length = t0->toEST-t0->frEST+1;
        t1->length = t1->toEST-t1->frEST+1;
      }
      Clist[ni] = c;
      Glist[ni] = g;
    } else {
      t0->ori = 'E';
      oris[ni] = 'E';
    }
    ni++;
    t0 = t1;
  }         
  
  st->orientation = BOTH;

  if ((numG==1) && (numC==1) && 
      (!Glist[0] || !Clist[0] || !Glist[1] || !Clist[1])) goto free_all;
  
  if (numG>=numC) {
    /* revisit all previous assignments that are inconsistent */
    for (i=0, t0=head->next_exon; i<ni; i++, t0=t1) {
      t1 = t0->next_exon;
      switch (oris[i]) {
        case 'G':
          break;
        case 'C': if (Glist[i]==NULL) { 
          /* compute the values for C */ 
          w1 = min(in_w, min(t0->length-1, t0->toGEN-t0->frGEN));
          w2 = min(in_w, min(t1->length-1, t1->toGEN-t1->frGEN));
          splice(_genSeq, t0->toGEN-w1, t0->toGEN+w1,
                 t1->frGEN-w2, t1->frGEN+w2, _estSeq,
                 t0->toEST-w1, t1->frEST+w2, &g, &c, FWD);
        } else {
          g = Glist[i];
        }
          
#ifdef SPLSCORE
          t0->splScore = g->score;
#endif

          t0->ori = 'G';
          t0->toGEN = g->xs;
          t0->toEST = g->ys;
          t1->frGEN = g->xe;
          t1->frEST = g->ye;
          t0->length = t0->toEST-t0->frEST+1;
          t1->length = t1->toEST-t1->frEST+1;

          break;
        case 'E': break;
        default : fatal("sim4b1.c: intron orientation not initialized.");
      }
      if (oris[i]!='E')
        wobble(t0, t1, "GT", "AG", _genSeq);
    }

    st->orientation = FWD;
  } else {
    /* analyze all assignments for consistency */
    for (i=0, t0=head->next_exon; i<ni; i++, t0=t1) {
      t1 = t0->next_exon;
      switch (oris[i]) {
        case 'C': break;
        case 'G': if (Clist[i]==NULL) {
          /* compute the values for C */
          w1 = min(in_w, min(t0->length-1, t0->toGEN-t0->frGEN));
          w2 = min(in_w, min(t1->length-1, t1->toGEN-t1->frGEN));
          splice(_genSeq, t0->toGEN-w1, t0->toGEN+w1,
                 t1->frGEN-w2, t1->frGEN+w2,
                 _estSeq, t0->toEST-w1, t1->frEST+w2, &g, &c, BWD);
        } else c = Clist[i];
                
#ifdef SPLSCORE
          t0->splScore = c->score;
#endif
          t0->ori = 'C';
          t0->toGEN = c->xs;
          t0->toEST = c->ys;
          t1->frGEN = c->xe;
          t1->frEST = c->ye;
          t0->length = t0->toEST-t0->frEST+1;
          t1->length = t1->toEST-t1->frEST+1;
          break;
        case 'E': break;
        default : fatal("sim4b1.c: intron orientation not initialized.");
      }
      if (oris[i]!='E')
        wobble(t0, t1, "CT", "AC", _genSeq);
    }

    st->orientation = BWD;
  }

  /* now free all memory allocated */
 free_all:
  for (i=0; i<ni; i++) {
    ckfree(Glist[i]);
    ckfree(Clist[i]);
  }

  ckfree(Glist);
  ckfree(Clist);
  ckfree(oris);
}

void 
Sim4::wobble(Exon *t0,
             Exon *t1,
             const char *donor,
             const char *acceptor,
             char *seq) {
  char *s = seq + t0->toGEN;      // first nt of donor
  char *q = seq + t1->frGEN - 3;  // first nt of acceptor

  if (DAcmp(s, donor[0], donor[1])) {
    /* match in place */
    if (DAcmp(q, acceptor[0], acceptor[1])) {
      return;
    } else if (DAcmp(q-1, acceptor[0], acceptor[1])) {
      t1->frGEN--;
      return;
    } else if (DAcmp(q+1, acceptor[0], acceptor[1])) {
      t1->frGEN++;
      return;
    }
  } else if (DAcmp(s-1, donor[0], donor[1])) {
    /* match is 1 off to the left */
    if (DAcmp(q, acceptor[0], acceptor[1])) {   
      t0->toGEN--;
      return;
    } else if (DAcmp(q-1, acceptor[0], acceptor[1])) {
      t0->toGEN--;
      t1->frGEN--;
      t0->toEST--;
      t1->frEST--;
      t0->length++;
      t1->length--;
      return;
    } else if (DAcmp(q+1, acceptor[0], acceptor[1])) {
      t0->toGEN--;
      t1->frGEN++;
      return;
    }
  } else if (DAcmp(s+1, donor[0], donor[1])) {
    /* match is 1 off to the right */
    if (DAcmp(q, acceptor[0], acceptor[1])) {
      t0->toGEN++;
      return;
    } else if (DAcmp(q-1, acceptor[0], acceptor[1])) {
      t0->toGEN++;
      t1->frGEN--;
      return;
    } else if (DAcmp(q+1, acceptor[0], acceptor[1])) {
      t0->toGEN++;
      t1->frGEN++;
      t0->toEST++;
      t1->frEST++;
      t0->length--;
      t1->length++;
      return;
    }
  } else if (DAcmp(q-1, acceptor[0], acceptor[1])) {
    /* match is 1 off to the left */
    t1->frGEN--;
    return;
  } else if (DAcmp(q+1, acceptor[0], acceptor[1])) {
    /* match is 1 off to the right */
    t1->frGEN++;
    return;
  }
}

void
Sim4::slide_intron(int in_w, Exon **lblock, sim4_stats_t *st) { 
  Exon *t0, *t1, *head = *lblock;
  splice_t *cell;
  char type;
  int numG=0, numC=0, numE=0, numN=0;

  t0 = head->next_exon;
  while (t0 && (t1=t0->next_exon) && t1->toGEN) {
    //fprintf(stderr, "Sim4::slide_intron()-- loop  t0=%p t1=%p\n", t0, t1);

    //  Don't do anything if the exon is short.  It bombs.
    //
    if (((t0->toGEN - t0->frGEN) <= 1) ||
        ((t1->toGEN - t1->frGEN) <= 1))
      fprintf(stderr, "Sim4::slide_intron()--  WARNING:  Short genomic range in exon!\n");

    splice_t *g = 0L;
    splice_t *c = 0L;

    if ((t1->frEST - t0->toEST - 1) == 0) {
      if (DAcmp(_genSeq+t0->toGEN, 'G', 'T') && 
          DAcmp(_genSeq+t1->frGEN-3, 'A', 'G')) {
        t0->ori = 'G';
        numG++;
#ifdef SPLSCORE
        t0->splScore = 999999;
#endif
      } else if (DAcmp(_genSeq+t0->toGEN, 'C', 'T') && 
                 DAcmp(_genSeq+t1->frGEN-3, 'A', 'C')) {
        t0->ori = 'C';
        numC++;
#ifdef SPLSCORE
        t0->splScore = 888888;
#endif
      } else {      
        int gtag=0, ctac=0;
        char *s;
        
        int w1 = min(in_w, min(t0->length - 2, t0->toGEN - t0->frGEN - 1));
        int w2 = min(in_w, min(t1->length - 2, t1->toGEN - t1->frGEN - 1));
        if (w1 < 0) w1 = 0;
        if (w2 < 0) w2 = 0;
        splice(_genSeq, t0->toGEN-w1, t0->toGEN+w1, t1->frGEN-w2, t1->frGEN+w2,
               _estSeq, t0->toEST-w1, t1->frEST+w2, &g, &c, BOTH);

        if (g->score>c->score) {
          cell = g;
          type = 'G';
        } else if (c->score>g->score) {
          cell = c;
          type = 'C';
        } else {
          cell = g;
          type = 'G';
        }
        
#ifdef SPLSCORE
        t0->splScore = cell->score;
#endif

        t0->toGEN  = cell->xs;
        t0->toEST  = cell->ys;
        t1->frGEN  = cell->xe;
        t1->frEST  = cell->ye;
        t0->length = t0->toEST - t0->frEST + 1;
        t1->length = t1->toEST - t1->frEST + 1;
        
        //wobble(&t0,&t1,(type=='G')? "GT":"CT",(type=='G')? "AG":"AC",_genSeq);
#if 0
        fprintf(stderr, "Sim4::slide_intron()-- before wobble: t0=%p t1=%p\n", t0, t1);
        fprintf(stderr, "Sim4::slide_intron()-- t0: GEN:%d-%d EST:%d-%d len:%d\n", t0->frGEN, t0->toGEN, t0->frEST, t0->toEST, t0->length);
        fprintf(stderr, "Sim4::slide_intron()-- t1: GEN:%d-%d EST:%d-%d len:%d\n", t1->frGEN, t1->toGEN, t1->frEST, t1->toEST, t1->length);
#endif

        if (type == 'G')
          wobble(t0, t1, "GT", "AG", _genSeq);
        else
          wobble(t0, t1, "CT", "AC", _genSeq);

#if 0
        fprintf(stderr, "Sim4::slide_intron()-- after wobble: t0=%p t1=%p\n", t0, t1);
        fprintf(stderr, "Sim4::slide_intron()-- t0: GEN:%d-%d EST:%d-%d len:%d\n", t0->frGEN, t0->toGEN, t0->frEST, t0->toEST, t0->length);
        fprintf(stderr, "Sim4::slide_intron()-- t1: GEN:%d-%d EST:%d-%d len:%d\n", t1->frGEN, t1->toGEN, t1->frEST, t1->toEST, t1->length);
#endif

        ckfree(g);
        ckfree(c);

        /* determine the type, based on the # matches w/ GT-AG (CT-AC) */

        s = _genSeq+t0->toGEN;
        if (*s=='G')
          gtag++;
        if (*s=='C')
          ctac++;
        s++;

        if (*s=='T') {
          gtag++;
          ctac++;
        } 
        s = _genSeq+t1->frGEN-3;

        if (*s=='A') {
          gtag++;
          ctac++;
        }
        ++s;

        if (*s=='G')
          gtag++;
        if (*s=='C')
          ctac++;

        if (gtag>ctac) {
          type = 'G';
          numG++;
        } else if (ctac>gtag) {
          type = 'C';
          numC++;
        } else {
          type = 'N';
          numN++;
        }

        t0->ori = type;
      }
    } else {
      t0->ori = 'E';
      numE++;
    }
    t0 = t1;
  }      
  
  st->orientation = BOTH;

  if ((numG > 0) && ((numC + numE + numN) == 0)) {
    st->orientation = FWD;
  } else if ((numC > 0) && ((numG + numE + numN) == 0)) {
    st->orientation = BWD;
  }

  if ((globalParams->_forceStrandPrediction) && (st->orientation == BOTH)) {
    if (numG > numC)
      st->orientation = FWD;
    if (numG < numC)
      st->orientation = BWD;

    //  otherwise, st->orientation = match orientation, but we
    //  don't know that here.  It's set in sim4string.C:run()
  }
}



bool
Sim4::get_match_quality(Exon *lblock, Exon *rblock, sim4_stats_t *st, int N)
{
  int  tcov;
  bool good_match;
  Exon *t;
  
  good_match = 1;
  st->icoverage = 0;
  t = lblock->next_exon;
  while (t->toGEN) {
    st->icoverage += t->toEST-t->frEST+1;
    if (100*t->edist>=5*(t->toEST-t->frEST+1)) {
      good_match = 0;
      break;
    }
    t = t->next_exon;
  }
  tcov = rblock->toEST-lblock->next_exon->frEST+1;
  if (lblock->next_exon->frEST>=.5*N &&
      tcov>=.75*(N-lblock->next_exon->frEST) &&
      st->icoverage>=max(.95*tcov,100))
    ;
  else if (rblock->toEST<=.5*N && tcov>=.75*rblock->toEST &&
           st->icoverage>=max(.95*tcov,100))
    ;
  else if ((tcov<.7*N) ||
           (st->icoverage<.9*tcov))
    good_match = 0;
  
  return good_match;
}


