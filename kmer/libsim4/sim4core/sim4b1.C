#include "sim4.H"

#define  DEFAULT_W         12
#define  DEFAULT_X         12    
#define  DEFAULT_K         16
#define  DEFAULT_C         12


void
Sim4::adjustBoundariesOfMarginalExons(Exon *Lblock) {
  coords *sig;
  char   tmp[50];
  Exon   *newthing;
  Exon   *tmp_block = Lblock->next_exon;

    /* condition for non-signal */
    if (tmp_block && tmp_block->toGEN &&
        (strncmp((char *)(seq1+tmp_block->frGEN-3), END_SIG, (size_t)2) || 
         (tmp_block->frEST!=1))) {
      sig = (G_score>=abs(C_score)) ? &last_AG : &last_AC; 
      if (sig->pos1 && (sig->pos2<=20)) {
        /* generated in extend_bw */
        assert(sig->pos2 > 1);
        (void)strcpy((char *)tmp,END_SIG);
        (void)strncpy((char *)(tmp+2),(char *)seq2,(size_t)sig->pos2-1);
        (void)strcpy((char *)(tmp+sig->pos2+1), START_SIG);
        newthing = bmatch(seq1,tmp,tmp_block->frGEN-3,sig->pos2+3,1,1); 
        if (newthing) {
          Lblock->next_exon->frGEN = sig->pos1;
          Lblock->next_exon->frEST = sig->pos2;
          Lblock->next_exon->length -= sig->pos2-1;
          newthing->next_exon = Lblock->next_exon;
          newthing->ori = (G_score>=abs(C_score)) ? 'G' : 'C';
          Lblock->next_exon = newthing; 
        }
      }
    }

    while (tmp_block && tmp_block->next_exon && tmp_block->next_exon->toGEN) 
      tmp_block = tmp_block->next_exon;
    if (tmp_block && tmp_block->toGEN &&
        (strncmp((char *)(seq1+tmp_block->toGEN),START_SIG,(size_t)2) || (tmp_block->toEST!=len2))) {   
      sig = (G_score>=abs(C_score)) ? &last_GT : &last_CT;
      if (sig->pos1 && (len2-sig->pos2<=20)) {
        assert(len2-sig->pos2 >= 0);
        (void)strcpy((char *)tmp,END_SIG);
        (void)strncpy((char *)(tmp+2),(char *)(seq2+sig->pos2),
                      (size_t)len2-sig->pos2);
        (void)strcpy((char *)(tmp+len2-sig->pos2+2),START_SIG);
        newthing = fmatch(seq1+sig->pos1-1,tmp,
                          len1-sig->pos1+1,len2-sig->pos2+4,
                          sig->pos1-1,sig->pos2+1);
        if (newthing) {
          tmp_block->toGEN = sig->pos1;
          tmp_block->toEST = sig->pos2;
          newthing->next_exon = tmp_block->next_exon;
          tmp_block->next_exon = newthing;
          tmp_block->ori = (G_score>=abs(C_score)) ? 'G' : 'C';
        }
      }
    } 
  }




void
Sim4::findLastAGandAC(Exon *tmp_block1) {
  int v;

  for (v=tmp_block1->frGEN-1; v<=tmp_block1->toGEN-3; v++) 
    if (!strncmp((char *)(seq1+v-2),"AG",(size_t)2)) {
      last_AG.pos1 = v+1;
      last_AG.pos2 = tmp_block1->frEST + (v-tmp_block1->frGEN)+1; 
      break;
    } 

  for (v=tmp_block1->frGEN-1; v<=tmp_block1->toGEN-3; v++)
    if (!strncmp((char *)(seq1+v-2),"AC",(size_t)2)) {
      last_AC.pos1 = v+1;
      last_AC.pos2 = tmp_block1->frEST + (v-tmp_block1->frGEN)+1;
      break;
    }
}




void
Sim4::findLastGTandCT(Exon *tmp_block) {
  int v;

  for (v=tmp_block->toGEN; v>=tmp_block->frGEN; v--)
    if (!strncmp((char *)(seq1+v),"GT",(size_t)2)) {
      last_GT.pos1 = v;
      last_GT.pos2 = tmp_block->toEST-(tmp_block->toGEN-v);
      break;
    }

  for (v=tmp_block->toGEN; v>=tmp_block->frGEN; v--)
    if (!strncmp((char *)(seq1+v),"CT",(size_t)2)) {
      last_CT.pos1 = v;
      last_CT.pos2 = tmp_block->toEST-(tmp_block->toGEN-v);
      break;
    }
}


inline
void
printExons(char *label, Exon *l) {
  fprintf(stdout, label);
  while (l) {
    fprintf(stdout, "%8d %8d %8d %8d %d\n",
            l->frGEN,
            l->toGEN,
            l->frEST,
            l->toEST,
            l->flag);
    
    l = l->next_exon;
  }
  fprintf(stdout, "----------------------------------------\n");
}

#define PRINTEXONS(S, L)
//#define PRINTEXONS(S, L)   printExons(S, L)


/* seq1 = genomic  DNA (text); seq2 = cDNA */
struct edit_script_list *
Sim4::SIM4(char *in_seq1, char *in_seq2, 
           int in_len1, int in_len2,
           int *dist_ptr, Exon **Exons,
           int *pA, int *pT,
           sim4_stats_t *st) {
  int    cflag, diff, cost, rollbflag;
  int    u, v, I, J;
  bool   good_match;
  Exon   *Lblock, *Rblock=NULL, *tmp_block, *last, *prev, *tmp_block1, *tmp_Lblock=NULL, *tmp_Rblock=NULL;

  struct edit_script_list *Script_head;

  seq1 = in_seq1;
  seq2 = in_seq2;
  len1 = in_len1;
  len2 = in_len2;

#if ABORT_EXPENSIVE
  st->tooManyMSPs = false;
  _mspManager.setLength(len2);
#endif

  //fprintf(stderr, "maskTail=%d percent=%f\n", globalParams->_ignorePolyTails, globalParams->_polyTailPercent);

  /* Compute the distance between two sequences A and B */

  *dist_ptr = 0;

  exon_cores(seq1-1, seq2-1, len1, len2, 1, 1, 0, wordSize, mspThreshold1, PERM);

  //  See if there are too many MSPs found.  If so, fail.
  //
#if ABORT_EXPENSIVE
  if (_mspManager.tooManyMSPs()) {
    st->tooManyMSPs     = true;
    st->numberOfMatches = _mspManager.numberOfMSPs();
    return(0L);
  }
#endif

  PRINTEXONS("0\n", exon_list);

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
    if (tmp_block->next_exon==NULL)
      Rblock = tmp_block;
    tmp_block = tmp_block->next_exon;
  }

  if (Lblock && 
      ((Lblock->frGEN>50000 && Lblock->frEST>100) || 
       ((len1-Rblock->toGEN>50000) && (len2-Rblock->toEST>100)))) {
    free_list(exon_list);

    exon_list = _mspManager.link(globalParams->_relinkWeight,
                                 DEFAULT_DRANGE,
                                 1,
                                 1,
                                 0,
                                 true,
                                 seq1, seq2);

    PRINTEXONS("0a\n", exon_list);

    tmp_block = Lblock = exon_list;
    while (tmp_block) {
      if (tmp_block->next_exon==NULL)
        Rblock = tmp_block;
      tmp_block = tmp_block->next_exon;
    }
  }
  _mspManager.clear();

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
    if (tmp_block->next_exon==NULL) 
      Rblock = tmp_block;
    tmp_block = tmp_block->next_exon; 
  }

  PRINTEXONS("0b\n", exon_list);

  /* enclose the current path in the (0,0,0,0) and (M+1,N+1,0,0) brackets */

  Lblock = new_exon (0,0,0,0,0,0,0,Lblock);
  if (Rblock == NULL)
    Rblock = Lblock;
  Rblock->next_exon = new_exon (len1+1,len2+1,0,0,0,0,0,NULL); 

  PRINTEXONS("0c\n", Lblock);

  /* compute current statistics */
  good_match = get_match_quality(Lblock, Rblock, st, len2);


  //  OK here
  PRINTEXONS("INIT\n", Lblock);



  tmp_block = Lblock;
  while ((tmp_block1 = tmp_block->next_exon)!=NULL) {

    PRINTEXONS("LOOP\n", Lblock);

    rollbflag = 0;
    diff = (int)(tmp_block1->frEST-tmp_block->toEST-1);
    if (diff) {
      if (diff<0) {
        int best_u;
        
        best_u = resolve_overlap(tmp_block,tmp_block1,seq1);

        tmp_block1->frGEN += best_u+1-tmp_block1->frEST;
        tmp_block1->frEST = best_u+1;
        if (((u=tmp_block1->toEST-tmp_block1->frEST+1)<=0) || (u<8) ||
            ((v=tmp_block1->toGEN-tmp_block1->frGEN+1)<=0) || (v<8)) { 
          /* remove exon associated with tmp_block1 */ 
          tmp_block->next_exon = tmp_block1->next_exon;
          tmp_block->flag = tmp_block1->flag; 
          rollbflag = 1;
          free(tmp_block1);
          //tmp_block1 = tmp_block; /* not necessary, just to keep it 'clean'*/
          tmp_block1 = NULL;
        }
        
        tmp_block->toGEN -= tmp_block->toEST-best_u; 
        tmp_block->toEST = best_u;
        if (((u=tmp_block->toEST-tmp_block->frEST+1)<=0) || (u<8) || 
            ((v=tmp_block->toGEN-tmp_block->frGEN+1)<=0) || (v<8)) {
          
          /* remove exon defined by tmp_block */
          prev = find_previous(Lblock,tmp_block);
          assert (prev!=NULL);
          prev->next_exon = tmp_block->next_exon;
          prev->flag = tmp_block->flag; 
          if (u>0) rollbflag = 1;
          free(tmp_block);
          tmp_block = prev;
        } 

        if (tmp_block->toGEN)
          tmp_block->length = tmp_block->toEST-tmp_block->frEST+1;
        if (tmp_block1 && tmp_block1->toGEN)
          tmp_block1->length = tmp_block1->toEST-tmp_block1->frEST+1;

      } else {
        /* bridge the gap */
        cflag = (tmp_block1->toEST && tmp_block->toEST) ? 0 : 1;
        if (diff && (tmp_block1->frGEN-tmp_block->toGEN-1>0)) {
          if (!cflag) {
            if (diff <= MAX_GRINIT) {
              cost = greedy(seq2+tmp_block->toEST,
                            seq1+tmp_block->toGEN, 
                            diff,
                            tmp_block1->frGEN-tmp_block->toGEN-1,
                            tmp_block->toEST,tmp_block->toGEN,
                            &tmp_Lblock, &tmp_Rblock);
            } else cost = max(wordSize,(int)(P*diff+1))+1;

            PRINTEXONS("greedy\n", tmp_Lblock);

            if (cost>max(wordSize,(int)(P*diff+1))) {
              if (!tmp_block->flag && !tmp_block1->flag) {
                exon_cores(seq1+tmp_block->toGEN-1,
                           seq2+tmp_block->toEST-1,
                           tmp_block1->frGEN-tmp_block->toGEN-1,
                           diff,
                           tmp_block->toGEN+1,
                           tmp_block->toEST+1,
                           1,
                           min(8,wordSize),
                           mspThreshold2,
                           TEMP);

                PRINTEXONS("1\n", exon_list);

                tmp_Lblock = tmp_Rblock = exon_list;
                while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL)) 
                  tmp_Rblock = tmp_Rblock->next_exon;

                if ((!tmp_Lblock && tmp_block1->frGEN-tmp_block->toGEN>50000) ||
                    (tmp_Lblock && (tmp_Lblock->frEST-tmp_block->toEST>100) && 
                     (tmp_Lblock->frGEN-tmp_block->frGEN>50000)) || 
                    (tmp_Lblock && (tmp_block1->frEST-tmp_Rblock->toEST>100) &&
                     (tmp_block1->frGEN-tmp_Rblock->frGEN>50000))) {
                  /* possible large intron; increase the score weight */
                  free_list(tmp_Lblock); 

                  exon_list = _mspManager.link(globalParams->_relinkWeight,
                                               DEFAULT_DRANGE, 
                                               tmp_block->toGEN + 1,
                                               tmp_block->toEST + 1,
                                               1,
                                               true,
                                               seq1, seq2); 
                  
                  PRINTEXONS("1a\n", exon_list);

                  tmp_Lblock = tmp_Rblock = exon_list;
                  while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                    tmp_Rblock = tmp_Rblock->next_exon;
                }    
                _mspManager.clear();

                if (tmp_Lblock) rollbflag = 1;
                else rollbflag = 0;   /* already 0 */
              } else 
                tmp_Lblock = tmp_Rblock = NULL;
            }
          } else if (tmp_block1->toGEN) {

            /* start of seq; find last_AG, last_AC */
            if (_accurateSequences)
              findLastAGandAC(tmp_block1);

            diff = (int)(min(diff,(int)(MAX_GRINIT/2)));
            u = min(4*diff,tmp_block1->frGEN-tmp_block->toGEN-1); 
            cost = EXTEND_BW(seq2+tmp_block->toEST+
                             (tmp_block1->frEST-tmp_block->toEST-1)-diff,
                             seq1+tmp_block->toGEN+
                             (tmp_block1->frGEN-tmp_block->toGEN-1)-u,
                             (int)diff, u,
                             tmp_block->toEST+
                             (tmp_block1->frEST-tmp_block->toEST-1)-diff,
                             tmp_block->toGEN+
                             (tmp_block1->frGEN-tmp_block->toGEN-1)-u,
                             &I, &J);
            if ((good_match==FALSE) || tmp_block->flag || (J==0) || (I==0)) {
              tmp_block1->frEST = I+1;
              tmp_block1->frGEN = J+1;
              tmp_block1->edist += cost;
              tmp_block1->length = tmp_block1->toEST-tmp_block1->frEST+1;
            }

            /* use blast if marginal gap still exists, and this is first scan */
            if (!(diff=(int)(tmp_block1->frEST-tmp_block->toEST-1)) ||
                tmp_block->flag) {
              /* blast-treated region or no gap */
              tmp_Rblock = tmp_Lblock = NULL;
            } else {
              exon_cores(seq1+tmp_block->toGEN-1,
                         seq2+tmp_block->toEST-1,
                         tmp_block1->frGEN-tmp_block->toGEN-1,
                         diff,
                         tmp_block->toGEN+1,
                         tmp_block->toEST+1,
                         1,
                         min(10,wordSize),
                         mspThreshold2,
                         TEMP);

              PRINTEXONS("2\n", exon_list);

              tmp_block -> flag = 1;
              tmp_Lblock = tmp_Rblock = exon_list;
              while (tmp_Rblock && tmp_Rblock->next_exon)
                tmp_Rblock = tmp_Rblock->next_exon;

              if ((!tmp_Lblock && tmp_block1->frGEN-tmp_block->toGEN>50000) ||
                  (tmp_Lblock && (tmp_Lblock->frEST-tmp_block->toEST>100) &&
                   (tmp_Lblock->frGEN-tmp_block->frGEN>50000)) ||
                  (tmp_Lblock && (tmp_block1->frEST-tmp_Rblock->toEST>100) &&
                   (tmp_block1->frGEN-tmp_Rblock->frGEN>50000))) {
                /* possible large intron; increase the score weight */
                free_list(tmp_Lblock);

                exon_list = _mspManager.link(globalParams->_relinkWeight,
                                             DEFAULT_DRANGE,
                                             tmp_block->toGEN + 1,
                                             tmp_block->toEST + 1,
                                             1,
                                             true,
                                             seq1, seq2);

                PRINTEXONS("2a\n", exon_list);

                tmp_Lblock = tmp_Rblock = exon_list;
                while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                  tmp_Rblock = tmp_Rblock->next_exon;
              }
              _mspManager.clear();

              if (tmp_Lblock) rollbflag = 1;
              else {
                tmp_block1->frEST = I+1;
                tmp_block1->frGEN = J+1;
                tmp_block1->edist += cost;
                tmp_block1->length = tmp_block1->toEST-tmp_block1->frEST+1;
              }
            }

          } else {

            if (_accurateSequences)
              findLastGTandCT(tmp_block);

            PRINTEXONS("tmp_block before extend\n", tmp_block);
            diff = (int)(min(diff,(int)(MAX_GRINIT/2)));
            cost = EXTEND_FW(seq2+tmp_block->toEST,
                             seq1+tmp_block->toGEN,
                             diff,
                             min(4*diff,tmp_block1->frGEN-tmp_block->toGEN-1),
                             tmp_block->toEST,tmp_block->toGEN,
                             &I, &J);
            PRINTEXONS("tmp_block after extend\n", tmp_block);
            if ((good_match==FALSE) || tmp_block1->flag || (I==len1) || (J==len2)) {
              if (tmp_block->toGEN) {
                tmp_block->toEST = I;
                tmp_block->toGEN = J;
                tmp_block->edist += cost;
                tmp_block->length = tmp_block->toEST-tmp_block->frEST+1;
                tmp_Rblock = tmp_Lblock = NULL;
              } else 
                /* special case: no initial exon */
                tmp_Lblock = tmp_Rblock = NULL;
            }
            PRINTEXONS("tmp_block after if\n", tmp_block);
            /* use blast if marginal gap still exists, and this is first scan */
            if (!(diff=(int)(tmp_block1->frEST-tmp_block->toEST-1)) ||
                tmp_block1->flag) {
              /* blast-treated region or no gap */
              tmp_Rblock = tmp_Lblock = NULL;
            } else {
              PRINTEXONS("tmp_block\n", tmp_block);
              PRINTEXONS("tmp_block1\n", tmp_block);
              exon_cores(seq1+tmp_block->toGEN-1,
                         seq2+tmp_block->toEST-1,
                         tmp_block1->frGEN-tmp_block->toGEN-1,
                         diff,
                         tmp_block->toGEN+1,
                         tmp_block->toEST+1,
                         1,
                         min(10,wordSize),
                         mspThreshold2,
                         TEMP);

              PRINTEXONS("3\n", exon_list);

              tmp_Lblock = tmp_Rblock = exon_list;
              while (tmp_Rblock && tmp_Rblock->next_exon)
                tmp_Rblock = tmp_Rblock->next_exon;

              if ((!tmp_Lblock && tmp_block1->frGEN-tmp_block->toGEN>50000) ||
                  (tmp_Lblock && (tmp_Lblock->frEST-tmp_block->toEST>100) &&
                   (tmp_Lblock->frGEN-tmp_block->frGEN>50000)) ||
                  (tmp_Lblock && (tmp_block1->frEST-tmp_Rblock->toEST>100) &&
                   (tmp_block1->frGEN-tmp_Rblock->frGEN>50000))) {
                /* possible large intron; increase the score weight */
                free_list(tmp_Lblock);

                exon_list = _mspManager.link(globalParams->_relinkWeight,
                                             DEFAULT_DRANGE, 
                                             tmp_block->toGEN + 1,
                                             tmp_block->toEST + 1,
                                             1,
                                             true,
                                             seq1, seq2);

                PRINTEXONS("3a\n", exon_list);

                tmp_Lblock = tmp_Rblock = exon_list;
                while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
                  tmp_Rblock = tmp_Rblock->next_exon;
              }
              _mspManager.clear();

              tmp_block1->flag = 1;
              if (tmp_Lblock) rollbflag = 1;
              else {
                if (tmp_block->toGEN) {
                  tmp_block->toEST = I;
                  tmp_block->toGEN = J;
                  tmp_block->edist += cost;
                  tmp_block->length = tmp_block->toEST-tmp_block->frEST+1;
                  tmp_Rblock = tmp_Lblock = NULL;
                } else
                  /* special case: no initial exon */
                  tmp_Lblock = tmp_Rblock = NULL;
              }
            }
          } 
        } else if (diff) {
          tmp_Rblock = tmp_Lblock = NULL;
        }
        
        /* merge block in the exon list; make connections 
           to the previous list of blocks; maintain 
           increasing order */

        if (tmp_Lblock) {       
          tmp_block->next_exon = tmp_Lblock;
          tmp_Rblock->next_exon = tmp_block1;

          PRINTEXONS("BEFORE MERGE tmp_block\n",  tmp_block);
          PRINTEXONS("BEFORE MERGE tmp_block1\n", tmp_block1);
          PRINTEXONS("BEFORE MERGE tmp_Lblock\n", tmp_Lblock);
          PRINTEXONS("BEFORE MERGE tmp_Rblock\n", tmp_Rblock);

          merge(&tmp_block,&tmp_block1);
        }
      }
    }  /* diff!=0 */
    if (!rollbflag) tmp_block = tmp_block1;
  }


  PRINTEXONS("FINAL\n", Lblock);


  /* compaction step; note: it resets the right end of the list to   */ 
  /* the last item in the block list                                 */

  compact_list(&(Lblock->next_exon), &Rblock);


  /* eliminate marginal small blocks at the start of the sequence;   */
  /* resets the empty alignment to one block (Lblock) only           */

  tmp_block = Lblock->next_exon;

  while ((tmp_block!=NULL) && (tmp_block->length<wordSize) && tmp_block->toGEN) {
    tmp_block1 = tmp_block; /* free */
    tmp_block = tmp_block->next_exon;
    free(tmp_block1); /* free */
  }
  Lblock->next_exon = tmp_block;

  /* eliminate marginal small blocks at the end of the sequence      */

  last = Lblock->next_exon;
  tmp_block = last;
  while (tmp_block!=NULL) {
    if (tmp_block->length>=wordSize)
      last = tmp_block;
    tmp_block = tmp_block->next_exon;
  }
  if (last && last->toGEN)
    last->next_exon = Rblock->next_exon;
  Rblock = last;

  /* if high accuracy requirement, adjust boundaries of marginal exons */
  if (_accurateSequences)
    adjustBoundariesOfMarginalExons(Lblock);

  /* Slide exon boundaries for optimal intron signals */
  if (get_sync_flag(Lblock, Rblock, 6) == TRUE)
    sync_slide_intron(6,&Lblock,st);
  else
    slide_intron(6,&Lblock,st);

  //fprintf(stdout, "init strand orientation = %d\n", st->orientation);

  /* decreasingly; script will be in reverse order */
  flip_list(&Lblock, &Rblock); 
  pluri_align(dist_ptr, Lblock, &Script_head, st); 
  flip_list(&Lblock, &Rblock);      /* increasingly */

  *pT = 0;
  *pA = 0;
  if (Script_head) {
    if (globalParams->_ignorePolyTails) {
      remove_polyT_front(&Script_head, Lblock, seq1, seq2, pT); 
      remove_polyA_back(&Script_head, Lblock, seq1, seq2, len2, pA);

      if (*pA || *pT)
        updateStatistics(Lblock, st);
    }

    //fprintf(stdout, "pA = %d\npT = %d\n", *pA, *pT);

    get_stats(Lblock, st);

    //fprintf(stdout, "final strand orientation = %d\n", st->orientation);

    *Exons = Lblock->next_exon;
  }
  free(Lblock);

  //  Memory leak when Script_head == 0L -- see pluri_align, too!

  return(Script_head);
}
