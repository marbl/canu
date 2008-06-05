#include "sim4.H"

#ifdef DEBUG_EXONS
#define PRINTEXONS(S, L)   (L)->printList(S)
#else
#define PRINTEXONS(S, L)
#endif

Sim4::edit_script_list *
Sim4::SIM4(int            *dist_ptr,
           Exon          **Exons,
           int            *pA,
           int            *pT,
           sim4_stats_t   *st) {

  int     rollbflag;
  Exon   *Lblock=0L, *tmp_Lblock=0L;
  Exon   *Rblock=0L, *tmp_Rblock=0L;
  Exon   *tmp_block=0L;
  Exon   *tmp_block1=0L;

  struct edit_script_list *Script_head;

  *dist_ptr = 0;
  *Exons    = 0L;
  *pA = 0;
  *pT = 0;

  //
  //  The call to exon_cores() that used to be here is now done in sim4string.
  //

  //  See if there are too many MSPs found.  If so, fail.
  //
  st->tooManyMSPs = false;
  if (_mspManager.tooManyMSPs()) {
    st->tooManyMSPs     = true;
    st->numberOfMatches = _mspManager.numberOfMSPs();
    return(0L);
  }

  PRINTEXONS("initial exon set\n", exon_list);

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
    if (tmp_block->next_exon==NULL)
      Rblock = tmp_block;
    tmp_block = tmp_block->next_exon;
  }

  if (Lblock && 
      ((Lblock->frGEN>50000 && Lblock->frEST>100) || 
       ((_genLen - Rblock->toGEN > 50000) && (_estLen - Rblock->toEST > 100)))) {
    //freeExonList(exon_list);  garbage collected

    exon_list = _mspManager.doLinking(globalParams->_relinkWeight,
                                      DEFAULT_DRANGE,
                                      1,
                                      1,
                                      0,
                                      true,
                                      _genSeq, _estSeq);

    PRINTEXONS("relink the initial stuff\n", exon_list);

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

  PRINTEXONS("initial exon set after possibly relinking\n", exon_list);

  /* enclose the current path in the (0,0,0,0) and (M+1,N+1,0,0) brackets */

#ifdef SHOW_PROGRESS
  fprintf(stderr, "exon bracket at start\n");
#endif
  Lblock = _exonManager.newExon(0,0,0,0,0,0,0,Lblock);
  if (Rblock == NULL)
    Rblock = Lblock;
#ifdef SHOW_PROGRESS
  fprintf(stderr, "exon bracket at end; Lblock = 0x%08lx, Rblock = 0x%08lx\n", Lblock, Rblock);
#endif
  Rblock->next_exon = _exonManager.newExon(_genLen+1,_estLen+1,0,0,0,0,0,NULL); 

  PRINTEXONS("initial exon set after inserting brackets\n", Lblock);

  /* compute current statistics */
  bool good_match = get_match_quality(Lblock, Rblock, st, _estLen);


  PRINTEXONS("after get_match_quality\n", Lblock);


#ifdef SHOW_PROGRESS
  fprintf(stderr, "before big nasty while loop\n");
#endif


  tmp_block = Lblock;
  while ((tmp_block1 = tmp_block->next_exon)!=NULL) {

    PRINTEXONS("start of loop to fill in missing pieces\n", Lblock);

    rollbflag = 0;

    //  This is the distance from this exon to the next exon
    //  in the EST
    //
    int diff = (int)(tmp_block1->frEST - tmp_block->toEST - 1);

#ifdef SHOW_PROGRESS
    fprintf(stdout, "tmp_block: %8d %8d %8d %8d %d diff=%d\n",
            tmp_block->frGEN,
            tmp_block->toGEN,
            tmp_block->frEST,
            tmp_block->toEST,
            tmp_block->flag,
            diff);
#endif

    if (diff) {

      if (diff < 0) {
        //  If the diff is less than zero, then there is an overlap in
        //  the EST.  Wobble the boundary using GTAG signals (so
        //  obviously, this won't work correctly if we are not cDNA).
        //
#ifdef SHOW_PROGRESS
        fprintf(stderr, "Called SIM4_block1() with diff=%d\n", diff);
#endif
        rollbflag = SIM4_block1(Lblock, tmp_block, tmp_block1);
      } else {

        //  Otherwise, there is a gap in the EST, and we need to fill
        //  it in.  This is done only if there is no overlap in the
        //  genomic.
        //
        if (tmp_block1->frGEN - tmp_block->toGEN - 1 > 0) {
          if (tmp_block1->toEST &&
              tmp_block->toEST) {
            //  We are not the first or last gap -- an interior gap
            //  between two exons.
            //
#ifdef SHOW_PROGRESS
            fprintf(stderr, "Called SIM4_block2()\n");
#endif
            rollbflag = SIM4_block2(tmp_Lblock,
                                    tmp_Rblock,
                                    tmp_block,
                                    tmp_block1);
          } else if (tmp_block1->toGEN) {
            //  Not the last gap, so must be the first gap.
            //
#ifdef SHOW_PROGRESS
            fprintf(stderr, "Called SIM4_block3()\n");
#endif
            rollbflag = SIM4_block3(good_match,
                                    tmp_Lblock,
                                    tmp_Rblock,
                                    tmp_block,
                                    tmp_block1);
          } else {
            //  By default, the last gap.
            //
#ifdef SHOW_PROGRESS
            fprintf(stderr, "Called SIM4_block4()\n");
#endif
            rollbflag = SIM4_block4(good_match,
                                    tmp_Lblock,
                                    tmp_Rblock,
                                    tmp_block,
                                    tmp_block1);
          } 
        } else {
          //  Overlapping genomic.  What these do when set to
          //  NULL is unknown.
          //
          tmp_Rblock = tmp_Lblock = NULL;
        }

        //  Merge block in the exon list; make connections to the
        //  previous list of blocks; maintain increasing order
        //
        if (tmp_Lblock) {       
          tmp_block->next_exon = tmp_Lblock;
          tmp_Rblock->next_exon = tmp_block1;

          PRINTEXONS("before merge tmp_block\n",  tmp_block);
          PRINTEXONS("before merge tmp_block1\n", tmp_block1);
          PRINTEXONS("before merge tmp_Lblock\n", tmp_Lblock);
          PRINTEXONS("before merge tmp_Rblock\n", tmp_Rblock);

          merge(&tmp_block,&tmp_block1);
        }
      }
    }

    //  If this exon block was not removed, move to the next.  If it was removed,
    //  we're already there.
    //
    if (rollbflag == 0)
      tmp_block = tmp_block1;
  }


  PRINTEXONS("all done -- final Lblock\n", Lblock);

#ifdef SHOW_PROGRESS
  fprintf(stderr, "sim4b1 -- before compact_list\n");
#endif

  /* compaction step; note: it resets the right end of the list to   */ 
  /* the last item in the block list                                 */

  compact_list(&(Lblock->next_exon), &Rblock);

#ifdef SHOW_PROGRESS
  fprintf(stderr, "sim4b1 -- before small block at start removal\n");
#endif

  /* eliminate marginal small blocks at the start of the sequence;   */
  /* resets the empty alignment to one block (Lblock) only           */

  tmp_block = Lblock->next_exon;

  while ((tmp_block!=NULL) && (tmp_block->length<wordSize) && tmp_block->toGEN) {
    tmp_block1 = tmp_block;
    tmp_block = tmp_block->next_exon;
    //freeExon(tmp_block1);  garbage collected
  }
  Lblock->next_exon = tmp_block;

  PRINTEXONS("all done -- after removing small blocks at the start\n", Lblock);

  //  eliminate marginal small blocks at the end of the sequence
  //  XXX:  Yes, there is a leak here.  That's why we garbage collect!

#ifdef SHOW_PROGRESS
  fprintf(stderr, "Rblock before end of list removal 0x%08lx\n", Rblock);
#endif

  Exon *last = Lblock->next_exon;
  tmp_block = last;
  while (tmp_block!=NULL) {
    if (tmp_block->length>=wordSize)
      last = tmp_block;
    tmp_block = tmp_block->next_exon;
  }
  if (last && last->toGEN)
    last->next_exon = Rblock->next_exon;
  Rblock = last;

#ifdef SHOW_PROGRESS
  fprintf(stderr, "Rblock after end of list removal 0x%08lx\n", Rblock);
#endif

  PRINTEXONS("all done -- after removing small blocks at the end\n", Lblock);

  /* if high accuracy requirement, adjust boundaries of marginal exons */
  if (_accurateSequences)
    adjustBoundariesOfMarginalExons(Lblock);

  /* Slide exon boundaries for optimal intron signals */
  if (globalParams->_slideIntrons) {
    if (get_sync_flag(Lblock, Rblock, 6) == 1)
      sync_slide_intron(6,&Lblock,st);
    else
      slide_intron(6,&Lblock,st);
  } else {
    //  Set orientation flag on introns to be unknown -- this has an
    //  undesired side effect of forcing the resulting match to have a
    //  strand orientation the same as the intron orientation (if one
    //  exon) instead of 'unknown'.
    Exon *t0 = Lblock->next_exon;
    Exon *t1 = NULL;

    while (t0 && (t1=t0->next_exon) && t1->toGEN) {
      t0->ori = 'E';
      t0 = t1;
    }
  }

  /* decreasingly; script will be in reverse order */
  flip_list(&Lblock, &Rblock); 
  pluri_align(dist_ptr, Lblock, &Script_head, st); 
  flip_list(&Lblock, &Rblock);      /* increasingly */

  *pT = 0;
  *pA = 0;
  if (Script_head) {
    if (globalParams->_ignorePolyTails) {
      remove_polyT_front(&Script_head, Lblock, _genSeq, _estSeq, pT); 
      remove_polyA_back(&Script_head, Lblock, _genSeq, _estSeq, _estLen, pA);

      if (*pA || *pT)
        updateStatistics(Lblock, st);
    }

    get_stats(Lblock, st);

    *Exons = Lblock->next_exon;
    //freeExon(Lblock);  garbage collected
  } else {
    *Exons = 0L;

    //freeExonList(Lblock);  garbage collected
  }

  //  Memory leak when Script_head == 0L -- see pluri_align, too!

  return(Script_head);
}
