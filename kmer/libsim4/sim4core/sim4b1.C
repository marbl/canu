#include "sim4.H"

//  SHOW_PROGRESS -- write the progress of Sim4::SIM4 to stderr
//  DEBUG_EXONS   -- dump the exons at various places
//
//#define SHOW_PROGRESS
//#define DEBUG_EXONS


#ifdef DEBUG_EXONS
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
  fflush(stdout);
}

#define PRINTEXONS(S, L)   printExons(S, L)
#else
#define PRINTEXONS(S, L)
#endif



/* seq1 = genomic  DNA (text); seq2 = cDNA */
struct edit_script_list *
Sim4::SIM4(int            *dist_ptr,
           Exon          **Exons,
           int            *pA,
           int            *pT,
           sim4_stats_t   *st) {

  int    rollbflag;
  Exon   *Lblock=0L, *tmp_Lblock=0L;
  Exon   *Rblock=0L, *tmp_Rblock=0L;
  Exon   *tmp_block=0L;
  Exon   *tmp_block1=0L;

  struct edit_script_list *Script_head;

  *dist_ptr = 0;
  *Exons    = 0L;
  *pA = 0;
  *pT = 0;

  //  Initialize the mspManager to fail if the match looks expensive.
  //
  st->tooManyMSPs = false;
  _mspManager.setLength(_estLen);


#ifdef SHOW_PROGRESS
  fprintf(stderr, "before exon_cores\n");
  double beforeExonCoresStartTime = getTime();
#endif

  exon_cores(_genSeq-1, _estSeq-1, _genLen, _estLen, 1, 1, 0, wordSize, mspThreshold1, PERM);

#ifdef SHOW_PROGRESS
  fprintf(stderr, "after exon_cores -- took %f seconds.\n", getTime() - beforeExonCoresStartTime);
#endif

  //  See if there are too many MSPs found.  If so, fail.
  //
  if (_mspManager.tooManyMSPs()) {
    st->tooManyMSPs     = true;
    st->numberOfMatches = _mspManager.numberOfMSPs();
    return(0L);
  }

  PRINTEXONS("0\n", exon_list);

  tmp_block = Lblock = exon_list;
  while (tmp_block) {
    if (tmp_block->next_exon==NULL)
      Rblock = tmp_block;
    tmp_block = tmp_block->next_exon;
  }

  if (Lblock && 
      ((Lblock->frGEN>50000 && Lblock->frEST>100) || 
       ((_genLen - Rblock->toGEN > 50000) && (_estLen - Rblock->toEST > 100)))) {
    free_list(exon_list);

    exon_list = _mspManager.doLinking(globalParams->_relinkWeight,
                                      DEFAULT_DRANGE,
                                      1,
                                      1,
                                      0,
                                      true,
                                      _genSeq, _estSeq);

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

#ifdef SHOW_PROGRESS
  fprintf(stderr, "exon bracket at start\n");
#endif
  Lblock = new_exon(0,0,0,0,0,0,0,Lblock);
  if (Rblock == NULL)
    Rblock = Lblock;
#ifdef SHOW_PROGRESS
  fprintf(stderr, "exon bracket at end; Lblock = 0x%08lx, Rblock = 0x%08lx\n", Lblock, Rblock);
#endif
  Rblock->next_exon = new_exon(_genLen+1,_estLen+1,0,0,0,0,0,NULL); 

  PRINTEXONS("0c\n", Lblock);

  /* compute current statistics */
  bool good_match = get_match_quality(Lblock, Rblock, st, _estLen);


  PRINTEXONS("INIT\n", Lblock);


#ifdef SHOW_PROGRESS
  fprintf(stderr, "before big nasty while loop\n");
#endif


  tmp_block = Lblock;
  while ((tmp_block1 = tmp_block->next_exon)!=NULL) {

    PRINTEXONS("LOOP\n", Lblock);

    rollbflag = 0;

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
#ifdef SHOW_PROGRESS
        fprintf(stderr, "Called SIM4_block1()\n");
#endif
        rollbflag = SIM4_block1(Lblock, tmp_block, tmp_block1);
      } else {
        /* bridge the gap */

        if (tmp_block1->frGEN - tmp_block->toGEN - 1 > 0) {
          if (tmp_block1->toEST &&
              tmp_block->toEST) {
#ifdef SHOW_PROGRESS
            fprintf(stderr, "Called SIM4_block2()\n");
#endif
            rollbflag = SIM4_block2(tmp_Lblock,
                                    tmp_Rblock,
                                    tmp_block,
                                    tmp_block1);
          } else if (tmp_block1->toGEN) {
#ifdef SHOW_PROGRESS
            fprintf(stderr, "Called SIM4_block3()\n");
#endif
            rollbflag = SIM4_block3(good_match,
                                    tmp_Lblock,
                                    tmp_Rblock,
                                    tmp_block,
                                    tmp_block1);
          } else {
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
    }
    if (rollbflag == 0)
      tmp_block = tmp_block1;
  }


  PRINTEXONS("FINAL\n", Lblock);

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
    ckfree(tmp_block1);
  }
  Lblock->next_exon = tmp_block;


  /* eliminate marginal small blocks at the end of the sequence      */

  //
  //  XXX:  Memory leak here??!?!
  //

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

  /* if high accuracy requirement, adjust boundaries of marginal exons */
  if (_accurateSequences)
    adjustBoundariesOfMarginalExons(Lblock);

  /* Slide exon boundaries for optimal intron signals */
  if (get_sync_flag(Lblock, Rblock, 6) == 1)
    sync_slide_intron(6,&Lblock,st);
  else
    slide_intron(6,&Lblock,st);

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
    ckfree(Lblock);
  } else {
    *Exons = 0L;

    free_list(Lblock);
  }

  //  Memory leak when Script_head == 0L -- see pluri_align, too!

  return(Script_head);
}
