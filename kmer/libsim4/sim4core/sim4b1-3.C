#include "sim4.H"

int
Sim4::SIM4_block3(bool     good_match,
                  Exon*   &tmp_Lblock,
                  Exon*   &tmp_Rblock,
                  Exon*   &tmp_block,
                  Exon*   &tmp_block1) {
  int  I, J;
  int rollbflag = 0;
  int  cost;

  //fprintf(stderr, "Called SIM4_block3()\n");

  /* start of seq; find last_AG, last_AC */
  if (_accurateSequences)
    findLastAGandAC(tmp_block1);

  //  These two blocks should do the same thing.  The first one isn't readable.

#if 0
  int diff = (int)(tmp_block1->frEST - tmp_block->toEST - 1);
  diff = (int)(min(diff,(int)(MAX_GRINIT/2)));

  int u = min(4*diff,tmp_block1->frGEN-tmp_block->toGEN-1); 
  cost = EXTEND_BW(_estSeq+tmp_block->toEST+
                   (tmp_block1->frEST-tmp_block->toEST-1)-diff,
                   _genSeq+tmp_block->toGEN+
                   (tmp_block1->frGEN-tmp_block->toGEN-1)-u,
                   (int)diff, u,
                   tmp_block->toEST+
                   (tmp_block1->frEST-tmp_block->toEST-1)-diff,
                   tmp_block->toGEN+
                   (tmp_block1->frGEN-tmp_block->toGEN-1)-u,
                   &I, &J);
#else
  int diff = min(tmp_block1->frEST - tmp_block->toEST - 1, MAX_GRINIT/2);
  int u    = min(4*diff, tmp_block1->frGEN - tmp_block->toGEN - 1); 

  cost = EXTEND_BW(_estSeq + tmp_block1->frEST - 1 - diff,
                   _genSeq + tmp_block1->frGEN - 1 - u,
                   diff,
                   u,
                   tmp_block1->frEST - 1 - diff,
                   tmp_block1->frGEN - 1 - u,
                   &I,
                   &J);
#endif


  if ((good_match==0) || tmp_block->flag || (J==0) || (I==0)) {
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
    exon_cores(_genSeq+tmp_block->toGEN-1,
               _estSeq+tmp_block->toEST-1,
               tmp_block1->frGEN-tmp_block->toGEN-1,
               diff,
               tmp_block->toGEN+1,
               tmp_block->toEST+1,
               1,
               spacedSeedExtMSS,
               mspThreshold2,
               TEMP);

    //PRINTEXONS("2\n", exon_list);

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
      //freeExonList(tmp_Lblock);  garbage collected

      exon_list = _mspManager.doLinking(globalParams->_relinkWeight,
                                        DEFAULT_DRANGE,
                                        tmp_block->toGEN + 1,
                                        tmp_block->toEST + 1,
                                        1,
                                        true,
                                        _genSeq, _estSeq);

      //PRINTEXONS("2a\n", exon_list);

      tmp_Lblock = tmp_Rblock = exon_list;
      while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
        tmp_Rblock = tmp_Rblock->next_exon;
    }
    _mspManager.clear();

    if (tmp_Lblock) {
      rollbflag = 1;
    } else {
      tmp_block1->frEST = I+1;
      tmp_block1->frGEN = J+1;
      tmp_block1->edist += cost;
      tmp_block1->length = tmp_block1->toEST-tmp_block1->frEST+1;
    }
  }
  return(rollbflag);
}
