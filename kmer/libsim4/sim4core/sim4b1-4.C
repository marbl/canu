#include "sim4.H"

int
Sim4::SIM4_block4(bool     good_match,
                  Exon*   &tmp_Lblock,
                  Exon*   &tmp_Rblock,
                  Exon*   &tmp_block,
                  Exon*   &tmp_block1) {
  int  I, J;
  int rollbflag = 0;
  int  cost;

  fprintf(stderr, "Called SIM4_block4()\n");

  if (_accurateSequences)
    findLastGTandCT(tmp_block);

  //PRINTEXONS("tmp_block before extend\n", tmp_block);

  int diff = (int)(tmp_block1->frEST - tmp_block->toEST - 1);
  diff = (int)(min(diff,(int)(MAX_GRINIT/2)));

  cost = EXTEND_FW(seq2+tmp_block->toEST,
                   seq1+tmp_block->toGEN,
                   diff,
                   min(4*diff,tmp_block1->frGEN-tmp_block->toGEN-1),
                   tmp_block->toEST,tmp_block->toGEN,
                   &I, &J);
  //PRINTEXONS("tmp_block after extend\n", tmp_block);
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
  //PRINTEXONS("tmp_block after if\n", tmp_block);
  /* use blast if marginal gap still exists, and this is first scan */
  if (!(diff=(int)(tmp_block1->frEST-tmp_block->toEST-1)) ||
      tmp_block1->flag) {
    /* blast-treated region or no gap */
    tmp_Rblock = tmp_Lblock = NULL;
  } else {
    //PRINTEXONS("tmp_block\n", tmp_block);
    //PRINTEXONS("tmp_block1\n", tmp_block);
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

    //PRINTEXONS("3\n", exon_list);

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

      //PRINTEXONS("3a\n", exon_list);

      tmp_Lblock = tmp_Rblock = exon_list;
      while ((tmp_Rblock!=NULL) && (tmp_Rblock->next_exon!=NULL))
        tmp_Rblock = tmp_Rblock->next_exon;
    }
    _mspManager.clear();

    tmp_block1->flag = 1;
    if (tmp_Lblock) {
      rollbflag = 1;
    } else {
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
  return(rollbflag);
}
