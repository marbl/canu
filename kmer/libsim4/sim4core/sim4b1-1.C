#include "sim4.H"


int
Sim4::SIM4_block1(Exon*  &Lblock,
                  Exon*  &tmp_block,
                  Exon*  &tmp_block1) {
  int rollbflag = 0;
  int best_u = resolve_overlap(tmp_block,tmp_block1,seq1);

  tmp_block1->frGEN += best_u + 1 - tmp_block1->frEST;
  tmp_block1->frEST  = best_u + 1;

  if (((tmp_block1->toEST - tmp_block1->frEST + 1) < 8) ||
      ((tmp_block1->toGEN - tmp_block1->frGEN + 1) < 8)) { 
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

  if (((tmp_block->toEST - tmp_block->frEST + 1) < 8) || 
      ((tmp_block->toGEN - tmp_block->frGEN + 1) < 8)) {
    /* remove exon defined by tmp_block */
    Exon *prev = find_previous(Lblock,tmp_block);
    assert (prev!=NULL);
    prev->next_exon = tmp_block->next_exon;
    prev->flag = tmp_block->flag; 
    if ((tmp_block->toEST - tmp_block->frEST + 1) > 0)
      rollbflag = 1;
    free(tmp_block);
    tmp_block = prev;
  } 

  if (tmp_block->toGEN)
    tmp_block->length = tmp_block->toEST - tmp_block->frEST + 1;
  if (tmp_block1 && tmp_block1->toGEN)
    tmp_block1->length = tmp_block1->toEST - tmp_block1->frEST + 1;

  return(rollbflag);
}
