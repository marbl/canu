#include "sim4.H"


int
Sim4::resolve_overlap(Exon *tmp_block, Exon *tmp_block1, char *seq) {          
  int   diff, best_u, l0, l1, u, cost;
  int    GTAG_score, CTAC_score;
  char *s1, *s2, *e1;
  
  diff = tmp_block1->frEST-tmp_block->toEST-1;
  if (diff>=0) return (tmp_block1->frEST-1);

  /* resolve overlap using the GT-AG criterion */
  /* u-1 = actual position in the sequence */
  
  l0 = tmp_block->length-diff;
  l1 = tmp_block1->length; 
  
  best_u = u = tmp_block1->frEST-1;
  s1 = seq+tmp_block->toGEN-(tmp_block->toEST-u);
  s2 = seq-2+tmp_block1->frGEN+u-tmp_block1->frEST;
  
  cost = 0;
  e1 = seq+tmp_block->toGEN;
  while (s1<=e1) {
    GTAG_score = CTAC_score = 0;
    GTAG_score += ((char)(*s1)=='G') ? 1 : 0;
    GTAG_score += ((char)(*(s1+1))=='T') ? 1 : 0;
    GTAG_score += ((char)(*s2)=='A') ? 1 : 0;
    GTAG_score += ((char)(*(s2+1))=='G') ? 1 : 0;
    
    if (GTAG_score > abs(cost) && ((l0>=8) || (l1>=8))) {
      cost = GTAG_score;
      best_u = u;
      if (cost == 4) break; 
    }
    
    CTAC_score += ((char)(*s1)=='C') ? 1 : 0;
    CTAC_score += ((char)(*(s1+1))=='T') ? 1 : 0;
    CTAC_score += ((char)(*s2)=='A') ? 1 : 0;
    CTAC_score += ((char)(*(s2+1))=='C') ? 1 : 0;
    
    if (CTAC_score > abs(cost)) {
      cost = -CTAC_score;
      best_u = u;
      if (cost == 4) break;
    }
    
    u++; s1++; s2++;
    l0++; l1--;
  }     
  
  return best_u;
}      


int
Sim4::SIM4_block1(Exon*  &Lblock,
                  Exon*  &tmp_block,
                  Exon*  &tmp_block1) {
  int rollbflag = 0;
  int best_u = resolve_overlap(tmp_block,tmp_block1,_genSeq);

  tmp_block1->frGEN += best_u + 1 - tmp_block1->frEST;
  tmp_block1->frEST  = best_u + 1;

  //fprintf(stderr, "sim4_block1()-- Lblock=%p tmp_block=%p tmp_block1=%p\n", Lblock, tmp_block, tmp_block1);

  if (((tmp_block1->toEST - tmp_block1->frEST + 1) < 8) ||
      ((tmp_block1->toGEN - tmp_block1->frGEN + 1) < 8)) { 

    // remove exon associated with tmp_block1

    tmp_block->next_exon = tmp_block1->next_exon;
    tmp_block->flag = tmp_block1->flag; 
    rollbflag = 1;
    ckfree(tmp_block1);
    tmp_block1 = NULL;
  }

  tmp_block->toGEN -= tmp_block->toEST-best_u; 
  tmp_block->toEST = best_u;

  if (((tmp_block->toEST - tmp_block->frEST + 1) < 8) || 
      ((tmp_block->toGEN - tmp_block->frGEN + 1) < 8)) {

    // remove exon defined by tmp_block

    Exon *prev = find_previous(Lblock, tmp_block);

    if (prev == 0L) {
      fprintf(stderr, "SIM4_block1(): Corrupted exon list, cannot find the previous exon.\n");
      for (; Lblock; Lblock = Lblock->next_exon)
        if (tmp_block == Lblock)
          fprintf(stderr, "  GEN f=%8d t=%8d  EST f=%8d t=%8d   flag=%d <- tried to find previous of this one\n",
                  Lblock->frGEN, Lblock->toGEN, Lblock->frEST, Lblock->toEST, Lblock->flag);
        else
          fprintf(stderr, "  GEN f=%8d t=%8d  EST f=%8d t=%8d   flag=%d\n",
                  Lblock->frGEN, Lblock->toGEN, Lblock->frEST, Lblock->toEST, Lblock->flag);
      kill(getpid(), SIGKILL);
    }

    prev->next_exon = tmp_block->next_exon;
    prev->flag = tmp_block->flag; 
    if ((tmp_block->toEST - tmp_block->frEST + 1) > 0)
      rollbflag = 1;
    ckfree(tmp_block);
    tmp_block = prev;
  } 

  if (tmp_block->toGEN)
    tmp_block->length = tmp_block->toEST - tmp_block->frEST + 1;
  if (tmp_block1 && tmp_block1->toGEN)
    tmp_block1->length = tmp_block1->toEST - tmp_block1->frEST + 1;

  return(rollbflag);
}
