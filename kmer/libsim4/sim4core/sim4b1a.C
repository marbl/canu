#include "sim4.H"


void
Sim4::adjustBoundariesOfMarginalExons(Exon *Lblock) {
  coords *sig;
  char   tmp[50];
  Exon   *newthing;
  Exon   *tmp_block = Lblock->next_exon;

  /* condition for non-signal */
  if (tmp_block && tmp_block->toGEN &&
      (strncmp((char *)(_genSeq+tmp_block->frGEN-3), END_SIG, (size_t)2) || 
       (tmp_block->frEST!=1))) {
    sig = (G_score>=abs(C_score)) ? &last_AG : &last_AC; 
    if (sig->pos1 && (sig->pos2<=20)) {
      /* generated in extend_bw */
      assert(sig->pos2 > 1);
      (void)strcpy((char *)tmp,END_SIG);
      (void)strncpy((char *)(tmp+2),(char *)_estSeq,(size_t)sig->pos2-1);
      (void)strcpy((char *)(tmp+sig->pos2+1), START_SIG);
      newthing = bmatch(_genSeq,tmp,tmp_block->frGEN-3,sig->pos2+3,1,1); 
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
      (strncmp((char *)(_genSeq+tmp_block->toGEN),START_SIG,(size_t)2) || (tmp_block->toEST!=_estLen))) {   
    sig = (G_score>=abs(C_score)) ? &last_GT : &last_CT;
    if (sig->pos1 && (_estLen-sig->pos2<=20)) {
      assert(_estLen-sig->pos2 >= 0);
      (void)strcpy((char *)tmp,END_SIG);
      (void)strncpy((char *)(tmp+2),(char *)(_estSeq+sig->pos2),
                    (size_t)_estLen-sig->pos2);
      (void)strcpy((char *)(tmp+_estLen-sig->pos2+2),START_SIG);
      newthing = fmatch(_genSeq+sig->pos1-1,tmp,
                        _genLen-sig->pos1+1,_estLen-sig->pos2+4,
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
    if (!strncmp((char *)(_genSeq+v-2),"AG",(size_t)2)) {
      last_AG.pos1 = v+1;
      last_AG.pos2 = tmp_block1->frEST + (v-tmp_block1->frGEN)+1; 
      break;
    } 

  for (v=tmp_block1->frGEN-1; v<=tmp_block1->toGEN-3; v++)
    if (!strncmp((char *)(_genSeq+v-2),"AC",(size_t)2)) {
      last_AC.pos1 = v+1;
      last_AC.pos2 = tmp_block1->frEST + (v-tmp_block1->frGEN)+1;
      break;
    }
}




void
Sim4::findLastGTandCT(Exon *tmp_block) {
  int v;

  for (v=tmp_block->toGEN; v>=tmp_block->frGEN; v--)
    if (!strncmp((char *)(_genSeq+v),"GT",(size_t)2)) {
      last_GT.pos1 = v;
      last_GT.pos2 = tmp_block->toEST-(tmp_block->toGEN-v);
      break;
    }

  for (v=tmp_block->toGEN; v>=tmp_block->frGEN; v--)
    if (!strncmp((char *)(_genSeq+v),"CT",(size_t)2)) {
      last_CT.pos1 = v;
      last_CT.pos2 = tmp_block->toEST-(tmp_block->toGEN-v);
      break;
    }
}

