#include "sim4polish.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "memory.h"

void
s4p_swapExons(sim4polish *p, int a, int b) {
  sim4polishExon  copyofa;

  memcpy(&copyofa,   p->exons+a, sizeof(sim4polishExon));
  memcpy(p->exons+a, p->exons+b, sizeof(sim4polishExon));
  memcpy(p->exons+b, &copyofa,   sizeof(sim4polishExon));
}

void
s4p_copyExon(sim4polishExon *copy, sim4polishExon *orig) {
  memcpy(copy, orig, sizeof(sim4polishExon));

  if (copy->estAlignment)
    free(copy->estAlignment);
  copy->estAlignment = 0L;
  if (orig->estAlignment)
    copy->estAlignment = (char *)memdup(orig->estAlignment, sizeof(char) * (strlen(orig->estAlignment) + 1));

  if (copy->genAlignment)
    free(copy->genAlignment); 
  copy->genAlignment = 0L;
 if (orig->genAlignment)
    copy->genAlignment = (char *)memdup(orig->genAlignment, sizeof(char) * (strlen(orig->genAlignment) + 1));
}


void
s4p_insertExon(sim4polish     *p,
               int             a, 
               u32bit          intronori,
               sim4polishExon *e) {
  int i;

  errno = 0;
  p->exons = (sim4polishExon *)realloc(p->exons, sizeof(sim4polishExon) * (p->numExons + 1));
  if (errno) {
    fprintf(stderr, "s4p_insertExon()-- Can't increase the number of exons!\n%s\n", strerror(errno));
    exit(1);
  }

  //  Move exons a and higher up one position
  for (i=p->numExons; i>a; i--)
    memcpy(p->exons+i, p->exons+i-1, sizeof(sim4polishExon));

  //  We trust that the user has set the intron orientation in the new
  //  exon, and that 'intronori' is the correct orientation for the
  //  previous intron.
  //
  if (a > 0)
    p->exons[a-1].intronOrientation = intronori;

  s4p_copyExon(p->exons+a, e);

  p->numExons++;
}


//  Inserts all the exons in e into the polish p at position a.
void
s4p_insertExons(sim4polish     *p,
                int             a, 
                u32bit          intronori,
                sim4polish     *e) {
  int i, j;

  errno = 0;
  p->exons = (sim4polishExon *)realloc(p->exons, sizeof(sim4polishExon) * (p->numExons + e->numExons));
  if (errno) {
    fprintf(stderr, "s4p_insertExon()-- Can't increase the number of exons!\n%s\n", strerror(errno));
    exit(1);
  }

  for (i=p->numExons-1, j=p->numExons + e->numExons - 1; i>=a; i--, j--)
    memcpy(p->exons+j, p->exons+i, sizeof(sim4polishExon));

  //  We trust that the user has set the intron orientation in the new
  //  exon, and that 'intronori' is the correct orientation for the
  //  previous intron.
  //
  if (a > 0)
    p->exons[a-1].intronOrientation = intronori;

  //  Insert the exons
  //
  for (i=0; i<e->numExons; i++, a++)
    s4p_copyExon(p->exons+a, e->exons+i);

  p->numExons += e->numExons;
}
