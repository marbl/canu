#ifndef SIM4_POLISH_H
#define SIM4_POLISH_H

//
//  Datastructures for writing, processing and reading the output of sim4
//

#include <stdio.h>
#include <stdlib.h>

#include "libbritypes.h"


#define SIM4_INTRON_ERROR        '?'  //  '??'
#define SIM4_INTRON_POSITIVE     '>'  //  '->'
#define SIM4_INTRON_NEGATIVE     '<'  //  '<-'
#define SIM4_INTRON_AMBIGUOUS    '-'  //  '--'
#define SIM4_INTRON_GAP          '='  //  '=='
#define SIM4_INTRON_NONE         '.'  //  '  '

#define SIM4_MATCH_ERROR         '?'
#define SIM4_MATCH_FORWARD       'f'
#define SIM4_MATCH_COMPLEMENT    'c'

#define SIM4_STRAND_ERROR        '?'
#define SIM4_STRAND_POSITIVE     'p'
#define SIM4_STRAND_NEGATIVE     'n'
#define SIM4_STRAND_UNKNOWN      'u'
#define SIM4_STRAND_INTRACTABLE  'I'
#define SIM4_STRAND_FAILED       'F'


typedef struct {
  u32bit           estFrom;
  u32bit           estTo;
  u32bit           genFrom;
  u32bit           genTo;
  u32bit           numMatches;
  u32bit           numMatchesN;
  u32bit           percentIdentity;
  u32bit           intronOrientation;
  char            *estAlignment;
  char            *genAlignment;
} sim4polishExon;


typedef struct {
  u32bit           estID;
  u32bit           estLen;
  u32bit           estPolyA;
  u32bit           estPolyT;

  u32bit           genID;
  u32bit           genLo;
  u32bit           genHi;

  u32bit           numMatches;
  u32bit           numMatchesN;
  u32bit           numCovered;          //  Number of bp covered in alignments
  u32bit           percentIdentity;
  u32bit           querySeqIdentity;    //  numMatches / (estLen - pA -pT)
  u32bit           matchOrientation;
  u32bit           strandOrientation;

  char            *comment;
  char            *estDefLine;
  char            *genDefLine;

  u32bit           numExons;
  sim4polishExon  *exons;
} sim4polish;



#ifdef __cplusplus
extern "C" {
#endif

int            s4p_genIDcompare(const void *a, const void *b);
int            s4p_estIDcompare(const void *a, const void *b);

sim4polish    *s4p_readPolish(FILE *F);
sim4polish    *s4p_copyPolish(sim4polish *orig);
void           s4p_destroyPolish(sim4polish *p);
void           s4p_printPolish(FILE *O, sim4polish *p);
void           s4p_printPolishColumn(FILE *O, sim4polish *p);
void           s4p_printPolishNormalized(FILE *O, sim4polish *p);

void           s4p_swapExons(sim4polish *p, int a, int b);
void           s4p_deleteExon(sim4polish *p, int a);
void           s4p_copyExon(sim4polish *p, int a, sim4polishExon *e);
void           s4p_overwriteExon(sim4polish *p, sim4polishExon *e, int a);
void           s4p_insertExon(sim4polish *p, int a, sim4polishExon *e);

sim4polish    *s4p_stringToPolish(char *s);
char          *s4p_polishToString(sim4polish *p);

#ifdef __cplusplus
}
#endif


#endif  //  SIM4_POLISH_H
