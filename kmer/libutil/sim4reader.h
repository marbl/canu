#ifndef SIM4READER_H
#define SIM4READER_H

#include <stdio.h>
#include <stdlib.h>

#define INTRON_POSITIVE    '>'  //  '->'
#define INTRON_NEGATIVE    '<'  //  '<-'
#define INTRON_AMBIGUOUS   '-'  //  '--'
#define INTRON_GAP         '='  //  '=='
#define INTRON_NONE        '.'  //  '  '

#define MATCH_FORWARD      'f'
#define MATCH_COMPLEMENT   'c'

#define STRAND_POSITIVE    'p'
#define STRAND_NEGATIVE    'n'
#define STRAND_UNKNOWN     'u'
#define STRAND_INTRACTABLE 'I'
#define STRAND_FAILED      'F'

#ifdef __cplusplus
extern "C" {
#endif

//  Utility for reading a whole line, safely, from a file.
//
typedef struct {
  unsigned int   l;
  unsigned int   a;
  char          *s;
  unsigned int   lineNumber;
} _line;

_line  *newLine(void);
void    deleteLine(_line *L);
void    readLine(FILE *F, _line *L);


typedef struct {
  int   estFrom;
  int   estTo;
  int   genFrom;
  int   genTo;

  int   numMatches;
  int   numMatchesN;
  int   percentIdentity;

  int   intronOrientation;

  char *estAlignment;
  char *genAlignment;
} sim4polishExon;


typedef struct {
  int              estID;
  int              estLen;
  int              estPolyA;
  int              estPolyT;

  int              genID;
  int              genLo;
  int              genHi;

  int              numMatches;
  int              numMatchesN;
  int              percentIdentity;
  int              querySeqIdentity;    //  numMatches / (estLen - pA -pT)
  int              matchOrientation;
  int              strandOrientation;

  char            *comment;
  char            *estDefLine;
  char            *genDefLine;

  int              numExons;
  sim4polishExon  *exons;
} sim4polish;


int genIDcompare(const void *a, const void *b);
int estIDcompare(const void *a, const void *b);


sim4polish *readPolish(FILE *F);
sim4polish *copyPolish(sim4polish *orig);
void        destroyPolish(sim4polish *p);
void        printPolish(FILE *O, sim4polish *p);
void        printPolishColumn(FILE *O, sim4polish *p);
void        printPolishNormalized(FILE *O, sim4polish *p);

//  Avoid a little namespace pollution
//
void        s4p_swapExons(sim4polish *p, int a, int b);
void        s4p_deleteExon(sim4polish *p, int a);
void        s4p_copyExon(sim4polish *p, int a, sim4polishExon *e);
void        s4p_overwriteExon(sim4polish *p, sim4polishExon *e, int a);
void        s4p_insertExon(sim4polish *p, int a, sim4polishExon *e);

#if 0
char       *s4p_polishToString(sim4polish *p);
#endif

sim4polish *s4p_stringToPolish(char *s);

#ifdef __cplusplus
}
#endif

#endif  //  SIM4READER_H

