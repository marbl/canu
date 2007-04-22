#ifndef SIM4_POLISH_H
#define SIM4_POLISH_H

//
//  Datastructures for writing, processing and reading the output of sim4
//

#include <stdio.h>
#include <stdlib.h>

#include "bio.h"


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
  u32bit           genRegionOffset;
  u32bit           genRegionLength;

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

int            s4p_genDEFcompare(const void *a, const void *b);
int            s4p_estDEFcompare(const void *a, const void *b);

sim4polish    *s4p_readPolish(FILE *F);
sim4polish    *s4p_copyPolish(sim4polish *orig);
sim4polish    *s4p_copyPolish_OneExon(sim4polish *orig, int exon);

void           s4p_removeAlignments(sim4polish *p);
void           s4p_removeDeflines(sim4polish *p);
void           s4p_destroyPolish(sim4polish *p);


//  Reverse complement an input polish
//
int            s4p_makeForward(sim4polish *p);
int            s4p_makeReverse(sim4polish *p);

//  Update the alignment scores based on the alignments that are
//  present.
//
void           s4p_updateAlignmentScores(sim4polish *p);

//  Approximate (integer) percent identity and coverage.
//
int            s4p_percentIdentityApprox(int numEdits, int alignmentLength);
int            s4p_percentCoverageApprox(sim4polish *p);

//  A very expensive and accurate calculation of the percent identity.
//
double         s4p_percentIdentityExact(sim4polish *p);
double         s4p_percentCoverageExact(sim4polish *p);


//  We allow the polish to be printed in different ways:
//
//  WARNING:  A flags value of zero MUST ALWAYS mean print the whole
//  polish.
//
#define S4P_PRINTPOLISH_FULL         0x0000
#define S4P_PRINTPOLISH_NOALIGNS     0x0001
#define S4P_PRINTPOLISH_NODEFS       0x0002

#define S4P_PRINTPOLISH_MINIMAL     (S4P_PRINTPOLISH_NOALIGNS | S4P_PRINTPOLISH_NODEFS)

//  The polish is not valuable, and the print routine is allowed to
//  modify it.  If this is not set, and you use NOALIGNS or NODEFS, a
//  copy of the polish is made at some non-zero expense.
//
#define S4P_PRINTPOLISH_NOTVALUABLE  0x8000

void           s4p_printPolish(FILE *O, sim4polish *p, u32bit flags);

void           s4p_swapExons(sim4polish *p, int a, int b);
void           s4p_deleteExon(sim4polish *p, int a);
void           s4p_copyExon(sim4polishExon *copy, sim4polishExon *orig);
void           s4p_insertExon(sim4polish      *p,
                              int              a,
                              u32bit           intronori,
                              sim4polishExon  *e);
void           s4p_insertExons(sim4polish     *p,
                               int             a, 
                               u32bit          intronori,
                               sim4polish     *e);

sim4polish    *s4p_stringToPolish(char *s);
char          *s4p_polishToString(sim4polish *p);

int            s4p_compatable(sim4polish *A, sim4polish *B);
int            s4p_IsSameRegion(sim4polish *A, sim4polish *B, int tolerance);
int            s4p_IsRegionOverlap(sim4polish *A, sim4polish *B);
int            s4p_IsSameExonModel(sim4polish *A, sim4polish *B, int tolerance);
void           s4p_compareExons_Overlap(sim4polish *A,
                                        sim4polish *B,
                                        double      overlapThreshold,
                                        u32bit     *numSame,
                                        u32bit     *numAOnly,
                                        u32bit     *numBOnly);
void           s4p_compareExons_Ends(sim4polish *A,
                                     sim4polish *B,
                                     s32bit      tolerance,
                                     u32bit     *numSame,
                                     u32bit     *numAOnly,
                                     u32bit     *numBOnly);


#ifdef __cplusplus
}
#endif


#endif  //  SIM4_POLISH_H
