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

//  We allow the polish to be printed in different ways:
//
//  WARNING: A flags value of zero MUST ALWAYS mean print the whole polish.
//
#define S4P_PRINTPOLISH_FULL         0x0000
#define S4P_PRINTPOLISH_NOALIGNS     0x0001
#define S4P_PRINTPOLISH_NODEFS       0x0002
#define S4P_PRINTPOLISH_MINIMAL     (S4P_PRINTPOLISH_NOALIGNS | S4P_PRINTPOLISH_NODEFS)

#define SIM4_FORMAT_UNKNOWN       0  //  No format specified, usually an error
#define SIM4_FORMAT_SIM4DB        1  //  Original sim4db, multilines, with alignments
#define SIM4_FORMAT_ATAC          2  //  Reduced sim4db, one record per line
#define SIM4_FORMAT_GFF3          3  //  GFF3, duh.


class sim4polishExon {
public:
  sim4polishExon() {
    _estFrom            = 0;
    _estTo              = 0;
    _genFrom            = 0;
    _genTo              = 0;
    _numMatches         = 0;
    _numMatchesN        = 0;
    _percentIdentity    = 0;
    _intronOrientation  = 0;
    _estAlignment       = NULL;
    _genAlignment       = NULL;
  };

  ~sim4polishExon() {
    delete [] _estAlignment;
    delete [] _genAlignment;
  };

  void   s4p_clearExon(void) {
    _estFrom            = 0;
    _estTo              = 0;
    _genFrom            = 0;
    _genTo              = 0;
    _numMatches         = 0;
    _numMatchesN        = 0;
    _percentIdentity    = 0;
    _intronOrientation  = 0;
    _estAlignment       = 0L;
    _genAlignment       = 0L;
  };

  void   s4p_copyExon(sim4polishExon *orig);

public:
  u32bit           _estFrom;
  u32bit           _estTo;
  u32bit           _genFrom;
  u32bit           _genTo;
  u32bit           _numMatches;
  u32bit           _numMatchesN;
  u32bit           _percentIdentity;
  u32bit           _intronOrientation;
  char            *_estAlignment;
  char            *_genAlignment;
};



class sim4polish {
private:
  void clear(void) {
    _estID             = 0;
    _estLen            = 0;
    _estPolyA          = 0;
    _estPolyT          = 0;

    _genID             = 0;
    _genRegionOffset   = 0;
    _genRegionLength   = 0;

    _numMatches        = 0;
    _numMatchesN       = 0;
    _numCovered        = 0;
    _percentIdentity   = 0;
    _querySeqIdentity  = 0;
    _matchOrientation  = 0;
    _strandOrientation = 0;

    _comment           = NULL;
    _estDefLine        = NULL;
    _genDefLine        = NULL;
    
    _numExons          = 0;
    _exons             = NULL;
  };

public:
  sim4polish() {
    clear();
  };

  sim4polish(FILE *F) {
    clear();
    s4p_readPolish(F);
  };

  sim4polish(sim4polish *orig) {
    clear();
    s4p_copyPolish(orig);
  };

  sim4polish(sim4polish *orig, u32bit exon) {
    clear();
    s4p_copyPolish(orig, exon);
  };

  ~sim4polish() {
    delete [] _comment;
    delete [] _estDefLine;
    delete [] _genDefLine;
    delete [] _exons;
  };

private:
  void   s4p_linesToPolish(u32bit lineNumber, u32bit maxLines, char **lines, u32bit *lengths);

public:
  void   s4p_readPolish(FILE *F);
  void   s4p_copyPolish(sim4polish *orig, u32bit exonNum=2147483648);

  void   s4p_stringToPolish(char *s);
  char  *s4p_polishToString(void);

  void   s4p_printPolish(FILE *O, u32bit flags);

  void   s4p_removeAlignments(void) {
    for (u32bit i=0; i<_numExons; i++) {
      delete [] _exons[i]._estAlignment;   _exons[i]._estAlignment = 0L;
      delete [] _exons[i]._genAlignment;   _exons[i]._genAlignment = 0L;
    }
  };
  void   s4p_removeDeflines(void) {
    delete [] _estDefLine;   _estDefLine = 0L;
    delete [] _genDefLine;   _genDefLine = 0L;
  };

  //  Reverse complement an input polish, returns true of it was reversed.
  //
  bool   s4p_makeForward(void);
  bool   s4p_makeReverse(void);

  //  Update the alignment scores based on the alignments that are present.
  //
  void   s4p_updateAlignmentScores(void);

  //  Approximate (integer) percent identity and coverage.
  //
  int    s4p_percentIdentityApprox(int numEdits, int alignmentLength);
  int    s4p_percentCoverageApprox(void);
  
  //  A very expensive and accurate calculation of the percent identity.
  //
  double s4p_percentIdentityExact(void);
  double s4p_percentCoverageExact(void);

  void   s4p_swapExons(u32bit a, u32bit b);
  void   s4p_deleteExon(u32bit a);
  void   s4p_insertExon(u32bit a, u32bit intronori, sim4polishExon  *e);
  void   s4p_insertExons(u32bit a, u32bit intronori, sim4polish     *e);

public:
  u32bit           _estID;
  u32bit           _estLen;
  u32bit           _estPolyA;
  u32bit           _estPolyT;

  u32bit           _genID;
  u32bit           _genRegionOffset;
  u32bit           _genRegionLength;

  u32bit           _numMatches;
  u32bit           _numMatchesN;
  u32bit           _numCovered;          //  Number of bp covered in alignments
  u32bit           _percentIdentity;
  u32bit           _querySeqIdentity;    //  numMatches / (estLen - pA -pT)
  u32bit           _matchOrientation;
  u32bit           _strandOrientation;

  char            *_comment;
  char            *_estDefLine;
  char            *_genDefLine;

  u32bit           _numExons;
  sim4polishExon  *_exons;
};

int            s4p_genIDcompare(const void *a, const void *b);
int            s4p_estIDcompare(const void *a, const void *b);

int            s4p_genDEFcompare(const void *a, const void *b);
int            s4p_estDEFcompare(const void *a, const void *b);

bool           s4p_compatable(sim4polish *A, sim4polish *B);
bool           s4p_IsSameRegion(sim4polish *A, sim4polish *B, int tolerance);
bool           s4p_IsRegionOverlap(sim4polish *A, sim4polish *B);
bool           s4p_IsSameExonModel(sim4polish *A, sim4polish *B, int tolerance);

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

#endif  //  SIM4_POLISH_H
