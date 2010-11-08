#ifndef SIM4_POLISH_H
#define SIM4_POLISH_H

//
//  Datastructures for writing, processing and reading the output of sim4
//

#include <stdio.h>
#include <stdlib.h>

#include "util++.H"
#include "bio++.H"
//#include "bio.h"

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

#define S4P_POLISH_S4DB          100
#define S4P_POLISH_GFF3          101
#define S4P_POLISH_ATAC          102

enum sim4polishStyle {
  sim4polishS4DB,
  sim4polishGFF3,
  sim4polishATAC
};

extern sim4polishStyle sim4polishStyleDefault;

//  sim4polishExon and sim4polish constructors should be private (with whatever builds them a
//  friend) but snapper2 needs to create empty objects and fill them manually.  We could provide a
//  constructor for that, but time is finite.

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
public:
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


  sim4polish() {
    clear();
  };

  friend class sim4polishBuilder;
  friend class sim4polishReader;
  friend class sim4polishWriter;
public:
  //  OBSOLETE
  //sim4polish(FILE *F) {
  //  fprintf(stderr, "OBSOLETE.\n");
  //  exit(1);
  //};

  sim4polish(readBuffer *rb, sim4polishStyle style) {
    clear();

    switch (style) {
      case sim4polishS4DB:
        s4p_readPolishS4DB(rb);
        break;
      case sim4polishGFF3:
        s4p_readPolishGFF3(rb);
        break;
      case sim4polishATAC:
        s4p_readPolishATAC(rb);
        break;
      default:
        fprintf(stderr, "sim4polish()-- ERROR: unknown style '%d'\n", style);
        exit(1);
    }
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
  void   s4p_readPolishS4DB(readBuffer *rb);
  void   s4p_readPolishGFF3(readBuffer *rb);
  void   s4p_readPolishATAC(readBuffer *rb);

  void   s4p_linesToPolishS4DB(u32bit lineNumber, u32bit maxLines, char **lines, u32bit *lengths);
  void   s4p_linesToPolishGFF3(u32bit lineNumber, u32bit maxLines, char **lines, u32bit *lengths);
  void   s4p_linesToPolishATAC(u32bit lineNumber, u32bit maxLines, char **lines, u32bit *lengths);

  void   s4p_copyPolish(sim4polish *orig, u32bit exonNum=2147483648);

public:

  //  Note that there is no (public) mechanism to convert these strings back to a sim4polish.  The
  //  only mechanism is through a readBuffer (aka, a file).
  //
  char  *s4p_polishToString(sim4polishStyle style);  // STYLE - add =sim4polishS4DB
private:
  char  *s4p_polishToStringS4DB(void);
  char  *s4p_polishToStringGFF3(void);
  char  *s4p_polishToStringATAC(void);

public:
  void   s4p_removeAlignments(void) {
    for (u32bit i=0; i<_numExons; i++) {
      delete [] _exons[i]._estAlignment;   _exons[i]._estAlignment = 0L;
      delete [] _exons[i]._genAlignment;   _exons[i]._genAlignment = 0L;
    }
  };
  void   s4p_removeDefLines(void) {
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
  u32bit           _querySeqIdentity;    //  numCovered / (estLen - pA -pT)
  u32bit           _matchOrientation;
  u32bit           _strandOrientation;

  char            *_comment;
  char            *_estDefLine;
  char            *_genDefLine;

  u32bit           _numExons;
  sim4polishExon  *_exons;
};

int  s4p_genIDcompare(const void *a, const void *b);
int  s4p_estIDcompare(const void *a, const void *b);

int  s4p_genDEFcompare(const void *a, const void *b);
int  s4p_estDEFcompare(const void *a, const void *b);

bool s4p_compatable(sim4polish *A, sim4polish *B);
bool s4p_IsSameRegion(sim4polish *A, sim4polish *B, int tolerance);
bool s4p_IsRegionOverlap(sim4polish *A, sim4polish *B);
bool s4p_IsSameExonModel(sim4polish *A, sim4polish *B, int tolerance);

void s4p_compareExons_Overlap(sim4polish *A,
                              sim4polish *B,
                              double      overlapThreshold,
                              u32bit     *numSame,
                              u32bit     *numAOnly,
                              u32bit     *numBOnly);

void s4p_compareExons_Ends(sim4polish *A,
                           sim4polish *B,
                           s32bit      tolerance,
                           u32bit     *numSame,
                           u32bit     *numAOnly,
                           u32bit     *numBOnly);

#endif  //  SIM4_POLISH_H
