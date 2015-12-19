
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_CNS/MultiAlignment_CNS.h
 *    src/AS_CNS/MultiAlignment_CNS_private.H
 *    src/AS_CNS/MultiAlignment_CNS_private.h
 *    src/utgcns/libcns/MultiAlignment_CNS_private.H
 *
 *  Modifications by:
 *
 *    Gennady Denisov from 2005-MAY-23 to 2007-OCT-25
 *      are Copyright 2005-2007 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUL-08 to 2013-AUG-01
 *      are Copyright 2005-2009,2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2006-FEB-13 to 2008-FEB-13
 *      are Copyright 2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-JAN-28 to 2009-SEP-25
 *      are Copyright 2008-2009 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2011-SEP-21
 *      are Copyright 2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-17 to 2015-JUL-23
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ABACUS_H
#define ABACUS_H

static const char *rcsid_ABACUS_H = "$Id:  $";

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

//  Probably can't change these

#define CNS_MIN_QV 0
#define CNS_MAX_QV 60

//  Tables and other static data - needs to be global so subclasses can access it

extern bool    DATAINITIALIZED;

extern double  EPROB[CNS_MAX_QV - CNS_MIN_QV + 1];  // prob of error for each quality value
extern double  PROB [CNS_MAX_QV - CNS_MIN_QV + 1];  // prob of correct call for each quality value (1-eprob)

#define CNS_NUM_SYMBOLS  6  //  -ACGTN

extern uint32   baseToIndex[256];
extern char     indexToBase[CNS_NUM_SYMBOLS];


typedef enum {
  abAbacus_Smooth  = 1, // only eliminate pairwise construction artifacts
  abAbacus_Poly_X  = 2, // align poly-X regions
  abAbacus_Indel   = 4  // push apart mushed block indels
} abAbacusRefineLevel;


#include "abBead.H"
#include "abColumn.H"
#include "abSequence.H"


class abAbacus {
public:
  abAbacus() {
    _sequencesLen = 0;
    _sequencesMax = 65536;
    _sequences    = new abSequence * [_sequencesMax];

    _columnsLen   = 0;
    _columnsMax   = 1024 * 1024;
    _columns      = new abColumn * [_columnsMax];

    _cnsBases     = new char  [_columnsMax];
    _cnsQuals     = new uint8 [_columnsMax];

    _firstColumn  = NULL;

    if (DATAINITIALIZED == false)
      initializeGlobals();
  };
  ~abAbacus() {
    delete [] _sequences;
    delete [] _columns;
    delete [] _cnsBases;
    delete [] _cnsQuals;
  };

private:
  void  initializeGlobals(void);

public:

  //  Adds gkpStore read 'readID' to the abacus; former AppendFragToLocalStore
  void            addRead(gkStore *gkpStore,
                          uint32 readID,
                          uint32 askip, uint32 bskip,
                          bool complemented,
                          map<uint32, gkRead *>     *inPackageRead,
                          map<uint32, gkReadData *> *inPackageReadData);

public:
  //abBeadID        unalignBeadFromColumn(abBeadID bid);
  //void            unalignTrailingGapBeads(abBeadID bid);

  //void            checkColumnBaseCounts();

  //void            lateralExchangeBead(abBeadID lid, abBeadID rid);

public:
  void    refreshMultiAlign(bool  recallBase  = false,
                            bool  highQuality = false);

  void    appendBases(uint32  bid,
                      uint32  bgn,
                      uint32  end);

  void    applyAlignment(uint32  bid,
                         int32   ahang,
                         int32   bhang,
                         int32  *trace, uint32 traceLen);


public:
  //  The iterators make calls to abAbacus, and to make them inline, they need to be defined after abAbacus.
  //  That means these constructors cannot be inlined.
  //abSeqBeadIterator  *createSeqBeadIterator(abSeqID sid);
  //abColBeadIterator  *createColBeadIterator(abColID cid);

  //
  //  Objects that we're trying to build multialigns for
  //

public:
  uint32             numberOfSequences(void) {
    return(_sequencesLen);
  };

  abSequence   *getSequence(uint32 sid) {
    assert(sid < _sequencesLen);
    return(_sequences[sid]);
  };
private:
  uint32            _sequencesLen;
  uint32            _sequencesMax;
  abSequence      **_sequences;

  //  The consensus sequence (actually, the multialign structure itself) Bases in columns get
  //  shifted, but columns themselves do not move.  However, we can insert columns.  This array is
  //  an in order list of the columns in the multialignment.

public:
  abColumn          *getFirstColumn(void) { return(_columns[0]);             };
  abColumn          *getLastColumn(void)  { return(_columns[_columnsLen-1]); };

  uint32             numberOfColumns(void) { return(_columnsLen); };

  abColumn          *getColumn(uint32 cid) {
    assert(cid < _columnsLen);
    return(_columns[cid]);
  };
  
  uint32            _columnsLen;
  uint32            _columnsMax;
  abColumn        **_columns;

  char             *bases(void) { return(_cnsBases); };
  uint8            *quals(void) { return(_cnsQuals); };

  char             *_cnsBases;
  uint8            *_cnsQuals;

  abColumn         *_firstColumn;

public:

  struct beadID {
    abColumn  *column;
    uint16     link;
  };

  map<beadID,uint32>  beadToRead;
  map<uint32,beadID>  readToBead;

  //  This is the former abMultiAlign
public:
  void                   mergeRefine(bool highQuality);

  void                   getConsensus(tgTig *tig);
  uint32                 getSequenceDeltas(uint32 sid, int32 *deltas);
  void                   getPositions(tgTig *tig);

  int32                  refineWindow(abColumn *bgnCol_column, abColumn *terCol);
  int32                  refine(abAbacusRefineLevel  level,
                                uint32               from = 0,
                                uint32               to   = UINT32_MAX);

  //  There are two multiAlign displays; this one, and one in tgTig.
  void                   display(FILE *F);
};

#endif
