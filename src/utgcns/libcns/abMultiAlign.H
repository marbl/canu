
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
 *    Brian P. Walenz beginning on 2014-NOV-17
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

typedef enum {
  abAbacus_Smooth  = 1, // only eliminate pairwise construction artifacts
  abAbacus_Poly_X  = 2, // align poly-X regions
  abAbacus_Indel   = 4  // push apart mushed block indels
} abAbacusRefineLevel;



//  This is the basic multialignment atom: A collection (possibly empty) of columns given by their
//  offsets in the global columnStore.
//
//  Former MANode
//
class abMultiAlign {
public:
  abMultiAlign() {
    lid   = abMultiAlignID();
    iid   = 0;
    first = abColID();
    last  = abColID();

    columnList.clear();
  };
  ~abMultiAlign() {
  };

  abMultiAlignID    const ident(void)          { return(lid); };
  void                    identSet(uint32 idx) { lid.set(idx); };  //  Use ONLY by abAbacus::addMultiAlign()

  abColID           const firstColumn(void) { return(first); };
  abColID           const lastColumn(void)  { return(last); };

  //  GetMANodeLength - also guarded against null nodes and returned -1
  uint32                  length(void) {
    return(columnList.size());
  };

  vector<abColID>        &columns(void) {  //  For unitigConsensus, to rebuild
    return(columnList);
  };

  void                   addColumnToMultiAlign(abColumn *column) {
    columnList.push_back(column->ident());

    if (column->prevID().isValid() == false)
      //  No previous, so this becomes the first column
      first = column->ident();

    if (column->nextID().isValid() == false)
      //  No next, so this becomes the last
      last = column->ident();
  };

  void                   getConsensus(abAbacus *abacus, char *bases, char *quals, uint32 &len, uint32 max);
  void                   getConsensus(abAbacus *abacus, char *bases, char *quals);
  void                   getConsensus(abAbacus *abacus, tgTig *tig);
  void                   getPositions(abAbacus *abacus, tgTig *tig);

  void                   mergeRefine(abAbacus *abacus, bool highQuality);
  int32                  refineWindow(abAbacus *abacus, abColumn *start_column, abColID stab_bgn);
  int32                  refine(abAbacus            *abacus,
                                abAbacusRefineLevel  level,
                                uint32               from = 0,
                                uint32               to   = UINT32_MAX);

  //  There are two multiAlign displays; this one, and one in tgTig.
  void                   display(abAbacus  *abacus,
                                 FILE      *F    = stdout,
                                 uint32     from = 0,
                                 uint32     to   = UINT32_MAX);

private:
  abMultiAlignID         lid;         //  MANode id in the manodeStore
  uint32                 iid;         //  tigStore tigID
  abColID                first;       //  first column in the multialignment
  abColID                last;
  vector<abColID>        columnList;  //  Used in AbacusRefine to get random access to specific column

  //friend class abAbacus;
};
