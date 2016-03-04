
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
 *    src/AS_CNS/AbacusRefine.C
 *    src/AS_CNS/AbacusRefine.c
 *    src/AS_CNS/MultiAlignment_CNS.c
 *    src/utgcns/libcns/AbacusRefine.C
 *
 *  Modifications by:
 *
 *    Michael Schatz on 2004-SEP-23
 *      are Copyright 2004 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2005-MAR-22
 *      are Copyright 2005 The Institute for Genomics Research, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter from 2005-MAR-30 to 2008-FEB-13
 *      are Copyright 2005-2006,2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Gennady Denisov from 2005-MAY-09 to 2008-JUN-06
 *      are Copyright 2005-2008 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-JUN-16 to 2013-AUG-01
 *      are Copyright 2005-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Aaron Halpern from 2005-SEP-29 to 2006-OCT-03
 *      are Copyright 2005-2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2008-FEB-27 to 2010-JUN-09
 *      are Copyright 2008-2010 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Jason Miller on 2011-SEP-21
 *      are Copyright 2011 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-30 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-10
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "abAbacus.H"


#define MAX_WINDOW_FOR_ABACUS_REFINE      100


enum ShiftStatus {
  LEFT_SHIFT  = (int) 'L', // Left Shifted (76)
  RIGHT_SHIFT = (int) 'R', // Right Shifted (82)
  UNSHIFTED   = (int) 'U', // Unshifted (85)
};


class abAbacusWork {
public:
  abAbacusWork() {
    start_column = NULL;
    end_column   = NULL;
    rows         = 0;
    columns      = 0;
    window_width = 0;
    shift        = UNSHIFTED;
    beads        = NULL;
    calls        = NULL;
  };

  abAbacusWork(abAbacus  *abacus,
               abColumn  *from,
               abColumn  *end);

  ~abAbacusWork() {
    delete [] beads;
    delete [] calls;
  };


  void   reset(void) {
    for (uint32 ii=0; ii<columns; ii++)
      calls[ii] = 'n';
  }


  void  setBase(int32 i, int32 j, char c) {
    if ((i < 0) || (rows    <= i) ||
        (j < 0) || (columns <= j))
      fprintf(stderr, "abAbacusWork::set()--  i=%d j=%d out of range of rows=%d columns=%d\n", i, j, rows, columns);

    assert( 0 <= i);
    assert( 0 <= j);

    assert(i < rows);
    assert(j < columns);

    beads[i * (columns + 2) + j + 1] = c;
  };

  char *getPtr(int32 i, int32 j) {
    if ((i < 0) || (rows    <= i) ||
        (j < 0) || (columns <= j))
      fprintf(stderr, "abAbacusWork::getPtr()-- i=%d j=%d out of range of rows=%d columns=%d\n", i, j, rows, columns);

    assert( 0 <= i);
    assert( 0 <= j);

    assert(i < rows);
    assert(j < columns);

    return(beads + i * (columns + 2) + j + 1);
  };

  char  getBase(int32 i, int32 j) {
    if ((i <  0) || (rows      <= i) ||
        (j < -1) || (columns+1 <= j))
      fprintf(stderr, "abAbacusWork::getCharr()-- i=%d j=%d out of range of rows=%d columns=%d\n", i, j, rows, columns);

    assert( 0 <= i);
    assert(-1 <= j);  //  Asking for j=-1 will return the 'n' on the border.

    assert(i <  rows);
    assert(j <= columns);  //  Asking for j=columns will return the 'n' on the border.

    return(beads[i * (columns + 2) + j + 1]);
  };


  abAbacusWork *clone(void) {
    abAbacusWork *clone = new abAbacusWork;

    clone->start_column   = start_column;
    clone->end_column     = end_column;
    clone->window_width   = window_width;
    clone->rows           = rows;
    clone->columns        = columns;
    clone->shift          = shift;
    clone->beads          = new char [rows * (columns+2)];
    clone->calls          = new char [columns];

    memcpy(clone->beads, beads, (rows * (columns+2)) * sizeof(char));
    memcpy(clone->calls, calls, (columns)            * sizeof(char));

    clone->abacus_indices = abacus_indices;

    return(clone);
  };


  void   show(void) {
    char form[32];

    sprintf(form, "'%%%d.%ds'\n", columns, columns);

    fprintf(stderr, "start_column  0x%16p %d\n", start_column, start_column->position());
    fprintf(stderr, "end_column    0x%16p %d\n", end_column,   end_column->position());
    fprintf(stderr, "rows          %d   columns  %d\n", rows, columns);
    fprintf(stderr, "window_width  %d\n", window_width);
    fprintf(stderr, "shift         %c\n", shift);

    //  Shows the 'n' border around each row.
    for (int32 i=0; i<rows; i++) {
      fprintf(stderr, "%03d - '", i);
      for (int32 j=-1; j<=columns; j++) {
        fprintf(stderr, "%c", beads[i * (columns + 2) + j + 1]);
      }
      fprintf(stderr, "'\n");
    }
  };


  uint32   scoreAbacus(int32 &cols);
  uint32   affineScoreAbacus(void);

  int32    merge(int32 merge_dir);

  void     refineOrigAbacus(void);

  uint32   leftShift(int32 &lcols);
  uint32   rightShift(int32 &rcols);

  //  The 'mixedShift', removed in r6675, seemed to shift the minor allele to the left, and the major allele to the right.

  void     applyAbacus(abAbacus *abacus);

public:
  uint32           *abacus_indices;

  abColumn         *start_column;
  abColumn         *end_column;
  int32             rows;
  int32             columns;
  int32             window_width;
  ShiftStatus       shift;

  char             *beads;
  char             *calls;
};







abAbacusWork::abAbacusWork(abAbacus  *abacus,
                           abColumn  *bgn,
                           abColumn  *end) {

  //  See the previous (to mid December 2015) for gigantic amounts of pain related to counting the
  //  number if reads that have any bases between the bgn and end columns.  It was scanning every
  //  64th (or so) column, and checking if the bead with that seqIdx was listed in abacus_indices.
  //  abacus_indices was a vector, initialized with push_back(), once for every read in abAbacus.
  //
  //  Instead, we can just test if each read intersects the columns, and then assign the read ID to
  //  an abAbacusWork-private index.

  assert(0);

  abacus_indices = new uint32 [abacus->numberOfSequences()];

  start_column   = bgn;
  end_column     = end;
  rows           = 0;                                       //  rows is incremented and set below
  window_width   = end->position() - bgn->position() + 1;
  columns        = 3 * window_width;
  shift          = UNSHIFTED;

  for (uint32 ii=0; ii<abacus->numberOfSequences(); ii++) {
    abacus_indices[ii] = UINT32_MAX;

    //  Needs to use abAbacus readTofBead and readTolBead; see unitigConsensus::refreshPositions
#if 0
    if (abacus->getSequence(ii)->lastColumn()->position() < bgn->position())
      continue;
    if (abacus->getSequence(ii)->firstColumn()->position() > end->position())
      continue;
#endif

    abacus_indices[ii] = rows++;
  }

  //  Now we can allocate beads and calls.

  beads          = new char [rows * (columns + 2)];         //  two extra gap columns, plus "null" borders
  calls          = new char [columns];

  //  Clear every base.

  for (uint32 ii=0; ii<rows * (columns+2); ii++)
    beads[ii] = 'n';

  //  Clear the calls[] to all 'n'

  reset();




  //  Fill the center third of abacus with bases.

  for (uint32 ss=0; ss<abacus->numberOfSequences(); ss++) {

    if (abacus_indices[ss] == UINT32_MAX)
      continue;

    //  Read is in the range.  We need to figure out a beadLink for some column, preferably close,
    //  so that we can start iterating over the beads.

    abColumn *bColumn = bgn;
    uint16    bLink   = UINT16_MAX;

    //  Find the earliest column with a link map.  The data for this is NOT IMPLEMENTED.

#if 0
    while (bColumn->beadReadIDsExist() == false)
      bColumn = bColumn->prev();

    //  Find the link, move back to the starting column

    for (uint32 ll=0; ll<bColumn->depth(); ll++)
      if (ss == bColumn->beadReadID(ll)) {
        bLink = ll;
        break;
      }
#endif

    assert(bLink != UINT16_MAX);

    while (bColumn != bgn) {
      bLink   = bColumn->bead(bLink)->nextOffset();
      bColumn = bColumn->next();
    }

    //  Now just add bases into abacus.

    uint32 columnIndex = 0;

    while (1) {
      setBase(abacus_indices[ss], window_width + columnIndex, bColumn->bead(bLink)->base());

      if (bColumn == end)
        break;

      bLink   = bColumn->bead(bLink)->nextOffset();
      bColumn = bColumn->next();
    }
  }

  //  Clear the border columns.  This is done last, so that we can initialize beads with no read
  //  coverage to 'n' (from the first initialization), and then clear the borders to '-' (here).

  for (uint32 i=0; i<rows; i++) {
    for (uint32 j=0; j<window_width; j++)
      setBase(i, j, '-');

    for (uint32 j=2 * window_width; j<columns; j++)
      setBase(i, j, '-');
  }
}






// cols is the number of "good" (non-null) columns found
// GD: This function counts the total number of bases which
//   - are different from column's "consensus" call and
//   - are not 'n'
//
uint32
abAbacusWork::scoreAbacus(int32 &cols) {
  uint32 score = 0;

  cols = 0;

  for (int32 cc=0; cc<columns; cc++) {
    uint32  counts[256] = { 0 };

    //  Sum the bases in this column.

    for (int32 rr=0; rr<rows; rr++) {
      char b = getBase(rr, cc);

      if ((b == '-' ) && (cc > 0) && (cc < columns - 1) &&
          ((getBase(rr, cc-1) == 'n')  ||
           (getBase(rr, cc+1) == 'n')))
        b = 'n';

      counts[b]++;
    }

    //  Pick the majority.

    char  best = 0;

    for (uint32 ii=0; ii<256; ii++)
      if (counts[best] < counts[ii])
        best = ii;

    //  Pick a call, and score the non-gap columns.

    if (counts['-'] + counts['n'] == rows) {
      calls[cc] = 'n';

    } else {
      cols++;

      calls[cc] = best;
      score    += rows - counts[best] - counts['n'];
    }
  }

  return(score);
}


uint32
abAbacusWork::affineScoreAbacus(void) {
  // This simply counts the number of opened gaps, to be used in tie breaker
  //   of edit scores.
  uint32 score=0;
  uint32 start_column;
  uint32 end_column;

  if (shift == LEFT_SHIFT) {
    start_column = 0;
    end_column   = columns/3;

  } else if (shift == RIGHT_SHIFT) {
    start_column = 2*columns/3;
    end_column   =   columns;

  } else {  //  shift == UNSHIFTED
    start_column =   columns/3;
    end_column   = 2*columns/3;
  }

  // Size of a gap does not matter, their number in a row does - GD

  for (uint32 i=0; i<rows; i++) {
    bool  in_gap = false;

    for (int32 j=start_column; j<end_column; j++) {
      if (getBase(i, j) != '-' ) {
        in_gap = false;

      } else if (in_gap == false) {
        in_gap = true;
        score++;
      }
    }
  }

  return(score);
}

int
abAbacusWork::merge(int32 merge_dir) {
  // sweep through abacus from left to right
  // testing for Level 1 (neighbor) merge compatibility of each column
  // with right neighbor and merge if compatible
  //
  //  GD: this code will merge practically any
  int32  mergeok, next_column_good, curr_column_good;
  char   prev, curr, next;
  int32  last_non_null = columns - 1;
  int32 first_non_null = 0;
  int32 columns_merged = 0;

  //fprintf(stderr, "merge()-- dir=%d columns=%d\n", merge_dir, columns);

  // determine the rightmost and leftmost columns
  // not totally composed of gaps
  for (int32 j=columns-1;j>0;j--)
    {
      int32 null_column = 1;
      for (int32 i=0; i<rows; i++) {
        curr = getBase(i,j);
        if (curr != '-') null_column = 0;
      }
      if (!null_column)
        break;
      last_non_null = j;
    }
  for (int32 j=0; j<columns;j++)
    {
      int32 null_column = 1;
      for (int32 i=0; i<rows; i++) {
        curr = getBase(i,j);
        if (curr != '-')
          null_column = 0;
      }
      if (!null_column)
        break;
      first_non_null = j;
    }

  //fprintf(stderr, "columns=%d first_non_null = %d last_non_null= %d\n",
  //        columns, first_non_null, last_non_null);

  if (merge_dir < 0)
    {
      for (int32 j=0;j<last_non_null;j++)
        {
          int32 num_gaps=0, num_ns=0;
          mergeok = 1;
          next_column_good = -1;
          for (int32 i=0;i<rows;i++)
            {
              curr = getBase(i,j);
              next = getBase(i,j+1);
              // at least in one column there should be a gap
              // or, alternatively, both should be 'n'
              if (curr != '-' && next != '-') {
                if (curr != 'n' || next != 'n') {
                  mergeok = 0;
                  break;
                }
                else
                  num_ns++;
              }
              else
                num_gaps++;

              // next column should contain at least one good base - a, c, g or t
              if (next != '-' && next != 'n') {
                next_column_good = i;
              }
            }
          //fprintf(stderr, "column= %d mergeok= %d next_column_good= %d\n", j, mergeok, next_column_good);
          if (mergeok && next_column_good >= 0 && num_gaps > num_ns)
            {
              columns_merged++;
              for (int32 i=0;i<rows;i++) {
                curr = getBase(i,j  );
                next = getBase(i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (next != '-' && next != 'n' )
                  {
                    setBase(i, j  , next);
                    setBase(i, j+1, curr);
                  }
              }
              // The entire j+1-th column now contains only gaps or n's
              // Remove it by shifting all the subsequent columns
              // one position to the left
              for (int32 i=0;i<rows;i++)
                {
                  curr = getBase(i,j  );
                  next = getBase(i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j+1; k<last_non_null; k++)
                    {
                      next= getBase(i,k+1);
                      setBase(i, k, next);
                    }
                  setBase(i, last_non_null, '-');
                }
              // Return to the previous coljumn to see if it can be merged again
              j--;
            }
        }
    }
  else /* merge_dir > 0 */
    {
      for (int32 j=last_non_null-1; j>first_non_null; j--)
        {
          int32 num_gaps=0, num_ns=0;
          mergeok = 1;
          curr_column_good = -1;
          for (int32 i=0;i<rows;i++)
            {
              curr = getBase(i,j);
              next = getBase(i,j+1);
              // in at least one column there should be a gap
              // or, alternatively, both should be 'n'
              if (curr != '-' && next != '-') {
                if (curr != 'n' || next != 'n') {
                  mergeok = 0;
                  break;
                }
                else
                  num_ns++;
              }
              else
                num_gaps++;
              // current column should contain at least one good base - a, c, g or t
              if (curr != '-' && curr != 'n')
                {
                  curr_column_good = i;
                }
            }
          //fprintf(stderr, "column= %d mergeok= %d next_column_good= %d\n", j, mergeok, next_column_good);
          if (mergeok && curr_column_good >= 0 && num_gaps > num_ns)
            {
              columns_merged++;
              for (int32 i=0;i<rows;i++) {
                curr = getBase(i,j  );
                next = getBase(i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (curr != '-' && curr != 'n' ) {
                  setBase(i, j  , next);
                  setBase(i, j+1, curr);
                }
              }
              // The entire j-th column contains gaps
              // Remove it by shifting all the previous columns
              // one position to the right
              for (int32 i=0;i<rows;i++)
                {
                  curr = getBase(i,j  );
                  next = getBase(i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j; k>first_non_null; k--)
                    {
                      prev = getBase(i,k-1);
                      setBase(i, k, prev);
                    }
                  setBase(i, first_non_null, '-');
                }
              // Return to the next column to see if it can be merged again
              j++;
            }
        }
    }

  return(columns_merged);
}




// lcols is the number of non-null columns in result

uint32
abAbacusWork::leftShift(int32 &lcols) {

  //fprintf(stderr, "leftShift()--\n");

  reset();
  //show();

  for (int32 j=window_width; j<2*window_width; j++) {
    for (int32 i=0; i<rows; i++) {

      //for (k=0; k<vreg.na; k++) {
      //  for (l=0; l<vreg.alleles[k].num_reads; l++) {
      //    i = vreg.alleles[k].read_ids[l];

      char  c = getBase(i, j);
      int32 ccol = j;

      if (c != '-' ) {
        //look to the left for a suitable placement
        // will be safe on left since abacus has 'n' border
        while (getBase(i, ccol-1) == '-')
          ccol--;

        // from ccol back up to j, look for column with matching call
        for (int32 pcol = ccol; pcol<j; pcol++) {
          char  call = calls[pcol];

          if (call != 'n' && call != c && c != 'n')
            // GD: consensus in a column == '-' ?
            continue;

          if (call == 'n') {
            // GD:
            // Found the leftmost column with non-gap consensus.
            // Now, reset its consensus "dynamically" to the
            // current base
            // Potential problem: the result will generally
            // depend on the order in which rows
            // are processed
            calls[pcol] = c;
#if 0
            fprintf(stderr, "j= %d i= %d calls[%d]= %c\n", j, i, pcol, c);
#endif
          }

          if (calls[pcol] == c || c == 'n') {
            // swap bases in columns pcol and j of row i
            setBase(i, j, '-');
            setBase(i, pcol, c);
            break;
          }
        }

        if (getBase(i, j) != '-')
          calls[j] = c;
      }
    }
  }

  //fprintf(stderr, "leftShift()-- PREMERGE\n");
  //show();

  merge(-1);

  //fprintf(stderr, "leftShift()-- POSTMERGE\n");
  //show();

  shift = LEFT_SHIFT;

  return(scoreAbacus(lcols));
}


uint32
abAbacusWork::rightShift(int32 &rcols) {
  // rcols is the number of non-null columns in result
  //int32 i, j, k, l, ccol, pcol;
  //char c, call;

  //fprintf(stderr, "rightShift()--\n");

  reset();
  //show();

  for (int32 j=2*window_width-1;j>window_width-1;j--) {
    for (int32 i=0; i<rows; i++) {
      //for (k=0; k<vreg.na; k++) {
      //for (l=0; l<vreg.alleles[k].num_reads; l++) {
      //i = vreg.alleles[k].read_ids[l];
      char  c = getBase(i,j);
      int32 ccol = j;

      if (c != '-' ) {
        //look to the right for a suitable placement
        // will be safe on right since abacus has 'n' border
        while (getBase(i,ccol+1) == '-' )
          ccol++;
        // now, from ccol back down to j, look for column with matching call
        for (int32 pcol = ccol;pcol>j;pcol--) {
          char call = calls[pcol];
          if (call != 'n' && call != c && c != 'n' )
            continue;
          if (call == 'n')
            calls[pcol] = c;
          if (calls[pcol] == c || c == 'n' ) {
            setBase(i,j,'-');
            setBase(i,pcol,c);
            break;
          }
        }
        if (getBase(i,j) != '-' )
          calls[j] = c;
      }
    }
  }

  //fprintf(stderr, "rightShift()-- PREMERGE\n");
  //show();

  merge(1);

  //fprintf(stderr, "rightShift()-- POSTMERGE\n");
  //show();

  shift = RIGHT_SHIFT;

  return(scoreAbacus(rcols));
}



#if 0
//  Relationship must be one of:
//
//  a) end gap moving left:
//
//     X > A > B > C > ... > -   becomes  X - A B C ...
//         ^________________/
//
//  b) non-gap moving left across only gap characters
//    (more efficient special case, since first gap and last
//     character can just be exchanged)
//
//     X > - > - > - > ... > A   becomes  X A - - - ...
//         ^________________/
static
void
leftEndShiftBead(abAbacus *abacus, abBeadID bid, abBeadID eid) {
  abBead    *shift = abacus->getBead(eid);
  abBeadID   aid   = abacus->getBead(bid)->prevID();

  assert(shift != NULL);

  if (abacus->getBase(shift->baseIdx()) != '-' ) {
    // assume first and internal characters are gaps
    abacus->lateralExchangeBead(bid, eid);
  }

  else {
    while (shift->prevID() != aid ) {
      abacus->lateralExchangeBead(shift->prevID(), shift->ident());
    }
  }
}


//  Relationship must be one of:
//
//  a) end gap moving right:
//
//      - > A > B > ... > C > X  becomes  A B ... C - X
//      \_________________^
//
//  b) non-gap moving right across only gap characters
//    (more efficient special case, since first gap and last
//     character can just be exchanged)
//
//      A > - > - > ... > - > X  becomes  - - - ... A X
//       \________________^
static
void
rightEndShiftBead(abAbacus *abacus, abBeadID bid, abBeadID eid) {
  abBead    *shift = abacus->getBead(bid);
  abBeadID   aid   = abacus->getBead(eid)->nextID();

  assert(shift != NULL);

  if (abacus->getBase(shift->baseIdx()) != '-' ) {
    // assume last and internal characters are gaps
    abacus->lateralExchangeBead(bid, eid);
  }

  else {
    while (shift->nextID() != aid ) {
      abacus->lateralExchangeBead(shift->ident(), shift->nextID());
    }
  }
}



void
abAbacusWork::applyAbacus(abAbacus *abacus) {

  beadID bid;  //  ALWAYS the id of bead
  beadID eid;  //  ALWAYS the id of exch

  //fprintf(stderr, "applyAbacus()--  shift=%c start=%d width=%d\n", shift, start_column.get(), window_width);
  //show();

  if (shift == LEFT_SHIFT) {
    abColumn *column = start_column;

    for (uint32 columnCount=0; columnCount < window_width; ) {
      bid.column = column;
      bid.link   = 0;

      //fprintf(stderr, "0; bid=%d eid=%d\n", bid.get(), eid.get());
      // Update all beads in a given column

      while (bid.link != UINT16_MAXisValid()) {
        bead = abacus->getBead(bid);
        char a_entry = getBase(abacus_indices[bead->seqIdx().get()] - 1, columnCount);

        //fprintf(stderr, "a_entry=%c bead=%c\n", a_entry, abacus->getBase(bead->baseIdx()));

        if (a_entry == 'n') {
          eid  = bead->upID();
          exch = abacus->getBead(eid);
          //fprintf(stderr, "1; bid=%d eid=%d\n", bid.get(), eid.get());
          abacus->unalignTrailingGapBeads(bid);

        } else if (a_entry != abacus->getBase(bead->baseIdx())) {
          //  Look for matching bead in frag and exchange
          eid  = bead->ident();
          exch = abacus->getBead(eid);

          //fprintf(stderr, "2; bid=%d eid=%d\n", bid.get(), eid.get());

          if (NULL == exch) {
            eid  = abacus->appendGapBead(bead->ident());
            bead = abacus->getBead(bid);
            abacus->alignBeadToColumn(abacus->getColumn(bead->colIdx())->nextID(),eid, "ApplyAbacus(1)");
            exch = abacus->getBead(eid);
          }

          //fprintf(stderr, "3; bid=%d eid=%d\n", bid.get(), eid.get());

          while (a_entry != abacus->getBase(exch->baseIdx())) {
            abBeadID eidp = exch->nextID();

            if (exch->nextID().isInvalid()) {
              eidp = abacus->appendGapBead(exch->ident());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);
              abacus->alignBeadToColumn(abacus->getColumn(exch->colIdx())->nextID(),eidp, "ApplyAbacus(2)");
              //fprintf(stderr, "4; bid=%d eid=%d\n", bid.get(), eid.get());

            } else if (exch->colIdx() == end_column) {
              eidp = abacus->appendGapBead(exch->ident());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);
              //fprintf(stderr, "5; bid=%d eid=%d\n", bid.get(), eid.get());

              abColumn *cid = column;
              abacus->appendColumn(exch->colIdx(),eidp);
              column = abacus->getColumn(cid);
            }

            //fprintf(stderr, "6; bid=%d eid=%d b col/frg=%d/%d e_col/frg=%d/%d\n",
            //        bid.get(), eid.get(),
            //        bead->colIdx().get(), bead->seqIdx().get(),
            //        exch->colIdx().get(), exch->seqIdx().get());

            eid  = eidp;
            exch = abacus->getBead(eid);
          }

          //fprintf(stderr,"LeftShifting bead %d (%c) with bead %d (%c).\n",
          //        bid.get(), abacus->getBase(bead->baseIdx()),
          //        eid.get(), abacus->getBase(exch->baseIdx()));

          leftEndShiftBead(abacus, bid, eid);
        } else {
          // no exchange necessary;
          eid  = bid;
          exch = bead;

          //fprintf(stderr, "7; bid=%d eid=%d\n", bid.get(), eid.get());
        }

        bid  = exch->downID();
        bead = NULL;

        //fprintf(stderr,"New bid is %d (%c), from %d down\n",
        //        bid.get(), (bid.isValid()) ? abacus->getBase(abacus->getBead(bid)->baseIdx()) : 'n',
        //        eid.get());
      }

      //  End of update; call base now.

      base = abacus->baseCall(column->ident(), true);

      column = column->next();
      columnCount++;
    }
  }


  if (shift == RIGHT_SHIFT) {
    abColumn *column = end_column;

    for (uint32 columnCount=0; columnCount < window_width; ) {
      bid = abacus->getBead(column->callID())->downID();

      while (bid.isValid()) {
        bead = abacus->getBead(bid);
        a_entry = getBase(abacus_indices[bead->seqIdx().get()] - 1, columns - columnCount - 1);

        if (a_entry == 'n') {
          eid  = bead->upID();
          exch = abacus->getBead(eid);
          abacus->unalignTrailingGapBeads(bid);
        } else if (a_entry != abacus->getBase(bead->baseIdx())) {
          //  Look for matching bead in frag and exchange
          eid  = bead->ident();
          exch = abacus->getBead(eid);

          if (NULL == exch) {
            eid  = abacus->prependGapBead(bead->ident());
            bead = abacus->getBead(bid);
            exch = abacus->getBead(eid);
            abacus->alignBeadToColumn(abacus->getColumn(bead->colIdx())->prevID(),eid, "ApplyAbacus(3)");
          }

          while (a_entry != abacus->getBase(exch->baseIdx())) {
            abBeadID eidp = exch->prevID();

            if (exch->prevID().isInvalid()) {
              eidp = abacus->prependGapBead(exch->ident());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);
              abacus->alignBeadToColumn(abacus->getColumn(exch->colIdx())->prevID(),eidp, "ApplyAbacus(4)");
            } else if (exch->colIdx() == start_column) {
              eidp = abacus->appendGapBead(exch->prevID());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);

              abColumn *cid = column;
              abacus->appendColumn(abacus->getColumn(exch->colIdx())->prevID(), eidp);
              column = abacus->getColumn(cid);
            }

            eid  = eidp;
            exch = abacus->getBead(eid);
          }

          //fprintf(stderr,"RightShifting bead %d (%c) with bead %d (%c).\n",
          //        eid.get(), abacus->getBase(exch->baseIdx()),
          //        bid.get(), abacus->getBase(bead->baseIdx()));

          rightEndShiftBead(abacus, eid, bid);
        } else {
          eid  = bid;
          exch = bead; // no exchange necessary;
        }

        bid  = exch->downID();
        bead = NULL;

        //fprintf(stderr,"New bid is %d (%c), from %d down\n",
        //        bid.get(), (bid.isValid()) ? abacus->getBase(abacus->getBead(bid)->baseIdx()) : 'n',
        //        eid.get());
      }

      base = abacus->baseCall(column->ident(), true);
      column = abacus->getColumn(column->prevID());
      columnCount++;
    }
  }
}

#endif


void
abAbacusWork::applyAbacus(abAbacus *abacus) {
}


//
//  In this case, we just look for a string of gaps in the consensus sequence
//
static
int32
IdentifyWindow_Smooth(abColumn *&bgnCol,
                      abColumn *&terCol) {

  terCol = bgnCol->next();

  //  Consensus not a gap, nothing to do.
  if (bgnCol->baseCall() != '-')
    return(0);

  int32  winLen = 1;

  //  Consensus is a gap.  Expand it to the maximum gap.

  while ((terCol->baseCall() == '-') &&
         (terCol->next() != NULL)) {
    terCol = terCol->next();
    winLen++;
  }

#ifdef DEBUG_IDENTIFY_WINDOW
  fprintf(stderr, "identifyWindow()-- gap at %d to %d  winLen=%d (return)\n",
          bgnCol->position(), ter->position(), winLen);
#endif

  return(winLen);
}




//
//  Here, we're looking for a string of the same character
//
static
int32
IdentifyWindow_Poly_X(abColumn *&bgnCol,
                      abColumn *&terCol) {

  terCol = bgnCol->next();

  int32 gapCount  = bgnCol->baseCount('-');
  int32 winLen    = 1;

  char  poly      = bgnCol->baseCall();

  if (poly == '-')
    return(0);

  while ((terCol->baseCall() == poly) || (terCol->baseCall() == '-'))  {
    if (terCol->next() == NULL)
      break;

    gapCount  += terCol->baseCount('-');
    winLen    += 1;

    terCol    = terCol->next();
  }

  if (winLen <= 2)
    return(0);

  // capture trailing gap-called columns

  while (terCol->baseCall() == '-' )  {
    if (terCol->majorityBase(true) != poly)
      break;

    if (terCol->next() == NULL)
      break;

    gapCount  += terCol->baseCount('-');
    winLen    += 1;

    terCol     = terCol->next();
  }

  // now that a poly run with trailing gaps is established, look for leading gaps

  abColumn *preCol = bgnCol->prev();

  while (preCol != NULL) {
    if ((preCol->baseCall() != '-') && (preCol->baseCall() != poly))
      break;

    bgnCol     = preCol;

    gapCount  += preCol->baseCount('-');
    winLen    += 1;

    preCol     = preCol->prev();
  }

  //fprintf(stderr,"POLYX candidate (%c) at column %d ter %d , width %d, gapcount %d\n",
  //        poly, bgnCol->position(), terCol.get(), winLen, gapCount);

  if ((bgnCol->prev() != NULL) &&
      (winLen   > 2) &&
      (gapCount > 0))
    return(winLen);

  return(0);
}





//
//  In this case, we look for a string mismatches, indicating a poor alignment region
//  which might benefit from Abacus refinement
//
//  heuristics:
//
//  > terle border on either side of window of width:  STABWIDTH
//  > fewer than STABMISMATCH in stable border
//
//  _              __              ___
//  SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
//  SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
//  SSSSS SSSSS => SSSSS .SSSS+ => SSSSS  .SSSS+
//  SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
//  SSSSS_SSSSS    SSSSS_.SSSS+    SSSSS__.SSSS+
//  |              |               |
//  |\_____________|_______________|____ growing 'gappy' window
//  |
//  bgnCol
//

#define STABWIDTH  6

static
int32
IdentifyWindow_Indel(abColumn *&bgnCol,
                     abColumn *&terCol) {

  terCol = bgnCol->next();

  int32 cum_mm     = bgnCol->mismatchCount();
  int32 stab_mm    = 0;
  int32 stab_gaps  = 0;
  int32 stab_width = 0;
  int32 stab_bases = 0;
  int32 winLen     = 0;

  if ((cum_mm == 0) || (bgnCol->baseCount('-') == 0))
    return(0);

  abColumn  *stableEnd = terCol;

  //  Compute the number of mismatches, gapes and bases in the next STABWIDTH columnd.

  while ((stableEnd->next() != NULL) && (stab_width < STABWIDTH)) {
    stab_mm    += stableEnd->mismatchCount();
    stab_gaps  += stableEnd->baseCount('-');
    stab_bases += stableEnd->depth();

    stableEnd   = stableEnd->next();

    stab_width++;
  }

  //  If no bases (how?) return an empty window.

  if (stab_bases == 0)
    return(0);

  //  While the number of mismatches is high, and the number of gaps is also high, shift the stable
  //  region ahead by one column.  Subtract out the values for the column we shift out (on the left)
  //  and add in the values for the column we shift in (on the right).

  while ((stab_mm   > 0.02 * stab_bases) ||  //  CNS_SEQUENCING_ERROR_EST
         (stab_gaps > 0.25 * stab_bases)) {

    int32 mm  = terCol->mismatchCount();
    int32 gp  = terCol->baseCount('-');
    int32 bps = terCol->depth();

    //  move terCol ahead

    if (stableEnd->next() != NULL) {
      stab_mm    += stableEnd->mismatchCount();
      stab_bases += stableEnd->depth();
      stab_gaps  += stableEnd->baseCount('-');

      stableEnd   = stableEnd->next();

      stab_mm    -= mm;
      stab_gaps  -= gp;
      stab_bases -= bps;
      cum_mm     += mm;

      terCol      = terCol->next();

      winLen++;
    } else {
      break;
    }
  }

  if (winLen > 1)
    return(winLen);

  return(0);
}




int32
abAbacus::refineWindow(abColumn     *bgnCol,
                       abColumn     *terCol) {

  abAbacusWork  *orig_abacus     = new abAbacusWork(this, bgnCol, terCol);
  int32          orig_columns    = 0;

  int32          orig_mm_score   = orig_abacus->scoreAbacus(orig_columns);

  abAbacusWork  *left_abacus     = orig_abacus->clone();
  int32          left_columns    = 0;
  abAbacusWork  *right_abacus    = orig_abacus->clone();
  int32          right_columns   = 0;

  int32          left_mm_score   = left_abacus->leftShift(left_columns);
  int32          right_mm_score  = right_abacus->rightShift(right_columns);

  // determine best score and apply abacus to real columns
  int32          orig_gap_score  = orig_abacus->affineScoreAbacus();
  int32          left_gap_score  = left_abacus->affineScoreAbacus();
  int32          right_gap_score = right_abacus->affineScoreAbacus();

  abAbacusWork  *best_abacus     = orig_abacus;
  int32          best_columns    = orig_columns;
  int32          best_gap_score  = orig_gap_score;
  int32          best_mm_score   = orig_mm_score;

  int32 orig_total_score  = orig_mm_score  + orig_columns  + orig_gap_score;
  int32 left_total_score  = left_mm_score  + left_columns  + left_gap_score;
  int32 right_total_score = right_mm_score + right_columns + right_gap_score;
  int32 best_total_score  = orig_total_score;

  int32 score_reduction   = 0;

  // Use the total score to refine the abacus
  if (left_total_score < orig_total_score || right_total_score < orig_total_score ) {
    if (left_total_score <= right_total_score ) {
      score_reduction += orig_total_score - left_total_score;
      //fprintf(stderr,"\nTry to apply LEFT abacus:\n");
      //ShowAbacus(left_abacus);
      best_abacus      = left_abacus;
      best_mm_score    = left_mm_score;
      best_columns     = left_columns;
      best_gap_score   = left_gap_score;
      best_total_score = left_total_score;
    } else {
      score_reduction += orig_total_score - right_total_score;
      //fprintf(stderr,"\nTry to apply RIGHT abacus:\n");
      //ShowAbacus(right_abacus);
      best_abacus      = right_abacus;
      best_mm_score    = right_mm_score;
      best_columns     = right_columns;
      best_gap_score   = right_gap_score;
      best_total_score = right_total_score;
    }
  }

  best_abacus->applyAbacus(this);

  delete orig_abacus;
  delete left_abacus;
  delete right_abacus;

  return(score_reduction);
}


//*********************************************************************************
// Abacus Refinement:
//   AbacusRefine contains the logic for sweeping through the multialignment,
//   and identifying candidate windows for refinement.
//   Each window is cast into an abacus, which is left and right shifted.
//   The best resulting base arrangement (if different from input) is then
//   applied to window of the MultiAlignment
//*********************************************************************************


//  AbacusRefine
//
//  ctgcns - from=0 to=-1 level=abAbacus_Smooth
//         - from=0 to=-1 level=abAbacus_Poly_X
//         - from=0 to=-1 level=abAbacus_Indel
//  utgcns - same
//
//  from,to are C-style.  Used to be INCLUSIVE, but never used anyway
//
int32
abAbacus::refine(abAbacusRefineLevel  level,
                 uint32               bgn,
                 uint32               end) {

#warning SKIPPING ALL REFINEMENTS
  return(0);

  if (end > numberOfColumns())
    end = numberOfColumns();

  abColumn *bgnCol = _columns[bgn];
  abColumn *endCol = _columns[end - 1];
  abColumn *terCol = NULL;

  int32    score_reduction = 0;

  while (bgnCol != endCol) {
    int32 window_width = 0;

    switch (level) {
      case abAbacus_Smooth:
        window_width = IdentifyWindow_Smooth(bgnCol, terCol);
        break;

      case abAbacus_Poly_X:
        window_width = IdentifyWindow_Poly_X(bgnCol, terCol);
        break;

      case abAbacus_Indel:
        window_width = IdentifyWindow_Indel(bgnCol, terCol);
        break;

      default:
        break;
    }


    //  If the window is too big, there's likely a polymorphism that won't respond well to abacus,
    //  so skip it.
    //
    if ((window_width > 0) &&
        (window_width < MAX_WINDOW_FOR_ABACUS_REFINE)) {

      //  Insert a gap column for maneuvering room if the window starts in the first.  This used to
      //  be logged, and BPW can't remember EVER seeing the message.

#if 0
      if (bgnCol->prev() == NULL) {
        abBeadID   firstbeadID = abacus->getBead(bgnCol->callID() )->downID();
        abBeadID   newbeadID   = abacus->appendGapBead(abacus->getBead(firstbeadID)->ident());

        fprintf(stderr, "abMultiAlign::refine()-- Adding gap bead "F_U32" after first bead "F_U32" to add abacus room for abutting left of multialignment\n",
                newbeadID.get(), firstbeadID.get());

        abacus->appendColumn(abacus->getBead(firstbeadID)->colIdx(), newbeadID);
      }
#endif

      //  Actually do the refinements.
      score_reduction += refineWindow(bgnCol, terCol);
    }

    //  Move to the column after the window we just examined.
    bgnCol = terCol;
  }

  //  WITH quality=1 make_v_list=1, all the rest defaults
  refreshColumns();
  recallBases(true);

  return(score_reduction);
}

