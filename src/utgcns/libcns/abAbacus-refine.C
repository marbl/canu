
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
 *    Brian P. Walenz beginning on 2014-DEC-30
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static char *rcsid = "$Id$";

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
    start_column = abColID();
    end_column   = abColID();
    rows         = 0;
    columns      = 0;
    window_width = 0;
    shift        = UNSHIFTED;
    beads        = NULL;
    calls        = NULL;
  };

  abAbacusWork(abAbacus        *abacus,
               abMultiAlignID   mid,
               abColID          from,
               abColID          end);

  ~abAbacusWork() {
    delete [] beads;
    delete [] calls;
  };


  void   reset(void) {
    for (uint32 j=0; j<columns; j++)
      calls[j] = 'n';
  }


  void  set(int32 i, int32 j, char c) {
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

  char  getChar(int32 i, int32 j) {
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

    fprintf(stderr, "start_column  %d\n", start_column.get());
    fprintf(stderr, "end_column    %d\n", end_column.get());
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


  void    GetBaseCount(abBaseCount &b) {
    b.clear();

    for (uint32 j=0; j<columns; j++)
      b.IncBaseCount(calls[j]);
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
  abColID           start_column;
  abColID           end_column;
  int32             rows;
  int32             columns;
  int32             window_width;
  ShiftStatus       shift;
  char             *beads;
  char             *calls;

  vector<uint32>   abacus_indices;
};



static
int
base2int(char b) {
  if (b == '-')             return 0;
  if (b == 'a' || b == 'A') return 1;
  if (b == 'c' || b == 'C') return 2;
  if (b == 'g' || b == 'G') return 3;
  if (b == 't' || b == 'T') return 4;
  if (b == 'n' || b == 'N') return 5;
  fprintf(stderr, "base2int b out of range");
  assert(0);
  return 255;
}


static
int
is_good_base(char b) {
  if (b == '-')             return 1;
  if (b == 'a' || b == 'A') return 1;
  if (b == 'c' || b == 'C') return 1;
  if (b == 'g' || b == 'G') return 1;
  if (b == 't' || b == 'T') return 1;
  if (b == 'n' || b == 'N') return 1;
  return 0;
}






abAbacusWork::abAbacusWork(abAbacus        *abacus,
                           abMultiAlignID   mid,
                           abColID          from,
                           abColID          end) {

  // from and end are ids of the first and last columns in the columnStore

  abMultiAlign    *ma = abacus->getMultiAlign(mid);


  //ResetIndex(abacus_indices,GetNumFragments(fragmentStore));
  abacus_indices.clear();
  for (uint32 ii=0; ii<abacus->numberOfSequences(); ii++)
    abacus_indices.push_back(0);
  //fprintf(stderr, "CREATE ABACUS abacus_indices size %d\n", abacus_indices.size());


  //  Count the number of rows and columns.  Also find the last column.
  //
  uint32        origNumColumns = 1;  //  We don't count the last column in the loop below.
  abColumn     *column  = abacus->getColumn(from);
  abColumn     *last    = abacus->getColumn(from);

  //fprintf(stderr, "column_ids= ");

  for (; (last->nextID() != end) && (last->nextID().isValid()); last = abacus->getColumn(last->nextID())) {
    origNumColumns++;
    //fprintf(stderr, "%u,", last->ident().get());
  }

  //fprintf(stderr, "%u\n", last->ident().get());
  //fprintf(stderr, "origNumColumns %d\n", origNumColumns);

  // GD: this is where base calling code should be called

  //
  //  For every bead in the first (last) column, build a map from the seqIdx
  //  to a row.

  // a little fragment may sneak in, have to ensure the abacus has a
  // row for it.  The introduction of late- and mid column was done
  // to eliminate a problem with a degenerate alignment consistenting
  // of essentially one long poly run.  (encountered in unitig
  // 1618966 of the NOV'01 human vanilla assembly) it happened that a
  // little fragment was caught even between the mid_column and end
  // column, so it's index wasn't in the index set...  which causes a
  // "SetAbacus" beyond row range error.  putting in two mid-columns
  // will hopefully catch all fragments in the abacus range.
  //

  //  Macaque, using overlap based trimming, needed mid_column points
  //  at rather small intervals to pass.  Even without OBT, macaque
  //  needed another point at 63 to pass.
  //
  //  This change was tested on macaque, and did not change the
  //  results (except for allowing one partition to finish....).  You
  //  can revert to the original behavior by undef'ing the following.

  //  For canu, we add new columns every 512 bases.  The original logic was convoluted and probably
  //  broken -- it added a new column every MAX_READ_LEN/2 (every 1k), which conflicts with the
  //  above comment.

  abColBeadIterator *bi         = NULL;
  uint32             rowcounter = 0;

  //  For the first column

  bi = abacus->createColBeadIterator(column->ident());

  for (abBeadID bid=bi->next(); bid.isValid(); bid=bi->next()) {
    uint32 si = abacus->getBead(bid)->seqIdx().get();

    abacus_indices[si] = ++rowcounter;
  }

  delete bi;

  //  For the last column

  bi = abacus->createColBeadIterator(last->ident());

  for (abBeadID bid=bi->next(); bid.isValid(); bid=bi->next()) {
    uint32 si = abacus->getBead(bid)->seqIdx().get();

    if (abacus_indices[si] == 0)
      abacus_indices[si] = ++rowcounter;
  }

  delete bi;

  //  For a bunch of intermediate columns

  uint32        colIdx  = 0;
  abColumn     *c = abacus->getColumn(from);

  for (; ((c->nextID() != end) && (c->nextID().isValid())); c = abacus->getColumn(c->nextID())) {
    colIdx++;

    if ((colIdx % 512) != 0)
      continue;

    bi = abacus->createColBeadIterator(c->ident());

    for (abBeadID bid = bi->next(); bid.isValid(); bid = bi->next()) {
      abBead *bead = abacus->getBead(bid);
      uint32  si   = bead->seqIdx().get();

      if (abacus_indices[si] == 0)
        abacus_indices[si] = ++rowcounter;
    }
    delete bi;
  }


  //  Fill out or class members.

  start_column  = from;
  end_column    = last->ident();
  rows          = rowcounter;
  window_width  = origNumColumns;
  columns       = 3 * origNumColumns;
  shift         = UNSHIFTED;
  beads         = new char [rows * (columns + 2)];    // two extra gap columns, plus "null" borders
  calls         = new char [columns];


  for (uint32 i=0; i<rows * (columns+2); i++)
    beads[i] = 'n'; // initialize to "null" code



  //  Fill the center third of abacus with chars from the columns

  uint32  colcounter = 0;

  while ((column->ident() != end) && (column->ident().isValid())) {
    abColBeadIterator *bi = abacus->createColBeadIterator(column->ident());

    for (abBeadID bid = bi->next(); bid.isValid(); bid=bi->next()) {
      abBead *bead = abacus->getBead(bid);

      set(abacus_indices[bead->seqIdx().get()] - 1,
          colcounter + origNumColumns,
          abacus->getBase(bead->baseIdx()));
    }

    delete bi;

    colcounter++;

    column = abacus->getColumn(column->nextID());
  }

  //  Clear the border columns

  for (uint32 i=0; i<rows; i++) {
    for (uint32 j=0; j<origNumColumns; j++)
      set(i, j, '-');

    for (uint32 j=2 * origNumColumns; j<columns; j++)
      set(i, j, '-');
  }

  reset();
}






// cols is the number of "good" (non-null) columns found
// GD: This function counts the total number of bases which
//   - are different from column's "consensus" call and
//   - are not 'n'
//
uint32
abAbacusWork::scoreAbacus(int32 &cols) {
  abBaseCount *counts = new abBaseCount [columns];

  //fprintf(stderr, "scoreAbacus()--\n");

#warning "this doesn't need to allocate an array of abBaseCount"

  uint32 score = 0;

  cols = 0;

  for (int32 i=0; i<rows; i++) {
    for (int32 j=0; j<columns; j++) {
      char b = getChar(i, j);

      if ((b == '-' ) && (j > 0) && (j < columns - 1) &&
          ((getChar(i, j-1) == 'n')  ||
           (getChar(i, j+1) == 'n')))
        b = 'n';

      counts[j].IncBaseCount(b);
    }
  }

  // now, for each column, generate the majority call
  for (int32 j=0; j<columns; j++) {
    if (counts[j].GetBaseCount('-') + counts[j].GetBaseCount('n') == counts[j].GetDepth()) {
      // null (all-gap) column. Flag with an 'n' basecall
      calls[j] = 'n';

    } else {
      cols++;

      calls[j] = counts[j].GetMaxBaseCountBase(false);

      // and then tally edit score
      score += counts[j].GetDepth() - counts[j].GetBaseCount(calls[j]) - counts[j].GetBaseCount('n');
    }
  }

  delete [] counts;

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
      if (getChar(i, j) != '-' ) {
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
        curr = getChar(i,j);
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
        curr = getChar(i,j);
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
              curr = getChar(i,j);
              next = getChar(i,j+1);
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
                curr = getChar(i,j  );
                next = getChar(i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (next != '-' && next != 'n' )
                  {
                    set(i, j  , next);
                    set(i, j+1, curr);
                  }
              }
              // The entire j+1-th column now contains only gaps or n's
              // Remove it by shifting all the subsequent columns
              // one position to the left
              for (int32 i=0;i<rows;i++)
                {
                  curr = getChar(i,j  );
                  next = getChar(i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j+1; k<last_non_null; k++)
                    {
                      next= getChar(i,k+1);
                      set(i, k, next);
                    }
                  set(i, last_non_null, '-');
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
              curr = getChar(i,j);
              next = getChar(i,j+1);
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
                curr = getChar(i,j  );
                next = getChar(i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (curr != '-' && curr != 'n' ) {
                  set(i, j  , next);
                  set(i, j+1, curr);
                }
              }
              // The entire j-th column contains gaps
              // Remove it by shifting all the previous columns
              // one position to the right
              for (int32 i=0;i<rows;i++)
                {
                  curr = getChar(i,j  );
                  next = getChar(i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j; k>first_non_null; k--)
                    {
                      prev = getChar(i,k-1);
                      set(i, k, prev);
                    }
                  set(i, first_non_null, '-');
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

      char  c = getChar(i, j);
      int32 ccol = j;

      if (c != '-' ) {
        //look to the left for a suitable placement
        // will be safe on left since abacus has 'n' border
        while (getChar(i, ccol-1) == '-')
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
            set(i, j, '-');
            set(i, pcol, c);
            break;
          }
        }

        if (getChar(i, j) != '-')
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
      char  c = getChar(i,j);
      int32 ccol = j;

      if (c != '-' ) {
        //look to the right for a suitable placement
        // will be safe on right since abacus has 'n' border
        while (getChar(i,ccol+1) == '-' )
          ccol++;
        // now, from ccol back down to j, look for column with matching call
        for (int32 pcol = ccol;pcol>j;pcol--) {
          char call = calls[pcol];
          if (call != 'n' && call != c && c != 'n' )
            continue;
          if (call == 'n')
            calls[pcol] = c;
          if (calls[pcol] == c || c == 'n' ) {
            set(i,j,'-');
            set(i,pcol,c);
            break;
          }
        }
        if (getChar(i,j) != '-' )
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
//  a) end gap moving left:
//
//      - > A > B > ... > C > X  becomes  A B ... C - X
//      \_________________^
//
//  b) non-gap moving left across only gap characters
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
  abColumn    *column;
  int32        columnCount=0;
  char         a_entry;

  abBeadID bid;  //  ALWAYS the id of bead
  abBeadID eid;  //  ALWAYS the id of exch

  abBead *bead = NULL;
  abBead *exch = NULL;

  char base;

  //fprintf(stderr, "applyAbacus()--  shift=%c start=%d width=%d\n", shift, start_column.get(), window_width);
  //show();

  if (shift == LEFT_SHIFT) {
    column = abacus->getColumn(start_column);
    assert(column != NULL);

    while (columnCount < window_width) {
      bid = abacus->getBead(column->callID())->downID();
      //fprintf(stderr, "0; bid=%d eid=%d\n", bid.get(), eid.get());
      // Update all beads in a given column

      while (bid.isValid()) {
        bead = abacus->getBead(bid);
        a_entry = getChar(abacus_indices[bead->seqIdx().get()] - 1, columnCount);

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

              abColID cid = column->ident();
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

      column = abacus->getColumn(column->nextID());
      columnCount++;
    }
  }


  if (shift == RIGHT_SHIFT) {
    column = abacus->getColumn(end_column);
    assert(column != NULL);

    while (columnCount < window_width) {
      bid = abacus->getBead(column->callID())->downID();

      while (bid.isValid()) {
        bead = abacus->getBead(bid);
        a_entry = getChar(abacus_indices[bead->seqIdx().get()] - 1, columns - columnCount - 1);

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

              abColID cid = column->ident();
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



//
//  In this case, we just look for a string of gaps in the consensus sequence
//
static
int32
IdentifyWindow_Smooth(abAbacus  *abacus,
                      abColumn *&start_column,
                      abColID   &stab_bgn) {

  stab_bgn = start_column->nextID();

  abColumn *stab = abacus->getColumn(stab_bgn);
  int32  win_length = 0;

#ifdef DEBUG_IDENTIFY_WINDOW
  fprintf(stderr, "identifyWindow()-- bgn=%d level=%d\n", stab_bgn.get(), level);
#endif

  if (abacus->getBase( start_column->callID() ) != '-')
    //  Consensus not a gap, nothing to do.
    return(win_length);

  //  Consensus is a gap.  Expand it to the maximum gap.

  while ((abacus->getBase( stab->callID()) == '-') &&
         (stab->nextID().isValid())) {
    stab_bgn = stab->nextID();
    stab     = abacus->getColumn(stab_bgn);

    win_length++;
  }
  
#ifdef DEBUG_IDENTIFY_WINDOW
  fprintf(stderr, "identifyWindow()-- gap at %d to %d  win_length=%d (return)\n",
          start_column->position(), stab->position(), win_length);
#endif

  return(win_length);
}




//
//  Here, we're looking for a string of the same character
//
static
int32
IdentifyWindow_Poly_X(abAbacus  *abacus,
                      abColumn *&start_column,
                      abColID   &stab_bgn) {

  stab_bgn = start_column->nextID();

  abColumn *stab = abacus->getColumn(stab_bgn);

#ifdef DEBUG_IDENTIFY_WINDOW
  fprintf(stderr, "identifyWindow()-- bgn=%d level=%d\n", stab_bgn.get(), level);
#endif

  int32 gap_count  = start_column->GetColumnBaseCount('-');
  int32 win_length = 1;

  char  poly  = abacus->getBase(start_column);
  char  cb    = abacus->getBase(stab);

  if (poly == '-')
    return(0);

  while ((cb == poly) || (cb == '-'))  {
    if (stab->nextID().isValid() == false)
      break;

    // move stab column ahead

    gap_count  += stab->GetColumnBaseCount('-');
    win_length += 1;

    stab_bgn    = stab->nextID();
    stab        = abacus->getColumn(stab_bgn);

    cb          = abacus->getBase(stab);
  }

  if (win_length <= 2)
    return(0);

  // capture trailing gap-called columns

  while (abacus->getBase(stab) == '-' )  {
    if (stab->GetMaxBaseCountBase(true) != poly )
      break;

    if (stab->nextID().isValid() == false)
      break;

    gap_count  += stab->GetColumnBaseCount('-');
    win_length += 1;

    stab_bgn    = stab->nextID();
    stab        = abacus->getColumn(stab_bgn);
  }

  // now that a poly run with trailing gaps is established, look for leading gaps

  abColumn *pre_start = start_column;

  while (pre_start->prevID().isValid()) {
    pre_start = abacus->getColumn(pre_start->prevID());
    cb        = abacus->getBase(pre_start);

    if ((cb != '-') && (cb != poly))
      break;

    start_column = pre_start;

    gap_count  += pre_start->GetColumnBaseCount('-');
    win_length += 1;
  }

  //fprintf(stderr,"POLYX candidate (%c) at column %d stab %d , width %d, gapcount %d\n",
  //        poly, start_column->position(), stab_bgn.get(), win_length, gap_count);

  if ((start_column->prevID().isValid() == true) &&
      (win_length > 2) &&
      (gap_count  > 0))
    return(win_length);

  return(0);
}





//
//  In this case, we look for a string mismatches, indicating a poor alignment region
//  which might benefit from Abacus refinement
//
//  heuristics:
//
//  > stable border on either side of window of width:  STABWIDTH
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
//  start_column
//

#define STABWIDTH  6

static
int32
IdentifyWindow_Indel(abAbacus  *abacus,
                     abColumn *&start_column,
                     abColID   &stab_bgn) {

  stab_bgn = start_column->nextID();

  abColumn *stab = abacus->getColumn(stab_bgn);

#ifdef DEBUG_IDENTIFY_WINDOW
  fprintf(stderr, "identifyWindow()-- bgn=%d level=%d\n", stab_bgn.get(), level);
#endif

  int32 cum_mm=0;
  int32 stab_mm=0;
  int32 stab_gaps=0;
  int32 stab_width=0;
  int32 stab_bases=0;
  abColumn *stab_end;
  int32 win_length = 0;

  cum_mm = start_column->mismatch();
  if (cum_mm > 0 && start_column->GetColumnBaseCount('-') > 0) {
    stab = start_column;
    stab = abacus->getColumn(start_column->nextID());
    stab_end = stab;
    while (stab_end->nextID().isValid() && stab_width < STABWIDTH) {
      stab_mm+=stab_end->mismatch();
      stab_gaps+=stab_end->GetColumnBaseCount('-');
      stab_bases+=stab_end->GetDepth();
      stab_end = abacus->getColumn(stab_end->nextID());
      stab_width++;
    }

    if (stab_bases == 0 )
      return(0);

    //  Floating point 'instability' here?
    while ((double)stab_mm/(double)stab_bases > 0.02 ||  //  CNS_SEQUENCING_ERROR_EST
           (double)stab_gaps/(double)stab_bases > .25  ){
      int32 mm=stab->mismatch();
      int32 gp=stab->GetColumnBaseCount('-');
      int32 bps=stab->GetDepth();
      // move stab column ahead
      if (stab_end->nextID().isValid() ) {
        stab_mm+=stab_end->mismatch();
        stab_bases+=stab_end->GetDepth();
        stab_gaps+=stab_end->GetColumnBaseCount('-');
        stab_end = abacus->getColumn(stab_end->nextID());
        stab_mm-=mm;
        stab_gaps-=gp;
        stab_bases-=bps;
        cum_mm+=mm;
        stab = abacus->getColumn(stab->nextID());
        win_length++;
      } else {
        break;
      }
    }
    stab_bgn = stab->ident();
  }

  if (win_length > 1)
    return(win_length);

  return(0);
}




int32
abMultiAlign::refineWindow(abAbacus     *abacus,
                           abColumn     *start_column,
                           abColID       stab_bgn) {

  abBaseCount    abacus_count;
  abAbacusWork  *orig_abacus     = new abAbacusWork(abacus, ident(), start_column->ident(), stab_bgn);
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
      left_abacus->GetBaseCount(abacus_count);
      best_abacus      = left_abacus;
      best_mm_score    = left_mm_score;
      best_columns     = left_columns;
      best_gap_score   = left_gap_score;
      best_total_score = left_total_score;
    } else {
      score_reduction += orig_total_score - right_total_score;
      //fprintf(stderr,"\nTry to apply RIGHT abacus:\n");
      //ShowAbacus(right_abacus);
      right_abacus->GetBaseCount(abacus_count);
      best_abacus      = right_abacus;
      best_mm_score    = right_mm_score;
      best_columns     = right_columns;
      best_gap_score   = right_gap_score;
      best_total_score = right_total_score;
    }
  }

  best_abacus->applyAbacus(abacus);

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
abMultiAlign::refine(abAbacus            *abacus,
                     abAbacusRefineLevel  level,
                     uint32               from,    // from and to are in ma's column coordinates
                     uint32               to) {

  if (to > length())
    to = length();

  abColID sid = columnList[from];      // id of the starting column
  abColID eid = columnList[to - 1];    // id of the ending column

  assert(sid == firstColumn());
  assert(eid == lastColumn());

  abColumn *start_column = abacus->getColumn(sid);
  abColID   stab_bgn;

  uint32    score_reduction = 0;


  while (start_column->ident() != eid) {
    int32 window_width = 0;

    switch (level) {
      case abAbacus_Smooth:
        IdentifyWindow_Smooth(abacus, start_column, stab_bgn);
        break;
      case abAbacus_Poly_X:
        IdentifyWindow_Poly_X(abacus, start_column, stab_bgn);
        break;
      case abAbacus_Indel:
        IdentifyWindow_Indel(abacus, start_column, stab_bgn);
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
      //
      if (start_column->prevID().isValid() == false) {
        abBeadID   firstbeadID = abacus->getBead(start_column->callID() )->downID();
        abBeadID   newbeadID   = abacus->appendGapBead(abacus->getBead(firstbeadID)->ident());

        fprintf(stderr, "abMultiAlign::refine()-- Adding gap bead "F_U32" after first bead "F_U32" to add abacus room for abutting left of multialignment\n",
                newbeadID.get(), firstbeadID.get());

        abacus->appendColumn(abacus->getBead(firstbeadID)->colIdx(), newbeadID);
      }

      //  Actually do the refinements.
      score_reduction += refineWindow(abacus, start_column, stab_bgn);
    }

    //  Move to the column after the window we just examined.
    start_column = abacus->getColumn(stab_bgn);
  }

  //  WITH quality=1 make_v_list=1, all the rest defaults
  abacus->refreshMultiAlign(ident(), true, true);

  return(score_reduction);
}
