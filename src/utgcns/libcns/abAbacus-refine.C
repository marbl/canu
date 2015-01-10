
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static char *rcsid = "$Id$";

#include "abAbacus.H"



#define MAX_WINDOW_FOR_ABACUS_REFINE      100

//  Needs to go...
#define  MAX_MID_COLUMN_NUM 100



enum ShiftStatus {
  LEFT_SHIFT  = (int) 'L', // Left Shifted
  RIGHT_SHIFT = (int) 'R', // Right Shifted
  UNSHIFTED   = (int) 'U', // Unshifted
  MIXED_SHIFT = (int) 'M'  // shifted in different directions
};



class abAbacusWork {
public:
  abAbacusWork() {
    start_column = abColID();
    end_column   = abColID();
    rows         = 0;
    columns      = 0;
    window_width = 0;
    shift        = MIXED_SHIFT;
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
    if ((i == -1) || (i >= rows))
      fprintf(stderr, "abAbacusWork::set()--  attempt to write beyond row range: i=%d j=%d rows=%d\n", i, j, rows);
    if ((j == -1) || (j >= columns))
      fprintf(stderr, "abAbacusWork::set()--  attempt to write beyond column range: i=%d j=%d columns=%d\n", i, j, columns);

    assert(i != -1);
    assert(i < rows);
    assert(j != -1);
    assert(j < columns);

    beads[i * (columns + 2) + j + 1] = c;
  };

  char *getPtr(int32 i, int32 j) {
    return(beads + i * (columns + 2) + j + 1);
  };
  char  getChar(int32 i, int32 j) {
    return(beads[i * (columns + 2) + j + 1]);
  };


  abAbacusWork *clone(void) {
    abAbacusWork *clone = new abAbacusWork;

    clone->start_column = start_column;
    clone->end_column   = end_column;
    clone->window_width = window_width;
    clone->rows         = rows;
    clone->columns      = columns;
    clone->shift        = shift;
    clone->beads        = new char [rows * (columns+2)];
    clone->calls        = new char [columns];

    memcpy(clone->beads, beads, (rows * (columns+2)) * sizeof(char));
    memcpy(clone->calls, calls, (columns)            * sizeof(char));

    return(clone);
  };


  void   show(void) {
    char form[10];

    sprintf(form, "%%%d.%ds\n", columns, columns);
    fprintf(stderr, "\nstart column: %d\n", start_column.get());
    for (uint32 i=0; i<rows; i++)
      fprintf(stderr, form, getPtr(i,0));
    fprintf(stderr,"\n");
    fprintf(stderr, form, calls);
  };


  void    GetBaseCount(abBaseCount &b) {
    b.clear();

    for (uint32 j=0; j<columns; j++)
      b.IncBaseCount(calls[j]);
  };


  uint32   scoreAbacus(int32 &cols);
  uint32   affineScoreAbacus(void);

  int32    merge(int32 merge_dir);

  void     refineOrigAbacus(abVarRegion &vreg);

  uint32   leftShift(abVarRegion &vreg, int32 &lcols);
  uint32   rightShift(abVarRegion &vreg, int32 &rcols);
  uint32   mixedShift(int32 &mcols, abVarRegion  vreg, int32 lpos, int32 rpos, char *tmpl, int32 long_allele, int32 short_allele);

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



  //  Count the number of rows and columns.  Also find the last column.
  //
  uint32        origNumColumns = 0;
  abColumn     *column  = abacus->getColumn(from);
  abColumn     *last    = abacus->getColumn(from);

  for (; (last->nextID() != end) && (last->nextID().isValid()); last = abacus->getColumn(last->nextID()))
    origNumColumns++;

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

  //  For ca3g, we add new columns every 512 bases.  The original logic was convoluted and probably
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
      cols = cols + 1;
      calls[j] = counts[j].GetMaxBaseCount(false);

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
#ifdef DEBUG_ABACUS
  fprintf(stderr, "columns=%d first_non_null = %d last_non_null= %d\n",
          columns, first_non_null, last_non_null);
#endif
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
#if 0
  fprintf(stderr, "Columns merged=%d\n", columns_merged);
#endif
  return columns_merged;
}


/* Refine the original abacus:
 *    - make sure that reads that belong to the same (confirmed) allele
 *      are aligned identically
 *    - this alignment will be dictated by the last read in the allele
 */

void
abAbacusWork::refineOrigAbacus(abVarRegion &vreg) {

  reset();

  for (uint32 j=window_width; j < 2*window_width; j++) {

    // Look through confirmed alleles only
    for (uint32 k=0; k<vreg.nca; k++) {

      // Line number corresponding to the last read in allele
      int32 il = vreg.alleles[k].read_ids[vreg.alleles[k].num_reads-1];
      char  c  = getChar(il, j);

      for (int32 l=0; l<vreg.alleles[k].num_reads-1; l++)
        set(vreg.alleles[k].read_ids[l], j, c);
    }
  }
}



uint32
abAbacusWork::leftShift(abVarRegion &vreg, int32 &lcols) {
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;

  reset();

#if 0
  fprintf(stderr, "Abacus region:\n");
  fprintf(stderr, "nr= %d na= %d nca= %d\n", vreg.nr, vreg.na, vreg.nca);
  fprintf(stderr, "Order of left-shifting alleles:\n");
  for (i=0; i<vreg.na; i++)
    {
      fprintf(stderr, "Allele %d uglen=%d weight= %d\n", i,
              vreg.alleles[i].uglen, vreg.alleles[i].weight);
      fprintf(stderr, "   Reads:\n");
      for (j=0; j<vreg.alleles[i].num_reads; j++)
        {
          int32 k, read_id = vreg.alleles[i].read_ids[j];
          int32 len = vreg.end-vreg.beg;
          fprintf(stderr, "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(stderr, "%c", vreg.reads[read_id].bases[k]);
          fprintf(stderr, "\n");
        }
    }
#endif

  for (j=window_width; j<2*window_width; j++) {
      for (k=0; k<vreg.na; k++) {
          for (l=0; l<vreg.alleles[k].num_reads; l++) {
              i = vreg.alleles[k].read_ids[l];
              c = getChar( i, j);
              ccol = j;
              if (c != '-' ) {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (getChar( i, ccol-1) == '-')
                    ccol--;

                  // from ccol back up to j, look for column with matching call
                  for (pcol = ccol; pcol<j; pcol++) {
                      call = calls[pcol];
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

                  if (getChar( i, j) != '-')
                    calls[j] = c;
                }
            }
        }
    }

#if 0
  fprintf(stderr, "Test calls=\n");
  for (j=0;j<columns;j++)
    fprintf(stderr, "%c", calls[j]);
  fprintf(stderr, "\n");
#endif

#ifdef DEBUG_ABACUS
  fprintf(stderr, "Abacus after LeftShift before Merge:\n");
  show();
#endif

  merge(-1);

#if 0
  fprintf(stderr, "Abacus after Merge:\n");
  show();
#endif

  shift = LEFT_SHIFT;

  return(scoreAbacus(lcols));
}


uint32
abAbacusWork::rightShift(abVarRegion &vreg, int32 &rcols) {
  // rcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;

  reset();

#if 0
  fprintf(stderr, "Abacus region:\n");
  fprintf(stderr, "nr= %d na= %d nca= %d\n", vreg.nr, vreg.na, vreg.nca);
  fprintf(stderr, "Order of left-shifting alleles:\n");
  for (i=0; i<vreg.na; i++)
    {
      fprintf(stderr, "Allele %d uglen=%d weight= %d\n", i,
              vreg.alleles[i].uglen, vreg.alleles[i].weight);
      fprintf(stderr, "   Reads:\n");
      for (j=0; j<vreg.alleles[i].num_reads; j++)
        {
          int32 k, read_id = vreg.alleles[i].read_ids[j];
          int32 len = vreg.end-vreg.beg;
          fprintf(stderr, "    %d   ", read_id);
          for (k=0; k<len; k++)
            fprintf(stderr, "%c", vreg.reads[read_id].bases[k]);
          fprintf(stderr, "\n");
        }
    }
#endif

  for (j=2*window_width-1;j>window_width-1;j--)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              c = getChar(i,j);
              ccol = j;
              if (c != '-' )
                {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while (getChar(i,ccol+1) == '-' )
                    ccol++;
                  // now, from ccol back down to j, look for column with matching call
                  for (pcol = ccol;pcol>j;pcol--)
                    {
                      call = calls[pcol];
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
    }

#ifdef DEBUG_ABACUS
  fprintf(stderr, "Abacus after RightShift before Merge:\n");
  show();
#endif

  merge(1);

  shift = RIGHT_SHIFT;

  return(scoreAbacus(rcols));
}

uint32
abAbacusWork::mixedShift(int32 &mcols, abVarRegion  vreg, int32 lpos, int32 rpos, char *tmpl, int32 long_allele, int32 short_allele) {
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  int32 window_beg, window_end;
  int32 shift =0;

  reset();

  if (shift == LEFT_SHIFT) {
    window_beg = 0;
    window_end = window_width;

  } else if (shift == UNSHIFTED) {
    window_beg = window_width;
    window_end = 2* window_width;

  } else {
    window_beg = 2*window_width;
    window_end = 3*window_width;
  }

  /* Populate calls */
  for (j=window_beg; j<window_end; j++)
    calls[j] = tmpl[j];

  /* Perform left shift */
  for (j=window_beg;j<=MIN(window_end, lpos);j++) {
      for (k=0; k<vreg.na; k++) {
          for (l=0; l<vreg.alleles[k].num_reads; l++) {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = getChar(i,j);
              ccol = j;
              if (c != '-' ) {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (( getChar(i,ccol-1) == '-' ) &&
                         (ccol > window_beg)) {
                    ccol--;
                  }
                  // now, from ccol back up to j, look for column with matching call
                  for (pcol = ccol;pcol<j;pcol++) {
                    call = calls[pcol];
                    if (call != 'n' && call != c && c != 'n')
                      // GD: consensus in a column == '-' ?
                      continue;

                    if (call == 'n') {
                      // GD: found the leftmost column with non-gap consensus =>
                      //     reset it consensus "dynamically" to the current base
                      //     Potential problem: this code is biased  in the sense that
                      //     the result will generally depend on the order in which
                      //     reads i(or rows) are processed
                      calls[pcol] = c;
                    }
                    if (calls[pcol] == c || c == 'n') {
                      // swap bases in columns pcol and j of row i
                      set(i,j,'-');
                      set(i,pcol,c);
                      break;
                    }
                  }
                  if (getChar(i,j) != '-' ) {
                    calls[j] = c;
                  }
                }
            }
        }
    }
#if 0
  fprintf(stderr, "In MixedShift: window_beg=%d lpos=%d rpos=%d  window_end=%d\n",
          window_beg, lpos, rpos, window_end);
  fprintf(stderr, "Abacus calls=\n");
  for (i=window_beg; i<window_end; i++)
    fprintf(stderr, "%c", calls[i]);
  fprintf(stderr, "\n");
#endif

  /* Perform right shift */
  for (j=window_end-1;j>(rpos>0?rpos:window_end);j--) {
      for (k=0; k<vreg.na; k++) {
          for (l=0; l<vreg.alleles[k].num_reads; l++) {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
#if 0
              fprintf(stderr, "i=%d vreg.alleles[i]=%d short_allele=%d\n", i, vreg.alleles[i], short_allele);
#endif
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = getChar(i,j);
              ccol = j;
              if (c != '-' ) {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while (( getChar(i,ccol+1) == '-') &&
                         (ccol+1<window_end) )
                    ccol++;
#if 0
                  fprintf(stderr, "ccol=%d\n", ccol);
#endif
                  // now, from ccol back down to j, look for column with matching call
                  for (pcol = ccol;pcol>j;pcol--) {
                    call = calls[pcol];
#if 0
                    fprintf(stderr, "i=%d j=%d c=%c pcol=%d call=%d \n", i, j, c, pcol, call);
#endif
                    if (call != 'n' && call != c && c != 'n' ) {
                      continue;
                    }

                    if (call == 'n') {
                      calls[pcol] = c;
                    }
#if 0
                    fprintf(stderr, "calls=%c c=%c\n", calls, c);
#endif
                    if (calls[pcol] == c || c == 'n' ) {
#if 0
                      fprintf(stderr, "Swapping elements (%d, %d)=%c  and (%d, %d)='-'\n",
                              i, j, c, i, pcol);
#endif
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
    }
  // MergeAbacus(abacus, 1);

  shift = MIXED_SHIFT;

  return(scoreAbacus(mcols));
}



static
void
leftEndShiftBead(abAbacus *abacus, abBeadID bid, abBeadID eid) {
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

  abBead *shift = abacus->getBead(eid);
  abBeadID aid = (abacus->getBead(bid))->prevID();

  assert(shift != NULL);

  if (abacus->getBase(shift->baseIdx()) != '-' ) {
    // assume first and internal characters are gaps
    abacus->lateralExchangeBead(bid, eid);
  } else {
    while (shift->prevID() != aid ) {
      abacus->lateralExchangeBead(shift->prevID(), shift->ident());
    }
  }
}

static
void
rightEndShiftBead(abAbacus *abacus, abBeadID bid, abBeadID eid) {
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

  abBead *shift = abacus->getBead(bid);
  abBeadID aid = abacus->getBead(eid)->nextID();
  abBeadID rid;
  assert(shift != NULL);
  if (abacus->getBase(shift->baseIdx()) != '-' ) {
    // assume last and internal characters are gaps
    abacus->lateralExchangeBead(bid, eid);
  } else {
    rid = shift->nextID();
    while (shift->nextID() != aid ) {
      abacus->lateralExchangeBead(shift->ident(), shift->nextID());
    }
  }
}




void
abAbacusWork::applyAbacus(abAbacus *abacus) {
  abColumn    *column;
  int32        columns=0;
  char         a_entry;
  double       fict_var;   // variation is a column
  abVarRegion  vreg;

  abBeadID bid;  //  ALWAYS the id of bead
  abBeadID eid;  //  ALWAYS the id of exch

  abBead *bead = NULL;
  abBead *exch = NULL;

  char base;

  vreg.nr = 0;

  if (shift == LEFT_SHIFT) {
    column = abacus->getColumn(start_column);
    assert(column != NULL);

    while (columns < window_width) {
      bid = abacus->getBead(column->callID())->downID();
#ifdef DEBUG_APPLYABACUS
      fprintf(stderr, "0; bid=%d eid=%d\n", bid, eid);
#endif

      // Update all beads in a given column

      while (bid.isValid()) {
        bead = abacus->getBead(bid);
        a_entry = getChar(abacus_indices[bead->seqIdx().get()] - 1, columns);

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr, "a_entry=%c bead=%c\n", a_entry, abacus->getBase(bead->baseIdx()));
#endif

        if (a_entry == 'n') {
          eid  = bead->upID();
          exch = abacus->getBead(eid);
#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "1; bid=%d eid=%d\n", bid, eid);
#endif
          abacus->unalignTrailingGapBeads(bid);

        } else if (a_entry != abacus->getBase(bead->baseIdx())) {
          //  Look for matching bead in frag and exchange
          eid  = bead->ident();
          exch = abacus->getBead(eid);

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "2; bid=%d eid=%d\n", bid, eid);
#endif

          if (NULL == exch) {
            eid  = abacus->appendGapBead(bead->ident());
            bead = abacus->getBead(bid);
            abacus->alignBeadToColumn(abacus->getColumn(bead->colIdx())->nextID(),eid, "ApplyAbacus(1)");
            exch = abacus->getBead(eid);
          }

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "3; bid=%d eid=%d\n", bid, eid);
#endif

          while (a_entry != abacus->getBase(exch->baseIdx())) {
            abBeadID eidp = exch->nextID();

            if (exch->nextID().isInvalid()) {
              eidp = abacus->appendGapBead(exch->ident());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);
              abacus->alignBeadToColumn(abacus->getColumn(exch->colIdx())->nextID(),eidp, "ApplyAbacus(2)");
#ifdef DEBUG_APPLYABACUS
              fprintf(stderr, "4; bid=%d eid=%d\n", bid, eid);
#endif

            } else if (exch->colIdx() == end_column) {
              eidp = abacus->appendGapBead(exch->ident());
              bead = abacus->getBead(bid);
              exch = abacus->getBead(eid);
#ifdef DEBUG_APPLYABACUS
              fprintf(stderr, "5; bid=%d eid=%d\n", bid, eid);
#endif

              abColID cid = column->ident();
              abacus->appendColumn(exch->colIdx(),eidp);
              column = abacus->getColumn(cid);
            }

#ifdef DEBUG_APPLYABACUS
            fprintf(stderr, "6; bid=%d eid=%d b col/frg=%d/%d e_col/frg=%d/%d\n",
                    bid, eid,
                    bead->colIdx(), bead->frag_index,
                    exch->colIdx(), exch->frag_index);
#endif

            eid  = eidp;
            exch = abacus->getBead(eid);
          }

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr,"LeftShifting bead %d (%c) with bead %d (%c).\n",
                  bid, abacus->getBase(bead->baseIdx()),
                  eid, abacus->getBase(exch->baseIdx()));
#endif

          leftEndShiftBead(abacus, bid, eid);
        } else {
          // no exchange necessary;
          eid  = bid;
          exch = bead;

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "7; bid=%d eid=%d\n", bid, eid);
#endif
        }

        bid  = exch->downID();
        bead = NULL;

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr,"New bid is %d (%c), from %d down\n",
                bid, (bid.isValid()) ? abacus->getBase(abacus->getBead(bid)->baseIdx()) : 'n',
                eid);
#endif
      }

      //  End of update; call base now.

      base = abacus->baseCall(vreg, column->ident(), true, fict_var, -1, 0, 0, 11);

      column = abacus->getColumn(column->nextID());
      columns++;
    }
  }


  if (shift == RIGHT_SHIFT) {
    column = abacus->getColumn(end_column);
    assert(column != NULL);

    while (columns<window_width) {
      bid = abacus->getBead(column->callID())->downID();

      while (bid.isValid()) {
        bead = abacus->getBead(bid);
        a_entry = getChar(abacus_indices[bead->seqIdx().get()] - 1, columns-columns-1);

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

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr,"RightShifting bead %d (%c) with bead %d (%c).\n",
                  eid, abacus->getBase(exch->baseIdx()),
                  bid, abacus->getBase(bead->baseIdx()));
#endif

          rightEndShiftBead(abacus, eid, bid);
        } else {
          eid  = bid;
          exch = bead; // no exchange necessary;
        }

        bid  = exch->downID();
        bead = NULL;

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr,"New bid is %d (%c), from %d down\n",
                bid, (bid>-1)?abacus->getBase(abacus->getBead(bid)->baseIdx()):'n',
                eid);
#endif
      }

      base = abacus->baseCall(vreg, column->ident(), true, fict_var, -1, 0, 0, 11);
      column = abacus->getColumn(column->prevID());
      columns++;
    }
  }
}

static
int32
IdentifyWindow(abAbacus               *abacus,
               abColumn              *&start_column,
               abColID                &stab_bgn,
               abAbacusRefineLevel     level) {
  int32   rc = 0;

#if 0
  abColumn *stab;
  abColumn *pre_start;
  char poly;
#endif

  uint32 win_length = 1;
  uint32 gap_count  = 0;

  stab_bgn = start_column->nextID();

  abColumn *stab = abacus->getColumn(stab_bgn);

  switch (level) {
  case abAbacus_Smooth:
    {
    // in this case, we just look for a string of gaps in the consensus sequence
    if (abacus->getBase( abacus->getBead(start_column->callID())->baseIdx() ) != '-' ) break;
    // here, there's a '-' in the consensus sequence, see if it expands
    while( abacus->getBase( abacus->getBead(stab->callID())->baseIdx()) == '-' )  {
      // move stab column ahead
      if (stab->nextID().isValid() ) {
        stab_bgn = stab->nextID();
        stab = abacus->getColumn(stab_bgn);
        win_length++;
      } else {
        break;
      }
    }
    if (win_length > 1)
      rc = win_length;
    }
    break;

  case abAbacus_Poly_X:
    {
    // here, we're looking for a string of the same character
    gap_count  =  start_column->GetColumnBaseCount('-');
    char poly  =  abacus->getBase(abacus->getBead(start_column->callID())->baseIdx());

    if (poly != '-' ) {
      char cb;

      while( (cb = abacus->getBase(abacus->getBead(stab->callID())->baseIdx())) == poly || cb == '-' )  {
        // move stab column ahead
        if (stab->nextID().isValid() ) {
          stab_bgn = stab->nextID();
          gap_count+=stab->GetColumnBaseCount('-');
          stab = abacus->getColumn(stab_bgn);
          win_length++;
        } else {
          break;
        }
      }
      // capture trailing gap-called columns
      if (win_length > 2 ) {
        while( abacus->getBase(abacus->getBead(stab->callID())->baseIdx()) == '-' )  {
          if (stab->GetMaxBaseCount(1) != poly ) break;
          if (stab->nextID().isValid() ) {
            stab_bgn = stab->nextID();
            gap_count+=stab->GetColumnBaseCount('-');
            stab = abacus->getColumn(stab_bgn);
            win_length++;
          } else {
            break;
          }
        }
        // now that a poly run with trailing gaps is established, look for leading gaps
        abColumn *pre_start = start_column;
        while (pre_start->prevID().isValid() ) {
          char cb;
          pre_start = abacus->getColumn(pre_start->prevID());
          if ((cb = abacus->getBase(abacus->getBead(pre_start->callID())->baseIdx())) != '-' && cb != poly ) break;
          start_column = pre_start;
          gap_count+=pre_start->GetColumnBaseCount('-');
          win_length++;
        }
      } else {
        break;
      }
    }
    if (start_column->prevID().isValid() && win_length > 2 && gap_count > 0) {
      //fprintf(stderr,"POLYX candidate (%c) at column %d, width %d, gapcount %d\n", poly,start_column->ma_index,win_length,gap_count);
      rc = win_length;
    }
    }
    break;
  case abAbacus_Indel:
    {
    /*
      in this case, we look for a string mismatches, indicating a poor alignment region
      which might benefit from Abacus refinement
      heuristics:
      > stable border on either side of window of width:  STABWIDTH
      > fewer than STABMISMATCH in stable border

      _              __              ___
      SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
      SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
      SSSSS SSSSS => SSSSS .SSSS+ => SSSSS  .SSSS+
      SSSSS SSSSS    SSSSS .SSSS+    SSSSS  .SSSS+
      SSSSS_SSSSS    SSSSS_.SSSS+    SSSSS__.SSSS+
      |               \               \
      |\_______________\_______________\______ growing 'gappy' window
      start_column
    */
    {
      int32 cum_mm=0;
      int32 stab_mm=0;
      int32 stab_gaps=0;
      int32 stab_width=0;
      int32 stab_bases=0;
      abColumn *stab_end;

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
        if (stab_bases == 0 ) break;
        //  Floating point 'instability' here?
        while( (double)stab_mm/(double)stab_bases > 0.02 ||  //  CNS_SEQUENCING_ERROR_EST
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
      if (win_length > 1 ) rc = win_length;
    }
    }
    break;
  default:
    break;
  }

  return(rc);
}



static
void
ShowCalls(abAbacusWork *abacus) {
  for (int32 j=0;j<abacus->columns;j++)
    fprintf(stderr, "%c", abacus->calls[j]);
  fprintf(stderr, "\n");
}



static
void
getConsensusForAbacus(abVarRegion     &vreg,
                      abVarRead       *reads,
                      abAbacusWork    *abacus,
                      char            *consensus[2]) {

  char bases[CNS_NALPHABET] = {'-', 'A', 'C', 'G', 'T', 'N'};

  // Allocate memory for consensus

  consensus[0] = new char [3 * abacus->window_width];
  consensus[1] = new char [3 * abacus->window_width];

  memset(consensus[0], '-', sizeof(char) * 3 * abacus->window_width);
  memset(consensus[1], '-', sizeof(char) * 3 * abacus->window_width);

  // Call consensus
  for (int32 i=0; i < 3*abacus->window_width; i++) {
    int32 bcount0[CNS_NALPHABET] = {0};
    int32 bcount1[CNS_NALPHABET] = {0};
    int32 best_count0=0, second_best_count0=0;
    int32 best_count1=0, second_best_count1=0;
    char cbase0=0, cbase1=0;
    for (int32 j=0; j<abacus->rows; j++) {

#if 0
      fprintf(stderr, " reads[%d][%d]= %c\n", j, i, reads[j].bases[i]);
#endif

      if (is_good_base(reads[j].bases[i])) {
        if   (vreg.alleles[j].id == 0)
          bcount0[base2int(reads[j].bases[i])]++;
        else
          bcount1[base2int(reads[j].bases[i])]++;
      }
    }
    for (int32 j=0; j<CNS_NALPHABET; j++) {
      if (best_count0 < bcount0[j]) {
        second_best_count0 = best_count0;
        best_count0 = bcount0[j];
        cbase0 = bases[j];
      } else if ( best_count0 >= bcount0[j] && second_best_count0 <  bcount0[j]) {
        second_best_count0  = bcount0[j];
      }
    }
    for (int32 j=0; j<CNS_NALPHABET; j++) {
      if (best_count1 < bcount1[j]) {
        second_best_count1 = best_count1;
        best_count1 = bcount1[j];
        cbase1 = bases[j];
      } else if ( best_count1 >= bcount1[j] && second_best_count1 <  bcount1[j]) {
        second_best_count1  = bcount1[j];
      }
    }

    if (best_count0 == second_best_count0)
      consensus[0][i] = 'N';
    else
      consensus[0][i] = cbase0;

    if (best_count1 == second_best_count1)
      consensus[1][i] = 'N';
    else
      consensus[1][i] = cbase1;
  }
}

/* Create ungapped consensus sequences and map them to gapped consensus sequences */
static
void
mapConsensus(int32 *imap[2], char *consensus[2], char *ugconsensus[2], int32 len, int32 uglen[2]) {

  uglen[0] = 0;
  uglen[1] = 0;

  ugconsensus[0] = new char [len];
  ugconsensus[1] = new char [len];

  imap[0] = new int32 [len];
  imap[1] = new int32 [len];

  for (int32 i=0; i<2; i++) {
    for (int32 j=0; j<len; j++)
      imap[i][j] = j;

    int32 k=0;

    for (int32 j=0; j<len; j++) {
      if (consensus[i][j] != '-') {
        ugconsensus[i][k] = consensus[i][j];
        imap[i][k] = j;
        k++;
      }
    }

    uglen[i] = k;
  }
}

/* Count gaps in the short and long consensus sequences */
static
void
CountGaps(char **consensus, int32 len, int32 *gapcount) {

  for (int32 i=0; i<2; i++)
    {
      int32 last_base = len-1;
      while ((last_base > 0) && (consensus[i][last_base] == '-'))
        last_base--;

      gapcount[i] = 0;
      int32 first_base = -1;
      for (int32 j=0; j<last_base + 1; j++)
        {
          if (consensus[i][j] != '-')
            first_base = j;

          if ((first_base >= 0) && (consensus[i][j] == '-'))
            gapcount[i]++;
        }
    }
}

/*
  Find an adjusted left boundary for long and short ungapped sequences (that is,
  the leftmost position starting from which the sequences will match):
  6.1 select a size of a "probe" k-mer, say, 3
  6.2 scan the short ungapped consensus from the left to the right
  6.3 for each position in the short ungapped consensus, get k-mer starting at
  this position
  6.4 find the leftmost position of this k-mer in the long ungapped sequence
  6.5 set adjusted left boundary of short ungapped consensus to the position of
  k-mer with leftmost occurrence in the long ungapped sequence
*/
static
void
FindAdjustedLeftBounds(int32 *adjleft, char **ugconsensus, int32 *uglen,
                       int32 short_allele, int32 long_allele) {
  int32   s, l;
  char *ps, *pl;

  adjleft[short_allele] = uglen[short_allele]-1;
  adjleft[long_allele]  = uglen[long_allele]-1;
  for (s=0; s<uglen[short_allele] - MSTRING_SIZE; s++)
    {
      ps = ugconsensus[short_allele] + s;
      for (l=0; l<uglen[long_allele] - MSTRING_SIZE; l++)
        {
          pl = ugconsensus[long_allele] + l;
          if (strncmp(pl, ps, MSTRING_SIZE) == 0)
            {
              if (adjleft[0] + adjleft[1] > s+l)
                {
                  adjleft[long_allele]  = l;
                  adjleft[short_allele] = s;
                }
            }
        }
    }
  if ((adjleft[long_allele]  == uglen[long_allele]-1) &&
      (adjleft[short_allele] == uglen[short_allele]-1))
    {
      adjleft[short_allele] = 0;
      adjleft[long_allele]  = 0;
    }
}

static
void
FindAdjustedRightBounds(int32 *adjright,  char **ugconsensus, int32 *uglen,
                        int32 short_allele, int32 long_allele) {
  int32   s, l;
  char *ps, *pl;

  adjright[short_allele] = uglen[short_allele]-1;
  adjright[long_allele]  = uglen[ long_allele]-1;
  for (s=uglen[short_allele] - MSTRING_SIZE-1; s>= 0; s--)
    {
      ps = ugconsensus[short_allele] + s;
      for (l=uglen[long_allele] - MSTRING_SIZE-1; l>=0; l--)
        {
          pl = ugconsensus[long_allele] + l;
          if (strncmp(pl, ps, MSTRING_SIZE) == 0)
            {
              if (adjright[0] + adjright[1] >
                  uglen[short_allele]-1 -(s+MSTRING_SIZE)
                  +uglen[long_allele] -1 -(l+MSTRING_SIZE))
                {
                  adjright[long_allele] =
                    uglen[long_allele]-1-(l+MSTRING_SIZE);
                  adjright[short_allele] =
                    uglen[short_allele]-1 -(s+MSTRING_SIZE);
                }
            }
        }
    }
  if ((adjright[long_allele]  == uglen[long_allele]-1) &&
      (adjright[short_allele] == uglen[short_allele]-1))
    {
      adjright[short_allele] = 0;
      adjright[long_allele]  = 0;
    }
}

static
void
GetLeftScore(char **ugconsensus, int32 *uglen, int32 **imap, int32 *adjleft,
             int32 short_allele, int32 long_allele, int32 *maxscore, int32 *maxpos) {
  int32 i, score = 0;

  *maxscore = 0;
  *maxpos   = adjleft[short_allele];
  i = 0;
  while ((i < uglen[short_allele] - adjleft[short_allele]) &&
         (i < uglen[ long_allele] - adjleft[ long_allele]))
    {
      int32 lpos = i + adjleft[long_allele];
      int32 spos = i + adjleft[short_allele];
      if (ugconsensus[short_allele][spos] == ugconsensus[long_allele][lpos])
        score++;
      else
        score--;
      if (*maxscore < score)
        {
          *maxscore = score;
          *maxpos   = spos;
        }
      i++;
    }
  /* Position in ungapped consensus */
  *maxpos = imap[short_allele][*maxpos];
}

static
void
GetRightScore(char **ugconsensus, int32 *uglen, int32 **imap, int32 *adjright,
              int32 short_allele, int32 long_allele, int32 *maxscore, int32 *maxpos) {
  int32 i, j, score = 0;

  *maxscore = 0;
  *maxpos   = uglen[short_allele]-1-adjright[short_allele];
  i = uglen[long_allele]-1;
  j = uglen[short_allele]-1;
  while ((j >= adjright[short_allele]) &&
         (i >= adjright[ long_allele]))
    {
      int32 lpos = i - adjright[long_allele];
      int32 spos = j - adjright[short_allele];
      if (ugconsensus[short_allele][spos] == ugconsensus[long_allele][lpos])
        score++;
      else
        score--;
      if (*maxscore < score)
        {
          *maxscore = score;
          *maxpos = spos;
        }
      i--;
      j--;
    }
  /* Position in ungapped consensus */
#if 0
  fprintf(stderr, "long_allele=%d  maxpos =%d \n", long_allele, *maxpos);
#endif
  *maxpos = imap[short_allele][*maxpos];
}

static
void
AdjustShiftingInterfaces(int32 *lpos, int32 *rpos, int32 lscore, int32 rscore,
                         int32 *adjleft, int32 *adjright, int32 long_allele, int32 short_allele) {
#if 0
  fprintf(stderr, "\nlpos=%d rpos=%d lscore=%d rscore=%d \n", *lpos, *rpos, lscore, rscore);
  fprintf(stderr, "adjleft = %d %d  adjright= %d %d \n", adjleft[0], adjleft[1],
          adjright[0], adjright[1]);
#endif
  if (adjleft[long_allele] > 5)
    {
      *lpos = -1;
      lscore = -1;
    }
  if (adjright[long_allele] > MAX_SIZE_OF_ADJUSTED_REGION)
    {
      *rpos = -1;
      rscore = -1;
    }

  /* Set teh posiytion of shifting interface */
  if (*lpos <= *rpos)
    {
    }
  else /* lpos >  rpos */
    {
      if ((lscore > 0) && (rscore > 0))
        {
          //         if (adjleft[short_allele] <= adjright[short_allele])
          if (lscore > rscore)
            *rpos = *lpos;
          else
            *lpos = *rpos;
        }
      else if ((lscore > 0) && (rscore <= 0))
        *rpos = -1;
      else
        *lpos = -1;
    }
}




static
void
GetTemplateForAbacus(char   *&tmpl,
                     char    *consensus[2],
                     int32   len,
                     char   *ugconsensus[2],
                     int32  *uglen,
                     int32   lpos,
                     int32   rpos,
                     int32  *imap[2],
                     int32  *adjleft,
                     int32  *adjright,
                     int32   short_allele,
                     int32   long_allele) {
  int32 i, j;

  tmpl = new char [len];

  for (i=0; i<len; i++)
    tmpl[i] = consensus[long_allele][i];

  /* Set Ns in the left part of the tmpl */
  i = 0;
  while ((imap[long_allele][i] <= lpos) &&
         (i < uglen[short_allele] - adjleft[short_allele]) &&
         (i < uglen[ long_allele] - adjleft[ long_allele]))
    {
      int32 lpos = i + adjleft[long_allele];
      int32 spos = i + adjleft[short_allele];
      if ((ugconsensus[short_allele][spos] != ugconsensus[long_allele][lpos]) &&
          (tmpl[imap[long_allele][lpos]] != '-'))
        tmpl[imap[long_allele][lpos]] = 'n';
      i++;
    }

  /* tmpl bases before adjusted left boundary should be 'N's */
  if (adjleft[long_allele] > 0 && lpos > 0)
    {
      i = imap[long_allele][adjleft[long_allele]-1];
      j = 0;
      while ((j < adjleft[short_allele]) && (i >= 0))
        {
          if (consensus[long_allele][i] != '-')
            {
              tmpl[i] = 'n';
              j++;
              i--;
            }
          else
            i--;
        }
    }

  /* Set Ns in the right part of the tmpl */
  i = uglen[long_allele]-1-adjright[long_allele];
  j = uglen[short_allele]-1-adjright[short_allele];
  while ((i >= adjleft[long_allele]) &&
         (j >= adjleft[short_allele]) &&
         (imap[long_allele][i] > rpos))
    {
      if ((ugconsensus[short_allele][j] != ugconsensus[long_allele][i]) &&
          (tmpl[imap[long_allele][i]] != '-'))
        {
          tmpl[imap[long_allele][i]] =  'n';
        }
      i--;
      j--;
    }

  /* tmpl bases after adjusted right boundary should be 'N's */
  if ((adjright[long_allele] > 0) && (rpos > 0))
    {
      for (i = uglen[long_allele]-adjright[long_allele];
           i < uglen[long_allele];
           i++)
        {
          j = imap[long_allele][i];
          if (consensus[long_allele][j] != '-')
            {
              tmpl[i] = 'n';
            }
        }
    }
}




int32
abMultiAlign::refineWindow(abAbacus     *abacus,
                           abColumn     *start_column,
                           abColID       stab_bgn) {

  int32  orig_columns  = 0;
  int32  left_columns  = 0;
  int32  right_columns = 0;
  int32  best_columns  = 0;

  // Mismatch, gap and total scores:
  int32   orig_mm_score     = 0;
  int32   left_mm_score     = 0;
  int32   right_mm_score    = 0;
  int32   best_mm_score     = 0;
  int32   orig_gap_score    = 0;
  int32   left_gap_score    = 0;
  int32   right_gap_score   = 0;
  int32   best_gap_score    = 0;
  int32   orig_total_score  = 0;
  int32   left_total_score  = 0;
  int32   right_total_score = 0;
  int32   best_total_score  = 0;
  int32   max_element       = 0;
  int32   score_reduction   = 0;

  abBaseCount    abacus_count;
  abAbacusWork  *orig_abacus   = new abAbacusWork(abacus, ident(), start_column->ident(), stab_bgn);
  abAbacusWork  *left_abacus   = NULL;
  abAbacusWork  *right_abacus  = NULL;
  abAbacusWork  *best_abacus   = NULL;

  abVarRegion      vreg;

  vreg.nr = orig_abacus->rows;

  //ShowAbacus(orig_abacus);

  // Process reads of the original abacus
  //AllocateDistMatrix(&vreg, 0);
  //AllocateMemoryForReads(&vreg.reads, orig_abacus->rows, orig_abacus->columns, QV_FOR_MULTI_GAP);

  vreg.allocateDistanceMatrix(0);

  vreg.beg = 0;
  vreg.end = orig_abacus->columns-1;

  orig_mm_score = orig_abacus->scoreAbacus(orig_columns);

  //GetReadsForAbacus(vreg.reads, orig_abacus);
  {
#if 0
    fprintf(stderr, "rows=%lu shift=%c window_width=%lu \nReads= \n",
            orig_abacus->rows, (char)orig_abacus->shift, orig_abacus->window_width);
#endif

    int32 shift = 0;

    if (orig_abacus->shift == UNSHIFTED)
      shift =     orig_abacus->columns;

    if (orig_abacus->shift == RIGHT_SHIFT)
      shift = 2 * orig_abacus->columns;

    for (uint32 i=0; i<orig_abacus->rows; i++) {
      vreg.reads[i].id = i;
      vreg.reads[i].allele_id = -1;
      vreg.reads[i].ave_qv = 20.;        // qvs are hardly available for abacus
      vreg.reads[i].uglen = 0;

      for (uint32 j=0; j<orig_abacus->columns; j++) {
        char  base = orig_abacus->getChar(i,j);

        if (is_good_base(base))
          vreg.reads[i].bases[j] = base;

        if (base != '-')
          vreg.reads[i].uglen++;
      }
    }
  }



#ifdef DEBUG_VAR_RECORDS
  OutputReads(stderr, vreg.reads, vreg.nr, orig_abacus->columns);
#endif

  //PopulateDistMatrix(vreg.reads, orig_abacus->columns, &vreg);
  vreg.populateDistanceMatrix();

#ifdef DEBUG_VAR_RECORDS
  OutputDistMatrix(stderr, &vreg);
#endif

  //AllocateMemoryForAlleles(&vreg.alleles, vreg.nr, &vreg.na);

  vreg.clusterReads();

  vreg.sortAllelesByLength();

  orig_abacus->refineOrigAbacus(vreg);
  //ShowAbacus(orig_abacus);
  orig_mm_score = orig_abacus->scoreAbacus(orig_columns);

#if 0
  fprintf(stderr, "vreg.alleles= ");
  for (i=0; i<vreg.nr; i++)
    fprintf(stderr, "%d", vreg.alleles[i]);

#endif

#ifdef DEBUG_ABACUS
  fprintf(stderr, "\n\nOrigCalls=\n");
  ShowCalls(orig_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(orig_abacus);
  fprintf(stderr, "\n");
#endif
  //ShowAbacus(orig_abacus);
  left_abacus = orig_abacus->clone();
  left_mm_score = left_abacus->leftShift(vreg, left_columns);
#ifdef DEBUG_ABACUS
  fprintf(stderr, "\n\nLeftShiftCalls=\n");
  ShowCalls(left_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(left_abacus);
  fprintf(stderr, "\n");
#endif
  right_abacus = orig_abacus->clone();
  right_mm_score = right_abacus->rightShift(vreg, right_columns);
#ifdef DEBUG_ABACUS
  fprintf(stderr, "\n\nRightShiftCalls=\n");
  ShowCalls(right_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(right_abacus);
  fprintf(stderr, "\n");
#endif
  //fprintf(stderr,"Abacus Report:\norig_mm_score: %d left_mm_score: %d right_mm_score: %d\n",
  //             orig_mm_score,left_mm_score,right_mm_score);
  //ShowAbacus(left_abacus);
  //ShowAbacus(right_abacus);
  // determine best score and apply abacus to real columns
  orig_gap_score  = orig_abacus->affineScoreAbacus();
  left_gap_score  = left_abacus->affineScoreAbacus();
  right_gap_score = right_abacus->affineScoreAbacus();
  best_abacus     = orig_abacus;
  best_columns    = orig_columns;
  best_gap_score  = orig_gap_score;
  best_mm_score   = orig_mm_score;
  orig_total_score  = orig_mm_score  + orig_columns  + orig_gap_score;
  left_total_score  = left_mm_score  + left_columns  + left_gap_score;
  right_total_score = right_mm_score + right_columns + right_gap_score;
  best_total_score  = orig_total_score;

#ifdef DEBUG_ABACUS
  fprintf(stderr, "In refineWindow: beg= %lu end= %d\n",
          start_column->ident(), stab_bgn);
  fprintf(stderr, "    abacus->columns= %d, abacus->rows= %d\n",
          orig_abacus->columns, orig_abacus->rows);
  fprintf(stderr, "    w_width left= %d orig= %d right= %d\n",
          left_abacus->window_width, orig_abacus->window_width,
          right_abacus->window_width);
  fprintf(stderr, "    mm_score left= %d orig= %d right= %d\n",
          left_mm_score, orig_mm_score, right_mm_score);
  fprintf(stderr, "     columns left= %d orig= %d right= %d\n",
          left_columns, orig_columns, right_columns);
  fprintf(stderr, "   gap_score left= %d orig= %d right= %d\n",
          left_gap_score, orig_gap_score, right_gap_score);
  fprintf(stderr, " total_score left= %d orig= %d right= %d\n",
          left_total_score, orig_total_score, right_total_score);
#endif

  // Use the total score to refine the abacus
  if (left_total_score < orig_total_score || right_total_score < orig_total_score ) {
    if (left_total_score <= right_total_score ) {
      score_reduction += orig_total_score - left_total_score;
      //fprintf(stderr,"\nTry to apply LEFT abacus:\n");
      //ShowAbacus(left_abacus);
      left_abacus->GetBaseCount(abacus_count);
#if 0
      fprintf(stderr, " Applying left abacus\n");
#endif
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
#if 0
      fprintf(stderr, " Applying right abacus\n");
#endif
      best_abacus      = right_abacus;
      best_mm_score    = right_mm_score;
      best_columns     = right_columns;
      best_gap_score   = right_gap_score;
      best_total_score = right_total_score;
    }
  }

#if 0
  fprintf(stderr, "Best Abacus Before MixedShift=\n");
  ShowAbacus(best_abacus);
#endif
  {
    int32 i, j;
    //char    **consensus=NULL;
    //char    **ugconsensus=NULL;
    int32     adjleft[2]={-1,-1};
    int32     adjright[2]={-1,-1};
    int32     gapcount[2];
    int32     short_allele=-1;
    int32     long_allele=-1;
    int32     lscore=0;
    int32     rscore=0;
    int32     lpos=-1;
    int32     rpos=-1;
    int32     mixed_columns=0;
    int32     mixed_mm_score=0;
    int32     mixed_gap_score=0;


    char   *consensus[2] = { NULL, NULL };

    getConsensusForAbacus(vreg, vreg.reads, best_abacus, consensus);

#if 0
    fprintf(stderr, "\nconsensus[0]=\n");
    for (i=0; i<3*best_abacus->window_width; i++)
      fprintf(stderr, "%c", consensus[0][i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "consensus[1]=\n");
    for (i=0; i<3*best_abacus->window_width; i++)
      fprintf(stderr, "%c", consensus[1][i]);
    fprintf(stderr, "\n\n");
#endif

    CountGaps(consensus, 3*best_abacus->window_width, gapcount);
    short_allele = (gapcount[0] >= gapcount[1]) ? 0 : 1;
    long_allele  = (gapcount[0] <  gapcount[1]) ? 0 : 1;

#if 0
    fprintf(stderr, "gapcounts[0, 1] = %d %d\n", gapcount[0], gapcount[1]);
#endif

    if (gapcount[short_allele] == 0) {
#if 0
      fprintf(stderr, "No MixedShift will be performed: gapcount[short_allele] = %d\n", gapcount[short_allele]);
#endif
      best_abacus->applyAbacus(abacus);

      delete    orig_abacus;
      delete    left_abacus;
      delete    right_abacus;

      delete [] consensus[0];
      delete [] consensus[1];

      return(score_reduction);
    }

    /* Now try the mixed consensus */

    char   *ugconsensus[2] = { NULL, NULL };
    int32  *imap[2]        = { NULL, NULL };
    int32   uglen[2]       = {    0,    0 };

    mapConsensus(imap, consensus, ugconsensus, 3 * best_abacus->window_width, uglen);

    if ((uglen[0] < MSTRING_SIZE) || (uglen[1] < MSTRING_SIZE)) {

#if 0
      fprintf(stderr, "No MixedShift will be performed: uglen = %d %d\n", uglen[0], uglen[1]);
#endif

      best_abacus->applyAbacus(abacus);

      delete    orig_abacus;
      delete    left_abacus;
      delete    right_abacus;

      delete [] consensus[0];
      delete [] consensus[1];

      delete [] ugconsensus[0];
      delete [] ugconsensus[1];

      delete [] imap[0];
      delete [] imap[1];

      return(score_reduction);
    }

#if 0
    fprintf(stderr, "\nuglen[0]=%d ugconsensus[0] =\n", uglen[0]);
    for (i=0; i<uglen[0]; i++)
      fprintf(stderr, "%c", ugconsensus[0][i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "uglen[1]=%d ugconsensus[1] =\n", uglen[1]);
    for (i=0; i<uglen[1]; i++)
      fprintf(stderr, "%c", ugconsensus[1][i]);
    fprintf(stderr, "\n\n");
#endif

    FindAdjustedLeftBounds(adjleft, ugconsensus, uglen, short_allele, long_allele);
    FindAdjustedRightBounds(adjright, ugconsensus, uglen, short_allele, long_allele);

#if 0
    fprintf(stderr, "Adjusted left bounds 0, 1= %d %d \n", adjleft[0], adjleft[1]);
    fprintf(stderr, "Adjusted right bounds 0, 1= %d %d \n", adjright[0], adjright[1]);
#endif

    GetLeftScore(ugconsensus, uglen, imap, adjleft, short_allele, long_allele, &lscore, &lpos);
    GetRightScore(ugconsensus, uglen, imap, adjright, short_allele, long_allele, &rscore, &rpos);
    AdjustShiftingInterfaces(&lpos, &rpos, lscore, rscore, adjleft, adjright, long_allele, short_allele);

    char *tmpl = NULL;
    GetTemplateForAbacus(tmpl, consensus, 3*best_abacus->window_width, ugconsensus, uglen, lpos, rpos, imap, adjleft, adjright, short_allele, long_allele);

    abAbacusWork *mixed_abacus = orig_abacus->clone();

#if 0
    {
      fprintf(stderr, "Template = \n");
      for (i=0; i<3*best_abacus->window_width; i++)
        fprintf(stderr, "%c", tmpl[i]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "Start calling MixedShift\n\n");
    fprintf(stderr, "   Final lpos=%d rpos=%d window_width=%d long_allele=%d short_allele=%d\n",
            lpos, rpos, best_abacus->window_width, long_allele, short_allele);
#endif

    mixed_mm_score = mixed_abacus->mixedShift(mixed_columns, vreg, lpos, rpos, tmpl, long_allele, short_allele);


#if 0
    fprintf(stderr, "Mixed abacus=\n");
    ShowAbacus(mixed_abacus);
    fprintf(stderr, "End calling MixedShift\n\n");
    fprintf(stderr, "mixed_mm_score=%d bast_score=%d\n", mixed_mm_score, best_mm_score);
    fprintf(stderr, "mixed_columns=%d best_columns=%d\n", mixed_columns, best_columns);
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d\n", mixed_gap_score, best_gap_score);
#endif
    mixed_gap_score  = mixed_abacus->affineScoreAbacus();
#if 0
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d mixed_columns=%d best_columns=%d mixed_mm_score=%d best_mm_score=%d\n", mixed_gap_score, best_gap_score, mixed_columns, best_columns, mixed_mm_score, best_mm_score);
#endif
    if ((mixed_gap_score <  best_gap_score) ||
        ((mixed_gap_score == best_gap_score) && (mixed_columns < best_columns))
        ||
        ((mixed_gap_score == best_gap_score) && (mixed_columns == best_columns) &&
         (mixed_mm_score < best_mm_score))) {
      best_abacus    = mixed_abacus;
      best_mm_score  = mixed_mm_score;
      best_columns   = mixed_columns;
      best_gap_score = mixed_gap_score;
    }
#if 0
    ShowCalls(best_abacus);
    fprintf(stderr, "Best Abacus after MixedShift =\n");
    ShowAbacus(best_abacus);
#endif

    //      OutputDistMatrix(stderr, &vreg);

#if 0
    fprintf(stderr, "Consensus0 =\n");
    for (int32 j=0; j<3*best_abacus->window_width; j++)
      fprintf(stderr, "%c", consensus[0][j]);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "Consensus1 =\n");
    for (int32 j=0; j<3*best_abacus->window_width; j++)
      fprintf(stderr, "%c", consensus[1][j]);
    fprintf(stderr, "\n\n");
#endif

    /* Otherwise, try to do a more sophisticated shift:
     * - only shifting reads of the shortest allele
     * - only within a subregion of abacus window where the alleles match
     */
#if 0
    fprintf(stderr, "Applying the Best abacus\n");
#endif
    best_abacus->applyAbacus(abacus);

    //      fprintf(stderr, "vreg.nr = %d\n", vreg.nr);

    delete orig_abacus;
    delete left_abacus;
    delete right_abacus;
    delete mixed_abacus;

    delete [] consensus[0];
    delete [] consensus[1];

    delete [] ugconsensus[0];
    delete [] ugconsensus[1];

    delete [] imap[0];
    delete [] imap[1];

    delete [] tmpl;
  }

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
                     uint32               from,
                     uint32               to) {

  // from and to are in ma's column coordinates
#if 0
  int32 sid, eid, stab_bgn;
  int32 ma_length = GetMANodeLength(ident());
  int32 score_reduction=0;
  int32 orig_length = ma_length;
  int32 refined_length = orig_length;
  abColumn *start_column;
  int32 i;
#endif

  if (length() < to)
    to = length();

  abColID sid = columnList[from];   // id of the starting column
  abColID eid = columnList[to];     // id of the ending column

  abColumn *start_column = abacus->getColumn(sid);
  abColID   stab_bgn;

  uint32    score_reduction = 0;

  while (start_column->ident() != eid) {
    int32 window_width = IdentifyWindow(abacus, start_column, stab_bgn, level);

    // start_column stands as the candidate for first column in window
    // look for window start and stop

    if (window_width > 0) {

      //  If the first column, insert a gap column for maneuvering room
      if (start_column->prevID().isValid() == false) {
        abBeadID   firstbeadID = abacus->getBead( start_column->callID() )->downID();
        abBeadID   newbeadID   = abacus->appendGapBead( abacus->getBead(firstbeadID)->ident() );

        fprintf(stderr, "Adding gapbead "F_U32" after "F_U32" to add abacus room for abutting left of multialignment\n",
                newbeadID.get(), firstbeadID.get());

        abacus->appendColumn(abacus->getBead(firstbeadID)->colIdx(), newbeadID);
      }

      //  if the window is too big, there's likely a polymorphism that won't respond well to abacus,
      //  so skip it.
      //
      //  BPW saw crashes with large window_width's (1333, 3252, 1858, 675, 855, 1563, 601, 1102).
      //  The longest window_width that worked was 573.  Previous versions used 100 here.  Not sure
      //  what it should be.
      //
      if (window_width < MAX_WINDOW_FOR_ABACUS_REFINE)
        score_reduction += refineWindow(abacus, start_column, stab_bgn);
      
      start_column = abacus->getColumn(stab_bgn);
    }

    start_column = abacus->getColumn(stab_bgn);
  }
 
  //  WITH quality=1 make_v_list=1, all the rest defaults
  abacus->refreshMultiAlign(ident(), 1, 1, 11, 1, NULL, NULL, 1, false);

  return(score_reduction);
}
