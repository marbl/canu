
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

static char *rcsid = "$Id: AbacusRefine.c,v 1.12 2011-12-09 02:59:56 brianwalenz Exp $";

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <map>

using namespace std;

#include "AS_global.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


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


static
char *
GetAbacus(AbacusDataStructure *a, int32 i, int32 j) {
  return (a->beads+i*(a->columns+2)+j+1);
}


static
void
SetAbacus(AbacusDataStructure *a, int32 i, int32 j, char c) {
  int32 offset = i * (a->columns+2) + j + 1;

  if ((i == -1) || (i > a->rows-1)) {
    fprintf(stderr, "i=%d j=%d a->rows=%d\n", i, j, a->rows);
    fprintf(stderr, "SetAbacus attempt to write beyond row range");
    assert(0);
  }

  if ((j == -1) || (j > a->columns-1)) {
    fprintf(stderr, "i=%d j=%d a->columns=%d\n", i, j, a->columns);
    fprintf(stderr, "SetAbacus attempt to write beyond column range");
    assert(0);
  }

  a->beads[offset] = c;
}


static
void
ResetCalls(AbacusDataStructure *a) {
  for (int32 j=0;j<a->columns;j++)
    a->calls[j] = 'n';
}


static
void
ResetIndex(VA_TYPE(int32) * indices, int32 n) {
  int32 value=0;
  for(int32 i=0;i<n;i++)
    SetVA_int32(indices,i,&value);
}


static
AbacusDataStructure *
CreateAbacus(int32 mid, int32 from, int32 end) {
  // from,and end are ids of the first and last columns in the columnStore
  AbacusDataStructure             *abacus;
  beadIdx              bid;
  int32               columns=1, rows=0, i, j, orig_columns, set_column;
  Column             *column,*last;
  ColumnBeadIterator  bi;
  Bead               *bead;
  MANode             *ma;
#define  MAX_MID_COLUMN_NUM 100
  static int32        mid_column_points[MAX_MID_COLUMN_NUM] = { 75, 150};
  Column             *mid_column[MAX_MID_COLUMN_NUM] = { NULL, NULL };
  int32               next_mid_column=0;
  int32               max_mid_columns = 0;

  //  Macaque, using overlap based trimming, needed mid_column points
  //  at rather small intervals to pass.  Even without OBT, macaque
  //  needed another point at 63 to pass.
  //
  //  This change was tested on macaque, and did not change the
  //  results (except for allowing one partition to finish....).  You
  //  can revert to the original behavior by undef'ing the following.
  //
#define EXTRA_MID_COLUMNS 1
#ifdef EXTRA_MID_COLUMNS
  if(mid_column_points[0]==75){
    for (i=0; i<MAX_MID_COLUMN_NUM; i++){
      mid_column_points[i] = i * ((AS_READ_MIN_LEN/2)-1) + 30;
      mid_column[i]=NULL;
    }
  }
  max_mid_columns=MAX_MID_COLUMN_NUM;

#else
  max_mid_columns=2;
#endif

  ma = GetMANode(manodeStore, mid);
  column = GetColumn(columnStore, from);

  assert(ma != NULL);
  assert(column != NULL);
  assert(abacus_indices != NULL);

  ResetIndex(abacus_indices,GetNumFragments(fragmentStore));

#if 0
  fprintf(stderr, "column_ids= %lu", from);
#endif
  // first, just determine requires number of rows and columns for Abacus
  while( column->next != end  && column->next != -1) { // this test looks subtly wrong: most of the loop could be done
    // as long as column!=NULL?  -- ALH
    columns++;

#if 0
    fprintf(stderr, ",%lu", column->next);
#endif

    if( next_mid_column<MAX_MID_COLUMN_NUM &&
        columns == mid_column_points[next_mid_column]){
      mid_column[next_mid_column] = GetColumn(columnStore,column->lid);
      next_mid_column++;
    }
    column = GetColumn(columnStore,column->next);
    // GD: this is where base calling code should be called
  }

#if 0
  fprintf(stderr, "\n");
#endif

  if(columns>MAX_MID_COLUMN_NUM*((AS_READ_MIN_LEN/2)-1)){
    fprintf(stderr,"WARNING: CreateAbacus called with such a large window that small fragments might slip through the cracks...\n");
  }

  max_mid_columns=next_mid_column;

  orig_columns = columns;
  last = column;
  column = GetColumn(columnStore, from);

  CreateColumnBeadIterator(column->lid, &bi);

  while ( (bid = NextColumnBead(&bi)).isValid() ) {
    bead = GetBead(beadStore,bid);
    rows++;
    SetVA_int32(abacus_indices,bead->frag_index,&rows);
  }

  CreateColumnBeadIterator(last->lid, &bi);

  while ( (bid = NextColumnBead(&bi)).isValid() ) {
    bead = GetBead(beadStore,bid);
    if ( *Getint32(abacus_indices,bead->frag_index) == 0 ) {
      rows++;
      SetVA_int32(abacus_indices,bead->frag_index,&rows);
    }
  }

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
  for (i=0; i<max_mid_columns; i++) {
    if (mid_column[i] != NULL) {
      CreateColumnBeadIterator(mid_column[i]->lid, &bi);

      while ((bid = NextColumnBead(&bi)).isValid()) {
        bead = GetBead(beadStore,bid);
        if ( *Getint32(abacus_indices,bead->frag_index) == 0 ) {
          rows++;
          SetVA_int32(abacus_indices,bead->frag_index,&rows);
        }
      }
    }
  }

  abacus = (AbacusDataStructure *) safe_malloc(sizeof(AbacusDataStructure));
  abacus->start_column = from;
  abacus->end_column = last->lid;
  abacus->rows = rows;
  abacus->window_width = orig_columns;
  abacus->columns = 3*orig_columns;
  abacus->shift = UNSHIFTED;
  abacus->beads = (char *) safe_calloc(rows*(abacus->columns+2),sizeof(char));
  abacus->calls = (char *) safe_calloc((abacus->columns),sizeof(char));
  // two extra gap columns, plus "null" borders

  // now, fill the center third of abacus with chars from the columns

  for (i=0;i<rows*(abacus->columns+2);i++) {
    abacus->beads[i] = 'n'; // initialize to "null" code
  }
  columns = 0;
  while( column->lid != end  && column->lid != -1) {
    CreateColumnBeadIterator(column->lid, &bi);

    set_column = columns+orig_columns;
    while ( (bid = NextColumnBead(&bi)).isValid() ) {
      bead = GetBead(beadStore,bid);
      SetAbacus(abacus, *Getint32(abacus_indices,bead->frag_index)-1,
                set_column, *Getchar(sequenceStore,bead->soffset));
    }
    columns++;
    column = GetColumn(columnStore,column->next);
  }

  for (i=0;i<rows;i++) {
    set_column = orig_columns;
    for (j=0;j<set_column;j++) {
      SetAbacus(abacus,i,j,'-');
    }
    set_column = 2*orig_columns-1;
    for (j=set_column+1;j<abacus->columns;j++) {
      SetAbacus(abacus,i,j,'-');
    }
  }
  ResetCalls(abacus);
  return abacus;
}

static
void
DeleteAbacus(AbacusDataStructure *abacus) {
  safe_free(abacus->beads);
  safe_free(abacus->calls);
  safe_free(abacus);
}

static
AbacusDataStructure *
CloneAbacus(AbacusDataStructure *abacus) {
  AbacusDataStructure *clone;
  int32 rows=abacus->rows;
  int32 columns=abacus->columns;
  clone = (AbacusDataStructure *) safe_malloc(sizeof(AbacusDataStructure));
  clone->beads = (char *) safe_calloc(rows*(columns+2),sizeof(char)); //
  clone->calls = (char *) safe_calloc((columns),sizeof(char));
  clone->rows = rows;
  clone->window_width = abacus->window_width;
  clone->columns = columns;
  clone->start_column = abacus->start_column;
  clone->end_column = abacus->end_column;
  clone->shift = abacus->shift;
  memcpy(clone->beads,abacus->beads,rows*(columns+2)*sizeof(char));
  memcpy(clone->calls,abacus->calls,columns*sizeof(char));
  return clone;
}


static
void
ShowAbacus(AbacusDataStructure *abacus) {
  char form[10];
  sprintf(form,"%%%d.%ds\n",abacus->columns,abacus->columns);
  fprintf(stderr,"\nstart column: %d\n",abacus->start_column);
  for (int32 i=0;i<abacus->rows;i++)
    fprintf(stderr,form,GetAbacus(abacus,i,0));
  fprintf(stderr,"\n");
  fprintf(stderr,form,abacus->calls);
}

  // cols is the number of "good" (non-null) columns found
  // GD: This function counts the total number of bases which
  //   - are different from column's "consensus" call and
  //   - are not 'n'
  //
static
int32
ScoreAbacus(AbacusDataStructure *abacus, int32 *cols) {
  BaseCount *counts;
  int32 score=0;
  char b;

  counts = (BaseCount *) safe_calloc(abacus->columns,sizeof(BaseCount));
  memset(counts,'\0',abacus->columns*sizeof(BaseCount));
  *cols=0;

  for (int32 i=0;i<abacus->rows;i++) {
    for (int32 j=0;j<abacus->columns;j++) {
      b = *GetAbacus(abacus,i,j);
      if ( b == '-' ) {
        if ( j>0 && j < abacus->columns-1) {
          if ((*GetAbacus(abacus,i,j-1) == 'n')  ||
              (*GetAbacus(abacus,i,j+1) == 'n') )
            {
              b = 'n';
            }
        }
      }
      IncBaseCount(&counts[j],b);
    }
  }
  // now, for each column, generate the majority call
  for (int32 j=0;j<abacus->columns;j++) {
    if ( GetBaseCount(&counts[j],'-') + GetBaseCount(&counts[j],'n') == counts[j].depth ) {
      // null (all-gap) column. Flag with an 'n' basecall
      abacus->calls[j] = 'n';
    }
    else {
      *cols=*cols+1;
      abacus->calls[j] = GetMaxBaseCount(&counts[j],0);
      // and then tally edit score
      score += counts[j].depth - counts[j].count[RINDEX[abacus->calls[j]]] - counts[j].count[CNS_NALPHABET-1]; // don't count 'n's
    }
  }

  safe_free(counts);
  return score;
}

static
int32
AffineScoreAbacus(AbacusDataStructure *abacus) {
  // This simply counts the number of opened gaps, to be used in tie breaker
  //   of edit scores.
  int32 score=0;
  char b;
  int32 start_column, end_column;

  if (abacus->shift == LEFT_SHIFT)
    {
      start_column = 0;
      end_column   = abacus->columns/3;
    }
  else if (abacus->shift == RIGHT_SHIFT)
    {
      start_column = 2*abacus->columns/3;
      end_column   =   abacus->columns;
    }
  else //  abacus->shift == UNSHIFTED
    {
      start_column =   abacus->columns/3;
      end_column   = 2*abacus->columns/3;
    }

  for (int32 i=0;i<abacus->rows;i++)
    {
      int32 in_gap=0;
      for (int32 j=start_column;j<end_column;j++)
        {
          b = *GetAbacus(abacus,i,j);
          //      if ( abacus->calls[j] != 'n')
          //      commented out in order to make gap_score
          //      of the orig_abacus non-zero - GD
          {// don't look at null columns
            if ( b != '-' )
              {
                in_gap=0;
              }
            else
              {
                // Size of a gap does not matter, their number in a row does - GD
                if ( ! in_gap )
                  {
                    in_gap = 1;
                    score++;
                  }
              }
          }
        }
    }
  return score;
}

static
int
MergeAbacus(AbacusDataStructure *abacus, int32 merge_dir) {
  // sweep through abacus from left to right
  // testing for Level 1 (neighbor) merge compatibility of each column
  // with right neighbor and merge if compatible
  //
  //  GD: this code will merge practically any
  int32  mergeok, next_column_good, curr_column_good;
  char   prev, curr, next;
  int32  last_non_null = abacus->columns-1;
  int32 first_non_null = 0;
  int32 columns_merged = 0;

  // determine the rightmost and leftmost columns
  // not totally composed of gaps
  for (int32 j=abacus->columns-1;j>0;j--)
    {
      int32 null_column = 1;
      for (int32 i=0; i<abacus->rows; i++) {
        curr = *GetAbacus(abacus,i,j);
        if (curr != '-') null_column = 0;
      }
      if (!null_column)
        break;
      last_non_null = j;
    }
  for (int32 j=0; j<abacus->columns;j++)
    {
      int32 null_column = 1;
      for (int32 i=0; i<abacus->rows; i++) {
        curr = *GetAbacus(abacus,i,j);
        if (curr != '-')
          null_column = 0;
      }
      if (!null_column)
        break;
      first_non_null = j;
    }
#ifdef DEBUG_ABACUS
  fprintf(stderr, "abacus->columns=%d first_non_null = %d last_non_null= %d\n",
          abacus->columns, first_non_null, last_non_null);
#endif
  if (merge_dir < 0)
    {
      for (int32 j=0;j<last_non_null;j++)
        {
          int32 num_gaps=0, num_ns=0;
          mergeok = 1;
          next_column_good = -1;
          for (int32 i=0;i<abacus->rows;i++)
            {
              curr = *GetAbacus(abacus,i,j);
              next = *GetAbacus(abacus,i,j+1);
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
              for (int32 i=0;i<abacus->rows;i++) {
                curr = *GetAbacus(abacus,i,j  );
                next = *GetAbacus(abacus,i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (next != '-' && next != 'n' )
                  {
                    SetAbacus(abacus, i, j  , next);
                    SetAbacus(abacus, i, j+1, curr);
                  }
              }
              // The entire j+1-th column now contains only gaps or n's
              // Remove it by shifting all the subsequent columns
              // one position to the left
              for (int32 i=0;i<abacus->rows;i++)
                {
                  curr = *GetAbacus(abacus,i,j  );
                  next = *GetAbacus(abacus,i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j+1; k<last_non_null; k++)
                    {
                      next= *GetAbacus(abacus,i,k+1);
                      SetAbacus(abacus, i, k, next);
                    }
                  SetAbacus(abacus, i, last_non_null, '-');
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
          for (int32 i=0;i<abacus->rows;i++)
            {
              curr = *GetAbacus(abacus,i,j);
              next = *GetAbacus(abacus,i,j+1);
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
              for (int32 i=0;i<abacus->rows;i++) {
                curr = *GetAbacus(abacus,i,j  );
                next = *GetAbacus(abacus,i,j+1);
                if (curr == 'n' && next == 'n')
                  {
                    continue;
                  }
                if (curr != '-' && curr != 'n' ) {
                  SetAbacus(abacus, i, j  , next);
                  SetAbacus(abacus, i, j+1, curr);
                }
              }
              // The entire j-th column contains gaps
              // Remove it by shifting all the previous columns
              // one position to the right
              for (int32 i=0;i<abacus->rows;i++)
                {
                  curr = *GetAbacus(abacus,i,j  );
                  next = *GetAbacus(abacus,i,j+1);
                  if (curr == 'n' && next == 'n')
                    continue;
                  for (int32 k=j; k>first_non_null; k--)
                    {
                      prev = *GetAbacus(abacus,i,k-1);
                      SetAbacus(abacus, i, k, prev);
                    }
                  SetAbacus(abacus, i, first_non_null, '-');
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

static
void
RefineOrigAbacus(AbacusDataStructure *abacus, VarRegion vreg) {

  ResetCalls(abacus);

  for (int32 j=abacus->window_width; j<2*abacus->window_width; j++)
    {
      // Look through confirmed alleles only
      for (int32 k=0; k<vreg.nca; k++)
        {
          // Line number corresponding to the last read in allele
          int32 il = vreg.alleles[k].read_ids[vreg.alleles[k].num_reads-1];
          char  c  = *GetAbacus(abacus, il, j);
          for (int32 l=0; l<vreg.alleles[k].num_reads-1; l++)
            {
              int32 i = vreg.alleles[k].read_ids[l];
              SetAbacus(abacus, i, j, c);
            }
        }
    }
}

static
int32
LeftShift(AbacusDataStructure *abacus, VarRegion vreg, int32 *lcols) {
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
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
  for (j=abacus->window_width; j<2*abacus->window_width; j++)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              c = *GetAbacus(abacus, i, j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (*GetAbacus(abacus, i, ccol-1) == '-') {
                    ccol--;
                  }
                  // from ccol back up to j, look for column with matching call
                  for (pcol = ccol; pcol<j; pcol++)
                    {
                      call = abacus->calls[pcol];
                      if ( call != 'n' && call != c && c != 'n') {
                        // GD: consensus in a column == '-' ?
                        continue;
                      }
                      if ( call == 'n') {
                        // GD:
                        // Found the leftmost column with non-gap consensus.
                        // Now, reset its consensus "dynamically" to the
                        // current base
                        // Potential problem: the result will generally
                        // depend on the order in which rows
                        // are processed
                        abacus->calls[pcol] = c;
#if 0
                        fprintf(stderr, "j= %d i= %d calls[%d]= %c\n", j, i, pcol, c);
#endif
                      }
                      if (abacus->calls[pcol] == c || c == 'n') {
                        // swap bases in columns pcol and j of row i
                        SetAbacus(abacus, i, j, '-');
                        SetAbacus(abacus, i, pcol, c);
                        break;
                      }
                    }
                  if (*GetAbacus(abacus, i, j) != '-')
                    abacus->calls[j] = c;
                }
            }
        }
    }
#if 0
  fprintf(stderr, "Test calls=\n");
  for (j=0;j<abacus->columns;j++)
    fprintf(stderr, "%c", abacus->calls[j]);
  fprintf(stderr, "\n");
#endif
#ifdef DEBUG_ABACUS
  fprintf(stderr, "Abacus after LeftShift before Merge:\n");
  ShowAbacus(abacus);
#endif
  MergeAbacus(abacus, -1);
#if 0
  fprintf(stderr, "Abacus after Merge:\n");
  ShowAbacus(abacus);
#endif
  abacus->shift = LEFT_SHIFT;
  return ScoreAbacus(abacus,lcols);
}

static
int32
RightShift(AbacusDataStructure *abacus, VarRegion vreg, int32 *rcols) {
 // rcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
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

  for (j=2*abacus->window_width-1;j>abacus->window_width-1;j--)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while ( *GetAbacus(abacus,i,ccol+1) == '-' )
                    ccol++;
                  // now, from ccol back down to j, look for column with matching call
                  for ( pcol = ccol;pcol>j;pcol--)
                    {
                      call = abacus->calls[pcol];
                      if ( call != 'n' && call != c && c != 'n' )
                        continue;
                      if ( call == 'n')
                        abacus->calls[pcol] = c;
                      if (abacus->calls[pcol] == c || c == 'n' ) {
                        SetAbacus(abacus,i,j,'-');
                        SetAbacus(abacus,i,pcol,c);
                        break;
                      }
                    }
                  if ( *GetAbacus(abacus,i,j) != '-' )
                    abacus->calls[j] = c;
                }
            }
        }
    }
#ifdef DEBUG_ABACUS
  fprintf(stderr, "Abacus after RightShift before Merge:\n");
  ShowAbacus(abacus);
#endif
  MergeAbacus(abacus, 1);
  abacus->shift = RIGHT_SHIFT;
  return ScoreAbacus(abacus,rcols);
}

static
int32
MixedShift(AbacusDataStructure *abacus, int32 *mcols, VarRegion  vreg, int32 lpos, int32 rpos,
                 char *tmpl, int32 long_allele, int32 short_allele) {
  // lcols is the number of non-null columns in result
  int32 i, j, k, l, ccol, pcol;
  char c, call;
  ResetCalls(abacus);
  int32 window_beg, window_end;
  int32 shift =0;

  if (abacus->shift == LEFT_SHIFT)
    {
      window_beg = 0;
      window_end = abacus->window_width;
    }
  else if (abacus->shift == UNSHIFTED)
    {
      window_beg = abacus->window_width;
      window_end = 2* abacus->window_width;
    }
  else
    {
      window_beg = 2*abacus->window_width;
      window_end = 3*abacus->window_width;
    }

  /* Populate calls */
  for (j=window_beg; j<window_end; j++)
    abacus->calls[j] = tmpl[j];

  /* Perform left shift */
  for (j=window_beg;j<=MIN(window_end, lpos);j++)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the left for a suitable placement
                  // will be safe on left since abacus has 'n' border
                  while (( *GetAbacus(abacus,i,ccol-1) == '-' ) &&
                         (ccol > window_beg)) {
                    ccol--;
                  }
                  // now, from ccol back up to j, look for column with matching call
                  for ( pcol = ccol;pcol<j;pcol++) {
                    call = abacus->calls[pcol];
                    if ( call != 'n' && call != c && c != 'n')
                      // GD: consensus in a column == '-' ?
                      continue;

                    if ( call == 'n') {
                      // GD: found the leftmost column with non-gap consensus =>
                      //     reset it consensus "dynamically" to the current base
                      //     Potential problem: this code is biased  in the sense that
                      //     the result will generally depend on the order in which
                      //     reads i(or rows) are processed
                      abacus->calls[pcol] = c;
                    }
                    if (abacus->calls[pcol] == c || c == 'n') {
                      // swap bases in columns pcol and j of row i
                      SetAbacus(abacus,i,j,'-');
                      SetAbacus(abacus,i,pcol,c);
                      break;
                    }
                  }
                  if ( *GetAbacus(abacus,i,j) != '-' ) {
                    abacus->calls[j] = c;
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
    fprintf(stderr, "%c", abacus->calls[i]);
  fprintf(stderr, "\n");
#endif

  /* Perform right shift */
  for (j=window_end-1;j>(rpos>0?rpos:window_end);j--)
    {
      for (k=0; k<vreg.na; k++)
        {
          for (l=0; l<vreg.alleles[k].num_reads; l++)
            {
              i = vreg.alleles[k].read_ids[l];
              // Only reads from short allele shouls be shifted
#if 0
              fprintf(stderr, "i=%d vreg.alleles[i]=%d short_allele=%d\n",
                      i, vreg.alleles[i], short_allele);
#endif
              if (vreg.alleles[i].id != short_allele)
                continue;

              c = *GetAbacus(abacus,i,j);
              ccol = j;
              if ( c != '-' )
                {
                  //look to the right for a suitable placement
                  // will be safe on right since abacus has 'n' border
                  while (( *GetAbacus(abacus,i,ccol+1) == '-') &&
                         (ccol+1<window_end) )
                    ccol++;
#if 0
                  fprintf(stderr, "ccol=%d\n", ccol);
#endif
                  // now, from ccol back down to j, look for column with matching call
                  for ( pcol = ccol;pcol>j;pcol--) {
                    call = abacus->calls[pcol];
#if 0
                    fprintf(stderr, "i=%d j=%d c=%c pcol=%d call=%d \n",
                            i, j, c, pcol, call);
#endif
                    if ( call != 'n' && call != c && c != 'n' ) {
                      continue;
                    }
                    if ( call == 'n') {
                      abacus->calls[pcol] = c;
                    }
#if 0
                    fprintf(stderr, "abacus->calls=%c c=%c\n", abacus->calls, c);
#endif
                    if (abacus->calls[pcol] == c || c == 'n' ) {
#if 0
                      fprintf(stderr, "Swapping elements (%d, %d)=%c  and (%d, %d)='-'\n",
                              i, j, c, i, pcol);
#endif
                      SetAbacus(abacus,i,j,'-');
                      SetAbacus(abacus,i,pcol,c);
                      break;
                    }
                  }
                  if ( *GetAbacus(abacus,i,j) != '-' )
                    abacus->calls[j] = c;
                }
            }
        }
    }
  // MergeAbacus(abacus, 1);
  abacus->shift = MIXED_SHIFT;
  return ScoreAbacus(abacus,mcols);
}



static
void
LeftEndShiftBead(beadIdx bid, beadIdx eid) {
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

  Bead *shift = GetBead(beadStore,eid);
  beadIdx aid = (GetBead(beadStore,bid))->prev;
  assert(shift != NULL);
  if ( *Getchar(sequenceStore,shift->soffset) != '-' ) {
    // assume first and internal characters are gaps
    LateralExchangeBead(bid, eid);
  } else {
    while ( shift->prev != aid ) {
      LateralExchangeBead(shift->prev, shift->boffset);
    }
  }
}

static
void
RightEndShiftBead(beadIdx bid, beadIdx eid) {
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

  Bead *shift = GetBead(beadStore,bid);
  beadIdx aid = (GetBead(beadStore,eid))->next;
  beadIdx rid;
  assert(shift != NULL);
  if ( *Getchar(sequenceStore,shift->soffset) != '-' ) {
    // assume last and internal characters are gaps
    LateralExchangeBead(bid, eid);
  } else {
    rid = shift->next;
    while ( shift->next != aid ) {
      LateralExchangeBead(shift->boffset, shift->next);
    }
  }
}



static
void
GetAbacusBaseCount(AbacusDataStructure *a, BaseCount *b) {
  int32 j;
  ResetBaseCount(b);
  for (j=0;j<a->columns;j++) {
    IncBaseCount(b,a->calls[j]);
  }
}

static
int
ColumnMismatch(Column *c) {
  char maxchar =  GetMaxBaseCount(&c->base_count,0);
  return c->base_count.depth - c->base_count.count[RINDEX[maxchar]];
}

static
char
GetBase(int32 s) {
  return *Getchar(sequenceStore,s);
}

static
void
ApplyAbacus(AbacusDataStructure *a, CNS_Options *opp) {
  Column    *column;
  int32        columns=0;
  char       a_entry;
  double     fict_var;   // variation is a column
  VarRegion  vreg;

  beadIdx bid;  //  ALWAYS the id of bead
  beadIdx eid;  //  ALWAYS the id of exch

  Bead *bead = NULL;
  Bead *exch = NULL;

  char base;

  vreg.nr = 0;

  if (a->shift == LEFT_SHIFT) {
    column = GetColumn(columnStore,a->start_column);
    assert(column != NULL);

    while (columns < a->window_width) {
      bid = GetBead(beadStore,column->call)->down;
#ifdef DEBUG_APPLYABACUS
      fprintf(stderr, "0; bid=%d eid=%d\n", bid, eid);
#endif

      // Update all beads in a given column

      while (bid.isValid()) {
        bead = GetBead(beadStore,bid);
        a_entry = *GetAbacus(a, *Getint32(abacus_indices,bead->frag_index) - 1, columns);

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr, "a_entry=%c bead=%c\n", a_entry, *Getchar(sequenceStore,bead->soffset));
#endif

        if (a_entry == 'n') {
          eid  = bead->up;
          exch = GetBead(beadStore, eid);
#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "1; bid=%d eid=%d\n", bid, eid);
#endif
          UnAlignTrailingGapBeads(bid);

        } else if (a_entry != *Getchar(sequenceStore,bead->soffset)) {
          //  Look for matching bead in frag and exchange
          eid  = bead->boffset;
          exch = GetBead(beadStore,eid);

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "2; bid=%d eid=%d\n", bid, eid);
#endif

          if (NULL == exch) {
            eid = AppendGapBead(bead->boffset);
            bead = GetBead(beadStore, bid);
            AlignBeadToColumn(GetColumn(columnStore,bead->column_index)->next,eid, "ApplyAbacus(1)");
            exch = GetBead(beadStore, eid);
          }

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "3; bid=%d eid=%d\n", bid, eid);
#endif

          while (a_entry != *Getchar(sequenceStore,exch->soffset)) {
            beadIdx eidp = exch->next;

            if (exch->next.isInvalid()) {
              eidp = AppendGapBead(exch->boffset);
              bead = GetBead(beadStore, bid);
              exch = GetBead(beadStore, eid);
              AlignBeadToColumn(GetColumn(columnStore,exch->column_index)->next,eidp, "ApplyAbacus(2)");
#ifdef DEBUG_APPLYABACUS
              fprintf(stderr, "4; bid=%d eid=%d\n", bid, eid);
#endif

            } else if (exch->column_index == a->end_column) {
              eidp = AppendGapBead(exch->boffset);
              bead = GetBead(beadStore, bid);
              exch = GetBead(beadStore, eid);
#ifdef DEBUG_APPLYABACUS
              fprintf(stderr, "5; bid=%d eid=%d\n", bid, eid);
#endif

              int32 cid = column->lid;
              ColumnAppend(exch->column_index,eidp);
              column = GetColumn(columnStore, cid);
            }

#ifdef DEBUG_APPLYABACUS
            fprintf(stderr, "6; bid=%d eid=%d b col/frg=%d/%d e_col/frg=%d/%d\n",
                    bid, eid,
                    bead->column_index, bead->frag_index,
                    exch->column_index, exch->frag_index);
#endif

            eid  = eidp;
            exch = GetBead(beadStore,eid);
          }

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr,"LeftShifting bead %d (%c) with bead %d (%c).\n",
                  bid, *Getchar(sequenceStore,bead->soffset),
                  eid, *Getchar(sequenceStore,exch->soffset));
#endif

          LeftEndShiftBead(bid, eid);
        } else {
          // no exchange necessary;
          eid  = bid;
          exch = bead;

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr, "7; bid=%d eid=%d\n", bid, eid);
#endif
        }

        bid  = exch->down;
        bead = NULL;

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr,"New bid is %d (%c), from %d down\n",
                bid, (bid.isValid()) ? *Getchar(sequenceStore,GetBead(beadStore,bid)->soffset) : 'n',
                eid);
#endif
      }

      //  End of update; call base now.

      BaseCall(column->lid, 1, fict_var, &vreg, -1, base, 0, opp);

      column = GetColumn(columnStore,column->next);
      columns++;
    }
  }


  if (a->shift == RIGHT_SHIFT) {
    column = GetColumn(columnStore,a->end_column);
    assert(column != NULL);

    while (columns<a->window_width) {
      bid = GetBead(beadStore,column->call)->down;

      while (bid.isValid()) {
        bead = GetBead(beadStore,bid);
        a_entry = *GetAbacus(a, *Getint32(abacus_indices,bead->frag_index) - 1, a->columns-columns-1);

        if (a_entry == 'n') {
          eid  = bead->up;
          exch = GetBead(beadStore, eid);
          UnAlignTrailingGapBeads(bid);
        } else if (a_entry != *Getchar(sequenceStore,bead->soffset)) {
          //  Look for matching bead in frag and exchange
          eid  = bead->boffset;
          exch = GetBead(beadStore, eid);

          if (NULL == exch) {
            eid  = PrependGapBead(bead->boffset);
            bead = GetBead(beadStore, bid);
            exch = GetBead(beadStore, eid);
            AlignBeadToColumn(GetColumn(columnStore,bead->column_index)->prev,eid, "ApplyAbacus(3)");
          }

          while (a_entry != *Getchar(sequenceStore,exch->soffset)) {
            beadIdx eidp = exch->prev;

            if (exch->prev.isInvalid()) {
              eidp = PrependGapBead(exch->boffset);
              bead = GetBead(beadStore, bid);
              exch = GetBead(beadStore, eid);
              AlignBeadToColumn(GetColumn(columnStore,exch->column_index)->prev,eidp, "ApplyAbacus(4)");
            } else if (exch->column_index == a->start_column) {
              eidp = AppendGapBead(exch->prev);
              bead = GetBead(beadStore, bid);
              exch = GetBead(beadStore, eid);

              int32 cid = column->lid;
              ColumnAppend(GetColumn(columnStore,exch->column_index)->prev,eidp);
              column = GetColumn(columnStore, cid);
            }

            eid  = eidp;
            exch = GetBead(beadStore, eid);
          }

#ifdef DEBUG_APPLYABACUS
          fprintf(stderr,"RightShifting bead %d (%c) with bead %d (%c).\n",
                  eid, *Getchar(sequenceStore,exch->soffset),
                  bid, *Getchar(sequenceStore,bead->soffset));
#endif

          RightEndShiftBead(eid, bid);
        } else {
          eid  = bid;
          exch = bead; // no exchange necessary;
        }

        bid  = exch->down;
        bead = NULL;

#ifdef DEBUG_APPLYABACUS
        fprintf(stderr,"New bid is %d (%c), from %d down\n",
                bid, (bid>-1)?*Getchar(sequenceStore,GetBead(beadStore,bid)->soffset):'n',
                eid);
#endif
      }

      BaseCall(column->lid, 1, fict_var, &vreg, -1, base, 0, opp);
      column = GetColumn(columnStore,column->prev);
      columns++;
    }
  }
}

static
int
IdentifyWindow(Column **start_column, int32 *stab_bgn, CNS_RefineLevel level) {
  Column *stab;
  Column *pre_start;
  int32 win_length=1;
  int32 rc=0;
  int32 gap_count=0;
  char poly;
  *stab_bgn = (*start_column)->next;
  stab = GetColumn(columnStore,*stab_bgn);
  switch (level) {
    case CNS_SMOOTH:
      // in this case, we just look for a string of gaps in the consensus sequence
      if ( GetBase( GetBead(beadStore,(*start_column)->call)->soffset ) != '-' ) break;
      // here, there's a '-' in the consensus sequence, see if it expands
      while( GetBase( GetBead(beadStore,stab->call)->soffset) == '-' )  {
        // move stab column ahead
        if ( stab->next != -1 ) {
          *stab_bgn = stab->next;
          stab = GetColumn(columnStore,*stab_bgn);
          win_length++;
        } else {
          break;
        }
      }
      if ( win_length > 1 ) rc = win_length;
      break;
    case CNS_POLYX:
      // here, we're looking for a string of the same character
      gap_count=GetColumnBaseCount(*start_column,'-');
      poly =  GetBase(GetBead(beadStore,(*start_column)->call)->soffset);
      if ( poly != '-' ) {
        char cb;

        while( (cb = GetBase(GetBead(beadStore,stab->call)->soffset)) == poly || cb == '-' )  {
          // move stab column ahead
          if ( stab->next != -1 ) {
            *stab_bgn = stab->next;
            gap_count+=GetColumnBaseCount(stab,'-');
            stab = GetColumn(columnStore,*stab_bgn);
            win_length++;
          } else {
            break;
          }
        }
        // capture trailing gap-called columns
        if ( win_length > 2 ) {
          while( GetBase(GetBead(beadStore,stab->call)->soffset) == '-' )  {
            if ( GetMaxBaseCount(&stab->base_count,1) != poly ) break;
            if ( stab->next != -1 ) {
              *stab_bgn = stab->next;
              gap_count+=GetColumnBaseCount(stab,'-');
              stab = GetColumn(columnStore,*stab_bgn);
              win_length++;
            } else {
              break;
            }
          }
          // now that a poly run with trailing gaps is established, look for leading gaps
          pre_start = *start_column;
          while ( pre_start->prev != -1 ) {
            char cb;
            pre_start = GetColumn(columnStore,pre_start->prev);
            if ( (cb = GetBase(GetBead(beadStore,pre_start->call)->soffset)) != '-' && cb != poly ) break;
            *start_column = pre_start;
            gap_count+=GetColumnBaseCount(pre_start,'-');
            win_length++;
          }
        } else {
          break;
        }
      }
      if ( (*start_column)->prev != -1 && win_length > 2 && gap_count > 0) {
        //fprintf(stderr,"POLYX candidate (%c) at column %d, width %d, gapcount %d\n", poly,(*start_column)->ma_index,win_length,gap_count);
        rc = win_length;
      }
      break;
    case CNS_INDEL:
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
        Column *stab_end;

        cum_mm = ColumnMismatch(*start_column);
        if ( cum_mm > 0 && GetColumnBaseCount(*start_column,'-') > 0) {
          stab = *start_column;
          stab = GetColumn(columnStore,(*start_column)->next);
          stab_end = stab;
          while ( stab_end->next != -1 && stab_width < STABWIDTH) {
            stab_mm+=ColumnMismatch(stab_end);
            stab_gaps+=GetColumnBaseCount(stab_end,'-');
            stab_bases+=GetDepth(stab_end);
            stab_end = GetColumn(columnStore,stab_end->next);
            stab_width++;
          }
          if ( stab_bases == 0 ) break;
          //  Floating point 'instability' here?
          while( (double)stab_mm/(double)stab_bases > 0.02 ||  //  CNS_SEQUENCING_ERROR_EST
                 (double)stab_gaps/(double)stab_bases > .25  ){
            int32 mm=ColumnMismatch(stab);
            int32 gp=GetColumnBaseCount(stab,'-');
            int32 bps=GetDepth(stab);
            // move stab column ahead
            if ( stab_end->next != -1 ) {
              stab_mm+=ColumnMismatch(stab_end);
              stab_bases+=GetDepth(stab_end);
              stab_gaps+=GetColumnBaseCount(stab_end,'-');
              stab_end = GetColumn(columnStore,stab_end->next);
              stab_mm-=mm;
              stab_gaps-=gp;
              stab_bases-=bps;
              cum_mm+=mm;
              stab = GetColumn(columnStore,stab->next);
              win_length++;
            } else {
              break;
            }
          }
          *stab_bgn = stab->lid;
        }
        if ( win_length > 1 ) rc = win_length;
      }
      break;
    default:
      break;
  }
  return rc;
}



static
void
ShowCalls(AbacusDataStructure *abacus) {
  for (int32 j=0;j<abacus->columns;j++)
    fprintf(stderr, "%c", abacus->calls[j]);
  fprintf(stderr, "\n");
}

static
 void
GetReadsForAbacus(Read *reads, AbacusDataStructure *abacus) {
  int32 i, j, shift=0;
  char base;

#if 0
  fprintf(stderr, "rows=%lu shift=%c window_width=%lu \nReads= \n",
          abacus->rows, (char)abacus->shift, abacus->window_width);
#endif
  if (abacus->shift == UNSHIFTED)
    shift = abacus->columns;
  else if (abacus->shift == RIGHT_SHIFT)
    shift = 2*abacus->columns;
  for (i=0; i<abacus->rows; i++) {
    reads[i].id = i;
    reads[i].allele_id = -1;
    reads[i].ave_qv = 20.;        // qvs are hardly available for abacus
    reads[i].uglen = 0;
    for (j=0; j<abacus->columns; j++) {
      base = *GetAbacus(abacus,i,j);
      if (is_good_base(base))
        reads[i].bases[j] = base;
      if (base != '-')
        reads[i].uglen++;
    }
  }
}



static
 void
GetConsensusForAbacus(VarRegion  *vreg, Read *reads, AbacusDataStructure *abacus,
                      char ***consensus) {
  char bases[CNS_NALPHABET] = {'-', 'A', 'C', 'G', 'T', 'N'};
  // Allocate memory for consensus
  *consensus = (char **)safe_malloc(2 * sizeof(char *));
  for (int32 i=0; i<2; i++) {
    (*consensus)[i] = (char *)safe_malloc(3*abacus->window_width * sizeof(char));
    for (int32 j=0; j<3*abacus->window_width; j++)
      (*consensus)[i][j] = '-';
  }

  // Call consensus
  for (int32 i=0; i<3*abacus->window_width; i++)
    {
      int32 bcount0[CNS_NALPHABET] = {0};
      int32 bcount1[CNS_NALPHABET] = {0};
      int32 best_count0=0, second_best_count0=0;
      int32 best_count1=0, second_best_count1=0;
      char cbase0=0, cbase1=0;
      for (int32 j=0; j<abacus->rows; j++) {
#if 0
        fprintf(stderr, " reads[%d][%d]= %c\n", j, i, reads[j].bases[i]);
#endif
        if (is_good_base(reads[j].bases[i]))
          {
            if   (vreg->alleles[j].id == 0)
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
        }
        else if (  best_count0 >= bcount0[j] &&
                   second_best_count0 <  bcount0[j]) {
          second_best_count0  = bcount0[j];
        }
      }
      for (int32 j=0; j<CNS_NALPHABET; j++) {
        if (best_count1 < bcount1[j]) {
          second_best_count1 = best_count1;
          best_count1 = bcount1[j];
          cbase1 = bases[j];
        }
        else if (  best_count1 >= bcount1[j] &&
                   second_best_count1 <  bcount1[j]) {
          second_best_count1  = bcount1[j];
        }
      }
      if (best_count0 == second_best_count0)
        (*consensus)[0][i] = 'N';
      else
        (*consensus)[0][i] = cbase0;
      if (best_count1 == second_best_count1)
        (*consensus)[1][i] = 'N';
      else
        (*consensus)[1][i] = cbase1;
    }
}

/* Create ungapped consensus sequences and map them to gapped consensus sequences */
static
 void
MapConsensus(int32 ***imap, char **consensus,  char ***ugconsensus,
             int32 len, int32 *uglen) {
  uglen[0] = uglen[1] = 0;
  *ugconsensus = (char **)safe_malloc(2*sizeof(char *));
  *imap        = (int32  **)safe_malloc(2*sizeof(int32  *));
  for (int32 i=0; i<2; i++)
    {
      (*ugconsensus)[i] = (char *)safe_malloc(len*sizeof(char));
      (*imap)[i]        = (int32  *)safe_malloc(len*sizeof(int32 ));
      for (int32 j=0; j<len; j++)
        (*imap)[i][j] = j;
      int32 k=0;
      for (int32 j=0; j<len; j++)
        {
          if (consensus[i][j] != '-')
            {
              (*ugconsensus)[i][k] = consensus[i][j];
              (*imap)[i][k] = j;
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
  GetTemplateForAbacus(char **tmpl, char **consensus, int32 len,
                           char **ugconsensus, int32 *uglen, int32 lpos, int32 rpos, int32 **imap,
                           int32 *adjleft, int32 *adjright, int32 short_allele, int32 long_allele) {
  int32 i, j;

  *tmpl = (char *)safe_malloc(len*sizeof(char));
  for (i=0; i<len; i++)
    (*tmpl)[i] = consensus[long_allele][i];

  /* Set Ns in the left part of the tmpl */
  i = 0;
  while ((imap[long_allele][i] <= lpos) &&
         (i < uglen[short_allele] - adjleft[short_allele]) &&
         (i < uglen[ long_allele] - adjleft[ long_allele]))
    {
      int32 lpos = i + adjleft[long_allele];
      int32 spos = i + adjleft[short_allele];
      if ((ugconsensus[short_allele][spos] != ugconsensus[long_allele][lpos]) &&
          ((*tmpl)[imap[long_allele][lpos]] != '-'))
        (*tmpl)[imap[long_allele][lpos]] = 'n';
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
              (*tmpl)[i] = 'n';
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
          ((*tmpl)[imap[long_allele][i]] != '-'))
        {
          (*tmpl)[imap[long_allele][i]] =  'n';
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
              (*tmpl)[i] = 'n';
            }
        }
    }
}




static
int
 RefineWindow(MANode *ma, Column *start_column, int32 stab_bgn,
                 CNS_Options *opp ) {
  int32 orig_columns=0, left_columns=0, right_columns=0, best_columns=0;
  // Mismatch, gap and total scores:
  int32 orig_mm_score=0, left_mm_score=0, right_mm_score=0, best_mm_score=0;
  int32 orig_gap_score=0, left_gap_score=0, right_gap_score=0, best_gap_score = 0;
  int32 orig_total_score, left_total_score, right_total_score, best_total_score;
  int32 max_element = 0, score_reduction = 0;
  BaseCount abacus_count;
  AbacusDataStructure *left_abacus, *orig_abacus, *right_abacus, *best_abacus;
  VarRegion  vreg;

  orig_abacus = CreateAbacus(ma->lid,start_column->lid,stab_bgn);
  vreg.nr = orig_abacus->rows;
  //ShowAbacus(orig_abacus);
  orig_mm_score = ScoreAbacus(orig_abacus,&orig_columns);

  // Process reads of the original abacus
  AllocateDistMatrix(&vreg, 0);
  AllocateMemoryForReads(&vreg.reads, orig_abacus->rows, orig_abacus->columns,
                         QV_FOR_MULTI_GAP);
  vreg.beg = 0;
  vreg.end = orig_abacus->columns-1;
  GetReadsForAbacus(vreg.reads, orig_abacus);
#ifdef DEBUG_VAR_RECORDS
  OutputReads(stderr, vreg.reads, vreg.nr, orig_abacus->columns);
#endif
  PopulateDistMatrix(vreg.reads, orig_abacus->columns, &vreg);
#ifdef DEBUG_VAR_RECORDS
  OutputDistMatrix(stderr, &vreg);
#endif
  AllocateMemoryForAlleles(&vreg.alleles, vreg.nr, &vreg.na);
  ClusterReads(vreg.reads, vreg.nr, vreg.alleles, &vreg.na,
               &vreg.nca, vreg.dist_matrix);
  SortAllelesByLength(vreg.alleles, vreg.nca, vreg.reads);

  RefineOrigAbacus(orig_abacus, vreg);
  //ShowAbacus(orig_abacus);
  orig_mm_score = ScoreAbacus(orig_abacus, &orig_columns);

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
  left_abacus = CloneAbacus(orig_abacus);
  left_mm_score = LeftShift(left_abacus, vreg, &left_columns);
#ifdef DEBUG_ABACUS
  fprintf(stderr, "\n\nLeftShiftCalls=\n");
  ShowCalls(left_abacus);
  fprintf(stderr, "Abacus=\n");
  ShowAbacus(left_abacus);
  fprintf(stderr, "\n");
#endif
  right_abacus = CloneAbacus(orig_abacus);
  right_mm_score = RightShift(right_abacus, vreg, &right_columns);
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
  orig_gap_score  = AffineScoreAbacus(orig_abacus);
  left_gap_score  = AffineScoreAbacus(left_abacus);
  right_gap_score = AffineScoreAbacus(right_abacus);
  best_abacus     = orig_abacus;
  best_columns    = orig_columns;
  best_gap_score  = orig_gap_score;
  best_mm_score   = orig_mm_score;
  orig_total_score  = orig_mm_score  + orig_columns  + orig_gap_score;
  left_total_score  = left_mm_score  + left_columns  + left_gap_score;
  right_total_score = right_mm_score + right_columns + right_gap_score;
  best_total_score  = orig_total_score;

#ifdef DEBUG_ABACUS
  fprintf(stderr, "In RefineWindow: beg= %lu end= %d\n",
          start_column->lid, stab_bgn);
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
  if ( left_total_score < orig_total_score ||
       right_total_score < orig_total_score )
    {
      if ( left_total_score <= right_total_score ) {
        score_reduction += orig_total_score - left_total_score;
        //fprintf(stderr,"\nTry to apply LEFT abacus:\n");
        //ShowAbacus(left_abacus);
        GetAbacusBaseCount(left_abacus, &abacus_count);
#if 0
        fprintf(stderr, " Applying left abacus\n");
#endif
        best_abacus      = left_abacus;
        best_mm_score    = left_mm_score;
        best_columns     = left_columns;
        best_gap_score   = left_gap_score;
        best_total_score = left_total_score;
      }
      else
        {
          score_reduction += orig_total_score - right_total_score;
          //fprintf(stderr,"\nTry to apply RIGHT abacus:\n");
          //ShowAbacus(right_abacus);
          GetAbacusBaseCount(right_abacus,&abacus_count);
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
    char  **consensus=NULL, **ugconsensus=NULL, *tmpl=NULL;
    int32   **imap=NULL, uglen[2]={0,0}, adjleft[2]={-1,-1}, adjright[2]={-1,-1};
    int32     gapcount[2], short_allele=-1, long_allele=-1;
    int32     lscore=0, rscore=0, lpos=-1, rpos=-1;
    int32     mixed_columns=0;
    int32   mixed_mm_score=0, mixed_gap_score=0;
    AbacusDataStructure *mixed_abacus=NULL;

    GetConsensusForAbacus(&vreg, vreg.reads, best_abacus, &consensus);
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
    if (gapcount[short_allele] == 0)
      {
#if 0
        fprintf(stderr, "No MixedShift will be performed: gapcount[short_allele] = %d\n", gapcount[short_allele]);
#endif
        ApplyAbacus(best_abacus, opp);

        DeleteAbacus(orig_abacus);
        DeleteAbacus(left_abacus);
        DeleteAbacus(right_abacus);
        if (vreg.nr > 0)
          {
            for (int32 j=0; j<vreg.nr; j++)
              {
                safe_free(vreg.alleles[j].read_ids);
                safe_free(vreg.alleles[j].read_iids);
                safe_free(vreg.dist_matrix[j]);
                safe_free(vreg.reads[j].bases);
                safe_free(vreg.reads[j].qvs);
              }
            safe_free(vreg.reads);
            safe_free(vreg.alleles);
            safe_free(vreg.dist_matrix);
          }
        safe_free(consensus[0]);
        safe_free(consensus[1]);
        safe_free(consensus);
        return score_reduction;
      }

    /* Now try the mixed consensus */
    MapConsensus(&imap, consensus, &ugconsensus, 3*best_abacus->window_width,
                 uglen);
    if ((uglen[0] < MSTRING_SIZE) || (uglen[1] < MSTRING_SIZE))
      {
#if 0
        fprintf(stderr, "No MixedShift will be performed: uglen = %d %d\n", uglen[0], uglen[1]);
#endif
        ApplyAbacus(best_abacus, opp);

        DeleteAbacus(orig_abacus);
        DeleteAbacus(left_abacus);
        DeleteAbacus(right_abacus);
        if (vreg.nr > 0)
          {
            for (int32 j=0; j<vreg.nr; j++)
              {
                safe_free(vreg.alleles[j].read_ids);
                safe_free(vreg.alleles[j].read_iids);
                safe_free(vreg.dist_matrix[j]);
                safe_free(vreg.reads[j].bases);
                safe_free(vreg.reads[j].qvs);
              }
            safe_free(vreg.reads);
            safe_free(vreg.alleles);
            safe_free(vreg.dist_matrix);
          }
        safe_free(consensus[0]);
        safe_free(consensus[1]);
        safe_free(consensus);
        for (i=0; i<2; i++)
          {
            safe_free(ugconsensus[i]);
            safe_free(imap[i]);
          }
        safe_free(ugconsensus);
        safe_free(imap);
        return score_reduction;
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
    FindAdjustedLeftBounds(adjleft, ugconsensus, uglen, short_allele,
                           long_allele);
    FindAdjustedRightBounds(adjright, ugconsensus, uglen, short_allele,
                            long_allele);
#if 0
    fprintf(stderr, "Adjusted left bounds 0, 1= %d %d \n", adjleft[0], adjleft[1]);
    fprintf(stderr, "Adjusted right bounds 0, 1= %d %d \n", adjright[0], adjright[1]);
#endif
    GetLeftScore(ugconsensus, uglen, imap, adjleft,
                 short_allele, long_allele, &lscore, &lpos);
    GetRightScore(ugconsensus, uglen, imap, adjright,
                  short_allele, long_allele, &rscore, &rpos);
    AdjustShiftingInterfaces(&lpos, &rpos, lscore, rscore,
                             adjleft, adjright, long_allele, short_allele);
    GetTemplateForAbacus(&tmpl, consensus, 3*best_abacus->window_width,
                         ugconsensus, uglen, lpos, rpos, imap, adjleft, adjright,
                         short_allele, long_allele);

    mixed_abacus = CloneAbacus(orig_abacus);
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
    mixed_mm_score = MixedShift(mixed_abacus, &mixed_columns, vreg, lpos, rpos,
                                tmpl, long_allele, short_allele);

#if 0
    fprintf(stderr, "Mixed abacus=\n");
    ShowAbacus(mixed_abacus);
    fprintf(stderr, "End calling MixedShift\n\n");
    fprintf(stderr, "mixed_mm_score=%d bast_score=%d\n", mixed_mm_score, best_mm_score);
    fprintf(stderr, "mixed_columns=%d best_columns=%d\n", mixed_columns, best_columns);
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d\n", mixed_gap_score, best_gap_score);
#endif
    mixed_gap_score  = AffineScoreAbacus(mixed_abacus);
#if 0
    fprintf(stderr, "mixed_gap_score=%d best_gap_score=%d mixed_columns=%d best_columns=%d mixed_mm_score=%d best_mm_score=%d\n", mixed_gap_score, best_gap_score, mixed_columns, best_columns, mixed_mm_score, best_mm_score);
#endif
    if ( (mixed_gap_score <  best_gap_score) ||
         ((mixed_gap_score == best_gap_score) && (mixed_columns < best_columns))
         ||
         ((mixed_gap_score == best_gap_score) && (mixed_columns == best_columns) &&
          (mixed_mm_score < best_mm_score)))
      {
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
    ApplyAbacus(best_abacus, opp);

    //      fprintf(stderr, "vreg.nr = %d\n", vreg.nr);

    DeleteAbacus(orig_abacus);
    DeleteAbacus(left_abacus);
    DeleteAbacus(right_abacus);
    DeleteAbacus(mixed_abacus);
    {
      safe_free(consensus[0]);
      safe_free(consensus[1]);
      safe_free(consensus);
      safe_free(ugconsensus[0]);
      safe_free(ugconsensus[1]);
      safe_free(ugconsensus);
      safe_free(imap[0]);
      safe_free(imap[1]);
      safe_free(imap);
      safe_free(tmpl);
    }
    if (vreg.nr > 0)
      {
        for (int32 j=0; j<vreg.nr; j++)
          {
            safe_free(vreg.alleles[j].read_ids);
            safe_free(vreg.alleles[j].read_iids);
            safe_free(vreg.dist_matrix[j]);
            safe_free(vreg.reads[j].bases);
            safe_free(vreg.reads[j].qvs);
          }
        safe_free(vreg.reads);
        safe_free(vreg.alleles);
        safe_free(vreg.dist_matrix);
      }
  }
  return score_reduction;
}


//*********************************************************************************
// Abacus Refinement:
//   AbacusRefine contains the logic for sweeping through the multialignment,
//   and identifying candidate windows for refinement.
//   Each window is cast into an abacus, which is left and right shifted.
//   The best resulting base arrangement (if different from input) is then
//   applied to window of the MultiAlignment
//*********************************************************************************

int
 AbacusRefine(MANode *ma, int32 from, int32 to, CNS_RefineLevel level,
                 CNS_Options *opp) {
  // from and to are in ma's column coordinates
  int32 sid, eid, stab_bgn;
  int32 ma_length = GetMANodeLength(ma->lid);
  int32 score_reduction=0;
  int32 orig_length = ma_length;
  int32 refined_length = orig_length;
  Column *start_column;
  int32 i;

  if(from < 0 || from > ma_length-1){
    fprintf(stderr, "AbacusRefine range (from) invalid");
    assert(0);
  }
  if ( to == -1 ) to = ma_length-1;
  if(to <= from || to > ma_length-1){
    fprintf(stderr, "AbacusRefine range (to) invalid");
    assert(0);
  }

  ResetIndex(abacus_indices,GetNumFragments(fragmentStore));
  sid = *Getint32(ma->columnList, from);   // id of the starting column
  eid = *Getint32(ma->columnList, to);     // id of the ending column
  start_column = GetColumn(columnStore,sid);

  while (start_column->lid != eid)
    {
      int32 window_width = IdentifyWindow(&start_column,&stab_bgn, level);
      // start_column stands as the candidate for first column in window
      // look for window start and stop

      if (window_width > 0)
        {
#ifdef DEBUG_ABACUS
          fprintf(stderr, "In AbacusRefine window_width= %d\n", window_width);
#endif
          //
          // refine in window
          if ( start_column->prev == -1 ) {
            // if start_column->prev == -1, insert a gap column for maneuvering room
            beadIdx newbead;
            Bead *firstbead;
            firstbead = GetBead(beadStore,GetBead(beadStore,start_column->call)->down);
            newbead   = AppendGapBead(firstbead->boffset);
            firstbead = GetBead(beadStore,GetBead(beadStore,start_column->call)->down);
            fprintf(stderr,"Adding gapbead "F_U64" after "F_U64" to add abacus room for abacus abutting left of multialignment\n",
                    (uint64)newbead.get(), (uint64)firstbead->boffset.get());
            ColumnAppend(firstbead->column_index,newbead);
          }

          //  if the window is too big, there's likely a polymorphism
          //  that won't respond well to abacus, so skip it
          //
          //  BPW saw crashes with large window_width's (1333, 3252,
          //  1858, 675, 855, 1563, 601, 1102).  The longest
          //  window_width that worked was 573.  Previous versions
          //  used 100 here.  Not sure what it should be.
          //
          if ( window_width < MAX_WINDOW_FOR_ABACUS_REFINE )
            score_reduction += RefineWindow(ma,start_column,stab_bgn, opp);
          
          start_column = GetColumn(columnStore, stab_bgn);
        }
      start_column = GetColumn(columnStore, stab_bgn);
    }
  RefreshMANode(ma->lid, 1, opp, NULL, NULL, 1, 0);
  return score_reduction;
}
