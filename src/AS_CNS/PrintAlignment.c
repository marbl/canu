
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

static char *rcsid = "$Id: PrintAlignment.c,v 1.2 2011-01-03 03:07:16 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


//  Width of PrintAlignment output
int32 ALNPAGEWIDTH=100;


void
PrintAlignment(FILE *print, int32 mid, int32 from, int32 to, CNS_PrintKey what) {
  /*
    Print the columns of MANode mid from column index "from" to column index "to"
    (use 0 and -1 to print all columns)
    here's the intent for the what values;
    CNS_QUIET      = (int)'Q', // quiet,  print nothing
    CNS_STATS_ONLY = (int)'S', // print only 1-line statistic summary
    CNS_ALIGNMENT  = (int)'A', // print the multialignment, sans CNS
    CNS_CONSENSUS  = (int)'C', // print the multialignment, with CNS
    CNS_DOTS       = (int)'D', // print the multialignment, dot format
    CNS_NODOTS     = (int)'N', // print the multialignment, nodot format
    CNS_EDIT_SCORE = (int)'E'  // print the edit score column by column
  */
  MANode *ma = GetMANode(manodeStore,mid);
  int32 ma_length=GetMANodeLength(mid);
  int32 i,num_frags;
#ifdef PRINTUIDS
  int64 *fids;
#else
  int32 *fids;
#endif
  char *types;
  int32 window_start, wi;
  VA_TYPE(char) *sequenceSpace,*qualitySpace;
  char *sequence, *quality;
  char pc;
  FragmentBeadIterator *read_it;
  beadIdx bid;
  Bead *bead;
  Fragment *fragment;
  SeqInterval *positions;
  int32 dots=0;

  if(what == CNS_VIEW_UNITIG)what=CNS_DOTS;
  if (what != CNS_CONSENSUS && what != CNS_DOTS && what != CNS_NODOTS && what != CNS_VERBOSE ) return;
  if (what == CNS_DOTS) dots = 1;
  if (what == CNS_NODOTS) dots = 2;
  if (to == -1 ) {
    to = ma_length;
  }
  if(from < 0 || from > to || to > ma_length){
    fprintf(stderr, "PrintAlignment column range invalid");
    assert(0);
  }
  // now, adjust from column so that start is divisible by 100
  // (purely for convenience in viewing)
  from = ((int) from/100)*100;
  if (((int) to/100)*100 != to ) {
    to = ((int) to/100 + 1)*100;
  } else {
    to = ((int) to/100)*100;
  }

#undef SHOW_MSA_ON_ONE_PAGE
#ifdef SHOW_MSA_ON_ONE_PAGE
  ALNPAGEWIDTH=to-from+1;
#endif

  sequenceSpace = CreateVA_char(ma_length);
  qualitySpace = CreateVA_char(ma_length);
  GetMANodeConsensus(mid,sequenceSpace,qualitySpace);
  sequence = Getchar(sequenceSpace,0);
  quality = Getchar(qualitySpace,0);
  num_frags = GetNumFragments(fragmentStore);
  read_it = (FragmentBeadIterator *) safe_calloc(num_frags,sizeof(FragmentBeadIterator));
#ifdef PRINTUIDS
  fids = (int64 *) safe_calloc(num_frags,sizeof(int64));
#else
  fids = (int32 *) safe_calloc(num_frags,sizeof(int32));
#endif
  types = (char *) safe_calloc(num_frags,sizeof(char));
  positions = (SeqInterval *) safe_calloc(num_frags,sizeof(SeqInterval));
  for (i=0;i<num_frags;i++) {
    int32 bgn_column;
    int32 end_column;
    fragment = GetFragment(fragmentStore,i);
    if ( fragment->deleted || fragment->manode != mid) {
      fids[i] = 0;
      continue;
    }
    bgn_column = (GetBead(beadStore,fragment->firstbead                         ))->column_index;
    end_column = (GetBead(beadStore,fragment->firstbead.get()+fragment->length-1))->column_index;
#ifdef PRINTUIDS
    if(fragment->type==AS_READ){
      fids[i] = fragment->uid;
    } else {
      fids[i] = fragment->iid;
    }
#else
    fids[i] = fragment->iid;
#endif
    types[i] = fragment->type;
    if ( bgn_column > -1 && end_column > -1 ) {
      positions[i].bgn = GetColumn(columnStore,bgn_column)->ma_index;
      positions[i].end = GetColumn(columnStore, end_column)->ma_index+1;
    }
    NullifyFragmentBeadIterator(&read_it[i]);
  }
  window_start = from;
  fprintf(print,"\n\n================  MultiAlignment ID %d ==================\n\n",ma->iid);
  while ( window_start < to ) {

    fprintf(print,"\n%d\n%-*.*s <<< consensus\n",window_start,ALNPAGEWIDTH,ALNPAGEWIDTH,&sequence[window_start]);
    fprintf(print,"%-*.*s <<< quality\n\n",ALNPAGEWIDTH,ALNPAGEWIDTH,&quality[window_start]);
    for (i=0;i<num_frags;i++) {
      if ( fids[i] == 0 ) continue;
      for (wi = window_start;wi<window_start+ALNPAGEWIDTH;wi++) {
        if ( IsNULLIterator(&read_it[i]) ) {
          if ( positions[i].bgn < wi && positions[i].end > wi ) {
            CreateFragmentBeadIterator(i,&read_it[i]);

            bid = NextFragmentBead(&read_it[i]);
            while ( GetColumn(columnStore,(bead=GetBead(beadStore,bid))->column_index)->ma_index < wi ) {
              bid = NextFragmentBead(&read_it[i]);
            }
            if (bid.isValid()) {
              pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);
              if (dots == 1) {
                // check whether matches consensus, and make it a dot if so
                if (pc == sequence[wi]) pc = '.';
              }
              if (dots == 2) {
                if (pc == sequence[wi]) pc = ' ';
              }
              fprintf(print,"%c",tolower(pc));
            }
          } else if ( positions[i].bgn ==  wi ) {
            CreateFragmentBeadIterator(i,&read_it[i]);
          } else if ( positions[i].bgn > window_start &&  positions[i].bgn < window_start+ALNPAGEWIDTH) {
            fprintf(print," ");
          } else if ( positions[i].end >= window_start &&  positions[i].end < window_start+ALNPAGEWIDTH) {
            fprintf(print," ");
          } else {
            break;
          }
        }
        if ( ! IsNULLIterator(&read_it[i]) ) {
          bid = NextFragmentBead(&read_it[i]);
          if (bid.isValid()) {
            pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);
            if (dots == 1 ) {
              // check whether matches consensus, and make it a dot if so
              if (pc == sequence[wi]) pc = '.';
            }
            if (dots == 2 ) {
              // check whether matches consensus, and make it a dot if so
              if (pc == sequence[wi]) pc = ' ';
            }
            fprintf(print,"%c",tolower(pc));
          } else {
            fprintf(print," ");
            NullifyFragmentBeadIterator(&read_it[i]);
          }
        }
#ifdef PRINTUIDS
        if ( wi == window_start+ALNPAGEWIDTH - 1 ) fprintf(print," <<< %ld (%c)\n",fids[i],types[i]);
#else
        if ( wi == window_start+ALNPAGEWIDTH - 1 ) fprintf(print," <<< %d (%c)\n",fids[i],types[i]);
#endif
      }
    }
    window_start+=ALNPAGEWIDTH;
  }
  safe_free(read_it);
  safe_free(fids);
  safe_free(types);
  safe_free(positions);
}

