
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

static char *rcsid = "$Id: PrintAlignment.c,v 1.4 2011-12-18 06:53:30 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"
#include "MicroHetREZ.h"
#include "AS_UTL_reverseComplement.h"


void
PrintAlignment(FILE *print, int32 mid, int32 from, int32 to) {
  MANode *ma        = GetMANode(manodeStore,mid);
  int32   ma_length = GetMANodeLength(mid);

  if (from < 0)        from = 0;
  if (to   < 0)        to   = ma_length;
  if (to > ma_length)  to   = ma_length;

  if (from > to)       return;

  int32 pageWidth = to - from + 1;

  VA_TYPE(char) *sequenceSpace = CreateVA_char(ma_length);
  VA_TYPE(char) *qualitySpace  = CreateVA_char(ma_length);

  GetMANodeConsensus(mid, sequenceSpace, qualitySpace);

  char *sequence = Getchar(sequenceSpace,0);
  char *quality  = Getchar(qualitySpace,0);

  int32 num_frags = GetNumFragments(fragmentStore);

  FragmentBeadIterator   *read_it   = (FragmentBeadIterator *) safe_calloc(num_frags, sizeof(FragmentBeadIterator));
  int32                  *fids      = (int32                *) safe_calloc(num_frags, sizeof(int32));
  char                   *types     = (char                 *) safe_calloc(num_frags, sizeof(char));
  SeqInterval            *positions = (SeqInterval          *) safe_calloc(num_frags, sizeof(SeqInterval));

  for (int32 i=0; i<num_frags; i++) {
    Fragment *fragment = GetFragment(fragmentStore,i);

    if ((fragment->deleted == true) || (fragment->manode != mid))
      continue;

    int32 bgn_column = (GetBead(beadStore,fragment->firstbead                         ))->column_index;
    int32 end_column = (GetBead(beadStore,fragment->firstbead.get()+fragment->length-1))->column_index;

    fids[i]  = fragment->iid;
    types[i] = fragment->type;

    if ((bgn_column > -1) &&
        (end_column > -1)) {
      positions[i].bgn = GetColumn(columnStore, bgn_column)->ma_index;
      positions[i].end = GetColumn(columnStore, end_column)->ma_index + 1;
    }

    NullifyFragmentBeadIterator(&read_it[i]);
  }

  int32 window_start = from;

  fprintf(print,"\n\n================  MultiAlignment ID %d ==================\n\n",ma->iid);

  while (window_start < to) {

    fprintf(print,"\n");
    fprintf(print,"%d\n", window_start);
    fprintf(print,"%-*.*s <<< consensus\n", pageWidth, pageWidth, sequence + window_start);
    fprintf(print,"%-*.*s <<< quality\n\n", pageWidth, pageWidth, quality  + window_start);

    for (int32 i=0; i<num_frags; i++) {
      if (fids[i] == 0)
        continue;

      for (int32 wi=window_start; wi<window_start + pageWidth; wi++) {

        if ( IsNULLIterator(&read_it[i]) ) {
          if ((positions[i].bgn < wi) &&
               (positions[i].end > wi)) {
            CreateFragmentBeadIterator(i, &read_it[i]);

            beadIdx   bid = NextFragmentBead(&read_it[i]);
            Bead     *bead;

            while (GetColumn(columnStore,(bead=GetBead(beadStore,bid))->column_index)->ma_index < wi )
              bid = NextFragmentBead(&read_it[i]);

            if (bid.isValid()) {
              char pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);
              
              if (pc == sequence[wi])
                pc = tolower(pc);
              else
                pc = toupper(pc);

              fprintf(print, "%c", pc);
            }

          } else if (positions[i].bgn ==  wi) {
            CreateFragmentBeadIterator(i, &read_it[i]);

          } else if ((positions[i].bgn > window_start) &&
                     (positions[i].bgn < window_start+pageWidth)) {
            fprintf(print," ");

          } else if ((positions[i].end >= window_start) &&
                     (positions[i].end < window_start+pageWidth)) {
            fprintf(print," ");

          } else {
            break;
          }
        }

        if ( ! IsNULLIterator(&read_it[i]) ) {
          beadIdx   bid = NextFragmentBead(&read_it[i]);
          Bead     *bead;

          if (bid.isValid()) {
            char pc = *Getchar(sequenceStore,(GetBead(beadStore,bid))->soffset);

            if (pc == sequence[wi])
              pc = tolower(pc);
            else
              pc = toupper(pc);

            fprintf(print, "%c", pc);

          } else {
            fprintf(print," ");
            NullifyFragmentBeadIterator(&read_it[i]);
          }
        }

        if (wi == window_start + pageWidth - 1)
          fprintf(print," <<< %d (%c)\n", fids[i],types[i]);
      }
    }

    window_start += pageWidth;
  }

  safe_free(read_it);
  safe_free(fids);
  safe_free(types);
  safe_free(positions);
}

