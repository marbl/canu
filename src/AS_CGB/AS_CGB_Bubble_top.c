
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
static char CM_ID[] 
= "$Id: AS_CGB_Bubble_top.c,v 1.1.1.1 2004-04-14 13:49:50 catmandew Exp $";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS_global.h"
#include "AS_CGB_all.h"
#include "AS_CGB_Bubble.h"
#include "AS_CGB_Bubble_log.h"
#include "AS_CGB_Bubble_Graph.h"
#include "AS_CGB_Bubble_GraphMethods.h"

IntFragment_ID
AS_CGB_Bubble_topo_sort(BubGraph_t bg, IntFragment_ID *out)
{
  IntFragment_ID q_start = 0, q_end = 0, f, opp_f, num_valid = 0;
  IntEdge_ID e;
  uint16 valid_and_unused = AS_CGB_BUBBLE_E_VALID | AS_CGB_BUBBLE_E_UNUSED;
  BG_E_Iter e_it;

  for (e = 0; e < GetNumEdges(BG_edges(bg)); ++e)
    if (BG_E_isSetFlag(bg, e, AS_CGB_BUBBLE_E_VALID)) {
      num_valid++;
      BG_E_setFlag(bg, e, AS_CGB_BUBBLE_E_UNUSED);
    }

  fprintf(BUB_LOG_G, "  * Found " F_IID " valid edges.\n", num_valid);

  num_valid = 0;
  for (f = 0; f < GetNumFragments(BG_vertices(bg)); ++f) 
    if (BG_V_isSetFlag(bg, f, AS_CGB_BUBBLE_V_VALID)) {
      num_valid++;
      if (BG_inDegree(bg, f, valid_and_unused) == 0) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
	fprintf(BUB_LOG_G, "Adding vertex " F_IID " (" F_IID ") as start.\n", f,
		get_iid_fragment(BG_vertices(bg), f));
#endif
	out[q_end++] = f;
      }
    }
  
  fprintf(BUB_LOG_G, "  * Found " F_IID " valid vertices.\n", num_valid);

  while (q_start < q_end) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
    fprintf(BUB_LOG_G, "PROCESSING vertex " F_IID " (" F_IID ").\n", out[q_start],
	    get_iid_fragment(BG_vertices(bg), out[q_start]));
#endif
    for (e = BGEI_bgn(bg, &e_it, out[q_start], bgeiOut, valid_and_unused);
	 !BGEI_end(&e_it);
	 e = BGEI_next(bg, &e_it, valid_and_unused)) {
      opp_f = BG_getOppositeVertex(bg, e, out[q_start]);
      BG_E_clearFlagSymmetric(bg, e, AS_CGB_BUBBLE_E_UNUSED);
      if (BG_inDegree(bg, opp_f, valid_and_unused) == 0) {
	out[q_end++] = opp_f;
#if AS_CGB_BUBBLE_VERY_VERBOSE
	fprintf(BUB_LOG_G, "Adding vertex " F_IID " (" F_IID ").\n", opp_f,
		get_iid_fragment(BG_vertices(bg), opp_f));
#endif
      }
    }
    q_start++;
  }
  
  if (q_end < num_valid) {
    fprintf(BUB_LOG_G, "  * WARNING: Only processed " F_IID " of " F_IID " vertices!  Cyclic graph!\n", q_end, num_valid);
    return 0;
  }
  else
    return num_valid;
}

