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

#include "eCR.h"
#include "ScaffoldGraph_CGW.h"


void
DumpContigMultiAlignInfo (char *label, MultiAlignT *cma, int contigID) {
  int           i, j;
  int           deletecma = 0;

  if (label) {
    fprintf(debug.diagnosticFP, "\n%s\n", label);
    fprintf(debug.diagnosticFP, "------------------------------------------------------------\n");
    fprintf(debug.diagnosticFP, "in DumpContigMultiAlignInfo, dumping info on contig %8d\n", contigID);
  }

  if (cma == NULL)
    cma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contigID, FALSE);

  fprintf(debug.diagnosticFP, "  contig %8d, strlen(consensus): %9ld\n",
          contigID, strlen(Getchar(cma->consensus, 0)));

  for (i = 0; i < GetNumIntMultiPoss(cma->f_list); i++) {
    IntMultiPos *pos = GetIntMultiPos(cma->f_list,i);
    fprintf(debug.diagnosticFP, "  fragment %8d, bgn: %10d, end: %10d, length: %10d\n",
            pos->ident, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end));
  }

  for (i = 0; i < GetNumIntUnitigPoss(cma->u_list); i++) {
    IntUnitigPos *pos = GetIntUnitigPos(cma->u_list, i);
    NodeCGW_T *unitig = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);

    MultiAlignT  *uma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);

    fprintf(debug.diagnosticFP, "  unitig %8d, bgn: %10d, end: %10d, length: %10d (consensus: %10d)\n",
            unitig->id, pos->position.bgn, pos->position.end,
            abs(pos->position.bgn - pos->position.end),
            strlen(Getchar(uma->consensus, 0)));

    for (j = 0; j < GetNumIntMultiPoss(uma->f_list); j++) {
      IntMultiPos *pos = GetIntMultiPos(uma->f_list, j);
      fprintf(debug.diagnosticFP, "  fragment %8d, bgn: %10d, end: %10d, length: %10d, source: %10d\n",
              pos->ident, pos->position.bgn, pos->position.end, abs(pos->position.bgn - pos->position.end),
              pos->sourceInt);
    }
  }
  fprintf(debug.diagnosticFP, "\n");
}


void
DumpUnitigInfo(char *label, NodeCGW_T *unitig) {
  int            i, j;
  MultiAlignT   *uma;
  IntUnitigPos  *pos;

  if (label) {
    fprintf(debug.diagnosticFP, "\n%s\n", label);
    fprintf(debug.diagnosticFP, "------------------------------------------------------------\n");
    fprintf(debug.diagnosticFP, "in DumpUnitigInfo, dumping info on unitig %8d\n", unitig->id);
  }

  uma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, unitig->id, TRUE);
  pos = GetIntUnitigPos(uma->u_list, i);

  fprintf(debug.diagnosticFP, "  unitig %8d, bgn: %10d, end: %10d, length: %10d (consensus: %10d)\n",
          unitig->id, pos->position.bgn, pos->position.end,
          abs(pos->position.bgn - pos->position.end),
          strlen(Getchar(uma->consensus, 0)));

  for (j = 0; j < GetNumIntMultiPoss(uma->f_list); j++) {
    IntMultiPos *mpos = GetIntMultiPos(uma->f_list, j);
    fprintf(debug.diagnosticFP, "  fragment %8d, bgn: %10d, end: %10d, length: %10d, source: %10d\n",
            mpos->ident, mpos->position.bgn, mpos->position.end, abs(mpos->position.bgn - mpos->position.end),
            mpos->sourceInt);
  }
  fprintf(debug.diagnosticFP, "\n");
}


void
DumpContigUngappedOffsets(char *label, int contigID) {
  int                    numCIs;
  int                    i;
  MultiAlignT           *cma;
  int32                 *offsets;

  static VA_TYPE(int32) *UngappedOffsets = NULL;

  if (label) {
    fprintf(debug.diagnosticFP, "\n%s\n", label);
    fprintf(debug.diagnosticFP, "------------------------------------------------------------\n");
    fprintf(debug.diagnosticFP, "in DumpContigUngappedOffsets, dumping info on contig %8d\n", contigID);
  }

  if(!UngappedOffsets) {
    UngappedOffsets = CreateVA_int32(1000);
  }

  cma = loadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, contigID, FALSE);
  GetMultiAlignUngappedOffsets(cma, UngappedOffsets);
  offsets = Getint32(UngappedOffsets, 0);

  numCIs = GetNumIntUnitigPoss(cma->u_list);

  for(i = 0; i < numCIs ; i++) {
    IntUnitigPos *pos = GetIntUnitigPos(cma->u_list, i);
    NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
    int bgn, end;

    // mp->position is an interval.  We need to subtract one from
    // the upper end of the interval
    if (pos->position.end < pos->position.bgn) {
      bgn = pos->position.bgn - 1;
      end = pos->position.end;
    } else {
      bgn = pos->position.bgn;
      end = pos->position.end - 1;
    }

    fprintf(debug.diagnosticFP, "in DCUO, unitig %d, (bgn, end): (%6d, %6d), offsets[bgn]: %10d, offsets[bgn]: %10d\n",
            node->id, bgn, end, offsets[bgn], offsets[end]);
    fprintf(debug.diagnosticFP, "in DCUO, unitig %d, pos->position.bgn: %10d, pos->position.end: %10d\n",
            node->id, pos->position.bgn, pos->position.end);
  }

  for(i = 0; i < GetNumIntMultiPoss(cma->f_list); i++) {
    IntMultiPos *mp = GetIntMultiPos(cma->f_list, i);
    int fragID = mp->sourceInt; // GetInfoByIID(ScaffoldGraph->iidToFragIndex, mp->ident)->fragIndex;
    // hack next line replaces above
    // int fragID = GetInfoByIID(ScaffoldGraph->iidToFragIndex, mp->ident)->fragIndex;
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
    int bgn, end;

    // mp->position is an interval.  We need to subtract one from
    // the upper end of the interval
    if (mp->position.end < mp->position.bgn) {
        bgn = mp->position.bgn - 1;
        end = mp->position.end;
      } else {
        bgn = mp->position.bgn;
        end = mp->position.end - 1;
      }

    fprintf(debug.diagnosticFP, "in DCUO, contig %8d, frag %10d, mp->position.bgn: %10d, mp->position.end: %10d, "
            "len: %10d, contained: %8d, source: %10d\n",
            contigID, frag->iid, mp->position.bgn, mp->position.end, abs(mp->position.bgn - mp->position.end),
            mp->contained, mp->sourceInt);
  }
}
