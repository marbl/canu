
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

static const char *rcsid = "$Id: AS_BAT_EvaluateMates.C,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_Unitig.H"
#include "AS_BAT_MateLocation.H"

void
evaluateMates(UnitigVector &unitigs, const char *output_prefix, const char *name) {

  //  [0] -- BOTH frag and mate are dovetail
  //  [1] -- ONE frag dovetail, ONE frag contained
  //  [2] -- BOTH frag and mate are contained

  int64   nunmated            = 0;

  int64   ngood[3]            = { 0, 0, 0 };

  int64   nbadFwd[3]          = { 0, 0, 0 };
  int64   nbadRev[3]          = { 0, 0, 0 };

  int64   ngoodExternalFwd[3] = { 0, 0, 0 };
  int64   ngoodExternalRev[3] = { 0, 0, 0 };

  int64   nbadExternalFwd[3]  = { 0, 0, 0 };
  int64   nbadExternalRev[3]  = { 0, 0, 0 };

  int64   nbadCompressed[3]   = { 0, 0, 0 };
  int64   nbadStretched[3]    = { 0, 0, 0 };

  int64   nbadNormal[3]       = { 0, 0, 0 };
  int64   nbadAnti[3]         = { 0, 0, 0 };
  int64   nbadOuttie[3]       = { 0, 0, 0 };

  if (IS)
    delete IS;
  IS = new InsertSizes(unitigs);

  for (uint32 ti=0; ti<unitigs.size(); ti++) {
    Unitig  *thisUtg = unitigs[ti];

    if ((thisUtg == NULL) ||
        (thisUtg->ufpath.size() < 2))
      continue;

    MateLocation          positions(unitigs, thisUtg);

    //positions.dumpHappiness(output_prefix, name);

    nunmated += positions.nunmated;

    for (int32 i=0; i<3; i++) {
      ngood[i] += positions.ngood[i];

      nbadFwd[i] += positions.nbadFwd[i];
      nbadRev[i] += positions.nbadRev[i];

      ngoodExternalFwd[i] += positions.ngoodExternalFwd[i];
      ngoodExternalRev[i] += positions.ngoodExternalRev[i];

      nbadExternalFwd[i] += positions.nbadExternalFwd[i];
      nbadExternalRev[i] += positions.nbadExternalRev[i];

      nbadCompressed[i] += positions.nbadCompressed[i];
      nbadStretched[i]  += positions.nbadStretched[i];

      nbadNormal[i] += positions.nbadNormal[i];
      nbadAnti[i]   += positions.nbadAnti[i];
      nbadOuttie[i] += positions.nbadOuttie[i];
    }
  }

  fprintf(logFile, "==> MATE HAPPINESS\n");
  fprintf(logFile, "\n");
  fprintf(logFile, "unmated          %10"F_S64P"\n", nunmated);
  fprintf(logFile, "\n");
  fprintf(logFile, "                  dove-dove  dove-cont  cont-cont\n");
  fprintf(logFile, "happy            %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngood[0],            ngood[1],            ngood[2]);
  fprintf(logFile, "\n");
  fprintf(logFile, "goodExternalFwd  %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngoodExternalFwd[0], ngoodExternalFwd[1], ngoodExternalFwd[2]);
  fprintf(logFile, "goodExternalRev  %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngoodExternalRev[0], ngoodExternalRev[1], ngoodExternalRev[2]);
  fprintf(logFile, "\n");
  fprintf(logFile, "badExternalFwd   %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadExternalFwd[0],  nbadExternalFwd[1],  nbadExternalFwd[2]);
  fprintf(logFile, "badExternalRev   %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadExternalRev[0],  nbadExternalRev[1],  nbadExternalRev[2]);
  fprintf(logFile, "badCompressed    %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadCompressed[0],   nbadCompressed[1],   nbadCompressed[2]);
  fprintf(logFile, "badStretched     %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadStretched[0],    nbadStretched[1],    nbadStretched[2]);
  fprintf(logFile, "badNormal        %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadNormal[0],       nbadNormal[1],       nbadNormal[2]);
  fprintf(logFile, "badAnti          %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadAnti[0],         nbadAnti[1],         nbadAnti[2]);
  fprintf(logFile, "badOuttie        %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadOuttie[0],       nbadOuttie[1],       nbadOuttie[2]);

  fprintf(logFile, "badFwd           %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadFwd[0],          nbadFwd[1],          nbadFwd[2]);
  fprintf(logFile, "badRev           %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadRev[0],          nbadRev[1],          nbadRev[2]);
}

