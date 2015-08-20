
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
 *    src/AS_BAT/AS_BAT_EvaluateMates.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz beginning on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

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

  writeLog("==> MATE HAPPINESS\n");
  writeLog("\n");
  writeLog("unmated          "F_S64"\n", nunmated);
  writeLog("\n");
  writeLog("                  dove-dove  dove-cont  cont-cont\n");
  writeLog("happy            %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngood[0],            ngood[1],            ngood[2]);
  writeLog("\n");
  writeLog("goodExternalFwd  %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngoodExternalFwd[0], ngoodExternalFwd[1], ngoodExternalFwd[2]);
  writeLog("goodExternalRev  %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", ngoodExternalRev[0], ngoodExternalRev[1], ngoodExternalRev[2]);
  writeLog("\n");
  writeLog("badExternalFwd   %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadExternalFwd[0],  nbadExternalFwd[1],  nbadExternalFwd[2]);
  writeLog("badExternalRev   %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadExternalRev[0],  nbadExternalRev[1],  nbadExternalRev[2]);
  writeLog("badCompressed    %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadCompressed[0],   nbadCompressed[1],   nbadCompressed[2]);
  writeLog("badStretched     %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadStretched[0],    nbadStretched[1],    nbadStretched[2]);
  writeLog("badNormal        %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadNormal[0],       nbadNormal[1],       nbadNormal[2]);
  writeLog("badAnti          %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadAnti[0],         nbadAnti[1],         nbadAnti[2]);
  writeLog("badOuttie        %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadOuttie[0],       nbadOuttie[1],       nbadOuttie[2]);

  writeLog("badFwd           %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadFwd[0],          nbadFwd[1],          nbadFwd[2]);
  writeLog("badRev           %10"F_S64P" %10"F_S64P" %10"F_S64P"\n", nbadRev[0],          nbadRev[1],          nbadRev[2]);
}

