
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: AS_OBT_acceptableOverlap.c,v 1.1 2008-10-14 03:05:37 brianwalenz Exp $";

#include "AS_OBT_acceptableOverlap.h"
#include "constants.H"

int
AS_OBT_acceptableOverlap(OVSoverlap ol) {
  int32   minq = AS_OVS_encodeQuality(OBT_MIN_ERATE / 100.0);

  //  Remember if the 5' ends are far apart -- but we only care if
  //  it's a forward match.
  //
  int far5prime = TRUE;
  if (ol.dat.obt.fwd) {
    if (ol.dat.obt.a_beg > ol.dat.obt.b_beg)
      far5prime = ((ol.dat.obt.a_beg - ol.dat.obt.b_beg) > OBT_FAR5PRIME);
    else
      far5prime = ((ol.dat.obt.b_beg - ol.dat.obt.a_beg) > OBT_FAR5PRIME);
  }

  int Adiff = ol.dat.obt.a_end - ol.dat.obt.a_beg;
  int Bdiff = ol.dat.obt.b_end - ol.dat.obt.b_beg;
  if (ol.dat.obt.b_end < ol.dat.obt.b_beg)
    Bdiff = ol.dat.obt.b_beg - ol.dat.obt.b_end;

  //  It's an acceptable overlap if the error is within tolerance,
  //  if it's long, and if the 5' ends are far apart.
  //
  return(((ol.dat.obt.erate < minq) ||
          ((Adiff > OBT_MIN_DIFF) &&
           (Bdiff > OBT_MIN_DIFF))) &&
         (far5prime));
}
