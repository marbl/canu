
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

#ifndef AS_OBT_ACCEPTABLEOVERLAP_H
#define AS_OBT_ACCEPTABLEOVERLAP_H

static const char *rcsid_AS_OBT_ACCEPTABLEOVERLAP_H = "$Id: AS_OBT_acceptableOverlap.h,v 1.1 2008-10-14 03:05:37 brianwalenz Exp $";

#include "AS_OVS_overlap.h"

int    AS_OBT_acceptableOverlap(OVSoverlap ol);

#endif
