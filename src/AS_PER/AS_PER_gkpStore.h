
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

#ifndef AS_PER_GKPFRGSTORE_H
#define AS_PER_GKPFRGSTORE_H

static const char *rcsid_AS_PER_GKPFRGSTORE_H = "$Id: AS_PER_gkpStore.h,v 1.62 2012-05-08 23:17:55 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_fileIO.h"

#define AS_IID_UNK     0
#define AS_IID_BAT     1
#define AS_IID_FRG     2
#define AS_IID_LIB     3

#define UID_NAMESPACE_AS 'U'


uint32 gkStore_decodeClearRegionLabel(const char *label);

class gkLibrary;
class gkFragment;
class gkStore;
class gkStream;
class gkClearRange;
class gkPlacement;

#include "gkStore.H"
#include "gkStream.H"

#include "gkClearRange.H"

#include "gkLibrary.H"
#include "gkFragment.H"
#include "gkPlacement.H"

#endif
