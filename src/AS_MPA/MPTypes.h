
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
/* $Id: MPTypes.h,v 1.5 2005-09-21 20:13:07 catmandew Exp $ */
#ifndef MPTYPES_H
#define MPTYPES_H 1

#include "cds.h"

#define ID_TYPE CDS_UID_t
#define F_MPID  F_UID
#define BOGUS_ID CDS_UINT64_MAX

#define UNIT_TYPE int64
#define F_MPUNIT  F_S64

#define COINCIDENT_THRESHOLD  50  // bp
#define CONFIRMATION_THRESHOLD 2  // 2 agreeing unsatisfieds confirm each other
#define STDDEVS_THRESHOLD 3.0     // mean +/- N stddevs

typedef enum
{
  ISS_ALL,
  ISS_HIGHER,
  ISS_LOWER,
  ISS_NUM_SCOPES
} InterSequenceScope_e;
#define INTER_S_DEFAULT_SCOPE  ISS_ALL

typedef enum
{
  CR_NATIVE,
  CR_COMPATIBLE
} CompressedRepresentation_e;

#endif //  MPTYPES_H
