
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
/*********************************************************************
   CVS_ID:  $Id: MicroHetInterfaceREZ.h,v 1.2 2004-09-23 20:25:28 mcschatz Exp $
 *********************************************************************/
#ifndef MicroHetInterfaceREZ_H
#define MicroHetInterfaceREZ_H

#include "MicroHetREZ.h"
#include "MicroHetScoreREZ.h"
#include "MicroHetPartitionsREZ.h"
#include "UtilsREZ.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"


Alignment_t *convert_array_to_alignment(char **ar, int c, int r);
Alignment_t *convert_IUM_to_alignment(IntUnitigMesg* ium, FragStoreHandle handle);
CGB_Type get_simulator_type(IntUnitigMesg* ium_mesg);

#endif








