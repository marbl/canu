
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
   CVS_ID:  $Id: MicroHetInterfaceREZ_test3.h,v 1.3 2005-03-22 19:07:42 jason_miller Exp $
 *********************************************************************/
#ifndef AS_REZ_MicroHetInterfaceREZ_H
#define AS_REZ_MicroHetInterfaceREZ_H

#include "MicroHetREZ_test3.h"
#include "MicroHetScoreREZ_test3.h"
#include "MicroHetPartitionsREZ_test3.h"
#include "UtilsREZ.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStorePartition.h"
#include "AS_PER_distStore.h"


Alignment_t *AS_REZ_convert_array_to_alignment(char **ar, int c, int r);
Alignment_t *AS_REZ_convert_IUM_to_alignment(IntUnitigMesg* ium, FragStoreHandle handle, tFragStorePartition *phandle,int compress);
CGB_Type AS_REZ_get_simulator_type(IntUnitigMesg* ium_mesg);


/* This function compresses shredded fragments from the same location 
   into basically a 1x coverage, such that there are no aritfical microhets;
   it also nulls out all but one position of a multibase gap, to reduce the
   effect of multibase indel polymorphisms on microhet detection. */
void compress_shreds_and_null_indels(int c, int r,  FragStoreHandle frag_store, tFragStorePartition *phandle,
			       char **array, int **id_array, int verbose);

#endif

