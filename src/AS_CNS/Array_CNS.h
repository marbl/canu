
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
   CVS_ID:  $Id: Array_CNS.h,v 1.4 2005-03-22 19:48:40 jason_miller Exp $
 *********************************************************************/
#ifndef AS_CNS_ARRAY_INCLUDE
#define AS_CNS_ARRAY_INCLUDE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStorePartition.h"
#include "AS_UTL_Var.h"
#define ESTDEPTH 32

int IMP2Array(IntMultiPos *frags, 
	      int num_frags, 
	      int length, 
	      FragStoreHandle frag_store, 
              tFragStorePartition *pfrag_store,
              FragStoreHandle bactig_store, 
	      int *depth, 
	      char ***multia, 
	      int ***id_array, 
	      int ***ori_array,
	      int show_cel_status,uint32 clrrng_flag);

int MultiAlignT2Array(MultiAlignT *ma, 
		      FragStoreHandle frag_store, 
		      tFragStorePartition *pfrag_store, 
		      FragStoreHandle bactig_store, 
		      int *depth, 
		      char ***multia, 
		      int ***id_array,
		      int ***ori_array,uint32 clrrng_flag);
#endif
