
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
/* 	$Id: AS_global.h,v 1.1.1.1 2004-04-14 13:47:07 catmandew Exp $	 */

/* This is the global include file that all C files in the AS subsystem should
   include.
*/
#ifndef AS_GLOBAL_H
#define AS_GLOBAL_H

#include "cds.h"
#include "AS_MSG_pmesg.h"
#include "PrimitiveVA.h"

// Constants that SHOULD be included

#define CGB_INVALID_CUTOFF           -12.0f
// A threshold value for Gene^s coverage statistic. Values BELOW this value
// have never been known to associated with unitigs with fragments that ARE
// not contiguous in the genome. They are guaranteed REPEATS.
#define CGB_UNIQUE_CUTOFF            10.0f
//#define CGB_UNIQUE_CUTOFF            12.0f
// A threshold value for Gene^s coverage statistic. Values above this value
// have never been known to associated with unitigs with fragments that are
// not contiguous in the genome.

#define CGB_TANDEM_REPEAT_THRESHOLD  50  
// A threshold distance in base pairs
// for the allowed slop between the minimum overlap and maximum overlap before
// calling an overlap definately a tandem repeat.

#define  AS_READ_ERROR_RATE         0.06
//#define  AS_READ_ERROR_RATE         0.10
    //  Errors per base allowed in matching regions between frag reads
#define  AS_GUIDE_ERROR_RATE        0.06
//#define  AS_GUIDE_ERROR_RATE        0.10
    //  Errors per base allowed in matching regions involving BAC ends
    //  or other guides.

#define AS_CGB_BPT_MIN_PREFIX   25
#define AS_CGB_BPT_MIN_SUFFIX   25
#define AS_CGB_BPT_ERROR_RATE    0.08
#define AS_CGB_BPT_PROB_THOLD    1e-6
// common parameters for calling BPnt_Seq_Comp_AS in AS_CGB & AS_URT

#ifndef max
# define max(a,b)		( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
# define min(a,b)		( ((a) < (b)) ? (a) : (b) )
#endif

// cgw and cns use NULLINDEX for a NULL index value
#define NULLINDEX (-1)

// A convenient assert for testing whether ptrs are null
// without bothering lint
#define AssertPtr(ptr) (assert((ptr) != NULL))


#endif
