
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
 * $Id: AS_CGB_miniunitigger.h,v 1.1.1.1 2004-04-14 13:49:42 catmandew Exp $
 *
 * Module: AS_CGB_miniunitigger.h
 * Description: A subroutine (object?) interface for the unitigger.
 * Assumptions: Too many to count.
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_MINIUNITIGGER_INCLUDE
#define AS_CGB_MINIUNITIGGER_INCLUDE

#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"

#if 0
VA_DEF(OFGMesg);
VA_DEF(OverlapMesg);
VA_DEF(IntMultiPos);
VA_DEF(IntUnitigMesg);
VA_DEF(char);
#endif

typedef void MiniUnitiggerObject;
// Hide the implementation of the Unitigger object.  I encapsulate the
// mini Unitigger (1) data with a void pointer and (2) methods by
// extern functions.

extern MiniUnitiggerObject * createMiniUnitigger( 
  int maxfrags,
  int maxedges,
  int maxtext
  );
// Create an "unitigger" object. The return value is non-NULL when
// successful.  The maxfrags is a pre-allocation amount for the
// maximum number of fragments.  The maxedges is a pre-allocation
// amount for the maximum number of directed edges in the graph. The
// maxtext is pre-allocation for the maximum amount of characters in
// the source test of the OFG messages.  You can specify zero for the
// pre-allocation amounts.

extern int set_cgb_unique_cutoff_MiniUnitigger
( MiniUnitiggerObject * self, float cgb_unique_cutoff);
// A typical value used for the Human genome is five. A return value
// of zero indicates success. 

extern int set_overlap_error_threshold_MiniUnitigger
( MiniUnitiggerObject * self, float overlap_error_threshold);
// An overlap_error_threshold of 0.06 means six percent.  This value is
// used to filter the overlaps on input.  To turn it off, set the
// value to one. A return value of zero indicates success.

extern int set_as_cgb_max_frag_iid_MiniUnitigger
( MiniUnitiggerObject * self, int as_cgb_max_frag_iid);
// Unfortunately I do not have a hash function that maps the fragment
// IIDs in the OFG messages to a dense interval of integers. Instead I
// just allocate an array of with range [0 .. as_cgb_max_frag_iid] to
// store the mapping. A return value of zero indicates success.

#if 0 // Dynamic binding
// If later we want to dynamicly bind a function inside the Mini
// Unitigger, then I propose to use a binding method similar to that
// used by the AS_CGB histogramming utility.  For example:
extern int extend_MiniUnitigger_with_YourDataType
(  
 MiniUnitiggerObject *self,
 size_t sizeof_YourDataType,
 YourDataType * (*indexdata)(YourDataType *b,int ib),
 void (*setdata)(YourDataType *a,int ib,YourDataType *b),
 void (*aggregate)(YourDataType *a,int ib,YourDataType *b),
 void (*printdata)(FILE *fout,
                   YourDataType *,
                   YourDataType *,
                   YourDataType *)
 );
// A return value of zero indicates success.
#endif  // Dynamic binding

typedef struct
{
  VA_TYPE(OFGMesg)     * the_ofg_messages;
  VA_TYPE(OverlapMesg) * the_ovl_messages;
  VA_TYPE(char)        * the_ofg_source;
  VA_TYPE(char)        * the_ovl_source;
} RunMiniUnitiggerParams;

typedef struct
{
  VA_TYPE(IntMultiPos)   * the_imp_messages;
  VA_TYPE(IntUnitigMesg) * the_ium_messages;
  VA_TYPE(char)          * the_imp_source;
  VA_TYPE(char)          * the_ium_source;
} RunMiniUnitiggerResults;

extern int run_MiniUnitigger
(
 MiniUnitiggerObject     * self,
 RunMiniUnitiggerParams  * params,
 RunMiniUnitiggerResults * results
 );
// We assume that the client owns the VAs in the params and results
// structures and is responsible for creating and destroying the VA
// arrays.

extern int destroyMiniUnitigger( MiniUnitiggerObject * self );
// Destroy an "unitigger" object. The return value is zero when successful.

#endif // AS_CGB_MINIUNITIGGER_INCLUDE
