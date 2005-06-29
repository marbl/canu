
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
/**************************************************************
 * AS_CNS/PublicAPI_CNS.h
 *
 * 'Public' functions that other subsystems are invited to link to.
 * CNS is the Consensus sybsystem of the Celer WGS assembler.
 *
 **************************************************************/
/*********************************************************************
 $Id: PublicAPI_CNS.h,v 1.5 2005-06-29 15:21:24 gdenisov Exp $
 *********************************************************************/

#ifndef PUBLICAPI_CNS_INCLUDE
#define PUBLICAPI_CNS_INCLUDE

#include "Globals_CNS.h"

void CNS_setExitStatus (int) ;

void CNS_setConsensusParametersToDefault() ;

void CNS_setConsensusParametersIndividually( int , int , float , float ,
    int , int , int , int , int ) ;

int MultiAlignUnitig(IntUnitigMesg *, FragStoreHandle , VA_TYPE(char) *,
    VA_TYPE(char) *, VA_TYPE(int32) *, CNS_PrintKey , int, 
    Overlap *(*)(COMPARE_ARGS), CNS_Options);

int MultiAlignContig(IntConConMesg *, VA_TYPE(char) *, VA_TYPE(char) *, 
    VA_TYPE(int32) *, CNS_PrintKey , Overlap *(*)(COMPARE_ARGS), CNS_Options);

#endif // PUBLICAPI_CNS_INCLUDE


