
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
 $Id: PublicAPI_CNS.h,v 1.4 2005-03-22 19:48:51 jason_miller Exp $
 *********************************************************************/

#ifndef PUBLICAPI_CNS_INCLUDE
#define PUBLICAPI_CNS_INCLUDE

#include "Globals_CNS.h"

void CNS_setExitStatus (int exitStatus) ;

void CNS_setConsensusParametersToDefault() ;

void CNS_setConsensusParametersIndividually(
    int useSDB,
    int usePartSDB,
    float priorSequencingErrorRate,
    float priorSNPRate,
    int numHaplotypes,
    int basecallingUsesPublicData,
    int basecallingFavorsPublicData,
    int outputDebugging,
    int isPartitioned) ;

int MultiAlignUnitig(IntUnitigMesg *iunitig,
                         FragStoreHandle frgStore,
			 VA_TYPE(char) *sequence,
			 VA_TYPE(char) *quality, 
			 VA_TYPE(int32) *deltas, 
			 CNS_PrintKey printwhat, 
			 int do_rez, 
			 Overlap *(*)(COMPARE_ARGS));
int MultiAlignContig(IntConConMesg *contig, 
			 VA_TYPE(char) *sequence,
			 VA_TYPE(char) *quality, 
			 VA_TYPE(int32) *deltas, 
			 CNS_PrintKey printwhat, 
			 Overlap *(*)(COMPARE_ARGS));

#endif // PUBLICAPI_CNS_INCLUDE


