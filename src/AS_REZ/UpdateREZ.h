
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
/**********************************************************************

        Module:  UpdateREZ.h

   Description:  Declaration of common data types used in UpdateREZ.c

    Programmer:  S. Lonardi (stelo@cs.purdue.edu)

       Written:  17 May 99
 **********************************************************************/

/*********************************************************************
   CVS_ID: $Id: UpdateREZ.h,v 1.2 2004-09-23 20:25:28 mcschatz Exp $
 *********************************************************************/

#ifndef UPDATE_REZ_H
#define UPDATE_REZ_H

#define DEBUG_UPDATE              0 // 0 no debug, 1 some, 2 more, 3 for dump scaffold contents

#define UPDATE_PARANOID_CHECKING  0 // set to 1 and it will check internal connectivity & trusted edges

#define COMPUTE_VARIANCE_ESTIMATE 0 // set this to 1 to assign an estimate to the variance to the
                                    // chunks we are inserting
//
// Trust_Edges() marks all the edges that the CI <cid> has to uniques
// and others CI in the fill_chunks as tentatively trusted: returns
// the number of RAW edge marked
//
int Trust_Edges(ScaffoldGraphT *,
		Scaffold_Fill_t *,
		Gap_Chunk_t * *,
		int32,
		int32);

//
// Update_Scaffold_Graph() updates the scaffolds by filling the
// appropriate gaps as indicated by the infos in the gapAssignment
// (only the chunks which has a "keep" flag TRUE will be
// considered). If the flag ChiSquare is TRUE then we will compute the
// ChiSquare test to re-mark internal edges
// The last boolean value indicates whether before insertinfg a chunk
// it should be split and the resulting surrogate gets inserted.
// Kind_Of_Fill_t  indicates whether the fill items are rocks, stones or walks
//
int Update_Scaffold_Graph(ScaffoldGraphT *,
			  Scaffold_Fill_t *,
			  int,
			  int,
			  int,
			  int copyAllOverlaps,
			  int,
                          Kind_Of_Fill_t);

#endif
