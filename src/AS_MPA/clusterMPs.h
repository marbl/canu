
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
/* $Id: clusterMPs.h,v 1.1.1.1 2004-04-14 13:51:58 catmandew Exp $ */
#ifndef CLUSTERMPS_H
#define CLUSTERMPS_H

#include <fstream>

#include <map>
#include <vector>

#include "CloneLibrary.h"
#include "MatePair.h"
#include "CompositeMPPolygon.h"

void ReadCloneLibs(map<ID_TYPE, CloneLibrary> & libs, ifstream & fin);
void ReadMatePairs(vector<MatePair> & mps, ifstream & fin);
void ClusterMPPs(const vector<CompositeMPPolygon<UNIT_TYPE> > & mpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & pmpps,
                 vector<CompositeMPPolygon<UNIT_TYPE> > & fmpps,
                 MatePairIndex_e mpii, unsigned int filterThresh);
void RefineInversions(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                      const vector<MatePair> & smpsv,
                      unsigned int filterThresh);
void RefineStretched(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                     const vector<MatePair> & smpsv);
void DetectTranspositions(const vector<CompositeMPPolygon<UNIT_TYPE> > & compressed,
                          const vector<CompositeMPPolygon<UNIT_TYPE> > & stretched,
                          const vector<CompositeMPPolygon<UNIT_TYPE> > & outties,
                          vector<CompositeMPPolygon<UNIT_TYPE> > & trans,
                          map<ID_TYPE, CloneLibrary> & libs,
                          double numStddevs,
                          unsigned int filterThresh);
void FilterMislabelledLibs(vector<CompositeMPPolygon<UNIT_TYPE> > & cmpps,
                           vector<CompositeMPPolygon<UNIT_TYPE> > & mmpps,
                           unsigned int minFilterCount);


#endif // #ifndef CLUSTERMPS_H
