/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005, J. Craig Venter Institute. All rights reserved.
 * Author: Brian Walenz
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

#ifndef ALIGNMENTDATAREADER_H
#define ALIGNMENTDATAREADER_H

static const char* rcsid_ALIGNMENTDATAREADER_H = "$Id$";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AlignmentDataFilter.h"
#include "AlignmentDataStats.h"
#include "FileUtils.h"
#include "IOException.h"
#include "ReadAlignment.h"

using namespace Utility;

namespace ReadAnalysis
{
	static const size_t DEFAULT_READ_ALIGNMENT_RESERVE_SIZE = 10240;
	
	class AlignmentDataReader
	{
	public:
		virtual void readData(string path);
		virtual void readData(FILE* stream);
		
		bool getVerbose()
		{
			return this->verbose;
		}
		
		void setVerbose(bool verbose)
		{
			this->verbose = verbose;
		}
		
		vector<AlignmentDataFilter*>& getFilters()
		{
			return this->filters;
		}
		
		vector<ReadAlignment*>& getData()
		{
			return this->data;
		}
		
		map<AS_IID, ReadAlignment*>& getIidMap()
		{
			return this->iidMap;
		}
		
		AlignmentDataStats& getDataStats()
		{
			return this->dataStats;
		}
		
	protected:
		bool verbose;
		vector<AlignmentDataFilter*> filters;
		FILE* stream;
		vector<ReadAlignment*> data;
		map<AS_IID, ReadAlignment*> iidMap;
		AlignmentDataStats dataStats;
		
		AlignmentDataReader();
		virtual ~AlignmentDataReader();
		
		virtual void processData();
		virtual void filterData();
		virtual bool filterReadAlign(ReadAlignment* readAlign);
	};
}

#endif
