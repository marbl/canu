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

static const char* rcsid_ALIGNMENTDATAREADER_H = "$Id: AlignmentDataReader.h,v 1.5 2011-09-05 21:23:26 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "AlignmentDataStats.h"
#include "AS_global.h"
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
		
		vector<ReadAlignment*>& getData()
		{
			return this->data;
		}
		
		AlignmentDataStats& getDataStats()
		{
			return this->dataStats;
		}
		
	protected:
		bool verbose;
		FILE* stream;
		vector<ReadAlignment*> data;
		AlignmentDataStats dataStats;
		
		AlignmentDataReader();
		virtual ~AlignmentDataReader();
		
		virtual void processData() = 0;
	};
}

#endif
