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

#ifndef READALIGNMENTPROFILER_H
#define READALIGNMENTPROFILER_H

static const char* rcsid_READALIGNMENTPROFILER_H = "$Id: ReadAlignmentProfiler.h,v 1.3 2011-09-03 01:29:50 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "AlignmentError.h"
#include "AlignmentErrorType.h"
#include "AS_global.h"
#include "BaseAlignment.h"
#include "ErrorUtils.h"
#include "FileUtils.h"
#include "ReadAlignment.h"
#include "StringUtils.h"

using namespace Utility;

namespace ReadAnalysis
{
	static const size_t DEFAULT_BASES_RESERVE_SIZE = 10240;
	static const size_t DEFAULT_BASES_RESIZE_FACTOR = 2;
	
	class ReadAlignmentProfiler
	{
	public:
		ReadAlignmentProfiler();
		
		inline static bool isBaseCall(const char base)
		{
			return (base == 'a') || (base == 'A') || 
				(base == 'c') || (base == 'C') || 
				(base == 'g') || (base == 'G') || 
				(base == 't') || (base == 'T');
		}
		
		inline static bool isBaseError(const char base)
		{
			return (base == 'A') || (base == 'C') || 
				(base == 'G') || (base == 'T');
		}
		
		inline static bool isBaseGap(const char base)
		{
			return base == '-';;
		}
		
		inline static bool isBaseUnknown(const char base)
		{
			return base == 'N';;
		}
		
		void profileData(vector<ReadAlignment*>& data);
		
		BaseAlignment* getBaseAlign(size_t index);
		
	protected:
		vector<ReadAlignment*> data;
		vector<BaseAlignment*> bases;
	};
}

#endif
