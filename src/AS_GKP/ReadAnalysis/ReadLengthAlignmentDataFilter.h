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

#ifndef READLENGTHALIGNMENTDATAFILTER_H
#define READLENGTHALIGNMENTDATAFILTER_H

static const char* rcsid_READLENGTHALIGNMENTDATAFILTER_H = "$Id: ReadLengthAlignmentDataFilter.h,v 1.1 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AlignmentDataFilter.h"
#include "AlignmentDataStats.h"
#include "ReadAlignment.h"
#include "StringUtils.h"

using namespace Utility;

namespace ReadAnalysis
{
	static const size_t DEFAULT_MIN_READ_LENGTH = 1;
	static const size_t DEFAULT_MAX_READ_LENGTH = UINT64_MAX;
	
	class ReadLengthAlignmentDataFilter : public AlignmentDataFilter
	{
	public:
		ReadLengthAlignmentDataFilter(size_t minLength = DEFAULT_MIN_READ_LENGTH, 
			size_t maxLength = DEFAULT_MAX_READ_LENGTH, AlignmentDataStats* dataStats = NULL);
		virtual ~ReadLengthAlignmentDataFilter();
		
		virtual bool filterReadAlign(ReadAlignment* readAlign);
		
		virtual string toString();
		
		size_t getMinLength()
		{
			return this->minLength;
		}
		
		void setMinLength(size_t minLength)
		{
			this->minLength = minLength;
		}
		
		size_t getMaxLength()
		{
			return this->maxLength;
		}
		
		void setMaxLength(size_t maxLength)
		{
			this->maxLength = maxLength;
		}
		
	protected:
		size_t minLength;
		size_t maxLength;
	};
}

#endif
