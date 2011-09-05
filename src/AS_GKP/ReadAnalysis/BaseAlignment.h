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

#ifndef BASEALIGNMENT_H
#define BASEALIGNMENT_H

static const char* rcsid_BASEALIGNMENT_H = "$Id: BaseAlignment.h,v 1.2 2011-09-05 16:49:44 mkotelbajcvi Exp $";

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

using namespace std;

#include "AlignmentError.h"
#include "AlignmentErrorType.h"

namespace ReadAnalysis
{
	class BaseAlignment
	{
	public:
		BaseAlignment(size_t position, size_t readsReserveSize, size_t errorTypeReserveSize);
		
		void addRead(AS_IID iid);
		void addError(AlignmentError* error);
		
		size_t getPosition()
		{
			return this->position;
		}
		
		void setPosition(size_t position)
		{
			this->position = position;
		}
		
		vector<AS_IID>& getReads()
		{
			return this->reads;
		}
		
		map<AlignmentErrorType, vector<AlignmentError*>* >& getErrors()
		{
			return this->errors;
		}
		
	protected:
		size_t position;
		vector<AS_IID> reads;
		map<AlignmentErrorType, vector<AlignmentError*>* > errors;
		
		void initErrorType(AlignmentErrorType type, size_t reserveSize);
	};
}

#endif
