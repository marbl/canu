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

#ifndef READALIGNMENT_H
#define READALIGNMENT_H

static const char* rcsid_READALIGNMENT_H = "$Id: ReadAlignment.h,v 1.6 2011-09-06 01:11:56 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AS_UTL_IID.h"

namespace ReadAnalysis
{
	static const size_t DEFAULT_READ_ALIGNMENT_SEQUENCE_RESERVE_SIZE = 512;
	
	class ReadAlignment
	{
	public:
		ReadAlignment(size_t index = 0);
		
		size_t getIndex()
		{
			return this->index;
		}
		
		void setIndex(size_t index)
		{
			this->index = index;
		}
		
		size_t getLength()
		{
			return this->length;
		}
		
		void setLength(size_t length)
		{
			this->length = length;
		}
		
		AS_IID getIid()
		{
			return this->iid;
		}
		
		void setIid(AS_IID iid)
		{
			this->iid = iid;
		}
		
		AS_IID getMateIid()
		{
			return this->mateIid;
		}
		
		void setMateIid(AS_IID mateIid)
		{
			this->mateIid = mateIid;
		}
		
		pair<size_t, size_t> getAlignedSegment()
		{
			return this->alignedSegment;
		}
		
		void setAlignedSegment(pair<size_t, size_t> alignedSegment)
		{
			this->alignedSegment = alignedSegment;
		}
		
		pair<size_t, size_t> getAlignment()
		{
			return this->alignment;
		}
		
		void setAlignment(pair<size_t, size_t> alignment)
		{
			this->alignment = alignment;
		}
		
		uint16 getIdentity()
		{
			return this->identity;
		}
		
		void setIdentity(uint16 identity)
		{
			this->identity = identity;
		}
		
		vector<char>& getSequence()
		{
			return this->sequence;
		}
		
		vector<char>& getGenomeSequence()
		{
			return this->sequence;
		}
		
	protected:
		size_t index;
		size_t length;
		AS_IID iid;
		AS_IID mateIid;
		pair<size_t, size_t> alignedSegment;
		pair<size_t, size_t> alignment;
		uint16 identity;
		vector<char> sequence;
		vector<char> genomeSequence;
	};
}

#endif
