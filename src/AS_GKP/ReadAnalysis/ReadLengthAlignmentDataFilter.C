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

static const char* rcsid = "$Id$";

#include "ReadLengthAlignmentDataFilter.h"

using namespace ReadAnalysis;

ReadLengthAlignmentDataFilter::ReadLengthAlignmentDataFilter(size_t minLength, size_t maxLength, AlignmentDataStats* dataStats)
	: AlignmentDataFilter(dataStats)
{
	this->minLength = minLength;
	this->maxLength = maxLength;
}

ReadLengthAlignmentDataFilter::~ReadLengthAlignmentDataFilter()
{
}

bool ReadLengthAlignmentDataFilter::filterReadAlign(ReadAlignment* readAlign)
{
	if ((readAlign->getLength() < this->minLength) || (readAlign->getLength() > this->maxLength))
	{
		this->recordFilteredRead(readAlign);
	
		return true;
	}
	else
	{
		return false;
	}
}

string ReadLengthAlignmentDataFilter::toString()
{
	string minLengthStr, maxLengthStr;
	
	return "ReadLengthAlignmentDataFilter (minLength=" + StringUtils::toString(this->minLength, minLengthStr) + 
		", maxLength=" + StringUtils::toString(this->maxLength, maxLengthStr) + ")";
}
