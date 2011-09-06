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

static const char* rcsid = "$Id: AlignmentDataFilter.C,v 1.1 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "AlignmentDataFilter.h"

using namespace ReadAnalysis;

AlignmentDataFilter::AlignmentDataFilter(AlignmentDataStats* dataStats)
{
	this->dataStats = dataStats;
}

AlignmentDataFilter::~AlignmentDataFilter()
{
}

bool AlignmentDataFilter::filterReadAlign(ReadAlignment* readAlign)
{
	return false;
}

void AlignmentDataFilter::filterData(vector<ReadAlignment*>& data)
{
}

string AlignmentDataFilter::toString()
{
	return "";
}

void AlignmentDataFilter::recordFilteredRead(ReadAlignment* readAlign)
{
	if (this->dataStats != NULL)
	{
		this->dataStats->addFilteredRead(this, readAlign->getIid());
	}
}
