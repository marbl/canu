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

static const char* rcsid = "$Id: ReadAlignmentProfiler.C,v 1.3 2011-09-03 01:29:50 mkotelbajcvi Exp $";

#include "ReadAlignmentProfiler.h"

using namespace ReadAnalysis;

ReadAlignmentProfiler::ReadAlignmentProfiler()
{
	this->bases.reserve(DEFAULT_BASES_RESERVE_SIZE);
}

void ReadAlignmentProfiler::profileData(vector<ReadAlignment*>& data)
{
	this->data = data;
	this->bases.clear();
	
	fprintf(stderr, "Generating profile using "F_U64" read alignment[s] ...\n", this->data.size());
	
	ReadAlignment* readAlign;
	BaseAlignment* baseAlign;
	char readBase, genomeBase;
	
	for (size_t a = 0; a < this->data.size(); a++)
	{
		readAlign = this->data[a];
		
		for (size_t b = 0; b < readAlign->getLength(); b++)
		{
			readBase = readAlign->getSequence()[b];
			genomeBase = readAlign->getGenomeSequence()[b];
			
			baseAlign = getBaseAlign(b);
			
			if (isBaseCall(readBase))
			{
				if (isBaseError(readBase))
				{
					if (isBaseError(genomeBase))
					{
						baseAlign->addError(AlignmentError(readAlign->getIid(), MISMATCH));
					}
					else if (isBaseGap(genomeBase))
					{
						baseAlign->addError(AlignmentError(readAlign->getIid(), INSERTION));
					}
				}
			}
			else if (isBaseGap(readBase))
			{
				baseAlign->addError(AlignmentError(readAlign->getIid(), DELETION));
			}
			else if (isBaseUnknown(readBase))
			{
				baseAlign->addError(AlignmentError(readAlign->getIid(), MISMATCH));
			}
			
			baseAlign->getReads().push_back(readAlign->getIid());
		}
	}
	
	fprintf(stderr, "Generated profile for "F_U64" base[s].\n", this->bases.size());
}

BaseAlignment* ReadAlignmentProfiler::getBaseAlign(size_t index)
{
	BaseAlignment* base = this->bases[index];
	
	if (base == NULL)
	{
		base = new BaseAlignment(index);
		
		if (index >= this->bases.size())
		{
			this->bases.reserve(index * DEFAULT_BASES_RESIZE_FACTOR);
			this->bases.resize(index);
		}
		
		this->bases[index] = base;
	}
	
	return base;
}
