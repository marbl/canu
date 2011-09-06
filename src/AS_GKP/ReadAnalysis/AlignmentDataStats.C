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

static const char* rcsid = "$Id: AlignmentDataStats.C,v 1.3 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "AlignmentDataStats.h"

using namespace ReadAnalysis;

AlignmentDataStats::AlignmentDataStats()
{
	this->numReads = 0;
	this->minReadLength = 0;
	this->maxReadLength = 0;
	this->meanReadLength = 0;
}

size_t AlignmentDataStats::getNumFilteredReads(AlignmentDataFilter* filter)
{
	if (filter != NULL)
	{
		return (this->filteredReadMap.count(filter) > 0) ? this->filteredReadMap[filter]->size() : 0;
	}
	else
	{
		size_t numFilteredReads = 0;
		
		for (map<AlignmentDataFilter*, vector<AS_IID>*>::iterator iterator = this->filteredReadMap.begin(); 
			iterator != this->filteredReadMap.end(); iterator++)
		{
			numFilteredReads += (*iterator).second->size();
		}
		
		return numFilteredReads;
	}
}

void AlignmentDataStats::addFilteredRead(AlignmentDataFilter* filter, AS_IID iid)
{
	vector<AS_IID>* filteredReads;
	
	filteredReads = (this->filteredReadMap.count(filter) > 0) ? this->filteredReadMap[filter] : NULL;
	
	if (filteredReads == NULL)
	{
		filteredReads = new vector<AS_IID>();
		
		this->filteredReadMap[filter] = filteredReads;
	}
	
	filteredReads->push_back(iid);
}

size_t AlignmentDataStats::getNumReads(size_t baseIndex)
{
	vector<AS_IID> reads;
	
	return this->getReads(baseIndex, reads).size();
}

vector<AS_IID>& AlignmentDataStats::getReads(size_t baseIndex, vector<AS_IID>& buffer)
{
	vector<AS_IID>* reads;
	
	for (map<size_t, vector<AS_IID>*>::iterator iterator = this->readLengthMap.begin(); iterator != this->readLengthMap.end(); iterator++)
	{
		if ((*iterator).first >= baseIndex)
		{
			reads = (*iterator).second;
			
			buffer.insert(buffer.end(), reads->begin(), reads->end());
		}
	}
	
	return buffer;
}

void AlignmentDataStats::addRead(ReadAlignment* readAlign)
{
	this->numReads++;
	
	this->minReadLength = (this->numReads > 1) ? min(this->minReadLength, readAlign->getLength()) : readAlign->getLength();
	this->maxReadLength = (this->numReads > 1) ? max(this->maxReadLength, readAlign->getLength()) : readAlign->getLength();
	this->meanReadLength = (this->numReads > 1) ? (double)(this->meanReadLength + readAlign->getLength()) / 2 : (double)readAlign->getLength();
	
	size_t readLengthIndex = readAlign->getLength() - 1;
	vector<AS_IID>* readLengthBucket = this->readLengthMap.count(readLengthIndex) ? this->readLengthMap[readLengthIndex] : NULL;
	
	if (readLengthBucket == NULL)
	{
		readLengthBucket = new vector<AS_IID>();
		
		this->readLengthMap[readLengthIndex] = readLengthBucket;
	}
	
	readLengthBucket->push_back(readAlign->getIid());
}
