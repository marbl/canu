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

static const char* rcsid = "$Id: SnapperAlignmentDataReader.C,v 1.7 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "SnapperAlignmentDataReader.h"

using namespace ReadAnalysis;

SnapperAlignmentDataReader::SnapperAlignmentDataReader()
	: AlignmentDataReader()
{
}

SnapperAlignmentDataReader::~SnapperAlignmentDataReader()
{
}

void SnapperAlignmentDataReader::processStart(string line, size_t lineNum)
{
	if (line != SNAPPER_READ_ALIGNMENT_START)
	{
		string lineStr;
		
		throw DataException(line, string("Read alignment start expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
}

void SnapperAlignmentDataReader::processReadInfo(ReadAlignment* readAlign, string line, size_t lineNum)
{
	size_t index, length;
	
	if (!sscanf(line.c_str(), SNAPPER_READ_INFO_FORMAT, &index, &length))
	{
		string lineStr;
		
		throw DataException(line, string("Read info expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
	
	readAlign->setIndex(index);
	readAlign->setLength(length);
}

void SnapperAlignmentDataReader::processReadDefinition(ReadAlignment* readAlign, string line, size_t lineNum)
{
	AS_IID iid, mateIid = 0;
	
	if (!sscanf(line.c_str(), SNAPPER_SANGER_READ_DEFINITION_FORMAT, &iid, &mateIid))
	{
		if (!sscanf(line.c_str(), SNAPPER_ILLUMINA_READ_DEFINITION_FORMAT, &iid))
		{
			string lineStr;
			
			throw DataException(line, string("Read definition expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
		}
	}
	
	readAlign->setIid(iid);
	readAlign->setMateIid(mateIid);
}

void SnapperAlignmentDataReader::processGenomeDefinition(ReadAlignment* readAlign, string line, size_t lineNum)
{
	// TODO: process, if there is anything useful
	/*
	if (!sscanf(line.c_str(), SNAPPER_GENOME_DEFINITION_FORMAT))
	{
		string lineStr;
		
		throw DataException(line, string("Genome definition expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
	*/
}

void SnapperAlignmentDataReader::processAlignmentInfo(ReadAlignment* readAlign, string line, size_t lineNum)
{
	pair<size_t, size_t> alignedSegment, alignment;
	uint16 identity;
	
	if (!sscanf(line.c_str(), SNAPPER_ALIGNMENT_INFO_FORMAT, &alignedSegment.first, &alignedSegment.second, 
		&alignment.first, &alignment.second, &identity))
	{
		string lineStr;
		
		throw DataException(line, string("Alignment info expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
	
	readAlign->setAlignedSegment(alignedSegment);
	readAlign->setAlignment(alignment);
	readAlign->setIdentity(identity);
}

void SnapperAlignmentDataReader::processReadSequence(ReadAlignment* readAlign, string line, size_t lineNum)
{
	if (StringUtils::isBlank(line))
	{
		string lineStr;
		
		throw DataException(line, string("Read sequence expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
	
	processSequence(readAlign->getSequence(), line, lineNum);
}

void SnapperAlignmentDataReader::processGenomeSequence(ReadAlignment* readAlign, string line, size_t lineNum)
{
	if (StringUtils::isBlank(line))
	{
		string lineStr;
		
		throw DataException(line, string("Genome sequence expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
	
	processSequence(readAlign->getGenomeSequence(), line, lineNum);
}

void SnapperAlignmentDataReader::processSequence(vector<char>& sequence, string line, size_t lineNum)
{
	copy(line.begin(), line.end(), sequence.begin());
}

void SnapperAlignmentDataReader::processEnd(string line, size_t lineNum)
{
	if (line != SNAPPER_READ_ALIGNMENT_END)
	{
		string lineStr;
		
		throw DataException(line, string("Read alignment end expected on line ") + StringUtils::toString(lineNum, lineStr) + ".");
	}
}

void SnapperAlignmentDataReader::processData()
{
	clock_t startClock = clock();
	
	fprintf(stderr, "Processing Snapper read alignment data:\n");
	
	size_t lineNum = 0;
	string line;
	
	ReadAlignment* readAlign;
	
	while (FileUtils::canRead(this->stream))
	{
		// Read alignment start
		line = readLine(this->stream, line, lineNum);
		
		if (line.empty())
		{
			break;
		}
		
		processStart(line, lineNum);
		
		// Read info
		line = readLine(this->stream, line, lineNum);
		readAlign = new ReadAlignment();
		processReadInfo(readAlign, line, lineNum);
		
		// Read definition
		line = readLine(this->stream, line, lineNum);
		processReadDefinition(readAlign, line, lineNum);
		
		// Genome definition
		line = readLine(this->stream, line, lineNum);
		processGenomeDefinition(readAlign, line, lineNum);
		
		// Alignment info
		line = readLine(this->stream, line, lineNum);
		processAlignmentInfo(readAlign, line, lineNum);
		
		// Read sequence
		line = readLine(this->stream, line, lineNum);
		processReadSequence(readAlign, line, lineNum);
		
		// Genome sequence
		line = readLine(this->stream, line, lineNum);
		processGenomeSequence(readAlign, line, lineNum);
		
		// Read alignment end
		line = readLine(this->stream, line, lineNum);
		processEnd(line, lineNum);
		
		if (!this->filterReadAlign(readAlign))
		{
			this->dataStats.addRead(readAlign);
			this->data.push_back(readAlign);
			this->iidMap[readAlign->getIid()] = readAlign;
		}
	}
	
	string timeStr;
	
	fprintf(stderr, "Processed "F_SIZE_T" Snapper read alignment[s] in: "F_STR"\n", this->data.size(), 
		StringUtils::toString(timeStr, startClock, true).c_str());
}
