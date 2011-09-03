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

static const char* rcsid = "$Id: SnapperAlignmentDataReader.C,v 1.4 2011-09-03 04:08:27 mkotelbajcvi Exp $";

#include "SnapperAlignmentDataReader.h"

using namespace ReadAnalysis;

SnapperAlignmentDataReader::SnapperAlignmentDataReader()
	: AlignmentDataReader()
{
}

SnapperAlignmentDataReader::~SnapperAlignmentDataReader()
{
}

void SnapperAlignmentDataReader::processData()
{
	fprintf(stderr, "Processing Snapper read alignment data ...\n");
	
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
		
		if (line != SNAPPER_READ_ALIGNMENT_START)
		{
			throw DataException(line, string("Read alignment start expected."));
		}
		
		// Read info
		line = readLine(this->stream, line, lineNum);
		
		readAlign = new ReadAlignment();
		
		if (!sscanf(line.c_str(), SNAPPER_READ_INFO_FORMAT, &readAlign->getIndex(), &readAlign->getLength()))
		{
			throw DataException(line, string("Read info expected."));
		}
		
		// Read definition
		line = readLine(this->stream, line, lineNum);
		
		if (!sscanf(line.c_str(), SNAPPER_READ_DEFINITION_FORMAT, &readAlign->getIid(), &readAlign->getMateIid()))
		{
			throw DataException(line, string("Read definition expected."));
		}
		
		// Genome definition
		line = readLine(this->stream, line, lineNum);
		
		/*
		if (!sscanf(line.c_str(), SNAPPER_GENOME_DEFINITION_FORMAT))
		{
			throw DataException(line, string("Genome definition expected."));
		}
		*/
		
		// Alignment info
		line = readLine(this->stream, line, lineNum);
		
		if (!sscanf(line.c_str(), SNAPPER_ALIGNMENT_INFO_FORMAT, &readAlign->getAlignedSegment().first, &readAlign->getAlignedSegment().second, 
			&readAlign->getAlignment().first, &readAlign->getAlignment().second, &readAlign->getIdentity()))
		{
			throw DataException(line, string("Read definition expected."));
		}
		
		// Read sequence
		line = readLine(this->stream, line, lineNum);
		
		copy(line.begin(), line.end(), readAlign->getSequence().begin());
		
		// Genome sequence
		line = readLine(this->stream, line, lineNum);
		
		copy(line.begin(), line.end(), readAlign->getGenomeSequence().begin());
		
		// Read alignment end
		line = readLine(this->stream, line, lineNum);
		
		if (line != SNAPPER_READ_ALIGNMENT_END)
		{
			throw DataException(line, string("Read alignment end expected."));
		}
		
		this->data.push_back(readAlign);
	}
	
	fprintf(stderr, "Processed "F_U64" Snapper read alignment[s].\n\n", this->data.size());
}
