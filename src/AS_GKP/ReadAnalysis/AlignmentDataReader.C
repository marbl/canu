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

static const char* rcsid = "$Id: AlignmentDataReader.C,v 1.4 2011-09-05 16:49:44 mkotelbajcvi Exp $";

#include "AlignmentDataReader.h"

using namespace ReadAnalysis;

AlignmentDataReader::AlignmentDataReader()
{
	this->verbose = false;
	this->stream = NULL;
	
	this->data.reserve(DEFAULT_READ_ALIGNMENT_RESERVE_SIZE);
}

AlignmentDataReader::~AlignmentDataReader()
{
	FileUtils::close(this->stream);
	
	this->stream = NULL;
}

void AlignmentDataReader::readData(string path)
{
	if (this->verbose)
	{
		fprintf(stderr, "Opening alignment data file stream: "F_STR"\n", path.c_str());
	}
	
	this->stream = FileUtils::openRead(path);
	
	this->readData(this->stream);
	
	if (this->verbose)
	{
		fprintf(stderr, "Closing alignment data file stream: "F_STR"\n", path.c_str());
	}
	
	FileUtils::close(this->stream);
	
	this->stream = NULL;
}

void AlignmentDataReader::readData(FILE* stream)
{
	this->data.clear();
	
	this->stream = stream;
	
	if (!FileUtils::canRead(stream))
	{
		string errorStr;
		
		throw IOException("Alignment data input stream cannot be read: " + ErrorUtils::getError(stream, errorStr));
	}
	
	this->processData();
	
	this->dataStats.setMeanReadLength(this->dataStats.getMeanReadLength() / this->data.size());
}

void AlignmentDataReader::processReadAlignmentStats(ReadAlignment* readAlign)
{
	this->dataStats.setMinReadLength(min(this->dataStats.getMinReadLength(), readAlign->getLength()));
	this->dataStats.setMaxReadLength(max(this->dataStats.getMaxReadLength(), readAlign->getLength()));
	this->dataStats.setMeanReadLength(this->dataStats.getMeanReadLength() + readAlign->getLength());
}
