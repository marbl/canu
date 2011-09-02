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

static const char* rcsid = "$Id: AlignmentDataReader.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "AlignmentDataReader.h"

using namespace ReadAnalysis;

AlignmentDataReader::AlignmentDataReader()
{
	this->stream = NULL;
	
	this->data.reserve(DEFAULT_READ_ALIGNMENT_RESERVE_SIZE);
}

AlignmentDataReader::~AlignmentDataReader()
{
	if (this->stream != NULL)
	{
		fclose(this->stream);
		
		this->stream = NULL;
	}
}

vector<ReadAlignment>& AlignmentDataReader::readData(const char* filePath)
{
	try
	{
		if ((filePath == NULL) || !AS_UTL_fileExists(filePath, 0, 1))
		{
			throw IOException(string("Alignment data file does not exist: ") + filePath);
		}
		
		FILE* stream = fopen(filePath, "r");
		
		ErrorUtils::throwIfError<IOException>(stream, string("Alignment data file cannot be opened for reading: ") + filePath);
		
		this->readData(stream, filePath);
		
		fclose(stream);
		
		return this->data;
	}
	catch (RuntimeException& e)
	{
		if (this->stream != NULL)
		{
			fclose(this->stream);
			
			this->stream = NULL;
		}
		
		throw e;
	}
}

vector<ReadAlignment>& AlignmentDataReader::readData(FILE* stream, const char* filePath)
{
	try
	{
		this->data.clear();
		
		this->stream = stream;
		
		if ((this->stream == NULL) || feof(this->stream))
		{
			throw new IOException((filePath != NULL) ? 
				string("Alignment data file cannot be read: ") + filePath : 
				string("Alignment data input stream cannot be read."));
		}
		
		this->processData(stream, filePath);
		
		return this->data;
	}
	catch (RuntimeException& e)
	{
		if (this->stream != NULL)
		{
			fclose(this->stream);
			
			this->stream = NULL;
		}
		
		throw e;
	}
}
