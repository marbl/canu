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

static const char* rcsid = "$Id: ReadAlignmentProfiler.C,v 1.4 2011-09-05 16:49:44 mkotelbajcvi Exp $";

#include "ReadAlignmentProfiler.h"

using namespace ReadAnalysis;

ReadAlignmentProfiler::ReadAlignmentProfiler()
{
	this->verbose = false;
	//this->bases.reserve(DEFAULT_BASES_RESERVE_SIZE);
}

ReadAlignmentProfiler::~ReadAlignmentProfiler()
{
	FileUtils::close(this->stream);
}

void ReadAlignmentProfiler::writeProfile(string path)
{
	if (this->bases.empty())
	{
		fprintf(stderr, "Read alignment data has not been processed.\n");
		
		return;
	}
	
	if (this->verbose)
	{
		fprintf(stderr, "Opening profile output file stream: "F_STR"\n", path.c_str());
	}
	
	this->stream = FileUtils::openWrite(path);
	
	this->writeProfile(this->stream);
	
	if (this->verbose)
	{
		fprintf(stderr, "Closing profile output file stream: "F_STR"\n", path.c_str());
	}
	
	FileUtils::close(this->stream);
	
	this->stream = NULL;
}

void ReadAlignmentProfiler::writeProfile(FILE* stream)
{
	if (this->bases.empty())
	{
		fprintf(stderr, "Read alignment data has not been processed.\n");
		
		return;
	}
	
	this->stream = stream;
	
	if (!FileUtils::canWrite(stream))
	{
		string errorStr;
		
		throw IOException("Profile output stream cannot be written to: " + ErrorUtils::getError(stream, errorStr));
	}
	
	clock_t startClock = clock();
	
	fprintf(stream, F_STR"numReadAlignments="F_U64"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX, this->data.size());
	fprintf(stream, F_STR"numBases="F_U64"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX, this->bases.size());
	
	fprintf(stream, F_STR"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX);
	fprintf(stream, F_STR"\"Base\"\t\"# of Reads\"\t\"# of Mismatch\"\t\"# of Insertion\"\t\"# of Deletion\"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX);
	
	BaseAlignment* baseAlign;
	
	for (size_t a = 0; a < this->bases.size(); a++)
	{
		baseAlign = this->bases[a];
		
		fprintf(stream, F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\t"F_U64"\n", a + 1, baseAlign->getReads().size(), 
			baseAlign->getErrors()[MISMATCH]->size(), baseAlign->getErrors()[INSERTION]->size(), 
			baseAlign->getErrors()[DELETION]->size());
	}
	
	fprintf(stderr, "Wrote profile of "F_U64" base[s]: time=%.2"F_F64P"s\n", this->bases.size(), 
		(double)(clock() - startClock) / (double)CLOCKS_PER_SEC);
}

void ReadAlignmentProfiler::profileData(vector<ReadAlignment*>& data, AlignmentDataStats& dataStats)
{
	clock_t startClock = clock();
	
	this->data = data;
	this->dataStats = dataStats;
	
	this->initBases();
	
	fprintf(stderr, "Generating profile of "F_U64" base[s]:\n", this->bases.size());
	
	ReadAlignment* readAlign;
	uint16 profiledPercent = 0, lastProfiledPercent = 0;
	
	for (size_t a = 0; a < this->data.size(); a++)
	{
		readAlign = this->data[a];
		
		for (size_t b = 0; b < readAlign->getLength(); b++)
		{
			this->profileBase(a, b);
		}
		
		profiledPercent = (uint16)(((double)a / (double)this->data.size()) * 100.0);
		
		if (profiledPercent >= (lastProfiledPercent + PROFILE_DATA_PERCENT_INCREMENT))
		{
			fprintf(stderr, "\t"F_U16"%% ("F_U64" read alignment[s]) completed ...\n", profiledPercent, a + 1);
			
			lastProfiledPercent = profiledPercent;
		}
	}
	
	fprintf(stderr, "Generated profile of "F_U64" base[s]: time=%.2"F_F64P"s\n", this->bases.size(), 
		(double)(clock() - startClock) / (double)CLOCKS_PER_SEC);
}

void ReadAlignmentProfiler::profileBase(size_t readIndex, size_t baseIndex)
{
	ReadAlignment* readAlign = this->data[readIndex];
	BaseAlignment* baseAlign = this->bases[baseIndex];
	
	char readBase = readAlign->getSequence()[baseIndex], genomeBase = readAlign->getGenomeSequence()[baseIndex];
	
	if (isBaseCall(readBase))
	{
		if (isBaseError(readBase))
		{
			if (isBaseError(genomeBase))
			{
				baseAlign->addError(new AlignmentError(readAlign->getIid(), MISMATCH));
			}
			else if (isBaseGap(genomeBase))
			{
				baseAlign->addError(new AlignmentError(readAlign->getIid(), INSERTION));
			}
		}
	}
	else if (isBaseGap(readBase))
	{
		baseAlign->addError(new AlignmentError(readAlign->getIid(), DELETION));
	}
	else if (isBaseUnknown(readBase))
	{
		baseAlign->addError(new AlignmentError(readAlign->getIid(), MISMATCH));
	}
	
	baseAlign->addRead(readAlign->getIid());
}

void ReadAlignmentProfiler::initBases()
{
	this->bases = vector<BaseAlignment*>();
	this->bases.reserve(this->dataStats.getMaxReadLength());
	
	for (size_t a = 0; a < this->bases.capacity(); a++)
	{
		// TODO: get reads per base stat to fill readsReserveSize
		BaseAlignment* baseAlign = new BaseAlignment(a, 1024, this->dataStats.getMaxReadLength() / 3);
		
		this->bases.push_back(baseAlign);
	}
}
