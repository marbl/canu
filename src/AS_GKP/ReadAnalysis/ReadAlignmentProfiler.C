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

static const char* rcsid = "$Id: ReadAlignmentProfiler.C,v 1.6 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "ReadAlignmentProfiler.h"

using namespace ReadAnalysis;

ReadAlignmentProfiler::ReadAlignmentProfiler()
{
	this->verbose = false;
	this->positionMode = DEFAULT;
	
	this->errorAssigner = &ReadAlignmentProfiler::defaultErrorAssigner;
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
	
	fprintf(stderr, "\nWriting profile of "F_SIZE_T" base[s] ...\n", this->bases.size());
	
	fprintf(stream, F_STR"numReadAlignments="F_SIZE_T"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX, this->data.size());
	fprintf(stream, F_STR"numFilteredReadAlignments="F_SIZE_T"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX, this->dataStats.getNumFilteredReads());
	fprintf(stream, F_STR"numBases="F_SIZE_T"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX, this->bases.size());
	
	fprintf(stream, F_STR"\n", PROFILE_DATA_OUTPUT_SUMMARY_PREFIX);
	
	fprintf(stream, F_STR"""\"Position\"\t\"# of Reads\"\t\"%% of Reads\"\t\"%% of Reads Error\"\t\"# of Mismatch\"\t\"# of Insertion\"\t\"# of Deletion\"\n", 
		PROFILE_DATA_OUTPUT_SUMMARY_PREFIX);
	
	BaseAlignment* baseAlign;
	size_t numBaseReads, numBaseMismatch, numBaseInsertion, numBaseDeletion;
	double percentBaseReads, percentBaseReadsError;
	
	for (size_t a = 0; a < this->bases.size(); a++)
	{
		baseAlign = this->bases[a];
		
		numBaseReads = this->dataStats.getNumReads(a);
		numBaseMismatch = baseAlign->getErrors()[MISMATCH]->size();
		numBaseInsertion = baseAlign->getErrors()[INSERTION]->size();
		numBaseDeletion = baseAlign->getErrors()[DELETION]->size();
		
		percentBaseReads = ((double)numBaseReads / (double)numBaseReads) * 100;
		percentBaseReadsError = ((double)(numBaseMismatch + numBaseInsertion + numBaseDeletion) / (double)numBaseReads) * 100;
		
		fprintf(stream, F_SIZE_T"\t"F_SIZE_T"\t%.3"F_F64P"\t%.3"F_F64P"\t"F_SIZE_T"\t"F_SIZE_T"\t"F_SIZE_T"\n", baseAlign->getPosition() + 1, numBaseReads, 
			percentBaseReads, percentBaseReadsError, numBaseMismatch, numBaseInsertion, numBaseDeletion);
	}
	
	string timeStr;
	
	fprintf(stderr, "Wrote profile of "F_SIZE_T" base[s] in: "F_STR"\n", this->bases.size(), 
		StringUtils::toString(timeStr, startClock, true).c_str());
}

void ReadAlignmentProfiler::profileData(vector<ReadAlignment*>& data, map<AS_IID, ReadAlignment*>& iidMap, 
	AlignmentDataStats& dataStats)
{
	clock_t startClock = clock();
	
	this->data = data;
	this->iidMap = iidMap;
	this->dataStats = dataStats;
	
	this->initBases();
	
	if (this->positionMode == END_DISTANCE)
	{
		this->errorAssigner = &ReadAlignmentProfiler::endDistanceErrorAssigner;
	}
	
	fprintf(stderr, "\nGenerating profile of "F_SIZE_T" base[s]:\n", this->bases.size());
	
	ReadAlignment* readAlign;
	double profiledPercent = 0, lastProfiledPercent = 0;
	
	for (size_t a = 0; a < this->data.size(); a++)
	{
		readAlign = this->data[a];
		
		for (size_t b = 0; b < readAlign->getLength(); b++)
		{
			this->profileBase(a, b);
		}
		
		profiledPercent = ((double)(a + 1) / (double)this->data.size()) * 100;
		
		if (profiledPercent >= (lastProfiledPercent + PROFILE_DATA_PERCENT_INCREMENT))
		{
			fprintf(stderr, "  %.0"F_F64P"%% ("F_SIZE_T" read alignment[s]) completed ...\n", profiledPercent, a + 1);
			
			lastProfiledPercent = profiledPercent;
		}
	}
	
	string timeStr;
	
	fprintf(stderr, "Generated profile of "F_SIZE_T" base[s] in: "F_STR"\n", this->bases.size(), 
		StringUtils::toString(timeStr, startClock, true).c_str());
}

void ReadAlignmentProfiler::profileBase(size_t readIndex, size_t baseIndex)
{
	ReadAlignment* readAlign = this->data[readIndex];
	
	char readBase = readAlign->getSequence()[baseIndex], genomeBase = readAlign->getGenomeSequence()[baseIndex];
	
	if (isBaseCall(readBase))
	{
		if (isBaseError(readBase))
		{
			if (isBaseError(genomeBase))
			{
				(this->*errorAssigner)(readAlign, new AlignmentError(readAlign->getIid(), MISMATCH), baseIndex);
			}
			else if (isBaseGap(genomeBase))
			{
				(this->*errorAssigner)(readAlign, new AlignmentError(readAlign->getIid(), INSERTION), baseIndex);
			}
		}
	}
	else if (isBaseGap(readBase))
	{
		(this->*errorAssigner)(readAlign, new AlignmentError(readAlign->getIid(), DELETION), baseIndex);
	}
	else if (isBaseUnknown(readBase))
	{
		(this->*errorAssigner)(readAlign, new AlignmentError(readAlign->getIid(), MISMATCH), baseIndex);
	}
}

void ReadAlignmentProfiler::initBases()
{
	this->bases = vector<BaseAlignment*>();
	this->bases.reserve(this->dataStats.getMaxReadLength());
	
	for (size_t a = 0; a < this->bases.capacity(); a++)
	{
		BaseAlignment* baseAlign = new BaseAlignment(a, this->dataStats.getMaxReadLength() / 3);
		
		this->bases.push_back(baseAlign);
	}
}

void ReadAlignmentProfiler::defaultErrorAssigner(ReadAlignment* readAlign, AlignmentError* error, 
	size_t baseIndex)
{
	this->bases[baseIndex]->addError(error);
}

void ReadAlignmentProfiler::endDistanceErrorAssigner(ReadAlignment* readAlign, AlignmentError* error, 
	size_t baseIndex)
{
	this->bases[readAlign->getLength() - baseIndex]->addError(error);
}
