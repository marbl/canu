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

#ifndef PROFILEREADALIGNMENT_H
#define PROFILEREADALIGNMENT_H

static const char* rcsid_PROFILEREADALIGNMENT_H = "$Id$";

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <set>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AlignmentDataFilter.h"
#include "AlignmentDataReader.h"
#include "ArgumentException.h"
#include "AS_UTL_fileIO.h"
#include "ErrorUtils.h"
#include "FileUtils.h"
#include "IllegalStateException.h"
#include "ReadAlignmentProfiler.h"
#include "ReadLengthAlignmentDataFilter.h"
#include "RuntimeException.h"
#include "SnapperAlignmentDataReader.h"

using namespace ReadAnalysis;
using namespace Utility;

static const char* STREAM_PATH = "-";
static const char* FILTER_PARAM_DELIMITER = ",";

class ProfileReadAlignmentOptions
{
public:
	ProfileReadAlignmentOptions()
	{
		this->verbose = false;
		this->positionMode = DEFAULT;
	}
	
	static bool isPathStream(string path)
	{
		return path == STREAM_PATH;
	}
	
	bool useStdin()
	{
		return isPathStream(this->inputPath);
	}
	
	bool useStdout()
	{
		return isPathStream(this->outputPath);
	}
	
	bool getVerbose()
	{
		return this->verbose;
	}
	
	void setVerbose(bool verbose)
	{
		this->verbose = verbose;
	}
	
	vector<string>& getFilterArgs()
	{
		return this->filterArgs;
	}
	
	BasePositionMode getPositionMode()
	{
		return this->positionMode;
	}
	
	void setPositionMode(BasePositionMode positionMode)
	{
		this->positionMode = positionMode;
	}
	
	string getInputPath()
	{
		return this->inputPath;
	}
	
	void setInputPath(string inputPath)
	{
		this->inputPath = inputPath;
	}
	
	string getOutputPath()
	{
		return this->outputPath;
	}
	
	void setOutputPath(string outputPath)
	{
		this->outputPath = outputPath;
	}
	
protected:
	bool verbose;
	vector<string> filterArgs;
	BasePositionMode positionMode;
	string inputPath;
	string outputPath;
};

void parseFilterArgs(ProfileReadAlignmentOptions& options, vector<AlignmentDataFilter*>& filters);
void parseCommandLine(ProfileReadAlignmentOptions& options, int numArgs, char** args);
void printUsage(const char* executableName);

int main(int numArgs, char** args);

#endif
