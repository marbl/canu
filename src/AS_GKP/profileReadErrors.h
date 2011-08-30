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

#ifndef PROFILEREADERRORS_H
#define PROFILEREADERRORS_H

static const char* RCSID_PROFILEREADERRORS_H = "$Id: profileReadErrors.h,v 1.1 2011-08-30 23:09:51 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <exception>
#include <map>
#include <string>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_IID.h"
#include "FileUtils.h"
#include "IllegalStateException.h"
#include "IOException.h"
#include "RuntimeException.h"
#include "StringUtils.h"

#define SNAPPER_FILE_LINE_BUFFER_SIZE 1536

#define INITIAL_ERROR_MATRIX_SIZE 1536
#define INITIAL_ERROR_MATRIX_BUCKET_SIZE 1024
#define ERROR_MATRIX_RESIZE_RESERVE_FACTOR 2

#define readSnapperFileLine(snapperFile, line, lineNum) \
	line = new char[SNAPPER_FILE_LINE_BUFFER_SIZE]; \
	line = FileUtils::readLine(snapperFile, line, SNAPPER_FILE_LINE_BUFFER_SIZE); \
	lineNum++;

#define isBaseCall(base) \
	(base == 'a') || (base == 'A') || \
	(base == 'c') || (base == 'C') || \
	(base == 'g') || (base == 'G') || \
	(base == 't') || (base == 'T')

#define isBaseError(base) \
	(base == 'A') || (base == 'C') || \
	(base == 'G') || (base == 'T')

#define isBaseGap(base) \
	base == '-'

#define isBaseUnknown(base) \
	base == 'N'

typedef enum AlignmentErrorType
{
	UNKNOWN, MISMATCH, INSERTION, DELETION
};

class AlignmentError
{
public:
	AlignmentError(AS_IID iid = 0, AlignmentErrorType type = UNKNOWN)
	{
		this->iid = iid;
		this->type = type;
	}
	
	AS_IID getIID()
	{
		return this->iid;
	}
	
	void setIID(AS_IID iid)
	{
		this->iid = iid;
	}
	
	AlignmentErrorType getType()
	{
		return this->type;
	}
	
	void setType(AlignmentErrorType type)
	{
		this->type = type;
	}
	
protected:
	AS_IID iid;
	AlignmentErrorType type;
};

void writeOutput(const char* outputFile, vector< vector<AlignmentError>* >& errorMatrix);

void processReadAlignment(AS_IID readIID, uint16 readLength, const char* readSequence, const char* genomeSequence, 
	vector< vector<AlignmentError>* >& errorMatrix);
void processSnapperFile(const char* snapperFile, map<AS_IID, uint16>& readMap, vector< vector<AlignmentError>* >& errorMatrix);

vector<AlignmentError>* getBaseErrorBucket(vector< vector<AlignmentError>* >& errorMatrix, size_t base);
AS_IID getReadIID(const char* readDefLine);
uint16 getReadLength(string readInfoLine);

void printUsage(const char* executableName);

int main(int argc, char** argv);

#endif
