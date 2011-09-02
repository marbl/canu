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

static const char* rcsid_PROFILEREADALIGNMENT_H = "$Id: profileReadAlignment.h,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace std;

#include "ArgumentException.h"
#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "ErrorUtils.h"
#include "IllegalStateException.h"
#include "ReadAlignmentProfiler.h"
#include "RuntimeException.h"

using namespace ReadAnalysis;
using namespace Utility;

typedef struct ProfileReadAlignmentOptions
{
	const char* inputFilePath;
	const char* outputFilePath;
};

void parseCommandLine(ProfileReadAlignmentOptions& options, int numArgs, char** args);
void printUsage(const char* executableName);

int main(int numArgs, char** args);

#endif
