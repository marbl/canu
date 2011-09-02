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

static const char* rcsid = "$Id: profileReadAlignment.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "profileReadAlignment.h"

void parseCommandLine(ProfileReadAlignmentOptions& options, int numArgs, char** args)
{
	if (numArgs > 1)
	{
		try
		{
			const char* arg, *argValue;
			int argIndex = 1;
			
			while (argIndex < numArgs)
			{
				arg = args[argIndex];
				argValue = ((argIndex + 1) < numArgs) ? args[argIndex + 1] : NULL;
				
				if (!StringUtils::isBlank(arg) && (strlen(arg) == 2))
				{
					switch (arg[1])
					{
						case 'i':
							options.inputFilePath = argValue;
							break;
						case 'o':
							options.outputFilePath = argValue;
							break;
						default:
							throw ArgumentException(string("Unknown option."), NULL, arg);
							break;
					}
				}
				
				argIndex += 2;
			}
			
			if ((options.inputFilePath == NULL) || (!StringUtils::areEqual(options.inputFilePath, "-") && !AS_UTL_fileExists(options.inputFilePath, 0, 1)))
			{
				throw ArgumentException(string("An existing alignment file or '-' for stdin must be specified: ") + options.inputFilePath, NULL, "-i");
			}
			
			if (options.outputFilePath == NULL)
			{
				throw ArgumentException(string("A writable output file or '-' for stdout must be specified: ") + options.outputFilePath, NULL, "-o");
			}
		}
		catch (RuntimeException& e)
		{
			fprintf(stderr, F_STR"\n", e.getMessage().c_str());
			
			exitFailure();
		}
	}
	else
	{
		printUsage(args[0]);
	}
}

void printUsage(const char* executableName)
{
	fprintf(stdout, "Creates a profile of read alignments by analysing alignment data.\n");
	fprintf(stdout, "Known alignment data formats: Snapper\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage: %s -i <alignment file> -o <output file>\n", executableName);
	fprintf(stdout, "\n");
	fprintf(stdout, "-i  path to an alignment file or '-' to read stdin\n");
	fprintf(stdout, "-o  path to an output file or '-' to write to stdout\n");
	
	exitSuccess();
}

int main(int numArgs, char** args)
{
	ErrorUtils::handleErrorSignals(ErrorUtils::printingSignalHandler);
	
	ProfileReadAlignmentOptions options;
	parseCommandLine(options, numArgs, args);
	
	fprintf(stderr, "Using options:\n  input="F_STR"\n  output="F_STR"\n\n", 
		StringUtils::areEqual(options.inputFilePath, "-") ? "<stdin>" : options.inputFilePath, 
		StringUtils::areEqual(options.outputFilePath, "-") ? "<stdout>" : options.outputFilePath);
}
