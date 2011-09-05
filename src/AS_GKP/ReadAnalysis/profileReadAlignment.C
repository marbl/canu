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

static const char* rcsid = "$Id: profileReadAlignment.C,v 1.4 2011-09-05 16:49:44 mkotelbajcvi Exp $";

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
							options.inputPath = argValue;
							argIndex++;
							break;
							
						case 'o':
							options.outputPath = argValue;
							argIndex++;
							break;
							
						case 'v':
							options.verbose = true;
							break;
							
						default:
							throw ArgumentException(arg, "Unknown option.");
							break;
					}
				}
				
				argIndex++;
			}
			
			options.useStdin = isPathStream(options.inputPath);
			options.useStdout = isPathStream(options.outputPath);
			
			if (!options.useStdin && !FileUtils::isReadable(options.inputPath))
			{
				throw ArgumentException("-i", string("A readable alignment data file or '-' for stdin must be specified: ") + options.inputPath);
			}
			
			if (!options.useStdout && !FileUtils::isWriteable(options.outputPath))
			{
				throw ArgumentException("-o", string("A writable output file or '-' for stdout must be specified: ") + options.outputPath);
			}
		}
		catch (RuntimeException& e)
		{
			fprintf(stderr, F_STR"\n", options.verbose ? e.what() : e.getMessage().c_str());
			
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
	fprintf(stdout, "Usage: %s {-v} -i <input path> -o <output path>\n", executableName);
	fprintf(stdout, "\n");
	fprintf(stdout, "-i  path to an alignment data file or '-' to read stdin\n");
	fprintf(stdout, "-o  path to an output file or '-' to write to stdout\n");
	fprintf(stdout, "-v  optional, whether to produce verbose log output\n");
	
	exitSuccess();
}

bool isPathStream(string path)
{
	return path == STREAM_PATH;
}

int main(int numArgs, char** args)
{
	ErrorUtils::handleErrorSignals(ErrorUtils::printingSignalHandler);
	
	ProfileReadAlignmentOptions options;
	parseCommandLine(options, numArgs, args);
	
	fprintf(stderr, "Using options:\n  verbose="F_STR"\n  input="F_STR"\n  output="F_STR"\n\n", 
		options.verbose ? "true" : "false", 
		options.useStdin ? "<stdin>" : options.inputPath.c_str(), 
		options.useStdout ? "<stdout>" : options.outputPath.c_str());
	
	try
	{
		SnapperAlignmentDataReader reader;
		reader.setVerbose(options.verbose);
		
		if (options.useStdin)
		{
			reader.readData(stdin);
		}
		else
		{
			reader.readData(options.inputPath);
		}
		
		ReadAlignmentProfiler profiler;
		profiler.setVerbose(options.verbose);
		
		profiler.profileData(reader.getData(), reader.getDataStats());
		
		if (options.useStdout)
		{
			profiler.writeProfile(stdout);
		}
		else
		{
			profiler.writeProfile(options.outputPath);
		}
	}
	catch (RuntimeException& e)
	{
		fprintf(stderr, F_STR"\n", options.verbose ? e.what() : e.getMessage().c_str());
		
		exitFailure();
	}
}
