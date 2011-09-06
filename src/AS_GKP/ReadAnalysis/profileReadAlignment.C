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

static const char* rcsid = "$Id: profileReadAlignment.C,v 1.5 2011-09-06 09:47:55 mkotelbajcvi Exp $";

#include "profileReadAlignment.h"

void parseFilterArgs(ProfileReadAlignmentOptions& options, vector<AlignmentDataFilter*>& filters)
{
	string filterArg;
	vector<string> params;
	
	for (size_t a = 0; a < options.getFilterArgs().size(); a++)
	{
		filterArg = options.getFilterArgs()[a];
		StringUtils::trim(filterArg, 2, " ", "\t");
		
		params = vector<string>();
		StringUtils::split(filterArg, params, FILTER_PARAM_DELIMITER);
		
		if (!params.empty() && (params[0] == "length"))
		{
			ReadLengthAlignmentDataFilter* filter = new ReadLengthAlignmentDataFilter();
			
			size_t minLength, maxLength;
			
			if ((params.size() >= 2) && sscanf(params[1].c_str(), F_SIZE_T, &minLength))
			{
				filter->setMinLength(minLength);
			}
			
			if ((params.size() >= 3) && sscanf(params[2].c_str(), F_SIZE_T, &maxLength))
			{
				filter->setMaxLength(maxLength);
			}
			
			filters.push_back(filter);
		}
	}
}

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
						case 'f':
							options.getFilterArgs().push_back(argValue);
							argIndex++;
							break;
							
						case 'm':
							if (StringUtils::areEqual(argValue, "enddist"))
							{
								options.setPositionMode(END_DISTANCE);
							}
							else if (StringUtils::areEqual(argValue, "default"))
							{
								options.setPositionMode(DEFAULT);
							}
							else
							{
								// TODO: add error
							}
							argIndex++;
							break;
							
						case 'i':
							options.setInputPath(argValue);
							argIndex++;
							break;
							
						case 'o':
							options.setOutputPath(argValue);
							argIndex++;
							break;
							
						case 'v':
							options.setVerbose(true);
							break;
							
						default:
							throw ArgumentException(arg, "Unknown option.");
							break;
					}
				}
				
				argIndex++;
			}
			
			if (!options.useStdin() && !FileUtils::isReadable(options.getInputPath()))
			{
				throw ArgumentException("-i", string("A readable alignment data file or '-' for stdin must be specified: ") + options.getInputPath());
			}
			
			if (!options.useStdout() && !FileUtils::isWriteable(options.getOutputPath()))
			{
				throw ArgumentException("-o", string("A writable output file or '-' for stdout must be specified: ") + options.getOutputPath());
			}
		}
		catch (RuntimeException& e)
		{
			fprintf(stderr, F_STR"\n", options.getVerbose() ? e.what() : e.getMessage().c_str());
			
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
	fprintf(stderr, "\n");
	fprintf(stderr, "Creates a profile of read alignments by analysing alignment data.\n");
	fprintf(stderr, "Known alignment data formats: Snapper\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: "F_STR" {-v} {-f <filter params>}* {-m <mode>} -i <input path> -o <output path>\n", executableName);
	fprintf(stderr, "\n");
	fprintf(stderr, "-v  optional, whether to produce verbose log output\n");
	fprintf(stderr, "-f  optional, multiple, specifies a filter to apply (see below)\n");
	fprintf(stderr, "-m  optional, specifies a positioning mode (see below)\n");
	fprintf(stderr, "-i  path to an alignment data file or '-' to read stdin\n");
	fprintf(stderr, "-o  path to an output file or '-' to write to stdout\n");
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Available Filters:\n");
	fprintf(stderr, "------------------\n");
	fprintf(stderr, "length{,minLength}{,maxLength}\n");
	fprintf(stderr, "    skips all read alignments with a length outside of a set of bounds\n");
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Positioning Modes:\n");
	fprintf(stderr, "------------------\n");
	fprintf(stderr, "'default'  positions are sequence indexes\n");
	fprintf(stderr, "'enddist'  positions are relative to their distance from the end of the read\n");
	fprintf(stderr, "\n");
	
	exitSuccess();
}

int main(int numArgs, char** args)
{
	ErrorUtils::handleErrorSignals(ErrorUtils::printingSignalHandler);
	
	ProfileReadAlignmentOptions options;
	parseCommandLine(options, numArgs, args);
	
	fprintf(stderr, "Using options:\n  verbose="F_STR"\n  positionMode="F_STR"\n  input="F_STR"\n  output="F_STR"\n\n", 
		options.getVerbose() ? "true" : "false", 
		(options.getPositionMode() == END_DISTANCE) ? "<end distance>" : "<default>", 
		options.useStdin() ? "<stdin>" : options.getInputPath().c_str(), 
		options.useStdout() ? "<stdout>" : options.getOutputPath().c_str());
	
	try
	{
		SnapperAlignmentDataReader reader;
		reader.setVerbose(options.getVerbose());
		
		vector<AlignmentDataFilter*> filters;
		parseFilterArgs(options, filters);
		
		AlignmentDataFilter* filter;
		
		if (!filters.empty())
		{
			fprintf(stderr, "Using filters:");
			
			for (size_t a = 0; a < filters.size(); a++)
			{
				filter = filters[a];
				filter->setDataStats(&reader.getDataStats());
				
				fprintf(stderr, "\n  "F_SIZE_T"="F_STR, a + 1, filter->toString().c_str());
				
				reader.getFilters().push_back(filter);
			}
			
			fprintf(stderr, "\n\n");
		}
		
		if (options.useStdin())
		{
			reader.readData(stdin);
		}
		else
		{
			reader.readData(options.getInputPath());
		}
		
		ReadAlignmentProfiler profiler;
		profiler.setVerbose(options.getVerbose());
		profiler.setPositionMode(options.getPositionMode());
		
		profiler.profileData(reader.getData(), reader.getIidMap(), reader.getDataStats());
		
		if (options.useStdout())
		{
			profiler.writeProfile(stdout);
		}
		else
		{
			profiler.writeProfile(options.getOutputPath());
		}
	}
	catch (RuntimeException& e)
	{
		fprintf(stderr, F_STR"\n", options.getVerbose() ? e.what() : e.getMessage().c_str());
		
		exitFailure();
	}
}
