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

#include "profileReadErrors.h"

static const char* rcsid = "$Id: profileReadErrors.C,v 1.4 2011-08-31 18:32:05 mkotelbajcvi Exp $";

void writeOutput(const char* outputFile, map<AS_IID, uint16>& readMap, vector<BasePosition*>& errorMatrix)
{
	FILE* outputFileHandle = NULL;
	
	try
	{
		outputFileHandle = fopen(outputFile, "w");
		
		if (outputFileHandle == NULL)
		{
			throw IOException(string("Unable to open output file: ") + outputFile);
		}
		
		fprintf(outputFileHandle, F_STR"\n", 
			"\"Base\"\t\"Relative Errors (%)\"\t\"Average Read Position (%)\"\t\"% Mismatch\"\t\"% Insertion\"\t\"% Deletion\"");
		
		for (size_t a = 0; a < errorMatrix.size(); a++)
		{
			BasePosition* baseErrorBucket = errorMatrix[a];
			
			if (baseErrorBucket != NULL)
			{
				size_t numMismatch = 0, numInsertion = 0, numDeletion = 0;
				
				for (size_t b = 0; b < baseErrorBucket->errors.size(); b++)
				{
					switch (baseErrorBucket->errors[b].type)
					{
						case MISMATCH:
							numMismatch++;
							break;
							
						case INSERTION:
							numInsertion++;
							break;
							
						case DELETION:
							numDeletion++;
							break;
							
						default:
							// TODO: add error
							break;
					}
				}
				
				fprintf(outputFileHandle, F_U64"\t"F_F64"\t"F_F64"\t"F_F64"\t"F_F64"\t"F_F64"\n", 
					a + 1, 
					(double)baseErrorBucket->errors.size() / baseErrorBucket->readsAtPosition, 
					baseErrorBucket->readLengthPercent, 
					(double)numMismatch / baseErrorBucket->errors.size(), 
					(double)numInsertion / baseErrorBucket->errors.size(), 
					(double)numDeletion / baseErrorBucket->errors.size());
			}
		}
		
		fclose(outputFileHandle);
		
		printf("Wrote "F_U64" base error counts to: "F_STR"\n", errorMatrix.size(), outputFile);
	}
	catch (RuntimeException& e)
	{
		if (outputFileHandle != NULL)
		{
			fclose(outputFileHandle);
		}
		
		fprintf(stderr, F_STR"\n", e.getMessage().c_str());
		
		exitFailure();
	}
}

void processReadAlignment(AS_IID readIID, uint16 readLength, const char* readSequence, const char* genomeSequence, 
	vector<BasePosition*>& errorMatrix)
{
	char readBase, genomeBase;
	BasePosition* baseErrorBucket;
	AlignmentError error;
	
	for (size_t a = 0; a < readLength; a++)
	{
		try
		{
			baseErrorBucket = getBaseErrorBucket(errorMatrix, a);
			readBase = readSequence[a];
			
			if (isBaseCall(readBase))
			{
				if (isBaseError(readBase))
				{
					error = AlignmentError(readIID);
					
					genomeBase = genomeSequence[a];
					
					if (isBaseError(genomeBase))
					{
						error.type = MISMATCH;
					}
					else if (isBaseGap(genomeBase))
					{
						error.type = INSERTION;
					}
					else
					{
						string baseStr, readIIDStr;
						
						throw IllegalStateException("Unknown genome base at " + StringUtils::toString(a, baseStr) + " in read " + 
							StringUtils::toString(readIID, readIIDStr) + ": " + genomeBase);
					}
					
					baseErrorBucket->addError(error, readLength);
				}
				
				baseErrorBucket->readsAtPosition++;
			}
			else if (isBaseGap(readBase))
			{
				baseErrorBucket->addError(AlignmentError(readIID, DELETION), readLength);
				
				baseErrorBucket->readsAtPosition++;
			}
			else if (isBaseUnknown(readBase))
			{
				baseErrorBucket->addError(AlignmentError(readIID, MISMATCH), readLength);
				
				baseErrorBucket->readsAtPosition++;
			}
			else if (!iscntrl(readBase))
			{
				string baseStr, readIIDStr;
				
				throw IllegalStateException("Unknown read base at " + StringUtils::toString(a, baseStr) + " in read " + 
					StringUtils::toString(readIID, readIIDStr) + ": " + readBase);
			}
		}
		catch (RuntimeException& e)
		{
			fprintf(stderr, F_STR"\n", e.getMessage().c_str());
		}
	}
}

void processSnapperFile(const char* snapperFile, map<AS_IID, uint16>& readMap, vector<BasePosition*>& errorMatrix)
{
	FILE* snapperFileHandle = NULL;
	
	try
	{
		snapperFileHandle = fopen(snapperFile, "r");
		
		if (snapperFileHandle == NULL)
		{
			throw IOException(string("Unable to open snapper alignment file: ") + snapperFile);
		}
		
		size_t numReads = 0;
		AS_IID readIID;
		uint16 readLength;
		const char* readSequence, *genomeSequence;
		
		size_t lineNum = 0;
		char* line;
		
		while (!feof(snapperFileHandle))
		{
			// Read section start
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			// Ignore trailing line
			if (StringUtils::isEmpty(line))
			{
				break;
			}
			
			if (!StringUtils::areEqual(line, "sim4begin"))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed read section start on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			// Read info
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (StringUtils::isBlank(line))
			{
				string lineNumStr;
				
				throw IllegalStateException("Missing expected content on line " + StringUtils::toString(lineNum, lineNumStr) + ".");
			}
			
			readLength = getReadLength(string(line));
			
			// Read definition
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (!StringUtils::startsWith(line, 1, "edef="))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed read definition on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			readIID = getReadIID(line);
			
			readMap[readIID] = readLength;
			
			// Genome definition
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (!StringUtils::startsWith(line, 1, "ddef="))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed genome definition on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			// Alignment info
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (StringUtils::isBlank(line))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed alignment information on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			// Read sequence
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (StringUtils::isBlank(line))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed read sequence on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			readSequence = line;
			
			// Genome sequence
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (StringUtils::isBlank(line))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed genome sequence on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			genomeSequence = line;
			
			processReadAlignment(readIID, readLength, readSequence, genomeSequence, errorMatrix);
			
			// Read section end
			readSnapperFileLine(snapperFileHandle, line, lineNum);
			
			if (!StringUtils::areEqual(line, "sim4end"))
			{
				string lineNumStr;
				
				throw IllegalStateException("Malformed read section end on line " + StringUtils::toString(lineNum, lineNumStr) + ": " + line);
			}
			
			numReads++;
		}
		
		fclose(snapperFileHandle);
		
		printf("Processed "F_U64" read alignment[s].\n", numReads);
	}
	catch (RuntimeException& e)
	{
		if (snapperFileHandle != NULL)
		{
			fclose(snapperFileHandle);
		}
		
		fprintf(stderr, F_STR"\n", e.getMessage().c_str());
		
		exitFailure();
	}
}

BasePosition* getBaseErrorBucket(vector<BasePosition*>& errorMatrix, size_t base)
{
	BasePosition* baseErrorBucket = errorMatrix[base];
	
	if (baseErrorBucket == NULL)
	{
		baseErrorBucket = new BasePosition(base);
		
		if (errorMatrix.size() <= base)
		{
			errorMatrix.reserve(base * ERROR_MATRIX_RESIZE_RESERVE_FACTOR);
			errorMatrix.resize(base);
		}
		
		errorMatrix[base] = baseErrorBucket;
	}
	
	return baseErrorBucket;
}

AS_IID getReadIID(const char* readDefLine)
{
	AS_IID readIID;
	
	sscanf(readDefLine, "edef="F_U32, &readIID);
	
	return readIID;
}

uint16 getReadLength(string readInfoLine)
{
	size_t lengthStartIndex = readInfoLine.find('[') + 1;
	uint16 readLength;
	
	sscanf(readInfoLine.substr(lengthStartIndex, readInfoLine.find('-') - lengthStartIndex).c_str(), F_U16, &readLength);
	
	return readLength;
}

void printUsage(const char* executableName)
{
	fprintf(stdout, "Generates histogram data for read alignment errors.\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "Usage: %s -s <snapper file> -o <output file>\n", executableName);
	fprintf(stdout, "\n");
	fprintf(stdout, "-s  path to a snapper alignment file\n");
	fprintf(stdout, "-o  path to an output file\n");
}

int main(int numArgs, char** args)
{
	if (numArgs > 1)
	{
		try
		{
			const char* arg, *argValue, *snapperFile = NULL, *outputFile = NULL;
			int argIndex = 1;
			
			while (argIndex < numArgs)
			{
				arg = args[argIndex];
				argValue = ((argIndex + 1) < numArgs) ? args[argIndex + 1] : NULL;
				
				if (!StringUtils::isBlank(arg) && (strlen(arg) == 2))
				{
					switch (arg[1])
					{
						case 's':
							snapperFile = argValue;
							break;
						case 'o':
							outputFile = argValue;
							break;
						default:
							fprintf(stderr, "Unknown option: "F_STR"\n", arg);
							
							exitFailure();
							break;
					}
				}
				
				argIndex += 2;
			}
			
			if ((snapperFile == NULL) && AS_UTL_fileExists(snapperFile, 0, 1))
			{
				throw ArgumentException(string("An existing snapper alignment file must be provided: ") + snapperFile, NULL, "-s");
			}
			
			if (outputFile == NULL)
			{
				throw ArgumentException(string("A writable output file must be provided: ") + outputFile, NULL, "-o");
			}
			
			printf("Using files:\n  snapper="F_STR"\n  output="F_STR"\n\n", snapperFile, outputFile);
			
			map<AS_IID, uint16> readMap;
			
			vector<BasePosition*> errorMatrix;
			errorMatrix.reserve(INITIAL_ERROR_MATRIX_SIZE);
			
			processSnapperFile(snapperFile, readMap, errorMatrix);
			
			writeOutput(outputFile, readMap, errorMatrix);
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
