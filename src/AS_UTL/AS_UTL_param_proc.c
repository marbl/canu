
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_param_proc.h"

static paramEntryT* firstParamEntry = NULL;
static paramEntryT* lastParamEntry = NULL;

static paramEntryT* allocParamEntry(void);
static void addParamEntry( paramEntryT * const newParam);

static paramEntryT* allocParamEntry(void)
{
  paramEntryT *newParam;
  
  newParam = (paramEntryT *) malloc( sizeof(paramEntryT));
  if (newParam != NULL)
	return newParam;
  else
  {
	fprintf( stderr, "Error allocating newParam in allocParamEntry()\n");
	assert(0);
  }
  return NULL; // Failure
}

static void addParamEntry( paramEntryT * const newParamEntry)
{
  if (firstParamEntry == NULL)
	firstParamEntry = newParamEntry;
  if (lastParamEntry != NULL)
	lastParamEntry->next = newParamEntry;
  lastParamEntry = newParamEntry;
}

char* getParam(char const * const paramNameIn)
{
  paramEntryT *currentParamEntry = firstParamEntry;
  int i, foundParamEntry = 0;
  char moduleName[INPUTMAX], paramName[INPUTMAX];
  int moduleNameLength, paramNameLength;
  char const *  paramNameInPtr = paramNameIn;

  moduleNameLength = 0;
  while (*paramNameInPtr != '.')
  {
	paramNameInPtr++;
	moduleNameLength++;
  }
  
  for ( i = 0; i < moduleNameLength; i++)
	moduleName[i] = paramNameIn[i];
  moduleName[i] = '\0';
  
  // move past '.'
  paramNameInPtr++;
  
  paramNameLength = 0;
  while (*paramNameInPtr != '\0')
  {
	paramNameInPtr++;
	paramNameLength++;
  }
  
  for ( i = 0; i < paramNameLength; i++)
	paramName[i] = paramNameIn[moduleNameLength + 1 + i];
  paramName[i] = '\0';
  
  while (currentParamEntry != NULL)
  {
	if (!strcmp( moduleName, currentParamEntry->paramModule))
	{
	  if (!strcmp( paramName, currentParamEntry->paramName))
	  {
		foundParamEntry = 1;
		break;
	  }
	}
	currentParamEntry = currentParamEntry->next;
  }
  
  if (foundParamEntry) {
	return (currentParamEntry->paramValue);
  }
  return (NULL);
}

int getAllParams(const char * const moduleNameIn, char **returnBuffer)
{
  paramEntryT *currentParamEntry = firstParamEntry;
  char *tempBufferPos;
  char tempBuffer[MAXRETURNBUFFERLENGTH];
  
  tempBufferPos = tempBuffer;
  
  while (currentParamEntry != NULL)
  {
	if (!strcmp(moduleNameIn, currentParamEntry->paramModule))
	{  
	  strncpy( tempBufferPos, currentParamEntry->paramModule, strlen(currentParamEntry->paramModule));
	  tempBufferPos += strlen(currentParamEntry->paramModule);
	  strncpy( tempBufferPos, ".", 1);
	  tempBufferPos++;
	  
	  strncpy( tempBufferPos, currentParamEntry->paramName, strlen(currentParamEntry->paramName));
	  tempBufferPos += strlen(currentParamEntry->paramName);
	  strncpy( tempBufferPos, " ", 1);
	  tempBufferPos++;
	  
	  strncpy( tempBufferPos, currentParamEntry->paramValue, strlen(currentParamEntry->paramValue));
	  tempBufferPos += strlen(currentParamEntry->paramValue);
	  strncpy( tempBufferPos, " ", 1);
	  tempBufferPos++;
	}
	currentParamEntry = currentParamEntry->next;
  }
  *returnBuffer = (char *) malloc(tempBufferPos - tempBuffer);
  strncpy( *returnBuffer, tempBuffer, tempBufferPos - tempBuffer);

  return( tempBufferPos - tempBuffer );
}

int loadParams(const char * const filename)
{
  FILE * paramFile;
  char buffer[ INPUTMAX ]; 
  char *currentChar;
  int i;
  int paramModuleLength, paramNameLength, paramValueLength;
  char *paramModuleStart, *paramNameStart, *paramValueStart;
  int verbose = 0;
  paramEntryT *newParamEntry;
  
  paramFile = fopen( filename, "r");
  if (paramFile == NULL)
  {
	fprintf(stderr, "Could not open %s\n", filename);
	return(PARAM_PROC_FAILURE);
  }
  
  while (fgets( buffer, INPUTMAX, paramFile) != NULL)
  {
	int done = 0;

	// fprintf( stderr, "raw: %s\n", buffer);
	currentChar = buffer;

	// ignore leading spaces
	while (isspace(*currentChar))
	  currentChar++;
	

	// check to see if line is a comment or empty (fgets sets null char past last char)
	if (*currentChar == '#' || *currentChar == 0)
	  done = 1;

	// now process rest of line
	// first grab everything up to =
	while (!done)
	{
	  ////// module //////

	  paramModuleLength = 0;
	  paramModuleStart = currentChar;
	  
	  // get module name
	  while ( isalnum(*currentChar))
	  {
		paramModuleLength++;
		currentChar++;
	  }
	  
	  if ( *currentChar != '.')
	  {
		fprintf( stderr, "Missing \'.\', line: %s\n", buffer);
		fprintf( stderr, "*currentChar: %c\n", *currentChar);
		assert(0);
	  }

	  if (verbose)
	  {
		fprintf( stderr, "paramModule: ");
		for ( i = 0; i < paramModuleLength; i++)
		  fprintf (stderr, "%c", paramModuleStart[i]);
		fprintf( stderr, "\n");
	  }
	  
	  // move past .
	  currentChar++;
	  
	  // check to make sure next character is not a space or any other baddie
	  if (!isalnum( *currentChar))
	  {
		fprintf( stderr, "Format error, line: %s\n", buffer);
		fprintf( stderr, "*currentChar: %c\n", *currentChar);
		assert(0);
	  }
		

	  ////// name //////
	  paramNameLength = 0;
	  paramNameStart = currentChar;
	
	  // get parameter name
	  while (isalnum(*currentChar) || ('_' == (*currentChar)))
	  {
		paramNameLength++;
		currentChar++;
	  }
	  
	  // ignore spaces
	  while (*currentChar == ' ')
		currentChar++;

	  // check for missing =
	  if (*currentChar != '=')
	  {
		fprintf( stderr, "Missing \'=\', line: %s\n", buffer);
		fprintf( stderr, "*currentChar: %c\n", *currentChar);
		assert(0);
	  }
	  
	  if (verbose)
	  {
		fprintf( stderr, "paramName: ");
		for ( i = 0; i < paramNameLength; i++)
		  fprintf (stderr, "%c", paramNameStart[i]);
		fprintf( stderr, "\n");
	  }
	  
	  // move past =
	  currentChar++;

	  // ignore spaces
	  while (*currentChar == ' ')
		currentChar++;


	  ////// value //////

	  paramValueLength = 0;
	  paramValueStart = currentChar;
	  
	  while (isalnum(*currentChar) || ('_' == (*currentChar)))
	  {
		paramValueLength++;
		currentChar++;
	  }

	  if (verbose)
	  {
		fprintf( stderr, "paramValue: ");
		for ( i = 0; i < paramValueLength; i++)
		  fprintf (stderr, "%c", paramValueStart[i]);
		fprintf( stderr, "\n");
	  }
	  
	  newParamEntry = allocParamEntry();
	  strncpy( newParamEntry->paramModule, paramModuleStart, paramModuleLength);
	  newParamEntry->paramModule[paramModuleLength] = '\0';

	  strncpy( newParamEntry->paramName, paramNameStart, paramNameLength);
	  newParamEntry->paramName[paramNameLength] = '\0';

	  strncpy( newParamEntry->paramValue, paramValueStart, paramValueLength);
	  newParamEntry->paramValue[paramValueLength] = '\0';
	  
	  addParamEntry( newParamEntry );
	  
	  done = 1;	  
	}
  }
  fclose (paramFile);
  return(PARAM_PROC_SUCCESS);
}

void main_test(void)
{
  char *allParams;
  int length;
  
  loadParams("/work/assembly/flanigmj/test_12_9/cds/AS/src/AS_UTL/mymod");
  loadParams("/work/assembly/flanigmj/test_12_9/cds/AS/src/AS_UTL/yourmod");

  fprintf( stderr, "\n Searching...\n");
  fprintf( stderr, "mymod.param1: %s\n", getParam("mymod.param1"));
  
  length = getAllParams( "mymod", &allParams );
  fprintf( stderr, "allParams: %s\n", allParams);
  fprintf( stderr, "length: %d\n", length);
  
  length = getAllParams( "yourmod", &allParams );
  fprintf( stderr, "allParams: %s\n", allParams);
  fprintf( stderr, "length: %d\n", length);
  
}
