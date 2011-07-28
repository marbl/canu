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

static const char *rcsid = "$Id: upgrade-v7-to-v8.C,v 1.2 2011-07-28 02:07:07 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <map>

using namespace std;

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_fileIO.h"

#define INF_STORE_FILENAME "inf"
#define LIB_STORE_FILENAME "lib"
#define STORE_RENAME_SUFFIX ".v7-old"

class gkLibraryOld
{
public:
	AS_UID libraryUID;

	double mean;
	double stddev;

	uint64 spareUN3;
	uint64 spareUN2;
	uint64 spareUN1;

	uint64 spareUTG:62;
	uint64 forceBOGunitigger:1;
	uint64 isNotRandom:1;

	uint64 spareALN:63;
	uint64 doNotTrustHomopolymerRuns:1;

	uint64 spareOBT:54;

	uint64 doTrim_initialNone:1;
	uint64 doTrim_initialMerBased:1;
	uint64 doTrim_initialFlowBased:1;
	uint64 doTrim_initialQualityBased:1;

	uint64 doRemoveDuplicateReads:1;

	uint64 doTrim_finalLargestCovered:1;
	uint64 doTrim_finalEvidenceBased:1;

	uint64 doRemoveSpurReads:1;
	uint64 doRemoveChimericReads:1;

	uint64 doConsensusCorrection:1;

	uint64 spareGKP:63;
	uint64 UNUSEDusePackedFragments:1;

	uint64 spareLIB:61;
	uint64 orientation:3;
};

map<AS_IID, const char*> dumpMap;

void lookupLibraryName(AS_IID iid, gkLibrary* lib)
{
	assert(dumpMap.count(iid) > 0);
	
	strcpy(lib->libraryName, dumpMap[iid]);
}

void processInfoStore(char* oldInfoStorePath, char* infoStorePath)
{
	gkStoreInfo info;
	
	errno = 0;
	
	FILE* oldInfoFile = fopen(oldInfoStorePath, "r");
	
	if (AS_UTL_safeRead(oldInfoFile, &info, "inf", sizeof(info), 1))
	{
		assert(info.gkVersion == 7);
		
		info.gkVersion = 8;
	}
	else
	{
		fprintf(stderr, "Unable to read information store:\nfile=%s\nerror=%s\n", oldInfoStorePath, strerror(errno));
		exit(1);
	}
	
	FILE* infoFile = fopen(infoStorePath, "w");
	
	AS_UTL_safeWrite(infoFile, &info, "inf", sizeof(info), 1);
	
	fprintf(stdout, "Information store upgraded to v8:\n%s\n", infoStorePath);
}

void processLibStore(char* oldLibStorePath, char* libStorePath)
{
	StoreStruct* oldLibStore = openStore(oldLibStorePath, "r"),
			*libStore = createIndexStore(libStorePath, "lib", sizeof(gkLibrary), 1);
	
	for (AS_IID a = getFirstElemStore(oldLibStore); a < getLastElemStore(oldLibStore); a++)
	{
		gkLibraryOld oldLib;
		gkLibrary lib;
		
		getIndexStore(oldLibStore, a, &oldLib);
		
		lib.libraryUID = oldLib.libraryUID;
		lib.mean = oldLib.mean;
		lib.stddev = oldLib.stddev;
		lib.spareUN3 = oldLib.spareUN3;
		lib.spareUN2 = oldLib.spareUN2;
		lib.spareUN1 = oldLib.spareUN1;
		lib.spareUTG = oldLib.spareUTG;
		lib.forceBOGunitigger = oldLib.forceBOGunitigger;
		lib.isNotRandom = oldLib.isNotRandom;
		lib.spareALN = oldLib.spareALN;
		lib.doNotTrustHomopolymerRuns = oldLib.doNotTrustHomopolymerRuns;
		lib.spareOBT = oldLib.spareOBT;
		lib.doTrim_initialNone = oldLib.doTrim_initialNone;
		lib.doTrim_initialMerBased = oldLib.doTrim_initialMerBased;
		lib.doTrim_initialFlowBased = oldLib.doTrim_initialFlowBased;
		lib.doTrim_initialQualityBased = oldLib.doTrim_initialQualityBased;
		lib.doRemoveDuplicateReads = oldLib.doRemoveDuplicateReads;
		lib.doTrim_finalLargestCovered = oldLib.doTrim_finalLargestCovered;
		lib.doTrim_finalEvidenceBased = oldLib.doTrim_finalEvidenceBased;
		lib.doRemoveSpurReads = oldLib.doRemoveSpurReads;
		lib.doRemoveChimericReads = oldLib.doRemoveChimericReads;
		lib.doConsensusCorrection = oldLib.doConsensusCorrection;
		lib.spareGKP = oldLib.spareGKP;
		lib.UNUSEDusePackedFragments = oldLib.UNUSEDusePackedFragments;
		lib.spareLIB = oldLib.spareLIB;
		lib.orientation = oldLib.orientation;
		
		lookupLibraryName(a, &lib);
		
		appendIndexStore(libStore, &lib);
		
		fprintf(stdout, "Appended v8 library:\nlibraryIID="F_U32"\nlibraryName=%s\n", a, lib.libraryName);
	}
	
	closeStore(oldLibStore);
	closeStore(libStore);
	
	fprintf(stdout, "Library store upgraded to v8:\n%s\n", libStorePath);
}

void getDumpMap(char* dumpPath)
{
	errno = 0;
	
	FILE* dumpFile = fopen(dumpPath, "r");
	
	if (errno != 0)
	{
		fprintf(stderr, "Unable to open dump file for reading:\nfile=%s\nerror=%s\n", dumpPath, strerror(errno));
		exit(1);
	}
	
	char* dumpLine;
	bool pastHeader = false;
	
	do
	{
		dumpLine = new char[10240];
		
		fgets(dumpLine, 256, dumpFile);
		
		if (dumpLine != NULL)
		{
			if (pastHeader)
			{
				if (strlen(dumpLine) > 0)
				{
					char* name = strtok(dumpLine, "\t");
					AS_IID iid = atoi(strtok(NULL, "\t"));
					
					dumpMap[iid] = name;
				}
			}
			else
			{
				pastHeader = true;
			}
		}
	}
	while ((dumpLine != NULL) && (strlen(dumpLine) > 0));
	
	fclose(dumpFile);
	
	fprintf(stdout, "Found %u library item[s] in dump.\n", dumpMap.size());
}

char* renameStore(char* storePath)
{
	char* oldStorePath = strdup(storePath);
	
	strcat(oldStorePath, STORE_RENAME_SUFFIX);
	
	if (!AS_UTL_fileExists(oldStorePath, 0, 1))
	{
		fprintf(stdout, "Renaming store:\n%s\n%s\n", storePath, oldStorePath);
	
		rename(storePath, oldStorePath);
	}
	else
	{
		fprintf(stdout, "Using existing old store backup:\n%s\n", oldStorePath);
	}
	
	return oldStorePath;
}

void printUsage(char* executableName)
{
	fprintf(stdout, "Usage: %s -s <store path> -d <dump path>\n", executableName);
	fprintf(stdout, "\n");
	fprintf(stdout, "-s  path to the gatekeeper store to use\n");
	fprintf(stdout, "-d  path to the tabular gatekeeper dump\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "Upgrades a library store to contain a library name.\n");
	fprintf(stdout, "Note: this is the main change between v7 and v8 of the store implementation.\n");
}

char* getSubStorePath(char* storePath, char* subStoreFileName)
{
	char* subStorePath = new char[strlen(storePath) + strlen(subStoreFileName) + 1];
	
	strcpy(subStorePath, storePath);
	strcat(subStorePath, "/");
	strcat(subStorePath, subStoreFileName);
	
	return subStorePath;
}

int main(int argc, char** argv)
{
	if (argc > 1)
	{
		char* storePath = NULL, *dumpPath = NULL;
		int arg = 1;
	
		while (arg < argc)
		{
			if (strcmp(argv[arg], "-s") == 0)
			{
				storePath = argv[++arg];
			}
			else if (strcmp(argv[arg], "-d") == 0)
			{
				dumpPath = argv[++arg];
			}
			else
			{
				fprintf(stderr, "Unknown argument: %s\n", argv[arg]);
				exit(1);
			}
	
			arg++;
		}
		
		if (storePath == NULL)
		{
			fprintf(stderr, "A store path must be provided.\n");
			exit(1);
		}
		
		if (!AS_UTL_fileExists(storePath, 1, 1))
		{
			fprintf(stderr, "Store path does not exist: %s\n", storePath);
			exit(1);
		}
		
		if (dumpPath == NULL)
		{
			fprintf(stderr, "A dump path must be provided.\n");
			exit(1);
		}
		
		if (!AS_UTL_fileExists(dumpPath, 0, 1))
		{
			fprintf(stderr, "Dump path does not exist: %s\n", dumpPath);
			exit(1);
		}
		
		char* libStorePath = getSubStorePath(storePath, LIB_STORE_FILENAME), 
				*infoStorePath = getSubStorePath(storePath, INF_STORE_FILENAME), 
				*oldLibStorePath = renameStore(libStorePath), 
				*oldInfoStorePath = renameStore(infoStorePath);
		
		getDumpMap(dumpPath);
		
		processLibStore(oldLibStorePath, libStorePath);
		processInfoStore(oldInfoStorePath, infoStorePath);
	}
	else
	{
		printUsage(argv[0]);
	}
	
	exit(0);
}
