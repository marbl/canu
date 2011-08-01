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

static const char *rcsid = "$Id: upgrade-v7-to-v8.C,v 1.4 2011-08-01 20:33:30 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <map>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_Hash.h"
#include "StringUtils.h"
#include "TestUtils.h"

#define INF_STORE_FILENAME "inf"
#define LIB_STORE_FILENAME "lib"
#define U2I_STORE_FILENAME "u2i"
#define UID_STORE_FILENAME "uid"
#define STORE_RENAME_SUFFIX ".v7-old"

#define AS_IID_LIB 3

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

void processInfoStore(char* oldInfoStorePath, char* infoStorePath)
{
	gkStoreInfo info;
	
	errno = 0;
	
	FILE* oldInfoFile = fopen(oldInfoStorePath, "r");
	
	if (AS_UTL_safeRead(oldInfoFile, &info, "inf", sizeof(info), 1))
	{
		assertTrue(info.gkVersion == 7, (string("Store is not version 7: version=") + StringUtils::toString(info.gkVersion)).c_str());
		
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

void processLibStore(map<AS_IID, const char*>& idMap, StoreStruct& oldLibStore, StoreStruct& libStore)
{
	size oldLibSize = sizeof(gkLibraryOld);
	
	for (map<AS_IID, const char*>::iterator iterator = idMap.begin(); iterator != idMap.end(); iterator++)
	{
		printf(F_U32"="F_STR"\n", (*iterator).first, (*iterator).second);
	}
	
	for (AS_IID a = 1; a <= idMap.size(); a++)
	{
		gkLibraryOld oldLib;
		gkLibrary lib;
		
		getIndexStore(&oldLibStore, a, &oldLib);
		
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
		
		assertTrue(idMap.count(a) > 0, StringUtils::concat(2, "Library does not contain an UID string: iid=", StringUtils::toString(a)));
		
		strcpy(lib.libraryName, idMap[a]);
		
		appendIndexStore(&libStore, &lib);
		
		fprintf(stdout, "Appended v8 library:\nlibraryIID="F_U32"\nlibraryName=%s\n", a, lib.libraryName);
	}
	
	closeStore(&oldLibStore);
	closeStore(&libStore);
	
	fprintf(stdout, "Library store upgraded to v8.\n");
}

void getIdMap(map<AS_IID, const char*>& idMap, char* uidStorePath, char* u2iStorePath)
{	
	HashTable_AS* uidToIidTable = LoadUIDtoIIDHashTable_AS(u2iStorePath);
	StoreStruct* uidStore = convertStoreToMemoryStore(openStore(uidStorePath, "r"));
	HashTable_Iterator_AS idIterator;
	
	InitializeHashTable_Iterator_AS(uidToIidTable, &idIterator);
	
	uint64 key = 0, value = 0;
	uint32 valueType = 0;
	
	while (NextHashTable_Iterator_AS(&idIterator, &key, &value, &valueType))
	{
		if (valueType == AS_IID_LIB)
		{
			uint32 actualLength = 0;
			int64 offset = 0;
			
			char* strTemp = getStringStorePtr(uidStore, value, &actualLength, &offset);
			
			assertTrue(actualLength > 0, StringUtils::concat(2, "UID string must not be empty: iid=", StringUtils::toString(value)));
			
			char* str = new char[actualLength];
			strcpy(str, strTemp);
			
			idMap[value] = str;
		}
	}
	
	closeStore(uidStore);
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
		char* storePath = NULL;
		int arg = 1;
	
		while (arg < argc)
		{
			if (strcmp(argv[arg], "-s") == 0)
			{
				storePath = argv[++arg];
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
		
		char* libStorePath = getSubStorePath(storePath, LIB_STORE_FILENAME), 
				*infoStorePath = getSubStorePath(storePath, INF_STORE_FILENAME), 
				*oldLibStorePath = renameStore(libStorePath), 
				*oldInfoStorePath = renameStore(infoStorePath);
		
		map<AS_IID, const char*> idMap;
		
		getIdMap(idMap, getSubStorePath(storePath, UID_STORE_FILENAME), 
			getSubStorePath(storePath, U2I_STORE_FILENAME));
		
		StoreStruct* oldLibStore = convertStoreToMemoryStore(openStore(oldLibStorePath, "r")),
			*libStore = createIndexStore(libStorePath, "lib", sizeof(gkLibrary), 1);
		
		processInfoStore(oldInfoStorePath, infoStorePath);
		processLibStore(idMap, *oldLibStore, *libStore);
	}
	else
	{
		printUsage(argv[0]);
	}
	
	exit(0);
}
