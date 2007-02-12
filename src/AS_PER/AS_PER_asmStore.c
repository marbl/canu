
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
static char CM_ID[] = "$Id: AS_PER_asmStore.c,v 1.5 2007-02-12 22:16:58 brianwalenz Exp $";

/*************************************************************************
 Module:  AS_PER_asmStore
 Description:
    A thin layer on top of the IndexStore supporing the storage and
 retrieval of records in assembly output file.
    The idea is to provide easier to use shortcuts for the common
 operations, and let the other operations be accessed through the
 generic Index Store API.

 Assumptions:
    Nothing special beyond genericStore.rtf

 Document:
      GenericStore.rtf

 *************************************************************************/

#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_asmStore.h"

char * ASM_Filenames[NUM_ASM_FILES] =
{
  "asm.mdi",
  "asm.bkt",
  
  "asm.afg",
  "asm.lkg",
  "asm.aci",
  "asm.asi"
  
  "asm.utg",
  "asm.utf",
  "asm.uci",
  "asm.usi",
  
  "asm.cco",
  "asm.ccf",
  "asm.ccu",
  
  "asm.dsc",
  
  "asm.scf",
  "asm.scg",
  "asm.scc",
  
  "asm.phash"
};

char * MAP_Filenames[NUM_MAP_FILES] =
{
  "map.chr",

  "map.fin",
  "map.phash",
};

void MakeDirectoryAsNeeded(char * path)
{
  DIR *dbDir;
  dbDir = opendir(path);
  if(dbDir == NULL)
  {
    char command[1024];
    int sysret;
    
    sprintf(command, "mkdir %s", path);
    sysret = system(command);
    assert(sysret == 0);
  }
  else
    closedir(dbDir);
}

int testOpenFile(char * path, char * name, char * mode)
{
  char fullName[FILENAME_MAX];
  FILE * fp;

  sprintf(fullName, "%s/%s", path, name);
  
  fp = fopen(fullName, mode);
  if(fp != NULL)
  {
    fclose(fp);
    return 1;
  }
  return 0;
}
int TestOpenAssemblyStore(AssemblyStore *asmStore)
{
  char * mode = "r+";
  int exists = 0;
  DIR *dbDir;
  int i;
  
  fprintf(stderr,"*** TestOpen %s\n", asmStore->storePath);

  dbDir = opendir(asmStore->storePath);

  if(dbDir != NULL)
  {
    int fileCount = 0;
    fprintf (stderr, "*** Directory exists %s... \n", asmStore->storePath);

    exists = -1;
    closedir (dbDir);

    for(i = 0; i < NUM_ASM_FILES; i++)
      fileCount += testOpenFile(asmStore->storePath, ASM_Filenames[i], mode);
    
    if(fileCount == NUM_ASM_FILES)
    {
      fprintf(stderr,"*  All files exist\n");
      exists = 1;
    }
    else
    {
      fprintf(stderr,"*  Directory exists -- %d files missing\n",
              NUM_ASM_FILES - fileCount);
    }
  }
  else
  {
    fprintf (stderr,
             "*** Directory DOES NOT exist %s... \n", asmStore->storePath);
  }
  return exists;
}


void removeFile(char * path, char * name)
{
  char command[FILENAME_MAX];
  int sysret;
  
  sprintf(command,"rm -f %s/%s", path, name);
  sysret = system(command);
  assert(sysret == 0);
}
int RemoveAssemblyStoreFiles(AssemblyStore *asmStore)
{
  int i;
  
  fprintf(stderr,"*** Remove %s\n", asmStore->storePath);
  for(i = 0; i < NUM_ASM_FILES; i++)
    removeFile(asmStore->storePath, ASM_Filenames[i]);
  
  return 0;
}
int RemoveMapStoreFiles(MapStore * mapStore)
{
  int i;
  
  fprintf(stderr,"*** Remove %s\n", mapStore->storePath);
  for(i = 0; i < NUM_MAP_FILES; i++)
    removeFile(mapStore->storePath, MAP_Filenames[i]);
  
  return 0;
}


void copyFile(char * filename, char * fromDir, char * toDir)
{
  char buffer[FILENAME_MAX];
  int sysret;
  
  sprintf(buffer,"cp %s/%s %s", fromDir, filename, toDir);
  sysret = system(buffer);
  assert(sysret == 0);
}
int CopyAssemblyStoreFiles(AssemblyStore *asmStore, char *path)
{
  int i;

  fprintf(stderr, "*** Copy %s//%s == > %s\n",
          getcwd(NULL,256), asmStore->storePath, path);

  MakeDirectoryAsNeeded(path);
  
  for(i = 0; i < NUM_ASM_FILES; i++)
    copyFile(ASM_Filenames[i], asmStore->storePath, path);
  
  return(0);
}
int CopyMapStoreFiles(MapStore *mapStore, char *path)
{
  int i;

  fprintf(stderr, "*** Copy %s//%s == > %s\n",
          getcwd(NULL,256), mapStore->storePath, path);

  MakeDirectoryAsNeeded(path);
  
  for(i = 0; i < NUM_MAP_FILES; i++)
    copyFile(MAP_Filenames[i], mapStore->storePath, path);
  
  return(0);
}


AssemblyStore * OpenAssemblyStoreCommon(char * path, char *mode)
{
  char name[FILENAME_MAX];
  AssemblyStore * asmStore;

  asmStore = (AssemblyStore *) calloc(1, sizeof(AssemblyStore));
  assert(asmStore != NULL);
  strcpy(asmStore->storePath, path);

  fprintf(stderr, "*** Open %s//%s\n", getcwd(NULL,256), asmStore->storePath);

  if(!strcmp(mode, "r"))
     MakeDirectoryAsNeeded(path);
     
  sprintf(name, "%s/asm.mdi", asmStore->storePath);
  asmStore->mdiStore = openASM_MDIStore(name, mode);
  sprintf(name, "%s/asm.bkt", asmStore->storePath);
  asmStore->bktStore = openASM_BucketStore(name, mode);
  
  sprintf(name, "%s/asm.afg", asmStore->storePath);
  asmStore->afgStore = openASM_AFGStore(name, mode);
  sprintf(name, "%s/asm.lkg", asmStore->storePath);
  asmStore->lkgStore = openASM_LKGStore(name, mode);
  sprintf(name, "%s/asm.aci", asmStore->storePath);
  asmStore->aciStore = openASM_InstanceStore(name, mode);
  sprintf(name, "%s/asm.asi", asmStore->storePath);
  asmStore->asiStore = openASM_InstanceStore(name, mode);
  
  sprintf(name, "%s/asm.utg", asmStore->storePath);
  asmStore->utgStore = openASM_UTGStore(name, mode);
  sprintf(name, "%s/asm.utf", asmStore->storePath);
  asmStore->utfStore = openASM_IIDStore(name, mode);
  sprintf(name, "%s/asm.uci", asmStore->storePath);
  asmStore->uciStore = openASM_InstanceStore(name, mode);
  sprintf(name, "%s/asm.usi", asmStore->storePath);
  asmStore->usiStore = openASM_InstanceStore(name, mode);
  
  sprintf(name, "%s/asm.cco", asmStore->storePath);
  asmStore->ccoStore = openASM_CCOStore(name, mode);
  sprintf(name, "%s/asm.ccf", asmStore->storePath);
  asmStore->ccfStore = openASM_IIDStore(name, mode);
  sprintf(name, "%s/asm.ccu", asmStore->storePath);
  asmStore->ccuStore = openASM_IIDStore(name, mode);
  
  sprintf(name, "%s/asm.dsc", asmStore->storePath);
  asmStore->dscStore = openASM_DSCStore(name, mode);
  
  sprintf(name, "%s/asm.scf", asmStore->storePath);
  asmStore->scfStore = openASM_SCFStore(name, mode);
  sprintf(name, "%s/asm.scg", asmStore->storePath);
  asmStore->scgStore = openASM_GapStore(name, mode);
  sprintf(name, "%s/asm.scc", asmStore->storePath);
  asmStore->sccStore = openASM_IIDStore(name, mode);
  
  if(NULLSTOREHANDLE == asmStore->mdiStore ||
     NULLSTOREHANDLE == asmStore->bktStore ||
     NULLSTOREHANDLE == asmStore->afgStore ||
     NULLSTOREHANDLE == asmStore->lkgStore ||
     NULLSTOREHANDLE == asmStore->aciStore ||
     NULLSTOREHANDLE == asmStore->asiStore ||
     NULLSTOREHANDLE == asmStore->utgStore ||
     NULLSTOREHANDLE == asmStore->utfStore ||
     NULLSTOREHANDLE == asmStore->uciStore ||
     NULLSTOREHANDLE == asmStore->usiStore ||
     NULLSTOREHANDLE == asmStore->ccoStore ||
     NULLSTOREHANDLE == asmStore->ccfStore ||
     NULLSTOREHANDLE == asmStore->ccuStore ||
     NULLSTOREHANDLE == asmStore->dscStore ||
     NULLSTOREHANDLE == asmStore->scfStore ||
     NULLSTOREHANDLE == asmStore->scgStore ||
     NULLSTOREHANDLE == asmStore->sccStore)
  {
    fprintf(stderr,"**** Failure to open Assembly Store ...\n");
    return NULL;
  }

  sprintf(name,"%s/asm.phash", asmStore->storePath);
  if(mode && *mode == 'r' && *(mode + 1) == '\0')
  {
    asmStore->hashTable = OpenReadOnlyPHashTable_AS(name);
  }
  else
  {
    asmStore->hashTable = OpenPHashTable_AS(name);
  }
  if(asmStore->hashTable == NULL)
  {
    fprintf(stderr,"**** Failed to open Assembly Persistent HashTable...\n");
    return NULL;
  }

  asmStore->gkpStore = NULL;

  return asmStore;
}
MapStore * OpenMapStoreCommon(char * path, char *mode)
{
  char name[FILENAME_MAX];
  MapStore * mapStore;

  mapStore = (MapStore *) calloc(1, sizeof(MapStore));
  assert(mapStore != NULL);
  strcpy(mapStore->storePath, path);

  fprintf(stderr, "*** Open %s//%s\n", getcwd(NULL,256), mapStore->storePath);

  if(!strcmp(mode, "r"))
     MakeDirectoryAsNeeded(path);

  sprintf(name, "%s/map.chr", mapStore->storePath);
  mapStore->chrStore = openASM_CHRStore(name, mode);
  
#ifdef NEVER
  sprintf(name, "%s/map.cfm", mapStore->storePath);
  mapStore->cfmStore = openASM_MemberStore(name, mode);
#endif
  
  sprintf(name, "%s/map.fin", mapStore->storePath);
  mapStore->finStore = openASM_InstanceStore(name, mode);
  
  if(NULLSTOREHANDLE == mapStore->chrStore ||
#ifdef NEVER
     NULLSTOREHANDLE == mapStore->cfmStore ||
#endif
     NULLSTOREHANDLE == mapStore->finStore)
  {
    fprintf(stderr,"**** Failure to open Map Store ...\n");
    return NULL;
  }

  sprintf(name,"%s/map.phash", mapStore->storePath);
  if(mode && *mode == 'r' && *(mode + 1) == '\0')
  {
    mapStore->hashTable = OpenReadOnlyPHashTable_AS(name);
  }
  else
  {
    mapStore->hashTable = OpenPHashTable_AS(name);
  }
  if(mapStore->hashTable == NULL)
  {
    fprintf(stderr,"**** Failed to open Map Store Persistent HashTable...\n");
    return NULL;
  }
  
  return mapStore;
}


AssemblyStore * OpenAssemblyStore(char * path)
{
  return OpenAssemblyStoreCommon(path, "r+");
}
AssemblyStore * OpenReadOnlyAssemblyStore(char * path)
{
  return OpenAssemblyStoreCommon(path, "r");
}
MapStore * OpenMapStore(char * path)
{
  return OpenMapStoreCommon(path, "r+");
}
MapStore * OpenReadOnlyMapStore(char * path)
{
  return OpenMapStoreCommon(path, "r");
}


int OpenGateKeeperStoreAssemblyStore(AssemblyStore * asmStore,
                                     char * gkpStorePath)
{
  asmStore->gkpStore = openGateKeeperStore(gkpStorePath, FALSE);
  assert(asmStore->gkpStore != NULL);
  return(1);
}


int OpenFragmentStoreAssemblyStore(AssemblyStore * asmStore,
                                   char * frgStorePath)
{
  return(1);
}


AssemblyStore * CreateAssemblyStore(char * path,
                                    char * gkpStorePath,
                                    char * frgStorePath)
{
  AssemblyStore * asmStore;
  char name[FILENAME_MAX];

  asmStore = (AssemblyStore *) calloc(1, sizeof(AssemblyStore));
  assert(asmStore != NULL);
  strcpy(asmStore->storePath, path);

  fprintf(stderr,"*** Create store %s at cwd %s\n",
          asmStore->storePath, getcwd(NULL, 256));

  MakeDirectoryAsNeeded(path);
  
  sprintf(name,"%s/asm.mdi", asmStore->storePath);
  asmStore->mdiStore = createASM_MDIStore(name, "mdi",1);
  sprintf(name,"%s/asm.bkt", asmStore->storePath);
  asmStore->bktStore = createASM_BucketStore(name, "bkt",1);
  
  sprintf(name,"%s/asm.afg", asmStore->storePath);
  asmStore->afgStore = createASM_AFGStore(name, "afg",1);
  sprintf(name,"%s/asm.lkg", asmStore->storePath);
  asmStore->lkgStore = createASM_LKGStore(name, "lnk",1);
  sprintf(name,"%s/asm.aci", asmStore->storePath);
  asmStore->aciStore = createASM_InstanceStore(name, "aci",1);
  sprintf(name,"%s/asm.asi", asmStore->storePath);
  asmStore->asiStore = createASM_InstanceStore(name, "asi",1);

  sprintf(name,"%s/asm.utg", asmStore->storePath);
  asmStore->utgStore = createASM_UTGStore(name, "utg",1);
  sprintf(name,"%s/asm.utf", asmStore->storePath);
  asmStore->utfStore = createASM_IIDStore(name, "utf",1);
  sprintf(name,"%s/asm.uci", asmStore->storePath);
  asmStore->uciStore = createASM_InstanceStore(name, "uci",1);
  sprintf(name,"%s/asm.usi", asmStore->storePath);
  asmStore->usiStore = createASM_InstanceStore(name, "usi",1);
  
  sprintf(name,"%s/asm.cco", asmStore->storePath);
  asmStore->ccoStore = createASM_CCOStore(name, "cco",1);
  sprintf(name,"%s/asm.ccf", asmStore->storePath);
  asmStore->ccfStore = createASM_IIDStore(name, "ccf",1);
  sprintf(name,"%s/asm.ccu", asmStore->storePath);
  asmStore->ccuStore = createASM_IIDStore(name, "ccu",1);

  sprintf(name,"%s/asm.dsc", asmStore->storePath);
  asmStore->dscStore = createASM_DSCStore(name, "dsc",1);

  sprintf(name,"%s/asm.scf", asmStore->storePath);
  asmStore->scfStore = createASM_SCFStore(name, "scf",1);
  sprintf(name,"%s/asm.scg", asmStore->storePath);
  asmStore->scgStore = createASM_GapStore(name, "scg",1);
  sprintf(name,"%s/asm.scc", asmStore->storePath);
  asmStore->sccStore = createASM_IIDStore(name, "scc",1);
  
  sprintf(name,"%s/asm.phash", asmStore->storePath);
  asmStore->hashTable = CreatePHashTable_AS(1024,name);

  if(gkpStorePath != NULL)
    OpenGateKeeperStoreAssemblyStore(asmStore, gkpStorePath);
  else
    asmStore->gkpStore = NULL;
  
  return asmStore;
}
MapStore * CreateMapStore(char * path)
{
  MapStore * mapStore;
  char name[FILENAME_MAX];

  mapStore = (MapStore *) calloc(1, sizeof(MapStore));
  assert(mapStore != NULL);
  strcpy(mapStore->storePath, path);

  fprintf(stderr,"*** Create store %s at cwd %s\n",
          mapStore->storePath, getcwd(NULL, 256));

  MakeDirectoryAsNeeded(path);
  
  sprintf(name,"%s/map.chr", mapStore->storePath);
  mapStore->chrStore = createASM_CHRStore(name, "chr",1);

#ifdef NEVER
  sprintf(name,"%s/map.cfm", mapStore->storePath);
  mapStore->cfmStore = createASM_MemberStore(name, "cfm",1);
#endif
  
  sprintf(name,"%s/map.fin", mapStore->storePath);
  mapStore->finStore = createASM_InstanceStore(name, "fin",1);
  
  sprintf(name,"%s/map.phash", mapStore->storePath);
  mapStore->hashTable = CreatePHashTable_AS(1024,name);

  return mapStore;
}


void CloseAssemblyStore(AssemblyStore *asmStore)
{
  fprintf(stderr,"*** Close directory %s\n", asmStore->storePath);

  if(asmStore->gkpStore != NULL)
  {
    closeGateKeeperStore(asmStore->gkpStore);
    free(asmStore->gkpStore);
  }
  
  if(asmStore->mdiStore != NULLSTOREHANDLE)
    closeStore(asmStore->mdiStore);
  if(asmStore->bktStore != NULLSTOREHANDLE)
    closeStore(asmStore->bktStore);
  
  if(asmStore->afgStore != NULLSTOREHANDLE)
    closeStore(asmStore->afgStore);
  if(asmStore->lkgStore != NULLSTOREHANDLE)
    closeStore(asmStore->lkgStore);
  if(asmStore->aciStore != NULLSTOREHANDLE)
    closeStore(asmStore->aciStore);
  if(asmStore->asiStore != NULLSTOREHANDLE)
    closeStore(asmStore->asiStore);
  
  if(asmStore->utgStore != NULLSTOREHANDLE)
    closeStore(asmStore->utgStore);
  if(asmStore->utfStore != NULLSTOREHANDLE)
    closeStore(asmStore->utfStore);
  if(asmStore->uciStore != NULLSTOREHANDLE)
    closeStore(asmStore->uciStore);
  if(asmStore->usiStore != NULLSTOREHANDLE)
    closeStore(asmStore->usiStore);
  
  if(asmStore->ccoStore != NULLSTOREHANDLE)
    closeStore(asmStore->ccoStore);
  if(asmStore->ccfStore != NULLSTOREHANDLE)
    closeStore(asmStore->ccfStore);
  if(asmStore->ccuStore != NULLSTOREHANDLE)
    closeStore(asmStore->ccuStore);
  
  if(asmStore->dscStore != NULLSTOREHANDLE)
    closeStore(asmStore->dscStore);
  
  if(asmStore->scfStore != NULLSTOREHANDLE)
    closeStore(asmStore->scfStore); 
  if(asmStore->scgStore != NULLSTOREHANDLE)
    closeStore(asmStore->scgStore); 
  if(asmStore->sccStore != NULLSTOREHANDLE)
    closeStore(asmStore->sccStore);
  
  if(asmStore->hashTable != NULL)
    ClosePHashTable_AS(asmStore->hashTable);
}
void CloseMapStore(MapStore *mapStore)
{
  fprintf(stderr,"*** Close directory %s\n", mapStore->storePath);

  if(mapStore->chrStore != NULLSTOREHANDLE)
    closeStore(mapStore->chrStore);

#ifdef NEVER
  if(mapStore->cfmStore != NULLSTOREHANDLE)
    closeStore(mapStore->cfmStore);
#endif

  if(mapStore->finStore != NULLSTOREHANDLE)
    closeStore(mapStore->finStore);
  
  if(mapStore->hashTable != NULL)
    ClosePHashTable_AS(mapStore->hashTable);
}
