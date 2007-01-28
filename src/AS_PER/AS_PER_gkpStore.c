
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
static char CM_ID[] = "$Id: AS_PER_gkpStore.c,v 1.10 2007-01-28 21:52:25 brianwalenz Exp $";

/*************************************************************************
 Module:  AS_PER_gkpfrgStore
 Description:
    A thin layer on top of the IndexStore supporing the storage and
 retrieval of records used by the gatekeeper records.
    The idea is to provide easier to use shortcuts for the common
 operations, and let the other operations be accessed through the
 generic Index Store API.

 Assumptions:
    Nothing special beyond genericStore.rtf

 Document:
      GenericStore.rtf

 *************************************************************************/

/* RCS Info
 * $Id: AS_PER_gkpStore.c,v 1.10 2007-01-28 21:52:25 brianwalenz Exp $
 * $Revision: 1.10 $
 *
 */
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
#include "AS_PER_gkpStore.h"

int CreateGateKeeperLinkRecordIterator(GateKeeperLinkStore store, uint32 startFromLink, 
				       uint32 followFrag, GateKeeperLinkRecordIterator *iterator){

  assert(startFromLink);

  if (startFromLink == 0)
    return(1);

  iterator->store          = store;
  iterator->prevLinkRecord = 0;
  iterator->linkRecord     = startFromLink;
  iterator->followFrag     = followFrag;

  return(0);
}

#if 0
//  UNUSED
int CreateGateKeeperLinkRecordFromFragmentIterator(GateKeeperLinkStore store, 
						   uint32 followFrag, GateKeeperLinkRecordIterator *iterator){
  GateKeeperFragmentRecord gkFrag;

  getGateKeeperFragmentStore(store, followFrag, &gkFrag);

 return CreateGateKeeperLinkRecordIterator(store, gkFrag.linkHead, followFrag, iterator);


}
#endif




int NextGateKeeperLinkRecordIterator(GateKeeperLinkRecordIterator *iterator, 
				     GateKeeperLinkRecord *link){
  int tprev = iterator->prevLinkRecord;

  if(iterator->linkRecord == 0){
#ifdef DEBUG
    fprintf(stderr,"*** Iterator bailing \n");
#endif
    return 0;
  }
#ifdef DEBUG
  fprintf(stderr,"*NextGateKeeperLinkRecord getting linkRecord %ld (%ld)\n", 
	  iterator->linkRecord, getLastElemStore(iterator->store));
#endif

  {
    int ret = getGateKeeperLinkStore(iterator->store, iterator->linkRecord, link);
    assert(ret == 0);
  }
  if(!(iterator->followFrag == link->frag1 || iterator->followFrag == link->frag2)){



  }

#ifdef DEBUG
    fprintf(stderr,"*** Got link (%d,%d) \n", link->frag1, link->frag2);
#endif

  iterator->prevLinkRecord = iterator->linkRecord;
  if(iterator->followFrag == link->frag1)
    iterator->linkRecord = link->frag1Next;
  else if(iterator->followFrag == link->frag2)
    iterator->linkRecord = link->frag2Next;
  else{
    fprintf(stderr,"* NextGateKeeperLinkRecordIterator Internal error!  iterator->followFrag = %d\n",
	    iterator->followFrag);
    fprintf(stderr,"* link is %d  prev = %d frag1 = %d frag2 = %d\n",
	    iterator->prevLinkRecord, tprev, link->frag1, link->frag2);
    assert(0);
  }

  return 1;

}



/***********************************************************************************/

/* Returns link index of link sought */
int findLink(GateKeeperLinkStore store, 
             uint32 frag,
             uint32 linkHead, 
             GateKeeperLinkRecord *searchlink,
             GateKeeperLinkRecord *foundlink){
       
  uint32 frag1IID     = searchlink->frag1;
  uint32 frag2IID     = searchlink->frag2;
  int linktype        = searchlink->type;
  uint64 distance     = searchlink->distance;
  int linkOrientation = searchlink->orientation;

  GateKeeperLinkRecordIterator iterator;
  GateKeeperLinkRecord link;

  if(linkHead == NULL_LINK)
    return NULL_LINK;

  CreateGateKeeperLinkRecordIterator(store, linkHead,
                                     frag, &iterator);

  while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
    if((link.deleted) ||
       (link.frag1        != frag1IID) ||
       (link.frag2        != frag2IID) ||
       (link.type         != linktype) || 
       ((linkOrientation  != AS_GKP_UNKNOWN) && (link.orientation != linkOrientation)) ||
       ((distance         != 0)              && (link.distance    != distance))){
      continue;
    }
#ifdef DEBUG_GKP
    fprintf(stderr,"* Found link (%d,%d) %d\n",
            frag1IID, frag2IID, link.type);
#endif
    if(foundlink)
      *foundlink = link;

    return iterator.prevLinkRecord;
  }
  return NULL_LINK;
}

				



/***********************************************************************************/
int unlinkLink_GKP(GateKeeperLinkStore gkplStore, 
		 GateKeeperFragmentStore     gkpStore, 
		 uint32 frag1,
		 uint32 frag2,
 	         GateKeeperFragmentRecord *gkf1, 
		 GateKeeperFragmentRecord *gkf2,
		 GateKeeperLinkRecord *newLink,
		 int deleteLinkIndex){

  GateKeeperLinkRecord link;
  int found = 0;

#ifdef DEBUG_GKP
 fprintf(stderr,"* unlink link %d (%d,%d) linkheads (%d)%d and (%d)%d\n",
	 deleteLinkIndex, newLink->frag1, newLink->frag2, frag1,gkf1->linkHead, frag2, 
	 gkf2->linkHead);
 verifyLink_GKP(gkplStore, gkpStore, deleteLinkIndex);
#endif

  if(gkf1->linkHead == deleteLinkIndex){
#ifdef DEBUG_GKP
    fprintf(stderr,"* Popped deleteLink %d from linkHead1\n", deleteLinkIndex);
#endif
    found = TRUE;
	   if(frag1 == newLink->frag1){
	       gkf1->linkHead = newLink->frag1Next;
	   }else{
	     assert(frag1 == newLink->frag2);
	     gkf1->linkHead = newLink->frag2Next;
	   }
	gkf1->numLinks--;
	setGateKeeperFragmentStore(gkpStore,newLink->frag1, gkf1);

  }else{
       GateKeeperLinkRecordIterator iterator;

       CreateGateKeeperLinkRecordIterator(gkplStore, gkf1->linkHead,
					  newLink->frag1, &iterator);

       while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
#ifdef DEBUG_GKP
	 fprintf(stderr,"* looking at link %d (%d,%d) next is %d\n",
		 iterator.prevLinkRecord, link.frag1, link.frag2, iterator.linkRecord);
#endif
	 if(iterator.linkRecord == deleteLinkIndex){ /* The next fetch gets us the deleteLink*/
	   found = TRUE;
	   if(link.frag1 == newLink->frag1){
	       link.frag1Next = newLink->frag1Next;
	   }else{
	     assert(link.frag2 == newLink->frag1);
	     link.frag2Next = newLink->frag1Next;
	   }
	   /* update the link */
	   setGateKeeperLinkStore(gkplStore, iterator.prevLinkRecord, &link);
	   break;
	 }
       }
       
      if(!found){
	fprintf(stderr,"* Failed to find link %d (%d,%d) starting from frag %d linkhead %d\n",
		deleteLinkIndex,newLink->frag1, newLink->frag2, newLink->frag1,gkf1->linkHead);
	assert(0);
      }else{
	   gkf1->numLinks--;
	   setGateKeeperFragmentStore(gkpStore,newLink->frag1, gkf1);
      }
  }
  found = FALSE;
  if(gkf2->linkHead == deleteLinkIndex){
#ifdef DEBUG_GKP
    fprintf(stderr,"* Popped deleteLink %d from linkHead2\n", deleteLinkIndex);
#endif
    found = TRUE;
	   if(frag2 == newLink->frag1){
	       gkf2->linkHead = newLink->frag1Next;
	   }else{
	     assert(frag2 == newLink->frag2);
	     gkf2->linkHead = newLink->frag2Next;
	   }
	gkf2->numLinks--;
	setGateKeeperFragmentStore(gkpStore,newLink->frag2, gkf2);
  }else{
       GateKeeperLinkRecordIterator iterator;

       CreateGateKeeperLinkRecordIterator(gkplStore, gkf2->linkHead,
					  newLink->frag2, &iterator);

       while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
#ifdef DEBUG_GKP
	 fprintf(stderr,"* looking at link %d (%d,%d) next is %d\n",
		 iterator.prevLinkRecord, link.frag1, link.frag2, iterator.linkRecord);
#endif
	 if(iterator.linkRecord == deleteLinkIndex){ /* The next fetch gets us the deleteLink*/
	   found = TRUE;
	   if(link.frag1 == newLink->frag2){
	       link.frag1Next = newLink->frag2Next;
	   }else{
	       assert(link.frag2 == newLink->frag2);
	       link.frag2Next = newLink->frag2Next;
	   }
	   /* update the link */
	   setGateKeeperLinkStore(gkplStore, iterator.prevLinkRecord, &link);
	   break;
	 }
       }
       
      if(!found){
	fprintf(stderr,"* Failed to find link %d (%d,%d) starting from frag %d linkhead %d\n",
		deleteLinkIndex,newLink->frag1, newLink->frag2, newLink->frag2, gkf2->linkHead);
	assert(0);
      }else{
	gkf2->numLinks--;
	setGateKeeperFragmentStore(gkpStore,newLink->frag2, gkf2);
      }

  }

  /* update the link */
  newLink->frag1Next = 0;
  newLink->frag2Next = 0;
  newLink->deleted = 1;
  setGateKeeperLinkStore(gkplStore, deleteLinkIndex, newLink);
  return(0);

}

/***********************************************************************************/
int linkLink_GKP(GateKeeperLinkStore gkplStore, 
		 GateKeeperFragmentStore     gkpStore, 
		 GateKeeperLinkRecord *newLink,
		 uint32 frag1,
		 uint32 frag2,
 	         GateKeeperFragmentRecord *gkf1, 
		 GateKeeperFragmentRecord *gkf2){

  int   newLinkIndex = getLastElemStore(gkplStore) + 1;


  /* Insert at head of lists */
  newLink->frag1Next = newLink->frag2Next = 0;
  newLink->deleted = 0;

  assert(frag1 == newLink->frag1 &&
	 frag2 == newLink->frag2);

    newLink->frag1Next = gkf1->linkHead;
    newLink->frag2Next = gkf2->linkHead;
    gkf1->linkHead = newLinkIndex;
    gkf2->linkHead = newLinkIndex;
    gkf1->numLinks++;
    gkf2->numLinks++;
#ifdef DEBUG_GKP
 fprintf(stderr,"* linkLink_GKP newLink %d f1:%d(U%lu) gkf1Head = %d f2:%d(U%lu) gkf2Head = %d\n",
	 newLinkIndex, newLink->frag1, gkf1->readUID, gkf1->linkHead, newLink->frag2, gkf2->readUID, gkf2->linkHead);

#endif

  /* First, append the new link record to the gkplStore, and remember it's index */
  setGateKeeperFragmentStore(gkpStore,newLink->frag1, gkf1);
  setGateKeeperFragmentStore(gkpStore,newLink->frag2, gkf2);
  appendGateKeeperLinkStore(gkplStore, newLink);

#ifdef DEBUG_GKP
  verifyLink_GKP(gkplStore, gkpStore, newLinkIndex);
#endif
  
  return(0);
}

void InitGateKeeperStore(GateKeeperStore *gkpStore, const char *path){
  AssertPtr(gkpStore);
  strcpy(gkpStore->storePath, path);
  gkpStore->hashTable = NULL;
  gkpStore->batStore = (GateKeeperBatchStore)0;
  gkpStore->frgStore = (GateKeeperFragmentStore)0;
  gkpStore->lnkStore = (GateKeeperLinkStore)0;
  gkpStore->dstStore = (GateKeeperDistanceStore)0;
  gkpStore->s_dstStore = (GateKeeperDistanceStore)0;
  gkpStore->btgStore = (GateKeeperBactigStore)0;
  gkpStore->locStore = (GateKeeperLocaleStore)0;
  gkpStore->s_locStore = (GateKeeperLocaleStore)0;
  gkpStore->seqStore = (GateKeeperSequenceStore)0;
}

int TestOpenGateKeeperStoreCommon(GateKeeperStore *gkpStore,const char *mode){
  char name[FILENAME_MAX];
  int exists = 0;

  DIR *dbDir;
  FILE *fp;
  fprintf(stderr,"*** TestOpen %s\n", gkpStore->storePath);

  dbDir = opendir(gkpStore->storePath);

  if  (dbDir != NULL)
    {
      int fileCount = 0;
      int upgrade_count = 0;
      fprintf (stderr,
	       "*** Directory exists %s... \n", gkpStore->storePath);

      exists = -1;
      closedir (dbDir);
      sprintf(name,"%s/gkp.bat", gkpStore->storePath);
      fp = fopen(name,mode);
      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.frg", gkpStore->storePath);
      fp = fopen(name,mode);
      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.lnk", gkpStore->storePath);
      fp = fopen(name,mode);
      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.loc", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.s_loc", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.seq", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.btg", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.dst", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      sprintf(name,"%s/gkp.s_dst", gkpStore->storePath);
      fp = fopen(name,mode);

      if(fp){
	fileCount++;
	fclose(fp);
      }
      
      sprintf(name,"%s/gkp.phash", gkpStore->storePath);

      fp = fopen(name,mode);
      if(fp){
	fileCount++;
	fclose(fp);
      }
      if(fileCount + upgrade_count == NUM_GKP_FILES){
	fprintf(stderr,"*  All files exist\n");
	exists = 1;
      }else{
        if(fileCount + 4 == NUM_GKP_FILES){
          fprintf(stderr,"*  Minimum set of files exists\n");
          fprintf(stderr,"*  Upgrade needed for new files\n");
          exists = 1;
        }else{
          fprintf(stderr,"*  Directory exists -- %d files missing\n",
                  NUM_GKP_FILES - fileCount);
        }
      }	
    } else {
      fprintf (stderr,
	       "*** Directory DOES NOT exist %s... \n", gkpStore->storePath);
    }

  return exists;

}

int TestOpenGateKeeperStore(GateKeeperStore *gkpStore){
  return TestOpenGateKeeperStoreCommon(gkpStore,"r+");
}

int TestOpenReadOnlyGateKeeperStore(GateKeeperStore *gkpStore){
  return TestOpenGateKeeperStoreCommon(gkpStore,"r");
}

int RemoveGateKeeperStoreFiles(GateKeeperStore *gkpStore){
  char buffer[FILENAME_MAX];

  fprintf(stderr,"*** Remove %s\n", gkpStore->storePath);

  sprintf(buffer,"rm -f %s/gkp.bat", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.frg", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.lnk", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.loc", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.s_loc", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.seq", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.btg", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.dst", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.s_dst", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"rm -f %s/gkp.phash", gkpStore->storePath);
  if(system(buffer) != 0) assert(0);

  return 0;
}

int CopyGateKeeperStoreFiles(GateKeeperStore *gkpStore, char *path){
  char buffer[FILENAME_MAX];

  fprintf(stderr,"*** Copy %s//%s == > %s\n", getcwd(NULL,256), gkpStore->storePath, path);

  sprintf(buffer,"cp %s/gkp.bat %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.frg %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.lnk %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.dst %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.s_dst %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.seq %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.btg %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.loc %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.s_loc %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  sprintf(buffer,"cp %s/gkp.phash %s", gkpStore->storePath, path);
  if(system(buffer) != 0) assert(0);
  
  return(0);
}



int OpenGateKeeperStoreCommon(GateKeeperStore *gkpStore, char *mode){
  char name[FILENAME_MAX];

  fprintf(stderr,"*** Open %s//%s\n", getcwd(NULL,256), gkpStore->storePath);

     sprintf(name,"%s/gkp.bat", gkpStore->storePath);
     gkpStore->batStore = openGateKeeperBatchStore(name,mode); 
     sprintf(name,"%s/gkp.frg", gkpStore->storePath);
     gkpStore->frgStore = openGateKeeperFragmentStore(name,mode); 
     sprintf(name,"%s/gkp.lnk", gkpStore->storePath);
     gkpStore->lnkStore = openGateKeeperLinkStore(name,mode); 
     sprintf(name,"%s/gkp.loc", gkpStore->storePath);
     gkpStore->locStore = openGateKeeperLocaleStore(name,mode); 
     sprintf(name,"%s/gkp.s_loc", gkpStore->storePath);
     gkpStore->s_locStore = openGateKeeperLocaleStore(name,mode); 
     sprintf(name,"%s/gkp.seq", gkpStore->storePath);
     gkpStore->seqStore = openGateKeeperSequenceStore(name,mode); 
     sprintf(name,"%s/gkp.dst", gkpStore->storePath);
     gkpStore->dstStore = openGateKeeperDistanceStore(name,mode); 
     sprintf(name,"%s/gkp.s_dst", gkpStore->storePath);
     gkpStore->s_dstStore = openGateKeeperDistanceStore(name,mode); 
     sprintf(name,"%s/gkp.btg", gkpStore->storePath);
     gkpStore->btgStore = openGateKeeperBactigStore(name,mode); 

     if(NULLSTOREHANDLE == gkpStore->batStore ||
	NULLSTOREHANDLE == gkpStore->frgStore ||
	NULLSTOREHANDLE == gkpStore->lnkStore ||
	NULLSTOREHANDLE == gkpStore->locStore ||
	NULLSTOREHANDLE == gkpStore->seqStore ||
	NULLSTOREHANDLE == gkpStore->dstStore ||
	NULLSTOREHANDLE == gkpStore->btgStore) {
       fprintf(stderr,"**** Failure to open Gatekeeper Store ...\n");
       return 1;
     }

     sprintf(name,"%s/gkp.phash", gkpStore->storePath);
     if(mode && *mode == 'r' && *(mode + 1) == '\0'){
       gkpStore->hashTable = OpenReadOnlyPHashTable_AS(name);
     }else{
       gkpStore->hashTable = OpenPHashTable_AS(name);
     }
     if(gkpStore->hashTable == NULL){
       fprintf(stderr,"**** Failed to open GateKeeper Persistent HashTable...\n");
       return 1;
     }
     return 0;
}

int OpenGateKeeperStore(GateKeeperStore *gkpStore){
  return OpenGateKeeperStoreCommon(gkpStore,"r+");

}
int OpenReadOnlyGateKeeperStore(GateKeeperStore *gkpStore){
  return OpenGateKeeperStoreCommon(gkpStore,"r");
}


int CreateGateKeeperStore(GateKeeperStore *gkpStore){
  char name[FILENAME_MAX];

  fprintf(stderr,"*** Create store %s at cwd %s\n", gkpStore->storePath, getcwd(NULL, 256));
     sprintf(name,"%s/gkp.bat", gkpStore->storePath);
     gkpStore->batStore = createGateKeeperBatchStore(name, "bat",1); 
     sprintf(name,"%s/gkp.frg", gkpStore->storePath);
     gkpStore->frgStore = createGateKeeperFragmentStore(name, "frg",1); 
     sprintf(name,"%s/gkp.lnk", gkpStore->storePath);
     gkpStore->lnkStore = createGateKeeperLinkStore(name,"lnk",1); 
     sprintf(name,"%s/gkp.loc", gkpStore->storePath);
     gkpStore->locStore = createGateKeeperLocaleStore(name,"loc",1); 
     sprintf(name,"%s/gkp.s_loc", gkpStore->storePath);
     gkpStore->s_locStore = createGateKeeperLocaleStore(name,"loc",1); 
     sprintf(name,"%s/gkp.seq", gkpStore->storePath);
     gkpStore->seqStore = createGateKeeperSequenceStore(name,"seq",1); 
     sprintf(name,"%s/gkp.dst", gkpStore->storePath);
     gkpStore->dstStore = createGateKeeperDistanceStore(name,"dst",1); 
     sprintf(name,"%s/gkp.s_dst", gkpStore->storePath);
     gkpStore->s_dstStore = createGateKeeperDistanceStore(name,"dst",1); 
     sprintf(name,"%s/gkp.btg", gkpStore->storePath);
     gkpStore->btgStore = createGateKeeperBactigStore(name,"btg",1); 
     sprintf(name,"%s/gkp.phash", gkpStore->storePath);
     gkpStore->hashTable = CreatePHashTable_AS(2048,name);
     return 0;
}



void CloseGateKeeperStore(GateKeeperStore *gkpStore){
  fprintf(stderr,"*** Close directory %s\n", gkpStore->storePath);

  if(gkpStore->batStore != NULLSTOREHANDLE)
    closeStore(gkpStore->batStore); 
  if(gkpStore->frgStore != NULLSTOREHANDLE)
    closeStore(gkpStore->frgStore); 
  if(gkpStore->lnkStore != NULLSTOREHANDLE)
    closeStore(gkpStore->lnkStore); 
  if(gkpStore->locStore != NULLSTOREHANDLE)
    closeStore(gkpStore->locStore); 
  if(gkpStore->s_locStore != NULLSTOREHANDLE)
    closeStore(gkpStore->s_locStore); 
  if(gkpStore->seqStore != NULLSTOREHANDLE)
    closeStore(gkpStore->seqStore); 
  if(gkpStore->btgStore != NULLSTOREHANDLE)
    closeStore(gkpStore->btgStore); 
  if(gkpStore->dstStore != NULLSTOREHANDLE)
    closeStore(gkpStore->dstStore); 
  if(gkpStore->s_dstStore != NULLSTOREHANDLE)
    closeStore(gkpStore->s_dstStore); 
  if(gkpStore->hashTable != NULL)
    ClosePHashTable_AS(gkpStore->hashTable);
}
