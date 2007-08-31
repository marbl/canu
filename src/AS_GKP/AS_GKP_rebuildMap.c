
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

static char const rcsid[] = "$Id: AS_GKP_rebuildMap.c,v 1.2 2007-08-31 21:06:16 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_GKP_include.h"


int
rebuildMap(char *hashFileName, 
           char *gkpStoreName) {
   HashTable_AS      *UIDtoIID   = CreateScalarHashTable_AS(32 * 1024);
   GateKeeperStore   *gkp        = openGateKeeperStore(gkpStoreName, FALSE);
   fragRecord        *fr;
   FragStream        *fs;
   StoreStat          stat;

   CDS_IID_t          i;
   int                success;
    
   if (gkp == NULL) {
      fprintf(stderr, "Failed to open %s\n", gkpStoreName);
      exit(1);
   }

   fr = new_fragRecord();
   fs = openFragStream(gkp, FRAG_S_ALL);

   statsStore(gkp->frg, &stat);
   resetFragStream(fs, stat.firstElem, stat.lastElem);

   // read all fragments from the store and enter them into our map
   while (nextFragStream(fs, fr)) {
      success = InsertInHashTable_AS(UIDtoIID, getFragRecordUID(fr), 0, (uint64)getFragRecordIID(fr), AS_IID_LIB);
      
      if (success == HASH_FAILURE) {  
         fprintf(stderr, "Error inserting UID"F_UID" and IID "F_IID"into hash table.\n", getFragRecordUID(fr), getFragRecordIID(fr));
         closeFragStream(fs);
         closeGateKeeperStore(gkp);
         
         exit(1);
      }
   }
   closeFragStream(fs);  
  
   // now add library info to map
   statsStore(gkp->lib, &stat);
   for (i = stat.firstElem; i <= stat.lastElem; i++) {    
      GateKeeperLibraryRecord *gkpl = getGateKeeperLibrary(gkp, i);
      LibraryMesg              lmesg;

      AS_PER_encodeLibraryFeatures(gkpl, &lmesg);
      success = InsertInHashTable_AS(UIDtoIID, gkpl->libraryUID, 0, (uint64)i, AS_IID_LIB);
      AS_PER_encodeLibraryFeaturesCleanup(&lmesg);
      
      if (success == HASH_FAILURE) {  
         fprintf(stderr, "Error inserting UID"F_UID" and IID "F_IID"into hash table.\n", gkpl->libraryUID, i);
         closeFragStream(fs);
         closeGateKeeperStore(gkp);
         
         exit(1);
      }
   }
  
   // now add the batch info to the map
   statsStore(gkp->bat, &stat);   
   for (i = stat.firstElem; i <= stat.lastElem; i++) {
      GateKeeperBatchRecord gkpb;
      getGateKeeperBatch(gkp, i, &gkpb);

      success = InsertInHashTable_AS(UIDtoIID, gkpb.batchUID, 0, (uint64)i, AS_IID_LIB);
      
      if (success == HASH_FAILURE) {  
         fprintf(stderr, "Error inserting UID"F_UID" and IID "F_IID"into hash table.\n", gkpb.batchUID, i);
         closeFragStream(fs);
         closeGateKeeperStore(gkp);
         
         exit(1);
      }
   }
      
   closeGateKeeperStore(gkp);
  
   // now dump the hash map we build
   SaveHashTable_AS(hashFileName, UIDtoIID);
   DeleteHashTable_AS(UIDtoIID);
   
   return 0;
}
