
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

static char const *rcsid = "$Id: AS_GKP_rebuildMap.c,v 1.8 2008-12-05 19:06:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_GKP_include.h"


int
rebuildMap(char *gkpStoreName) {
   GateKeeperStore   *gkp        = openGateKeeperStore(gkpStoreName, FALSE);

   if (gkp == NULL)
     fprintf(stderr, "Failed to open %s\n", gkpStoreName), exit(1);

   {
     char name[FILENAME_MAX];
     sprintf(name,"%s/map", gkpStoreName);
     gkp->UIDtoIID = CreateScalarHashTable_AS();
     SaveHashTable_AS(name, gkp->UIDtoIID);
   }

   //  Insert batch info
   //
   {
     AS_IID                i    = 0;
     GateKeeperBatchRecord gkpb = {0};

     int64  firstElem = getFirstElemStore(gkp->bat);
     int64  lastElem  = getLastElemStore(gkp->bat);

     for (i=firstElem; i<=lastElem; i++) {
       getGateKeeperBatch(gkp, i, &gkpb);

       if (InsertInHashTable_AS(gkp->UIDtoIID, AS_UID_toInteger(gkpb.batchUID), 0, (uint64)i, AS_IID_BAT) == HASH_FAILURE)
         fprintf(stderr, "Error inserting batch %s,"F_IID" into hash table.\n",
                 AS_UID_toString(gkpb.batchUID), i);
     }
   }

   //  Insert library info
   //
   {
     AS_IID                   i    = 0;
     GateKeeperLibraryRecord *gkpl = NULL;

     int64  firstElem = getFirstElemStore(gkp->lib);
     int64  lastElem  = getLastElemStore(gkp->lib);

     for (i=firstElem; i<=lastElem; i++) {
       gkpl = getGateKeeperLibrary(gkp, i);

       if (InsertInHashTable_AS(gkp->UIDtoIID, AS_UID_toInteger(gkpl->libraryUID), 0, (uint64)i, AS_IID_LIB) == HASH_FAILURE)
         fprintf(stderr, "Error inserting library %s,"F_IID" into hash table.\n",
                 AS_UID_toString(gkpl->libraryUID), i);
     }
   }

   //  Insert fragment info
   {
     FragStream *fs = openFragStream(gkp, FRAG_S_ALL);
     fragRecord  fr = {0};

     while (nextFragStream(fs, &fr)) {
       if (InsertInHashTable_AS(gkp->UIDtoIID, AS_UID_toInteger(getFragRecordUID(&fr)), 0, (uint64)getFragRecordIID(&fr), AS_IID_FRG) == HASH_FAILURE)
         fprintf(stderr, "Error inserting UID %s and IID "F_IID"into hash table.\n",
                 AS_UID_toString(getFragRecordUID(&fr)), getFragRecordIID(&fr));
     }

     closeFragStream(fs);
   }

   //  This saves the updated hash table for us.
   //
   closeGateKeeperStore(gkp);

   return 0;
}
