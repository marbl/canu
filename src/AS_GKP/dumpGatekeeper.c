
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
static char CM_ID[] = "$Id: dumpGatekeeper.c,v 1.14 2007-02-12 22:16:57 brianwalenz Exp $";

/* Dump the gatekeeper stores for debug */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"

int
main(int argc, char **argv) {
   GateKeeperStore   *gkpStore;
   char              *gkpStoreName = NULL;

   int                dumpBatchInfo    = 0;
   int                dumpLibraryInfo  = 0;
   int                dumpFragmentInfo = 0;

   int                firstFrag = 0;
   int                lastFrag  = 0;

   int arg = 1;
   int err = 0;
   while (arg < argc) {
     if        (strcmp(argv[arg], "-b") == 0) {
       firstFrag        = atoi(argv[++arg]);
       dumpFragmentInfo = 1;
     } else if (strcmp(argv[arg], "-e") == 0) {
       lastFrag         = atoi(argv[++arg]);
       dumpFragmentInfo = 1;
     } else if (strcmp(argv[arg], "-F") == 0) {
       dumpFragmentInfo = 1;
     } else if (strcmp(argv[arg], "-L") == 0) {
       dumpLibraryInfo  = 1;
     } else if (strcmp(argv[arg], "-B") == 0) {
       dumpBatchInfo    = 1;
     } else if (strcmp(argv[arg], "-g") == 0) {
       gkpStoreName = argv[++arg];
     } else {
       fprintf(stderr, "unknown option '%s'\n", argv[arg]);
     }

     arg++;
   }

   if ((err) || (gkpStoreName == NULL)) {
     fprintf(stderr, "usage: %s [-b begin] [-e end] [-B] [-L] [-F] -g gkpStore]\n", argv[0]);
     fprintf(stderr, "  -b begin   dump frags from iid begin to end\n");
     fprintf(stderr, "  -e end     dump frags from iid begin to end\n");
     fprintf(stderr, "  -B         dump batch records\n");
     fprintf(stderr, "  -L         dump library records\n");
     fprintf(stderr, "  -F         dump fragment records\n");
     exit(1);
   }

   gkpStore = openGateKeeperStore(gkpStoreName, FALSE);


   if (dumpBatchInfo) {
     StoreStat stat;
     int       i;

     statsStore(gkpStore->bat, &stat);
     
     for (i=stat.firstElem; i<=stat.lastElem; i++) {
       GateKeeperBatchRecord gkpb;
       time_t                createdtime;

       getGateKeeperBatchStore(gkpStore->bat, i, &gkpb);

       createdtime = (time_t)gkpb.created;

       fprintf(stdout, "<batch>\n");
       fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpb.UID, i);
       fprintf(stdout, "<name>%s</name>\n", gkpb.name);
       fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpb.created, ctime(&createdtime));
       if (gkpb.comment[0])
         fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpb.comment);
       fprintf(stdout, "<nFrags>"F_S32"</nFrags>\n", gkpb.numFragments);
       fprintf(stdout, "<nLibs>"F_S32"</nLibs)\n", gkpb.numLibraries);
       fprintf(stdout, "<nLibsS>"F_S32"</nLibsS>\n", gkpb.numLibraries_s);
       fprintf(stdout, "</batch>\n");
     }

     fprintf(stdout, "<total>\n");
     fprintf(stdout, "<nFragsTotal>"F_S32"</nFragsTotal>\n", getNumGateKeeperFragments(gkpStore->frg));
     fprintf(stdout, "<nLibsTotal>"F_S32"</nLibsTotal>\n", getNumGateKeeperLibrarys(gkpStore->lib));
     fprintf(stdout, "<nLibsSTotal>"F_S32"</nLibsSTotal>\n", getNumGateKeeperLibrarys(gkpStore->lis));
     fprintf(stdout, "</total>\n");
   }





   if (dumpLibraryInfo) {
     StoreStat stat;
     int       i;

     statsStore(gkpStore->lib, &stat);

     for (i=stat.firstElem; i<=stat.lastElem; i++) {
       GateKeeperLibraryRecord gkpl;
       time_t                  createdtime;

       getGateKeeperLibraryStore(gkpStore->lib, i, &gkpl);

       createdtime = (time_t)gkpl.created;

       fprintf(stdout, "<library>\n");
       fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpl.UID, i);
       fprintf(stdout, "<name>%s</name>\n", gkpl.name);
       fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpl.created, ctime(&createdtime));
       if (gkpl.comment[0])
         fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpl.comment);
       fprintf(stdout, "<deleted>%d</deleted>\n", gkpl.deleted);
       fprintf(stdout, "<redefined>%d</redefined>\n", gkpl.redefined);
       fprintf(stdout, "<orientation>%d</orientation>\n", gkpl.orientation);
       fprintf(stdout, "<mean>%f</mean>\n", gkpl.mean);
       fprintf(stdout, "<stddev>%f</stddev>\n", gkpl.stddev);
       fprintf(stdout, "<numFeatures>%d</numFeatures>\n", gkpl.numFeatures);
       fprintf(stdout, "<prevInstanceID>"F_IID"</prevInstanceID>\n", gkpl.prevInstanceID);
       fprintf(stdout, "<prevID>"F_IID"</prevID>\n", gkpl.prevID);
       fprintf(stdout, "<birthBatch>%d</birthBatch>\n", gkpl.birthBatch);
       fprintf(stdout, "<deathBatch>%d</deathBatch>\n", gkpl.deathBatch);
       fprintf(stdout, "</library>\n");
     }

       
     statsStore(gkpStore->lis, &stat);

     for (i=stat.firstElem; i<=stat.lastElem; i++) {
       GateKeeperLibraryRecord gkpl;
       time_t                  createdtime;

       getGateKeeperLibraryStore(gkpStore->lis, i, &gkpl);

       createdtime = (time_t)gkpl.created;

       fprintf(stdout, "<shadow_library>\n");
       fprintf(stdout, "<ident>"F_UID","F_IID"</ident>\n", gkpl.UID, i);
       fprintf(stdout, "<name>%s</name>\n", gkpl.name);
       fprintf(stdout, "<created>("F_U64") %s</created>\n", gkpl.created, ctime(&createdtime));
       if (gkpl.comment[0])
         fprintf(stdout, "<comment>\n%s\n</comment>\n", gkpl.comment);
       fprintf(stdout, "<deleted>%d</deleted>\n", gkpl.deleted);
       fprintf(stdout, "<redefined>%d</redefined>\n", gkpl.redefined);
       fprintf(stdout, "<orientation>%d</orientation>\n", gkpl.orientation);
       fprintf(stdout, "<mean>%f</mean>\n", gkpl.mean);
       fprintf(stdout, "<stddev>%f</stddev>\n", gkpl.stddev);
       fprintf(stdout, "<numFeatures>%d</numFeatures>\n", gkpl.numFeatures);
       fprintf(stdout, "<prevInstanceID>"F_IID"</prevInstanceID>\n", gkpl.prevInstanceID);
       fprintf(stdout, "<prevID>"F_IID"</prevID>\n", gkpl.prevID);
       fprintf(stdout, "<birthBatch>%d</birthBatch>\n", gkpl.birthBatch);
       fprintf(stdout, "<deathBatch>%d</deathBatch>\n", gkpl.deathBatch);
       fprintf(stdout, "</shadow_library>\n");
     }
   }
     


   if (dumpFragmentInfo) {
     StreamHandle             frags = openStream(gkpStore->frg, NULL, 0);
     GateKeeperFragmentRecord gkpf;
     StoreStat                stat;
     int                      i;
     PHashValue_AS            value;

     statsStore(gkpStore->frg, &stat);

     if (firstFrag < stat.firstElem)
       firstFrag = stat.firstElem;
     if ((lastFrag > stat.lastElem) || (lastFrag == 0))
       lastFrag = stat.lastElem;

     resetStream(frags, firstFrag, lastFrag);
     i = firstFrag;

     while (nextStream(frags, &gkpf)) {
       if (gkpf.deleted) {
         fprintf(stdout," * Deleted Fragment ("F_UID","F_IID") refs:%d batch("F_U16","F_U16")\n",
                 gkpf.UID,
                 i,
                 value.refCount,
                 gkpf.birthBatch,
                 gkpf.deathBatch);
       } else {
         fprintf(stdout,"* Fragment ("F_UID","F_IID") refs:%d  batch("F_U16","F_U16")\n",
                 gkpf.UID, 
                 i,
                 value.refCount,
                 gkpf.birthBatch,
                 gkpf.deathBatch);
       }

#if 0
       fprintf(stdout,"\tLink (" F_IID "," F_IID ") dist:" F_IID " type:%c ori:%c\n",
               link.frag1, link.frag2, link.distance, link.type,
               getLinkOrientation( &link ));
#endif

       i++;
     }
   }

   exit(0);
}
