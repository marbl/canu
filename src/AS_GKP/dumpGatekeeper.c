
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
static char CM_ID[] = "$Id: dumpGatekeeper.c,v 1.8 2007-01-25 09:02:12 brianwalenz Exp $";

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

int  nerrs = 0;   // Number of errors in current run
int maxerrs = 10; // Number of errors allowed before we punt

static MesgReader reader;
static MesgWriter writer;


int
main(int argc, char * argv []) {

   int status = 0;
   int  summary;
   int  quiet;
   int fragID = -1;
   char *gatekeeperStorePath;
   GateKeeperStore gkpStore;

   summary = quiet = 0;

   {
     int ch,errflg=0;
     optarg = NULL;
     while (!errflg && ((ch = getopt(argc, argv, "sqF:")) != EOF))
       switch(ch) {
       case 'F':
	 fragID = atoi(optarg);
	 break;
       case 's':
	 summary = TRUE;
	 break;
       case 'q':
	 quiet = TRUE;
	 break;
       case '?':
	 fprintf(stdout,"Unrecognized option -%c",optopt);
       default :
	 errflg++;
       }

     if(argc - optind != 1 ) {
       fprintf (stdout, "USAGE:  dumpGatekeeper [-qs] <gatekeeperStorePath>\n");
       exit (EXIT_FAILURE);
     }

     gatekeeperStorePath = argv[optind++];
   }
   


   InitGateKeeperStore(&gkpStore, gatekeeperStorePath);
   OpenReadOnlyGateKeeperStore(&gkpStore);


   /**************** DUMP Batches  *************/
   if(fragID == -1)
   {
     GateKeeperBatchRecord gkpb;
     StoreStat stat;
     int64 i;
     statsStore(gkpStore.batStore, &stat);
     fprintf(stdout,"* Stats for Batch Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     fprintf(stdout,"* Printing Batches\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperBatchStore(gkpStore.batStore,i,&gkpb);
       
       fprintf(stdout,"* Batch " F_S64 " UID:" F_UID " name:%s comment:%s created:(" F_TIME_T ") %s \n",
	       i,gkpb.UID, gkpb.name, gkpb.comment,
               gkpb.created, ctime(&gkpb.created));
       fprintf(stdout,"\t Num Fragments " F_S32 "\n", gkpb.numFragments);
       fprintf(stdout,"\t Num Locales " F_S32 "\n", gkpb.numLocales);
       fprintf(stdout,"\t Num s_Locales " F_S32 "\n", gkpb.num_s_Locales);
       fprintf(stdout,"\t Num Bactigs " F_S32 "\n", gkpb.numBactigs);
       fprintf(stdout,"\t Num Distances " F_S32 "\n", gkpb.numDistances);
       fprintf(stdout,"\t Num s_Distances " F_S32 "\n", gkpb.num_s_Distances);
       fprintf(stdout,"\t Num Links " F_S32 "\n", gkpb.numLinks);
       fprintf(stdout,"\t Num Sequences " F_S32 "\n", gkpb.numSequences);
       fprintf(stdout,"\t Num Screens " F_S32 "\n", gkpb.numScreens);
     }
     gkpb.numFragments = getNumGateKeeperFragments(gkpStore.frgStore);
     gkpb.numLocales = getNumGateKeeperLocales(gkpStore.locStore);
     gkpb.num_s_Locales = getNumGateKeeperLocales(gkpStore.s_locStore);
     gkpb.numBactigs = getNumGateKeeperBactigs(gkpStore.btgStore);
     gkpb.numDistances = getNumGateKeeperDistances(gkpStore.dstStore);
     gkpb.num_s_Distances = getNumGateKeeperDistances(gkpStore.s_dstStore);
     gkpb.numScreens = getNumGateKeeperScreens(gkpStore.scnStore);
     gkpb.numRepeats = getNumGateKeeperRepeats(gkpStore.rptStore);
     gkpb.numLinks = getNumGateKeeperLinks(gkpStore.lnkStore);
     gkpb.numSequences = getNumGateKeeperSequences(gkpStore.seqStore);
     fprintf(stdout,"* Final Stats\n");
     fprintf(stdout,"\t Num Fragments " F_S32 "\n",gkpb.numFragments );
     fprintf(stdout,"\t Num Locales " F_S32 "\n", gkpb.numLocales);
     fprintf(stdout,"\t Num s_Locales " F_S32 "\n", gkpb.num_s_Locales);
     fprintf(stdout,"\t Num Bactigs " F_S32 "\n", gkpb.numBactigs);
     fprintf(stdout,"\t Num Distances " F_S32 "\n", gkpb.numDistances);
     fprintf(stdout,"\t Num s_Distances " F_S32 "\n", gkpb.num_s_Distances);
     fprintf(stdout,"\t Num Links " F_S32 "\n", gkpb.numLinks);
     fprintf(stdout,"\t Num Sequences " F_S32 "\n", gkpb.numSequences);
     fprintf(stdout,"\t Num Screens " F_S32 "\n", gkpb.numScreens);
   }

   if(summary)
     exit(0);

   /**************** DUMP Dists  *************/
   if(fragID == -1)
   {
     GateKeeperDistanceRecord gkpd;
     StoreStat stat;
     StoreStat shadow_stat;
     int64 i;
     int64 j;
     statsStore(gkpStore.dstStore, &stat);
     fprintf(stdout,"* Stats for Dist Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stdout,"* Printing Dists\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperDistanceStore(gkpStore.dstStore,i,&gkpd);
       if (!quiet)
         fprintf(stdout,"* Dist " F_S64 " UID:" F_UID " del:%d red:%d mean:%f std:%f batch(" F_U16 "," F_U16 ") prevID:" F_IID " prevInstanceID: " F_IID "\n",
                 i,gkpd.UID, gkpd.deleted, gkpd.redefined, gkpd.mean, gkpd.stddev,
                 gkpd.birthBatch, gkpd.deathBatch, gkpd.prevID, gkpd.prevInstanceID);
     }

     if(!quiet)
       fprintf(stdout,"* Printing Shadowed Dists\n");
       
     statsStore(gkpStore.s_dstStore, &shadow_stat);
     fprintf(stdout,"* Stats for s_Dist Store are first:" F_S64 " last :" F_S64 "\n",
             shadow_stat.firstElem, shadow_stat.lastElem);
       
     j = shadow_stat.firstElem;
       
     for(j = 1; j <= shadow_stat.lastElem; j++){
       getGateKeeperDistanceStore(gkpStore.s_dstStore,j,&gkpd);
         
       if(!quiet)
         fprintf(stdout,"* Dist " F_S64 " UID:" F_UID " del:%d red:%d mean:%f std:%f batch(" F_U16 "," F_U16 ") prevID: " F_IID " prevInstanceID:" F_IID "\n",
                 j,gkpd.UID, gkpd.deleted, gkpd.redefined, gkpd.mean, gkpd.stddev,
                 gkpd.birthBatch, gkpd.deathBatch, gkpd.prevID, gkpd.prevInstanceID);
     }
   }
     
     /**************** DUMP RPTs   *************/
   if(fragID == -1)
   {
     GateKeeperRepeatRecord gkpr;
     StoreStat stat;
     int64 i;
     statsStore(gkpStore.rptStore, &stat);
     fprintf(stdout,"* Stats for Repeat Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stdout,"* Printing Repeats\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperRepeatStore(gkpStore.rptStore,i,&gkpr);
       
       if(!quiet)
         fprintf(stdout,"* Repeat " F_S64 " UID:" F_UID " which:%s \n",
                 i, gkpr.UID, gkpr.which);
     }
   }


   /**************** DUMP SCNs   *************/
   if(fragID == -1)
   {
     GateKeeperScreenRecord gkps;
     StoreStat stat;
     int64 i;
     statsStore(gkpStore.scnStore, &stat);
     fprintf(stdout,"* Stats for Screen Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stdout,"* Printing Screens\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperScreenStore(gkpStore.scnStore,i,&gkps);
       
       if(!quiet)
         fprintf(stdout,"* Screen " F_S64 " UID:" F_UID " repeatID: " F_IID " batch:(" F_U16 "," F_U16 ")\n",
                 i, gkps.UID, gkps.repeatID, gkps.birthBatch, gkps.deathBatch);
     }
   }


   /**************** DUMP BACs and BACtigs *************/
   if(fragID == -1)
   {
     StreamHandle frags = openStream(gkpStore.locStore,NULL,0);
     GateKeeperLocaleRecord gkpl;
     StoreStat stat;
     int64 i = 1;
     
     statsStore(gkpStore.locStore, &stat);
     fprintf(stdout,"* Stats for Locale Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
     
     if(!quiet)
       fprintf(stdout,"* Printing Locales\n");
     
     while(nextStream(frags, &gkpl)){

       if(!quiet){
         
	 fprintf(stdout,"* Locale " F_S64 ": typ:%c UID: " F_UID " sid:" F_IID " bac:%d len:" F_IID " del:%d red:%d has:%d prev:" F_IID " btgs: " F_S32 " first:" F_S32 " batch(" F_U16 "," F_U16 ") \n",
                 i,gkpl.type, gkpl.UID, gkpl.sequenceID, 
                 gkpl.isBac,
                 gkpl.lengthID,
                 gkpl.deleted, gkpl.redefined, gkpl.hasSequence, gkpl.prevInstanceID, gkpl.numBactigs, gkpl.firstBactig,
                 gkpl.birthBatch, gkpl.deathBatch);
	 if(gkpl.redefined){
	   fprintf(stdout,"\tREDEFINED: prevID " F_IID "\n", gkpl.prevID);
	 }
       }
       if(gkpl.type == AS_UNFINISHED ||
          gkpl.type == AS_FINISHED){
         GateKeeperSequenceRecord gkpseq;
         getGateKeeperSequenceStore(gkpStore.seqStore, gkpl.sequenceID, &gkpseq);
         assert(gkpseq.localeID == i);
         if(!quiet)
           fprintf(stdout,"\tSequence " F_UID "\n", gkpseq.UID);
       }
       if(gkpl.numBactigs > 0){
         int32 j;
         int cnt;
         
         fprintf(stdout,"\tBactigs\n");
         for(cnt = 0, j = gkpl.firstBactig; cnt < gkpl.numBactigs; cnt++, j++){
           GateKeeperBactigRecord gkpbtg;
           
           getGateKeeperBactigStore(gkpStore.btgStore, j, &gkpbtg);
           assert(gkpbtg.bacID == i);
           assert(gkpbtg.seqID == gkpl.sequenceID);
           if(!quiet)
             fprintf(stdout,"\t\t id:" F_S32 " UID:" F_UID " length:" F_IID " del:%d has:%d\n",
                     j, gkpbtg.UID, gkpbtg.length, gkpbtg.deleted, gkpbtg.hasSequence);
         }
         
       }
       i++;
     }
   }

   
   if(fragID == -1)
   {
     StreamHandle frags = openStream(gkpStore.s_locStore,NULL,0);
     GateKeeperLocaleRecord gkpl;
     StoreStat stat;
     int64 i = 1;
     
     statsStore(gkpStore.s_locStore, &stat);
     fprintf(stdout,"* Stats for Shadow Locale Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
     
     fprintf(stdout,"* Printing Shadowed Locales\n");
     
     while(nextStream(frags, &gkpl)){
       
       if(!quiet)
         fprintf(stdout,"* Locale " F_S64 ": typ:%c UID: " F_UID " sid:" F_IID " bac:%d len:" F_IID " del:%d red:%d has:%d prev:" F_IID " btgs:" F_S32 " first:" F_S32 " batch(" F_U16 "," F_U16 ")\n",
                 i,gkpl.type, gkpl.UID, gkpl.sequenceID, 
                 gkpl.isBac,
                 gkpl.lengthID,
                 gkpl.deleted, gkpl.redefined, gkpl.hasSequence, gkpl.prevInstanceID, gkpl.numBactigs, gkpl.firstBactig,
                 gkpl.birthBatch, gkpl.deathBatch);
       if(gkpl.redefined){
	 if(!quiet)
	   fprintf(stdout,"\tREDEFINED: prevID " F_IID "\n", gkpl.prevID);
       }
       
       if(gkpl.type == AS_UNFINISHED ||
          gkpl.type == AS_FINISHED){
         GateKeeperSequenceRecord gkpseq;
         getGateKeeperSequenceStore(gkpStore.seqStore, gkpl.sequenceID, &gkpseq);
	 if(!quiet)
	   fprintf(stdout,"\tSequence " F_UID "\n", gkpseq.UID);
       }
       if(gkpl.numBactigs > 0){
         int64 j;
         int cnt;
         
	 if(!quiet)
	   fprintf(stdout,"\tBactigs\n");
         for(cnt = 0, j = gkpl.firstBactig; cnt < gkpl.numBactigs; cnt++, j++){
           GateKeeperBactigRecord gkpbtg;
           
           getGateKeeperBactigStore(gkpStore.btgStore, j, &gkpbtg);
           if(!quiet)
             fprintf(stdout,"\t\t id:" F_S64 " UID:" F_UID " length:%d del:%d has:%d\n",
                     j, gkpbtg.UID, gkpbtg.length, gkpbtg.deleted, gkpbtg.hasSequence);
         }
         
       }
       i++;
     }
   }

   /**************** DUMP Fragments and Links *************/
   {
     StreamHandle frags = openStream(gkpStore.frgStore,NULL,0);
     GateKeeperFragmentRecord gkpf;
     StoreStat stat;
     int64 i = 1;
     PHashValue_AS value;
     
     statsStore(gkpStore.frgStore, &stat);
     fprintf(stdout,"* Stats are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     if(fragID != -1){
       resetStream(frags,fragID, fragID + 1);
       i = fragID;
       if(!quiet)
         fprintf(stdout,"* Printing fragments %d-%d\n", fragID, fragID + 1);
     }else{
       resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
       if(!quiet)
         fprintf(stdout,"* Printing ALL fragments \n");
     }
     
     while(nextStream(frags, &gkpf)){
       GateKeeperLinkRecordIterator iterator;
       GateKeeperLinkRecord link;

       if(HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore.hashTable, 
                                                    UID_NAMESPACE_AS,
                                                    gkpf.readUID, 
                                                    (gkpf.type == AS_BACTIG?AS_IID_BTG:AS_IID_FRAG), 
                                                    FALSE,
                                                    stdout,
                                                    &value)){
         
         if(!quiet){
           if(!gkpf.deleted)
             fprintf(stdout,"# *****ERROR******");
           else
             fprintf(stdout,"# Deleted Fragment ");
           
           fprintf(stdout,F_S64 "(" F_IID "): UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                   i, value.IID, gkpf.readUID, 
                   gkpf.type,
                   value.refCount, gkpf.numLinks, gkpf.linkHead,
                   gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
         }
       }else{
         if(!quiet){
           if(!gkpf.deleted){
             fprintf(stdout,"* Fragment " F_S64 ": UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                     i, gkpf.readUID, 
                     gkpf.type,
                     value.refCount, gkpf.numLinks, gkpf.linkHead,
                     gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
           }else{
             fprintf(stdout,"* Redefined Fragment " F_S64 " (now " F_IID "): UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                     i, value.IID, gkpf.readUID, 
                     gkpf.type,
                     value.refCount, gkpf.numLinks, gkpf.linkHead,
                     gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
           }
         }
         if(gkpf.numLinks > 0){
           CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpf.linkHead,
                                              i, &iterator);
           while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
             if(!quiet)
               fprintf(stdout,"\tLink (" F_IID "," F_IID ") dist:" F_IID " type:%c ori:%c\n",
                       link.frag1, link.frag2, link.distance, link.type,
                       getLinkOrientation( &link ));
             
           }
         }
       }
       i++;
     }
   }

   exit(status != GATEKEEPER_SUCCESS);
}
