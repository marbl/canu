
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
static char CM_ID[] = "$Id: dumpGatekeeper.c,v 1.1.1.1 2004-04-14 13:51:39 catmandew Exp $";

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


int  main(int argc, char * argv [])

{
   int status = 0;
#if 0
   char  Output_File_Name [FILENAME_MAX];
   char cmd[FILENAME_MAX * 4];
   char tmpFilePath[FILENAME_MAX];
#endif
   int  summary;
   int  quiet;
   int fragID = -1;
   char *gatekeeperStorePath;
   GateKeeperStore gkpStore;

   summary = quiet = 0;
   /**************** Process Command Line Arguments *********************/
   { /* Parse the argument list using "man 3 getopt". */ 
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
	 fprintf(stderr,"Unrecognized option -%c",optopt);
       default :
	 errflg++;
       }

     
     if(argc - optind != 1 )
       {
	 fprintf (stderr, "USAGE:  dumpGatekeeper [-qs] <gatekeeperStorePath>\n");
	 exit (EXIT_FAILURE);
       }

     gatekeeperStorePath = argv[optind++];

     /* End of command line parsing */
   }
   


   /**************** Open or Create Files *********************/
   fprintf(stderr,"* GatekeeperStorePath is %s\n",
	   gatekeeperStorePath);

   InitGateKeeperStore(&gkpStore, gatekeeperStorePath);
   OpenReadOnlyGateKeeperStore(&gkpStore);


   /**************** DUMP Batches  *************/
   if(fragID == -1)
   {
     GateKeeperBatchRecord gkpb;
     StoreStat stat;
     int64 i;
     statsStore(gkpStore.batStore, &stat);
     fprintf(stderr,"* Stats for Batch Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     fprintf(stderr,"* Printing Batches\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperBatchStore(gkpStore.batStore,i,&gkpb);
       
       fprintf(stderr,"* Batch " F_S64 " UID:" F_UID " name:%s comment:%s created:(" F_TIME_T ") %s \n",
	       i,gkpb.UID, gkpb.name, gkpb.comment,
               gkpb.created, ctime(&gkpb.created));
       fprintf(stderr,"\t Num Fragments " F_S32 "\n", gkpb.numFragments);
       fprintf(stderr,"\t Num Locales " F_S32 "\n", gkpb.numLocales);
       fprintf(stderr,"\t Num s_Locales " F_S32 "\n", gkpb.num_s_Locales);
       fprintf(stderr,"\t Num Bactigs " F_S32 "\n", gkpb.numBactigs);
       fprintf(stderr,"\t Num Distances " F_S32 "\n", gkpb.numDistances);
       fprintf(stderr,"\t Num s_Distances " F_S32 "\n", gkpb.num_s_Distances);
       fprintf(stderr,"\t Num Links " F_S32 "\n", gkpb.numLinks);
       fprintf(stderr,"\t Num Sequences " F_S32 "\n", gkpb.numSequences);
       fprintf(stderr,"\t Num Screens " F_S32 "\n", gkpb.numScreens);
       fprintf(stderr,"\t Num Plates " F_S32 "\n", gkpb.numPlates);
       fprintf(stderr,"\t Num Wells " F_S32 "\n", gkpb.numWells);
     }
     gkpb.numFragments = getNumGateKeeperFragments(gkpStore.frgStore);
     gkpb.numLocales = getNumGateKeeperLocales(gkpStore.locStore);
     gkpb.num_s_Locales = getNumGateKeeperLocales(gkpStore.s_locStore);
     gkpb.numBactigs = getNumGateKeeperBactigs(gkpStore.btgStore);
     gkpb.numDistances = getNumGateKeeperDistances(gkpStore.dstStore);
     gkpb.num_s_Distances = getNumGateKeeperDistances(gkpStore.s_dstStore);
     gkpb.numScreens = getNumGateKeeperScreens(gkpStore.scnStore);
     gkpb.numRepeats = getNumGateKeeperRepeats(gkpStore.rptStore);
     gkpb.numPlates = getNumGateKeeperSequencePlates(gkpStore.sqpStore);
     gkpb.numWells = getNumGateKeeperWells(gkpStore.welStore);
     gkpb.numLinks = getNumGateKeeperLinks(gkpStore.lnkStore);
     gkpb.numSequences = getNumGateKeeperSequences(gkpStore.seqStore);
     fprintf(stderr,"* Final Stats\n");
     fprintf(stderr,"\t Num Fragments " F_S32 "\n",gkpb.numFragments );
     fprintf(stderr,"\t Num Locales " F_S32 "\n", gkpb.numLocales);
     fprintf(stderr,"\t Num s_Locales " F_S32 "\n", gkpb.num_s_Locales);
     fprintf(stderr,"\t Num Bactigs " F_S32 "\n", gkpb.numBactigs);
     fprintf(stderr,"\t Num Distances " F_S32 "\n", gkpb.numDistances);
     fprintf(stderr,"\t Num s_Distances " F_S32 "\n", gkpb.num_s_Distances);
     fprintf(stderr,"\t Num Links " F_S32 "\n", gkpb.numLinks);
     fprintf(stderr,"\t Num Sequences " F_S32 "\n", gkpb.numSequences);
     fprintf(stderr,"\t Num Screens " F_S32 "\n", gkpb.numScreens);
     fprintf(stderr,"\t Num Plates " F_S32 "\n", gkpb.numPlates);
     fprintf(stderr,"\t Num Wells " F_S32 "\n", gkpb.numWells);
     
     
   }

   if(summary)
     exit(0);

   /**************** DUMP Dists  *************/
   if(fragID == -1)
   {
     GateKeeperDistanceRecord gkpd;
     GateKeeperLibDonorRecord gkpldr;
     StoreStat stat;
     StoreStat shadow_stat;
     int64 i;
     int64 j;
     statsStore(gkpStore.dstStore, &stat);
     fprintf(stderr,"* Stats for Dist Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stderr,"* Printing Dists\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperDistanceStore(gkpStore.dstStore,i,&gkpd);
       getGateKeeperLibDonorStore(gkpStore.libStore,i,&gkpldr);
       
       if(!quiet){
         fprintf(stderr,"* Dist " F_S64 " UID:" F_UID " del:%d red:%d mean:%f std:%f batch(" F_U16 "," F_U16 ") prevID:" F_IID " prevInstanceID: " F_IID "\n",
                 i,gkpd.UID, gkpd.deleted, gkpd.redefined, gkpd.mean, gkpd.stddev,
                 gkpd.birthBatch, gkpd.deathBatch, gkpd.prevID, gkpd.prevInstanceID);
         if(gkpldr.set){
           fprintf(stderr,"\t* Donor: %u\n",gkpldr.idonor);
         }else{
           fprintf(stderr,"\t* No Donor set\n");
         }
       }
       
       if(!quiet)
         fprintf(stderr,"* Printing Shadowed Dists\n");
       
       statsStore(gkpStore.s_dstStore, &shadow_stat);
       fprintf(stderr,"* Stats for s_Dist Store are first:" F_S64 " last :" F_S64 "\n",
               shadow_stat.firstElem, shadow_stat.lastElem);
       
       j = shadow_stat.firstElem;
       
       for(j = 1; j <= shadow_stat.lastElem; j++){
         getGateKeeperDistanceStore(gkpStore.s_dstStore,j,&gkpd);
         
         if(!quiet)
           fprintf(stderr,"* Dist " F_S64 " UID:" F_UID " del:%d red:%d mean:%f std:%f batch(" F_U16 "," F_U16 ") prevID: " F_IID " prevInstanceID:" F_IID "\n",
                   j,gkpd.UID, gkpd.deleted, gkpd.redefined, gkpd.mean, gkpd.stddev,
                   gkpd.birthBatch, gkpd.deathBatch, gkpd.prevID, gkpd.prevInstanceID);
       }
     }
   }
     
     /**************** DUMP RPTs   *************/
   if(fragID == -1)
   {
     GateKeeperRepeatRecord gkpr;
     StoreStat stat;
     int64 i;
     statsStore(gkpStore.rptStore, &stat);
     fprintf(stderr,"* Stats for Repeat Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stderr,"* Printing Repeats\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperRepeatStore(gkpStore.rptStore,i,&gkpr);
       
       if(!quiet)
         fprintf(stderr,"* Repeat " F_S64 " UID:" F_UID " which:%s \n",
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
     fprintf(stderr,"* Stats for Screen Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     i = stat.firstElem;
     
     if(!quiet)
       fprintf(stderr,"* Printing Screens\n");
     
     for(i = 1; i <= stat.lastElem; i++){
       getGateKeeperScreenStore(gkpStore.scnStore,i,&gkps);
       
       if(!quiet)
         fprintf(stderr,"* Screen " F_S64 " UID:" F_UID " repeatID: " F_IID " batch:(" F_U16 "," F_U16 ")\n",
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
     fprintf(stderr,"* Stats for Locale Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
     
     if(!quiet)
       fprintf(stderr,"* Printing Locales\n");
     
     while(nextStream(frags, &gkpl)){

       if(!quiet){
         
	 fprintf(stderr,"* Locale " F_S64 ": typ:%c UID: " F_UID " sid:" F_IID " bac:%d len:" F_IID " del:%d red:%d has:%d prev:" F_IID " btgs: " F_S32 " first:" F_S32 " batch(" F_U16 "," F_U16 ") \n",
                 i,gkpl.type, gkpl.UID, gkpl.sequenceID, 
                 gkpl.isBac,
                 gkpl.lengthID,
                 gkpl.deleted, gkpl.redefined, gkpl.hasSequence, gkpl.prevInstanceID, gkpl.numBactigs, gkpl.firstBactig,
                 gkpl.birthBatch, gkpl.deathBatch);
	 if(gkpl.redefined){
	   fprintf(stderr,"\tREDEFINED: prevID " F_IID "\n", gkpl.prevID);
	 }
       }
       if(gkpl.type == AS_UNFINISHED ||
          gkpl.type == AS_FINISHED){
         GateKeeperSequenceRecord gkpseq;
         getGateKeeperSequenceStore(gkpStore.seqStore, gkpl.sequenceID, &gkpseq);
         assert(gkpseq.localeID == i);
         if(!quiet)
           fprintf(stderr,"\tSequence " F_UID "\n", gkpseq.UID);
       }
       if(gkpl.numBactigs > 0){
         int32 j;
         int cnt;
         
         fprintf(stderr,"\tBactigs\n");
         for(cnt = 0, j = gkpl.firstBactig; cnt < gkpl.numBactigs; cnt++, j++){
           GateKeeperBactigRecord gkpbtg;
           
           getGateKeeperBactigStore(gkpStore.btgStore, j, &gkpbtg);
           assert(gkpbtg.bacID == i);
           assert(gkpbtg.seqID == gkpl.sequenceID);
           if(!quiet)
             fprintf(stderr,"\t\t id:" F_S32 " UID:" F_UID " length:" F_IID " del:%d has:%d\n",
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
     fprintf(stderr,"* Stats for Shadow Locale Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
     
     fprintf(stderr,"* Printing Shadowed Locales\n");
     
     while(nextStream(frags, &gkpl)){
       
       if(!quiet)
         fprintf(stderr,"* Locale " F_S64 ": typ:%c UID: " F_UID " sid:" F_IID " bac:%d len:" F_IID " del:%d red:%d has:%d prev:" F_IID " btgs:" F_S32 " first:" F_S32 " batch(" F_U16 "," F_U16 ")\n",
                 i,gkpl.type, gkpl.UID, gkpl.sequenceID, 
                 gkpl.isBac,
                 gkpl.lengthID,
                 gkpl.deleted, gkpl.redefined, gkpl.hasSequence, gkpl.prevInstanceID, gkpl.numBactigs, gkpl.firstBactig,
                 gkpl.birthBatch, gkpl.deathBatch);
       if(gkpl.redefined){
	 if(!quiet)
	   fprintf(stderr,"\tREDEFINED: prevID " F_IID "\n", gkpl.prevID);
       }
       
       if(gkpl.type == AS_UNFINISHED ||
          gkpl.type == AS_FINISHED){
         GateKeeperSequenceRecord gkpseq;
         getGateKeeperSequenceStore(gkpStore.seqStore, gkpl.sequenceID, &gkpseq);
	 if(!quiet)
	   fprintf(stderr,"\tSequence " F_UID "\n", gkpseq.UID);
       }
       if(gkpl.numBactigs > 0){
         int64 j;
         int cnt;
         
	 if(!quiet)
	   fprintf(stderr,"\tBactigs\n");
         for(cnt = 0, j = gkpl.firstBactig; cnt < gkpl.numBactigs; cnt++, j++){
           GateKeeperBactigRecord gkpbtg;
           
           getGateKeeperBactigStore(gkpStore.btgStore, j, &gkpbtg);
           if(!quiet)
             fprintf(stderr,"\t\t id:" F_S64 " UID:" F_UID " length:%d del:%d has:%d\n",
                     j, gkpbtg.UID, gkpbtg.length, gkpbtg.deleted, gkpbtg.hasSequence);
         }
         
       }
       i++;
     }
   }

   /**************** DUMP Fragments and Links *************/
   /**************** Also dump aux Frags ******************/
   {
     StreamHandle frags = openStream(gkpStore.frgStore,NULL,0);
     GateKeeperFragmentRecord gkpf;
     GateKeeperAuxFragRecord gkpafr;
     StoreStat stat;
     int64 i = 1;
     PHashValue_AS value;
     
     statsStore(gkpStore.frgStore, &stat);
     fprintf(stderr,"* Stats are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     if(fragID != -1){
       resetStream(frags,fragID, fragID + 1);
       i = fragID;
       if(!quiet)
         fprintf(stderr,"* Printing fragments %d-%d\n", fragID, fragID + 1);
     }else{
       resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
       if(!quiet)
         fprintf(stderr,"* Printing ALL fragments \n");
     }
     
     while(nextStream(frags, &gkpf)){
       GateKeeperLinkRecordIterator iterator;
       GateKeeperLinkRecord link;
       
       getGateKeeperAuxFragStore(gkpStore.auxStore, i, &gkpafr);
       
       if(HASH_SUCCESS != LookupTypeInPHashTable_AS(gkpStore.hashTable, 
                                                    UID_NAMESPACE_AS,
                                                    gkpf.readUID, 
                                                    (gkpf.type == AS_BACTIG?AS_IID_BTG:AS_IID_FRAG), 
                                                    FALSE,
                                                    stderr,
                                                    &value)){
         
         if(!quiet){
           if(!gkpf.deleted)
             fprintf(stderr,"# *****ERROR******");
           else
             fprintf(stderr,"# Deleted Fragment ");
           
           fprintf(stderr,F_S64 "(" F_IID "): UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                   i, value.IID, gkpf.readUID, 
                   gkpf.type,
                   value.refCount, gkpf.numLinks, gkpf.linkHead,
                   gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
           if(gkpafr.set){
             if(!gkpafr.deleted)
               fprintf(stderr,"# *****ERROR******");
             else
               fprintf(stderr,"# Deleted Auxiliary Fragment ");
             fprintf(stderr,"\tplate:" F_IID ", well:%u, lib:" F_IID "\n",
                     gkpafr.iplate, gkpafr.iwell, gkpafr.ilib);
           }else{
             fprintf(stderr,"# No Auxiliary Fragment\n");
           }
         }
       }else{
         if(!quiet){
           if(!gkpf.deleted){
             fprintf(stderr,"* Fragment " F_S64 ": UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                     i, gkpf.readUID, 
                     gkpf.type,
                     value.refCount, gkpf.numLinks, gkpf.linkHead,
                     gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
           }else{
             fprintf(stderr,"* Redefined Fragment " F_S64 " (now " F_IID "): UID:" F_UID " type:%c refs:%d links:%u(" F_IID ") lID:" F_IID " sID:" F_IID " bID:" F_IID " batch(" F_U16 "," F_U16 ")\n",
                     i, value.IID, gkpf.readUID, 
                     gkpf.type,
                     value.refCount, gkpf.numLinks, gkpf.linkHead,
                     gkpf.localeID, gkpf.seqID, gkpf.bactigID, gkpf.birthBatch, gkpf.deathBatch);
           }
           if(gkpafr.set){
             if(gkpafr.deleted)
               fprintf(stderr,"# *****ERROR****** - Deleted Auxiliary Fragment");
             else
               fprintf(stderr,"# Auxiliary Fragment ");
             fprintf(stderr,"\tplate:" F_IID ", well:%u, lib:" F_IID "\n",
                     gkpafr.iplate, gkpafr.iwell, gkpafr.ilib);
           }else{
             fprintf(stderr,"# No Auxiliary Fragment\n");
           }
         }
         if(gkpf.numLinks > 0){
           CreateGateKeeperLinkRecordIterator(gkpStore.lnkStore, gkpf.linkHead,
                                              i, &iterator);
           while(NextGateKeeperLinkRecordIterator(&iterator, &link)){
             if(!quiet)
               fprintf(stderr,"\tLink (" F_IID "," F_IID ") dist:" F_IID " type:%c ori:%c\n",
                       link.frag1, link.frag2, link.distance, link.type,
                       getLinkOrientation( &link ));
             
           }
         }
       }
       i++;
     }
   }


   /**************** DUMP Plates *************/
   if(fragID == -1){
     StreamHandle frags = openStream(gkpStore.sqpStore,NULL,0);
     GateKeeperSequencePlateRecord gkpsqp;
     StoreStat stat;
     int i = 1;
     
     statsStore(gkpStore.sqpStore, &stat);
     fprintf(stderr,"* Stats for Sequence Plate Store are first:" F_S64 " last :" F_S64 "\n",
	     stat.firstElem, stat.lastElem);
     
     resetStream(frags,STREAM_FROMSTART, STREAM_UNTILEND);
     
     if(!quiet)
       fprintf(stderr,"* Printing Sequence Plates\n");
     
     while(nextStream(frags, &gkpsqp)){

       if(!quiet){
         
	 fprintf(stderr,"* Plate %d: UID: " F_UID " deleted:%d firstWell:" F_IID " numWells:%u mate:" F_IID "\n",
                 i, gkpsqp.UID, gkpsqp.deleted,
                 gkpsqp.firstWell, gkpsqp.numWells, gkpsqp.mate);
       }
       
       if(gkpsqp.numWells > 0){
         int j;
         int cnt;
         
         fprintf(stderr,"\tWells\n");
         for(cnt = 0, j = gkpsqp.firstWell; cnt < gkpsqp.numWells; cnt++, j++){
           GateKeeperWellRecord gkpwel;
           
           getGateKeeperWellStore(gkpStore.welStore, j, &gkpwel);
           if(!quiet)
             fprintf(stderr,"\t\t well: %u, deleted: %u, frag: " F_IID ", lib: " F_IID "\n",
                     gkpwel.ewell, gkpwel.deleted, gkpwel.ifrag, gkpwel.ilib);
         }
       }
       i++;
     }
   }
   
#if 0
   /**************** Close files   *********************/
   if(status == GATEKEEPER_SUCCESS){ //  OK to update persistent data and generate output
     /* Remove temporary files if any */
     fprintf(stderr,"#  Successful run with %d errors < %d maxerrs ..removing temp backup files\n",
	     nerrs , maxerrs);
     sprintf(cmd,"rm -rf %s", tmpFilePath);
     if(system(cmd) != 0) assert(0);
   }else{
     fprintf(stderr,"# Too Many Errors -- removing output and exiting \n");
     sprintf(cmd,"rm -f %s", Output_File_Name);
     if(system(cmd) != 0) assert(0);
     sprintf(cmd,"rm -rf %s", gatekeeperStorePath);
     if(system(cmd) != 0) assert(0);
     if(append)
       rename(tmpFilePath,gatekeeperStorePath);
   }
#endif
   
   
   exit(status != GATEKEEPER_SUCCESS);
}


/***********************************************************************************/
FILE *  File_Open
    (const char * Filename, const char * Mode, int exitOnFailure)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL && exitOnFailure)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }

