
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
#include <stdlib.h> 
#include <stdio.h> 
#include <assert.h>
#include <unistd.h> /* man 3 getopt */

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"

int main(int argc, char *argv[]){

#if 0
 FragStreamHandle jane;
#endif
 FragStoreHandle source;
 ReadStructp myRead;
 int load = FALSE;
 int32 begin = -1, end = -1;
 int i;
 int ch;
 if(argc  < 2 ){
  fprintf(stderr,"Usage: %s [-b <firstElem>] [-e <lastElem>] [-l] <StorePath1> \n",
	  argv[0]);
  fprintf(stderr,"   -l option causes frag store to be loaded into memory, rather than opened\n");
  exit(1);
 }

    while ((ch = getopt(argc, argv, "b:e:l")) != EOF){
      switch(ch) {
      case 'l':
	load = TRUE;
	break;
      case 'e':
	end = atoi(optarg);
	fprintf(stderr,"* end = %d\n", end);
	break;
      case 'b':
	begin = atoi(optarg);
	fprintf(stderr,"* begin = %d\n", begin);
	break;
      default:
	fprintf(stderr,"* Unknown option %s\n", optarg);
	break;
      }
    }

 if(load){
 fprintf(stdout,"* LOADing FragStore %s\n", argv[optind]);
    source = loadFragStorePartial(argv[optind],STREAM_FROMSTART, STREAM_UNTILEND);
 }else{
 fprintf(stdout,"* Opening FragStore %s\n", argv[optind]);
    source = openFragStore(argv[optind],"r");
 }
 if(source == NULLSTOREHANDLE){
   exit(1);
 }

#if 0
   {
     DistStore distStore;
     char distStoreName[1024];
     DistRecord distRecord;
     StreamHandle distStream;
     StoreStat stats;
     int i, start;

     sprintf(distStoreName,"%s/db.dst", argv[1]);
     fprintf(stdout,"* Opening DistStore %s\n", distStoreName);
     distStore = openDistStore(distStoreName,"r"); 
     distStream = openStream(distStore,NULL,0);
     if(0 != statsStore(distStore, &stats))
       assert(0);
     
     fprintf(stdout,"* Distance Records (" F_S64 "-" F_S64 ")\n\n",
	     stats.firstElem, stats.lastElem);
#if 0
     i = getStartIndexStream(distStream);
     while(nextStream(distStream,&distRecord) ){
       fprintf(stdout,"* Dist Record %d d:%d iid:" F_IID " "
	       "uid:" F_UID " %f +- %f\n",
	       i, distRecord.deleted, distRecord.IID, distRecord.UID,
	       distRecord.mean, distRecord.stddev);
     }
#else   // Just an example
     {
       for(i = stats.firstElem;i <= stats.lastElem ; i++){
	int ret =  getDistStore(distStore,i,&distRecord);
	
       fprintf(stdout,"* Dist Record %d d:%d iid:" F_IID " uid:" F_UID " %f +- %f\n",
	       i, distRecord.deleted, distRecord.IID, distRecord.UID,
	       distRecord.mean, distRecord.stddev);
       }
     }	 
#endif
     closeDistStore(distStore);
   }

 fprintf(stdout,"* Closed DistStore \n\n\n");
#endif

 myRead =  new_ReadStruct();

 fprintf(stdout,"* Dumping fragStore %s (%d,%d) of (" F_S64 "," F_S64 ")\n",
	 argv[optind],begin,end,
         getFirstElemFragStore(source), getLastElemFragStore(source));

#if 0
 jane = openFragStream(source,NULL,0);

   while(nextFragStream(jane, myRead, FRAG_S_ALL)){
     dump_ReadStruct(myRead, stdout);
   }
   closeFragStream(jane);
   fprintf(stdout,"* Closing jane, source\n");
#else
   if(begin < getFirstElemFragStore(source) || begin > getLastElemFragStore(source)){
     begin = getFirstElemFragStore(source);
   }
   if(end > getLastElemFragStore(source) || end < getFirstElemFragStore(source)){
     end = getLastElemFragStore(source);
   }
   for(i = begin; i <= end; i++){
     getFragStore(source, i, FRAG_S_ALL, myRead);
     dump_ReadStruct(myRead, stdout, FALSE);
   }

#endif
   fprintf(stdout,"* Bye Bye\n");

     
     


   exit(0);
}
