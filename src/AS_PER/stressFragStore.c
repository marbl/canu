
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
#include <unistd.h>

#include "AS_global.h" 
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"
#include "AS_UTL_timer.h"
#include "AS_UTL_rand.h"
#include "AS_PER_fragDB.h"

#define BLOCKSIZE 1
#define NUMREADS 10000

int main(int argc, char *argv[]){

  int seed = 37;
  int useDB = 0;
  int numReads = NUMREADS;
  int blockSize = BLOCKSIZE;
  int random = FALSE;

 FragStoreHandle source;
 FragStreamHandle jane;
 ReadStructp myRead;
 int firstElem, lastElem;
 int i;
 TimerT timer;
 time_t t;
   double total_time;
   int64 cycles;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
									"dDrRn:b:s:")) != EOF)){
      switch(ch) {
		case 's':
		  seed = atoi(optarg);
		  fprintf(stderr,"* Seed set to %d \n",seed);
		  break;
		case 'd':
		  useDB = 0;
		  fprintf(stderr,"* Using fragStore (default)\n");
		  break;
		
		case 'D':
		  useDB = 1;
		  fprintf(stderr,"* Using DB\n");
		  break;
		case 'r':
		  random = 0;
		  fprintf(stderr,"* Sequential test\n");
		  break;
		case 'R':
		  random = 1;
		  fprintf(stderr,"* Random test (default)\n");
		  break;
		case 'n':
		  numReads = atoi(optarg);
		  fprintf(stderr,"* Number of reads set to %d\n", numReads);
		  break;
		case 'b':
		  blockSize = atoi(optarg);
		  fprintf(stderr,"* Blocksize (for sequential reads) set to %d\n", blockSize);
		  break;
  
		default:
		  fprintf(stderr, "Usage: %s -[rRn:b:] <StorePath1> \n",  argv[0]);
		  exit(1);

	
	  }
	  
	}
	
  }
  
 if(argc  < 2 ){
  exit(1);
 }
 t = time(0);
 
 if(!useDB)
 {
   fprintf(stderr,"* Testing using fragStore\n");
   source = openFragStore(argv[optind],"r");
   
   firstElem = getFirstElemFragStore(source);
   lastElem = getLastElemFragStore(source);
 }else
 {
   fprintf(stderr,"* Testing using db\n");
   firstElem = 1;
   lastElem = 3600000;
 }
 
 
 if(random){
 fprintf(stderr, "* Random Stress Test  started at %s\n", ctime(&t));
 }else{
 fprintf(stderr, "* Ordered Stress Test started at %s\n", ctime(&t));
 }


 fprintf(stderr, "about to InitTimerT in main\n");

 InitTimerT(&timer);
 InitRandom_AS((long)seed);
 myRead =  new_ReadStruct();
 
 if(random){
   for(i = 0; i < numReads; i++){
     int index = GetRand_AS(firstElem, lastElem, TRUE);
     StartTimerT(&timer);
	 if(useDB)
	 {
	   getFragDB(index,FRAG_S_ALL, myRead, blockSize);
	   //dump_ReadStruct(myRead, stderr);
	 }
	 else
	 {
	   getFragStore(source, index, FRAG_S_ALL, myRead);
	 }
	 
     StopTimerT(&timer);
     
     if(i%(numReads/10) == 0){
       t = time(0);
       fprintf(stderr,"* After %d reads at time %s \n", i, ctime(&t));
     }
   }
   total_time = TotalTimerT(&timer, &cycles);
   fprintf(stdout, "stressFragStore did %ld random reads in %g seconds (avg %g)\n",
	   cycles, total_time, total_time/cycles);

 }else{
   int numUnitigs = numReads/blockSize;
   fprintf(stdout,"*Simulated read of %d unitigs of size %d\n",
	   numUnitigs, blockSize);
   for(i = 0; i < numUnitigs; i++){
     int index = GetRand_AS(firstElem, lastElem, TRUE);
     int j;

	 if (seed == -1)
	   index = 1;
	 
     for(j = 0; j < blockSize; j++){
       
     StartTimerT(&timer);
	 if(useDB)
	 {
	   getFragDB(index+j,FRAG_S_ALL, myRead, blockSize);
//	   dump_ReadStruct(myRead, stderr);
	 }
	 else
	 {
	   getFragStore(source, index+j, FRAG_S_ALL, myRead);
	 }
	 
     StopTimerT(&timer);
     if((i*blockSize + j)%(numReads/10) == 0){
       t = time(0);
       fprintf(stderr,"* After %d reads at time %s \n", i*blockSize + j, ctime(&t));
     }
     }
   }
   total_time = TotalTimerT(&timer, &cycles);
   fprintf(stdout, "stressFragStore did %ld block reads in %g seconds (avg %g)\n",
	   cycles, total_time, total_time/numReads);

 }
   fprintf(stdout,"* Bye Bye\n");

   exit(0);
}
