
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
#include "strings.h"
#include "AS_PER_ReadStruct.h" 
#include "AS_PER_fragStore.h" 
#include "AS_PER_genericStore.h"
#include "AS_PER_distStore.h"

char InputName[2048];

int main(int argc, char *argv[]){


 FragStoreHandle source;
 ReadStructp myRead;
 int i, j, k = 0;
 int numFiles = 1;
 int64 numFragmentsPerFile;
 int64 first, last;
 char buffer[2048];
 // char dummy[1000];  // Don't know why, but this program doesn't work if you delete this!

 FILE *outfp;

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
		       "n:")) != EOF)){
      switch(ch) {
      case 'n':
	numFiles = atoi(optarg);
	fprintf(stderr,"* numFiles set to %d\n", numFiles);
	break;
      case '?':
	fprintf(stderr,"Unrecognized option -%c",optopt);
      default :
	errflg++;
      }
    }
  }

  fprintf(stderr,"* Opening frag store %s\n", argv[optind]);

  source = -1;
  source = openFragStore(argv[optind],"r");

  strcpy(InputName, argv[optind]);
 sprintf(buffer,"%s.dump.0", InputName);
 outfp = fopen(buffer,"w");
 AssertPtr(outfp);

 myRead =  new_ReadStruct();

 first = getFirstElemFragStore(source);
 last = getLastElemFragStore(source);

 numFragmentsPerFile = (last - first)/numFiles + 1;
 fprintf(stdout,"* Dumping fragStore %s (" F_S64 "," F_S64 ") in " F_S64 " frag batches\n",
	 InputName,first, last, numFragmentsPerFile);



 for(j = 0, i = first; i <= last; i++, j++){
#if 1
   if(j == numFragmentsPerFile){
     //     fprintf(stderr,"* j == " F_S64 "\n", numFragmentsPerFile);
     fflush(NULL);
     j = 0;
#if 1
     fclose(outfp);
     sprintf(buffer,"%s.dump.%d", InputName,++k);
     outfp = fopen(buffer,"w");
     AssertPtr(outfp);
     fprintf(stderr,"* Just opened %s\n", buffer);
#endif
   }
#endif
   getNDumpFragStore(source, i, &myRead, outfp);
   if(j%1000 == 0)
     fflush(NULL);
 }


 fclose(outfp);

   exit(0);
  }
