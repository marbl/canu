
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
/* CountMessages
 *    Extracts messages of a given set of types
 *
 * $Id: ExtractMessages.c,v 1.3 2005-03-22 19:06:20 jason_miller Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "AS_global.h"
#include <assert.h>

#define MAX_MESG NUM_OF_REC_TYPES
int main(int argc, char **argv)
{ 
  GenericMesg *pmesg;
  int include[MAX_MESG + 1];
  int i;
  int Include = TRUE;
  int numTypesExtracted = 0;
  MesgReader reader = InputFileType_AS(stdin);

  if(argc < 2){
    fprintf(stderr,"* usage: extractmessages [-x] [list of message types] < <input file>\n");
    fprintf(stderr,"*        e.g., extractmessages IUM < a004.cgb \n");
    exit(1);
  }

  if(argv[1][0] == '-' ){
    if( argv[1][1] == 'x'){
      Include = FALSE;
    }else{
      fprintf(stderr,"* Illegal option %s...ignored\n",
	      argv[1]);
    }
  }
  for(i = 0; i <= MAX_MESG; i++)
    include[i] = !Include;

  for(i = 1; i < argc; i++){
    int type;

    if(argv[i][0] == '-')
      continue;

      type = GetMessageType(argv[i]);
    if(type >= 1 && type <= MAX_MESG){
      include[type] = Include;
      fprintf(stderr,"* %s messages of type %d %s\n",
	      (Include?"Including":"Excluding"),type, GetMessageName(type));
      numTypesExtracted++;
    }else{
      fprintf(stderr,"* Unknown message type %s...ignoring\n",
	      argv[i]);
    }
  }
  if(numTypesExtracted == 0){
    fprintf(stderr,"*** No message types specified...exiting\n");
    exit(1);
  }
 while (reader(stdin,&pmesg) != EOF){
   //   fprintf(stderr,"* MAX_MESG %d pmesg->t %d\n",
   //	   MAX_MESG, pmesg->t);
   assert(pmesg->t <= MAX_MESG );
   if(include[pmesg->t]){
    WriteProtoMesg_AS(stdout,pmesg);
    fflush(stdout);
   }
 }
 exit(0);
}

