
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

/* Read a file from stdin, throw away all but the ILK messages.
   Sort teh ILK messages.
   Output them to stdout.
*/


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
#include "AS_UTL_Var.h"

VA_DEF(InternalLinkMesg)

VA_TYPE(InternalLinkMesg) *Links;

static int CompareLinks (const void *c1, const void *c2){
  InternalLinkMesg *l1 = (InternalLinkMesg *)c1;
  InternalLinkMesg *l2 = (InternalLinkMesg *)c2;
  int diff;

  diff = l1->ifrag1 - l2->ifrag1;
  if(diff)
    return diff;

  diff = l1->ifrag2 - l2->ifrag2;
  if(diff)
    return diff;


  return 0;
}


int main(void){

  MesgReader reader = InputFileType_AS(stdin);
  MesgWriter writer = OutputFileType_AS(AS_PROTO_OUTPUT);
  GenericMesg *pmesg;
  InternalLinkMesg *links;
  CDS_CID_t tmp;

  Links = CreateVA_InternalLinkMesg(100);
  while(  (EOF != (reader)(stdin, &pmesg))){
    switch(pmesg->t){
    case MESG_ILK:
      links = pmesg->m;
      if(links->ifrag1 > links->ifrag2){
	tmp = links->ifrag1;
	links->ifrag1 = links->ifrag2;
	links->ifrag2 = tmp;
      }
      AppendInternalLinkMesg(Links, (InternalLinkMesg *)pmesg->m);
      break;
    default:
      break;
    }
  }

  fprintf(stderr,"* Read %d ILK messages \n",
	  (int) GetNumInternalLinkMesgs(Links));

  links = GetInternalLinkMesg(Links,0);
  qsort(links, GetNumInternalLinkMesgs(Links), sizeof(InternalLinkMesg), CompareLinks);

  {
    int i;
    GenericMesg outMesg;

    outMesg.t = MESG_ILK;
    for(i = 0; i < GetNumInternalLinkMesgs(Links); i++){
      outMesg.m = GetInternalLinkMesg(Links,i);
      writer(stdout, &outMesg);
    }
  }
  return 0;
}
