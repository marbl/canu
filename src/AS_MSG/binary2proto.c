
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
/* $Id: binary2proto.c,v 1.2 2004-09-23 20:25:24 mcschatz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include "AS_global.h"

int main(int argc, char * argv[])
{
  GenericMesg *pmesg;
  MesgReader reader = InputFileType_AS(stdin);
  ProtoIOMode mode = AS_HUMAN_MODE;
  
  if(argc == 2){
    char c = argv[1][1];
    fprintf(stderr,"* argv[1] = %s %c\n", argv[1], c);
    if(c == 'H'){
      mode = AS_HUMAN_MODE;
      fprintf(stderr,"* Saw a H\n");
    }else if(c == 'D'){
      fprintf(stderr,"* Saw a D\n");
      mode = AS_DROS_MODE;
    }
  } 
  fprintf(stderr,"* mode is %c\n", mode);
  SetProtoMode_AS(mode);

  while(reader(stdin,&pmesg) != EOF) {
    WriteProtoMesg_AS(stdout,pmesg);
    fflush(stdout);
  }
  exit (0);
}
