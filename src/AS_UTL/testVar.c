
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

#include <assert.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"

typedef struct{
  int x;
  int y;
  int z;
}GorkT;

VA_DEF(GorkT)


int main(int argc, char **argv){
  int i;
  GorkT gork;
  long int length;
  VA_TYPE(GorkT) *gorks;
  if(argc < 2)
    length = 10000;
  else
    length = atoi(argv[1]);

  gorks = CreateVA_GorkT(length);

  gork.x = 0;
  gork.y = 0;
  gork.z = 0;

  SetGorkT(gorks,length/2, &gork);

  for(i = 0; i < GetNumGorkTs(gorks); i++){
    GorkT *gp = GetGorkT(gorks, i);
    assert(gp->x == 0 && gp->y == 0 && gp->z == 0);
  }

  SetGorkT(gorks,length, &gork);

  for(i = 0; i < GetNumGorkTs(gorks); i++){
    GorkT *gp = GetGorkT(gorks, i);
    assert(gp->x == 0 && gp->y == 0 && gp->z == 0);
  }
  return 0;
}

