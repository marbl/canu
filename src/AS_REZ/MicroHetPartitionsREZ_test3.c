
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
static char CM_ID[] = "$Id: MicroHetPartitionsREZ_test3.c,v 1.3 2005-03-22 19:07:44 jason_miller Exp $";

#include "AS_UTL_skiplist.h"
#include "UtilsREZ.h"
#include "MicroHetScoreREZ_test3.h"
#include "MicroHetPartitionsREZ_test3.h"

#define DEBUG 2

SL_DEF(Partition_t)

void AS_REZ_free_marker(Marker_t *m)
{
  free(m->set);
  free(m);
}

/* functions to allocate and free a marker of size l */
Marker_t *AS_REZ_allocate_marker(int l)
{
  int i;
  Marker_t* m = (Marker_t*) safe_malloc(sizeof(Marker_t));
  m->set      = (int*) safe_calloc(sizeof(int),l);
  for(i=0; i<l; i++)
    m->set[i] = TRUE;
  m->len      = l;
  return m;
}

void AS_REZ_print_marker(Marker_t *m){
  int i;
  for(i=0; i<m->len; i++)
    printf("|%d|",m->set[i]);
  printf("\n");
}

