
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

static char CM_ID[] = "$Id: AS_UTL_skiplist_test.c,v 1.6 2008-06-27 06:29:21 brianwalenz Exp $";


#include "AS_UTL_skiplist.h"


typedef struct{
  int hi;
  double how;
  char are;
  int* you;
} Complicated;


Complicated* alloc_Complicated(int len)
{
  Complicated *c= (Complicated*)safe_malloc(sizeof(Complicated));
  c->you = (int*)safe_calloc(sizeof(int),len);
  return c;
}

void free_Complicated(Complicated* c){
  safe_free(c->you);
  safe_free(c);
}

SL_DEF(int)
SL_DEF(Complicated)


int main(void){
  int i;
  sl_item it;
  double x;
  Complicated *com;
  SL_TYPE(int) *sl = CreateSL_int(FALSE,NULL);
  SL_TYPE(Complicated) *sl2 = CreateSL_Complicated(TRUE,(SLF) free_Complicated);
  int dummy=4711;
  double lu = (double)GetRand_AS(0,1000,TRUE);



  for(i=0; i<=4; i++){
    x = (double)GetRand_AS(0,1000,TRUE);
    InsertSL_int(x,&dummy,sl);
  }
  PrintSL_int(sl);


  printf("looking up %lf \n",lu);
  it = LookupSL_int(lu,sl);
  printf("result of lookup = %lf\n",it->key);

  printf("inserting %lf \n",lu);
  it = InsertSL_int(lu,&dummy,sl);
  PrintSL_int(sl);

  printf("deleting %lf \n",lu);
  DeleteSL_int(lu,sl);
  PrintSL_int(sl);
  FreeSL_int(sl);



  for(i=0; i<=4; i++){
    x = (double)GetRand_AS(0,1000,TRUE);
    com = (Complicated*) alloc_Complicated(10);
    InsertSL_Complicated(x,com,sl2);
  }
  PrintSL_Complicated(sl2);


  printf("looking up %lf \n",lu);
  it = LookupSL_Complicated(lu,sl2);
  printf("result of lookup = %lf\n",it->key);

  printf("inserting %lf \n",lu);
  com = (Complicated*)safe_malloc(sizeof(Complicated));
  it = InsertSL_Complicated(lu,com,sl2);
  PrintSL_Complicated(sl2);

  printf("deleting %lf \n",lu);
  DeleteSL_Complicated(lu,sl2);
  PrintSL_Complicated(sl2);
  FreeSL_Complicated(sl2);

  return 0;
}
