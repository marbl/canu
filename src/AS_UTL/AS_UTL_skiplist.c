
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
/**********************************************************************
 Module:      AS_UTL_skiplist.c
 Description: Implementation of a sorted sequence data structure that
              supports lookup, insert and delete a (key,value) pair. 
	      Using access macros, one can define sorted sequences
	      with different value types. The key is always a double.
	      THe data structure is implemented a skiplist (Pugh 91)
	      which supports the above three operations in exspected time O(log n)
	      with an exponential tail estimate.
	      The creation of the skiplist takes a boolean parameter which determines
	      whether the skiplist frees the information its value field points to,
	      or whether the user does this.
 Assumptions: - All inserted keys are bigger than minf and smaller than pinf   
**********************************************************************/

static char CM_ID[] = "$Id: AS_UTL_skiplist.c,v 1.1.1.1 2004-04-14 13:53:45 catmandew Exp $";


/* 
   Implementation of dynamic skiplists supporting the operations
   INSERT, DELETE und LOOKUP.
*/


#include "AS_UTL_skiplist.h"
#include "unistd.h"
#include <float.h>

/* throws a fair coin */
static int coin(void)
{
  return(random()&01);
} 


/* initializes the random number generator */
static void init_random(int seed)
{ 
  time_t l = seed;
  if (l==0) time(&l);
  srandom(l);
}



static void clear_item(sl_item it)
{
  it->up   = NULL;
  it->down = NULL;
  it->succ = NULL;
  it->pred = NULL;
}


static sl_item add_item(sl_item x, SkipList *sl, int i)
{
  sl_item a,b,y,it1,it2;
  
  if(i == sl->no_of_levels)
    {
      a=new_item;
      assert(a != NULL);
      b=new_item;  
      assert(b != NULL);
	
      clear_item(a);
      clear_item(b);
      
      a->down = sl->head;
      b->down = sl->tail;
      sl->head->up = a;
      sl->tail->up = b;
      a->succ = b;
      b->pred = a;
      a->key = minf;
      b->key = pinf;
      
      sl->head = a;
      sl->tail = b;

      sl->no_of_levels++;		
    }
  
  y = new_item;      
  assert(y != NULL);

  clear_item(y);
  y->key  = x->key;
  y->down = x;
  x->up = y;
  
  it1 = x->pred;
  while(it1->up == NULL)
    it1 = it1->pred;
  it1 = it1->up;
  it2 = x->succ;
  while(it2->up == NULL)
    it2=it2->succ;
  it2=it2->up;
  
  it1->succ=y;
  y->pred=it1;
  it2->pred=y;
  y->succ=it2;
  return(y);
}		


/* main functions */
/******************/
/* deallocates the skiplist */

void Free_SL(SkipList *sl)
{
  sl_item it1,it2;
  
  it1=sl->head;
  sl->head=sl->head->down;
  while(it1 != NULL)
    {	
      it2=it1->succ;
      while(it2 != NULL)
	{
	  if( sl->freeData && it1->down == NULL && it1->value != NULL )
	    sl->free_value((void*) it1->value);
	  free(it1);
	  it1=it2;
	  it2=it1->succ;
	}
      free(it1);
      it1=sl->head;
      if(sl->head != NULL)
	sl->head=sl->head->down;
    }
  free(sl);
}


/* print function for small lists */

void Print_SL(SkipList *sl)
{
  sl_item it1,it2;
  
  printf("\n");
  it1 = sl->head;
  while(it1 != NULL)
    {	
      it2=it1;
      while(it2 != NULL )
	{	
	  if( it2->key != pinf && it2->key != minf)
	    printf("%3.2lf\t",it2->key);
	  it2=it2->succ;
	}
      printf("\n");
      it1 = it1->down;
    }
  printf("\nskiplist with %d elements and %d levels\n",sl->no_of_elements,sl->no_of_levels); 
}
 
 


SkipList *Create_SL(int fd, SLF free_value)
{
  SkipList *sl = (SkipList*) malloc(sizeof(struct sl));

  sl_item a,b;
  assert( sl != NULL);
  init_random(0);
  sl->free_value = free_value;

  sl->head       = new_item;
  assert(sl->head != NULL);

  clear_item(sl->head);
  sl->freeData   = fd;
  sl->no_of_elements = 0;
  sl->no_of_levels   = 2;
  sl->head->key = minf;

  sl->tail       = new_item;
  assert(sl->tail != NULL);

  clear_item(sl->tail);
  sl->tail->key = pinf;
  sl->tail->pred = sl->head;
  sl->head->succ = sl->tail;

  a = new_item;
  assert(a != NULL);
  b = new_item;
  assert(b != NULL);
  clear_item(a);
  clear_item(b);
  
  a->down = sl->head;
  b->down = sl->tail;
  a->key = minf;
  b->key = pinf;
  sl->head->up = a;
  sl->tail->up = b;
  a->succ = b;
  b->pred = a;

  sl->head = a;
  sl->tail = b;

  return(sl);
}




/* returns the element with the biggest key less than or 
   equal to key  */

sl_item Lookup_SL(keyType key, SkipList *sl)
{	
  sl_item it1;
  int found,fertig;
  
  it1=sl->head;
  fertig = FALSE;
  while(fertig == FALSE)
    {
      found = FALSE;
      while(found == FALSE)
	{	
	if(it1->succ->key > key)
	  found = TRUE;
	else
	  it1=it1->succ;
	}
      if(it1->down == NULL)
	fertig = TRUE;
      else
	it1=it1->down;
    }
  return(it1);
}	
 
 

/* inserts an element in the skiplist if it is not already present
   returns the freshly inserted or present element */

 
sl_item Insert_SL(keyType key, valueType value, SkipList * sl)
{	
  sl_item it1,a,newi;
  int m;
  int i;	
  
  a = Lookup_SL(key,sl);
  if(a->key == key)	
    return(a);
  //else
    {
      newi=new_item;
      assert(newi != NULL);

      clear_item(newi);
      sl->no_of_elements++;
      
      newi->key=key;
      newi->value=value;
      
      newi->succ=a->succ;
      a->succ->pred=newi;
      newi->pred=a;
      a->succ=newi;
      
      it1 = newi;
      i = 2;
      m = coin();
      while(m == TRUE)
	{
	  it1=add_item(it1,sl,i++);
	  m = coin();
	}
      return(newi);
    }
}
 



/* deletes the element with key key from the skiplist  */

void Delete_SL(keyType key, SkipList *sl)
{	
  sl_item a,it;
  
  a = Lookup_SL(key,sl);
  if(a->key == key)
    {
      while(a->up != NULL)
	a=a->up;
      while(a != NULL)
	{	 
	  a->pred->succ = a->succ;
	  a->succ->pred = a->pred;

	  if(a->pred->key == minf  && a->succ->key == pinf)
	    {	
	      sl->head = sl->head->down;
	      sl->tail = sl->tail->down;
	      //    sl->head->succ = sl->tail;
	      //	      sl->tail->pred = sl->head;
	      free(sl->tail->up);
	      free(sl->head->up);
	      sl->head->up = NULL;
	      sl->tail->up = NULL;
	      sl->no_of_levels--;
	    }

	  it = a->down; 
 	  if( sl->freeData && a->down == NULL )
	    sl->free_value( (void*) a->value);
	  free(a);
	  a = it;
	}
      sl->no_of_elements--;
    }
}	


sl_item Min_SL(SkipList *sl)
{
  sl_item it = sl->head;
  while( it->down != NULL )
    it = it->down;
  return it->succ;
}


sl_item Max_SL(SkipList *sl)
{
  sl_item it = sl->tail;
  while( it->down != NULL )
    it = it->down;
  return it->pred;
}


int GetNum_SL(SkipList* sl){
  return sl->no_of_elements;
}

 
 
