
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
/* $Id: hash.c,v 1.1.1.1 2004-04-14 13:53:57 catmandew Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Htag
{ struct Htag *next;
  char        *value;
  int          id;
} HashEntry;

typedef struct
{ int         size, count;
  char       *lasthash;
  HashEntry **vector;
} HashTable;

static int HashKey(char *entry)
{ int i, key, glob;

  key = 0;
  glob = 0;
  for (i = 0; entry[i] != '\0'; i++)
    { glob = (glob << 8) | entry[i];
      if (i % 3 == 2)
        { key = key ^ glob;
          glob = 0;
        }
    }
  key = key ^ glob;
  return (key);
}

void PrintHashTable(HashTable *table)
{ int i;
  HashEntry *c;

  printf("\nHASH TABLE (%d,%d):\n",table->count,table->size);
  for (i = 0; i < table->size; i++)
    if (table->vector[i] != NULL)
      { printf("  Vector %4d:\n",i);
        for (c = table->vector[i]; c != NULL; c = c->next)
          printf("    %4d: '%s'\n",c->id,c->value);
      }
}

HashTable *NewHashTable(int size)
{ HashTable *table;
  int i;

  table = (HashTable *) malloc(sizeof(HashTable));
  table->vector   = (HashEntry **) malloc(sizeof(HashEntry *)*size);
  table->size     = size;
  table->count    = 0;
  table->lasthash = "";
  for (i = 0; i < size; i++)
    table->vector[i] = NULL;
  return (table);
}

static void DoubleHashTable(HashTable *table)
{ HashTable *newone;
  HashEntry *c, *d;
  int i, key;

  newone = NewHashTable(2*table->size-1);
  for (i = 0; i < table->size; i++)
    for (c = table->vector[i]; c != NULL; c = d)
      { key = HashKey(c->value) % newone->size;
        d = c->next;
        c->next = newone->vector[key];
        newone->vector[key] = c;
      }
  free(table->vector);
  table->size = newone->size;
  table->vector = newone->vector;
  free(newone);
}

char *GetLastHash(HashTable *table)
{ return (table->lasthash); }

int HashLookup(char *entry, HashTable *table, int *isdef)
{ int key;
  HashEntry *chain;

  key = HashKey(entry) % table->size;
  chain = table->vector[key];
  while (chain != NULL)
    { if (strcmp(chain->value,entry) == 0)
        { table->lasthash = chain->value;
          if (chain->id >= 0)
            { if (*isdef)
                *isdef = -1;
              else
                *isdef = 1;
              return (chain->id);
            }
          else if (*isdef)
            { chain->id = -(chain->id+1);
              return (chain->id);
            }
          else
            return (-(chain->id+1));
        }
      chain = chain->next;
    }
  return (-1);
}
    
int HashAdd(char *entry, HashTable *table, int isdef)
{ int key;
  HashEntry *chain;

  key = HashKey(entry) % table->size;
  chain = table->vector[key];
  while (chain != NULL)
    { if (strcmp(chain->value,entry) == 0)
        return (-1);
      chain = chain->next;
    }
  chain = (HashEntry *) malloc(sizeof(HashEntry));
  chain->next = table->vector[key];
  table->vector[key] = chain;
  if (isdef)
    chain->id = table->count;
  else
    chain->id = -(table->count+1);
  table->count += 1;
  chain->value = strdup(entry);
  table->lasthash = chain->value;
  if (chain->id > table->size*.25)
    DoubleHashTable(table);
  return (table->count-1);
}

char *HashRefButNoDef(HashTable *table)
{ int i;
  HashEntry *c;

  for (i = 0; i < table->size; i++)
    for (c = table->vector[i]; c != NULL; c = c->next)
      if (c->id < 0)
        return (c->value);
  return (NULL);
}
