
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/fraguid2iid.c,v $
$Revision: 1.1.1.1 $
**********************************************************************/

//  Convert a list of fragment uids to iids using a fragment store

/*********************************************************************/
// headers
// standard headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include  <assert.h>
#include  <ctype.h>

// project headers
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_distStore.h"
#include "AS_ORA_fragments.h"
#include "AS_ORA_overlaps.h"
#include "AS_ORA_statistics.h"
#include "AS_ORA_inlines.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_PHash.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"




/*********************************************************************/


/*********************************************************************/
// defines
#ifndef MAX_SEQUENCE_LENGTH
#define MAX_SEQUENCE_LENGTH AS_READ_MAX_LEN
#endif

#ifndef MAX_SOURCE_LENGTH
#define MAX_SOURCE_LENGTH 512
#endif

#define PRINTUID /*output two fields: uid iid*/
/*********************************************************************/


/*********************************************************************/
// structures
/*********************************************************************/


/*********************************************************************/

#define  HASH_LOAD_FACTOR     0.3333
    //  Portion of hash table to be occupied by keys.
#define  MAX_LIST_LEN             1000000
    //  Length of insert list
#define  SKIP_MODULUS         11
    //  Modulus for secondary hash function to determine skip distance
    //  for linear probing.


static int  Cmp
    (const void * a, const void * b)
  {
   int64  a_val, b_val;

   a_val = * ((int64 *) a);
   b_val = * ((int64 *) b);

   if  (a_val < b_val)
       return  -1;
   else if  (b_val < a_val)
       return  1;
     else
       return  0;
  }
static int  Hash_Find
    (int64 key, int64 tab [], int64 size);
static int64  Next_Odd_Prime
    (int64 N);
static void  Read_Frag_IDs
    (FILE * fidfile, int64 * * hash_table, int64 * hash_table_size);
void  Read_Fragments
    (char * fragstore_name, char * gkpstore_name, int64 * hash_table, int64 hash_size);

/*********************************************************************/

int printUID=0;
int printSTATUS=0;

int main( int argc, char ** argv )
{
  char              * fragstore_name = NULL;
  char              * gkpstore_name = NULL;
  char              * uid_filename = NULL;
  int64  * hash_table, hash_table_size;
  FILE  * fidfile;
  
  // parse the command line parameters
  // use getopt(): see "man 3 getopt"
  {
    int ch, errflg = 0;
    optarg = NULL;
    while( !errflg && ((ch = getopt( argc, argv, "s:i:USg:" )) != EOF) )
    {
      switch( ch )
      {
        case 's':
          fragstore_name = optarg;
          break;
        case 'i':
          uid_filename = optarg;
          break;
        case 'U':
	  printUID=1;
	  break;
        case 'S':
          printSTATUS=1;
	  break;
        case 'g':
          gkpstore_name = optarg;
          break;
        case '?':
          fprintf( stderr, "Unrecognized option -%c\n", optopt );
        default:
          errflg++;
          break;
      }
    }

    // need fragstore_name & min_overlap and one or both of
    // input and output ovl filenames
    if( errflg != 0 || fragstore_name == NULL||
	  (gkpstore_name!=NULL&&printSTATUS==0)||
	(printSTATUS && gkpstore_name==NULL))
    {
      fprintf( stderr, "Usage: %s\n"
               "       -s fragstorename\n"
               "       -i uid-filename\n"
               "       [-U]                          (gives UID in output)\n"
               "       [-S -g gatekeeperstorename]   (gives status [1=deleted] of frag)\n",
               argv[0] );
    return 1;
    }
  }

  fprintf(stderr,"%s: WARNING:\n"
                 "   MAY RETURN MULTIPLE IIDs PER UID\n",argv[0]);

  if(!printSTATUS)
    fprintf(stderr,"  IT IS RECOMMENDED TO RUN WITH THE -S OPTION SO THAT\n"
                   "  DELETION STATUS OF FRAGMENTS IS OUTPUT!\n");
  else 
    assert(gkpstore_name!=NULL);


   fidfile = fopen (uid_filename, "r");
   if  (fidfile == NULL)
       {
        fprintf (stderr, "ERROR:  Can't open file \"%s\"\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }

   Read_Frag_IDs (fidfile, & hash_table, & hash_table_size);
   fprintf (stderr, "Hash table size = " F_S64 "\n", hash_table_size);
   fclose (fidfile);

  // read in the fragment store
  fprintf( stderr, "reading fragments...\n" );
  Read_Fragments (fragstore_name, gkpstore_name, hash_table, hash_table_size);

  fprintf( stderr, "Done.\n" );
  return 0;
}



void  Read_Fragments
    (char * fragstore_name, char * gkpstore_name, int64 * hash_table, int64 hash_size)

//  Read fragments from fragstore name  fragstore_name  and
//  output iid's of frags whose uid's are in the hash table.

{
  FragStoreHandle store_handle;
  FragStreamHandle stream_handle;
  ReadStructp read_struct;
  GateKeeperStore gkpStore;
  GateKeeperFragmentRecord gkpFrag;

  // open the frag store
  store_handle = openFragStore( fragstore_name, "r" );
  // NOTE: is this a valid check?
  if( store_handle < 0 )
  {
    fprintf( stderr,
             "Failed to open fragment store %s.\n",
             fragstore_name );
    return;
  }

  // get a stream handle to iterate through fragments
  stream_handle = openFragStream( store_handle, NULL, 0 );
  // NOTE: is this a valid check?
  if( stream_handle < 0 )
  {
    fprintf( stderr, "Failed to open fragment stream on store %s.\n",
             fragstore_name );
    closeFragStore( store_handle );
    return;
  }

  if(printSTATUS){
    InitGateKeeperStore(&gkpStore,gkpstore_name);
    OpenReadOnlyGateKeeperStore(&gkpStore);
    assert(TestOpenGateKeeperStore(&gkpStore) == TRUE);
  }

  // get a reusable read structure
  read_struct = new_ReadStruct();
  if( read_struct == NULL )
  {
    fprintf( stderr, "Failed to allocate read structure.\n" );
    closeFragStream( stream_handle );
    closeFragStore( store_handle );
    return;
  }

  // read in the fragments & copy pertinent data to the fragments array
  while( nextFragStream( stream_handle, read_struct, FRAG_S_FIXED ) )
  {
   uint64  frag_uid;

   getAccID_ReadStruct (read_struct, & frag_uid);
   if  (Hash_Find (frag_uid, hash_table, hash_size))
       {
        uint32  frag_iid;

        getReadIndex_ReadStruct(read_struct, & frag_iid);
	
	if(printUID){
	  printf (F_UID " ",frag_uid);
	}
	printf("%d", frag_iid);
	if(printSTATUS){
	  if(getGateKeeperFragmentStore(gkpStore.frgStore,frag_iid,&gkpFrag)!=0)
            assert(0);
	  printf(" %d",(gkpFrag).deleted);

	}
	printf("\n");
       }
  }
  // clean up
  delete_ReadStruct( read_struct );
  closeFragStream( stream_handle );
  closeFragStore( store_handle );

  return;
}



static int  Hash_Find
    (int64 key, int64 tab [], int64 size)

//  Return  TRUE  iff  key  occurs in hash table  tab [0 .. (size - 1)] .

  {
   int  skip, sub;

   sub = key % size;
   skip = 1 + (key % SKIP_MODULUS);

   while  (tab [sub] != -1 && tab [sub] != key)
     sub = (sub + skip) % size;

   return  (tab [sub] == key);
  }



static int64  Next_Odd_Prime
    (int64 N)

//  Return the first odd prime  >= N .

  {
   int64  Div, Last;

   if  (N % 2 == 0)
       N ++;
   while  (TRUE)
     {
      Last = (int64) (sqrt ((double) N));
      for  (Div = 3;  Div <= Last;  Div += 2)
        if  (N % Div == 0)
            break;
      if  (Div > Last)
          return  N;
      N += 2;
     }
  }



static void  Read_Frag_IDs
    (FILE * fp, int64 * * tab, int64 * size)

//  Read integers from  fp  and insert them into a hash table  (* tab) .
//  Set  (* size)  to the size of the hash table.  The hash function is
//  just  key % size .

  {
   int64  i, ct, key, skip, sub;

   for  (ct = 0;  fscanf (fp, F_S64, & key) != EOF;  ct ++)
     ;

   //fprintf (stderr, "ct = %d\n", ct);
   (* size) = Next_Odd_Prime ((int64) (ct / HASH_LOAD_FACTOR));
   (* tab) = (int64 *) malloc ((* size) * sizeof (int64));
   if  (* tab == NULL)
       {
        fprintf (stderr, "ERROR:  malloc failed\n");
        exit (EXIT_FAILURE);
       }

// fprintf (stderr, "(* size) = %ld\n", (* size));
   for  (i = 0;  i < (* size);  i ++)
     (* tab) [i] = -1;

   rewind (fp);

   while  (fscanf (fp, F_S64, & key) != EOF)
     {

      sub = key % (* size);
      skip = 1 + (key % SKIP_MODULUS);

      while  ((* tab) [sub] != -1 && (* tab) [sub] != key)
{
// fprintf (stderr, "try sub = %d\n", sub);
        sub = (sub + skip) % (* size);
}

// fprintf (stderr, "key = %ld  sub = %d\n", key, sub);
      if  ((* tab) [sub] == -1)
          (* tab) [sub] = key;
     }

   return;
  }



