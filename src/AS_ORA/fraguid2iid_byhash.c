
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
$Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_ORA/Attic/fraguid2iid_byhash.c,v $
$Revision: 1.1 $
**********************************************************************/

//  Convert a list of fragment uids to iids using a gkp store

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

static  GateKeeperStore my_gkp_store ; // See AS_PER_gkpStore.h

int uid2iid(uint64 uid){
  PHashValue_AS value;
  static firstFailure=1;
  if(HASH_FAILURE == LookupInPHashTable_AS(my_gkp_store.hashTable, 
					       UID_NAMESPACE_AS,
					       uid,
					       &value)){
    if(firstFailure){
      fprintf(stderr,"Tried to look up iid of unknown uid: " F_UID "; this may reflect trying to use a deleted fragment; further instances will not be reported.\n",uid);
      firstFailure=0;
    }
    return (-1);
  }
  return (value.IID);
}


/*********************************************************************/

int main( int argc, char ** argv )
{
  char              * gkpstore_name = NULL;
  char              * uid_filename = NULL;
  FILE  * fidfile;
  uint64 uid;

  int printUID=0;
  int printSTATUS=0;

  // parse the command line parameters
  // use getopt(): see "man 3 getopt"
  {
    int ch, errflg = 0;
    optarg = NULL;
    while( !errflg && ((ch = getopt( argc, argv, "g:i:SU" )) != EOF) )
    {
      switch( ch )
      {
        case 'i':
          uid_filename = optarg;
          break;
        case 'U':
	  printUID=1;
	  break;
        case 'g':
          gkpstore_name = optarg;
          break;
        case 'S':
          printSTATUS=1;
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
    if( errflg != 0 ||
	(printSTATUS && gkpstore_name==NULL))
    {
      fprintf( stderr, "Usage: %s\n"
               "       -g gatekeeperstorename\n"
               "       -i uidlist-filename\n"
               "       [-U]                          (gives UID in output)\n"
               "       [-S]   (gives status [1=deleted] of frag)\n",
               argv[0] );
    return 1;
    }
  }

  assert(gkpstore_name!=NULL);
  InitGateKeeperStore(&my_gkp_store, gkpstore_name);
  OpenReadOnlyGateKeeperStore(&my_gkp_store);

  if(!printSTATUS)
    fprintf(stderr,"  IT IS RECOMMENDED TO RUN WITH THE -S OPTION SO THAT\n"
                   "  DELETION STATUS OF FRAGMENTS IS OUTPUT!\n");

   fidfile = fopen (uid_filename, "r");
   if  (fidfile == NULL)
       {
        fprintf (stderr, "ERROR:  Can't open file \"%s\"\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }

   while(fscanf(fidfile,F_UID,&uid)==1){
     int32 iid = uid2iid(uid);
     if(printUID){
       printf(F_UID "\t",uid);
     }
     printf(F_IID,iid);
     if(printSTATUS&&iid>0){
       GateKeeperFragmentRecord gkpFrag;
       if(getGateKeeperFragmentStore(my_gkp_store.frgStore,iid,&gkpFrag)!=0)
	 assert(0);
       printf(" %d",(gkpFrag).deleted);
     }
     printf("\n");
   }

  fprintf( stderr, "Done.\n" );
   fclose (fidfile);
  return 0;
}




