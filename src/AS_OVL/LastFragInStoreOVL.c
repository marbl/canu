
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
/*************************************************
* Module:  LastFragInStore.c
* Description:
*   Prints the internal ID of the last fragment in the fragment store
*   specified on the command line.
* 
*    Programmer:  A. Delcher
*       Written:  10 Aug 2000
*  Last Revised:  
* 
*************************************************/

/* RCS info
 * $Id: LastFragInStoreOVL.c,v 1.4 2005-03-22 19:49:19 jason_miller Exp $
 * $Revision: 1.4 $
*/

static char  CM_ID []
    = "$Id: LastFragInStoreOVL.c,v 1.4 2005-03-22 19:49:19 jason_miller Exp $";


/*************************************************************************/
/* System include files */
/*************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>


/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"


/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/

static void  Usage
    (char * executable_name);



int  main
    (int argc, char * argv [])

  {
   char  * fragstore_path = NULL;
   int64  last_stored_frag;


   // Parse the argument list using "man 3 getopt".

   {
    int  ch, errflg = 0;
    optarg = NULL;

    while  (! errflg && ((ch = getopt (argc, argv, "h")) != EOF))
      switch  (ch)
        {
         case  'h' :
           Usage (argv [0]);
           exit (EXIT_FAILURE);
           break;
         default :
           errflg ++;
        }

    if  (argc - optind != 1)
        {
         Usage (argv [0]);
         exit (EXIT_FAILURE);
        }

    fragstore_path = argv [optind ++];
   }

   if  (existsFragStore (fragstore_path))
       {
        FragStore  frag_store = openFragStore (fragstore_path, "r");

        last_stored_frag = getLastElemFragStore (frag_store);
        closeFragStore (frag_store);
       }
     else
       {
        fprintf (stderr, "ERROR:  Could not open fragstore \"%s\"\n",
                 fragstore_path);
        exit (EXIT_FAILURE);
       }

   printf ("Last frag in store is iid = " F_S64 "\n", last_stored_frag);

   return  0;
  }



static void  Usage
    (char * executable_name)

//  Print parameters and options to run this program, named
//   executable_name .

  {
   fprintf (stderr,
           "USAGE:  %s [-h] <FragStore>\n"
           "Prints iid of last frag in store\n"
           "Options:\n"
           "-h   print this message\n",
           executable_name);
  }

