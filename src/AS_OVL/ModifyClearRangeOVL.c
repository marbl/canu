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
* Module:  ModifyClearRangeOVL.c
* Description:
*   Change the clear range value of reads in specified frag Store
*   to the values specified in the input file
* 
*    Programmer:  A. Delcher
*       Started:  21 Sep 2004
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: ModifyClearRangeOVL.c,v 1.2 2005-03-22 19:49:19 jason_miller Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: ModifyClearRangeOVL.c,v 1.2 2005-03-22 19:49:19 jason_miller Exp $";


//  System include files

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>


//  Local include files

#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_UTL_version.h"


//  Type definitions



//  Static Globals

static char  * Frag_Store_Path = NULL;
    // Location of fragement store to be updated
static char  * Clear_Mod_Path = NULL;
    // File containing new clear ranges



//  Static Functions

static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FragStoreHandle  frag_store;
     // fragment store to be modified
   FILE  * fp;
   ReadStructp  fsread;
   time_t  now;
   unsigned int  fiid, old_lo, old_hi, new_lo, new_hi;
   int  lo, hi;
   int  status;

   now = time (NULL);
   fprintf (stderr, "### Starting at %s\n", ctime (& now));

   Parse_Command_Line  (argc, argv);

   frag_store = openFragStore (Frag_Store_Path, "rw+");
   fp = File_Open (Clear_Mod_Path, "r");
   fsread = new_ReadStruct ();

   while  (fscanf (fp, "%u %d %d", & fiid, & lo, & hi) == 3)
     {
      if  (lo == 0 && hi == 0)
          continue;

      if  (lo < 0)
          lo = 0;
      getFragStore (frag_store, fiid, FRAG_S_ALL, fsread);   
      getClearRegion_ReadStruct (fsread, & old_lo, & old_hi,
           READSTRUCT_ORIGINAL);
      new_lo = old_lo + lo;
      new_hi = old_lo + hi;
      setClearRegion_ReadStruct (fsread, new_lo, new_hi, READSTRUCT_OVL);
      status = setFragStore (frag_store, fiid, fsread);
     }

   fclose (fp);
   closeFragStore (frag_store);

   now = time (NULL);
   fprintf (stderr, "### Finished at %s\n", ctime (& now));
   
   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "h")) != EOF))
     switch  (ch)
       {
        case  'h' :
          Usage (argv [0]);
          exit (EXIT_FAILURE);

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);
          // fall through

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 2)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Frag_Store_Path = argv [optind ++];

   Clear_Mod_Path = argv [optind ++];

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s <FrgStore> <ClrMods>\n"
       "\n"
       "Change the clear range value of reads in <FrgStore> to\n"
       "to the values specified in <ClrMods>\n"
       "Format of entries in <ClrMods> is:\n"
       "<iid> <clr-begin> <clr-end>\n"
       "\n"
       "Options:\n"
       "-h  print this message\n",
       command);

   return;
  }
