
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
* Module:  iid2fastaOVL.c
* Description:
*   Take a list of fragment iid's (from stdin) and print their
*   clear range sequences in multifasta format (to stdout)
* 
*    Programmer:  A. Delcher
*       Written:  13 Jun 2000
* 
* 
* Assumptions:
*  argv [1] is the fragment store in which fragments are looked up
*
*************************************************/

/* RCS info
 * $Id: iid2fastaOVL.c,v 1.3 2005-03-22 19:07:13 jason_miller Exp $
 * $Revision: 1.3 $
*/


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
#include  "AS_UTL_version.h"




int  main
    (int argc, char * argv [])

  {
   FragStore  store;
   static char sequence [AS_READ_MAX_LEN + 1];
   static char quality [AS_READ_MAX_LEN + 1];
   ReadStructp  fsread = new_ReadStruct ();
   uint32  clear_begin, clear_end;
   int  iid, len;
   char  * p, * s;

   if  (argc < 2)
       {
        fprintf (stderr, "USAGE:  %s <fragstore>\n",
                 argv [0]);
        exit (EXIT_FAILURE);
       }

   store = openFragStore (argv [1], "r");
   if  (store == NULLSTOREHANDLE)
       {
        fprintf (stderr, "ERROR:  Can't open store \"%s\"\n",
                 argv [1]);
        exit (EXIT_FAILURE);
       }

   while  (scanf ("%d", & iid) == 1)
     {
      getFragStore (store, iid, FRAG_S_SEQUENCE, fsread);
      getSequence_ReadStruct (fsread, sequence, quality, AS_READ_MAX_LEN);
      getClearRegion_ReadStruct (fsread, & clear_begin, & clear_end, 
				 READSTRUCT_LATEST);

      len = clear_end - clear_begin;
      assert (len <= AS_READ_MAX_LEN);

      p = sequence + clear_begin;
      p [len] = '\0';

      printf (">iid=%d\n", iid);
      s = p;
      {
       int  ct;

       ct = 0;
       while  (* s != '\0')
         {
          if  (ct == 60)
              {
               putchar ('\n');
               ct = 0;
              }
          putchar (* s);
          ct ++;
          s ++;
         }
       if  (ct > 0)
           putchar ('\n');
      }
     }

   delete_ReadStruct(fsread);

   return  0;
  }


