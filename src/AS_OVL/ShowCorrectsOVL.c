
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
* Module:  ShowCorrectsOVL.c
* Description:
*   Based on overlaps between DNA fragment sequences, make corrections
*   to single bases in the sequences.
* 
*    Programmer:  A. Delcher
*       Started:   11 Dec 2000
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: ShowCorrectsOVL.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/

static char CM_ID[] = "$Id: ShowCorrectsOVL.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $";


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
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "FragCorrectOVL.h"



//  Type definitions



//  Static Functions

static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * fp;
   Correction_Output_t  msg;

   if  (argc < 2)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   fp = File_Open (argv [1], "rb");

   while  (fread (& msg, sizeof (Correction_Output_t), 1, fp) == 1)
     {
      if  (msg . frag . is_ID)
          printf (">%8u  %s  %s\n",
                  msg . frag . iid,
                  msg . frag . keep_left ? "KL" : "  ",
                  msg . frag . keep_right ? "KR" : "  ");
        else
          {
           printf (" %4u ", msg . corr . pos);
           switch  ((Vote_Value_t) msg . corr . type)
             {
              case  EXTENSION :
                printf ("extend\n");
                break;
              case  DELETE :
                printf ("d\n");
                break;
              case  A_SUBST :
                printf ("x a\n");
                break;
              case  C_SUBST :
                printf ("x c\n");
                break;
              case  G_SUBST :
                printf ("x g\n");
                break;
              case  T_SUBST :
                printf ("x t\n");
                break;
              case  A_INSERT :
                printf ("i a\n");
                break;
              case  C_INSERT :
                printf ("i c\n");
                break;
              case  G_INSERT :
                printf ("i g\n");
                break;
              case  T_INSERT :
                printf ("i t\n");
                break;
              default :
                fprintf (stderr, "ERROR:  Bad correction type = %d\n",
                         (int) msg . corr . type);
             }
          }
     }

   fclose (fp);

   return  0;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s <filename>\n"
           "\n"
           "Dump ASCII version of fragment corrections listed in\n"
           " <filename>  to stdout\n"
           "\n"
           "Options:\n",
           command);

   return;
  }



