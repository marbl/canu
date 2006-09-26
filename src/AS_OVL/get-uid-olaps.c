
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
* Module:  get-uid-olaps.c
* Description:
*   Reads the stream of messages produced by the overlapper
*   and prints overlaps with uids
*
*    Programmer:  A. Delcher
*       Written:  3 Sep 99
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: get-uid-olaps.c,v 1.6 2006-09-26 21:07:45 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char fileID[] = "$Id: get-uid-olaps.c,v 1.6 2006-09-26 21:07:45 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include  "AS_OVL_delcher.h"


typedef  struct
  {
   int32 iid;
   int64 uid;
  }  ID_Pair_t;


static int  Cmp
    (const void * a, const void * b)
  {
   ID_Pair_t  * ap, * bp;

   ap = (ID_Pair_t *) a;
   bp = (ID_Pair_t *) b;

   if  (ap -> iid < bp -> iid)
       return  -1;
   else if  (bp -> iid < ap -> iid)
       return  1;
     else
       return  0;
  }



int main  (int argc, char * argv [])

  {
   FILE  * ovlfile;
   char  * infile_name;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   GenericMesg  * pmesg;
   ID_Pair_t  * list;
   int  list_len, frag_ct;
   int  ch, error_flag;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "")) != EOF))
     switch  (ch)
       {
        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind >= argc)
       {
        fprintf (stderr, 
                 "USAGE:  %s <ovl-file>\n", 
                 argv [0]);
        exit (EXIT_FAILURE);
       }

   infile_name = strdup (argv [optind ++]);
   assert (infile_name != NULL);
   fprintf (stderr, "Input Overlap File = %s\n", infile_name);

   ovlfile = File_Open (infile_name, "r");

   read_msg_fn = (MesgReader)InputFileType_AS (ovlfile);

   pmesg = (GenericMesg *) safe_malloc (sizeof (GenericMesg));
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) safe_malloc (sizeof (AuditMesg));
      
   list_len = 10000;
   list = (ID_Pair_t *) safe_calloc (list_len, sizeof (ID_Pair_t));
   frag_ct = 0;

   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_OFG :
          {
           OFGMesg  * ofg_mesg = gmesg -> m;
          
           if  (frag_ct >= list_len - 1)
               {list_len *= 2;
                list = (ID_Pair_t *) safe_realloc
                           (list, list_len * sizeof (ID_Pair_t));
               }

           list [frag_ct] . uid = ofg_mesg -> eaccession;
           list [frag_ct] . iid = ofg_mesg -> iaccession;
           frag_ct ++;

           break;
          }

        default :
          break;
       }

   rewind (ovlfile);

   qsort ((void *) list, (size_t) frag_ct, sizeof (ID_Pair_t), Cmp);

#if  0
   for  (i = 0;  i < frag_ct;  i ++)
     printf ("%7d %15ld\n", list [i] . iid, list [i] . uid);
   exit (-1);
#endif

   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_OVL :
          {
           OverlapMesg  * ovl_mesg = gmesg -> m;
           ID_Pair_t  * ap, * bp;

           ap = bsearch ((void *) (& (ovl_mesg -> aifrag)),
                         (void *) list, (size_t) frag_ct,
                         sizeof (ID_Pair_t), Cmp);
           bp = bsearch ((void *) (& (ovl_mesg -> bifrag)),
                         (void *) list, (size_t) frag_ct,
                         sizeof (ID_Pair_t), Cmp);
           
           if  (ap == NULL)
               {
                fprintf (stderr, "ERROR:  aid = %d not found\n",
                         ovl_mesg -> aifrag);
                exit (-1);
               }
           if  (bp == NULL)
               {
                fprintf (stderr, "ERROR:  bid = %d not found\n",
                         ovl_mesg -> bifrag);
                exit (-1);
               }
           printf ("%15" F_UIDP " %15" F_UIDP "  %4d  %4d  %c  %4.2f\n",
                    ap -> uid, bp -> uid,
                    ovl_mesg -> ahg,
                    ovl_mesg -> bhg,
                    (char) (ovl_mesg -> orientation),
                    100.0 * ovl_mesg -> quality);
           break;
          }

        default :
          break;
       }

   fclose (ovlfile);

   return  0;
  }



