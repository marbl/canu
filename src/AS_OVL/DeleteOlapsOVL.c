
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
* Module:  DeleteOlapsOVL.c
* Description:
*   Based on overlaps between DNA fragment sequences, make corrections
*   to single bases in the sequences.
* 
*    Programmer:  A. Delcher
*       Started:  19 Dec 2000
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: DeleteOlapsOVL.c,v 1.4 2005-03-22 19:49:18 jason_miller Exp $
 * $Revision: 1.4 $
*/

static char CM_ID[] = "$Id: DeleteOlapsOVL.c,v 1.4 2005-03-22 19:49:18 jason_miller Exp $";


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



//  Constants

#define  EXPANSION_FACTOR            1.4
    // Factor by which to grow memory in olap array when reading it



//  Type definitions

typedef  struct
  {
   int32  lo_iid, hi_iid;
  }  Frag_Pair_t;



//  Static Globals

static char  * Delete_Path;
    // Name of file with list of overlaps to delete
static Frag_Pair_t  * Delete_List = NULL;
    // Sorted list of overlaps to delete
static char  * New_OVL_Path;
    // Name of file to which to output OVLs after deletions
static int  Num_Deletes;
    // Number of entries in  Delete_List
static char  * Old_OVL_Path;
    // Name of file from which to read existing OVLs
static MesgWriter  Write_Msg_Fn;
    // Pointer to function to write OVL messages.



//  Static Functions

static void  Build_Delete_List
    (void);
static int  By_Lo_IID
    (const void * a, const void * b);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * in_stream, * out_stream;
   MesgReader  read_msg_fn;
   GenericMesg  * gmesg = NULL;
   OverlapMesg  * olm = NULL;

   Write_Msg_Fn = OutputFileType_AS (AS_BINARY_OUTPUT);

   Parse_Command_Line  (argc, argv);

   Build_Delete_List ();

   in_stream = File_Open (Old_OVL_Path, "r");
   read_msg_fn = InputFileType_AS (in_stream);

   out_stream = File_Open (New_OVL_Path, "w");

   while  (read_msg_fn (in_stream, & gmesg) != EOF)
     {
      if  (gmesg != NULL)
          {
           if  (gmesg -> t != MESG_OVL)
               Write_Msg_Fn (out_stream, gmesg);
             else
               {
                Frag_Pair_t  pair;
                
                olm = (OverlapMesg *) gmesg -> m;
                if  (olm -> aifrag < olm -> bifrag)
                    {
                     pair . lo_iid = olm -> aifrag;
                     pair . hi_iid = olm -> bifrag;
                    }
                  else
                    {
                     pair . lo_iid = olm -> bifrag;
                     pair . hi_iid = olm -> aifrag;
                    }

                if  (bsearch (& pair, Delete_List, Num_Deletes,
                              sizeof (Frag_Pair_t), By_Lo_IID) == NULL)
                    Write_Msg_Fn (out_stream, gmesg);
               }
          }
     }

   fclose (in_stream);
   fclose (out_stream);

   return  0;
  }



static void  Build_Delete_List
    (void)

//  Read list of overlaps from  Delete_Path  and save them
//  in  Delete_List .

  {
   FILE  * fp;
   int  list_size = 10000;
   int  lo, hi;

   Delete_List = (Frag_Pair_t *) Safe_malloc (list_size * sizeof (Frag_Pair_t));
   Num_Deletes = 0;

   fp = File_Open (Delete_Path, "r");

   while  (fscanf (fp, "%d %d", & lo, & hi) == 2)
     {
      assert (lo < hi);
      if  (Num_Deletes >= list_size)
          {
           list_size *= EXPANSION_FACTOR;
           Delete_List = (Frag_Pair_t *) Safe_realloc (Delete_List,
                          list_size * sizeof (Frag_Pair_t));
           assert (Num_Deletes < list_size);
          }

      Delete_List [Num_Deletes] . lo_iid = lo;
      Delete_List [Num_Deletes] . hi_iid = hi;
      Num_Deletes ++;
     }

   Delete_List = (Frag_Pair_t *) Safe_realloc (Delete_List,
                  Num_Deletes * sizeof (Frag_Pair_t));

   fclose (fp);

   return;
  }



static int  By_Lo_IID
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* Frag_Pair_t) 's,
//  first by  lo_iid , then by  hi_iid.
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Frag_Pair_t  * x, * y;

   x = (Frag_Pair_t *) a;
   y = (Frag_Pair_t *) b;

   if  (x -> lo_iid < y -> lo_iid)
       return  -1;
   else if  (x -> lo_iid > y -> lo_iid)
       return  1;
   else if  (x -> hi_iid < y -> hi_iid)
       return  -1;
   else if  (x -> hi_iid > y -> hi_iid)
       return  1;

   return  0;
  }



static int  OVL_Max_int
    (int a, int b)

//  Return the larger of  a  and  b .

  {
   if  (a < b)
       return  b;

   return  a;
  }



static int  OVL_Min_int
    (int a, int b)

//  Return the smaller of  a  and  b .

  {
   if  (a < b)
       return  a;

   return  b;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "-P")) != EOF))
     switch  (ch)
       {
        case  'P' :
          Write_Msg_Fn = OutputFileType_AS (AS_PROTO_OUTPUT);
          break;

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 3)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Old_OVL_Path = argv [optind ++];

   Delete_Path = argv [optind ++];

   New_OVL_Path = argv [optind ++];

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s [-P] <OldOVLFile> <DeleteListFile> <NewOVLFile>\n"
       "\n"
       "Copy contents of  <OldOVLFile>  to  <NewOVLFile>  but leave out\n"
       "any OVL messages that occur in  <DeleteListFile>\n"
       "\n"
       "Options:\n"
       "-P             Use ASCII (Proto) output\n",
       command);

   return;
  }

