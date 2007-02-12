
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
* Module:  FilterOlapsOVL.c
* Description:
*   Read an .ovl file and create a new one containing only those
*   overlaps that join unitigs in the input list of unitig pairs.
* 
*    Programmer:  A. Delcher
*       Started:  11 Jan 2000
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: FilterOlapsOVL.c,v 1.8 2007-02-12 22:16:57 brianwalenz Exp $
 * $Revision: 1.8 $
*/

static char CM_ID[] = "$Id: FilterOlapsOVL.c,v 1.8 2007-02-12 22:16:57 brianwalenz Exp $";


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
   int32  lo_id, hi_id;
  }  ID_Pair_t;



//  Static Globals

static int  IUM_List_Len = 0;
    // Number of entries in  IUM_List
static char  * IUM_List_Path;
    // Name of file with list unitig IDs of fragments
static int  * IUM_List = NULL;
    // Array of unitig IDs indexed by fragment IDs
static char  * New_OVL_Path;
    // Name of file to which to output OVLs after deletions
static int  Num_Unitig_Pairs = 0;
    // Number of entries in  Unitig_Pair_List
static char  * Old_OVL_Path;
    // Name of file from which to read existing OVLs
static ID_Pair_t  * Unitig_Pair_List = NULL;
    // Sorted list of overlaps to delete
static char  * Unitig_Pair_Path;
    // Name of file from which list of unitig pairs is read

//  Static Functions

static int  By_Lo_ID
    (const void * a, const void * b);
static int  OVL_Max_int
    (int a, int b);
static int  OVL_Min_int
    (int a, int b);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Read_IUM_List
    (void);
static void  Read_Unitig_Pairs
    (void);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * in_stream, * out_stream;
   GenericMesg  * gmesg = NULL;
   OverlapMesg  * olm = NULL;

   Parse_Command_Line  (argc, argv);

   Read_IUM_List ();

   Read_Unitig_Pairs ();

   in_stream = File_Open (Old_OVL_Path, "r");

   out_stream = File_Open (New_OVL_Path, "w");

   while  (ReadProtoMesg_AS (in_stream, & gmesg) != EOF)
     {
      if  (gmesg != NULL)
          {
           if  (gmesg -> t != MESG_OVL)
               WriteProtoMesg_AS (out_stream, gmesg);
             else
               {
                ID_Pair_t  pair;
                int  a_uni, b_uni;
                
                olm = (OverlapMesg *) gmesg -> m;
                a_uni = IUM_List [olm -> aifrag];
                b_uni = IUM_List [olm -> bifrag];

                if  (a_uni == b_uni)
                    WriteProtoMesg_AS (out_stream, gmesg);
                  else
                    {
                     void  * found;

                     if  (a_uni < b_uni)
                         {
                          pair . lo_id = a_uni;
                          pair . hi_id = b_uni;
                         }
                       else
                         {
                          pair . lo_id = b_uni;
                          pair . hi_id = a_uni;
                         }
                     found = bsearch (& pair, Unitig_Pair_List, Num_Unitig_Pairs,
                                   sizeof (ID_Pair_t), By_Lo_ID);
                     if  (found != NULL)
                         WriteProtoMesg_AS (out_stream, gmesg);
                    }
               }
          }
     }

   fclose (in_stream);
   fclose (out_stream);

   return  0;
  }



static void  Read_Unitig_Pairs
    (void)

//  Read list of unitig pairs from  Unitig_Pair_Path  and save them
//  in  Unitig_Pair_List .

  {
   FILE  * fp;
   int  list_size = 10000;
   int  lo, hi;

   Unitig_Pair_List = (ID_Pair_t *) safe_malloc (list_size * sizeof (ID_Pair_t));
   Num_Unitig_Pairs = 0;

   fp = File_Open (Unitig_Pair_Path, "r");

   while  (fscanf (fp, "%d %d", & lo, & hi) == 2)
     {
      assert (lo < hi);
      if  (Num_Unitig_Pairs >= list_size)
          {
           list_size *= EXPANSION_FACTOR;
           Unitig_Pair_List = (ID_Pair_t *) safe_realloc (Unitig_Pair_List,
                          list_size * sizeof (ID_Pair_t));
           assert (Num_Unitig_Pairs < list_size);
          }

      Unitig_Pair_List [Num_Unitig_Pairs] . lo_id = lo;
      Unitig_Pair_List [Num_Unitig_Pairs] . hi_id = hi;
      Num_Unitig_Pairs ++;
     }

   Unitig_Pair_List = (ID_Pair_t *) safe_realloc (Unitig_Pair_List,
                        Num_Unitig_Pairs * sizeof (ID_Pair_t));

   fclose (fp);

   return;
  }



static int  By_Lo_ID
    (const void * a, const void * b)

//  Compare the values in  a  and  b  as  (* ID_Pair_t) 's,
//  first by  lo_id , then by  hi_id.
//  Return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   ID_Pair_t  * x, * y;

   x = (ID_Pair_t *) a;
   y = (ID_Pair_t *) b;

   if  (x -> lo_id < y -> lo_id)
       return  -1;
   else if  (x -> lo_id > y -> lo_id)
       return  1;
   else if  (x -> hi_id < y -> hi_id)
       return  -1;
   else if  (x -> hi_id > y -> hi_id)
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
          fprintf(stderr, "-P depricated; protoIO is default.\n");
          break;

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 4)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Old_OVL_Path = argv [optind ++];

   IUM_List_Path = argv [optind ++];

   Unitig_Pair_Path = argv [optind ++];

   New_OVL_Path = argv [optind ++];

   return;
  }



static void  Read_IUM_List
    (void)

//  Read list of unitig IDs for each frament from file
//  IUM_List_Path  and save it in global  IUM_List .

  {
   FILE  * fp;

   fp = File_Open (IUM_List_Path, "rb");

   fread (& IUM_List_Len, sizeof (int), 1, fp);

   IUM_List = (int *) safe_malloc (IUM_List_Len * sizeof (int));

   fread (IUM_List, sizeof (int), IUM_List_Len, fp);

   fclose (fp);

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s [-P] <OldOVLFile> <IUMListFile> <UnitigPairFile> <NewOVLFile>\n"
       "\n"
       "Copy contents of  <OldOVLFile>  to  <NewOVLFile>  but leave out\n"
       "any OVL messages with frags not in same unitig or in  <UnitigPairFile>\n"
       "\n"
       "Options:\n"
       "-P             Use ASCII (Proto) output\n",
       command);

   return;
  }

