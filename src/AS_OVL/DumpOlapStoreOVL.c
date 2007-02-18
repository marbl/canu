
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
* Module:  DumpOlapStoreOVL.c
* Description:
*   Dump (to stdout) the contents of an overlap store
* 
*    Programmer:  A. Delcher
*       Started:  19 Feb 2001
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: DumpOlapStoreOVL.c,v 1.6 2007-02-18 14:04:49 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char CM_ID[] = "$Id: DumpOlapStoreOVL.c,v 1.6 2007-02-18 14:04:49 brianwalenz Exp $";


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
#include  "AS_PER_gkpStore.h"
#include  "AS_PER_genericStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "OlapStoreOVL.h"


//  Constants

#define  DEFAULT_FRAGS_PER_FILE  20000       // Make bigger after testing
    //  Default value for global  Frags_Per_File
#define  VERBOSE                 0
    //  If  1  will print lots of extra output



//  Type definitions



//  Static Globals

static uint32  Frags_Per_File = DEFAULT_FRAGS_PER_FILE;
    // Number of fragments who overlaps are stored in each
    // (except maybe the last) data file in the store.
static char  * Input_Store_Path = NULL;
    // Path to store being added to
static uint32  Max_Old_Frag = 0;
    // The highest fragment ID in the previous store
static uint32  * Old_Offset = NULL;
    // Holds offsets to old overlaps in the store files
static char  * OVL_Store_Path = NULL;
    // Holds path to overlap store to be dumped
static uint32  Start_Frag = 1;
    // First fragment iid to dump
static uint32  Stop_Frag = UINT_MAX;
    // Last fragment iid to dump



//  Static Functions

static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
#if  0
   FILE  * fp;
   Short_Olap_Data_t  olap;
   uint32  header [3];
   int  file_index, num_data_files;
   char  filename [MAX_FILENAME_LEN];
   int  i, j;

   Input_Store_Path = argv [1];

   sprintf (filename, "%s/offset.olap", Input_Store_Path);
   fp = File_Open (filename, "rb");

   Safe_fread (header, sizeof (uint32), 3, fp);
   Max_Old_Frag = header [0];
   assert (header [1] == sizeof (Short_Olap_Data_t));
   Frags_Per_File = header [2];
   printf ("Hi Frag ID = %u\n", Max_Old_Frag);
   printf ("sizeof (Short_Olap_Data_t) = %u\n", header [1]);
   printf ("Frags per File = %u\n", Frags_Per_File);

   Old_Offset = (uint32 *) safe_malloc ((Max_Old_Frag + 2) * sizeof (uint32));
   Safe_fread (Old_Offset, sizeof (uint32), Max_Old_Frag + 2, fp);

   fclose (fp);

   file_index = 0;
   for  (i = 0;  i <= Max_Old_Frag;  i ++)
     {
      if  (i % Frags_Per_File == 1)
          {
           file_index ++;
           sprintf (filename, "%s/data%02d.olap", Input_Store_Path, file_index);
           fp = File_Open (filename, "rb");
          }

      printf ("> %8d  %9u:\n", i, Old_Offset [i]);
      for  (j = Old_Offset [i];  j < Old_Offset [i + 1];  j ++)
        {
         Safe_fread (& olap, sizeof (Short_Olap_Data_t), 1, fp);
         printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
                 i, olap . b_iid,
                 olap . flipped ? 'I' : 'N',
                 olap . a_hang, olap . b_hang,
                 olap . orig_erate / 10.0, olap . corr_erate / 10.0);
        }

      if  (i > 0 && i % Frags_Per_File == 0)
          fclose (fp);
     }
#else
   OVL_Store_t  * my_store;
   OVL_Stream_t  * my_stream;
   Long_Olap_Data_t  olap;
   uint32  last;

   Parse_Command_Line  (argc, argv);

   my_store = New_OVL_Store ();
   my_stream = New_OVL_Stream ();

   Open_OVL_Store (my_store, OVL_Store_Path);
   last = Last_Frag_In_OVL_Store (my_store);
   if  (Stop_Frag > last)
       Stop_Frag = last;
   Init_OVL_Stream (my_stream, Start_Frag, Stop_Frag, my_store);

   while  (Next_From_OVL_Stream (& olap, my_stream))
     {
      printf ("    %8d %8d %c %5d %5d %4.1f %4.1f\n",
              olap . a_iid, olap . b_iid,
              olap . flipped ? 'I' : 'N',
              olap . a_hang, olap . b_hang,
              Expand_Quality (olap . orig_erate) * 100.0,
              Expand_Quality (olap . corr_erate) * 100.0);
     }

   Free_OVL_Stream (my_stream);
   Free_OVL_Store (my_store);
#endif

   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "b:e:")) != EOF))
     switch  (ch)
       {
        case  'b' :
          Start_Frag = (uint32) strtol (optarg, & p, 10);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal start fragment \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  'e' :
          Stop_Frag = (uint32) strtol (optarg, & p, 10);
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal stop fragment \"%s\"\n",
                        optarg);
               errflg = TRUE;
              }
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);
          // fall through

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 1)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   OVL_Store_Path = argv [optind ++];

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [-b <Start_Frag>] [-e <End_Frag>]\n"
           "           <OVL_Store>\n"
           "\n"
           "Dumps overlap information in  <OVL_Store>  for\n"
           "fragments from  <Start_Frag>  to  <End_Frag> , inclusive.\n"
           "If  <Start_Frag>  is omitted will start at  1 .\n"
           "If  <End_Frag>  is omitted will end at last fragment .\n",
           command);

   return;
  }

