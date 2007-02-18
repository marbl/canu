
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
* Module:  CatCorrectsOVL.c
* Description:
*   Concatenate fragment correction files on command line to produce
*   a single file.
* 
*    Programmer:  A. Delcher
*       Started:   24 Jan 2001
* 
* Assumptions:
*   Files on the command line *MUST* be in order
*   by fragment id.
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: CatCorrectsOVL.c,v 1.6 2007-02-18 14:04:49 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char CM_ID[] = "$Id: CatCorrectsOVL.c,v 1.6 2007-02-18 14:04:49 brianwalenz Exp $";


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
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "OlapStoreOVL.h"
#include  "FragCorrectOVL.h"


//  Constants

#define  INITIAL_FILE_LIST_LEN   100
    //  Number of entries in initial list of correction files




typedef  char *  char_Ptr_t;

VA_DEF (char_Ptr_t)



//  Global Variables

static VA_TYPE (char_Ptr_t)  * File_List;
    // Array of names of correction files to process
static char  * File_List_Path = NULL;
    // Path to file with names of correction files to process
static char  * Outfile_Path = NULL;



//  Static Functions

static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   FILE  * fp, * outfile;
   Correction_Output_t  msg;
   int32  prev_id = -1;
   int  i, num_files;

   Parse_Command_Line (argc, argv);

   outfile = File_Open (Outfile_Path, "wb");

   num_files = GetNumVA_char_Ptr_t (File_List);
   for  (i = 0;  i < num_files;  i ++)
     {
      char  * filename;

      filename = * GetVA_char_Ptr_t (File_List, i);
      fp = File_Open (filename, "rb");

      while  (fread (& msg, sizeof (Correction_Output_t), 1, fp) == 1)
        {
         if  (msg . frag . is_ID)
             {
              if  ((int) msg . frag . iid <= prev_id)
                  {
                   fprintf (stderr,
                            "ERROR:  frag IDs out of order\n"
                            "Hit frag  %d  in file  %s  preceeding frag was %d\n",
                            msg . frag . iid, argv [i], prev_id);
                   exit (EXIT_FAILURE);
                  }
              prev_id = (int) msg . frag . iid;
             }

         fwrite (& msg, sizeof (Correction_Output_t), 1, outfile);
        }

      fclose (fp);
     }

   fclose (outfile);

   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   char  buffer [MAX_FILENAME_LEN];
   int  ch, errflg = FALSE;
   int  i, n;

   File_List = CreateVA_char_Ptr_t (INITIAL_FILE_LIST_LEN);
   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "L:o:")) != EOF))
     switch  (ch)
       {
        case  'L' :
          File_List_Path = optarg;
          break;

        case  'o' :
          Outfile_Path = optarg;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   if  (Outfile_Path == NULL)
       {
        fprintf (stderr, "No output file specified with -o\n");
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   if  (File_List_Path != NULL)
       {
        FILE  * fp;
        char_Ptr_t  name;

        fp = File_Open (File_List_Path, "r");
        while  (fscanf (fp, "%s", buffer) == 1)
          {
           name = strdup (buffer);
           Appendchar_Ptr_t (File_List, & name);
          }

        fclose (fp);
       }

   while  (optind < argc)
     {
      Appendchar_Ptr_t (File_List, argv + optind);
      optind ++;
     }

   n = GetNumVA_char_Ptr_t (File_List);
   if  (n <= 0)
       {
        fprintf (stderr, "ERROR:  No correction files specified\n");
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }
   for  (i = 0;  i < n;  i ++)
     {
      FILE  * fp;

      fp = File_Open (* GetVA_char_Ptr_t (File_List, i), "r");
      fclose (fp);
     }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s -o <outfile> [-L <listfile>]  <infile-1> ... <infile-n>\n"
           "\n"
           "Concatenate fragment corrections in <infile-1> ... <infile-n>\n"
           "to a single file <outfile>\n"
           "\n"
           "Options:\n"
           "  -L <listfile> Specify a file containing names of correction files\n"
           "     (useful if too many to fit on command line)\n"
           "  -o <outfile>  Specify output, REQUIRED\n",
           command);

   return;
  }



