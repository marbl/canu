
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
* Module:  CatEratesOVL.c
* Description:
*   Concatenate erate correction files on command line to produce
*   a single file.
* 
*    Programmer:  A. Halpern 
*       Started:   23 Dec 2003
* 
* Assumptions:
*   Files on the command line *MUST* be in order
*   by fragment id of first frag in overlap.
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: CatEratesOVL.c,v 1.1 2004-09-23 20:32:57 mcschatz Exp $
 * $Revision: 1.1 $
*/

static char CM_ID[] = "$Id: CatEratesOVL.c,v 1.1 2004-09-23 20:32:57 mcschatz Exp $";


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
#include  "OlapStoreOVL.h"



//  Constants

#define  INITIAL_FILE_LIST_LEN   100
    //  Number of entries in initial list of correction files



//  Type definitions

typedef  struct
  {
   unsigned  is_ID : 1;
   unsigned  keep_left : 1;     // set true if left overlap degree is low
   unsigned  keep_right : 1;    // set true if right overlap degree is low
   unsigned  iid : 29;
  }  Frag_ID_t;

typedef  struct
  {
   unsigned  is_ID : 1;
   unsigned  pos : 20;    // position in fragment
   unsigned  type : 11;
  }  Correction_t;

typedef  union
  {
   Frag_ID_t  frag;
   Correction_t  corr;
  }  Correction_Output_t;

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
   int16 *erate;
   int32 header[3];
   int32 low=999999999;
   int32 high=-1;
   int32 ttlNum=0;

   Parse_Command_Line (argc, argv);

   outfile = File_Open (Outfile_Path, "wb");

   num_files = GetNumVA_char_Ptr_t (File_List);
   
   for  (i = 0;  i < num_files;  i ++)
     {
      char  * filename;
      int32 num,lo,hi;

      filename = * GetVA_char_Ptr_t (File_List, i);
      fp = File_Open (filename, "rb");
      Safe_fread (header, sizeof (int32), 3, fp);
      lo = header[0];
      hi = header[1];
      num = header [2];
      if(lo<low)low=lo;
      if(hi>high)high=hi;
      ttlNum+=num;
      fclose(fp);
      fprintf(stderr,"File %s lo %d hi %d num %d\n",
	      filename,lo,hi,num);
     }

   header[0]=low;
   header[1]=high;
   header[2]=ttlNum;

   fprintf(stderr,"total lo %d hi %d num %d\n",
	   low,high,ttlNum);

   Safe_fwrite(header,sizeof(int32),3,outfile);

   for  (i = 0;  i < num_files;  i ++)
     {
      char  * filename;
      int32 num,lo,hi;

      filename = * GetVA_char_Ptr_t (File_List, i);
      fp = File_Open (filename, "rb");
      Safe_fread (header, sizeof (int32), 3, fp);
      num = header[2];
      erate = (int16 *) Safe_malloc (num * sizeof (int16));
      Safe_fread (erate, sizeof (int16), num, fp);
      Safe_fwrite(erate, sizeof(int16),num,outfile);
      fclose (fp);
      free(erate);
     }

   fprintf (stderr, "Finished\n");

   return  0;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   char  buffer [MAX_FILENAME_LEN];
   int  ch, errflg = FALSE, len;
   int  i, n;
   char  * p;

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
        fprintf (stderr, "ERROR:  No erate files specified\n");
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
           "Concatenate erate corrections in <infile-1> ... <infile-n>\n"
           "to a single file <outfile>\n"
           "\n"
           "Options:\n"
           "  -L <listfile> Specify a file containing names of erate files\n"
           "     (useful if too many to fit on command line)\n"
           "  -o <outfile>  Specify output, REQUIRED\n",
           command);

   return;
  }



