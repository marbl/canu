
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
* Module:  UpdateEratesOVL.c
* Description:
*   Add binary list of overlap error rates to corrected error-rate
*   field of appropriate entries in an overlap store.
* 
*    Programmer:  A. Delcher
*       Started:  16 Mar 2001
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: UpdateEratesOVL.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/

static char CM_ID[] = "$Id: UpdateEratesOVL.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $";


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
#include  "OlapStoreOVL.h"


//  Type definitions



//  Static Globals

static char  * OVL_Store_Path = NULL;
    // Location of overlap store to be updated
static char  * Erate_File_Path = NULL;
    // File containing corrected error rates to be stored
    // in the overlap store



//  Static Functions

static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Set_Corrected_Erate
    (const char * const path, int32 lo_id, int32 hi_id, uint16 * erate, int num);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   int32  header [3];
   uint16  * erate = NULL;
   int  num;

   Parse_Command_Line  (argc, argv);

   {
     FILE * fp = File_Open (Erate_File_Path, "rb");
     Safe_fread (header, sizeof (int32), 3, fp);
     num = header [2];
     erate = (uint16 *) Safe_malloc (num * sizeof (uint16));
     Safe_fread (erate, sizeof (uint16), num, fp);
     fclose (fp);
   }
   assert(NULL != erate);
   Set_Corrected_Erate (OVL_Store_Path, header [0], header [1], erate, num);

   free(erate);
   fprintf (stderr, "Finished\n");
   
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

   OVL_Store_Path = argv [optind ++];

   Erate_File_Path = argv [optind ++];

   return;
  }



static void  Set_Corrected_Erate
    (const char * const path, int32 lo_id, int32 hi_id, uint16 * erate, int num)

//  Set the  corr_erate  field in the overlap store in  path .
//  Set it for fragments  lo_id .. hi_id  to the values in
//   erate [0 .. (num - 1)]  which are in the same order as the
//  overlaps in the store.

  {
    FILE * fp = NULL;
    Short_Olap_Data_t  buff = {0};
    Short_Olap_Data_t  io_buff [IO_BUFF_SIZE] = {{0}};
    uint32  * olap_offset = NULL, max_frag = 0, frags_per_file = 0;
    char  filename [MAX_FILENAME_LEN] = "\0";
    size_t  file_position = 0, io_buff_start = 0;
    int  num_frags, io_buff_ct = 0;
    int  i, j, ct, first, file_index;
    
    assert(NULL != path);
    {
     uint32 header [3];
     sprintf (filename, "%s/offset.olap", path);
     assert(NULL == fp);
     fp = File_Open (filename, "rb");
     
     Safe_fread (header, sizeof (uint32), 3, fp);
     max_frag = header [0];
     assert (header [1] == sizeof (Short_Olap_Data_t));
     frags_per_file = header [2];
     fprintf (stderr, "sizeof (Short_Olap_Data_t) = %u\n", header [1]);
     fprintf (stderr, "Frags per File = %u\n", frags_per_file);
     
     assert (1 <= lo_id && lo_id <= hi_id);
     if  (max_frag < lo_id)
         {
          fprintf (stderr,
                   "ERROR:  No overlaps for frags in this range--nothing to do\n");
          return;
         }
     if  (max_frag < hi_id)
         {
          fprintf (stderr, "Hi frag %d past last ovlStore frag %d\n",
                   hi_id, max_frag);
          hi_id = max_frag;
         }

     num_frags = 2 + hi_id - lo_id;   // go 1 past the end
     olap_offset = (uint32 *) Safe_malloc (num_frags * sizeof (uint32));
     CDS_FSEEK (fp, (off_t) (lo_id * sizeof (uint32)), SEEK_CUR);
     Safe_fread (olap_offset, sizeof (uint32), num_frags, fp);
     
     assert(NULL != fp);
     fclose (fp);
     fp = NULL; 
   }

   ct = 0;
   first = TRUE;
   file_index = (int) ceil ((double) lo_id / frags_per_file);

   for  (i = lo_id;  i <= hi_id;  i ++) {
     if  (first || i % frags_per_file == 1)
       {
	   assert(file_index < 100);
           sprintf (filename, "%s/data%02d.olap", path, file_index);
	   assert(NULL == fp);
           fp = File_Open (filename, "rb+");
	   assert(NULL != fp);

           file_position = olap_offset [i - lo_id] * sizeof (Short_Olap_Data_t);
           io_buff_start = file_position;
           io_buff_ct = 0;
           CDS_FSEEK(fp, (off_t) file_position, SEEK_SET);
           first = FALSE;
          }

     if  (i % frags_per_file == 0){
	assert(NULL != fp);
	while  (fread (& buff, sizeof (Short_Olap_Data_t), 1, fp) == 1)
            {
             buff . corr_erate = erate [ct];
             io_buff [io_buff_ct ++] = buff;
             ct ++;
             file_position += sizeof (Short_Olap_Data_t);
             if  (io_buff_ct >= IO_BUFF_SIZE)
                 {
                  assert (io_buff_ct == IO_BUFF_SIZE);
                  CDS_FSEEK(fp, (off_t) io_buff_start, SEEK_SET);
                  Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
                  io_buff_ct = 0;
                  io_buff_start = file_position;
                  CDS_FSEEK (fp, (off_t) file_position, SEEK_SET);
                 }
            }
     } else {
	for  (j = olap_offset [i - lo_id];  j < olap_offset [i + 1 - lo_id];  j ++)
	  {
	     assert(NULL != fp);
             Safe_fread (& buff, sizeof (Short_Olap_Data_t), 1, fp);
             buff . corr_erate = erate [ct];
             io_buff [io_buff_ct ++] = buff;
             ct ++;
             file_position += sizeof (Short_Olap_Data_t);
             if  (io_buff_ct >= IO_BUFF_SIZE)
                 {
                  assert (io_buff_ct == IO_BUFF_SIZE);
                  CDS_FSEEK (fp, (off_t) io_buff_start, SEEK_SET);
                  Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
                  io_buff_ct = 0;
                  io_buff_start = file_position;
                  CDS_FSEEK (fp, (off_t) file_position, SEEK_SET);
                 }
            }
      }

      if  (i % frags_per_file == 0 || i == hi_id)
          {
	   assert(NULL != fp);
           if  (io_buff_ct > 0)
               {
                CDS_FSEEK (fp, (off_t) io_buff_start, SEEK_SET);
                Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
               }
           fclose (fp); 
	   fp = NULL;
           file_index ++;
          }
   }

   if  (ct != num)
       {
        fprintf (stderr, "Number of corrected erates = %d\n", num);
        fprintf (stderr, "OVLStore entries for frags %d .. %d = %d\n",
                 lo_id, hi_id, ct);
        fprintf (stderr, "ERROR:  They *must* match\n");
        exit (EXIT_FAILURE);
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
       "USAGE:  %s <OVLStore> <ErateFile>\n"
       "\n"
       "Adds new error rates in <ErateFile> to overlap store <OVLStore>\n"
       "as corrected error rates\n"
       "\n"
       "Options:\n"
       "-h  print this message\n",
       command);

   return;
  }
