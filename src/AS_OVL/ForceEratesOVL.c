
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

static char CM_ID[] = "$Id: ForceEratesOVL.c,v 1.3 2005-03-22 19:06:52 jason_miller Exp $";


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
#include  "OlapStoreOVL.h"



//  Constants


//  Type definitions



//  Static Globals

static uint16  Modified_Erate;
    // Value to set erate of specified overlaps to
static char  * OVL_Store_Path = NULL;
    // Location of overlap store to be updated
static char  * Olap_List_File_Path = NULL;
static FILE  * Olap_List_fp;
    // Name and fcb for file containing list of overlaps to set erate for
    // in the overlap store
static int32  Curr_A_ID;
    // First id of next batch of overlaps to change
static int32  * B_ID_List = NULL;
static int  B_ID_List_Len = 0;
static int  B_ID_List_Size = 0;
    // Holds B ids of next batch of overlaps t change
    // All have the same A id, which is  Curr_A_ID;
static int32  Next_A_ID, Next_B_ID;
    // Hold next line in  Olap_List_fp .


//  Static Functions

static int  Is_On_List
    (int32 x, int32 * list, int n);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Set_New_Erates
    (const char * const path);
static void  Update_ID_List
    (void);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   time_t  now;

   Parse_Command_Line  (argc, argv);

   Olap_List_fp = File_Open (Olap_List_File_Path, "r");

   {
    struct stat buffer ;
    fstat(fileno(Olap_List_fp),&buffer);
    if(buffer.st_size == 0){
      fprintf(stderr,"%s: Empty list of overlaps with erates to force (not necessarily a problem!)\n",argv[0]);
      exit(0);
    }
   }

   B_ID_List_Size = 100;
   B_ID_List = (int32 *) Safe_calloc (B_ID_List_Size, sizeof (int32));

   if  (fscanf (Olap_List_fp, "%d %d", & Next_A_ID, & Next_B_ID) != 2)
       {
        fprintf (stderr, "ERROR:  Can't read %s\n", Olap_List_File_Path);
        exit (EXIT_FAILURE);
       }
   Update_ID_List ();

   Set_New_Erates (OVL_Store_Path);

   now = time (NULL);
   fprintf (stderr, "Finished at %s\n", ctime (& now));
   
   return  0;
  }



static int  Is_On_List
    (int32 x, int32 * list, int n)

//  Return  TRUE  iff  x  is in  list [0 .. (n - 1)] .

  {
   int  i;

   for  (i = 0;  i < n;  i ++)
     if  (x == list [i])
         return  TRUE;

   return  FALSE;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   double  x;
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

   if  (errflg || optind != argc - 3)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   OVL_Store_Path = argv [optind ++];

   Olap_List_File_Path = argv [optind ++];

   x = strtod (argv [optind ++], NULL) / 100.0;
   Modified_Erate = (int16) Shrink_Quality ((float32) x);
   fprintf (stderr, "Setting error rate to %.3f%%\n",
        100.0 * Expand_Quality (Modified_Erate));
     // Convert and re-convert in case max is exceeded or there's another problem.

   return;
  }



static void  Set_New_Erates
    (const char * const path)

//  Set the  corr_erate  field in the overlap store in  path
//  to the value  Modified_Erate  for all overlaps that match
//  one in  Olap_List_fp .

  {
    FILE * fp = NULL;
    Short_Olap_Data_t  buff = {0};
    Short_Olap_Data_t  io_buff [IO_BUFF_SIZE] = {{0}};
    uint32  * olap_offset = NULL, max_frag = 0, frags_per_file = 0;
    int32  lo_id, hi_id;
    char  filename [MAX_FILENAME_LEN] = {0};
    size_t  file_position = 0, io_buff_start = 0;
    int  num_frags, io_buff_ct = 0, change_ct = 0;
    int  i, j, ct, first, file_index;
    
    assert(NULL != path);
    lo_id = 1;

    {
     uint32 header [3];
     sprintf (filename, "%s/offset.olap", path);
     assert(NULL == fp);
     fp = File_Open (filename, "rb");
     
     Safe_fread (header, sizeof (uint32), 3, fp);
     hi_id = max_frag = header [0];
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
     CDS_FSEEK (fp, lo_id * sizeof (uint32), SEEK_CUR);
     Safe_fread (olap_offset, sizeof (uint32), num_frags, fp);
     
     assert(NULL != fp);
     fclose (fp);
     fp = NULL; 
   }

   ct = 0;
   first = TRUE;
   file_index = (int) ceil ((double) lo_id / frags_per_file);

   for  (i = lo_id;  i <= hi_id;  i ++)
     {
      if  (Curr_A_ID != -1 && i > Curr_A_ID)
          {
           Update_ID_List ();
           assert  (Curr_A_ID == -1 || i <= Curr_A_ID);
          }

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
           CDS_FSEEK (fp, file_position, SEEK_SET);
           first = FALSE;
          }

      if  (i % frags_per_file == 0)
          {
           assert(NULL != fp);
           while  (fread (& buff, sizeof (Short_Olap_Data_t), 1, fp) == 1)
             {
              if  (i == Curr_A_ID
                     && Is_On_List (buff . b_iid, B_ID_List, B_ID_List_Len))
                  {
                   buff . corr_erate = Modified_Erate;
                   change_ct ++;
                  }
              io_buff [io_buff_ct ++] = buff;
              ct ++;
              file_position += sizeof (Short_Olap_Data_t);
              if  (io_buff_ct >= IO_BUFF_SIZE)
                  {
                   assert (io_buff_ct == IO_BUFF_SIZE);
                   CDS_FSEEK (fp, io_buff_start, SEEK_SET);
                   Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
                   io_buff_ct = 0;
                   io_buff_start = file_position;
                   CDS_FSEEK (fp, file_position, SEEK_SET);
                  }
             }
          }
        else
          {
           for  (j = olap_offset [i - lo_id];  j < olap_offset [i + 1 - lo_id];  j ++)
             {
              assert(NULL != fp);
              Safe_fread (& buff, sizeof (Short_Olap_Data_t), 1, fp);
              if  (i == Curr_A_ID
                     && Is_On_List (buff . b_iid, B_ID_List, B_ID_List_Len))
                  {
                   buff . corr_erate = Modified_Erate;
                   change_ct ++;
                  }
              io_buff [io_buff_ct ++] = buff;
              ct ++;
              file_position += sizeof (Short_Olap_Data_t);
              if  (io_buff_ct >= IO_BUFF_SIZE)
                  {
                   assert (io_buff_ct == IO_BUFF_SIZE);
                   CDS_FSEEK (fp, io_buff_start, SEEK_SET);
                   Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
                   io_buff_ct = 0;
                   io_buff_start = file_position;
                   CDS_FSEEK (fp, file_position, SEEK_SET);
                  }
             }
          }

      if  (i % frags_per_file == 0 || i == hi_id)
          {
           assert(NULL != fp);
           if  (io_buff_ct > 0)
               {
                CDS_FSEEK (fp, io_buff_start, SEEK_SET);
                Safe_fwrite (io_buff, sizeof (Short_Olap_Data_t), io_buff_ct, fp);
               }
           fclose (fp); 
           fp = NULL;
           file_index ++;
          }
     }

   fprintf (stderr, "%d changes made\n", change_ct);

   return;
  }



static void  Update_ID_List
    (void)

//  Put  Next_A_ID  into  Curr_A_ID  and read from
//   Olap_List_fp  putting B id's into  B_ID_List
//  until find a different  A id.  Put that A id into
//   Next_A_ID  and its corresponding B id into
//   Next_B_ID .  If hit end of file, set  Next_A_ID
//  to -1.

  {

   Curr_A_ID = Next_A_ID;
   B_ID_List [0] = Next_B_ID;
   B_ID_List_Len = 1;

   if  (Next_A_ID == -1)
       return;

   while  (fscanf (Olap_List_fp, "%d %d", & Next_A_ID, & Next_B_ID) == 2)
     {
      if  (Next_A_ID < Curr_A_ID)
          {
           fprintf (stderr, "ERROR:  IDs out of order  %d > %d\n",
                Curr_A_ID, Next_A_ID);
           exit (EXIT_FAILURE);
          }
      if  (Curr_A_ID < Next_A_ID)
          return;
      // Same A id.  Add B id to list
      if  (B_ID_List_Len == B_ID_List_Size)
          {
           B_ID_List_Size *= 2;
           B_ID_List = (int32 *) Safe_realloc (B_ID_List,
                B_ID_List_Size * sizeof (int32));
          }
      B_ID_List [B_ID_List_Len ++] = Next_B_ID;
     }

   Next_A_ID = Next_B_ID = -1;

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s <OVLStore> <OlapListFile> <NewErate>\n"
       "\n"
       "Sets corrected error rates in <OvlStore> to <NewErate>\n"
       "for overlaps in <OlapListFile>\n"
       "<OlapListFile> must be sorted by id of first entry\n"
       "<NewErate> is expressed as a percent, e.g., 9.9 for 9.9%% error\n"
       "\n"
       "Options:\n"
       "-h  print this message\n",
       command);

   return;
  }
