
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
* Module:  GrowOlapStoreOVL.c
* Description:
*   Create or add to a compact, binary-file, sorted list of
*   fragment overlaps.
* 
*    Programmer:  A. Delcher
*       Started:  8 Feb 2001
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: GrowOlapStoreOVL.c,v 1.2 2004-09-23 20:25:25 mcschatz Exp $
 * $Revision: 1.2 $
*/

static char CM_ID[] = "$Id: GrowOlapStoreOVL.c,v 1.2 2004-09-23 20:25:25 mcschatz Exp $";


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
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_UTL_version.h"
#include  "OlapStoreOVL.h"



//  Constants

#define  DEFAULT_FRAGS_PER_FILE  2000000
    //  Default value for global  Frags_Per_File
#define  INITIAL_FILE_LIST_LEN   100
    //  Number of entries in initial list of .ovl files
#define  MAX_BACKUPS             10
    //  Most backup versions of an overlap store that are allowed
    //  Uses about 3GB of memory
    //  The most overlaps that can be stored in memory at a time



//  Type definitions

typedef  struct
  {
   signed int  needs_update : 2;
   signed int  old_exists : 2;
  }  File_Status_t;

typedef  struct
  {
   FILE  * fp;
   int  ct;
  }  Tmp_File_Info_t;

typedef  char *  char_Ptr_t;

VA_DEF (char_Ptr_t)



//  Static Globals

static int  Delete_Old = FALSE;
    // Set true if appending without backup.  Indicates
    // old files should be deleted after the new are safely
    // written.
static char  * Dump_Format_Input_Path = NULL;
    // Name of file holding overlaps in format produced by  dump-olap-store
static File_Status_t  * File_Status = NULL;
    // Holds information about files that store overlap data
static int  File_Status_Size = 0;
    // Number of entries allocated for  File_Status;
static uint32  Frags_Per_File = DEFAULT_FRAGS_PER_FILE;
    // Number of fragments whose overlaps are stored in each
    // (except maybe the last) data file in the store.
static uint32  Highest_Frag;
    // Highest frag ID either in old store or new overlaps
static char  * Input_Store_Path = NULL;
    // Path to store being added to
static uint32  Max_ID = 0;
    // The highest fragment ID seen so far
static uint32  Max_Old_Frag = 0;
    // The highest fragment ID in the previous store
static uint32  * New_Offset = NULL;
    // Holds offsets to  New_Olap  array
static int32  New_Offset_Size = 0;
    // Number of entries in  New_Offset  array
static Long_Olap_Data_t  New_Olap [MAX_OLAP_BATCH];
    // Holds current batch of overlaps being processed
static int  Num_Files;
    // The number of data files holding overlaps
static int32  Num_Olaps;
    // Number of entries in  New_Olap  array
static uint32  * Old_Offset = NULL;
    // Holds offsets to old overlaps in the store files
static char  * Output_Store_Path = NULL;
    // Path to resulting store
static VA_TYPE (char_Ptr_t)  * OVL_File_List;
    // Array of names of .ovl files to process
static char  * OVL_File_List_Path = NULL;
    // Path to file with names of .ovl files to process
static int  Save_Space = FALSE;
    // If  TRUE  will delete and rename .tmp files as they
    // are produced rather than waiting to end of entire batch
static Short_Olap_Data_t  Sorted_Olap [MAX_OLAP_BATCH];
    // Holds sorted list of current batch of overlaps
static Tmp_File_Info_t  * Tmp_File_Info = NULL;
    // Holds information about temporary files of overlaps
static int  Tmp_File_Info_Size = 0;
    // Number of entries allocated for  Tmp_File_Info;
static int  Verbose = 0;
    // Determines level of extra printouts



//  Static Functions

static void  Output_New_Batch
    (void);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Process_New_Batch
    (Long_Olap_Data_t new_olap [], int num_olaps);
static void  Process_Olaps_From_Dump_File
    (void);
static void  Process_OVL_Files
    (void);
static void  Process_Tmp_File
    (int i);
static void  Tmp_Output_Olap
    (const Long_Olap_Data_t * olap);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   time_t  now;
   int  i;

   now = time (NULL);
   fprintf (stderr, "### Starting at  %s\n", ctime (& now));

   Parse_Command_Line  (argc, argv);

   File_Status_Size = INITIAL_FILE_LIST_LEN;
   File_Status
       = (File_Status_t *) Safe_malloc (File_Status_Size * sizeof (File_Status_t));

   Tmp_File_Info_Size = INITIAL_FILE_LIST_LEN;
   Tmp_File_Info
       = (Tmp_File_Info_t *) Safe_malloc
             (Tmp_File_Info_Size * sizeof (Tmp_File_Info_t));

   Num_Files = 0;
   Num_Olaps = 0;

   if  (Dump_Format_Input_Path != NULL)
       Process_Olaps_From_Dump_File ();
     else
       Process_OVL_Files ();


   // Now add overlaps from any non-empty temp files to the store

   for  (i = 1;  i <= Num_Files;  i ++)
     {
      if  (Tmp_File_Info [i] . ct > 0)
          Process_Tmp_File (i);
      if  (Tmp_File_Info [i] . fp != NULL)
          {
           char  file_name [MAX_FILENAME_LEN];

           sprintf (file_name, "%s/new%02d.tmp", Output_Store_Path, i);
           Safe_remove (file_name);
          }
     }

   now = time (NULL);
   fprintf (stderr, "### Finished at  %s\n", ctime (& now));

   return  0;
  }



static void  Output_New_Batch
    (void)

//  Merge overlaps in  Sorted_Olap  with previous overlaps
//  (if any) and output them to the store in  Output_Store_Path .

  {
   FILE  * new_fp, * old_fp = NULL;
   Short_Olap_Data_t  olap;
   char  input_name [MAX_FILENAME_LEN], output_name [MAX_FILENAME_LEN];
   char  offset_name [MAX_FILENAME_LEN];
   int  file_index = 0;
   uint32  header [3];
   uint32  * combined_offset;
   int  i, j, where = 0;

   fprintf (stderr, "Outputting new overlaps  Highest_Frag = %d\n", Highest_Frag);

   combined_offset = (uint32 *) Safe_malloc ((Highest_Frag + 2) * sizeof (uint32));
   combined_offset [0] = 0;

   for  (file_index = 1;  file_index <= Num_Files;  file_index ++)
     {
      uint32  lo_frag, hi_frag;

      where = 0;
      hi_frag = file_index * Frags_Per_File;
      lo_frag = 1 + hi_frag - Frags_Per_File;
      if  (hi_frag > Highest_Frag)
          hi_frag = Highest_Frag;

      sprintf (output_name, "%s/data%02d.tmp", Output_Store_Path, file_index);
      if  (Verbose > 0)
          {
           fprintf (stderr, "file_index = %d  Save_Space = %c  old_exists = %c"
                    "  needs_update = %c\n", file_index, Save_Space ? 'T' : 'F',
                    File_Status [file_index] . old_exists ? 'T' : 'F',
                    File_Status [file_index] . needs_update ? 'T' : 'F');
           fprintf (stderr, "lo_frag = %d  hi_frag = %d  Max_Old_Frag = %d\n",
                    lo_frag, hi_frag, Max_Old_Frag);
           if  (File_Status [file_index] . old_exists)
               {
                int  i;

                fprintf (stderr, "Old_Offset [%d] = %d   Old_Offset [%d] = %d\n",
                         Max_Old_Frag, Old_Offset [Max_Old_Frag], Max_Old_Frag + 1,
                         Old_Offset [Max_Old_Frag + 1]);

                for  (i = 0;  i < 10;  i ++)
                  fprintf (stderr, "Old_Offset [%2d] = %7d\n", i, Old_Offset [i]);
               }
          }

      if  (Save_Space && ! File_Status [file_index] . needs_update)
          {
           if  (File_Status [file_index] . old_exists)
               {
                uint32  prev = 0;

                for  (j = lo_frag;  j <= hi_frag;  j ++)
                  if  (j <= Max_Old_Frag + 1)
                      prev = combined_offset [j] = Old_Offset [j];
                    else
                      combined_offset [j] = prev;
                if  (hi_frag == Max_Old_Frag)
                    where = Old_Offset [hi_frag + 1];
               }
             else
               {
                for  (j = lo_frag;  j <= hi_frag;  j ++)
                  combined_offset [j] = 0;
                new_fp = File_Open (output_name, "wb");
                fclose (new_fp);
                File_Status [file_index] . needs_update = TRUE;
                where = 0;
               }
           continue;
          }

      new_fp = File_Open (output_name, "wb");

      if  (lo_frag <= Max_Old_Frag)
          {
           assert (Input_Store_Path != NULL);
           sprintf (input_name, "%s/data%02d.olap", Input_Store_Path, file_index);
           old_fp = File_Open (input_name, "rb");
          }

      for  (i = lo_frag;  i <= hi_frag;  i ++)
        {
         combined_offset [i] = where;

         if  (i <= Max_Old_Frag)
             {
              // If this is the last fragment in the data file the next offset
              // will be zero in the next file, and can't be used to determine
              // how many overlaps to read
              if  (i % Frags_Per_File == 0)
                  while  (fread (& olap, sizeof (Short_Olap_Data_t), 1, old_fp) == 1)
                    {
                     Safe_fwrite (& olap, sizeof (Short_Olap_Data_t), 1, new_fp);
                     where ++;
                    }
                else
                  for  (j = Old_Offset [i];  j < Old_Offset [i + 1];  j ++)
                    {
                     Safe_fread (& olap, sizeof (Short_Olap_Data_t), 1, old_fp);
                     Safe_fwrite (& olap, sizeof (Short_Olap_Data_t), 1, new_fp);
                     where ++;
                    }
             }
         if  (i <= Max_ID)
             for  (j = New_Offset [i];  j < New_Offset [i + 1];  j ++)
               {
                Safe_fwrite (Sorted_Olap + j, sizeof (Short_Olap_Data_t), 1, new_fp);
                where ++;
               }
        }

      fclose (new_fp);
      if  (lo_frag <= Max_Old_Frag)
          {
           fclose (old_fp);
           if  (Save_Space)
               Safe_remove (input_name);
          }

      File_Status [file_index] . old_exists = TRUE;
      if  (Verbose > 0)
          {
           int  i;

           for  (i = 0;  i < 10;  i ++)
             fprintf (stderr, "combined_offset [%2d] = %7d\n", i, combined_offset [i]);
          }
     }

   if  (Highest_Frag % Frags_Per_File == 0)
       where = 0;
   combined_offset [Highest_Frag + 1] = where;

   header [0] = Highest_Frag;
   header [1] = sizeof (Short_Olap_Data_t);
   header [2] = Frags_Per_File;
   sprintf (offset_name, "%s/offset.tmp", Output_Store_Path);
   new_fp = File_Open (offset_name, "wb");
   Safe_fwrite (header, sizeof (uint32), 3, new_fp);
   Safe_fwrite (combined_offset, sizeof (uint32), Highest_Frag + 2, new_fp);
   if  (Verbose > 0)
       {
        int  i;

        fprintf (stderr, "where = %d\n", where);
        fprintf (stderr, "                      Highest_Frag = %9d\n", Highest_Frag);
        for  (i = Highest_Frag - 9;  i <= Highest_Frag + 1;  i ++)
          if  (i >= 0)
              fprintf (stderr, "    combined_offset [%7d] = %9d\n",
                       i, combined_offset [i]);
       }
   fclose (new_fp);

   if  (Delete_Old)
       {
        int  ct = (int) ceil ((double) Max_Old_Frag / Frags_Per_File);

        for  (i = 1;  i <= ct && ! Save_Space;  i ++)
          {
           sprintf (input_name, "%s/data%02d.olap", Input_Store_Path, i);
           Safe_remove (input_name);
          }
        sprintf (input_name, "%s/offset.olap", Input_Store_Path);
        Safe_remove (input_name);
       }

   fprintf (stderr, "Renaming tmp files & writing offset file\n");
   for  (i = 1;  i <= Num_Files;  i ++)
     if  (! Save_Space || File_Status [i] . needs_update)
         {
          sprintf (input_name, "%s/data%02d.tmp", Output_Store_Path, i);
          sprintf (output_name, "%s/data%02d.olap", Output_Store_Path, i);
          Safe_rename (input_name, output_name);
         }
   sprintf (input_name, "%s/offset.tmp", Output_Store_Path);
   sprintf (output_name, "%s/offset.olap", Output_Store_Path);
   Safe_rename (input_name, output_name);


   // Now set up for next batch in case there is one

   free (Old_Offset);
   Old_Offset = combined_offset;
   Max_Old_Frag = Highest_Frag;

   if  (Output_Store_Path != NULL)
       {
        Input_Store_Path = Output_Store_Path;
        Delete_Old = TRUE;
       }

   return;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   char  buffer [MAX_FILENAME_LEN], command [MAX_FILENAME_LEN];
   int  ch, errflg = FALSE;
   int  append = FALSE;
   int  backup = FALSE;
   int  create = FALSE;
   int  force = FALSE;
   int  i, n;

   fprintf (stderr, "* Working directory is %s\n", getcwd (NULL, 256));
   OVL_File_List = CreateVA_char_Ptr_t (INITIAL_FILE_LIST_LEN);

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "aAcD:fi:L:o:hSv:")) != EOF))
     switch  (ch)
       {
        case  'a' :
          append = TRUE;
          backup = TRUE;
          break;

        case  'A' :
          append = TRUE;
          backup = FALSE;
          break;

        case  'c' :
          create = TRUE;
          break;

        case  'D' :
          Dump_Format_Input_Path = optarg;
          break;

        case  'f' :
          force = TRUE;
          break;

        case  'i' :
          Input_Store_Path = optarg;
          break;

        case  'L' :
          OVL_File_List_Path = optarg;
          break;

        case  'o' :
          Output_Store_Path = optarg;
          break;

        case  'S' :
          Save_Space = TRUE;
          break;

        case  'v' :
          Verbose = strtol (optarg, NULL, 10);
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

   if  (OVL_File_List_Path != NULL)
       {
        FILE  * fp;
        char_Ptr_t  name;

        fp = File_Open (OVL_File_List_Path, "r");
        while  (fscanf (fp, "%s", buffer) == 1)
          {
           name = strdup (buffer);
           Appendchar_Ptr_t (OVL_File_List, & name);
          }

        fclose (fp);
       }

   while  (optind < argc)
     {
      Appendchar_Ptr_t (OVL_File_List, argv + optind);
      optind ++;
     }

   if  (backup && ! append)
       {
        fprintf (stderr, "ERROR:  Backup only with append\n");
        exit (EXIT_FAILURE);
       }
   if  (force && append)
       {
        fprintf (stderr, "ERROR:  Cannot force and append\n");
        exit (EXIT_FAILURE);
       }
   if  (create && append)
       {
        fprintf (stderr, "ERROR:  Cannot create and append\n");
        exit (EXIT_FAILURE);
       }
   if  (append && (Input_Store_Path == NULL || Output_Store_Path != NULL))
       {
        fprintf (stderr, "ERROR:  On append specify input store only (-i)\n");
        exit (EXIT_FAILURE);
       }
   if  (! append && (Output_Store_Path == NULL
          || (! create && Input_Store_Path == NULL)))
       {
        fprintf (stderr,
                 "ERROR:  Must specify output store (-o) and\n"
                 "        input store (-i) if not create (-c)\n");
        exit (EXIT_FAILURE);
       }

   Delete_Old = (append ||
                   (Input_Store_Path != NULL && Output_Store_Path != NULL
                      && strcmp (Input_Store_Path, Output_Store_Path) == 0));

   n = GetNumVA_char_Ptr_t (OVL_File_List);
   if  (n <= 0 && Dump_Format_Input_Path == NULL)
       {
        fprintf (stderr, "ERROR:  No .ovl files specified\n");
        exit (EXIT_FAILURE);
       }
   for  (i = 0;  i < n;  i ++)
     {
      FILE  * fp;

      fp = File_Open (* GetVA_char_Ptr_t (OVL_File_List, i), "r");
      fclose (fp);
     }

   if  (backup)
       {
        for  (i = 1;  i <= MAX_BACKUPS;  i ++)
          {
           sprintf (buffer, "%s.back%03d", Input_Store_Path, i);
           if  (mkdir (buffer, S_IRWXU | S_IRWXG | S_IROTH) == 0)
               break;
           if  (errno != EEXIST)
               {
                fprintf (stderr, "ERROR:  Can't create directory \"%s\"\n",
                         buffer);
                exit (EXIT_FAILURE);
               }
          }
        if  (i > MAX_BACKUPS)
            {
             fprintf (stderr, "ERROR:  Can't create backup directory\n");
             exit (EXIT_FAILURE);
            }
        sprintf (command, "cp %s/offset.olap %s/data*.olap %s",
                 Input_Store_Path, Input_Store_Path, buffer);
        if  (system (command) != 0)
            {
             fprintf (stderr, "ERROR:  Failed to backup store\n");
             exit (EXIT_FAILURE);
            }
       }

   if  (Output_Store_Path == NULL)
       Output_Store_Path = Input_Store_Path;
     else
       {
        if  (mkdir (Output_Store_Path, S_IRWXU | S_IRWXG | S_IROTH) != 0)
            {
             if  (errno == EEXIST)
                 {
                  if  (force)
                      {
                       sprintf (command, "rm -f %s/*.olap", Output_Store_Path);
                       if  (system (command) != 0)
                           {
                            fprintf (stderr, "ERROR:  Failed to clear existing store\n");
                            exit (EXIT_FAILURE);
                           }
                      }
                    else
                      {
                       fprintf (stderr,
                       "ERROR:  Directory \"%s\" exists.  Use -f to overwrite it.\n",
                                Output_Store_Path);
                       exit (EXIT_FAILURE);
                      }
                 }
               else
                 {
                  fprintf (stderr,
                  "ERROR:  Failed to create directory \"%s\".\n",
                           Output_Store_Path);
                  exit (EXIT_FAILURE);
                 }
            }
       }

   if  (Input_Store_Path != NULL)
       {
        FILE  * fp;
        uint32  header [3];
        int  num_data_files;

        sprintf (buffer, "%s/offset.olap", Input_Store_Path);
        fp = File_Open (buffer, "rb");

        Safe_fread (header, sizeof (uint32), 3, fp);
        Max_Old_Frag = header [0];
        assert (header [1] == sizeof (Short_Olap_Data_t));
        Frags_Per_File = header [2];

        Old_Offset = (uint32 *) Safe_malloc ((Max_Old_Frag + 2) * sizeof (uint32));
        Safe_fread (Old_Offset, sizeof (uint32), Max_Old_Frag + 2, fp);

        fclose (fp);

        if  (Verbose > 0)
            {
             fprintf (stderr, "                 Max_Old_Frag = %9d\n",
                      Max_Old_Frag);
             fprintf (stderr, "    Old_Offset [Max_Old_Frag] = %9d\n",
                      Old_Offset [Max_Old_Frag]);
             fprintf (stderr, "Old_Offset [Max_Old_Frag + 1] = %9d\n",
                      Old_Offset [Max_Old_Frag + 1]);
            }

        // Check that can open all data files.
        num_data_files = (int) ceil ((double) Max_Old_Frag / Frags_Per_File);
        for  (i = 1;  i <= num_data_files;  i ++)
          {
           sprintf (buffer, "%s/data%02d.olap", Input_Store_Path, i);
           fp = File_Open (buffer, "rb");
           fclose (fp);
          }
       }

   return;
  }



static void  Process_New_Batch
    (Long_Olap_Data_t new_olap [], int num_olaps)

//  Process the new overlaps in  new_olap [0 .. (num_olaps - 1)]
//  by sorting them and then combining them with any old overlaps and
//  writing them to the appropriate output files.

  {
   int  i, old_num_files;

   // First sort the new overlaps by  a_iid  using a bucket sort
   fprintf (stderr, "Sorting new overlaps\n");

   if  (New_Offset_Size < Max_ID + 2)
       {
        New_Offset_Size = Max_ID + 2;
        New_Offset = (uint32 *) Safe_realloc (New_Offset,
                                             New_Offset_Size * sizeof (uint32));
       }

   for  (i = 0;  i <= Max_ID + 1;  i ++)
     New_Offset [i] = 0;

   for  (i = 0;  i < num_olaps;  i ++)
     New_Offset [new_olap [i] . a_iid] ++;

   // Determine which files have changes
   old_num_files = Max_Old_Frag / Frags_Per_File;
   if  (Max_Old_Frag % Frags_Per_File > 0)
       old_num_files ++;
   for  (i = 1;  i <= old_num_files;  i ++)
     File_Status [i] . old_exists = TRUE;

   if  (Max_ID <= Max_Old_Frag)
       {
        Highest_Frag = Max_Old_Frag;
        Num_Files = old_num_files;
       }
     else
       {
        Highest_Frag = Max_ID;
        Num_Files = Highest_Frag / Frags_Per_File;
        if  (Highest_Frag % Frags_Per_File > 0)
            Num_Files ++;
       }
   if  (Num_Files >= File_Status_Size)
       {
        File_Status_Size *= 2;
        if  (Num_Files >= File_Status_Size)
            File_Status_Size = Num_Files + 1;
        File_Status = (File_Status_t *) Safe_realloc (File_Status,
                          File_Status_Size * sizeof (File_Status_t));
       }
   for  (i = old_num_files + 1;  i <= Num_Files;  i ++)
     File_Status [i] . old_exists = FALSE;

   for  (i = 1;  i <= Num_Files;  i ++)
     File_Status [i] . needs_update = FALSE;
   for  (i = 1;  i <= Max_ID;  i ++)
     if  (New_Offset [i] > 0)
         File_Status [1 + (i - 1) / Frags_Per_File] . needs_update = TRUE;

   for  (i = 1;  i <= Max_ID + 1;  i ++)
     New_Offset [i] += New_Offset [i - 1];

   for  (i = num_olaps - 1;  i >= 0;  i --)
     {
      int  j;

      j = -- New_Offset [new_olap [i] . a_iid];
      Sorted_Olap [j] . b_iid = new_olap [i] . b_iid;
      Sorted_Olap [j] . flipped = new_olap [i] . flipped;
      Sorted_Olap [j] . a_hang = new_olap [i] . a_hang;
      Sorted_Olap [j] . b_hang = new_olap [i] . b_hang;
      Sorted_Olap [j] . orig_erate = new_olap [i] . orig_erate;
      Sorted_Olap [j] . corr_erate = new_olap [i] . corr_erate;
     }

   Output_New_Batch ();

   return;
  }



static void  Process_Olaps_From_Dump_File
    (void)

//  Open file in  Dump_Format_Input_Path  and read its overlaps
//  and add them to the store files.  The data in the file is
//  assumed to have no duplicates and each overlap should appear
//  only once (i.e., A overlapping B should not also be reported
//  as B overlapping A.

  {
   FILE  * fp;
   Long_Olap_Data_t  fwd_olap, rev_olap;
   int  a_iid, b_iid, a_hang, b_hang;
   double  orig_erate, corr_erate;
   char  orient [100];

   fp = File_Open (Dump_Format_Input_Path, "r");

   while  (fscanf (fp, "%d %d %s %d %d %lf %lf",
                   & a_iid, & b_iid, orient, & a_hang, & b_hang,
                   & orig_erate, & corr_erate) != EOF)
     {
      fwd_olap . a_iid = a_iid;
      fwd_olap . b_iid = b_iid;
      fwd_olap . a_hang = a_hang;
      fwd_olap . b_hang = b_hang;
      fwd_olap . flipped = (tolower (orient [0]) != 'n');
      fwd_olap . orig_erate = Shrink_Quality (0.01 * orig_erate);
      fwd_olap . corr_erate = Shrink_Quality (0.01 * corr_erate);
      

      rev_olap . a_iid = b_iid;
      rev_olap . b_iid = a_iid;
      if  (fwd_olap . flipped)
          {
           rev_olap . a_hang = b_hang;
           rev_olap . b_hang = a_hang;
          }
        else
          {
           rev_olap . a_hang = - a_hang;
           rev_olap . b_hang = - b_hang;
          }
      rev_olap . flipped = fwd_olap . flipped;
      rev_olap . orig_erate = fwd_olap . orig_erate;
      rev_olap . corr_erate = fwd_olap . corr_erate;

      if  (a_iid > Max_ID)
          Max_ID = a_iid;
      if  (b_iid > Max_ID)
          Max_ID = b_iid;
      Num_Olaps += 2;

      Tmp_Output_Olap (& fwd_olap);
      Tmp_Output_Olap (& rev_olap);
     }

   return;
  }



static void  Process_OVL_Files
    (void)

//  Read the OVL files in global  OVL_File_List  and save them to
//  temp files.  When a temp file becomes full, read it, sort it,
//  and add it to the store.

  {
   time_t  now;
   int  i, num_ovl_files;

   num_ovl_files = GetNumVA_char_Ptr_t (OVL_File_List);
   for  (i = 0;  i < num_ovl_files;  i ++)
     {
      FILE  * fp;
      char  * filename;
      MesgReader  read_msg_fn;
      MessageType  imesgtype;
      GenericMesg  * pmesg;
      OverlapMesg  * ovl_mesg;

      filename = * GetVA_char_Ptr_t (OVL_File_List, i);
      fp = File_Open (filename, "r");
      read_msg_fn = InputFileType_AS (fp);

      now = time (NULL);
      fprintf (stderr, "Starting file \"%s\"\n", filename);
      fprintf (stderr, "  at %s\n", ctime (& now));


      while  (EOF != read_msg_fn (fp, & pmesg))
        {
         Long_Olap_Data_t  fwd_olap, rev_olap;

         imesgtype = pmesg -> t;
         switch  (imesgtype)
           {
            case MESG_OVL:
              ovl_mesg = pmesg -> m;
              switch (ovl_mesg -> orientation)
                {
                 case  AS_NORMAL :
                   fwd_olap . a_iid = ovl_mesg -> aifrag;
                   fwd_olap . b_iid = ovl_mesg -> bifrag;
                   fwd_olap . a_hang = ovl_mesg -> ahg;
                   fwd_olap . b_hang = ovl_mesg -> bhg;
                   fwd_olap . flipped = FALSE;
                   fwd_olap . orig_erate
                       = fwd_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   rev_olap . a_iid = ovl_mesg -> bifrag;
                   rev_olap . b_iid = ovl_mesg -> aifrag;
                   rev_olap . a_hang = - ovl_mesg -> ahg;
                   rev_olap . b_hang = - ovl_mesg -> bhg;
                   rev_olap . flipped = FALSE;
                   rev_olap . orig_erate
                       = rev_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   break;

                 case  AS_INNIE :
                   fwd_olap . a_iid = ovl_mesg -> aifrag;
                   fwd_olap . b_iid = ovl_mesg -> bifrag;
                   fwd_olap . a_hang = ovl_mesg -> ahg;
                   fwd_olap . b_hang = ovl_mesg -> bhg;
                   fwd_olap . flipped = TRUE;
                   fwd_olap . orig_erate
                       = fwd_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   rev_olap . a_iid = ovl_mesg -> bifrag;
                   rev_olap . b_iid = ovl_mesg -> aifrag;
                   rev_olap . a_hang = ovl_mesg -> bhg;
                   rev_olap . b_hang = ovl_mesg -> ahg;
                   rev_olap . flipped = TRUE;
                   rev_olap . orig_erate
                       = rev_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   break;

                 case  AS_OUTTIE :
                   fwd_olap . a_iid = ovl_mesg -> bifrag;
                   fwd_olap . b_iid = ovl_mesg -> aifrag;
                   fwd_olap . a_hang = - ovl_mesg -> ahg;
                   fwd_olap . b_hang = - ovl_mesg -> bhg;
                   fwd_olap . flipped = TRUE;
                   fwd_olap . orig_erate
                       = fwd_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   rev_olap . a_iid = ovl_mesg -> aifrag;
                   rev_olap . b_iid = ovl_mesg -> bifrag;
                   rev_olap . a_hang = - ovl_mesg -> bhg;
                   rev_olap . b_hang = - ovl_mesg -> ahg;
                   rev_olap . flipped = TRUE;
                   rev_olap . orig_erate
                       = rev_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   break;

                 case  AS_ANTI :

                   fwd_olap . a_iid = ovl_mesg -> aifrag;
                   fwd_olap . b_iid = ovl_mesg -> bifrag;
                   fwd_olap . a_hang = - ovl_mesg -> bhg;
                   fwd_olap . b_hang = - ovl_mesg -> ahg;
                   fwd_olap . flipped = FALSE;
                   fwd_olap . orig_erate
                       = fwd_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   rev_olap . a_iid = ovl_mesg -> bifrag;
                   rev_olap . b_iid = ovl_mesg -> aifrag;
                   rev_olap . a_hang = ovl_mesg -> bhg;
                   rev_olap . b_hang = ovl_mesg -> ahg;
                   rev_olap . flipped = FALSE;
                   rev_olap . orig_erate
                       = rev_olap . corr_erate = Shrink_Quality (ovl_mesg -> quality);

                   break;

                 case  AS_UNKNOWN :
                 default :
                   fprintf (stderr,
                   "YIKES:  Bad overlap orientation = %d for a = %d  b = %d\n",
                            (int) ovl_mesg -> orientation,
                            ovl_mesg -> aifrag, ovl_mesg -> bifrag);
                   continue;
                }

              if  (ovl_mesg -> aifrag > Max_ID)
                  Max_ID = ovl_mesg -> aifrag;
              if  (ovl_mesg -> bifrag > Max_ID)
                  Max_ID = ovl_mesg -> bifrag;
              Num_Olaps += 2;

              Tmp_Output_Olap (& fwd_olap);
              Tmp_Output_Olap (& rev_olap);

              break;

            default :
              ;  // Ignore other types of messages
           }

        }

      fclose (fp);
     }

   return;
  }


static void  Process_Tmp_File
    (int idx)

//  Process the new overlaps saved in the temporary file referenced
//  in  Tmp_File_Info [idx] .  More specifically:  read them; sort them;
//  and merge them with any overlaps previous saved in the overlap store.

  {
   Tmp_File_Info_t  * tfi;

   tfi = Tmp_File_Info + idx;

   rewind (tfi -> fp);
   assert (tfi -> ct <= MAX_OLAP_BATCH);

   if  (Verbose > 0)
       fprintf (stderr, "Process_Tmp_File #%d  ct = %d\n", idx, tfi -> ct);
   Safe_fread (New_Olap, sizeof (Long_Olap_Data_t), tfi -> ct, tfi -> fp);

   Process_New_Batch (New_Olap, tfi -> ct);

   rewind (tfi -> fp);
   tfi -> ct = 0;

   return;
  }



static void  Tmp_Output_Olap
    (const Long_Olap_Data_t * olap)

//  Write  (* olap)  to its appropriate temporary file, kept
//  track of in global  Tmp_File_Info  arrray.  If that file
//  is full, then sort and merge all the overlaps in it together
//  with any prior overlaps for that range, and clear the temporary
//  file.

  {
   char  file_name [MAX_FILENAME_LEN];
   Tmp_File_Info_t  * tfi;
   int  file_index;
   int  i;

   file_index = (olap -> a_iid + Frags_Per_File - 1) / Frags_Per_File;

   if  (file_index >= Tmp_File_Info_Size)
       {
        Tmp_File_Info_Size *= 2;
        if  (file_index >= Tmp_File_Info_Size)
            Tmp_File_Info_Size = file_index + 1;
        Tmp_File_Info = (Tmp_File_Info_t *) Safe_realloc (Tmp_File_Info,
                          Tmp_File_Info_Size * sizeof (Tmp_File_Info_t));
       }

   if  (file_index > Num_Files)
       {
        for  (i = Num_Files + 1;  i <= file_index;  i ++)
          {
           Tmp_File_Info [i] . fp = NULL;
           Tmp_File_Info [i] . ct = 0;
          }
        Num_Files = file_index;
       }

   tfi = Tmp_File_Info + file_index;

   if  (tfi -> fp == NULL)
       {
        sprintf (file_name, "%s/new%02d.tmp", Output_Store_Path, file_index);
        assert (strlen (file_name) < MAX_FILENAME_LEN - 1);
        tfi -> fp = File_Open (file_name, "wb+");
       }

   Safe_fwrite (olap, sizeof (Long_Olap_Data_t), 1, tfi -> fp);
   tfi -> ct ++;

   if  (tfi -> ct >= MAX_OLAP_BATCH)
       Process_Tmp_File (file_index);

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s [-aAcfhS] [-i <InputStore> [-o <OutputStore>]\n"
       "        [-L <ListFile>] [-D <DumpFile>] [-v <VerboseLevel>]\n"
       "        <OlapFile1> [<OlapFile2> ...]\n"
       "\n"
       "Create or add to a compact, sorted, binary-file of fragment\n"
       "overlap information\n"
       "\n"
       "Options:\n"
       "-i   specifies input olap store\n"
       "-o   specifies output olap store\n"
       "-a   append with backup (requires -i)\n"
       "-A   append without backup (requires -i)\n"
       "-c   create (requires -o)\n"
       "-D   specify a file containing overlaps in dump-olap-store format\n"
       "-f   force, if output exists delete it\n"
       "-h   help, print this message\n"
       "-L   specify a file containing names of overlap files\n"
       "     (useful if too many to fit on command line)\n"
       "-S   save space by renaming each tmp file as it is finished\n"
       "-v   set verbose level for extra printouts\n",
       command);

   return;
  }

