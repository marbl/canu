
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
* Module:  MakeScriptOVL.c
* Description:
*   Create a batch script of LSF commands to compute the overlaps
*   specified on the command line.
* 
*    Programmer:  A. Delcher
*************************************************/

/* RCS info
 * $Id: MakeScriptOVL.c,v 1.6 2007-01-29 20:41:17 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char  CM_ID []
    = "$Id: MakeScriptOVL.c,v 1.6 2007-01-29 20:41:17 brianwalenz Exp $";


/*************************************************************************/
/* System include files */
/*************************************************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>

/*************************************************************************/
/* Local include files */
/*************************************************************************/

#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"

# define THE_SHELL "/usr/bin/csh"

#define  AS_BIN_VARIABLE        "AS_BIN"
    //  Standard environment variable name for location of binaries
#define  ASCII_APPEND_COMMAND   "proto_append"
    //  Command to use to append overlap files if ascii mode
#define  BINARY_APPEND_COMMAND  "binary_append"
    //  Command to use to append overlap files if binary mode
#define  DEFAULT_LSF_PROJECT_NAME   "IR:OPERATIONS:M"
    //  Tag used for project name and as prefixes for job names and
    //  output logs.
#define  DEFAULT_SCRIPTFILE_NAME  "fragovl.script"
    //  Name of scriptfile produced if not specified on the
    //  command line
#define  JOB_GROUP                "/Assembly"
    //  Prefix used for jobname in LSF commands
#define  POPULATOR_CMD       "PopulateFragStore"
    //  Command for first job which populates the fragment store and
    //  computes no overlaps
#define  LSF_QUEUE_NAME       "assembly"
    //  LSF queue to which jobs are submitted
//#define  MACHINE_SPEC_STRING  "select[mem>=1100 && ncpus>=4] span[hosts=1]"
#define  MACHINE_SPEC_STRING  "rusage[physmem=2800]"
    //  LSF option to specify machine characteristics needed to run job
#define  MAX_LSF_ARRAY_JOBS       2000
    //  Most LSF jobs allowed in an LSF array.  Will adjust batch size
    //  to ensure that at most this many array jobs are created.
//#define  PROCESSOR_SPEC_OPTION  "-x -n 2,4"
#define  PROCESSOR_SPEC_OPTION  "-x"
    //  LSF option to specify number of available processors needed to
    //  run job
#define  TAG_STRING             "lsfovl"
    //  Tag to use on lsf script file names and jobids

#if  1
#define  MIN_BATCH_SIZE       800000
    //  Make new batch plus old frags at least this big if possible
#define  NEW_BATCH_SIZE       400000
    //  Ideal number of new fragments to process in each LSF job
#define  OLD_BATCH_SIZE       2000000
    //  Ideal number of old fragments to process in each LSF job
#else
#define  MIN_BATCH_SIZE       30000
    //  Make new batch plus old frags at least this big if possible
#define  NEW_BATCH_SIZE       10000
    //  Ideal number of new fragments to process in each LSF job
#define  OLD_BATCH_SIZE       50000
    //  Ideal number of old fragments to process in each LSF job
#endif

/*************************************************************************/
/* Type definitions */
/*************************************************************************/



//  Constants

const char *  DEFAULT_OLAP_CMD = "overlap";
    //  Command for subsequent jobs which read the fragment store and
    //  produce overlaps


/*************************************************************************/
/* Static Globals */
/*************************************************************************/

static char  * Append_Command = BINARY_APPEND_COMMAND;
static char  * Input_Store_Name = NULL;
static char  * Output_Store_Name = NULL;

char  AS_BIN_Path [MAX_NAME_LEN];
    //  Holds default path to overlap command if commend not specified
    //  by  -C  option
static Batch_ID  Batch_Msg_UID = 0;
static char  * E_Option = NULL;
    //  Indicates path to check with -E option of bsub command
static char  * LSF_Project_Name = DEFAULT_LSF_PROJECT_NAME;
static int  Concat_OVL_Files = FALSE;
    //  If set true then create script to concatenate all
    //  .ovl files into a single .ovl file
static int  Even_Batch_Sizes = FALSE;
    //  If set true by  -e  option, make sizes even multiples
    //  of  OLD_BATCH_SIZE  and  NEW_BATCH_SIZE  as much as possible,
    //  i.e., instead of making sizes equal.
static int  First_Is_Rename = FALSE;
static const char  * Olap_Cmd;
    //  Command to do overlapping; can be set by  -C  option
static int  Use_LSF = TRUE;



/*************************************************************************/
/* External Global Definitions */
/*************************************************************************/



/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/

static int  Count_Frag_Messages
    (Input_Stream in_stream);
static void  Emit_Script
    (char * scriptfile_name, char * lsf_queue_name, char * outfile_prefix,
     char * outfile_name, int old_ct, int new_ct,
     char * store_path, char * options);
static int  Input_File_Is_Binary
    (FILE * fp);
static void  Usage
    (char * executable_name);

int  main  (int argc, char * argv [])

  {
   FILE  * fp, * range_file;
   char  * outfile_name, * outfile_prefix, * scriptfile_name;
   char  * fragstore_path = NULL, * lsf_queue_name;
   char  option_string [MAX_NAME_LEN] = {'\0'};
   int64  first_frag, last_frag, last_stored_frag;
   int  use_ascii_output = FALSE;
   int  new_frag_ct, old_frag_ct;

   Olap_Cmd = DEFAULT_OLAP_CMD;

   // Parse the argument list using "man 3 getopt".

   outfile_prefix = NULL;
   scriptfile_name = NULL;
   lsf_queue_name  = LSF_QUEUE_NAME;
   {
    int  ch, errflg = 0;
    optarg = NULL;

    while  (! errflg
              && ((ch = getopt (argc, argv, "1C:eE:hk:K:l:L:o:Pq:Qs:t:w")) != EOF))
      switch  (ch)
        {
         case  '1' :
           Concat_OVL_Files = TRUE;
           break;
         case  'C' :
           Olap_Cmd = strdup (optarg);
           break;
         case  'e' :
           Even_Batch_Sizes = TRUE;
           break;
         case  'E' :
           E_Option = optarg;
           break;
         case  'h' :
           Usage (argv [0]);
           exit (EXIT_FAILURE);
           break;
         case  'k' :
           strcat (option_string, " -k ");
           strcat (option_string, optarg);
           break;
         case  'K' :
           strcat (option_string, " -K \"");
           strcat (option_string, optarg);
           strcat (option_string, "\"");
           break;
         case  'l' :
           strcat (option_string, " -l ");
           strcat (option_string, optarg);
           break;
         case  'L' :
           LSF_Project_Name = strdup (optarg);
           break;
         case  'o' :
           outfile_prefix = strdup (optarg);
           break;
         case  'P' :
           strcat (option_string, " -P");
           Append_Command = ASCII_APPEND_COMMAND;
           use_ascii_output = TRUE;
           break;
         case  'q' :
           lsf_queue_name = strdup (optarg);
           break;
         case  'Q' :
           Use_LSF = FALSE;
           break;
         case  's' :
           scriptfile_name = strdup (optarg);
           break;
         case  't' :
           strcat (option_string, " -t ");
           strcat (option_string, optarg);
           break;
         case  'w' :
           strcat (option_string, " -w");
           break;
         case  '?' :
           fprintf (stderr, "Unrecognized option -%c", optopt);
         default :
           errflg++;
        }

    if  (argc - optind != 2)
        {
         Usage (argv [0]);
         exit (EXIT_FAILURE);
        }


    fragstore_path = argv [optind ++];
    range_file = File_Open (argv [optind], "r");
    if  (fscanf (range_file, "batch:" F_U64 "\n", & Batch_Msg_UID) != 1
           || fscanf (range_file, "lo:" F_S64 "  hi:" F_S64 "\n", & first_frag, & last_frag)
                  != 2)
        {
         fprintf (stderr, "ERROR in range file \"%s\"\n", argv [optind]);
         exit (EXIT_FAILURE);
        }
    fclose (range_file);

    if  (outfile_prefix == NULL)
        {
         char  * p;

         outfile_prefix = strdup (argv [optind]);
         p = strrchr (outfile_prefix, '.');
         if  (p != NULL)
             * p = '\0';
        }
   }

   assert (outfile_prefix != NULL);
   outfile_name = (char *) malloc (strlen (outfile_prefix) + 5);
   assert (outfile_name != NULL);
#if 0 // A memory leak.
   outfile_name = strdup (outfile_prefix);
   strcat (outfile_name, ".ovl");
#else
   sprintf(outfile_name, "%s.ovl", outfile_prefix);
#endif
   fp = fopen (outfile_name, "r");

   if  (fp != NULL)
       {
        // check if binary or proto
        if  (Input_File_Is_Binary (fp))
            {
             if  (use_ascii_output)
                 {
                  fprintf (stderr, "ERROR: specified -P when %s is binary\n",
                           outfile_name);
                  exit (EXIT_FAILURE);
                 }
            }
          else
            {
             if  (! use_ascii_output)
                 {
                  fprintf (stderr, "ERROR: need -P when %s is ASCII (proto)\n",
                           outfile_name);
                  exit (EXIT_FAILURE);
                 }
            }
        fclose (fp);
       }
     else
       {
        // do rename instead of append first time
        First_Is_Rename = TRUE;
       }


   if  (existsFragStore (fragstore_path))
       {
        FragStore  frag_store = openFragStore (fragstore_path, "r");

        last_stored_frag = getLastElemFragStore (frag_store);
        closeFragStore (frag_store);
       }
     else
       {
        fprintf (stderr, "ERROR:  Could not open fragstore \"%s\"\n",
                 fragstore_path);
        exit (EXIT_FAILURE);
       }

   if  (last_stored_frag < last_frag)
       {
        fprintf (stderr,
                 "ERROR:  fragment range " F_S64 ".." F_S64 " goes beyond last frag " F_S64 "\n",
                 first_frag, last_frag, last_stored_frag);
        exit (EXIT_FAILURE);
       }

   if  (scriptfile_name == NULL)
       scriptfile_name = DEFAULT_SCRIPTFILE_NAME;


   old_frag_ct = first_frag - 1;
   new_frag_ct = last_frag - old_frag_ct;

   Emit_Script (scriptfile_name, lsf_queue_name, outfile_prefix,
                outfile_name, old_frag_ct, new_frag_ct,
                fragstore_path, option_string);

   printf ("Overlap command is \"%s%s\"\n", AS_BIN_Path, Olap_Cmd);
   printf ("fragment range = " F_S64 ".. " F_S64 "\n", first_frag, last_frag);
   printf ("option_string = \"%s\"\n", option_string);
   printf ("fragstore_path = \"%s\"\n", fragstore_path);
   printf ("outfile_prefix = \"%s\"\n", outfile_prefix);
   printf ("scriptfile_name = \"%s\"\n", scriptfile_name);

   fprintf (stderr, "### Return from main\n");

   return  0;
  }



static int  Count_Frag_Messages
    (Input_Stream in_stream)

//  Return the number of  AS_ADD  fragment messages in  in_stream .

  {
   GenericMesg  * pmesg;
   int  ct = 0;

   while  (ReadProtoMesg_AS (in_stream, & pmesg) != EOF)
     {
      MessageType  imesgtype = pmesg -> t;

      switch  (imesgtype)
        {
         case  MESG_IFG :
           {
            InternalFragMesg  * ifg = (InternalFragMesg*) pmesg -> m;

            if  (ifg -> action == AS_ADD)
                ct ++;
 
            break;
           }

         default :
           ;  // Nothing
        }
     }

   fclose (in_stream);

   return  ct;
  }



static void  Emit_Script
    (char * scriptfile_name, char * lsf_queue_name, char * outfile_prefix,
     char * outfile_name, int old_ct, int new_ct,
     char * store_path, char * options)

//  Output to file  scriptfile_name  a set of LSF commands
//  to compute overlaps between the first  old_ct  fragments in the fragment store
//  store_path  and the following  new_ct  fragments after them
//  in the same fragment store.
//  The output overlap file will have prefix  outfile_prefix .
//  The name of the single output file (if  Concat_OVL_Files  is true) is
//  outfile_name .
//  The LSF queue to which the jobs are
//  submitted is  lsf_queue_name .
//  Options to overlap commands are in  options .

  {
   FILE  * scriptfile, * appendfile;
   char  id_tag [MAX_NAME_LEN];
   char  job_name [MAX_NAME_LEN];
   char  working_dir [MAX_NAME_LEN];
   char  append_filename [MAX_NAME_LEN];
   char  olap_script_name [MAX_NAME_LEN];
   char  command [MAX_NAME_LEN];
   char  * as_bin_env;
   int  new_batches, tag_ct;
   double  new_incr;
   int  last, i, j;

   tag_ct = 0;
   sprintf (id_tag, "%s" F_PID_T, TAG_STRING, getpid ());
   sprintf (working_dir, "OlapDir" F_PID_T, getpid ());
   sprintf (command, "mkdir %s", working_dir);
   system (command);

   scriptfile = File_Open (scriptfile_name, "w");
   fprintf (scriptfile, "#!" THE_SHELL "\n");
   fprintf (scriptfile, "set QUEUE = %s\n", lsf_queue_name);
   fprintf (scriptfile, "set MACHINE_SPEC = \"%s\"\n", MACHINE_SPEC_STRING);
   if  (E_Option != NULL)
       fprintf (scriptfile, "set CHECK_PATH = \"%s\"\n", E_Option);

   // identify the AS_BIN path, if set in the environment
   if  ( (as_bin_env = getenv( AS_BIN_VARIABLE )) != NULL
           && Olap_Cmd == DEFAULT_OLAP_CMD)
       sprintf (AS_BIN_Path, "%s/", as_bin_env);
   else
       AS_BIN_Path[0] = '\0';
   
   if  (new_ct <= 0)
       {
        fprintf (stderr, "No new fragments.  Nothing to do\n");
        fprintf (scriptfile, "echo \"No new fragments.  Nothing to do\"");
        fclose (scriptfile);

        sprintf (command, "chmod a+x %s", scriptfile_name);
        system (command);

        return;
       }

   if  (! Use_LSF)
       {
        fprintf (scriptfile, "%s%s -h %7d-%-7d -r %7d-%-7d %s"
                 " -o %s.%dTMP.ovl %s\n",
                 AS_BIN_Path,
                 Olap_Cmd, old_ct + 1, old_ct + new_ct, 1, old_ct + new_ct,
                 options, outfile_prefix, 1, store_path);
        if  (! Concat_OVL_Files)
            fprintf (scriptfile, "mv %s.%dTMP.ovl %s.r%d-%dh%d-%d.ovl\n",
                     outfile_prefix, 1, outfile_prefix,
                     1, old_ct + new_ct, old_ct + 1, old_ct + new_ct);
          else
            {
             if  (First_Is_Rename)
                 fprintf (scriptfile, "mv %s.%dTMP.ovl %s\n",
                          outfile_prefix, 1, outfile_name);
               else
                 {
                  fprintf (scriptfile, "%s %s %s.%dTMP.ovl\n",
                           Append_Command, outfile_name, outfile_prefix, 1);
                  fprintf (scriptfile, "rm %s.%dTMP.ovl\n", outfile_prefix, 1);
                 }
            }
        fclose (scriptfile);

        sprintf (command, "chmod a+x %s", scriptfile_name);
        system (command);

        return;
       }

   new_batches = 1 + (new_ct - 1) / NEW_BATCH_SIZE;    // Clark's version
   if  (new_batches > MAX_LSF_ARRAY_JOBS)
       new_batches = MAX_LSF_ARRAY_JOBS;
   if  (Even_Batch_Sizes)
       {
        new_incr = ceil ((double) new_ct / (new_batches * NEW_BATCH_SIZE))
                     * NEW_BATCH_SIZE;
        new_batches = ceil (new_ct / new_incr);
       }
     else
       new_incr = ((double) new_ct) / new_batches;

   last = old_ct;
   for  (i = 0;  i < new_batches;  i ++)
     {
      int  old_batches;
      double  old_incr;
      int  first, lo;

      first = last + 1;
      last = rint (old_ct + (i + 1) * new_incr);

      // "Double up" first batch if too few old fragments
      while  (last < MIN_BATCH_SIZE && i < new_batches - 1)
        {
         i ++;
         last = rint (old_ct + (i + 1) * new_incr);
        }

      // Set last directly for last batch to eliminate round-off errors
      if  (i >= new_batches - 1)
          last = old_ct + new_ct;

      old_batches = 1 + (last - 1) / OLD_BATCH_SIZE;    // Clark's version
      if  (old_batches > MAX_LSF_ARRAY_JOBS / new_batches)
          old_batches = MAX_LSF_ARRAY_JOBS / new_batches;

      if  (Even_Batch_Sizes)
          {
           old_incr = ceil ((double) last / (old_batches * OLD_BATCH_SIZE))
                        * OLD_BATCH_SIZE;
           old_batches = ceil (last / old_incr);
          }
        else
          old_incr = ((double) last) / old_batches;

      lo = 1;
      for  (j = 0;  j < old_batches;  j ++)
        {
         FILE  * olap_script_fp;
         int  hi;

         // Set hi directly for last batch to eliminate round-off errors
         if  (j >= old_batches - 1)
             hi = last;
           else
             hi = rint ((j + 1) * old_incr);

         sprintf (job_name, "%s.%d", id_tag, ++ tag_ct);
         sprintf (olap_script_name, "%s/%s.sh", working_dir, job_name);

         if  (first < lo)
             first = lo;

         olap_script_fp = File_Open (olap_script_name, "w");
         fprintf (olap_script_fp, "#!" THE_SHELL "\n");
         fprintf (olap_script_fp, "echo \"path = \" $path\n");
         fprintf (olap_script_fp, "set QUEUE = %s\n", lsf_queue_name);
         fprintf (olap_script_fp, "set MACHINE_SPEC = \"%s\"\n",
                  MACHINE_SPEC_STRING);
         fprintf (olap_script_fp, "%s%s -h %7d-%-7d -r %7d-%-7d %s"
                  " -o %s.%dTMP.ovl %s\n",
                  AS_BIN_Path,
                  Olap_Cmd, first, last, lo, hi, options,
                  outfile_prefix, tag_ct, store_path);
         fprintf (olap_script_fp, "if  (${status} == 0)  then\n");
         fprintf (olap_script_fp, "echo \"olap job %s successful\"\n",
                  job_name);
         if  (! Concat_OVL_Files)
             fprintf (olap_script_fp, "mv %s.%dTMP.ovl %s.r%d-%dh%d-%d.ovl\n",
                      outfile_prefix, tag_ct, outfile_prefix,
                      lo, hi, first, last);
           else
             fprintf (olap_script_fp, "mv %s.%dTMP.ovl %s.%d.ovl\n",
                      outfile_prefix, tag_ct, outfile_prefix, tag_ct);
         fprintf (olap_script_fp, "bdel -J \"%s/%sFAIL\"\n",
                  JOB_GROUP, job_name);
         fprintf (olap_script_fp, "else\n");
         fprintf (olap_script_fp, "echo \"olap job %s failed  status = \" ${status}\n",
                  job_name);
         fprintf (olap_script_fp,
                  "echo \"This email was sent automatically by the LSF overlapper\""
                  " >! %s/%s.FAILED\n",
                  working_dir, job_name);
         fprintf (olap_script_fp,
                  "echo \"olap job %s failed\" >> %s/%s.FAILED\n",
                  job_name, working_dir, job_name);
         fprintf (olap_script_fp, "mail ${user} < %s/%s.FAILED\n",
                  working_dir, job_name);
         fprintf (olap_script_fp, "bdel -J \"%s/%sSUCCESS\"\n",
                  JOB_GROUP, id_tag);
         fprintf (olap_script_fp, "bdel -J \"%s/%sFAIL\"\n",
                  JOB_GROUP, job_name);
         fprintf (olap_script_fp, "exit -1\n");
         fprintf (olap_script_fp, "endif\n");
         fclose (olap_script_fp);
         sprintf (command, "chmod a+x %s", olap_script_name);
         system (command);


         sprintf (olap_script_name, "%s/%sFAIL.sh", working_dir, job_name);
         olap_script_fp = File_Open (olap_script_name, "w");
         fprintf (olap_script_fp, "#!" THE_SHELL "\n");
         fprintf (olap_script_fp, "set QUEUE = %s\n", lsf_queue_name);
         fprintf (olap_script_fp, "set MACHINE_SPEC = \"%s\"\n",
                  MACHINE_SPEC_STRING);
         fprintf (olap_script_fp, "echo \"olap job %s failed\"\n",
                  job_name);
         fprintf (olap_script_fp,
                  "echo \"This email was sent automatically by the LSF overlapper\""
                  " >! %s/%s.FAILED\n",
                  working_dir, job_name);
         fprintf (olap_script_fp,
                  "echo \"olap job %s failed\" >> %s/%s.FAILED\n",
                  job_name, working_dir, job_name);
         fprintf (olap_script_fp, "mail ${user} < %s/%s.FAILED\n",
                  working_dir, job_name);
         fprintf (olap_script_fp, "bdel -J \"%s/%sSUCCESS\"\n",
                  JOB_GROUP, id_tag);
         fclose (olap_script_fp);
         sprintf (command, "chmod a+x %s", olap_script_name);
         system (command);


         lo = hi + 1;
        }
     }
            
   if  (tag_ct > 0)
       {
        char  * option_string = "";

        if  (E_Option != NULL)
            option_string = "-E \"ls -l $CHECK_PATH\" ";

        fprintf (scriptfile,
                 "set bsubReturn = `bsub %s-q $QUEUE %s -R \"$MACHINE_SPEC\""
                 " -J \"%s/%s[1-%d:1]\" -o \"%s.%%I.log\""
                 " -P \"%s\" -i \"%s/%s.%%I.sh\" " THE_SHELL " -s`\n",
                 option_string,
                 PROCESSOR_SPEC_OPTION, JOB_GROUP, id_tag, tag_ct, id_tag,
                 LSF_Project_Name, working_dir, id_tag);
        fprintf (scriptfile,
                 "set olapJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Overlap LSF ID = $olapJobID\"\n");
        for  (i = 1;  i <= tag_ct;  i ++)
          fprintf (scriptfile,
                   "bsub -q $QUEUE -J \"%s/%s.%dFAIL\" -w \"exit(${olapJobID}[%d])\""
                   " -o \"%s.%dFAIL.log\" -P \"%s\" %s/%s.%dFAIL.sh\n",
                   JOB_GROUP, id_tag, i,
                   i, id_tag, i, LSF_Project_Name, working_dir, id_tag, i);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -K -q $QUEUE -J \"%s/%sSUCCESS\""
                 " -w \"done(${olapJobID}[1-%d])\""
//                 " -o \"%sSUCCESS.log\" \"echo Success; rm -rf %s\"`\n",
                 " -o \"%sSUCCESS.log\" -P \"%s\" \"echo Success\"`\n",
                 JOB_GROUP, id_tag,
//                 tag_ct, id_tag, working_dir);
                 tag_ct, id_tag, LSF_Project_Name);
        fprintf (scriptfile,
                 "set successJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Success LSF ID = $successJobID\"\n");

        // Kludge because LSF leaves a pending job
        fprintf (scriptfile,
                 "bdel ${successJobID}\n");
       }

   if  (Concat_OVL_Files)
       {
        sprintf (append_filename, "%s.append", id_tag);
        appendfile = File_Open (append_filename, "w");
        for  (i = 1; i <= tag_ct;  i ++)
          {
           if  (First_Is_Rename && i == 1)
               fprintf (appendfile, "mv %s.%d.ovl %s\n",
                        outfile_prefix, i, outfile_name);
             else
               {
                fprintf (appendfile, "%s %s %s.%d.ovl\n",
                         Append_Command, outfile_name, outfile_prefix, i);
                fprintf (appendfile, "rm %s.%d.ovl\n", outfile_prefix, i);
               }
          }
        fclose (appendfile);
        sprintf (command, "chmod a+x %s", append_filename);
        system (command);

#if  0
        sprintf (job_name, "%sappend", id_tag);
        fprintf (scriptfile,
                 "bsub -K -q $QUEUE -J \"/%s\" -o \"%s.log\" -P \"%s\""
                 " -w \"%s/%sSUCCESS\" %s\n",
                 job_name, job_name, LSF_Project_Name, JOB_GROUP, id_tag,
                 append_filename);

        // Kludge because LSF leaves a pending job
        fprintf (scriptfile,
                 "bdel -J \"/%s\"\n", job_name);
#else
        // Just run the append without LSF  -K on success job
        // should ensure that this will not run before it.
        fprintf (scriptfile,
                 "%s\n", append_filename);
#endif
       }

   fclose (scriptfile);

   sprintf (command, "chmod a+x %s", scriptfile_name);
   system (command);

   return;
  }



static int  Input_File_Is_Binary
    (FILE * fp)

//  Returns  TRUE  iff  fp  is a binary-form proto I/O file.
//  Assumes  fp  is at the beginning of a valid message file.

  {
   int  c;

   c = fgetc (fp);
   ungetc (c, fp);
   return  (c == 0);
  }



static void  Usage
    (char * executable_name)

//  Print parameters and options to run this program, named
//   executable_name .

  {
   fprintf (stderr,
           "USAGE:  %s [-1ehPQw] [-o <ovlfiletag>]\n"
           "        [-K <kmerhitlimit>] [-q <LSFqueue>]\n"
           "        [-L <LSFprojectname>] [-C <overlapper>]\n"
           "        [-s <scriptfile>] [-t <numThreads>]\n"
           "        <FragStore> <rangefile>\n"
           "Creates LSF script to compute overlaps\n"
           "If  -o  is not specified, output from script will be  <tag>.ovl\n"
           "where  <tag>  is the part of  <rangefile>  preceding the\n"
           "last dot\n"
           "Options:\n"
           "-1   concatenate all .ovl files together\n"
           "-C   specify overlap command to use, *MUST* specify full path\n"
           "     default is  $AS_BIN/overlap\n"
           "-e   use even-multiple batch sizes\n"
           "-h   print this message\n"
           "-k   specify file of kmers to skip in hash table\n"
           "-K   specify min hash-table kmer frequency to ignore\n"
           "     can use  [g,c,p]  format to set K value automatically\n"
           "     so that the probability of K or more hits in a genome\n"
           "     of length g containing c copies of the kmer is\n"
           "     less than p  (Put in quotes to get past the shell.)\n"
           "     e.g.,  -K \"[3e9,10,1e-9]\"\n"
           "-l   specifies max olaps per old frag end per hash batch\n"
           "-L   specifies LSF project name\n"
           "-o   specifies prefix of ovl output file names\n"
           "-P   proto (ASCII) output (default is binary)\n"
           "-q   specifies LSF queue for jobs\n"
           "-Q   *DON'T* use LSF\n"
           "-s   specifies name of LSF scriptfile\n"
           "-t   specifies number of parallel threads for overlap jobs\n"
           "-w   filter overlaps with too many errors in a window\n",
           executable_name);
  }

