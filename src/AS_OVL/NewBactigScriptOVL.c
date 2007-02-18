
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
* Module:  NewBactigScriptOVL.c
* Description:
*   Create a batch script of LSF commands to compute overlaps
*   between all fragments in the specified frag store and the
*   range of BACtigs in the specified BACtig store.  The range
*   is specified in the range-file named on the command line.
* 
*    Programmer:  A. Delcher
*       Written:  15 Feb 00
*  Last Revised:  15 Feb 00
* 
*************************************************/

/* RCS info
 * $Id: NewBactigScriptOVL.c,v 1.6 2007-02-18 14:04:50 brianwalenz Exp $
 * $Revision: 1.6 $
*/

static char  CM_ID []
    = "$Id: NewBactigScriptOVL.c,v 1.6 2007-02-18 14:04:50 brianwalenz Exp $";


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
#include  "AS_PER_gkpStore.h"
#include  "AS_PER_genericStore.h"
#include  "AS_PER_fragStore.h"
#include  "AS_PER_distStore.h"
#include  "AS_UTL_PHash.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_overlap.h"


#define  DEFAULT_LSF_PROJECT_NAME  "IR:HUM_ASM:M"
    //  Name of LSF project if none specified with  -L  option
#define  DEFAULT_SCRIPTFILE_NAME  "bactigovl.script"
    //  Name of scriptfile produced if not specified on the
    //  command line
#define  JOB_GROUP                "/Assembly"
    //  Prefix used for jobname in LSF commands
#define  OLAP_CMD                 "$OVERLAY/bin/overlap_ca -c -m"
    //  Command for subsequent jobs which read the fragment store and
    //  produce overlaps
#define  PROJECT_TAG             "newbac"
    //  Tag used for prefixes for job names and output logs.
#define  LSF_QUEUE_NAME           "assembly"
    //  LSF queue to which jobs are submitted
//#define  MACHINE_SPEC_STRING      "select[mem>=1100 && ncpus>=4] span[hosts=1]"
//#define  MACHINE_SPEC_STRING      "rusage[physmem=2500] select[ncpus>=4] span[hosts=1]"
#define  MACHINE_SPEC_STRING      "rusage[physmem=2500]"
    //  LSF option to specify machine characteristics needed to run job
#define  MAX_LSF_ARRAY_JOBS       2000
    //  Most LSF jobs allowed in an LSF array.  Will adjust batch size
    //  to ensure that at most this many array jobs are created.
//#define  PROCESSOR_SPEC_OPTION    "-n 2,4"
//#define  PROCESSOR_SPEC_OPTION    " -x -n 2,4"
#define  PROCESSOR_SPEC_OPTION    " -x"
    //  LSF option to specify number of available processors needed to
    //  run job

#if  1
#define  NEW_BATCH_SIZE           2000
    //  Ideal number of new BACtigs to process in each LSF job
#define  OLD_BATCH_SIZE           500000
    //  Ideal number of old fragments to process in each LSF job
#else
#define  NEW_BATCH_SIZE           2000
    //  Ideal number of new BACtigs to process in each LSF job
#define  OLD_BATCH_SIZE           100000
    //  Ideal number of old fragments to process in each LSF job
#endif

/*************************************************************************/
/* Type definitions */
/*************************************************************************/



/*************************************************************************/
/* Static Globals */
/*************************************************************************/

static Batch_ID  Batch_Msg_UID = 0;
static char  * LSF_Project_Name = DEFAULT_LSF_PROJECT_NAME;


/*************************************************************************/
/* External Global Definitions */
/*************************************************************************/


/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/

static void  Emit_Script
    (char * scriptfile_name, char * lsf_queue_name,
     int old_ct, int64 first_new, int64 last_new,
     char * fragstore_path, char * bactigstore_path, char * bol_path,
     char * options);




int  main  (int argc, char * argv [])

  {
   FILE  * range_file;
   char  * outfile_name, * scriptfile_name, * lsf_queue_name;
   char  * bactigstore_path, * fragstore_path, * bol_path;
   char  option_string [MAX_NAME_LEN] = {'\0'};
   int  exists;
   int64  first_bactig, last_bactig, last_stored_bactig;
   int  old_frag_ct;

   exists = 0;


   // Parse the argument list using "man 3 getopt".

   outfile_name = NULL;
   scriptfile_name = NULL;
   lsf_queue_name  = LSF_QUEUE_NAME;
   {
    int  ch, errflg = 0;
    optarg = NULL;

    while  (! errflg && ((ch = getopt (argc, argv, "L:Ps:t:q:")) != EOF))
      switch  (ch)
        {
         case  'L' :
           LSF_Project_Name = strdup (optarg);
           break;
         case  'P' :
           strcat (option_string, " -P");
           break;
         case  'q' :
           lsf_queue_name = strdup (optarg);
           break;
         case  's' :
           scriptfile_name = strdup (optarg);
           break;
         case  't' :
           strcat (option_string, " -t ");
           strcat (option_string, optarg);
           break;
         case  '?' :
           fprintf (stderr, "Unrecognized option -%c", optopt);
         default :
           errflg++;
        }

    if  (argc - optind != 4)
        {
         fprintf (stderr,
                 "USAGE:  %s [options] <FragStore> <BactigStore>\n"
                 "          <BOLPath> <RangeFile>\n"
                 "Creates script to overlap all of FragStore with specified\n"
                 "  range of BactigStore\n"
                 "Output .bol files go to <BOLPath>\n"
                 "Use -L to specify LSF project name\n"
                 "Use -P to force ASCII output\n"
                 "Use -q to specify name of LSF queue\n"
                 "Use -s to specify name of script file\n"
                 "Use -t to designate number of parallel threads\n",
                 argv [0]);
         exit (EXIT_FAILURE);
        }

    fragstore_path = argv [optind ++];
    bactigstore_path = argv [optind ++];
    bol_path = argv [optind ++];

    range_file = File_Open (argv [optind], "r");
    if  (fscanf (range_file, "batch:" F_UID "\n", & Batch_Msg_UID) != 1
           || fscanf (range_file, "lo:" F_S64 "  hi:" F_S64 "\n",
                      & first_bactig, & last_bactig)
                  != 2)
        {
         fprintf (stderr, "ERROR in range file \"%s\"\n", argv [optind]);
         exit (EXIT_FAILURE);
        }
    fclose (range_file);
   }

   if  (existsFragStore (fragstore_path))
       {
        FragStore  oldfrag_store = openFragStore (fragstore_path, "r");

        old_frag_ct = getLastElemFragStore (oldfrag_store);
        closeFragStore (oldfrag_store);
       }
     else
       {
        fprintf (stderr, "ERROR:  Could not open fragstore \"%s\"\n",
                 fragstore_path);
        exit (EXIT_FAILURE);
       }

   if  (existsFragStore (bactigstore_path))
       {
        FragStore  bactig_store = openFragStore (bactigstore_path, "r");

        last_stored_bactig = getLastElemFragStore (bactig_store);
        closeFragStore (bactig_store);
       }
     else
       {
        fprintf (stderr, "ERROR:  Could not open bactigstore \"%s\"\n",
                 bactigstore_path);
        exit (EXIT_FAILURE);
       }

   if  (last_stored_bactig < last_bactig)
       {
        fprintf (stderr,
                 "ERROR:  bactig range " F_S64 ".." F_S64 " goes beyond last bactig " F_S64 "\n",
                 first_bactig, last_bactig, last_stored_bactig);
        exit (EXIT_FAILURE);
       }

   if  (scriptfile_name == NULL)
       scriptfile_name = DEFAULT_SCRIPTFILE_NAME;


   Emit_Script (scriptfile_name, lsf_queue_name,
                old_frag_ct, first_bactig, last_bactig,
                fragstore_path, bactigstore_path, bol_path,
                option_string);


   printf ("old_frag_ct = %d\n", old_frag_ct);
   printf ("bactig range = " F_S64 " .. " F_S64 "\n", first_bactig, last_bactig);
   printf ("option_string = \"%s\"\n", option_string);
   printf ("fragstore_path = \"%s\"\n", fragstore_path);
   printf ("bactigstore_path = \"%s\"\n", bactigstore_path);
   printf ("scriptfile_name = \"%s\"\n", scriptfile_name);

   fprintf (stderr, "### Return from main\n");

   return  0;
  }



static void  Emit_Script
    (char * scriptfile_name, char * lsf_queue_name,
     int old_ct, int64 first_new, int64 last_new,
     char * fragstore_path, char * bactigstore_path, char * bol_path,
     char * options)

//  Output to file  scriptfile_name  a set of LSF commands
//  to overlap all fragments already in the old fragment store
//  with fragments numbered  first_new .. last_new  in the
//  BACtig store.  There are  old_ct  fragments  in the fragment store.
//  lsf_queue_name  is the queue to which the LSF jobs are submitted.
//  fragstore_path ,  bactigstore_path  and  bol_path  are the locations
//  of the old fragment store, BACtig store and .bol output files,
//  respectively.  Options to overlap commands
//  are in  options .

  {
   FILE  * scriptfile;
   char  id_tag [MAX_NAME_LEN];
   char  job_name [MAX_NAME_LEN];
   char  working_dir [MAX_NAME_LEN];
   char  olap_script_name [MAX_NAME_LEN];
   char  command [MAX_NAME_LEN];
   int64  new_ct;
   int  new_batches, new_incr, tag_ct;
   int  last, i, j;

   tag_ct = 0;
   sprintf (id_tag, "%s" F_PID_T, PROJECT_TAG, getpid ());
   sprintf (working_dir, "OlapDir" F_PID_T, getpid ());
   sprintf (command, "mkdir %s", working_dir);
   system (command);

   scriptfile = File_Open (scriptfile_name, "w");
   fprintf (scriptfile, "#!/usr/bin/csh\n");
   fprintf (scriptfile, "set QUEUE = %s\n", lsf_queue_name);
   fprintf (scriptfile, "set MACHINE_SPEC = \"%s\"\n", MACHINE_SPEC_STRING);

    if  (first_new < 1 || last_new < first_new)
        {
         fprintf (stderr, "No new BACtigs.  Nothing to do\n");
         fprintf (scriptfile, "echo \"No new BACtigs.  Nothing to do\"");
         fclose (scriptfile);

         return;
        }

   new_ct = 1 + last_new - first_new;
   new_batches = ceil ((double) new_ct / NEW_BATCH_SIZE);
   if  (new_batches > MAX_LSF_ARRAY_JOBS)
       new_batches = MAX_LSF_ARRAY_JOBS;
   new_incr = (int) (new_ct / new_batches);

   last = first_new - 1;
   for  (i = 0;  i < new_batches;  i ++)
     {
      int  old_batches, old_incr;
      int  first, lo;

      first = last + 1;
      last += new_incr;

      if  (i >= new_batches - 1)
          last = last_new;

      old_batches = ceil ((double) old_ct / OLD_BATCH_SIZE);
      if  (old_batches > MAX_LSF_ARRAY_JOBS / new_batches)
          old_batches = MAX_LSF_ARRAY_JOBS / new_batches;
      old_incr = (int) (old_ct / old_batches);
      lo = 1;
      for  (j = 0;  j < old_batches;  j ++)
        {
         FILE  * olap_script_fp;
         int  hi;

         if  (j == old_batches - 1)
             hi = old_ct;
           else
             hi = lo + old_incr - 1;
         sprintf (job_name, "%s.%d", id_tag, ++ tag_ct);
         sprintf (olap_script_name, "%s/%s.sh", working_dir, job_name);

         olap_script_fp = File_Open (olap_script_name, "w");
         fprintf (olap_script_fp, "#!/usr/bin/csh\n");
         fprintf (olap_script_fp, "set QUEUE = %s\n", lsf_queue_name);
         fprintf (olap_script_fp, "set MACHINE_SPEC = \"%s\"\n",
                  MACHINE_SPEC_STRING);
         fprintf (olap_script_fp, "%s -h %7d-%-7d -r %7d-%-7d %s"
                  " -b %s/" F_UID ".%d.bol %s %s\n",
                  OLAP_CMD, first, last, lo, hi, options,
                  bol_path, Batch_Msg_UID, tag_ct, fragstore_path,
                  bactigstore_path);
         fprintf (olap_script_fp, "if  (${status} == 0)  then\n");
         fprintf (olap_script_fp, "echo \"olap job %s successful\"\n",
                  job_name);
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
         sprintf (command, "chmod u+x %s", olap_script_name);
         system (command);


         sprintf (olap_script_name, "%s/%sFAIL.sh", working_dir, job_name);
         olap_script_fp = File_Open (olap_script_name, "w");
         fprintf (olap_script_fp, "#!/usr/bin/csh\n");
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
         sprintf (command, "chmod u+x %s", olap_script_name);
         system (command);

         lo = hi + 1;
        }
     }

   if  (tag_ct > 0)
       {
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -q $QUEUE %s -R \"$MACHINE_SPEC\""
                 " -J \"%s/%s[1-%d:1]\" -o \"%s.%%I.log\""
                 " -P \"%s\" -i \"%s/%s.%%I.sh\" csh -s`\n",
                 PROCESSOR_SPEC_OPTION, JOB_GROUP, id_tag, tag_ct, id_tag,
                 LSF_Project_Name, working_dir, id_tag);
        fprintf (scriptfile,
                 "set olapJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Overlap LSF ID = $olapJobID\"\n");
        for  (i = 1;  i <= tag_ct;  i ++)
          fprintf (scriptfile,
                   "bsub -q $QUEUE -P \"%s\" -J \"%s/%s.%dFAIL\""
                   " -w \"exit(${olapJobID}[%d])\""
                   " -o \"%s.%dFAIL.log\" %s/%s.%dFAIL.sh\n",
                   LSF_Project_Name, JOB_GROUP, id_tag, i,
                   i, id_tag, i, working_dir, id_tag, i);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -K -q $QUEUE -P \"%s\" -J \"%s/%sSUCCESS\""
                 " -w \"done(${olapJobID}[1-%d])\""
                 " -o \"%sSUCCESS.log\" \"echo Success; rm -rf %s\"`\n",
                 LSF_Project_Name, JOB_GROUP, id_tag,
                 tag_ct, id_tag, working_dir);
        fprintf (scriptfile,
                 "set successJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Success LSF ID = $successJobID\"\n");
        fprintf (scriptfile, "bjobs  ${successJobID}\n");
        fprintf (scriptfile, "if  (${status} != 0)  then\n");
        fprintf (scriptfile, "exit -1\n");
        fprintf (scriptfile, "else\n");
        fprintf (scriptfile, "bdel ${successJobID}\n");
        fprintf (scriptfile, "endif\n");
       }

            
   fclose (scriptfile);

   return;
  }

