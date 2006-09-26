
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
* Module:  CorrectScriptOVL.c
* Description:
*   Create a script of commands (LSF or not by option) to compute
*   fragment corrections and corrected overlaps and save the
*   revised overlap error rates in the overlap store.
* 
*    Programmer:  A. Delcher
*       Written:   6 Mar 2001
*  Last Revised:  
* 
*************************************************/

/* RCS info
 * $Id: CorrectScriptOVL.c,v 1.5 2006-09-26 21:07:45 brianwalenz Exp $
 * $Revision: 1.5 $
*/

static char  CM_ID []
    = "$Id: CorrectScriptOVL.c,v 1.5 2006-09-26 21:07:45 brianwalenz Exp $";


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
#include  "AS_UTL_version.h"
#include  "AS_OVL_overlap.h"
#include  "OlapStoreOVL.h"

# define THE_SHELL "/usr/bin/csh"

#define  AVAIL_SPACE             2.7e9
    //  Number of bytes available for data structures in a batch
    //  of  correct-frags
#define  AVG_FRAG_LEN            700.0
    //  (Over-)Estimate of average clear-range bytes per fragment
#define  CAT_CORRECT_CMD         "cat-corrects"
    //  Command to concatenate correction files from separate fragment
    //  ranges into a single correction file
#define  DEFAULT_BINARY_PATH   "$AS_BIN"
    //  Standard environment variable name for location of binaries
#define  DEFAULT_SCRIPTFILE_NAME  "fragovl.script"
    //  Name of scriptfile produced if not specified on the
    //  command line
#define  ERATE_UPDATE_CMD         "update-erates"
    //  Command to add binary array of revised overlap error rates
    //  to the overlap store
#define  DEFAULT_FRAG_BATCH_SIZE       250000
    //  Number of fragments to process in each batch of  correct-frags
#define  FRAG_CORRECT_CMD         "correct-frags"
    //  Command to read overlaps and determine corrections for a range
    //  of fragments
#define  INITIAL_FILE_LIST_LEN   100
    //  Number of entries in initial list of correction files
#define  JOB_GROUP             "/Assembly"
    //  Prefix used for jobname in LSF commands
#define  MAX_STRING_LEN        5000
    //  Length of longest possible command or name
#define  LSF_QUEUE_NAME        "assembly"
    //  LSF queue to which jobs are submitted
#define  MACHINE_SPEC_STRING   "select[physmem>2400]rusage[physmem=2400]"
    //  LSF option to specify machine characteristics needed to run job
#define  MAX_LSF_ARRAY_JOBS    2000
    //  Most LSF jobs allowed in an LSF array.  Will fail if exceeded.
#define  OLAP_BATCH_SIZE       1000000
    //  Number of fragments to process in each batch of  correct-olaps
#define  OLAP_CORRECT_CMD      "correct-olaps"
    //  Command to read corrections for fragments and re-do overlaps
    //  using the corrected fragments
#define  PROCESSOR_SPEC_OPTION  "-x"
    //  LSF option to specify number of available processors needed to
    //  run job
#define  SPACE_PER_FRAG_IN_FC  13 * AVG_FRAG_LEN
    //  Number of bytes used for each fragment in a batch of
    //  correct-frags.  13 is the 3 bytes per base to store frag data/quality
    //  plus 10 bytes for  Vote_Tally_t .
#define  SPACE_PER_OLAP_IN_FC  12.0
    //  Number of bytes used for each overlap in a batch of
    //  correct-frags
#define  SPACE_PER_FRAG_IN_OC  3 * AVG_FRAG_LEN
    //  Number of bytes used for each fragment in a batch of
    //  correct-olaps.  4 is the 3 bytes per base to store frag data/quality
    //  plus 1 bytes for indel/adjustments
#define  SPACE_PER_OLAP_IN_OC  20.0
    //  Number of bytes used for each overlap in a batch of
    //  correct-olaps
#define  TAG_STRING            "corr"
    //  Tag to use on lsf script file names and jobids


/*************************************************************************/
/* Type definitions */
/*************************************************************************/

typedef  struct
  {
   int32  lo, hi;
  }  Range_t;



/*************************************************************************/
/* Static Globals */
/*************************************************************************/

static int  Adapt_Batch_Sizes = FALSE;
  //  If set true by  -A  option will make batch sizes based on
  //  number of overlaps in the overlap store.  Must be used in
  //  conjunction with the  -S  option.
static char  * Binary_Path = DEFAULT_BINARY_PATH;
static char  * c_Option = NULL;
static char  * Correct_File_Path = NULL;
static char  * d_Option = NULL;
static char  * E_Option = NULL;
static char  * F_Option = NULL;
static int  Frag_Batch_Size = DEFAULT_FRAG_BATCH_SIZE;
    // Number of fragments to correct in each batch
char  * Frag_Store_Path = NULL;
static char  * k_Option = NULL;
static char  * LSF_Project_Name = NULL;
static char  * LSF_Queue_Name = NULL;
static char  * o_Option = NULL;
static int  P_Option = FALSE;
static char  * q_Option = NULL;
static char  * Scriptfile_Name = DEFAULT_SCRIPTFILE_NAME;
static char  * S_Option = NULL;
static int  Use_LSF = TRUE;
static char  * v_Option = NULL;
static char  * x_Option = NULL;
static char  * X_Option = NULL;


/*************************************************************************/
/* External Global Definitions */
/*************************************************************************/



/*************************************************************************/
/* Function prototypes for internal static functions */
/*************************************************************************/

static void  Calculate_Batches
    (int num_frags, int * num_batches, Range_t * * batch_range,
     double space_per_frag, double space_per_olap, int fixed_batch_size);
static void  Emit_Script
    (void);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main  (int argc, char * argv [])

  {
   Parse_Command_Line  (argc, argv);

   Emit_Script ();

   fprintf (stderr, "Finished\n");

   return  0;
  }



static void  Calculate_Batches
    (int num_frags, int * num_batches, Range_t * * batch_range,
     double space_per_frag, double space_per_olap, int fixed_batch_size)

//  Allocate and fill array  (* batch_range)  with the ranges of fragments
//  to do in each batch.  If global  Adapt_Batch_Sizes  is true,
//  use the actual number of overlaps in the overlap store for this.
//  Allow  space_per_frag  bytes for each fragment and  space_per_olap
//  bytes for each overlap.  Otherwise, just make fixed-size batches of
//  size  fixed_batch_size.  Set  (* num_batches)
//  to the number of batches, i.e., the number of entries in the array.
//  num_frags is the total number of fragments that need to be
//  processed.

  {
   OVL_Store_t  * ovl_store;
   Range_t  * bp;
   int32  batch_lo, batch_hi;
   uint32  last_ovl_frag;
   uint32  * offset;
   int  i;

   if  (Adapt_Batch_Sizes)
       {
        uint32  frag_ct, olap_ct, max_olaps, frag_olaps;
        int  batch_ct;

        ovl_store = New_OVL_Store ();
        Open_OVL_Store (ovl_store, S_Option);
        last_ovl_frag = Last_Frag_In_OVL_Store (ovl_store);
        if  (last_ovl_frag != num_frags)
            {
             fprintf (stderr,
                      "WARNING:  ovl store has %u frags  frag store has %d\n",
                      last_ovl_frag, num_frags);
            }

        offset = Load_Frag_Offsets (ovl_store);

        frag_ct = 0;
        olap_ct = 0;
        max_olaps = 0;
        batch_ct = 1;
        (* batch_range) = (Range_t *) safe_malloc
                                  (batch_ct * sizeof (Range_t));
        batch_lo = 1;

        for  (i = 1;  i <= last_ovl_frag;  i ++)
          {
           if  (offset [i + 1] >= offset [i])
               frag_olaps = offset [i + 1] - offset [i];
             else
               frag_olaps = max_olaps;
           if  (frag_olaps > max_olaps)
               max_olaps = frag_olaps;

           if  ((frag_ct + 1) * space_per_frag
                   + (olap_ct + frag_olaps) * space_per_olap
                      <= AVAIL_SPACE)
               {
                frag_ct ++;
                olap_ct += frag_olaps;
               }
             else
               {
                (* batch_range) [batch_ct - 1] . lo = batch_lo;
                (* batch_range) [batch_ct - 1] . hi = i - 1;
fprintf (stderr, "Frag space = %.0f   Olap space = %.0f\n",
         frag_ct * space_per_frag, olap_ct * space_per_olap);
                batch_ct ++;
                (* batch_range) = (Range_t *) safe_realloc
                     ((* batch_range), batch_ct * sizeof (Range_t));
                batch_lo = i;
                frag_ct = 1;
                olap_ct = frag_olaps;
               }
          }
        (* batch_range) [batch_ct - 1] . lo = batch_lo;
        (* batch_range) [batch_ct - 1] . hi = last_ovl_frag;
        (* num_batches) = batch_ct;

for  (i = 0;  i < batch_ct;  i ++)
  fprintf (stderr, "%8d %8d %9d\n", (* batch_range) [i] . lo,
           (* batch_range) [i] . hi,
           1 + (* batch_range) [i] . hi - (* batch_range) [i] . lo);
       }
     else
       {
        (* num_batches) = (int) ceil ((double) num_frags / fixed_batch_size);
        bp = (* batch_range) = (Range_t *) safe_malloc
                                  ((* num_batches) * sizeof (Range_t));

        batch_lo = 1;
        for  (i = 0;  i < (* num_batches);  i ++)
          {
           batch_hi = batch_lo + fixed_batch_size - 1;
           if  (batch_hi > num_frags)
               batch_hi = num_frags;
           bp -> lo = batch_lo;
           bp -> hi = batch_hi;
           batch_lo = batch_hi + 1;
           bp ++;
          }
       }

   return;
  }



static void  Emit_Script
    (void)

//  Output to file  Scriptfile_Name  a set of commands
//  to correct fragments and overlaps.  Fragments are in
//  the fragment store in  Frag_Store_Path .  Overlaps are in
//  the overlap store in  S_option .

  {
   int64  lo, hi, num_frags;
   int  num_fc_batches, num_oc_batches;
   FILE  * listfile, * scriptfile;
   Range_t  * fc_batch_range, * oc_batch_range;
   char  id_tag [MAX_STRING_LEN];
   char  job_name [MAX_STRING_LEN];
   char  listfile_name [MAX_STRING_LEN];
   char  script_name [MAX_STRING_LEN];
   char  working_dir [MAX_STRING_LEN];
   char  command [MAX_STRING_LEN];
   char  sys_command [MAX_STRING_LEN];
   int  i, tag_ct;

   if  (existsFragStore (Frag_Store_Path))
       {
        FragStore  frag_store = openFragStore (Frag_Store_Path, "r");

        num_frags = getLastElemFragStore (frag_store);
        closeFragStore (frag_store);
       }
     else
       {
        fprintf (stderr, "ERROR:  Could not open fragstore \"%s\"\n",
                 Frag_Store_Path);
        exit (EXIT_FAILURE);
       }

   scriptfile = File_Open (Scriptfile_Name, "w");
   fprintf (scriptfile, "#!" THE_SHELL "\n");

   Calculate_Batches (num_frags, & num_fc_batches, & fc_batch_range,
                      SPACE_PER_FRAG_IN_FC, SPACE_PER_OLAP_IN_FC,
                      Frag_Batch_Size);

   sprintf (command, "%s/%s", Binary_Path, FRAG_CORRECT_CMD);
   if  (d_Option != NULL)
       {
        strcat (command, " -d ");
        strcat (command, d_Option);
       }
   if  (F_Option != NULL)
       {
        strcat (command, " -F ");
        strcat (command, F_Option);
       }
   if  (k_Option != NULL)
       {
        strcat (command, " -k ");
        strcat (command, k_Option);
       }
   if  (S_Option != NULL)
       {
        strcat (command, " -S ");
        strcat (command, S_Option);
       }
   if  (v_Option != NULL)
       {
        strcat (command, " -v ");
        strcat (command, v_Option);
       }
   if  (x_Option != NULL)
       {
        strcat (command, " -x ");
        strcat (command, x_Option);
       }

   if  (! Use_LSF)
       {
        if  (num_fc_batches == 1)
            {
             fprintf (scriptfile, "%s -o %s %s 1 " F_S64 "\n",
                      command, Correct_File_Path, Frag_Store_Path, num_frags);
             fprintf (scriptfile, "if  (${status} != 0)  then\n");
             fprintf (scriptfile, "    echo \"%s failed\"\n", FRAG_CORRECT_CMD);
             fprintf (scriptfile, "    exit -1\n");
             fprintf (scriptfile, "endif\n");
            }
          else
            {
             sprintf (listfile_name, "%s" F_PID_T ".cor.list",
                      TAG_STRING, getpid ());
             listfile = File_Open (listfile_name, "w");

             for  (i = 0;  i < num_fc_batches;  i ++)
               {
                fprintf (scriptfile, "%s -o %s.%03d %s %d %d\n",
                         command, Correct_File_Path, i + 1, Frag_Store_Path,
                         fc_batch_range [i] . lo, fc_batch_range [i] . hi);
                fprintf (scriptfile, "if  (${status} != 0)  then\n");
                fprintf (scriptfile, "    echo \"%s %d %d (.%03d) failed\"\n",
                         FRAG_CORRECT_CMD, fc_batch_range [i] . lo,
                         fc_batch_range [i] . hi, i + 1);
                fprintf (scriptfile, "    exit -1\n");
                fprintf (scriptfile, "endif\n");

                fprintf (listfile, "%s.%03d\n", Correct_File_Path, i + 1);
               }
             fclose (listfile);

             fprintf (scriptfile, "%s/%s -L %s -o %s\n",
                      Binary_Path, CAT_CORRECT_CMD, listfile_name, Correct_File_Path);
             fprintf (scriptfile, "if  (${status} != 0)  then\n");
             fprintf (scriptfile, "    echo \"%s failed\"\n", CAT_CORRECT_CMD);
             fprintf (scriptfile, "    exit -1\n");
             fprintf (scriptfile, "endif\n");
            }
       }
     else
       {
        sprintf (id_tag, "%s%s" F_PID_T, TAG_STRING, "frag", getpid ());
        sprintf (working_dir, "CorrectDir" F_PID_T, getpid ());
        sprintf (sys_command, "mkdir %s", working_dir);
        system (sys_command);

        fprintf (scriptfile, "set QUEUE = %s\n", LSF_Queue_Name);
        fprintf (scriptfile, "set MACHINE_SPEC = \"%s\"\n", MACHINE_SPEC_STRING);

        if  (num_fc_batches > MAX_LSF_ARRAY_JOBS)
            {
             fprintf (stderr, "ERROR:  Exceeded LSF max array size\n");
             exit (EXIT_FAILURE);
            }

        sprintf (listfile_name, "%s" F_PID_T ".cor.list", TAG_STRING, getpid ());
        listfile = File_Open (listfile_name, "w");

        tag_ct = 0;
        for  (i = 0;  i < num_fc_batches;  i ++)
          {
           FILE  * script_fp;

           sprintf (job_name, "%s.%d", id_tag, ++ tag_ct);
           sprintf (script_name, "%s/%s.sh", working_dir, job_name);

           script_fp = File_Open (script_name, "w");
           fprintf (script_fp, "#!" THE_SHELL "\n");
           fprintf (script_fp, "ls -l %s/db.frg\n", Frag_Store_Path);
           fprintf (script_fp, "%s -o %s.%03d %s %d %d\n",
                    command, Correct_File_Path, i + 1, Frag_Store_Path,
                    fc_batch_range [i] . lo, fc_batch_range [i] . hi);
           fprintf (script_fp, "if  (${status} == 0)  then\n");
           fprintf (script_fp, "    bdel -J \"%s/%sFAIL\"\n",
                    JOB_GROUP, job_name);
           fprintf (script_fp, "  else\n");
           fprintf (script_fp,
                    "    echo \"%s %d %d (.%03d, job %s) failed"
                    "  status = \" ${status}\n",
                    FRAG_CORRECT_CMD, fc_batch_range [i] . lo,
                    fc_batch_range [i] . hi, i + 1, job_name);
           fprintf (script_fp,
                    "    echo \"This email was sent automatically by"
                    " LSF fragment correction\" >! %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp,
                    "    echo \"%s %d %d (.%03d, job %s) failed  status = \" ${status}"
                    " >> %s/%s.FAILED\n",
                    FRAG_CORRECT_CMD, fc_batch_range [i] . lo,
                    fc_batch_range [i] . hi, i + 1, job_name, working_dir, job_name);
           fprintf (script_fp, "    mail ${user} < %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp, "    bdel -J \"%s/%sSUCCESS\"\n",
                    JOB_GROUP, id_tag);
           fprintf (script_fp, "    bdel -J \"%s/%sFAIL\"\n",
                    JOB_GROUP, job_name);
           fprintf (script_fp, "    exit -1\n");
           fprintf (script_fp, "endif\n");
           fclose (script_fp);
           sprintf (sys_command, "chmod a+x %s", script_name);
           system (sys_command);

           sprintf (script_name, "%s/%sFAIL.sh", working_dir, job_name);
           script_fp = File_Open (script_name, "w");
           fprintf (script_fp, "#!" THE_SHELL "\n");
           fprintf (script_fp,
                    "echo \"%s %d %d (.%03d, job %s) failed  status = \" ${status}\n",
                    FRAG_CORRECT_CMD, fc_batch_range [i] . lo,
                    fc_batch_range [i] . hi, i + 1, job_name);
           fprintf (script_fp,
                    "echo \"This email was sent automatically by LSF fragment"
                    " correction\" >! %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp,
                    "echo \"%s %d %d (.%03d, job %s) failed  status = \" ${status}"
                    " >> %s/%s.FAILED\n",
                    FRAG_CORRECT_CMD, fc_batch_range [i] . lo,
                    fc_batch_range [i] . hi, i + 1, job_name, working_dir, job_name);
           fprintf (script_fp, "    mail ${user} < %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp, "    bdel -J \"%s/%sSUCCESS\"\n",
                    JOB_GROUP, id_tag);
           fclose (script_fp);
           sprintf (sys_command, "chmod a+x %s", script_name);
           system (sys_command);

           fprintf (listfile, "%s.%03d\n", Correct_File_Path, i + 1);
          }
        fclose (listfile);

        assert (tag_ct > 0);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -q $QUEUE %s -R \"$MACHINE_SPEC\""
                 " -E \"ls -l %s/db.frg\""
                 " -J \"%s/%s[1-%d:1]\" -o \"%s.%%I.log\""
                 " -P \"%s\" -i \"%s/%s.%%I.sh\" " THE_SHELL " -s`\n",
                 PROCESSOR_SPEC_OPTION,
                 Frag_Store_Path, JOB_GROUP, id_tag, tag_ct, id_tag,
                 LSF_Project_Name, working_dir, id_tag);
        fprintf (scriptfile,
                 "set CorrectJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"FragCorrect LSF ID = $CorrectJobID\"\n");
        for  (i = 1;  i <= tag_ct;  i ++)
          fprintf (scriptfile,
                   "bsub -q $QUEUE -J \"%s/%s.%dFAIL\" -w"
                   " \"exit(${CorrectJobID}[%d])\""
                   " -o \"%s.%dFAIL.log\" -P \"%s\" %s/%s.%dFAIL.sh\n",
                   JOB_GROUP, id_tag, i,
                   i, id_tag, i, LSF_Project_Name, working_dir, id_tag, i);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -K -q $QUEUE -J \"%s/%sSUCCESS\""
                 " -w \"done(${CorrectJobID}[1-%d])\""
                 " -o \"%sSUCCESS.log\" -P \"%s\" \"echo Success\"`\n",
                 JOB_GROUP, id_tag,
                 tag_ct, id_tag, LSF_Project_Name);
        fprintf (scriptfile,
                 "set successJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Success LSF ID = $successJobID\"\n");

        // Kludge because LSF leaves a pending job
        fprintf (scriptfile,
                 "bdel ${successJobID}\n");

        fprintf (scriptfile, "%s/%s -L %s -o %s\n",
                 Binary_Path, CAT_CORRECT_CMD, listfile_name,
                 Correct_File_Path);
        fprintf (scriptfile, "if  (${status} != 0)  then\n");
        fprintf (scriptfile, "    echo \"%s failed\"\n", CAT_CORRECT_CMD);
        fprintf (scriptfile, "    exit -1\n");
        fprintf (scriptfile, "endif\n");
       }

   // Now do batches for  correct-olaps

   Calculate_Batches (num_frags, & num_oc_batches, & oc_batch_range,
                      SPACE_PER_FRAG_IN_OC, SPACE_PER_OLAP_IN_OC,
                      OLAP_BATCH_SIZE);

   sprintf (command, "%s/%s", Binary_Path, OLAP_CORRECT_CMD);

   if  (c_Option != NULL)
       {
        strcat (command, " -c ");
        strcat (command, c_Option);
       }
   if  (F_Option != NULL)
       {
        strcat (command, " -F ");
        strcat (command, F_Option);
       }
   if  (P_Option)
       {
        strcat (command, " -P");
       }
   if  (q_Option != NULL)
       {
        strcat (command, " -q ");
        strcat (command, q_Option);
       }
   if  (S_Option != NULL)
       {
        strcat (command, " -S ");
        strcat (command, S_Option);
       }
   if  (v_Option != NULL)
       {
        strcat (command, " -v ");
        strcat (command, v_Option);
       }
   if  (X_Option != NULL)
       {
        strcat (command, " -X ");
        strcat (command, X_Option);
       }

   sprintf (id_tag, "%s%s" F_PID_T, TAG_STRING, "olap", getpid ());
   if  (num_oc_batches == 1)
       {
        if  (o_Option != NULL)
            {
             strcat (command, " -o ");
             strcat (command, o_Option);
            }
        if  (S_Option == NULL)
            fprintf (scriptfile, "%s %s %s 1 " F_S64 "\n",
                     command, Frag_Store_Path, Correct_File_Path, num_frags);
          else
            fprintf (scriptfile, "%s -e %s.erate %s %s 1 " F_S64 "\n",
                     command, id_tag, Frag_Store_Path, Correct_File_Path, num_frags);
        fprintf (scriptfile, "if  (${status} != 0)  then\n");
        fprintf (scriptfile, "    echo \"%s failed\"\n", OLAP_CORRECT_CMD);
        fprintf (scriptfile, "    exit -1\n");
        fprintf (scriptfile, "endif\n");
        if  (S_Option != NULL)
            {
             fprintf (scriptfile, "%s/%s %s %s.erate\n",
                      Binary_Path, ERATE_UPDATE_CMD, S_Option, id_tag);
             fprintf (scriptfile, "if  (${status} != 0)  then\n");
             fprintf (scriptfile, "    echo \"%s %s %s.erate failed\"\n",
                      ERATE_UPDATE_CMD, S_Option, id_tag);
             fprintf (scriptfile, "    exit -1\n");
             fprintf (scriptfile, "endif\n");
            }
       }
   else if  (! Use_LSF || S_Option == NULL)
       {
        for  (i = 0;  i < num_oc_batches;  i ++)
          {
           char  option [MAX_STRING_LEN] = "";

           if  (o_Option != NULL)
               sprintf (option, "-o %s.%03d ", o_Option, i + 1);
           lo = oc_batch_range [i] . lo;
           hi = oc_batch_range [i] . hi;
           assert  (hi <= num_frags);
           if  (S_Option == NULL)
               fprintf (scriptfile, "%s %s%s %s " F_S64 " " F_S64 "\n",
                        command, option, Frag_Store_Path, Correct_File_Path, lo, hi);
             else
               fprintf (scriptfile, "%s -e %s.erate %s%s %s " F_S64 " " F_S64 "\n",
                        command, id_tag, option, Frag_Store_Path, Correct_File_Path,
                        lo, hi);
           fprintf (scriptfile, "if  (${status} != 0)  then\n");
           fprintf (scriptfile, "    echo \"%s " F_S64 " " F_S64 " failed\"\n",
                    OLAP_CORRECT_CMD, lo, hi);
           fprintf (scriptfile, "    exit -1\n");
           fprintf (scriptfile, "endif\n");
           if  (S_Option != NULL)
               {
                fprintf (scriptfile, "%s/%s %s %s.erate\n",
                         Binary_Path, ERATE_UPDATE_CMD, S_Option, id_tag);
                fprintf (scriptfile, "if  (${status} != 0)  then\n");
                fprintf (scriptfile, "    echo \"%s %s %s.erate failed\"\n",
                         ERATE_UPDATE_CMD, S_Option, id_tag);
                fprintf (scriptfile, "    exit -1\n");
                fprintf (scriptfile, "endif\n");
               }
          }
       }
     else
       {
        if  (num_oc_batches > MAX_LSF_ARRAY_JOBS)
            {
             fprintf (stderr, "ERROR:  Exceeded LSF max array size\n");
             exit (EXIT_FAILURE);
            }

        tag_ct = 0;
        for  (i = 0;  i < num_oc_batches;  i ++)
          {
           FILE  * script_fp;
           char  option [MAX_STRING_LEN] = "";

           if  (o_Option != NULL)
               sprintf (option, "-o %s.%03d ", o_Option, i + 1);
           lo = oc_batch_range [i] . lo;
           hi = oc_batch_range [i] . hi;
           assert  (hi <= num_frags);

           sprintf (job_name, "%s.%d", id_tag, ++ tag_ct);
           sprintf (script_name, "%s/%s.sh", working_dir, job_name);

           script_fp = File_Open (script_name, "w");
           fprintf (script_fp, "#!" THE_SHELL "\n");
           fprintf (script_fp, "ls -l %s/db.frg\n", Frag_Store_Path);
           fprintf (script_fp, "%s -e %s.erate.%03d %s%s %s " F_S64 " " F_S64 "\n",
                    command, id_tag, i + 1, option, Frag_Store_Path,
                    Correct_File_Path, lo, hi);
           fprintf (script_fp, "if  (${status} == 0)  then\n");
           fprintf (script_fp, "    bdel -J \"%s/%sFAIL\"\n",
                    JOB_GROUP, job_name);
           fprintf (script_fp, "  else\n");
           fprintf (script_fp,
                    "    echo \"%s " F_S64 " " F_S64 " (.%03d, job %s) failed"
                    "  status = \" ${status}\n",
                    OLAP_CORRECT_CMD, lo, hi, i + 1, job_name);
           fprintf (script_fp,
                    "    echo \"This email was sent automatically by"
                    " LSF fragment correction\" >! %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp,
                    "    echo \"%s " F_S64 " " F_S64 " (.%03d, job %s) failed  status = \" ${status}"
                    " >> %s/%s.FAILED\n",
                    OLAP_CORRECT_CMD, lo, hi, i + 1, job_name, working_dir, job_name);
           fprintf (script_fp, "    mail ${user} < %s/%s.FAILED\n",
                    working_dir, job_name);
           fprintf (script_fp, "    bdel -J \"%s/%sSUCCESS\"\n",
                    JOB_GROUP, id_tag);
           fprintf (script_fp, "    bdel -J \"%s/%sFAIL\"\n",
                    JOB_GROUP, job_name);
           fprintf (script_fp, "    exit -1\n");
           fprintf (script_fp, "endif\n");
           fclose (script_fp);
           sprintf (sys_command, "chmod a+x %s", script_name);
           system (sys_command);
          }

        assert (tag_ct > 0);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -q $QUEUE %s -R \"$MACHINE_SPEC\""
                 " -E \"ls -l %s/db.frg\""
                 " -J \"%s/%s[1-%d:1]\" -o \"%s.%%I.log\""
                 " -P \"%s\" -i \"%s/%s.%%I.sh\" " THE_SHELL " -s`\n",
                 PROCESSOR_SPEC_OPTION,
                 Frag_Store_Path, JOB_GROUP, id_tag, tag_ct, id_tag,
                 LSF_Project_Name, working_dir, id_tag);
        fprintf (scriptfile,
                 "set CorrectJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"OlapCorrect LSF ID = $CorrectJobID\"\n");
        for  (i = 1;  i <= tag_ct;  i ++)
          fprintf (scriptfile,
                   "bsub -q $QUEUE -J \"%s/%s.%dFAIL\" -w"
                   " \"exit(${CorrectJobID}[%d])\""
                   " -o \"%s.%dFAIL.log\" -P \"%s\" %s/%s.%dFAIL.sh\n",
                   JOB_GROUP, id_tag, i,
                   i, id_tag, i, LSF_Project_Name, working_dir, id_tag, i);
        fprintf (scriptfile,
                 "set bsubReturn = `bsub -K -q $QUEUE -J \"%s/%sSUCCESS\""
                 " -w \"done(${CorrectJobID}[1-%d])\""
                 " -o \"%sSUCCESS.log\" -P \"%s\" \"echo Success\"`\n",
                 JOB_GROUP, id_tag,
                 tag_ct, id_tag, LSF_Project_Name);
        fprintf (scriptfile,
                 "set successJobID = `echo $bsubReturn | cut -f 1 -d \\>"
                 " | cut -f 2 -d \\<`\n");
        fprintf (scriptfile,
                 "echo \"Success LSF ID = $successJobID\"\n");

        // Kludge because LSF leaves a pending job
        fprintf (scriptfile,
                 "bdel ${successJobID}\n");

        for  (i = 1;  i <= num_oc_batches;  i ++)
          {
           fprintf (scriptfile, "%s %s %s.erate.%03d\n",
                    ERATE_UPDATE_CMD, S_Option, id_tag, i);
           fprintf (scriptfile, "if  (${status} != 0)  then\n");
           fprintf (scriptfile, "    echo \"%s %s %s.erate.%03d failed\"\n",
                    ERATE_UPDATE_CMD, S_Option, id_tag, i);
           fprintf (scriptfile, "    exit -1\n");
           fprintf (scriptfile, "endif\n");
          }
       }

   fclose (scriptfile);

   sprintf (sys_command, "chmod a+x %s", Scriptfile_Name);
   system (sys_command);

   return;
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
             && ((ch = getopt (argc, argv,
                          "AbB:c:d:E:F:g:hk:L:o:Pq:Q:s:S:v:x:X:")) != EOF))
     switch  (ch)
       {
        case  'A' :
          Adapt_Batch_Sizes = TRUE;
          break;

        case  'b' :
          Use_LSF = FALSE;
          break;

        case  'B' :
          Binary_Path = optarg;
          break;

        case  'c' :
          c_Option = optarg;
          break;

        case  'd' :
          d_Option = optarg;
          break;

        case  'E' :
          E_Option = optarg;
          break;

        case  'F' :
          F_Option = optarg;
          break;

        case  'g' :
          Frag_Batch_Size = strtol (optarg, & p, 10);
          if  (p == optarg || p == NULL || Frag_Batch_Size < 1)
              {
               fprintf (stderr, "ERROR:  Bad Fragment Batch Size = \"%s\"\n",
                        optarg);
               exit (-1);
              }
          fprintf (stderr, "Fragment Batch Size set to %d\n", Frag_Batch_Size);
          break;

        case  'h' :
          Usage (argv [0]);
          exit (EXIT_SUCCESS);

        case  'k' :
          k_Option = optarg;
          break;

        case  'L' :
          LSF_Project_Name = optarg;
          break;

        case  'o' :
          o_Option = optarg;
          break;

        case  'P' :
          P_Option = TRUE;
          break;

        case  'q' :
          q_Option = optarg;
          break;

        case  'Q' :
          LSF_Queue_Name = optarg;
          break;

        case  's' :
          Scriptfile_Name = optarg;
          break;

        case  'S' :
          S_Option = optarg;
          break;

        case  'v' :
          v_Option = optarg;
          break;

        case  'x' :
          x_Option = optarg;
          break;

        case  'X' :
          X_Option = optarg;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc - 2)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   if  (F_Option == NULL && S_Option == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify overlaps with -F or -S\n");
        exit (EXIT_FAILURE);
       }

   if  (Use_LSF && LSF_Project_Name == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify LSF project name with -L\n");
        exit (EXIT_FAILURE);
       }

   if  (Use_LSF && LSF_Queue_Name == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify LSF queue name with -Q\n");
        exit (EXIT_FAILURE);
       }

   if  (Adapt_Batch_Sizes && S_Option == NULL)
       {
        fprintf (stderr, "ERROR:  Must specify overlap store with -S\n");
        exit (EXIT_FAILURE);
       }

   Frag_Store_Path = argv [optind ++];

   Correct_File_Path = argv [optind ++];

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
       "USAGE:  %s [-AbhP] [-B <bin-dir>] [-c <cgb-file>] [-d <degr-thresh>]\n"
       "           [-g <fragbatchsize>] [-k <kmer-len>] [-o <ovl_file>]\n"
       "           [-q <qual_val>] [-s <scriptname>] [-v <num>] [-x <exclude-end>]\n"
       "           [-F OlapFile | -S OlapStore]\n"
       "           -L <proj-name> -Q <queue>\n"
       "           <FragStore> <CorrectFile>\n"
       "\n"
       "Generates script to correct all fragments in  <FragStore>\n"
       "and recalculate their overlaps, outputting the ones better than\n"
       "the quality threshold to a new <ovl_file>\n"
       "\n"
       "Options:\n"
       "-A         adapt batch sizes to number of overlaps in the overlap\n"
       "           store.  Must use -S option also\n"
       "-b         create a simple sequential batch file that doesn't\n"
       "           use LSF\n"
       "-B <path>  specify path to executables, default is  $AS_BIN \n"
       "-c <file>  specify CGB file which is used to determine a\n"
       "           list of overlaps that do not break existing unitigs.\n"
       "           Files  ium.id  and  unitig.pair  are generated listing\n"
       "           the unitig of each fragment and pairs of unitigs that\n"
       "           have a corrected overlap between them.  Only use if\n"
       "           generating a single batch.\n"
       "-d <num>   specify degree threshold for  correct-frags  program\n"
       "-F <file>  specify file of sorted overlaps to use (in the format\n"
       "           produced by  get-olaps\n"
       "-g <num>   specify number of fragments to correct in each batch\n"
       "-h         output this message\n"
       "-k <num>   specify minimum exact match region to confirm a fragment\n"
       "           in  correct-frags  program\n"
       "-L <name>  specify LSF project name, REQUIRED\n"
       "-o <file>  specifies name of file to which OVL messages go\n"
       "-P         make output .ovl file ASCII (proto) format\n"
       "-q <num>   overlaps less than this error rate are\n"
       "           automatically output\n"
       "-Q <name>  specify LSF queue, REQUIRED\n"
       "-s <name>  specify output scriptfile name, default is %s\n"
       "-S <store> specify the binary overlap store containing overlaps to use\n"
       "-v <num>   specify level of verbose outputs, higher is more\n"
       "-x <num>   specify number of bases on end of exact match to exclude\n"
       "           when confirming in  correct-frags  program\n",
       command, DEFAULT_SCRIPTFILE_NAME);

   return;
  }

