
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
* Module:  batchify.c
* Description:
*   Reads a stream of fragment messages (i.e., a  .frg ,  .inp  or
*   .urc  file) and splits it into separate batches, each (except the
*   last) containing the specified number of fragments.
*
*    Programmer:  A. Delcher
*       Written:  1 Jul 99
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: batchify.c,v 1.4 2005-03-22 19:49:26 jason_miller Exp $
 * $Revision: 1.4 $
*/

static char fileID[] = "$Id: batchify.c,v 1.4 2005-03-22 19:49:26 jason_miller Exp $";

#define FILE_NAME_FORMAT "b%03d."

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"


#define  DEFAULT_BATCH_SIZE         1000000


int main  (int argc, char * argv [])

  {
   FILE  * infile, * outfile;
   char  * infile_name, * outfile_name;
   char  * p, * file_suffix, * file_tail;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   MesgWriter  write_msg_fn;
   GenericMesg  * pmesg;
   AuditLine  audit_line;
   AuditMesg  * adt_mesg, * new_adt_mesg;
   char  batch_line [100];
   int64  batch_size = DEFAULT_BATCH_SIZE, frag_ct;
   int  use_binary_output = TRUE;
   int  batch_num;
   int  ch, error_flag, len, i, j;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "Pb:")) != EOF))
     switch  (ch)
       {
        case  'P':
          use_binary_output = FALSE;
          break;
        case  'b':
          batch_size = STR_TO_INT64(optarg, & p, 10);
          if  (batch_size <= 0 || p == NULL)
              {
               fprintf (stderr, "ERROR:  Illegal batch size = \"%s\"\n",
                        optarg);
               exit (EXIT_FAILURE);
              }
          break;
        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind >= argc)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-P] [-b <batch size>] <dataset_name>\n", 
                 argv [0]);
        exit (EXIT_FAILURE);
       }

   infile_name = strdup (argv [optind ++]);
   assert (infile_name != NULL);
   fprintf (stderr, "Input File_Name = %s\n", infile_name);
   len = strlen (infile_name);
   outfile_name = (char *) malloc (len + 20);
   for  (i = len - 1;  i >= 0;  i --)
     if  (infile_name [i] == '.')
         break;
   for  (j = 0;  j <= i;  j ++)
     outfile_name [j] = infile_name [j];
   file_suffix = infile_name + i + 1;
   file_tail = outfile_name + j;

   infile = fopen (infile_name, "r");
   if  (infile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file %s\n",
                 infile_name);
        exit (EXIT_FAILURE);
       }
   read_msg_fn = InputFileType_AS (infile);

   if  (use_binary_output)
       write_msg_fn = OutputFileType_AS (AS_BINARY_OUTPUT);
     else
       write_msg_fn = OutputFileType_AS (AS_PROTO_OUTPUT);

   frag_ct = batch_num = 0;
   sprintf (file_tail, FILE_NAME_FORMAT, batch_num);
   strcat (file_tail, file_suffix);
   outfile = fopen (outfile_name, "w");
   if  (outfile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file %s\n",
                 outfile_name);
        exit (EXIT_FAILURE);
       }

   pmesg = (GenericMesg *) malloc (sizeof (GenericMesg));
   assert (pmesg != NULL);
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) malloc (sizeof (AuditMesg));
   assert (pmesg -> m != NULL);
   new_adt_mesg = (AuditMesg *) pmesg -> m;
   new_adt_mesg -> list = & audit_line;
      

   while  (read_msg_fn (infile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_ADT :
          adt_mesg = (AuditMesg *) gmesg -> m;
          sprintf (batch_line, "Batch %d", batch_num);
          AppendAuditLine_AS (adt_mesg, & audit_line, time (0), "batchify",
                              "$Revision: 1.4 $", batch_line);
          write_msg_fn (outfile, gmesg);
          break;

        case  MESG_FRG :
        case  MESG_IFG :
        case  MESG_SFG :
          frag_ct ++;
          if  (frag_ct > batch_size)
              {

               fclose (outfile);
               batch_num ++;
               frag_ct = 1;

               sprintf (file_tail, FILE_NAME_FORMAT, batch_num);
               strcat (file_tail, file_suffix);
               outfile = fopen (outfile_name, "w");
               if  (outfile == NULL)
                   {
                    fprintf (stderr, "ERROR:  Could not open file %s\n",
                             outfile_name);
                    exit (EXIT_FAILURE);
                   }

               audit_line . next = NULL;
               audit_line . name = "batchify";
               audit_line . complete = time (0);
               audit_line . version = "$Revision: 1.4 $";
               sprintf (batch_line, "Batch %d", batch_num);
               audit_line . comment = batch_line;
               write_msg_fn (outfile, pmesg);
              }
          // fall through

        default :
          write_msg_fn (outfile, gmesg);
       }

   fclose (infile);
   fclose (outfile);

   return  0;
  }



