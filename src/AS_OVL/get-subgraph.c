
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
* Module:  get-subgraph.c
* Description:
*   Reads the stream of messages produced by the overlapper
*   and filters those that do not refer to fragments contained
*   int the list of internal fragment ID's specified on the
*   the command line.
*
*    Programmer:  A. Delcher
*       Written:  15 Jul 99
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: get-subgraph.c,v 1.4 2005-03-22 19:49:19 jason_miller Exp $
 * $Revision: 1.4 $
*/

static char fileID[] = "$Id: get-subgraph.c,v 1.4 2005-03-22 19:49:19 jason_miller Exp $";

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
#include  "AS_OVL_delcher.h"


#define  HASH_LOAD_FACTOR     0.3333
    //  Portion of hash table to be occupied by keys.
#define  SKIP_MODULUS         11
    //  Modulus for secondary hash function to determine skip distance
    //  for linear probing.


static int  Hash_Find
    (int32 key, int32 tab [], int32 size);
static int32  Next_Odd_Prime
    (int32 N);
static void  Read_Frag_IDs
    (FILE * fidfile, int32 * * hash_table, int32 * hash_table_size);



int main  (int argc, char * argv [])

  {
   FILE  * ovlfile, * outfile, * fidfile;
   char  * infile_name, * outfile_name;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   MesgWriter  write_msg_fn;
   GenericMesg  * pmesg;
   AuditLine  audit_line;
   AuditMesg  * new_adt_mesg;
   int32  * hash_table, hash_table_size;
   char  label_line [1000];
   int  use_binary_output = TRUE;
   int  ch, error_flag, len, i, j;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "P")) != EOF))
     switch  (ch)
       {
        case  'P':
          use_binary_output = FALSE;
          break;
        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind >= argc)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-P] <ovl-file> <frag-id-file>\n", 
                 argv [0]);
        exit (EXIT_FAILURE);
       }

   infile_name = strdup (argv [optind ++]);
   assert (infile_name != NULL);
   fprintf (stderr, "Input Overlap File = %s\n", infile_name);
   len = strlen (infile_name);
   outfile_name = (char *) Safe_malloc (len + 20);
   for  (i = len - 1;  i >= 0;  i --)
     if  (infile_name [i] == '.')
         break;
   for  (j = 0;  j <= i;  j ++)
     outfile_name [j] = infile_name [j];
   strcpy (outfile_name + j, "sub.");
   strcat (outfile_name, infile_name + i + 1);
   fprintf (stderr, "Output Overlap File = %s\n", outfile_name);

   ovlfile = File_Open (infile_name, "r");
   outfile = File_Open (outfile_name, "w");
   fidfile = File_Open (argv [optind], "r");

   Read_Frag_IDs (fidfile, & hash_table, & hash_table_size);
fprintf (stderr, "Hash table size = %d\n", hash_table_size);

   read_msg_fn = InputFileType_AS (ovlfile);
   if  (use_binary_output)
       write_msg_fn = OutputFileType_AS (AS_BINARY_OUTPUT);
     else
       write_msg_fn = OutputFileType_AS (AS_PROTO_OUTPUT);

   pmesg = (GenericMesg *) Safe_malloc (sizeof (GenericMesg));
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) Safe_malloc (sizeof (AuditMesg));
   new_adt_mesg = pmesg -> m;
   new_adt_mesg -> list = & audit_line;
      

   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_ADT :
          {
           AuditMesg  * adt_mesg = gmesg -> m;

           sprintf (label_line, "%s %s %s", argv [0], infile_name,
                    argv [optind]);
           AppendAuditLine_AS (adt_mesg, & audit_line, time (0), "get-subgraph",
                               "$Revision: 1.4 $", label_line);
           write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_ILK :
          {
           InternalLinkMesg  * ilk_mesg = gmesg -> m;
          
           if  (Hash_Find (ilk_mesg -> ifrag1, hash_table, hash_table_size)
                  && Hash_Find (ilk_mesg -> ifrag2, hash_table, hash_table_size))
               write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_OFG :
          {
           OFGMesg  * ofg_mesg = gmesg -> m;
          
           if  (Hash_Find (ofg_mesg -> iaccession, hash_table,
                           hash_table_size))
               write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_OVL :
          {
           OverlapMesg  * ovl_mesg = gmesg -> m;

           if  (Hash_Find (ovl_mesg -> aifrag, hash_table, hash_table_size)
                  && Hash_Find (ovl_mesg -> aifrag, hash_table, hash_table_size))
               write_msg_fn (outfile, gmesg);
           break;
          }

        default :
          write_msg_fn (outfile, gmesg);
       }

   fclose (ovlfile);
   fclose (outfile);

   return  0;
  }



static int  Hash_Find
    (int32 key, int32 tab [], int32 size)

//  Return  TRUE  iff  key  occurs in hash table  tab [0 .. (size - 1)] .

  {
   int  skip, sub;

   sub = key % size;
   skip = 1 + (key % SKIP_MODULUS);

   while  (tab [sub] != -1 && tab [sub] != key)
     sub = (sub + skip) % size;

   return  (tab [sub] == key);
  }



static int32  Next_Odd_Prime
    (int32 N)

//  Return the first odd prime  >= N .  Return  0  if can't find
//  one.

  {
   int32  Div, Last;

   if  (N % 2 == 0)
       N ++;
   while  (N < INT_MAX)
     {
      Last = (int64) (sqrt ((double) N));
      for  (Div = 3;  Div <= Last;  Div += 2)
        if  (N % Div == 0)
            break;
      if  (Div > Last)
          return  N;
      N += 2;
     }

   return  0;
  }



static void  Read_Frag_IDs
    (FILE * fp, int32 * * tab, int32 * size)

//  Read integers from  fp  and insert them into a hash table  (* tab) .
//  Set  (* size)  to the size of the hash table.  The hash function is
//  just  key % size .

  {
   int32  i, ct, key, skip, sub;

   for  (ct = 0;  fscanf (fp, "%d", & key) != EOF;  ct ++)
     ;

fprintf (stderr, "ct = %d\n", ct);
   (* size) = Next_Odd_Prime ((int32) (ct / HASH_LOAD_FACTOR));
   assert (* size > 0);
   (* tab) = (int32 *) Safe_malloc ((* size) * sizeof (int32));

fprintf (stderr, "(* size) = %d\n", (* size));
   for  (i = 0;  i < (* size);  i ++)
     (* tab) [i] = -1;

   rewind (fp);

   while  (fscanf (fp, "%d", & key) != EOF)
     {
      sub = key % (* size);
      skip = 1 + (key % SKIP_MODULUS);

      while  ((* tab) [sub] != -1 && (* tab) [sub] != key)
{
fprintf (stderr, "try sub = %d\n", sub);
        sub = (sub + skip) % (* size);
}

fprintf (stderr, "key = %d  sub = %d\n", key, sub);
      if  ((* tab) [sub] == -1)
          (* tab) [sub] = key;
     }

   return;
  }
