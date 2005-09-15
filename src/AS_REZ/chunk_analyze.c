
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
* Module:  chunk_analyze.c
* Description:
*   Read output from the unitigger and list invalid chunks, i.e.,
*   those with fragments from separate copies of the genome,
*   based on celsim coordinates.  Also output a  chunk_analyze.cam
*   file to show all the genome positions from which such chunks
*   came.
* 
*    Programmer:  A. Delcher
*       Written:  11 Jun 99
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: chunk_analyze.c,v 1.5 2005-09-15 15:20:16 eliv Exp $
 * $Revision: 1.5 $
*/

static char fileID[] = "$Id: chunk_analyze.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";

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


#define  ALLOWED_GAP           100     // gaps <= this are ignored
#define  BUFF_LEN              100
#define  COVERAGE_THRESHOLD    12.0
#define  MAX_LINE              1000


typedef struct
  {
   int  id;
   int  lo, hi;
   int  frags;
  }  Range_t;


#define  NUM_COLOURS   6
#define  UNIQUE_COLOUR         1
#define  CONSISTENT_COLOUR     2
#define  ONEFRAG_COLOUR        3
#define  REPEAT_COLOUR         4
#define  BADUNIQUE_COLOUR      5
char  * Colour_String [NUM_COLOURS]
    = {
       "C000000 T2 S  # Unused",
       "CFFFF00 T2 S  # Unique",
       "C00FF00 T2 S  # Consistent",
       "C0000FF T2 S  # OneFrag",
       "CFF0000 T2 S  # Repeat",
       "CFF9A11 T2 S  # BadUnique"   // Orange
      };
FILE  * Cam_File;


int  Cam_ID = 1;
int  Chunk_Left, Chunk_Right;
int  Num_Frags;
int  Passes_Discriminator;
Range_t  * Range;
int  Range_Ct = 0;
int  Range_Size;
int  Repeat_Chunk_Ct = 0;


void  Add_Range_Entry
    (int fid, int a_end, int b_end);
void  Analyze_Range
    (int chunk_id);
void  Commatize
    (long int n, char buff [], int buff_len);



int main  (int argc, char * argv [])

  {
   FILE  * infile;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   char  buff [BUFF_LEN];
   int  a_end, b_end, chunk_id, fid;
   int  has_r_flag = FALSE;
   int  i, ct;

   if  (argc < 2)
       {
        fprintf (stderr, "USAGE:  %s <cgb-file>\n", argv [0]);
        exit (-1);
       }

   infile = fopen (argv [1], "r");
   if  (infile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file %s\n",
                 argv [1]);
        exit (-1);
       }

   Cam_File = fopen ("chunk_analyze.cam", "w");
   assert (Cam_File != NULL);

   for  (i = 0;  i < NUM_COLOURS;  i ++)
     fprintf (Cam_File, "%dCKA: %s\n", i, Colour_String [i]);


   Range_Size = 10000;
   Range = (Range_t *) malloc (Range_Size * sizeof (Range_t));
   assert (Range != NULL);

   read_msg_fn = (MesgReader)InputFileType_AS (infile);

   while  (read_msg_fn (infile, & gmesg) != EOF)
     {
      if  (gmesg && gmesg -> t == MESG_IUM)
          {
           IntUnitigMesg  * msg;
           char  * p;

           // deal with the generic message as an IUM message
           msg = (IntUnitigMesg *) gmesg -> m;

           chunk_id = msg -> iaccession;
           Passes_Discriminator
               = (msg -> coverage_stat >= COVERAGE_THRESHOLD);
           Num_Frags = msg -> num_frags;
           p = strstr (msg -> source, "gen>");
           if  (p == NULL)
               {
                fprintf (stderr, "ERROR:  Missing  gen>  in source entry\n");
                has_r_flag = FALSE;
                a_end = b_end = -1;
               }
             else
               {
                char  tag [MAX_LINE];

                ct = sscanf (p + 4, "%s [%d,%d]", tag, & a_end, & b_end);
                assert (ct == 3);
                has_r_flag = (tag [0] == 'r');
                if  (a_end <= b_end)
                    {
                     Chunk_Left = a_end;
                     Chunk_Right = b_end;
                    }
                  else
                    {
                     Chunk_Left = b_end;
                     Chunk_Right = a_end;
                    }
               }

           Range_Ct = 0;

           for  (i = 0;  i < Num_Frags;  i ++)
             {
              fid = msg -> f_list [i] . ident;
              p = strstr (msg -> f_list [i] . source, "\n[");
              if  (p != NULL
                     && sscanf (p + 2, "%d,%d", & a_end, & b_end) == 2
                     && has_r_flag)
                  Add_Range_Entry (fid, a_end, b_end);
             }

           Analyze_Range (chunk_id);
          }
     }

   fclose (Cam_File);

   Commatize (Repeat_Chunk_Ct, buff, BUFF_LEN);
   printf ("\nRepeat Chunks = %11s\n", buff);

   return  0;
  }



void  Add_Range_Entry    (int fid, int a_end, int b_end)

//  Append a new entry for  fid  with celsim coordinates
//  [a_end,b_end]  onto global array  Range  and increment
//  Range_Ct .  Expand array and  Range_Size  if necessary.

  {
   if  (Range_Ct >= Range_Size - 1)
       {
        Range_Size *= 2;
        Range = (Range_t *) realloc (Range, Range_Size * sizeof (Range_t));
        assert (Range != NULL);
       }

   Range [Range_Ct] . id = fid;
   if  (a_end <= b_end)
       {
        Range [Range_Ct] . lo = a_end;
        Range [Range_Ct] . hi = b_end;
       }
     else
       {
        Range [Range_Ct] . lo = b_end;
        Range [Range_Ct] . hi = a_end;
       }
   Range_Ct ++;

   return;
  }



void  Analyze_Range
    (int chunk_id)

//  Output information from global  Range  array for the
//  chunk with  chunk_id .

  {
   char  buff1 [100], buff2 [100], buff3 [100], buff4 [100];
   int  i, j, hi, lo;
   int  copy = 0, frags;

   if  (Range_Ct == 0)
       {
        if  (Num_Frags == 1)
            fprintf (Cam_File, "%dCKAID: %d A%dCKA %d R3 # chunk %d\n",
                     Cam_ID ++, Chunk_Left, ONEFRAG_COLOUR, Chunk_Right,
                     chunk_id);
        else if  (Passes_Discriminator)
            fprintf (Cam_File, "%dCKAID: %d A%dCKA %d R1 # chunk %d\n",
                     Cam_ID ++, Chunk_Left, UNIQUE_COLOUR, Chunk_Right,
                     chunk_id);
          else
            fprintf (Cam_File, "%dCKAID: %d A%dCKA %d R3 # chunk %d\n",
                     Cam_ID ++, Chunk_Left, CONSISTENT_COLOUR, Chunk_Right,
                     chunk_id);
        return;
       }

   Repeat_Chunk_Ct ++;

   // Sort according to low celsim coord

   for  (i = 0;  i < Range_Ct - 1;  i ++)
     for  (j = i + 1;  j < Range_Ct;  j ++)
       if  (Range [i] . lo > Range [j] . lo)
           {
            Range_t  save;

            save = Range [i];
            Range [i] = Range [j];
            Range [j] = save;
           }

   printf ("\nChunk #%d:\n", chunk_id);

   // Compress overlapping regions to front of  Range

   lo = Range [0] . lo;
   hi = Range [0] . hi;
   frags = 1;
   copy = 0;

   for  (i = 1;  i < Range_Ct;  i ++)
     if  (Range [i] . lo <= hi + ALLOWED_GAP)   // overlap
         {
          if  (Range [i] . hi > hi)
              hi = Range [i] . hi;
          frags ++;
         }
       else
         {
          Range [copy] . lo = lo;
          Range [copy] . hi = hi;
          Range [copy] . frags = frags;
          copy ++;
          lo = Range [i] . lo;
          hi = Range [i] . hi;
          frags = 1;
         }

   Range [copy] . lo = lo;
   Range [copy] . hi = hi;
   Range [copy] . frags = frags;


   copy ++;

   for  (i = 0;  i < copy;  i ++)
     {
      int  colour;

      Commatize (Range [i] . lo, buff1, 100);
      Commatize (Range [i] . hi, buff2, 100);
      Commatize (Range [i] . hi - Range [i] . lo, buff3, 100);
      printf ("  frags = %3d  %11s .. %11s   len = %9s",
              Range [i] . frags, buff1, buff2, buff3);
      if  (i == 0)
          printf ("\n");
        else
          {
           Commatize (Range [i] . lo - Range [i - 1] . hi, buff4, 100);
           printf ("   gap = %10s\n", buff4);
          }
      if  (Passes_Discriminator)
          colour = BADUNIQUE_COLOUR;
        else
          colour = REPEAT_COLOUR;

      fprintf (Cam_File,
               "%dCKAID: %d A%dCKA %d R1 # chunk %d copy %d of %d  frags = %d (%d others)",
               Cam_ID ++, Range [i] . lo, colour,
               Range [i] . hi, chunk_id, i + 1, copy, Range [i] . frags,
               Range_Ct - Range [i] . frags);
      if  (i > 0)
          {
           Commatize (Range [i - 1] . lo, buff1, 100);
           Commatize (Range [i - 1] . hi, buff2, 100);
           fprintf (Cam_File, "  prev = [%s .. %s]",
                    buff1, buff2);
          }
      if  (i < copy - 1)
          {
           Commatize (Range [i + 1] . lo, buff1, 100);
           Commatize (Range [i + 1] . hi, buff2, 100);
           fprintf (Cam_File, "  next = [%s .. %s]",
                    buff1, buff2);
          }
      fputc ('\n', Cam_File);
     }

   return;
  }



void  Commatize
    (long int n, char buff [], int buff_len)

//  Put printable representation of  n  into  buff  if it fits
//  within length  buff_len , with commas every 3 digits.
//  If it doesn't fit, fill  buff  with stars.

  {
   char  tmp [1000];
   int  i, len;

   len = 0;
   do
     {
      if  (len % 4 == 3)
          tmp [len ++] = ',';
      tmp [len ++] = '0' + n % 10;
      n /= 10;
     } while  (n > 0);

   if  (len >= buff_len)
       {
        memset (buff, '*', buff_len - 1);
        buff [buff_len - 1] = '\0';
        return;
       }

   for  (i = 0;  i < len;  i ++)
     buff [i] = tmp [len - i - 1];
   buff [len] = '\0';

   return;
  }
