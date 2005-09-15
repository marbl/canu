
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
* Module:  screen_analyze.c
* Description:
*   Reads a stream of screened fragment messages (i.e., a
*   .urc  file) and outputs a list of which fragments are screened
*   by which screen items, and summary statistics on the number of
*   fragments and bases screened by each screen element and each
*   screen-element class.
*   Output goes to the standard output.
*
*    Programmer:  A. Delcher
*       Written:  6 Jul 99
*  Last Revised:  26 Aug 99  Convert to proto I/O instead of ASCII
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: screen_analyze.c,v 1.5 2005-09-15 15:20:16 eliv Exp $
 * $Revision: 1.5 $
*/

static char fileID[] = "$Id: screen_analyze.c,v 1.5 2005-09-15 15:20:16 eliv Exp $";

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
#include <ctype.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"


#define  BUFF_LEN              1000
#define  DEFAULT_LONG_LISTING  FALSE
#define  HARD_RELEVANCE_MASK   AS_OVL_HEED_RPT
#define  MAX_SCREEN_ITEMS      10000
#define  REPORT_INTERVAL       10000


typedef struct
  {
   int  screen_id;
   int  screen_class;
   int  lo, hi;
  }  Screen_Item_t;

typedef struct
  {
   int  id;
   int  hits;
   int  bases;
  }  Stat_Item_t;


static int  Class_Entries = 0;
static int  Ignore_Relevance = FALSE;
static int  Item_Entries = 0;
static Stat_Item_t  Item_Stats [MAX_SCREEN_ITEMS], Class_Stats [MAX_SCREEN_ITEMS];
static int  Long_Listing = DEFAULT_LONG_LISTING;
static unsigned  Relevance_Mask = HARD_RELEVANCE_MASK;
static int  Screened_Frag_Ct = 0;
static int  Total_Frags = 0;
static int  Total_Screened_Bases = 0;


static void  Add_Entry
    (Stat_Item_t * stat, int * ct, int id, int bases);
static void  Analyze_Screen_Items
    (int fid, Screen_Item_t * screen, int screen_ct);
static void  Commatize
    (long int n, char buff [], int buff_len);



int main  (int argc, char * argv [])

  {
   FILE  * urcfile;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   char  buff [BUFF_LEN];
   long int  clear_sum = 0;
   int  i, ch, error_flag;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "al")) != EOF))
     switch  (ch)
       {
        case  'a':
          // List all screen matches regardless of relevance

          Ignore_Relevance = TRUE;
          break;

        case  'l':
          // List screen matches of each fragment

          Long_Listing = TRUE;
          break;

        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind > argc - 1)
       {
        fprintf (stderr, 
                 "USAGE:  %s <urc-file>\n", 
                 argv [0]);
        fprintf (stderr, "Options:\n");
        fprintf (stderr, " -a  Use *all* matches regardless of relevance\n");
        fprintf (stderr, " -l  Long list of screen items for each fragment\n");
        exit (EXIT_FAILURE);
       }

   urcfile = fopen (argv [optind], "r");
   if  (urcfile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file %s\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }
   read_msg_fn = (MesgReader)InputFileType_AS (urcfile);

   if  (! Long_Listing)
       fprintf (stderr, "\n");
   while  (read_msg_fn (urcfile, & gmesg) != EOF && gmesg != NULL)
     {
      if  (! Long_Listing && Total_Frags % REPORT_INTERVAL == 0)
          {
           char  buff [BUFF_LEN];

           Commatize (Total_Frags, buff, BUFF_LEN);
           fprintf (stderr, "\rFrags read = %s", buff);
          }

      if  (gmesg -> t == MESG_SFG)
          {
           IntScreenMatch  * p;
           ScreenedFragMesg  * sfg_mesg;
           Screen_Item_t  screen [MAX_SCREEN_ITEMS];
           int  screen_ct = 0;
           int  frag_id, frag_len;

           Total_Frags ++;

           sfg_mesg = (ScreenedFragMesg*) gmesg -> m;
           frag_id = sfg_mesg -> iaccession;
           frag_len = sfg_mesg -> clear_rng . end - sfg_mesg -> clear_rng . bgn;
           clear_sum += frag_len;

           for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
             {
              if  ((p -> relevance & Relevance_Mask) || Ignore_Relevance)
                  {
                   screen [screen_ct] . screen_id = p -> iwhat;
                   screen [screen_ct] . screen_class = p -> repeat_id;
                   screen [screen_ct] . lo = p -> where . bgn;
                   screen [screen_ct] . hi = p -> where . end;
                   screen_ct ++;
                  }
             }

           Analyze_Screen_Items (frag_id, screen, screen_ct);
          }
     }

   if  (! Long_Listing)
       {
        Commatize (Total_Frags, buff, BUFF_LEN);
        fprintf (stderr, "\rFrags read = %s\n", buff);
       }

   if  (Total_Frags == 0)
       {
        fprintf (stderr, "\a\a\aERROR:  0 fragments read\n");
        exit (EXIT_FAILURE);
       }
   if  (clear_sum == 0)
       {
        fprintf (stderr, "\a\a\aERROR:  0 clear bases read\n");
        exit (EXIT_FAILURE);
       }

   printf ("\n%8s %15s %5s %16s %5s\n",
           "Screen", "Frags", "", "Total", "");
   printf ("%8s %15s %5s %16s %5s\n",
           "Item", "Hit", "", "Bases", "");
   for  (i = 0;  i < Item_Entries;  i ++)
     {
      char  buff1 [BUFF_LEN], buff2 [BUFF_LEN];

      Commatize (Item_Stats [i] . hits, buff1, BUFF_LEN);
      Commatize (Item_Stats [i] . bases, buff2, BUFF_LEN);
      printf ("%8d %12s (%5.2f%%) %13s (%5.2f%%)\n",
              Item_Stats [i] . id,
              buff1, (100.0 * Item_Stats [i] . hits) / Total_Frags,
              buff2, (100.0 * Item_Stats [i] . bases) / clear_sum);
     }

   printf ("\n%8s %15s %5s %16s %5s\n",
           "Screen", "Frags", "", "Total", "");
   printf ("%8s %15s %5s %16s %5s\n",
           "Class", "Hit", "", "Bases", "");
   for  (i = 0;  i < Class_Entries;  i ++)
     {
      char  buff1 [BUFF_LEN], buff2 [BUFF_LEN];

      Commatize (Class_Stats [i] . hits, buff1, BUFF_LEN);
      Commatize (Class_Stats [i] . bases, buff2, BUFF_LEN);
      printf ("%8d %12s (%5.2f%%) %13s (%5.2f%%)\n",
              Class_Stats [i] . id,
              buff1, (100.0 * Class_Stats [i] . hits) / Total_Frags,
              buff2, (100.0 * Class_Stats [i] . bases) / clear_sum);
     }

   Commatize (Screened_Frag_Ct, buff, BUFF_LEN);
   printf ("\nTotal Screened Frags = %11s (%4.1f%%)\n",
           buff, (100.0 * Screened_Frag_Ct) / Total_Frags);
   Commatize (Total_Screened_Bases, buff, BUFF_LEN);
   printf ("Total Screened Bases = %11s (%4.1f%%)\n",
           buff, (100.0 * Total_Screened_Bases) / clear_sum);
   Commatize (clear_sum, buff, BUFF_LEN);
   printf ("   Total Clear Bases = %11s\n", buff);
   Commatize (Total_Frags, buff, BUFF_LEN);
   printf ("     Total Fragments = %11s\n", buff);

   return  0;
  }



static void  Add_Entry
    (Stat_Item_t * stat, int * ct, int id, int bases)

//  Increment entry for  id  in  stat [0 .. (ct - 1)]  by
//  adding  bases  to its field and adding  1  to  hits  field.
//  If necessary create a new entry and increment  (* ct) .

  {
   int  i, j;

   for  (i = (* ct) - 1;  i >= 0 && stat [i] . id >= id;  i --)
     if  (stat [i] . id == id)
         {
          stat [i] . hits ++;
          stat [i] . bases += bases;
          return;
         }

   for  (j = (* ct) - 1;  j > i;  j --)
     stat [j + 1] = stat [j];

   stat [i + 1] . id = id;
   stat [i + 1] . hits = 1;
   stat [i + 1] . bases = bases;
   (* ct) ++;

   return;
  }



static void  Analyze_Screen_Items
    (int fid, Screen_Item_t * screen, int screen_ct)

//  Update global counters and output messages regarding the
//  screen items in  screen [0 .. (screen_ct - 1)]  for
//  fragment  fid .

  {
   int  i, j, hi, lo;
   int  bases, class_ct, item_ct;

   if  (screen_ct == 0)
       return;

   Screened_Frag_Ct ++;

   for  (i = 0;  i < screen_ct - 1;  i ++)
     for  (j = i + 1;  j < screen_ct;  j ++)
       if  (screen [i] . screen_id > screen [j] . screen_id
              || (screen [i] . screen_id == screen [j] . screen_id
                  && screen [i] . lo > screen [j] . lo))
        {
         Screen_Item_t  save;

         save = screen [i];
         screen [i] = screen [j];
         screen [j] = save;
        }

   screen [screen_ct] . screen_id = INT_MAX;
   lo = screen [0] . lo;
   hi = screen [0] . hi;
   item_ct = 0;
   bases = 0;

   for  (i = 1;  i <= screen_ct;  i ++)
     if  (screen [i] . screen_id != screen [i - 1] . screen_id)
         {
          bases += hi - lo;
          Add_Entry (Item_Stats, & Item_Entries, screen [i - 1] . screen_id,
                     bases);
          item_ct ++;
          lo = screen [i] . lo;
          hi = screen [i] . hi;
          bases = 0;
         }
     else if  (screen [i] . lo <= hi)   // overlap
         {
          if  (screen [i] . hi > hi)
              hi = screen [i] . hi;
         }
       else
         {
          bases += hi - lo;
          lo = screen [i] . lo;
          hi = screen [i] . hi;
         }

   if  (Long_Listing && item_ct > 1)
       {
        printf ("fid:%6d  screened by ", fid);
        for  (i = 0;  i < screen_ct;  i ++)
          if  (screen [i] . screen_id != screen [i + 1] . screen_id)
              printf (" %6d", screen [i] . screen_id);
        printf ("\n");
       }


   // Redo by class instead of item

   for  (i = 0;  i < screen_ct - 1;  i ++)
     for  (j = i + 1;  j < screen_ct;  j ++)
       if  (screen [i] . screen_class > screen [j] . screen_class
              || (screen [i] . screen_class == screen [j] . screen_class
                  && screen [i] . lo > screen [j] . lo))
        {
         Screen_Item_t  save;

         save = screen [i];
         screen [i] = screen [j];
         screen [j] = save;
        }

   screen [screen_ct] . screen_class = INT_MAX;
   lo = screen [0] . lo;
   hi = screen [0] . hi;
   class_ct = 0;
   bases = 0;

   for  (i = 1;  i <= screen_ct;  i ++)
     if  (screen [i] . screen_class != screen [i - 1] . screen_class)
         {
          bases += hi - lo;
          Add_Entry (Class_Stats, & Class_Entries, screen [i - 1] . screen_class,
                     bases);
          class_ct ++;
          lo = screen [i] . lo;
          hi = screen [i] . hi;
          bases = 0;
         }
     else if  (screen [i] . lo <= hi)   // overlap
         {
          if  (screen [i] . hi > hi)
              hi = screen [i] . hi;
         }
       else
         {
          bases += hi - lo;
          lo = screen [i] . lo;
          hi = screen [i] . hi;
         }


   // Redo for everything

   for  (i = 0;  i < screen_ct - 1;  i ++)
     for  (j = i + 1;  j < screen_ct;  j ++)
       if  (screen [i] . lo > screen [j] . lo)
        {
         Screen_Item_t  save;

         save = screen [i];
         screen [i] = screen [j];
         screen [j] = save;
        }

   lo = screen [0] . lo;
   hi = screen [0] . hi;
   bases = 0;

   for  (i = 1;  i < screen_ct;  i ++)
     if  (screen [i] . lo <= hi)   // overlap
         {
          if  (screen [i] . hi > hi)
              hi = screen [i] . hi;
         }
       else
         {
          bases += hi - lo;
          lo = screen [i] . lo;
          hi = screen [i] . hi;
         }

   bases += hi - lo;
   Total_Screened_Bases += bases;

   return;
  }



static void  Commatize
    (long int n, char buff [], int buff_len)

//  Put printable representation of  n  into  buff  if it fits
//  within length  buff_len , with commas every 3 digits.
//  If it doesn't fit, fill  buff  with stars.

  {
   char  tmp [1000];
   int  i, len;

   if  (n == 0)
       {
        strcpy (buff, "0");
        return;
       }

   len = 0;
   while  (n > 0)
     {
      if  (len % 4 == 3)
          tmp [len ++] = ',';
      tmp [len ++] = '0' + n % 10;
      n /= 10;
     }

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



