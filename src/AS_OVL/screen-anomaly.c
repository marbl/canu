
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
* Module:  screen-anomaly.c
* Description:
*   Reads a stream of screened fragment messages (i.e., a
*   .urc  file) and a list of frag IDs and prints the screen information
*   of the fragments on the list.
*
*    Programmer:  A. Delcher
*       Written:  26 Aug 99
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: screen-anomaly.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/

static char fileID[] = "$Id: screen-anomaly.c,v 1.1.1.1 2004-04-14 13:52:45 catmandew Exp $";

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
#define  DEFAULT_LONG_LISTING  TRUE
#define  HARD_RELEVANCE_MASK   AS_OVL_HEED_RPT
#define  MAX_SCREEN_ITEMS      10000
#define  REPORT_INTERVAL       10000
#define  MAX_SCREEN            10000
  // Maximum number of screen items a fragment can have


typedef struct
  {
   int  screen_id;
   int  screen_class;
   int  lo, hi;
  }  Screen_Item_t;

typedef  struct
  {
   unsigned  bgn : 16;
   unsigned  end : 16;
  }  Screen_Range_t;

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
static void  Coalesce_Screen_Info
    (Screen_Range_t space [], int * top);
static void  Commatize
    (long int n, char buff [], int buff_len);



int main  (int argc, char * argv [])

  {
   FILE  * urcfile, * idlist_file, * uid_file;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   char  line [BUFF_LEN], buff [BUFF_LEN];
   long int  clear_sum = 0;
   int  curr_id;
   int  ch, error_flag;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "")) != EOF))
     switch  (ch)
       {
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
   read_msg_fn = InputFileType_AS (urcfile);

   optind ++;
   idlist_file = fopen (argv [optind], "r");
   if  (idlist_file == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file %s\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }
   curr_id = -1;

   uid_file = fopen ("uids.list", "w");
   assert (uid_file != NULL);

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
           int  frag_id,  frag_len;

           Total_Frags ++;

           sfg_mesg = (ScreenedFragMesg*) gmesg -> m;
           frag_id = sfg_mesg -> iaccession;
           frag_len = sfg_mesg -> clear_rng . end - sfg_mesg -> clear_rng . bgn;
           clear_sum += frag_len;

           while  (curr_id < frag_id)
             {
              if  (fgets (line, BUFF_LEN, idlist_file) == NULL)
                  curr_id = INT_MAX;
                else
                  sscanf (line, "%d", & curr_id);
             }

           if  (curr_id == frag_id)
               {
                Screen_Range_t  hard_screen [MAX_SCREEN], all_screen [MAX_SCREEN];
                int  hard_ct = 0, all_ct = 0;
                int  all_len, hard_len;
                int  i;

                fprintf (uid_file, F_UID "\n", sfg_mesg -> eaccession);
                for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
                  {
                   all_screen [all_ct] . bgn = p -> where . bgn;
                   all_screen [all_ct] . end = p -> where . end;
                   all_ct ++;
                   if  (p -> relevance & AS_OVL_HEED_RPT)
                       {
                        hard_screen [hard_ct] . bgn = p -> where . bgn;
                        hard_screen [hard_ct] . end = p -> where . end;
                        hard_ct ++;
                       }
                  }
                Coalesce_Screen_Info (all_screen, & all_ct);
                all_len = 0;
                for  (i = 0;  i < all_ct;  i ++)
                  all_len += all_screen [i] . end - all_screen [i] . bgn;
                hard_len = 0;
                Coalesce_Screen_Info (hard_screen, & hard_ct);
                  hard_len += hard_screen [i] . end - hard_screen [i] . bgn;

                printf (
  "frag %7d  len = %4d   hard_screen = %3.0f%%   all_screen = %3.0f%%\n",
                        frag_id, frag_len,
                        (100.0 * hard_len) / frag_len,
                        (100.0 * all_len) / frag_len);
                for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
                  printf ("  %7d  whr = %4d,%4d  pof = %6d,%6d  rel = %d\n",
                          p -> iwhat,
                          p -> where . bgn,
                          p -> where . end,
                          p -> portion_of . bgn,
                          p -> portion_of . end,
                          p -> relevance);
               }
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

   fclose (urcfile);
   fclose (idlist_file);
   fclose (uid_file);

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



static void  Coalesce_Screen_Info
    (Screen_Range_t space [], int * top)

//  Sort entries in  space [0 .. (* top) - 1]
//  Then combine overlapping adjacent entries.  Reduce  (* top)
//  by the number of combined entries.

  {
   int  i, j, hi;

   hi = (* top - 1);

   //  Sort ascending by  bgn  value
   for  (i = 0;  i < hi;  i ++)
     for  (j = i + 1;  j <= hi;  j ++)
       if  (space [i] . bgn > space [j] . bgn)
           {
            Screen_Range_t  save = space [i];

            space [i] = space [j];
            space [j] = save;
           }

   //  Combine overlapping entries
   for  (i = 0, j = i + 1;  j <= hi;  j ++)
     {
      if  (space [j] . bgn > space [i] . end)
          space [++ i] = space [j];
        else if  (space [j] . end
                    > space [i] . end)
          space [i] . end = space [j] . end;
     }

   (* top) = i + 1;

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



