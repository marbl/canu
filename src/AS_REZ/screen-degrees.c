
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
* Module:  screen-degrees.c
* Description:
*   Reads a stream of screened fragment messages (i.e., a
*   .urc  file) and a list of the numbers of overlaps for each
*   fragment, and puts them together.  Output goes to the standard
*   output.
*
*    Programmer:  A. Delcher
*       Written:  6 Jul 99
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: screen-degrees.c,v 1.4 2005-03-22 19:49:26 jason_miller Exp $
 * $Revision: 1.4 $
*/

static char fileID[] = "$Id: screen-degrees.c,v 1.4 2005-03-22 19:49:26 jason_miller Exp $";

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
#define  END_SLUSH             40
  // A screen item that matches within this many bases of the end of the
  // fragment counts as matching that end of the fragment
#define  HARD_RELEVANCE_MASK   AS_OVL_HEED_RPT
#define  MAX_SCREEN            10000
#define  REPORT_INTERVAL       10000
#define  TAG_LEN               1000


typedef  struct
  {
   char  * tag;
   int  min_degree, max_degree;
   double  degree_sum;
   int  degree_ct;
  }  Match_Entry_t;

typedef  struct
  {
   int  what;
   int  min_lo, max_hi;
  }  Range_t;


static void  Add_Match_Entry
    (char * tag, int degree, Match_Entry_t a [], int * n);
static void  Commatize
    (long int n, char buff [], int buff_len);
void  Insert_Sorted
    (int x, int a [], int * n);
static void  Make_Lowercase
    (char * s);
static void  Rev_Complement
    (char * s);
static void  Update_Range
    (int what, int lo, int hi, Range_t a [], int * n);


static int  Ignore_Relevance = FALSE;
static int  Long_Listing = DEFAULT_LONG_LISTING;
static unsigned  Relevance_Mask = HARD_RELEVANCE_MASK;
static int  Total_Frags = 0;
static FILE  * Unscreened_Degree_File;
static int  Unscreened_Histo = FALSE;


int main  (int argc, char * argv [])

  {
   FILE  * urcfile, * degfile;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   Range_t  screen_range [MAX_SCREEN];
   Match_Entry_t  match_entry [2 * MAX_SCREEN];
   char  buff [BUFF_LEN];
   long int  total_overlaps = 0;
   int  ch, error_flag, range_len = 0, entry_len = 0;
   int  i;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "alu")) != EOF))
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

        case  'u':
          // Create data file of unscreened degrees for histogram

          Unscreened_Histo = TRUE;
          break;

        case  '?':
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind > argc - 2)
       {
        fprintf (stderr, 
                 "USAGE:  %s <urc-file> <olap-degree-file>\n", 
                 argv [0]);
        fprintf (stderr, "Options:\n");
        fprintf (stderr, " -a  Use *all* matches regardless of relevance\n");
        exit (EXIT_FAILURE);
       }

   urcfile = fopen (argv [optind], "r");
   if  (urcfile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file \"%s\"\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }
   read_msg_fn = InputFileType_AS (urcfile);

   optind ++;
   degfile = fopen (argv [optind], "r");
   if  (degfile == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file \"%s\"\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }

   if  (Unscreened_Histo)
       {
        Unscreened_Degree_File = fopen ("unscreened.histo", "w");
        if  (Unscreened_Degree_File == NULL)
            {
             fprintf (stderr, "ERROR:  Could not open file \"unscreened.histo\"\n");
             exit (EXIT_FAILURE);
            }
        fprintf (Unscreened_Degree_File, "%s Unscreened Degrees\n",
                 argv [optind - 1]);
       }

   if  (! Long_Listing)
       fprintf (stderr, "\n");
   while  (read_msg_fn (urcfile, & gmesg) != EOF && gmesg != NULL)
     if  (gmesg -> t == MESG_SFG)
         {
          IntScreenMatch  * p;
          ScreenedFragMesg  * sfg_mesg;
          char  tag [TAG_LEN];
          int  a_screen_list [MAX_SCREEN], b_screen_list [MAX_SCREEN];
          int  a_list_len, b_list_len;
          int  frag_id, id, a_degree, b_degree, frag_len, scan_result;
          int  j, len;

	  sfg_mesg = (ScreenedFragMesg *) gmesg -> m;
          frag_id = sfg_mesg -> iaccession;
          frag_len = sfg_mesg -> clear_rng . end - sfg_mesg -> clear_rng . bgn;

          scan_result = fscanf (degfile, "%d %d %d", & id, & a_degree, & b_degree);
          if  (scan_result == EOF)
              {
               id = frag_id;
               a_degree = b_degree = 0;
              }
          else if  (scan_result != 3)
              {
               fprintf (stderr, "ERROR:  Only %d entries in degree file\n",
                        scan_result);
               exit (EXIT_FAILURE);
              }
          if  (id != frag_id)
              {
               fprintf (stderr, "ID error:  urc frag id = %d  degree frag id = %d\n",
                        frag_id, id);
               exit (EXIT_FAILURE);
              }

          if  (! Long_Listing && Total_Frags % REPORT_INTERVAL == 0)
              {
               Commatize (Total_Frags, buff, BUFF_LEN);
               fprintf (stderr, "\rFrags read = %s", buff);
              }

          Total_Frags ++;
          total_overlaps += a_degree + b_degree;

          a_list_len = b_list_len = 0;
          for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
            if  ((p -> relevance & Relevance_Mask) || Ignore_Relevance)
                {
                 if  (p -> where . bgn < END_SLUSH)
                     Insert_Sorted (p -> iwhat, a_screen_list, & a_list_len);
                 if  (p -> where . end > frag_len - END_SLUSH)
                     Insert_Sorted (p -> iwhat, b_screen_list, & b_list_len);
                 Update_Range (p -> iwhat, p -> portion_of . bgn, p -> portion_of . end,
                               screen_range, & range_len);
                }

          j = 0;
          for  (i = 0;  i < a_list_len;  i ++)
            {
             if  (i > 0)
                 tag [j ++] = ',';
             sprintf (tag + j, "%d%n", a_screen_list [i], & len);
             j += len;
            }
          tag [j] = '\0';
          assert (j < TAG_LEN);
          if  (Long_Listing)
              printf ("%8d %6d %-25s : ",
                      frag_id, a_degree, tag);
          Add_Match_Entry (tag, a_degree, match_entry, & entry_len);

          j = 0;
          for  (i = 0;  i < b_list_len;  i ++)
            {
             if  (i > 0)
                 tag [j ++] = ',';
             sprintf (tag + j, "%d%n", b_screen_list [i], & len);
             j += len;
            }
          tag [j] = '\0';
          assert (j < TAG_LEN);
          if  (Long_Listing)
              printf ("%6d %-25s\n",
                      b_degree, tag);
          Add_Match_Entry (tag, b_degree, match_entry, & entry_len);
         }

   if  (! Long_Listing)
       {
        Commatize (Total_Frags, buff, BUFF_LEN);
        fprintf (stderr, "\rFrags read = %s\n", buff);
       }

   fclose (urcfile);
   fclose (degfile);
   if  (Unscreened_Histo)
       fclose (Unscreened_Degree_File);

   if  (strcmp (match_entry [0] . tag, "") == 0)
       match_entry [0] . tag = strdup ("<Unscreened>");
   printf ("\n%-24s  %9s  %9s  %9s  %9s %12s\n",
           "Screened By", "Frag Ends", "MinOlaps", "MaxOlaps",
           "AvgOlaps", "TotalOlaps");
   for  (i = 0;  i < entry_len;  i ++)
     {
      int  tag_len;

      tag_len = strlen (match_entry [i] . tag);
      if  (tag_len <= 24)
          printf ("%-24s  %9d  %9d  %9d  %9.1f %12.0f\n",
                  match_entry [i] . tag, 
                  match_entry [i] . degree_ct,
                  match_entry [i] . min_degree,
                  match_entry [i] . max_degree,
                  match_entry [i] . degree_sum / match_entry [i] . degree_ct,
                  match_entry [i] . degree_sum);
        else
          printf ("%-24s\n%-24s  %9d  %9d  %9d  %9.1f %12.0f\n",
                  match_entry [i] . tag, "",
                  match_entry [i] . degree_ct,
                  match_entry [i] . min_degree,
                  match_entry [i] . max_degree,
                  match_entry [i] . degree_sum / match_entry [i] . degree_ct,
                  match_entry [i] . degree_sum);
     }
    printf ("\n%24s  %9d  %9s  %9s  %9.1f %12ld\n",
            "Total:", 2 * Total_Frags,
            "", "", 
            total_overlaps / (2.0 * Total_Frags),
            total_overlaps);

   Commatize (total_overlaps / 2, buff, BUFF_LEN);
   printf ("\nTotal_Overlaps = %s\n", buff);

   printf ("\n\nRange of Screen Item Matched\n");
   printf ("%9s  %9s  %9s\n",
           "Screen ID", "Lo", "Hi");
   for  (i = 0;  i < range_len;  i ++)
     printf ("%9d  %9d  %9d\n",
             screen_range [i] . what,
             screen_range [i] . min_lo,
             screen_range [i] . max_hi);

   return  0;
  }



static void  Add_Match_Entry
    (char * tag, int degree, Match_Entry_t a [], int * n)

//  Add  degree  info  for match identified by  tag  to
//  a [0 .. (*n - 1)] .

  {
   int  i, left, right, mid;

   if  (Unscreened_Histo && strcmp (tag, "") == 0)
       fprintf (Unscreened_Degree_File, "%d\n", degree);

   left = 0;
   right = (* n) - 1;

   while  (left <= right)
     {
      mid = (left + right) / 2;
      if  (strcmp (a [mid] . tag, tag) == 0)
          {
           if  (degree < a [mid] . min_degree)
               a [mid] . min_degree = degree;
           if  (degree > a [mid] . max_degree)
               a [mid] . max_degree = degree;
           a [mid] . degree_sum += degree;
           a [mid] . degree_ct ++;
           return;
          }
      if  (strcmp (tag, a [mid] . tag) < 0)
          right = mid - 1;
        else
          left = mid + 1;
     }

   for  (i = (* n);  i > left;  i --)
     a [i] = a [i - 1];

   a [left] . tag = strdup (tag);;
   a [left] . min_degree = degree;
   a [left] . max_degree = degree;
   a [left] . degree_sum = degree;
   a [left] . degree_ct = 1;
   (* n) ++;

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



void  Insert_Sorted
    (int x, int a [], int * n)

//  Insert value  x  into sorted array  a [] .   n  is the number
//  of entries in  a  and will be incremented unless  x  is alread in
//  a .

  {
   int  i, lo, hi, mid;

   lo = 0;
   hi = (* n) - 1;

   while  (lo <= hi)
     {
      mid = (lo + hi) / 2;
      if  (a [mid] == x)
          return;
      if  (x < a [mid])
          hi = mid - 1;
        else
          lo = mid + 1;
     }

   for  (i = (* n);  i > lo;  i --)
     a [i] = a [i - 1];

   a [lo] = x;
   (* n) ++;

   return;
  }



static char  Complement
    (char ch)

/*  Return the DNA complement of  ch . */

  {
   switch  (tolower ((int) ch))
     {
      case  'a' :
        return  't';
      case  'c' :
        return  'g';
      case  'g' :
        return  'c';
      case  't' :
        return  'a';
      default :
        return  tolower ((int) ch);
     }
  }



static void  Make_Lowercase
    (char * s)

/* Set all characters in string  s  to lowercase */

  {
   int  i;

   for  (i = 0;  s [i] != '\0';  i ++)
     s [i] = tolower (s [i]);

   return;
  }



static void  Rev_Complement
    (char * s)

/* Set string  s  to its DNA reverse complement. */

  {
   char  ch;
   int  i, j, len;

   len = strlen (s);
   for  (i = 0, j = len - 1;  i < j;  i ++, j --)
     {
      ch = Complement (s [i]);
      s [i] = Complement (s [j]);
      s [j] = ch;
     }

   if  (i == j)
       s [i] = Complement (s [i]);

   return;
  }



static void  Update_Range
    (int what, int lo, int hi, Range_t a [], int * n)

//  Search for  what  in  a [0 .. (*n - 1)] .  If
//  present, extend  min_lo  and  max_hi  to include  lo  and
//  hi , if needed.  If not present, add a new entry (in order)
//  and set  min_lo, max_hi  to  lo, hi .

  {
   int  i, left, right, mid;

   left = 0;
   right = (* n) - 1;

   while  (left <= right)
     {
      mid = (left + right) / 2;
      if  (a [mid] . what == what)
          {
           if  (lo < a [mid] . min_lo)
               a [mid] . min_lo = lo;
           if  (hi > a [mid] . max_hi)
               a [mid] . max_hi = hi;
           return;
          }
      if  (what < a [mid] . what)
          right = mid - 1;
        else
          left = mid + 1;
     }

   for  (i = (* n);  i > left;  i --)
     a [i] = a [i - 1];

   a [left] . what = what;
   a [left] . min_lo = lo;
   a [left] . max_hi = hi;
   (* n) ++;

   return;
  }
