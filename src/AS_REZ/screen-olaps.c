
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
* Module:  screen-olaps.c
* Description:
*   Reads a stream of screened fragment messages (i.e., a
*   .urc  file) and a prints patterns at boundaries of screen
*   regions depending on how various #if's are modified.
*   Not production-quality code.
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
 * $Id: screen-olaps.c,v 1.3 2005-03-22 19:08:04 jason_miller Exp $
 * $Revision: 1.3 $
*/

static char fileID[] = "$Id: screen-olaps.c,v 1.3 2005-03-22 19:08:04 jason_miller Exp $";

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


#define  MAX_SCREEN            1000
#define  ONLY_UNSCREENED       0


typedef  struct
  {
   int  lo, hi;
  }  Screen_Info;


static int  Start_Ct = 0;
static int  Iso_Ct = 0;


static void  Add_Coverage
    (Screen_Info screen_list [], int lo, int hi, int * n);
static void  Make_Lowercase
    (char * s);
static void  Rev_Complement
    (char * s);
static void  stripWhiteSpace
    (char * target, char * source, int maxlen);



int main  (int argc, char * argv [])

  {
   FILE  * urcfile;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   Screen_Info  screen_list [MAX_SCREEN];
   int  list_len;
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

   while  (read_msg_fn (urcfile, & gmesg) != EOF && gmesg != NULL)
     if  (gmesg -> t == MESG_SFG)
         {
          IntScreenMatch  * p;
          ScreenedFragMesg  * sfg_mesg;
          char  seq [AS_READ_MAX_LEN * 2];
          char  buff [AS_READ_MAX_LEN * 2];
          int  frag_id, frag_len;
          IntScreenMatch  * other_matches [1000];
          int  other_len, found_100001, found_700031, found_start;
          

	  sfg_mesg = gmesg -> m;
          frag_id = sfg_mesg -> iaccession;

          stripWhiteSpace (seq, sfg_mesg -> sequence, AS_READ_MAX_LEN * 2);
          frag_len = sfg_mesg -> clear_rng . end - sfg_mesg -> clear_rng . bgn;
          if  (sfg_mesg -> clear_rng . bgn > 0)
              memmove (seq, seq + sfg_mesg -> clear_rng . bgn, frag_len);
          seq [frag_len] = '\0';

          list_len = 0;
          other_len = 0;
          found_100001 = FALSE;
          found_700031 = FALSE;
          found_start = FALSE;
          for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
            {
             Add_Coverage (screen_list, p -> where . bgn, p -> where . end,
                           & list_len);
             if  (p -> iwhat == 100001
                    || p -> iwhat == 700031)
                 other_matches [other_len ++] = p;
             if  (p -> iwhat == 100001)
                 found_100001 = TRUE;
             if  (p -> iwhat == 700031)
                 found_700031 = TRUE;
if  (p -> iwhat == 700031 && p -> portion_of . bgn == 0
       && ((p -> direction == AS_FORWARD && p -> where . bgn >= 50)
          || (p -> direction == AS_REVERSE && p -> where . end <= frag_len - 50)))
    {
     Start_Ct ++;
     found_start = TRUE;
    }
            }
#if  0
if  (found_start && ! found_100001)
    {
     int  i;

     printf ("frag = %7d ", frag_id);
     for  (i = 0;  i < other_len;  i ++)
          printf ("  [%6d: %4d %4d %6d %6d %c]",
                  other_matches [i] -> iwhat,
                  other_matches [i] -> where . bgn,
                  other_matches [i] -> where . end,
                  other_matches [i] -> portion_of . bgn,
                  other_matches [i] -> portion_of . end,
                  other_matches [i] -> direction == AS_FORWARD ? 'f' : 'r');
     putchar ('\n');

     if  (other_len == 1
            && ((other_matches [0] -> direction == AS_FORWARD
                   && other_matches [0] -> where . bgn >= 50)
                 || (other_matches [0] -> direction == AS_REVERSE
                       && other_matches [0] -> where . end <= frag_len - 50)))
         {
          printf ("frag = %7d ", frag_id);
          printf ("  [%6d: %4d %4d %6d %6d %c]",
                  other_matches [0] -> iwhat,
                  other_matches [0] -> where . bgn,
                  other_matches [0] -> where . end,
                  other_matches [0] -> portion_of . bgn,
                  other_matches [0] -> portion_of . end,
                  other_matches [0] -> direction == AS_FORWARD ? 'f' : 'r');
          putchar ('\n');
          Iso_Ct ++;
         }
    }
continue;
          if  (found_100001 && found_700031)
              {
               int  i, j, top;
               IntScreenMatch  * stack [100];

               for  (i = 0;  i < other_len - 1;  i ++)
                 for  (j = i + 1;  j < other_len;  j ++)
                   if  (other_matches [i] -> where . bgn >
                          other_matches [j] -> where . bgn)
                       {
                        IntScreenMatch  * save;

                        save = other_matches [i];
                        other_matches [i] = other_matches [j];
                        other_matches [j] = save;
                       }

               j = -1;
               top = 0;
               for  (i = 1;  i < other_len;  i ++)
                 if  (other_matches [i - 1] -> iwhat
                        != other_matches [i] -> iwhat)
                     {
                      if  (j != i - 1)
                          stack [top ++] = other_matches [i - 1];
                      stack [top ++] = other_matches [i];
                      j = i;
                     }

               printf ("frag = %7d ", frag_id);
               if  (stack [0] -> iwhat == 100001)
                   {
                    for  (i = 0;  i < top;  i ++)
                      printf ("  [%6d: %4d %4d %6d %6d %c]",
                              stack [i] -> iwhat,
                              stack [i] -> where . bgn,
                              stack [i] -> where . end,
                              stack [i] -> portion_of . bgn,
                              stack [i] -> portion_of . end,
                              stack [i] -> direction == AS_FORWARD ? 'f' : 'r');
                   }
                 else
                   {
                    for  (i = top - 1;  i >= 0;  i --)
                      printf ("  [%6d: %4d %4d %6d %6d %c]",
                              stack [i] -> iwhat,
                              stack [i] -> where . bgn,
                              stack [i] -> where . end,
                              stack [i] -> portion_of . bgn,
                              stack [i] -> portion_of . end,
                              stack [i] -> direction == AS_FORWARD ? 'f' : 'r');
                   }
               putchar ('\n');
              }
          continue;
#endif


          for  (p = sfg_mesg -> screened;  p != NULL;  p = p -> next)
            if  (p -> iwhat == 700031)
                {
                 int  i, len, start, stop;
                 int  lo = p -> portion_of . bgn;
#if  0
                 int  hi = p -> portion_of . end;
printf ("%4d %6d\n", lo, hi);
#else
#if  1
                 if  (lo == 0)
                     {
                      if  (p -> direction == AS_FORWARD)
                          {
                           int  match_bgn = p -> where . bgn;

                           start = 0;
#if  ONLY_UNSCREENED
                           for  (i = 0;
                                   i < list_len
                                     && screen_list [i] . hi < match_bgn;
                                   i ++)
                             start = screen_list [i] . hi;
                           if  (match_bgn > screen_list [i] . lo
                                  || match_bgn - start < 20)
                               continue;
                           if  (match_bgn - start > 40)
                               start = match_bgn - 40;
#else
                           if  (start < match_bgn - 40)
                               start = match_bgn - 40;
#endif

                           len = match_bgn - start;
                           strncpy (buff, seq + start, len);
                           buff [len] = '\0';
                           Make_Lowercase (buff);
                           printf ("%40s | ", buff);

                           stop = match_bgn + 60;
                           if  (stop > frag_len)
                               stop = frag_len;
                           seq [stop] = '\0';
                           strcpy (buff, seq + match_bgn);
                           Make_Lowercase (buff);
                           printf ("%-60s %6d . %3d F", buff,
                                   frag_id, match_bgn);
                          }
                      else if  (p -> direction == AS_REVERSE)
                          {
                           int  match_end = p -> where . end;

                           stop = frag_len;
#if  ONLY_UNSCREENED
                           for  (i = list_len - 1;
                                   i >= 0
                                     && screen_list [i] . lo > match_end;
                                   i --)
                             stop = screen_list [i] . lo;
                           if  (match_end < screen_list [i] . hi
                                  || stop - match_end < 20)
                               continue;
                           if  (stop - match_end > 40)
                               stop = match_end + 40;
#else
                           if  (match_end + 40 < stop)
                               stop = match_end + 40;
#endif

                           len = stop - match_end;
                           seq [stop] = '\0';
                           strcpy (buff, seq + match_end);
                           Make_Lowercase (buff);
                           Rev_Complement (buff);
                           printf ("%40s | ", buff);

                           start = match_end - 60;
                           if  (start < 0)
                               start = 0;
                           strncpy (buff, seq + start, match_end - start);
                           buff [match_end - start] = '\0';
                           Make_Lowercase (buff);
                           Rev_Complement (buff);
                           printf ("%-60s %6d . %3d R", buff,
                                   frag_id, p -> where . end);
                          }
                      for  (i = 0;  i < other_len;  i ++)
                        if  (other_matches [i] -> iwhat != 700031
                             || (other_matches [i] -> iwhat == 700031
                                 && other_matches [i] -> portion_of . end >= 5000))
                            printf ("  [%6d: %4d %4d %6d %6d]",
                                    other_matches [i] -> iwhat,
                                    other_matches [i] -> where . bgn,
                                    other_matches [i] -> where . end,
                                    other_matches [i] -> portion_of . bgn,
                                    other_matches [i] -> portion_of . end);
                      putchar ('\n');
                     }
#else
                 if  (hi >= 5354)
                     {
                      if  (p -> direction == AS_REVERSE)
                          {
                           int  match_bgn = p -> where . bgn;

                           start = 0;
                           for  (i = 0;
                                   i < list_len
                                     && screen_list [i] . hi < match_bgn;
                                   i ++)
                             start = screen_list [i] . hi;
                           if  (match_bgn > screen_list [i] . lo
                                  || match_bgn - start < 20)
                               continue;
                           if  (match_bgn - start > 40)
                               start = match_bgn - 40;

                           len = match_bgn - start;
                           strncpy (buff, seq + start, len);
                           strcpy (buff + len, " | ");
                           seq [match_bgn + 20] = '\0';
                           strcat (buff, seq + match_bgn);
                           Make_Lowercase (buff);
                           Rev_Complement (buff);
                           printf ("%-63s %6d . %3d F", buff,
                                   frag_id, match_bgn);
                          }
                      else if  (p -> direction == AS_FORWARD)
                          {
                           int  match_end = p -> where . end;

                           stop = frag_len;
                           for  (i = list_len - 1;
                                   i >= 0
                                     && screen_list [i] . lo > match_end;
                                   i --)
                             stop = screen_list [i] . lo;
                           if  (match_end < screen_list [i] . hi
                                  || stop - match_end < 20)
                               continue;
                           if  (stop - match_end > 40)
                               stop = match_end + 40;

                           len = stop - match_end;
                           strncpy (buff, seq + match_end - 20, 20);
                           strcpy (buff + 20, " | ");
                           seq [stop] = '\0';
                           strcat (buff, seq + match_end);
                           Make_Lowercase (buff);
                           printf ("%-63s %6d . %3d R", buff,
                                   frag_id, p -> where . end);
                          }
                      for  (i = 0;  i < other_len;  i ++)
                        if  (other_matches [i] -> iwhat != 700031
                              || other_matches [i] -> iwhat == 700031
                                   && other_matches [i] -> portion_of . bgn <= 50)
                            printf ("  [%6d: %4d %4d %6d %6d]",
                                    other_matches [i] -> iwhat,
                                    other_matches [i] -> where . bgn,
                                    other_matches [i] -> where . end,
                                    other_matches [i] -> portion_of . bgn,
                                    other_matches [i] -> portion_of . end);
                      putchar ('\n');
                     }
#endif
#endif
                }
         }

   fclose (urcfile);

printf ("Start_Ct = %d   Iso_Ct = %d\n", Start_Ct, Iso_Ct);

   return  0;
  }



static void  Add_Coverage
    (Screen_Info list [], int lo, int hi, int * n)

//  Add interval  lo .. hi  to intervals in  list  so that  list
//  contains the union of the intervals.  Entries in  list  are
//  in ascending order and do not overlap.   (* n)  is the number
//  of entries in  list  and is either increased or decreased.

  {
   int  i, j;

   for  (i = 0;  i < (* n);  i ++)
     if  (lo <= list [i] . hi)
         break;

   if  (i == (* n))
       {
        assert ((* n) < MAX_SCREEN);
        list [i] . lo = lo;
        list [i] . hi = hi;
        (* n) ++;
        return;
       }

   if  (hi < list [i] . lo)
       {
        assert ((* n) < MAX_SCREEN);
        for  (j = (* n);  j > i;  j --)
          list [j] = list [j - 1];
        list [i] . lo = lo;
        list [i] . hi = hi;
        (* n) ++;
        return;
       }

   if  (lo < list [i] . lo)
       list [i] . lo = lo;
   if  (hi > list [i] . hi)
       list [i] . hi = hi;
   for  (j = i + 1;  j < (* n) && list [j] . lo <= hi;  j ++)
     if  (list [i] . hi < list [j] . hi)
         list [i] . hi = list [j] . hi;
   while  (j < (* n))
     list [++ i] = list [j ++];
   (* n) = i + 1;

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



/* Function stripWhiteSpace:
   Input:  source   string of maximum length maxlen
           maxlen   maximum length of source
   Output: target
   
   Description:
     Copy non-space characters from source to target.
*/

static void  stripWhiteSpace
    (char *target, char *source, int maxlen)

  {
  int i = 0;
  *target = '\0';
  while(i < maxlen){
    if(!isspace(*source)){
      *target++ = *source;
      i++;
    }
    if(*source == '\0')
      break;
    source++;
  }

}
