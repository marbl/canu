
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
*   Collect and output the number of overlaps for each fragment.
*   Input is from the standard input and is the format produced by
*   get-olaps
*
*    Programmer:  A. Delcher
*       Written:  6 Jul 99
*  Last Revised:  26 Aug 99  Get separate degrees for each end of frag
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: get-degrees.c,v 1.2 2004-09-23 20:25:26 mcschatz Exp $
 * $Revision: 1.2 $
*/

static char fileID[] = "$Id: get-degrees.c,v 1.2 2004-09-23 20:25:26 mcschatz Exp $";

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
#include  <limits.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "OlapStoreOVL.h"


#define  INITIAL_NUM_FRAGS  1000000
#define  MAX_LINE  10000
#define  MAX_RECORD_LENGTH 200


typedef struct
  {
   int  a_degree, b_degree;
  }  Degree_t;


static char *  OVL_Store_Path = NULL;
  // If using overlap store, will be set to directory
  // where it is (by the -S option)


static void  Get_Degrees_From_Store
    (void);
static int  Get_Olap
    (int * from, int * to, int * a_hang, int * b_hang, char * orient);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Usage
    (char * command);



int  main
    (int argc, char * argv [])

  {
   int  num_frags = INITIAL_NUM_FRAGS;
   int  min_frag = INT_MAX, max_frag = INT_MIN;
   int  from, to, a_hang, b_hang;
   int  i;
   char  orient [MAX_LINE];
   Degree_t  * degree;

   Parse_Command_Line  (argc, argv);

   degree = (Degree_t *) calloc (num_frags, sizeof (Degree_t));
   assert (degree != NULL);

   if  (OVL_Store_Path != NULL)
       {
        Get_Degrees_From_Store ();

        return  0;
       }

   while  (Get_Olap (& from, & to, & a_hang, & b_hang, orient))
     {
      if  (from < min_frag)
          min_frag = from;
      if  (from > max_frag)
          max_frag = from;
      if  (to < min_frag)
          min_frag = to;
      if  (to > max_frag)
          max_frag = to;

      if  (from >= num_frags)
          {
           degree = (Degree_t *) realloc (degree, (from + 1) * sizeof (Degree_t));
           assert (degree != NULL);
           for  (i = num_frags;  i <= from;  i ++)
             degree [i] . a_degree = degree [i] . b_degree = 0;
           num_frags = from + 1;
          }

      if  (to >= num_frags)
          {
           degree = (Degree_t *) realloc (degree, (to + 1) * sizeof (Degree_t));
           assert (degree != NULL);
           for  (i = num_frags;  i <= to;  i ++)
             degree [i] . a_degree = degree [i] . b_degree = 0;
           num_frags = to + 1;
          }

      assert (0 <= from && from < num_frags && 0 <= to && to < num_frags);

      switch  (orient [0])
        {
         case  'I' :
           degree [to] . b_degree ++;
           if  (b_hang <= 0)
               degree [to] . a_degree ++;
           if  (b_hang >= 0)
               degree [from] . b_degree ++;
           break;
         case  'O' :
           degree [to] . a_degree ++;
           if  (b_hang <= 0)
               degree [to] . b_degree ++;
           if  (b_hang >= 0)
               degree [from] . a_degree ++;
           break;
         case  'N' :
           degree [to] . a_degree ++;
           if  (b_hang <= 0)
               degree [to] . b_degree ++;
           if  (b_hang >= 0)
               degree [from] . b_degree ++;
           break;
         default :
           fprintf (stderr,
                    "\a\a\aERROR:  Bad orientation = \'%c\' for a = %d  b = %d\n",
                    orient [0], from, to);
        }
     }

   if  (1 < min_frag)
       min_frag = 1;
   for  (i = min_frag;  i <= max_frag;  i ++)
     printf ("%8d %6d %6d\n", i, degree [i] . a_degree, degree [i] . b_degree);

   return  0;
  }



static void  Get_Degrees_From_Store
    (void)

//  Read overlaps in  OVL_Store_Path  and print their degrees
//  to standard output

  {
   OVL_Store_t  * ovl_store;
   OVL_Stream_t  * ovl_stream;
   Long_Olap_Data_t  olap;
   int  prev_a_iid, left_degr, right_degr;
   
   ovl_store = New_OVL_Store ();
   if  (Open_OVL_Store (ovl_store, OVL_Store_Path) != 0)
       {
        fprintf (stderr, "ERROR:  Failed to open overlap store \"%s\"\n",
                 OVL_Store_Path);
        exit (EXIT_FAILURE);
       }
   ovl_stream = New_OVL_Stream ();
   Init_OVL_Stream (ovl_stream, 1, Last_Frag_In_OVL_Store (ovl_store),
                    ovl_store);

   prev_a_iid = 1;
   left_degr = right_degr = 0;
   while  (Next_From_OVL_Stream (& olap, ovl_stream))
     {
      if  (prev_a_iid < olap . a_iid)
          {
           printf ("%8d %6d %6d\n", prev_a_iid, left_degr, right_degr);
           left_degr = right_degr = 0;
           while  (++ prev_a_iid < olap . a_iid)
             printf ("%8d %6d %6d\n", prev_a_iid, 0, 0);
          }
      assert (prev_a_iid == olap . a_iid);
      if  (olap . a_hang <= 0)
          left_degr ++;
      if  (olap . b_hang >= 0)
          right_degr ++;
     }
   printf ("%8d %6d %6d\n", prev_a_iid, left_degr, right_degr);

   Free_OVL_Stream (ovl_stream);
   Free_OVL_Store (ovl_store);

   return;
  }



static int  Get_Olap
    (int * from, int * to, int * a_hang, int * b_hang, char * orient)

//  Read next overlap and put its information into  (* from) , (* to) ,
//  (* a_hang) , (* b_hang)  and  orient .  Return  TRUE  if an
//  overlap is read;  FALSE , otherwise.

  {
   char  line [MAX_LINE];

   if  (fgets (line, MAX_RECORD_LENGTH, stdin) == NULL)
       return  FALSE;
   if  (sscanf (line, "%d %d %d %d %s", from, to, a_hang, b_hang, orient) != 5)
       return  FALSE;

   return  TRUE;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "S:")) != EOF))
     switch  (ch)
       {
        case  'S' :
          OVL_Store_Path = optarg;
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE: %s [-S OVLStore]\n"
           "\n"
           "Reads overlaps from standard input, or store specified\n"
           "by  -S  option, and produces a list of overlap degrees\n"
           "on both ends of each fragement.  Output goes to standard\n"
           "output\n"
           "\n"
           "Options:\n"
           "-S   Specify overlap store from which to get overlaps\n",
           command);

   return;
  }



