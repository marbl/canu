
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
* Module:  frag-anomaly.c
* Description:
*   Reads the stream of messages produced by the overlapper
*   and identifies fragments that have unusual overlap patterns.
*   Eventually will generate a revised version of the overlap
*   file with overlaps involving those fragments removed.
*
*    Programmer:  A. Delcher
*       Written:  30 Aug 99
*  Last Revised:  
* 
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: frag-anomaly.c,v 1.3 2005-03-22 19:07:12 jason_miller Exp $
 * $Revision: 1.3 $
*/

static char fileID[] = "$Id: frag-anomaly.c,v 1.3 2005-03-22 19:07:12 jason_miller Exp $";

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <assert.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>

#include  "AS_global.h"
#include  "AS_MSG_pmesg.h"
#include  "AS_OVL_delcher.h"


#define  MIN_BP_OLAPS             5
  // minimum number of overlaps needed on other end of fragment with
  // 0 overlaps on this end in order to be "weak"
#define  OLAP_SLUSH               30
  // amout of overlap tolerated between fragments from each end of
  // a chimeric-like fragment
#define  SCORE_THRESHOLD          -2.0
  // log probability of overlaps must be less than this on both sides
  // of a chimeric-like fragment
#define  SHOW_SPAN                0
  // if true show overlap that supports span


typedef  struct
  {
   char  end;
   char  orient;
   int16 len;
   int  frag_id;
  }  Olap_Entry_t;

typedef  struct
  {
   Olap_Entry_t  * olap;
   int16  num_olaps;
   signed int  len : 15;
   unsigned  reject : 1;
  }  Frag_Entry_t;


static Frag_Entry_t  * Frag;


static void  Add_Entry
    (Frag_Entry_t * f, char end, char orient, int len, int frag_id);



int main  (int argc, char * argv [])

  {
   FILE  * ovlfile, * outfile;
   char  * infile_name, * outfile_name;
   GenericMesg  * gmesg = NULL;
   MesgReader  read_msg_fn;
   MesgWriter  write_msg_fn;
   GenericMesg  * pmesg;
   AuditLine  audit_line;
   AuditMesg  * new_adt_mesg;
   int  frag_size, max_fragid, reject_ct = 0;
   int64  olaps_kept = 0, olaps_removed = 0, olap_ct = 0;
   char  label_line [1000];
   int  use_binary_output = TRUE, show_progress = FALSE;
   int  ch, error_flag, len, i, j;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "pP")) != EOF))
     switch  (ch)
       {
        case  'p' :
          show_progress = TRUE;
          break;
        case  'P' :
          use_binary_output = FALSE;
          break;
        case  '?' :
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
        default :
          error_flag ++;
        }

   if  (optind >= argc)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-P] <ovl-file> <frag-id-file>\n", 
                 argv [0]);
        fprintf (stderr,
                 "  Options:  -P  generate proto-format (ASCII) output\n");
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
   strcpy (outfile_name + j, "filt.");
   strcat (outfile_name, infile_name + i + 1);
   fprintf (stderr, "Output Overlap File = %s\n", outfile_name);

   ovlfile = File_Open (infile_name, "r");
   outfile = File_Open (outfile_name, "w");

   read_msg_fn = InputFileType_AS (ovlfile);
   if  (use_binary_output)
       write_msg_fn = OutputFileType_AS (AS_BINARY_OUTPUT);
     else
       write_msg_fn = OutputFileType_AS (AS_PROTO_OUTPUT);


   // Get fragment count and allocate array to hold adjacency lists

   frag_size = 1000000;
   Frag = (Frag_Entry_t *) Safe_calloc (frag_size, sizeof (Frag_Entry_t));

   if  (show_progress)
       fprintf (stderr, "\n");

   max_fragid = 0;
   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_OFG :
          {
           OFGMesg  * ofg_mesg = (OFGMesg*) gmesg -> m;
           int  id = ofg_mesg -> iaccession;
           int  i;

           if  (id >= frag_size)
               {
                int  new_size = 2 * frag_size;

                Frag = (Frag_Entry_t *) Safe_realloc
                           (Frag, new_size * sizeof (Frag_Entry_t));
                for  (i = frag_size;  i < new_size;  i ++)
                  memset (Frag + i, 0, sizeof (Frag_Entry_t));
                frag_size = new_size;
               }
           if  (id > max_fragid)
               max_fragid = id;

           Frag [id] . len = ofg_mesg -> clear_rng . end - ofg_mesg -> clear_rng . bgn;

           if  (show_progress && max_fragid % 10000 == 0)
               fprintf (stderr, "\r Frags read = %d", max_fragid);

           break;
          }

        default :
          break;
       }

   frag_size = 1 + max_fragid;
   Frag = (Frag_Entry_t *) Safe_realloc (Frag, frag_size * sizeof (Frag_Entry_t));

   if  (show_progress)
       fprintf (stderr, "\r Frags read = %d\n", frag_size);


   // Add overlaps to adjacency lists

   rewind (ovlfile);

   olap_ct = 0;
   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_OVL :
          {
           int  a_frag, b_frag;
           char  a_end, b_end, orient;
           OverlapMesg  * ovl_mesg = (OverlapMesg*) gmesg -> m;

           a_frag = ovl_mesg -> aifrag;
           b_frag = ovl_mesg -> bifrag;

           switch  (ovl_mesg -> orientation)
             {
              case  AS_NORMAL :
                a_end = 'B';
                b_end = 'A';
                orient = 'S';    // same
                break;
              case  AS_INNIE :
                a_end = 'B';
                b_end = 'B';
                orient = 'F';    // flipped
                break;
              case  AS_OUTTIE :
                a_end = 'A';
                b_end = 'A';
                orient = 'F';    // flipped
                break;
              case  AS_ANTI :
                a_end = 'A';
                b_end = 'B';
                orient = 'S';    // same
                break;
              case  AS_UNKNOWN :
                continue;
              default :
                fprintf (stderr, "ERROR:  Unknown overlap orientation\n");
                assert (FALSE);
             }

           if  (ovl_mesg -> bhg >= 0)
               Add_Entry (Frag + a_frag, a_end, orient,
                          Frag [a_frag] . len - ovl_mesg -> ahg, b_frag);
           Add_Entry (Frag + b_frag, b_end, orient,
                      Frag [b_frag] . len - max (0, ovl_mesg -> bhg), a_frag);

           olap_ct ++;
           if  (show_progress && olap_ct % 10000 == 0)
               fprintf (stderr, "\r Olaps read = " F_S64, olap_ct);

           break;
          }

        default :
          break;
       }

   if  (show_progress)
       fprintf (stderr, "\r Olaps read = " F_S64 "\n", olap_ct);


   printf ("\n%7s %4s  %6s %6s %7s   %6s %6s %7s\n",
           "Frag", "Len", "A max", "A degr", "Score", "B min", "B degr", "Score");
   for  (i = 0;  i < frag_size;  i ++)
     {
      int  a_max, a_degree, b_max, b_degree, b_min;
      int  spanned, hang;
      double  a_score, b_score;
      int  j;

      a_max = b_max = 0;
      a_degree = b_degree = 0;

      for  (j = 0;  j < Frag [i] . num_olaps;  j ++)
        if  (Frag [i] . olap [j] . end == 'A')
            {
             if  (Frag [i] . olap [j] . len > a_max)
                 a_max = Frag [i] . olap [j] . len;
             a_degree ++;
            }
          else
            {
             if  (Frag [i] . olap [j] . len > b_max)
                 b_max = Frag [i] . olap [j] . len;
             b_degree ++;
            }

      if  (a_degree == 0)
          a_score = 0.0;
        else
          a_score = log10 (a_max / ((double) (Frag [i] . len)))
                      * a_degree;
      if  (b_degree == 0)
          b_score = 0.0;
        else
          b_score = log10 (b_max / ((double) (Frag [i] . len)))
                      * b_degree;
      b_min = Frag [i] . len - b_max;

      spanned = FALSE;
      if  (a_degree == 0 && b_degree >= MIN_BP_OLAPS && b_max < Frag [i] . len)
          {
           printf ("%7d %4d  %6s %6s %7s   %6d %6d %7.2f",
                   i, Frag [i] . len, "", "zero", "",
                   b_min, b_degree, b_score);

           for  (j = 0;  j < Frag [i] . num_olaps && ! spanned;  j ++)
             if  (Frag [i] . olap [j] . orient == 'S')
                 {
                  int  id = Frag [i] . olap [j] . frag_id;
                  int  k;

                  hang = Frag [i] . len - Frag [i] . olap [j] . len;
                  for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                    if  (Frag [id] . olap [k] . end == 'A'
                           && Frag [id] . olap [k] . frag_id != i
                           && Frag [Frag [id] . olap [k] . frag_id] . len
                                - Frag [id] . olap [k] . len >= hang)
                        {
                         spanned = TRUE;
#if  SHOW_SPAN
printf ("\n olap1 = %d %c %c %d  olap2 = %d %c %c %d",
        Frag [i] . olap [j] . frag_id,
        Frag [i] . olap [j] . end,
        Frag [i] . olap [j] . orient,
        Frag [i] . olap [j] . len,
        Frag [id] . olap [k] . frag_id,
        Frag [id] . olap [k] . end,
        Frag [id] . olap [k] . orient,
        Frag [id] . olap [k] . len);
#endif
                        }
                 }
               else
                 {
                  int  id = Frag [i] . olap [j] . frag_id;
                  int  k;

                  hang = Frag [i] . len - Frag [i] . olap [j] . len;
                  for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                    if  (Frag [id] . olap [k] . end == 'B'
                           && Frag [id] . olap [k] . frag_id != i
                           && Frag [Frag [id] . olap [k] . frag_id] . len
                                - Frag [id] . olap [k] . len >= hang)
                        {
                         spanned = TRUE;
#if  SHOW_SPAN
printf ("\n olap1 = %d %c %c %d  olap2 = %d %c %c %d",
        Frag [i] . olap [j] . frag_id,
        Frag [i] . olap [j] . end,
        Frag [i] . olap [j] . orient,
        Frag [i] . olap [j] . len,
        Frag [id] . olap [k] . frag_id,
        Frag [id] . olap [k] . end,
        Frag [id] . olap [k] . orient,
        Frag [id] . olap [k] . len);
#endif
                        }
                 }

           printf ("  %-7s\n", spanned ? "spanned" : "");
          }
      else if  (b_degree == 0 && a_degree >= 5 && a_max < Frag [i] . len)
          {
           printf ("%7d %4d  %6d %6d %7.2f   %6s %6s %7s",
                   i, Frag [i] . len, a_max, a_degree, a_score,
                   "", "zero", "");

           spanned = FALSE;
           for  (j = 0;  j < Frag [i] . num_olaps && ! spanned;  j ++)
             if  (Frag [i] . olap [j] . orient == 'S')
                 {
                  int  id = Frag [i] . olap [j] . frag_id;
                  int  k;

                  hang = Frag [i] . len - Frag [i] . olap [j] . len;
                  for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                    if  (Frag [id] . olap [k] . end == 'B'
                           && Frag [id] . olap [k] . frag_id != i
                           && Frag [Frag [id] . olap [k] . frag_id] . len
                                - Frag [id] . olap [k] . len >= hang)
                        {
                         spanned = TRUE;
#if  SHOW_SPAN
printf ("\n olap1 = %d %c %c %d  olap2 = %d %c %c %d",
        Frag [i] . olap [j] . frag_id,
        Frag [i] . olap [j] . end,
        Frag [i] . olap [j] . orient,
        Frag [i] . olap [j] . len,
        Frag [id] . olap [k] . frag_id,
        Frag [id] . olap [k] . end,
        Frag [id] . olap [k] . orient,
        Frag [id] . olap [k] . len);
#endif
                        }
                 }
               else
                 {
                  int  id = Frag [i] . olap [j] . frag_id;
                  int  k;

                  hang = Frag [i] . len - Frag [i] . olap [j] . len;
                  for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                    if  (Frag [id] . olap [k] . end == 'A'
                           && Frag [id] . olap [k] . frag_id != i
                           && Frag [Frag [id] . olap [k] . frag_id] . len
                                - Frag [id] . olap [k] . len >= hang)
                        {
                         spanned = TRUE;
#if  SHOW_SPAN
printf ("\n olap1 = %d %c %c %d  olap2 = %d %c %c %d",
        Frag [i] . olap [j] . frag_id,
        Frag [i] . olap [j] . end,
        Frag [i] . olap [j] . orient,
        Frag [i] . olap [j] . len,
        Frag [id] . olap [k] . frag_id,
        Frag [id] . olap [k] . end,
        Frag [id] . olap [k] . orient,
        Frag [id] . olap [k] . len);
#endif
                        }
                 }

           printf ("  %-7s\n", spanned ? "spanned" : "");
          }
      else if  (a_score <= SCORE_THRESHOLD && b_score <= SCORE_THRESHOLD
                  && a_max - b_min <= OLAP_SLUSH)
          {
           printf ("%7d %4d  %6d %6d %7.2f   %6d %6d %7.2f",
                   i, Frag [i] . len, a_max, a_degree, a_score,
                   b_min, b_degree, b_score);

           spanned = FALSE;
           for  (j = 0;  j < Frag [i] . num_olaps && ! spanned;  j ++)
             if  (Frag [i] . olap [j] . end == 'A')
                 {
                  int  id = Frag [i] . olap [j] . frag_id;
                  int  k;

                  if  (Frag [i] . olap [j] . orient == 'S')
                      {
                       for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                         if  (Frag [id] . olap [k] . end == 'A')
                             {
                              int  m, find_id;

                              find_id = Frag [id] . olap [k] . frag_id;
                              for  (m = 0;  m < Frag [i] . num_olaps && ! spanned;
                                      m ++)
                                if  (Frag [i] . olap [m] . frag_id == find_id
                                       && Frag [i] . olap [m] . end == 'B')
                                    spanned = TRUE;
                             }
                      }
                    else
                      {
                       for  (k = 0;  k < Frag [id] . num_olaps && ! spanned;  k ++)
                         if  (Frag [id] . olap [k] . end == 'B')
                             {
                              int  m, find_id;

                              find_id = Frag [id] . olap [k] . frag_id;
                              for  (m = 0;  m < Frag [i] . num_olaps && ! spanned;
                                      m ++)
                                if  (Frag [i] . olap [m] . frag_id == find_id
                                       && Frag [i] . olap [m] . end == 'B')
                                    spanned = TRUE;
                             }
                      }
                 }

           printf ("  %-7s\n", spanned ? "spanned" : "");
          }
      if  (spanned)
          {
           Frag [i] . reject = TRUE;
           reject_ct ++;
          }

      if  (show_progress && i % 10000 == 0)
          fprintf (stderr, "\r Frags processed = %d", i);
     }

   if  (show_progress)
       fprintf (stderr, "\r Frags processed = %d\n", frag_size);


   rewind (ovlfile);

   pmesg = (GenericMesg *) Safe_malloc (sizeof (GenericMesg));
   pmesg -> t = MESG_ADT;
   pmesg -> m = (AuditMesg *) Safe_malloc (sizeof (AuditMesg));
   new_adt_mesg = (AuditMesg*) pmesg -> m;
   new_adt_mesg -> list = & audit_line;

   fprintf (stderr, "\n Fragments rejected = %9d\n", reject_ct);


   olap_ct = 0;
   while  (read_msg_fn (ovlfile, & gmesg) != EOF && gmesg != NULL)
     switch  (gmesg -> t)
       {
        case  MESG_ADT :
          {
           AuditMesg  * adt_mesg = (AuditMesg*) gmesg -> m;

           sprintf (label_line, "%s %s %s", argv [0], infile_name,
                    argv [optind]);
           AppendAuditLine_AS (adt_mesg, & audit_line, time (0), "frag-anomaly screen",
                               "$Revision: 1.3 $", label_line);
           write_msg_fn (outfile, gmesg);
           break;
          }

        case  MESG_OVL :
          {
           OverlapMesg  * ovl_mesg = (OverlapMesg*) gmesg -> m;

           if  (! Frag [ovl_mesg -> aifrag] . reject
                  && ! Frag [ovl_mesg -> bifrag] . reject)
               {
                write_msg_fn (outfile, gmesg);
                olaps_kept ++;
               }
             else
               olaps_removed ++;

           olap_ct ++;
           if  (show_progress && olap_ct % 10000 == 0)
               fprintf (stderr, "\r Olaps output = " F_S64 , olap_ct);

           break;
          }

        default :
          write_msg_fn (outfile, gmesg);
       }

   if  (show_progress)
       fprintf (stderr, "\r Olaps output = " F_S64 "\n", olap_ct);

   fprintf (stderr, "   Overlaps removed = %9" F_S64P "\n", olaps_removed);
   fprintf (stderr, "      Overlaps kept = %9" F_S64P "\n", olaps_kept);

   fclose (ovlfile);
   fclose (outfile);

   return  0;
  }



static void  Add_Entry
    (Frag_Entry_t * f, char end, char orient, int len, int frag_id)

//  Add entry to  (* f)  indicating an overlap over  end  of this
//  fragment, with orientation  orient , of length  len , to other
//  fragment  frag_id .

  {
   f -> olap = (Olap_Entry_t *) Safe_realloc (f -> olap,
                 1 + f -> num_olaps * sizeof (Olap_Entry_t));
   f -> olap [f -> num_olaps] . end = end;
   f -> olap [f -> num_olaps] . orient = orient;
   f -> olap [f -> num_olaps] . len = len;
   f -> olap [f -> num_olaps] . frag_id = frag_id;

   f -> num_olaps ++;

   return;
  }



