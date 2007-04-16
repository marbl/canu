
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
//  fgbovlhits.c
//
//  Read output of  ls -l  *.ovl  and check ranges in filenames to make
//  sure all are present.  Also check file size versus "area" to be
//  overlapped and report anomalies.  Then count how many overlaps
//  (in the format produced by  get-olaps ) hit each file.


#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <errno.h>
#include  <fcntl.h>
#include  <string.h>
#include  <unistd.h>

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

typedef  struct Node
  {
   int  rlo, rhi, hlo, hhi;
   long int  filesize;
   double  area, ratio;
   long int  ovl_ct;
  }  Node_t;


static int  By_rlo
    (const void * a, const void * b)

//  Compare the  rlo  fields in  a  and  b  as  (Node_t *) 's and
//  return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Node_t  * x, * y;

   x = (Node_t *) a;
   y = (Node_t *) b;

   if  (x -> rlo < y -> rlo)
       return  -1;
   else if  (x -> rlo > y -> rlo)
       return  1;
   else if  (x -> rhi < y -> rhi)
       return  -1;
   else if  (x -> rhi > y -> rhi)
       return  1;
   else if  (x -> hlo < y -> hlo)
       return  -1;
   else if  (x -> hlo > y -> hlo)
       return  1;
   else if  (x -> hhi < y -> hhi)
       return  -1;
   else if  (x -> hhi > y -> hhi)
       return  1;

   return  0;
  }



static int  By_hlo
    (const void * a, const void * b)

//  Compare the  rlo  fields in  a  and  b  as  (Node_t *) 's and
//  return  -1  if  a < b ,  0  if  a == b , and  1  if  a > b .
//  Used for  qsort .

  {
   Node_t  * x, * y;

   x = (Node_t *) a;
   y = (Node_t *) b;

   if  (x -> hlo < y -> hlo)
       return  -1;
   else if  (x -> hlo > y -> hlo)
       return  1;
   else if  (x -> hhi < y -> hhi)
       return  -1;
   else if  (x -> hhi > y -> hhi)
       return  1;
   else if  (x -> rlo < y -> rlo)
       return  -1;
   else if  (x -> rlo > y -> rlo)
       return  1;
   else if  (x -> rhi < y -> rhi)
       return  -1;
   else if  (x -> rhi > y -> rhi)
       return  1;

   return  0;
  }



int  main
    (int argc, char * argv [])

  {
   FILE  * olap_file;
   Node_t  * list = NULL;
   int  list_size, ct = 0;
   int  basic_version = FALSE;
   char  line [10000];
   int  ch, error_flag;
   int  old_rlo, old_rhi, old_hhi, old_max;
   int  i;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "b")) != EOF))
     switch  (ch)
       {
        case  'b' :
          basic_version = TRUE;
          break;
        case  'h' :
          error_flag = TRUE;
          break;
        default :
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
          error_flag = TRUE;
        }

   if  (optind != argc - 1 || error_flag)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-b] <olapfile>\n"
                 "  Read output of  ls -l  of  .ovl  files (from stdin)\n"
                 "  with canonical names and verify coverage of fragment\n"
                 "  ranges and give ratio of  filesize/fragmentarea\n"
                 "  Also count number of olaps in  <olapfile>  that hit\n"
                 "  each file.  Output (to stdout) is same as  olapcoverage\n"
                 "  program followed by hit counts.  Format of  <olapfile>\n"
                 "  just needs one olap per line, with iids as first 2 entries\n"
                 "  separated by spaces.  This is produced by  get-olaps \n"
                 "  for example.\n",
                 argv [0]);
        fprintf (stderr,
                 "  Options:  -b  only use canonical name at end of line\n"
                 "                no filesize info used\n");
        exit (EXIT_FAILURE);
       }

   olap_file = fopen (argv [optind], "r");
   if  (olap_file == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file \"%s\"\n",
                 argv [optind]);
        exit (EXIT_FAILURE);
       }

   list_size = 1000;
   list = (Node_t *) malloc (list_size * sizeof (Node_t));
   assert (list != NULL);

   while  (fgets (line, 10000, stdin) != NULL)
     {
      char  permission [100];
      char  number [100];
      char  owner [100];
      char  group [100];
      char  month [100];
      char  day [100];
      char  time [100];
      char  name [5000];
      long int  size;
      char  * p;

      if  (basic_version)
          {
           if  (ct >= list_size - 1)
               {
                list_size += 1000;
                list = (Node_t *) realloc (list, list_size * sizeof (Node_t));
                assert (list != NULL);
               }
           list [ct] . filesize = 0;
           list [ct] . ovl_ct = 0;
           p = strrchr (line, '.');
           if (p == NULL)
              {
               fprintf (stderr, "bad filename:  %s", line);
              }
           for  (p --;  p >= line && * p != '.';  p --)
             ;
           if  (sscanf (p + 1, "r%d-%dh%d-%d.ovl",
                        & list [ct] . rlo, & list [ct] . rhi,
                        & list [ct] . hlo, & list [ct] . hhi) == 4)
               ct ++;
          }
        else
          {
           if  (sscanf (line, "%s %s %s %s %ld %s %s %s %s",
                        permission, number, owner, group, & size, month, day,
                        time, name)
                  == 9)
               {
                if  (ct >= list_size - 1)
                    {
                     list_size += 1000;
                     list = (Node_t *) realloc (list, list_size * sizeof (Node_t));
                     assert (list != NULL);
                    }
                list [ct] . filesize = size;
                list [ct] . ovl_ct = 0;
                p = strrchr (name, '.');
                if (p == NULL)
                   {
                    fprintf (stderr, "bad filename:  %s", line);
                   }
                for  (p --;  p >= name && * p != '.';  p --)
                  ;
                if  (sscanf (p + 1, "r%d-%dh%d-%d.ovl",
                             & list [ct] . rlo, & list [ct] . rhi,
                             & list [ct] . hlo, & list [ct] . hhi) == 4)
                    ct ++;
               }
          }
     }

   qsort (list, ct, sizeof (Node_t), By_rlo);

   old_rlo = -1;
   old_rhi = 0;
   old_hhi = 0;
   old_max = 0;
   for  (i = 0;  i < ct;  i ++)
     {
      int  effective_hlo;

      if  (list [i] . hlo < list [i] . rlo)
          effective_hlo = list [i] . rlo;
        else
          effective_hlo = list [i] . hlo;

      list [i] . area = (double) (1 + list [i] . hhi - effective_hlo)
                          * (1 + list [i] . rhi - list [i] . rlo);

      if  (effective_hlo <= list [i] . rhi)
          {
           double  olap;

           olap = 1 + list [i] . rhi - effective_hlo;
           list [i] . area -= (olap * (olap + 1)) / 2;
          }

      list [i] . ratio = list [i] . filesize / list [i] . area;

      printf ("%9d - %9d  %9d - %9d  %9ld  %5.1f  %5.1f\n",
              list [i] . rlo, list [i] . rhi,
              list [i] . hlo, list [i] . hhi,
              list [i] . filesize,
              list [i] . area * 1e-9,
              list [i] . ratio * 1e6);

#if  1
      if  (list [i] . rlo != old_rlo)
          {
           if  (list [i] . rlo != old_rhi + 1)
               {
                printf ("ERROR:  Gap in r coverage\n");
                exit (-1);
               }
           if  (list [i] . hlo != list [i] . rlo)
               {
                printf ("WARNING:  hlo should be rlo\n");
               }
           if  (old_hhi != old_max)
               {
                printf ("ERROR:  Missing hi coverage\n");
                exit (-1);
               }

           old_rlo = list [i] . rlo;
           old_rhi = list [i] . rhi;
           old_hhi = list [i] . hhi;
          }
      else if  (list [i] . rhi != old_rhi)
          {
           if  (list [i] . rhi <= old_hhi)
               {
                printf ("ERROR:  Unexpected jump in rhi\n");
                exit (-1);
               }
           if  (list [i] . hlo != old_hhi + 1)
               {
                printf ("ERROR:  Gap in h coverage\n");
                exit (-1);
               }

           old_rhi = list [i] . rhi;
           old_hhi = list [i] . hhi;
           if  (old_hhi > old_max)
               old_max = old_hhi;
          }
        else
          {
           if  (list [i] . hlo != old_hhi + 1)
               {
                printf ("ERROR:  Gap in h coverage\n");
                exit (-1);
               }

           old_hhi = list [i] . hhi;
           if  (old_hhi > old_max)
               old_max = old_hhi;
          }
#endif
     }

#if  1
   if  (old_hhi != old_max || old_rhi != old_max)
       {
        printf ("ERROR:  Missing hi coverage\n");
        exit (-1);
       }
#endif

//  Now count the number of overlaps in the file specified on the
//  command line that hit each entry in list.

   while  (fgets (line, 10000, olap_file) != NULL)
     {
      int  a, b;

      if  (sscanf (line, "%d %d", & a, & b) != 2)
          fprintf (stderr, "Bad format on line: %s", line);
        else
          {
           int  lo, hi, mid;

           if  (b < a)
               {
                int  save = a;

                a = b;
                b = save;
               }

           lo = 0;
           hi = ct - 1;
           while  (lo <= hi)
             {
              mid = (lo + hi) / 2;
              if  (a < list [mid] . rlo
                     || (a <= list [mid] . rhi && b < list [mid] . hlo))
                  hi = mid - 1;
              else if  (a > list [mid] . rhi
                     || (a >= list [mid] . rlo && b > list [mid] . hhi))
                  lo = mid + 1;
                else
                  {
                   assert (list [mid] . rlo <= a
                             && a <= list [mid] . rhi
                             && list [mid] . hlo <= b
                             && b <= list [mid] . hhi);
                   list [mid] . ovl_ct ++;
                   break;
                  }
             }
           if  (hi < lo)
               {
                fprintf (stderr, "Did not find olap %8d v %8d\n", a, b);
               }
          }
     }

   printf ("\nOlap hits:\n");
   for  (i = 0;  i < ct;  i ++)
     printf ("r %8d-%8d  h %8d - %8d  %8ld  %5.1f  %8.1f %7.1f\n",
             list [i] . rlo, list [i] . rhi,
             list [i] . hlo, list [i] . hhi,
             list [i] . ovl_ct,
             list [i] . area * 1e-9,
             (1e10 * list [i] . ovl_ct) / list [i] . area,
             list [i] . ratio * 1e6);

   return  0;
  }






