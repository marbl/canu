
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
//  olapcoverage.c
//
//  Read output of  ls -l  *.ovl  and check ranges in filenames to make
//  sure all are present.  Also check file size versus "area" to be
//  overlapped and report anomalies.


#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <errno.h>
#include  <fcntl.h>
#include  <sys/types.h>
#include  <string.h>
#include  <dirent.h>
#include  <sys/stat.h>
#include  <unistd.h>

static int  Error_Ct = 0;

#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE  (1)
#endif

struct  Node_t
  {
   int  rlo, rhi, hlo, hhi;
   long int  filesize;
   double  area, ratio;
   char  * message;

   void  Print  (void)
     {
      printf ("%9d - %9d  %9d - %9d  %9ld  %6.1f  %5.1f  %s\n",
              rlo, rhi, hlo, hhi, filesize,
              area * 1e-9, ratio * 1e6, message);
     }

   void  Set_Area_And_Ratio  (void)
     {
      double  base, a, b;
      int  b_ref;

      if  (rlo > hhi)
          {
           area = ratio = 0.0;
           message = "No overlaps possible";
           return;
          }
      if  (rhi <= hhi)
          a = 1 + rhi - rlo;
        else
          {
           a = 1 + hhi - rlo;
           message = "Extraneous range";
          }
      if  (rlo <= hlo)
          {
           base = 1 + hhi - hlo;
           b_ref = hlo;
          }
        else
          {
           base = 1 + hhi - rlo;
           message = "Extraneous range";
           b_ref = rlo;
          }
      if  (b_ref < rhi)
          b = 1 + b_ref - rlo;
        else
          b = 1 + rhi - rlo;

      // Trapezoid with base  base  and sides  a  and  b .
      area = base * (a + b) / 2.0;

      assert (area > 0.0);
      ratio = filesize / area;
     }
  };


struct  Point
  {
   int  x, y;
  };

class  Skyline_t
  {
  private:
   Point  * pt;
   int  n;

  public:
   Skyline_t  (int size)
     {
      assert (size > 0);
      pt = (Point *) malloc (size * sizeof (Point));
      assert (pt != NULL);
      pt [0] . x = 1;
      pt [0] . y = 0;
      n = 1;
     }
   ~ Skyline_t  ()
     {
      free (pt);
     }

   void  Add_Block  (Node_t * p)
   //  Add the region specified in  p  to this skyline
     {
      int  i, j, end, top;

      for  (i = 0;  i < n && p -> rlo > pt [i] . x;  i ++)
        ;

      if  (p -> rlo != pt [i] . x
             || p -> hlo != pt [i] . y + 1)
          {
           printf ("ERROR:  Unexpected new block p rlo=%d hlo=%d pt[%d]x=%d y=%d\n",
                   p->rlo, p->hlo, i, pt[i].x, pt[i].y);
           Error_Ct ++;
           // exit (EXIT_FAILURE);
          }

      if  (i > 0 && pt [i - 1] . y == p -> hhi)   // Extends previous entry
          {
           if  (i == n - 1)   // Last entry
               {
                if  (p -> rhi < p -> hhi)
                    pt [i] . y = p -> rhi;
                  else
                    pt [i] . y = p -> hhi;
                pt [i] . x = pt [i] . y + 1;
                return;
               }
           if  (i == n - 2)   // Next to last entry
               {
                assert (pt [i] . y == pt [i + 1] . y);
                if  (p -> rhi < pt [i + 1] . x - 1)
                    pt [i] . x = p -> rhi + 1;
                else if  (p -> rhi == pt [i + 1] . x - 1)   // Exact fit
                    {
                     pt [i] = pt [i + 1];
                     n --;
                    }
                  else
                    {
                     if  (p -> rhi > p -> hhi)
                         end = p -> hhi;
                       else
                         end = p -> rhi;
                     pt [i] . y = end;
                     pt [i] . x = end + 1;
                     n --;
                    }
                return;
               }
           assert (pt [i] . y < pt [i + 1] . y);
           if  (p -> rhi >= pt [i + 1] . x)   // Duplicate coverage
               {
                if  (p -> rhi > pt [i + 2] . x)
                    end = pt [i + 2] . x;
                  else
                    end = p -> rhi;
                if  (p -> hhi <= pt [i + 1] . y)
                    top = p -> hhi;
                  else
                    top = pt [i + 1] . y;
                printf ("ERROR:  Duplicate coverage of r%d-%d vs h%d-%d\n",
                        pt [i + 1] . x, end, p -> hlo, top);
                Error_Ct ++;
                exit (EXIT_FAILURE);
               }
           if  (p -> rhi < pt [i + 1] . x - 1)  // Shrinks this entry
               {
                pt [i] . x = p -> rhi + 1;
                return;
               }
           if  (p -> hhi == pt [i + 1] . y)   // Combine i-1, i & i+1
               {
                for  (j = i + 2;  j < n;  j ++)
                  pt [j - 2] = pt [j];
                n -= 2;
                return;
               }
           // Combine i-1 & i
           for  (j = i + 1;  j < n;  j ++)
             pt [j - 1] = pt [j];
           n --;
           return;
          }

      // Can ignore previous entry if there is one
      if  (i == n - 1)   // Last entry
          {
           pt [i] . y = p -> hhi;
           if  (p -> rhi < p -> hhi)
               pt [i + 1] . y = p -> rhi;
             else
               pt [i + 1] . y = p -> hhi;
           pt [i + 1] . x = pt [i] . y + 1;
           n ++;
           return;
          }
      if  (i == n - 2)   // Next to last entry
          {
           assert (pt [i] . y == pt [i + 1] . y);
           if  (p -> rhi >= pt [i + 1] . x)   // Extends past previous last
               {
                pt [i] . y = p -> hhi;
                if  (p -> rhi > p -> hhi)
                    end = p -> hhi;
                  else
                    end = p -> rhi;
                pt [i + 1] . y = end;
                pt [i + 1] . x = end + 1;
               }
           else if  (p -> rhi == pt [i + 1] . x - 1)   // Exact fit
               {
                pt [i] . y = p -> hhi;
               }
             else   // Split this entry
               {
                for  (j = n;  j > i;  j --)
                  pt [j] = pt [j - 1];
                n ++;
                pt [i + 1] . x = p -> rhi + 1;
                pt [i] . y = p -> hhi;
               }
           return;
          }
      assert (pt [i] . y < pt [i + 1] . y);
      if  (p -> rhi >= pt [i + 1] . x)   // Duplicate coverage
          {
           if  (p -> rhi > pt [i + 2] . x)
               end = pt [i + 2] . x;
             else
               end = p -> rhi;
           if  (p -> hhi <= pt [i + 1] . y)
               top = p -> hhi;
             else
               top = pt [i + 1] . y;
           printf ("ERROR:  Duplicate coverage of r%d-%d vs h%d-%d\n",
                   pt [i + 1] . x, end, p -> hlo, top);
           Error_Ct ++;
           exit (EXIT_FAILURE);
          }
      if  (p -> rhi < pt [i + 1] . x - 1)  // Split this entry
          {
           for  (j = n;  j > i;  j --)
             pt [j] = pt [j - 1];
           n ++;
           pt [i + 1] . x = p -> rhi + 1;
           pt [i] . y = p -> hhi;
           return;
          }
      if  (p -> hhi == pt [i + 1] . y)   // Combine i & i+1
          {
           for  (j = i + 2;  j < n;  j ++)
             pt [j - 1] = pt [j];
           n --;
           pt [i] . y = p -> hhi;
           return;
          }
      // Exact fit
      pt [i] . y = p -> hhi;

      return;
     }

   void  Check_Unfilled  (int x, int y)
   // Print unfilled areas below block with lower left corner
   // at  (x, y)  and then modify skyline to fill them in.

     {
      int  i, j, next;

      for  (i = 0;  i < n;  i ++)
        if  (pt [i] . x < x && pt [i] . y < y)
            {
             if  (i < n - 1)
                 next = pt [i + 1] . x;
             else if  (x <= y)
                 next = y;
               else
                 next = y + 1;
             if  (next <= x)
                 {
                  printf ("ERROR:  Uncovered range r%d-%d vs h%d-%d\n",
                          pt [i] . x, next - 1, pt [i] . y + 1, y);
                  Error_Ct ++;
                  pt [i] . y = y;
                 }
               else
                 {
                  printf ("ERROR:  Uncovered range r%d-%d vs h%d-%d\n",
                          pt [i] . x, x - 1, pt [i] . y + 1, y);
                  Error_Ct ++;
                  if  (pt [i] . y < y - 1)
                      {
                       printf ("ERROR:  Uncovered range r%d-%d vs h%d-%d\n",
                               x, next - 1, pt [i] . y + 1, y - 1);
                       Error_Ct ++;
                      }
                  for  (j = n;  j > i;  j --)
                    pt [j] = pt [j - 1];
                  pt [i] . y = y;
                  pt [i + 1] . x = x;
                  pt [i + 1] . y = y - 1;
                  i ++;
                  n ++;
                 }
             if  (i == n - 1)
                 {
                  pt [n] . x = next;
                  pt [n] . y = (next == y) ? y - 1 : y;
                  i ++;
                  n ++;
                 }
            }
        else if  (pt [i] . y < y - 1)
            {
             if  (i < n - 1)
                 next = pt [i + 1] . x;
               else
                 next = y;
             printf ("ERROR:  Uncovered range r%d-%d vs h%d-%d\n",
                     pt [i] . x, next - 1, pt [i] . y + 1, y - 1);
             Error_Ct ++;
             pt [i] . y = y - 1;
             if  (i == n - 1)
                 {
                  pt [n] . x = next;
                  pt [n] . y = y - 1;
                  i ++;
                  n ++;
                 }
            }

      // Merge duplicates
      j = 0;
      for  (i = 1;  i < n - 1;  i ++)
        if  (pt [i] . y != pt [j] . y)
            pt [++ j] = pt [i];
      if  (n > 1)
          pt [++ j] = pt [n - 1];
      n = j + 1;

      return;
     }

   int  Get_Top  (void)
   // Return the highest y coordinate in the skyline

     {
      int  i, top;

      top = 0;
      for  (i = 1;  i < n;  i ++)
        if  (pt [i] . y > pt [top] . y)
            top = i;

      return  pt [top] . y;
     }

   Point  Lowest_Left  (void)
     {
      Point  result;
      int  i, lo = 0;

      assert (n > 0);
      for  (i = 1;  i < n;  i ++)
        if  (pt [i] . y < pt [lo] . y
               || (pt [i] . y == pt [lo] . y && pt [i] . x < pt [lo] . x))
            lo = i;
      result . x = pt [lo] . x;
      result . y = pt [lo] . y + 1;

      return  result;
     }

   void  Print  (void)
     {
      int  i;
      printf ("Skyline with n = %d:\n", n);

      for  (i = 0;  i < n;  i ++)
        printf ("  (%8d %8d)\n", pt [i] . x, pt [i] . y);

      return;
     }

  };



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



static int  By_hlo_rlo
    (const void * a, const void * b)

//  Compare the  hlo  fields in  a  and  b  as  (Node_t *) 's and
//  return  -1  if  a < b  and  1  if  a > b .  If equal then
//  break tie using  rlo .  Used for  qsort .

  {
   Node_t  * x, * y;

   x = (Node_t *) a;
   y = (Node_t *) b;

   if  (x -> hlo < y -> hlo)
       return  -1;
   else if  (x -> hlo > y -> hlo)
       return  1;
   else if  (x -> rlo < y -> rlo)
       return  -1;
   else if  (x -> rlo > y -> rlo)
       return  1;

   return  0;
  }



int  main
    (int argc, char * argv [])

  {
   Node_t  * list = NULL;
   int  list_size, ct = 0;
   int  basic_version = FALSE;
   char  line [10000];
   int  ch, error_flag;
   int  old_rlo, old_rhi, old_hhi, old_max;
   int  i, top;

   optarg = NULL;
   error_flag = FALSE;
   while  (! error_flag && ((ch = getopt (argc, argv, "b")) != EOF))
     switch  (ch)
       {
        case  'b' :
          basic_version = TRUE;
          break;
        default :
          fprintf (stderr, "Unrecognized option \"-%c\"\n", optopt);
          error_flag = TRUE;
        }

   if  (optind < argc || error_flag)
       {
        fprintf (stderr, 
                 "USAGE:  %s [-b]\n"
                 "  Read output of  ls -l  of  .ovl  files with canonical names\n"
                 "  and verify coverage of fragment ranges and give ratio\n"
                 "  of  filesize/fragmentarea\n",
                 argv [0]);
        fprintf (stderr,
                 "  Options:  -b  only use canonical name at end of line\n"
                 "                no filesize info used\n");
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

      if  (ct >= list_size - 1)
          {
           list_size += 1000;
           list = (Node_t *) realloc (list, list_size * sizeof (Node_t));
           assert (list != NULL);
          }
      if  (basic_version)
          {
           list [ct] . filesize = 0;
           list [ct] . message = "";
           p = strrchr (line, '.');
           assert (p != NULL);
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
                list [ct] . filesize = size;
                list [ct] . message = "";
                p = strrchr (name, '.');
                assert (p != NULL);
                for  (p --;  p >= name && * p != '.';  p --)
                  ;
                if  (sscanf (p + 1, "r%d-%dh%d-%d.ovl",
                             & list [ct] . rlo, & list [ct] . rhi,
                             & list [ct] . hlo, & list [ct] . hhi) == 4)
                    ct ++;
               }
          }
     }

   qsort (list, ct, sizeof (Node_t), By_hlo_rlo);

   Skyline_t  skyline (1 + ct);

   printf ("%21s  %21s  %9s  %6s  %5s  %s\n",
           "Old Frag Range", "New Frag Range", "FileSize", "Area", "Ratio", "Note");
   for  (i = 0;  i < ct;  i ++)
     {
      Point  low_left;

      list [i] . Set_Area_And_Ratio ();

      low_left = skyline . Lowest_Left ();
      if  (list [i] . rlo != low_left . x
             || list [i] . hlo != low_left . y)
          {
           skyline . Check_Unfilled (list [i] . rlo, list [i] . hlo);
          }

      list [i] . Print ();

      skyline . Add_Block (list + i);
     }

   top = skyline . Get_Top ();
   skyline . Check_Unfilled (top + 1, top);
   printf ("%d error%s\n", Error_Ct, Error_Ct == 1 ? "" : "s");

#if  0
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
#endif

   return  0;
  }





