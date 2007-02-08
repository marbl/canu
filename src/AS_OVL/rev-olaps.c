
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
* Module:  rev-olaps.c
* Description:
*   Reads file of overlaps (in the format produced by  get-olaps )
*   and prints each one out in reverse orientation (in the same format).
*   Input is from  stdin .  Output is to  stdout .
* 
*    Programmer:  A. Delcher
*       Started:   5 Dec 2000
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: rev-olaps.c,v 1.5 2007-02-08 02:46:00 brianwalenz Exp $
 * $Revision: 1.5 $
*/

static char CM_ID[] = "$Id: rev-olaps.c,v 1.5 2007-02-08 02:46:00 brianwalenz Exp $";


//  System include files

#include  <stdlib.h>
#include  <stdio.h>
#include  <ctype.h>



int  main
    (int argc, char **argv)

  {
   long int  a_id, b_id;
   int  a_hang, b_hang, hang1 = 0, hang2 = 0;
   char  orient [10];
   double  error_rate;

   while  (scanf ("%ld %ld %d %d %s %lf",
           & a_id, & b_id, & a_hang, & b_hang, orient, & error_rate) == 6)
     {
      orient [0] = toupper (orient [0]);
      switch  (orient [0])
        {
         case  'I' :
         case  'O' :
           hang1 = b_hang;
           hang2 = a_hang;
           break;
           
         case  'N' :
           hang1 = -a_hang;
           hang2 = -b_hang;
           break;
        }
      printf ("%8ld %8ld %5d %5d  %c %5.2f\n",
              b_id, a_id, hang1, hang2, orient [0], error_rate);
     }

   return  0;
  }

