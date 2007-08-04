
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

//  $Id: AS_global.c,v 1.2 2007-08-04 22:27:35 brianwalenz Exp $

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//  Nonsense values, mostly for making sure everybody that uses an
//  error rate calls AS_configure() at startup.

double AS_OVL_ERROR_RATE = -90.0;
double AS_CGW_ERROR_RATE = -90.0;
double AS_CNS_ERROR_RATE = -90.0;
double AS_MAX_ERROR_RATE =   0.25;

//  We take argc and argv, so, maybe, eventually, we'll want to parse
//  something out of there.  We return argc in case what we parse we
//  want to remove.
//
int
AS_configure(int argc, char **argv) {
  char *p = NULL;

  AS_OVL_ERROR_RATE = 0.06;
  AS_CGW_ERROR_RATE = 0.10;
  AS_CNS_ERROR_RATE = 0.06;

  p = getenv("AS_OVL_ERROR_RATE");
  if (p) {
    AS_OVL_ERROR_RATE = atof(p);
    if ((AS_OVL_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_OVL_ERROR_RATE)) {
      fprintf(stderr, "%s: ERROR:  Invalid AS_OVL_ERROR_RATE ('%s'); should be between 0.0 and %0.2f\n", argv[0], p, AS_MAX_ERROR_RATE);
      exit(1);
    }
    fprintf(stderr, "%s: AS_configure()-- AS_OVL_ERROR_RATE set to %0.2f\n", argv[0], AS_OVL_ERROR_RATE);
  }

  p = getenv("AS_CGW_ERROR_RATE");
  if (p) {
    AS_CGW_ERROR_RATE = atof(p);
    if ((AS_CGW_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_CGW_ERROR_RATE)) {
      fprintf(stderr, "%s: ERROR:  Invalid AS_CGW_ERROR_RATE ('%s'); should be between 0.0 and %0.2f\n", argv[0], p, AS_MAX_ERROR_RATE);
      exit(1);
    }
    fprintf(stderr, "%s: AS_configure()-- AS_CGW_ERROR_RATE set to %0.2f\n", argv[0], AS_CGW_ERROR_RATE);
  }

  p = getenv("AS_CNS_ERROR_RATE");
  if (p) {
    AS_CNS_ERROR_RATE = atof(p);
    if ((AS_CNS_ERROR_RATE < 0.0) || (AS_MAX_ERROR_RATE < AS_CNS_ERROR_RATE)) {
      fprintf(stderr, "%s: ERROR:  Invalid AS_CNS_ERROR_RATE ('%s'); should be between 0.0 and %0.2f\n", argv[0], p, AS_MAX_ERROR_RATE);
      exit(1);
    }
    fprintf(stderr, "%s: AS_configure()-- AS_CNS_ERROR_RATE set to %0.2f\n", argv[0], AS_CNS_ERROR_RATE);
  }

  return(argc);
}
