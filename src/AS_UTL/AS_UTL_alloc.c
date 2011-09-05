
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

static const char *rcsid = "$Id: AS_UTL_alloc.c,v 1.16 2011-09-05 16:49:44 mkotelbajcvi Exp $";

//  We explicitly do not include AS_UTL_alloc.h here, because it
//  redefines malloc(), calloc(), realloc() and free() to be errors.
//  We want everyone to use the safe_*() versions supplied here.
//
//#include "AS_UTL_alloc.h"

//  We want to include AS_global.h to get the F_SIZE_T definition, but that
//  includes AS_UTL_alloc.h, so we pretend we've already included it.
//
#define AS_UTL_ALLOC_H

#include <assert.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

#include "AS_global.h"


void *
safe_calloc(size_t num, size_t len) {

  if ((num == 0) || (len == 0))
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = calloc(num, len);

  if (p == NULL)
    fprintf(stderr, "Could not calloc memory ("F_SIZE_T" * "F_SIZE_T" bytes = "F_SIZE_T")\n",
            num, len, num*len);
  assert(p != NULL);

  return(p);
}


void *
safe_malloc(size_t len) {

  if (len == 0)
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = malloc(len);

  if (p == NULL)
    fprintf(stderr, "Could not malloc memory ("F_SIZE_T" bytes)\n", len);
  assert(p != NULL);

  return(p);
}


void *
safe_realloc(void *q, size_t len) {

  if (len == 0)
    return(NULL);    //  Bail, user didn't request anything.

  void  *p = realloc(q, len);

  if (p == NULL)
    fprintf(stderr, "Could not realloc memory ("F_SIZE_T" bytes)\n", len);
  assert(p != NULL);

  return(p);
}


void
safe_free2(void *q) {
  free(q);
}
