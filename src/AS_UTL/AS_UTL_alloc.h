
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

#ifndef AS_UTL_ALLOC_H
#define AS_UTL_ALLOC_H

static const char *rcsid_AS_UTL_ALLOC_H = "$Id: AS_UTL_alloc.h,v 1.7 2008-10-08 22:03:00 brianwalenz Exp $";

//
//  The safe_*alloc routines are the same as the normal routines,
//  except they print a standard message and assert if memory cannot
//  be allocated.
//
//  safe_free is somewhat useless.  It frees the memory, then sets the
//  pointer to NULL.  It's useless because many times the pointer we
//  set to NULL is a local copy, e.g.:
//
//     void freeSomeStructure(SomeStructure *p) {
//       safe_free(p->data);
//       safe_free(p);
//     }
//

void *safe_calloc(size_t num, size_t len);
void *safe_malloc(size_t len);
void *safe_realloc(void *q, size_t len);
void  safe_free2(void *);

#define safe_free(Q) { safe_free2(Q); Q = NULL; }

//  And, thanks, GNU.  strdup() (and, sigh, probably lots others) are
//  implemented as a macro that calls free(), which then gets expanded
//  into our bogus function.
//
#ifndef X86_GCC_LINUX

#define malloc(X)     use_safe_malloc_instead(X)
#define calloc(X,Y)   use_safe_calloc_instead(X,Y)
#define realloc(X,Y)  use_safe_realloc_instead(X,Y)
#define free(X)       use_safe_free_instead(X)

#endif

#endif // AS_UTL_ALLOC_H
