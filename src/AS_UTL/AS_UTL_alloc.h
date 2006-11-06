
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

//  Memory allocation

void *safe_calloc(size_t num, size_t len);
void *safe_malloc(size_t len);
void *safe_realloc(void *q, size_t len);

#define safe_free(Q) { free(Q); Q = NULL; }

#if 0
#define malloc(X)  ERROR_MALLOC(X)
#define calloc(X)  ERROR_CALLOC(X)
#define realloc(X) ERROR_REALLOC(X)
#define free(X)    ERROR_FREE(X)
#endif

#endif // AS_UTL_ALLOC_H
