/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005, J. Craig Venter Institute. All rights reserved.
 * Author: Brian Walenz
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

#ifndef EXCEPTIONUTILS_H
#define EXCEPTIONUTILS_H

static const char* rcsid_EXCEPTIONUTILS_H = "$Id: ExceptionUtils.h,v 1.3 2011-08-04 18:18:56 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <string>

using namespace std;

#include "AS_global.h"

#if defined(__GLIBC__)

#include <execinfo.h>

#else

static  int    backtrace(void **buffer, int size) { return(0); }
static  char **backtrace_symbols(void *const *buffer, int size) { return(NULL); }

#endif

#define DEFAULT_STACK_TRACE_DEPTH 15
#define DEFAULT_STACK_TRACE_CALLER "ExceptionUtils"

typedef struct
{
	size depth;
	char** lines;
} StackTrace;

class ExceptionUtils
{
public:
	static StackTrace& getStackTrace(const char* caller = DEFAULT_STACK_TRACE_CALLER, size depth = DEFAULT_STACK_TRACE_DEPTH);
};

#endif
