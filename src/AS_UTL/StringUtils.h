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

#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <exception>
#include <string>

#include "ArgumentException.h"
#include "AS_global.h"
#include "AS_UTL_alloc.h"

static const char *rcsid_STRING_UTILS_H = "$Id: StringUtils.h,v 1.2 2011-07-28 04:35:33 mkotelbajcvi Exp $";

static const char* WHITESPACE_CHARS = " \t\n\r";

class StringUtils
{
public:
	static bool isBlank(const char* str);
	static bool isEmpty(const char* str);
	static bool isEmpty(std::string str);

	static bool areEqual(const char* str1, const char* str2);

	static const char* trim(const char* str, size num, ...);
	static const char* trim(const char* str, size num, const char** toTrim);
	static const char* trimStart(const char* str, size num, ...);
	static const char* trimStart(const char* str, size num, const char** toTrim);
	static const char* trimEnd(const char* str, size num, ...);
	static const char* trimEnd(const char* str, size num, const char** toTrim);
};

template<typename T>
T* getArgs(size num, T* args, va_list& argsList);

#endif
