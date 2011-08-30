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

#ifndef VARUTILS_H
#define VARUTILS_H

static const char* rcsid_VARUTILS_H = "$Id: VarUtils.h,v 1.4 2011-08-30 23:09:51 mkotelbajcvi Exp $";

#include <stdarg.h>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>

using namespace std;

#include "AS_global.h"

#define initArgs(lastStaticArg) \
	va_list argsList; \
	va_start(argsList, lastStaticArg)

class VarUtils
{
public:
	template<class T>
	inline static T* getArgs(size_t num, va_list& argsList)
	{
		T* args = new T[num];
		T arg;
		
		for (size_t a = 0; (a < num) && ((arg = va_arg(argsList, T)) != NULL); a++)
		{
			args[a] = arg;
		}
		
		va_end(argsList);
		
		return args;
	}
};

#endif
