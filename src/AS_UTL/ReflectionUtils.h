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

#ifndef REFLECTIONUTILS_H
#define REFLECTIONUTILS_H

static const char* rcsid_REFLECTIONUTILS_H = "$Id: ReflectionUtils.h,v 1.1 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <typeinfo>

using namespace std;

#if defined(__GLIBCXX__)
#include <cxxabi.h>

using namespace abi;
#endif

#include "AS_global.h"
#include "AS_UTL_alloc.h"

namespace Utility
{
	static const char STACK_ENTRY_MANGLED_NAME_PREFIX = '(';
	static const char STACK_ENTRY_MANGLED_NAME_SUFFIX = '+';
	static const char* STACK_ENTRY_OBJECT_DELIMITER = "::";
	static const char* STACK_ENTRY_FUNCTION_SUFFIX = "()";
	
	static const char* STACK_ENTRY_MAIN_NAME = "main";
	static const char* STACK_ENTRY_LIBC_START_MAIN_NAME = "__libc_start_main";
	static const char* STACK_ENTRY_GXX_PERSONALITY_NAME = "__gxx_personality_v0";
	
	static const size_t INITIAL_DEMANGLING_NAME_BUFFER_SIZE = 128;
	
	class ReflectionUtils
	{
	public:
		static bool demangleStackEntry(string& entry, bool includePath = false);
		static bool demangle(string& name);
		
	private:
		ReflectionUtils()
		{
		}
	};
}

#endif
