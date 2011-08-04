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

#ifndef ASSERT_H
#define ASSERT_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>

using namespace std;

#include "AS_global.h"
#include "AssertionException.h"
#include "AssertionType.h"

static const char* RCSID_ASSERT_H = "$Id: Assert.h,v 1.1 2011-08-04 14:34:41 mkotelbajcvi Exp $";

#ifndef __ASSERT_FUNCTION
#define __ASSERT_FUNCTION __func__
#endif

class Assert
{
public:
	inline static void assertFalse(bool condition, const char* message = NULL)
	{
		if (condition)
		{
			throw AssertionException(ASSERT_FALSE, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertTrue(bool condition, const char* message = NULL)
	{
		if (!condition)
		{
			throw AssertionException(ASSERT_TRUE, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertNull(void* obj, const char* message = NULL)
	{
		if (obj != NULL)
		{
			throw AssertionException(ASSERT_NULL, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertNotNull(void* obj, const char* message = NULL)
	{
		if (obj == NULL)
		{
			throw AssertionException(ASSERT_NOT_NULL, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertEmpty(const char* str, const char* message = NULL)
	{
		if (str != NULL)
		{
			throw AssertionException(ASSERT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertEmpty(string str, const char* message = NULL)
	{
		if (!str.empty())
		{
			throw AssertionException(ASSERT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertNotEmpty(const char* str, const char* message = NULL)
	{
		if (str == NULL)
		{
			throw AssertionException(ASSERT_NOT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertNotEmpty(string str, const char* message = NULL)
	{
		if (str.empty())
		{
			throw AssertionException(ASSERT_NOT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertEquals(void* obj1, void* obj2, const char* message = NULL)
	{
		if (obj1 != obj2)
		{
			throw AssertionException(ASSERT_EQUALS, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
	
	inline static void assertNotEquals(void* obj1, void* obj2, const char* message = NULL)
	{
		if (obj1 == obj2)
		{
			throw AssertionException(ASSERT_NOT_EQUALS, message, __FILE__, __LINE__, __ASSERT_FUNCTION);
		}
	}
};

#endif
