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

#ifndef ASSERTS_H
#define ASSERTS_H

static const char* rcsid_ASSERTS_H = "$Id: Asserts.h,v 1.1 2011-08-11 17:34:34 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>

using namespace std;

#include "AS_global.h"
#include "AssertionException.h"
#include "AssertionType.h"

class Asserts
{
public:
	inline static void assertFalse(bool condition, const char* message = NULL)
	{
		if (condition)
		{
			throw AssertionException(message, NULL, ASSERT_FALSE);
		}
	}
	
	inline static void assertTrue(bool condition, const char* message = NULL)
	{
		if (!condition)
		{
			throw AssertionException(message, NULL, ASSERT_TRUE);
		}
	}
	
	inline static void assertNull(void* obj, const char* message = NULL)
	{
		if (obj != NULL)
		{
			throw AssertionException(message, NULL, ASSERT_NULL);
		}
	}
	
	inline static void assertNotNull(void* obj, const char* message = NULL)
	{
		if (obj == NULL)
		{
			throw AssertionException(message, NULL, ASSERT_NOT_NULL);
		}
	}
	
	inline static void assertEmpty(const char* str, const char* message = NULL)
	{
		if (str != NULL)
		{
			throw AssertionException(message, NULL, ASSERT_EMPTY);
		}
	}
	
	inline static void assertEmpty(string str, const char* message = NULL)
	{
		if (!str.empty())
		{
			throw AssertionException(message, NULL, ASSERT_EMPTY);
		}
	}
	
	inline static void assertNotEmpty(const char* str, const char* message = NULL)
	{
		if (str == NULL)
		{
			throw AssertionException(message, NULL, ASSERT_NOT_EMPTY);
		}
	}
	
	inline static void assertNotEmpty(string str, const char* message = NULL)
	{
		if (str.empty())
		{
			throw AssertionException(message, NULL, ASSERT_NOT_EMPTY);
		}
	}
	
	inline static void assertEquals(void* obj1, void* obj2, const char* message = NULL)
	{
		if (obj1 != obj2)
		{
			throw AssertionException(message, NULL, ASSERT_EQUALS);
		}
	}
	
	inline static void assertNotEquals(void* obj1, void* obj2, const char* message = NULL)
	{
		if (obj1 == obj2)
		{
			throw AssertionException(message, NULL, ASSERT_NOT_EQUALS);
		}
	}
};

#endif
