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

static const char* rcsid_ASSERTS_H = "$Id: Asserts.h,v 1.4 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using namespace std;

#include "AS_global.h"
#include "AssertionException.h"
#include "AssertionType.h"

namespace Utility
{
	class Asserts
	{
	public:
		inline static void assertFalse(bool condition, string message = "")
		{
			if (condition)
			{
				throw AssertionException(ASSERT_FALSE, message);
			}
		}
		
		inline static void assertTrue(bool condition, string message = "")
		{
			if (!condition)
			{
				throw AssertionException(ASSERT_TRUE, message);
			}
		}
		
		inline static void assertNull(void* obj, string message = "")
		{
			if (obj != NULL)
			{
				throw AssertionException(ASSERT_NULL, message);
			}
		}
		
		inline static void assertNotNull(void* obj, const char* message = "")
		{
			if (obj == NULL)
			{
				throw AssertionException(ASSERT_NOT_NULL, message);
			}
		}
		
		inline static void assertEmpty(const char* str, string message = "")
		{
			if (!StringUtils::isEmpty(str))
			{
				throw AssertionException(ASSERT_EMPTY, message);
			}
		}
		
		inline static void assertEmpty(string str, string message = "")
		{
			if (!str.empty())
			{
				throw AssertionException(ASSERT_EMPTY, message);
			}
		}
		
		inline static void assertNotEmpty(const char* str, string message = "")
		{
			if (StringUtils::isEmpty(str))
			{
				throw AssertionException(ASSERT_NOT_EMPTY, message);
			}
		}
		
		inline static void assertNotEmpty(string str, string message = "")
		{
			if (str.empty())
			{
				throw AssertionException(ASSERT_NOT_EMPTY, message);
			}
		}
		
		template<class T>
		inline static void assertEquals(T obj1, T obj2, string message = "")
		{
			if (obj1 != obj2)
			{
				throw AssertionException(ASSERT_EQUALS, message);
			}
		}
		
		template<class T>
		inline static void assertNotEquals(T obj1, T obj2, string message = "")
		{
			if (obj1 == obj2)
			{
				throw AssertionException(ASSERT_NOT_EQUALS, message);
			}
		}
		
		inline static void assertSame(void* obj1, void* obj2, string message = "")
		{
			if (obj1 != obj2)
			{
				throw AssertionException(ASSERT_EQUALS, message);
			}
		}
		
		inline static void assertNotSame(void* obj1, void* obj2, string message = "")
		{
			if (obj1 == obj2)
			{
				throw AssertionException(ASSERT_NOT_EQUALS, message);
			}
		}
		
	private:
		Asserts()
		{
		}
	};
}

#endif
