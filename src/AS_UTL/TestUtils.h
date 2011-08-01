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

#ifndef TESTUTILS_H
#define TESTUTILS_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vector>

using namespace std;

#include "AS_global.h"
#include "AssertionException.h"
#include "AssertionType.h"
#include "StringUtils.h"

static const char* RCSID_TESTUTILS_H = "$Id: TestUtils.h,v 1.5 2011-08-01 20:33:35 mkotelbajcvi Exp $";

typedef void (*TestFunction)();

#define assertFalse(condition, message) \
	if (condition) \
	{ \
		throw AssertionException(ASSERT_FALSE, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertTrue(condition, message) \
	if (!condition) \
	{ \
		throw AssertionException(ASSERT_TRUE, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertEquals(obj1, obj2, message) \
	if (obj1 != obj2) \
	{ \
		throw AssertionException(ASSERT_EQUALS, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertNotEquals(obj1, obj2, message) \
	if (obj1 == obj2) \
	{ \
		throw AssertionException(ASSERT_NOT_EQUALS, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertNull(obj, message) \
	if (obj != NULL) \
	{ \
		throw AssertionException(ASSERT_NULL, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertNotNull(obj, message) \
	if (obj == NULL) \
	{ \
		throw AssertionException(ASSERT_NOT_NULL, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertEmpty(str, message) \
	if (strlen(str) != 0) \
	{ \
		throw AssertionException(ASSERT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertNotEmpty(str, message) \
	if (strlen(str) == 0) \
	{ \
		throw AssertionException(ASSERT_NOT_EMPTY, message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

class TestUtils
{
public:
	static void runTests(vector<TestFunction>& tests);
};

#endif
