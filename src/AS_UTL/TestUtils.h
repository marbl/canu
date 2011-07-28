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

#ifndef TEST_UTILS_H_
#define TEST_UTILS_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include "AssertionException.h"

static const char* RCSID_TEST_UTILS_H_ = "$Id: TestUtils.h,v 1.3 2011-07-28 11:31:00 mkotelbajcvi Exp $";

#define assertFalse(condition, message) \
	if (condition) \
	{ \
		throw AssertionException(message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertTrue(condition, message) \
	if (!condition) \
	{ \
		throw AssertionException(message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertEquals(obj1, obj2, message) \
	if (obj1 != obj2) \
	{ \
		throw AssertionException(message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#define assertNotEquals(obj1, obj2, message) \
	if (obj1 == obj2) \
	{ \
		throw AssertionException(message, __FILE__, __LINE__, __ASSERT_FUNCTION); \
	}

#endif
