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

#ifndef AS_UTL_TEST_H
#define AS_UTL_TEST_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <iostream>

static const char* rcsid = "$Id: TestUtils.h,v 1.1 2011-07-28 04:23:50 mkotelbajcvi Exp $";

#define assertFalse(condition, message) \
	if (condition) \
	{ \
		if (message != NULL) \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert false failed: " << message; \
		} \
		else \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert false failed."; \
		} \
		\
		exit(EXIT_FAILURE); \
	}

#define assertTrue(condition, message) \
	if (!condition) \
	{ \
		if (message != NULL) \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert true failed: " << message; \
		} \
		else \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert true failed."; \
		} \
		\
		exit(EXIT_FAILURE); \
	}

#define assertEquals(obj1, obj2, message) \
	if (obj1 != obj2) \
	{ \
		if (message != NULL) \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert equals failed: " << message << "\nobj1=" \
				<< obj1 << "\nobj2=" << obj2; \
		} \
		else \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert equals failed." << "\nobj1=" \
				<< obj1 << "\nobj2=" << obj2; \
		} \
		\
		exit(EXIT_FAILURE); \
	}

#define assertNotEquals(obj1, obj2, message) \
	if (obj1 != obj2) \
	{ \
		if (message != NULL) \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert not equals failed: " << message << "\nobj1=" \
				<< obj1 << "\nobj2=" << obj2; \
		} \
		else \
		{ \
			std::cerr << __FILE__ << ":" << __LINE__ << " " << __ASSERT_FUNCTION << " Assert not equals failed." << "\nobj1=" \
				<< obj1 << "\nobj2=" << obj2; \
		} \
		\
		exit(EXIT_FAILURE); \
	}

#endif
