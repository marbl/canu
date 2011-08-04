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

static const char* rcsid = "$Id: testExceptions.C,v 1.1 2011-08-04 19:10:59 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>
#include <vector>

using namespace std;

#include "ArgumentException.h"
#include "Assert.h"
#include "AssertionException.h"
#include "IllegalStateException.h"
#include "StringUtils.h"
#include "TestUtils.h"

void testCauseDepth()
{
	Assert::assertTrue(StringUtils::findAll(ArgumentException("exception1", 
		new ArgumentException("exception2", new ArgumentException("exception3", NULL, "arg3"), "arg2"), "arg1"), 
		1, "Caused by: ").size() == (MAX_CAUSE_DEPTH - 1), "cause depth failed");
}

int main(int argc, char** argv)
{
	vector<TestFunction> tests;
	tests.push_back(&testCauseDepth);
	
	TestUtils::runTests(tests);
}
