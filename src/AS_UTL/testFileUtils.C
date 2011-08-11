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

static const char* rcsid = "$Id: testFileUtils.C,v 1.4 2011-08-11 17:34:34 mkotelbajcvi Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>
#include <vector>

using namespace std;

#include "Asserts.h"
#include "FileUtils.h"
#include "TestUtils.h"

void testGetPath()
{
	Asserts::assertTrue(string(FileUtils::getPath(3, "/test1", "//test2", "test3/")) == "/test1/test2/test3/", "get path was different");
}

int main(int argc, char** argv)
{
	vector<TestFunction> tests;
	tests.push_back(&testGetPath);
	
	TestUtils::runTests(tests);
}
