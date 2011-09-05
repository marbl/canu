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

static const char* rcsid = "$Id: testReflectionUtils.C,v 1.1 2011-09-05 16:49:45 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <typeinfo>
#include <vector>

using namespace std;

#include "Asserts.h"
#include "ReflectionUtils.h"
#include "TestUtils.h"

using namespace Utility;

void testDemangleTypeName()
{
	string str = typeid(ReflectionUtils).name();
	
	ReflectionUtils::demangle(str);
	
	Asserts::assertEquals(str, string("Utility::ReflectionUtils"), "demangling type name failed");
}

int main(int argc, char** argv)
{
	vector<TestFunction> tests;
	tests.push_back(&testDemangleTypeName);
	
	TestUtils::runTests(tests);
}
