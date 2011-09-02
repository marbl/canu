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

static const char* rcsid = "$Id: TestUtils.C,v 1.4 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "TestUtils.h"

using namespace Utility;

void TestUtils::runTests(vector<TestFunction>& tests)
{
	ErrorUtils::handleErrorSignals(ErrorUtils::printingSignalHandler);
	atexit(printResults);
	
	runningTests = tests;
	successful = 0, errors = 0;
	
	for (size_t a = 0; a < runningTests.size(); a++)
	{
		try
		{
			runningTests[a]();
			
			successful++;
		}
		catch (exception& e)
		{
			fprintf(stderr, "TEST ERROR => %s\n", e.what());
			
			errors++;
		}
	}
}

void TestUtils::printResults()
{
	fprintf(stderr, "TESTS FINISHED => "F_U64" total, "F_U64" successful, "F_U64" errors\n", runningTests.size(), successful, errors);
}
