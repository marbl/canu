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

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string>

#include "StringUtils.h"
#include "TestUtils.h"

void testIsBlankBlankString()
{
	assertTrue(StringUtils::isBlank(" \t\n\r"), "string blank");
}

void testIsBlankNotBlankString()
{
	assertFalse(StringUtils::isBlank(" test\t\n\r"), "string not blank");
}

void testAreEqual()
{
	assertTrue(StringUtils::areEqual("test", "test"), "strings are equal");
}

void testTrimStart()
{
	assertEquals(std::string(StringUtils::trimStart("trim1trim2|test", 2, "trim1", "trim2")), std::string("|test"), "start trimmed string is different");
}

void testTrimEnd()
{
	assertEquals(std::string(StringUtils::trimEnd("test|trim2trim1", 2, "trim1", "trim2")), std::string("test|"), "end trimmed string is different");
}

void testTrim()
{
	assertEquals(std::string(StringUtils::trim("trim1trim2|test|trim2trim1", 2, "trim1", "trim2")), std::string("|test|"), "trimmed string is different");
}

int main(int argc, char** argv)
{
	testIsBlankBlankString();
	testIsBlankNotBlankString();
	
	testAreEqual();
	
	testTrimStart();
	testTrimEnd();
	testTrim();
}
