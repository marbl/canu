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
#include <exception>
#include <string>
#include <vector>

using namespace std;

#include "StringUtils.h"
#include "TestUtils.h"

static const char* RCSID = "$Id: testStringUtils.C,v 1.4 2011-08-01 20:33:36 mkotelbajcvi Exp $";

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

void testStartsWith()
{
	assertTrue(StringUtils::startsWith("test1test2test3", 2, "test1", "test2"), "starts with failed");
}

void testEndsWith()
{
	assertTrue(StringUtils::endsWith("test1test2test3", 2, "test3", "test2"), "ends with failed");
}

void testTrimStart()
{
	assertEquals(string(StringUtils::trimStart("trim1trim2|test", 2, "trim1", "trim2")), string("|test"), "start trimmed string is different");
}

void testTrimEnd()
{
	assertEquals(string(StringUtils::trimEnd("test|trim2trim1", 2, "trim1", "trim2")), string("test|"), "end trimmed string is different");
}

void testTrim()
{
	assertEquals(string(StringUtils::trim("trim1trim2|test|trim2trim1", 2, "trim1", "trim2")), string("|test|"), "trimmed string is different");
}

void testToStringUnsigned()
{
	assertEquals(string(StringUtils::toString((unsigned)1)), string("1"), "string of unsigned is different");
}

void testToStringUnsignedLong()
{
	assertEquals(string(StringUtils::toString((unsigned long)1)), string("1"), "string of unsigned long is different");
}

void testToStringInt()
{
	assertEquals(string(StringUtils::toString(-1)), string("-1"), "string of integer is different");
}

void testToStringLong()
{
	assertEquals(string(StringUtils::toString(-1L)), string("-1"), "string of long is different");
}

void testToStringFloat()
{
	assertEquals(string(StringUtils::toString(-1.0F)), string("-1"), "string of float is different");
}

void testToStringDouble()
{
	assertEquals(string(StringUtils::toString(-1.0)), string("-1"), "string of double is different");
}

void testToStringChar()
{
	assertEquals(string(StringUtils::toString('a')), string("a"), "string of char is different");
}

int main(int argc, char** argv)
{
	vector<TestFunction> tests;
	tests.push_back(&testIsBlankBlankString);
	tests.push_back(&testIsBlankNotBlankString);
	tests.push_back(&testAreEqual);
	tests.push_back(&testStartsWith);
	tests.push_back(&testEndsWith);
	tests.push_back(&testTrimStart);
	tests.push_back(&testTrimEnd);
	tests.push_back(&testTrim);
	tests.push_back(&testToStringUnsigned);
	tests.push_back(&testToStringUnsignedLong);
	tests.push_back(&testToStringInt);
	tests.push_back(&testToStringLong);
	tests.push_back(&testToStringFloat);
	tests.push_back(&testToStringDouble);
	tests.push_back(&testToStringChar);
	
	TestUtils::runTests(tests);
}
