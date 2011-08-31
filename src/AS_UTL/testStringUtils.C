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

static const char* rcsid = "$Id: testStringUtils.C,v 1.11 2011-08-31 06:49:27 mkotelbajcvi Exp $";

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

using namespace std;

#include "Asserts.h"
#include "StringUtils.h"
#include "TestUtils.h"

void testConcat()
{
	string result;
	result.resize(256);
	
	Asserts::assertTrue(StringUtils::concat(result, F_U32""F_STR, 0, "|test1|") == "0|test1|", "concat failed");
}

void testFindAll()
{
	string str("|test1|test2|test3|");
	vector<size_t> indexes;
	size_t expectedIndexes[3] = { 1, 7, 13 };
	
	Asserts::assertTrue(StringUtils::findAll(str, indexes, 1, "test") == vector<size_t>(expectedIndexes, expectedIndexes + 3), "find all failed");
}

void testSplit()
{
	string str("|test1%test2!test3$");
	vector<string> tokens;
	string expectedTokens[5] = { string(), string("test1"), string("test2"), string("test3"), string() };
	
	Asserts::assertTrue(StringUtils::split(str, tokens, 4, "|", "%", "!", "$") == vector<string>(expectedTokens, expectedTokens + 5), "split failed");
}

void testIsBlankBlankString()
{
	Asserts::assertTrue(StringUtils::isBlank(" \t\n\r"), "string blank");
}

void testIsBlankNotBlankString()
{
	Asserts::assertFalse(StringUtils::isBlank(" test\t\n\r"), "string not blank");
}

void testAreEqual()
{
	Asserts::assertTrue(StringUtils::areEqual("test", "test"), "strings are equal");
}

void testStartsWith()
{
	Asserts::assertTrue(StringUtils::startsWith("test1test2test3", 2, "test1", "test2"), "starts with failed");
}

void testEndsWith()
{
	Asserts::assertTrue(StringUtils::endsWith("test1test2test3", 2, "test3", "test2"), "ends with failed");
}

void testTrimStart()
{
	string str("trim1trim2|test");
	
	Asserts::assertTrue(StringUtils::trimStart(str, 2, "trim1", "trim2") == "|test", "start trimmed string is different");
}

void testTrimEnd()
{
	string str("test|trim2trim1");
	
	Asserts::assertTrue(StringUtils::trimEnd(str, 2, "trim1", "trim2") == "test|", "end trimmed string is different");
}

void testTrim()
{
	string str("trim1trim2|test|trim2trim1");
	
	Asserts::assertTrue(StringUtils::trim(str, 2, "trim1", "trim2") == "|test|", "trimmed string is different");
}

void testToStringUnsigned()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString((unsigned)1, str) == "1", "string of unsigned is different");
}

void testToStringUnsignedLong()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString((unsigned long)1, str) == "1", "string of unsigned long is different");
}

void testToStringInt()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString(-1, str) == "-1", "string of integer is different");
}

void testToStringLong()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString(-1L, str) == "-1", "string of long is different");
}

void testToStringFloat()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString(-1.0F, str) == "-1", "string of float is different");
}

void testToStringDouble()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString(-1.0, str) == "-1", "string of double is different");
}

void testToStringChar()
{
	string str;
	
	Asserts::assertTrue(StringUtils::toString('a', str) == "a", "string of char is different");
}

int main(int argc, char** argv)
{
	vector<TestFunction> tests;
	tests.push_back(&testConcat);
	tests.push_back(&testFindAll);
	tests.push_back(&testSplit);
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
