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

void testToStringUnsigned()
{
	assertEquals(std::string(StringUtils::toString((unsigned)1)), std::string("1"), "string of unsigned is different");
}

void testToStringUnsignedLong()
{
	assertEquals(std::string(StringUtils::toString((unsigned long)1)), std::string("1"), "string of unsigned long is different");
}

void testToStringInt()
{
	assertEquals(std::string(StringUtils::toString(-1)), std::string("-1"), "string of integer is different");
}

void testToStringLong()
{
	assertEquals(std::string(StringUtils::toString(-1L)), std::string("-1"), "string of long is different");
}

void testToStringFloat()
{
	assertEquals(std::string(StringUtils::toString(-1.0F)), std::string("-1"), "string of float is different");
}

void testToStringDouble()
{
	assertEquals(std::string(StringUtils::toString(-1.0)), std::string("-1"), "string of double is different");
}

void testToStringChar()
{
	assertEquals(std::string(StringUtils::toString('a')), std::string("a"), "string of char is different");
}

int main(int argc, char** argv)
{
	testIsBlankBlankString();
	testIsBlankNotBlankString();
	
	testAreEqual();
	
	testTrimStart();
	testTrimEnd();
	testTrim();
	
	testToStringUnsigned();
	testToStringUnsignedLong();
	testToStringInt();
	testToStringLong();
	testToStringFloat();
	testToStringDouble();
	testToStringChar();
}
