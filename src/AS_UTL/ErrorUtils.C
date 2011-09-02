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

static const char* rcsid = "$Id: ErrorUtils.C,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include "ErrorUtils.h"

using namespace Utility;

void ErrorUtils::exceptionSignalHandler(int signalNum)
{
	throw SignalException(signalNum);
}

void ErrorUtils::printingSignalHandler(int signalNum)
{
	fprintf(stderr, "Signal ("F_STR"): "F_STR"\n", getSignalName(signalNum).c_str(), getSignalMessage(signalNum).c_str());
	
	ExceptionUtils::printStackTrace(stderr, DEFAULT_STACK_TRACE_LINE_DELIMITER, "\t", "ErrorUtils");
	
	signal(signalNum, SIG_DFL);
	raise(signalNum);
}

void ErrorUtils::handleErrorSignals(sighandler_t signalHandler)
{
	handleSignals(ERROR_SIGNALS, signalHandler);
}

void ErrorUtils::handleExitSignals(sighandler_t signalHandler)
{
	handleSignals(EXIT_SIGNALS, signalHandler);
}

void ErrorUtils::handleSignals(const int* signalNums, sighandler_t signalHandler)
{
	for (size_t a = 0; signalNums[a] != NO_SIGNAL; a++)
	{
		handleSignal(signalNums[a], signalHandler);
	}
}

sighandler_t ErrorUtils::handleSignal(int signalNum, sighandler_t signalHandler)
{
	return signal(signalNum, signalHandler);
}

string ErrorUtils::getSignalMessage(int signalNum, string message)
{
	return isSignal(signalNum) ? string(strsignal(signalNum)) + 
		(!message.empty() ? (StringUtils::startsWith(message, 1, ": ") ? ": " : "") + message : "") : string();
}

string ErrorUtils::getSignalName(int signalNum)
{
	return isSignal(signalNum) ? SIGNAL_NAMES[signalNum] : string();
}

bool ErrorUtils::isSignalIgnored(int signalNum)
{
	return getSignalHandler(signalNum) == SIG_IGN;
}

bool ErrorUtils::isSignalHandled(int signalNum)
{
	sighandler_t signalHandler = getSignalHandler(signalNum);
	
	return (signalHandler != SIG_ERR) && (signalHandler != SIG_IGN) && (signalHandler != SIG_DFL);
}

bool ErrorUtils::isSignal(int signalNum)
{
	return SIGNAL_NAMES.count(signalNum) != 0;
}

sighandler_t ErrorUtils::getSignalHandler(int signalNum)
{
	sighandler_t signalHandler = signal(signalNum, tempSignalHandler);
	signal(signalNum, signalHandler);
	
	return signalHandler;
}

string ErrorUtils::getError(string message)
{
	return getError(errno, message);
}

string ErrorUtils::getError(FILE* stream, string message)
{
	return getError(ferror(stream));
}

string ErrorUtils::getError(int errorNum, string message)
{
	return isError(errorNum) ? string(strerror(errorNum)) + 
		(!message.empty() ? (StringUtils::startsWith(message, 1, ": ") ? ": " : "") + message : "") : string();
}

bool ErrorUtils::hasError(FILE* stream)
{
	return isError(ferror(stream));
}

bool ErrorUtils::hasError()
{
	return isError(errno);
}

bool ErrorUtils::isError(int errorNum)
{
	return errorNum > NO_ERROR;
}
