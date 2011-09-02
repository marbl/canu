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

#ifndef ERRORUTILS_H
#define ERRORUTILS_H

static const char* rcsid_ERRORUTILS_H = "$Id: ErrorUtils.h,v 1.1 2011-09-02 14:59:27 mkotelbajcvi Exp $";

#include <unistd.h>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

using namespace std;

#include "AS_global.h"
#include "ExceptionUtils.h"
#include "RuntimeException.h"
#include "SignalException.h"
#include "StreamUtils.h"
#include "StringUtils.h"

namespace Utility
{
	static map<int, string> getSignalNameMap()
	{
		map<int, string> signalNameMap;
		signalNameMap[SIGHUP] = string("SIGHUP");
		signalNameMap[SIGINT] = string("SIGINT");
		signalNameMap[SIGQUIT] = string("SIGQUIT");
		signalNameMap[SIGILL] = string("SIGILL");
		signalNameMap[SIGTRAP] = string("SIGTRAP");
		signalNameMap[SIGABRT] = string("SIGABRT");
		signalNameMap[SIGIOT] = string("SIGIOT");
		signalNameMap[SIGBUS] = string("SIGBUS");
		signalNameMap[SIGFPE] = string("SIGFPE");
		signalNameMap[SIGKILL] = string("SIGKILL");
		signalNameMap[SIGUSR1] = string("SIGUSR1");
		signalNameMap[SIGSEGV] = string("SIGSEGV");
		signalNameMap[SIGUSR2] = string("SIGUSR2");
		signalNameMap[SIGPIPE] = string("SIGPIPE");
		signalNameMap[SIGALRM] = string("SIGALRM");
		signalNameMap[SIGTERM] = string("SIGTERM");
		signalNameMap[SIGSTKFLT] = string("SIGSTKFLT");
		signalNameMap[SIGCLD] = string("SIGCLD");
		signalNameMap[SIGCHLD] = string("SIGCHLD");
		signalNameMap[SIGCONT] = string("SIGCONT");
		signalNameMap[SIGSTOP] = string("SIGSTOP");
		signalNameMap[SIGTSTP] = string("SIGTSTP");
		signalNameMap[SIGTTIN] = string("SIGTTIN");
		signalNameMap[SIGTTOU] = string("SIGTTOU");
		signalNameMap[SIGURG] = string("SIGURG");
		signalNameMap[SIGXCPU] = string("SIGXCPU");
		signalNameMap[SIGXFSZ] = string("SIGXFSZ");
		signalNameMap[SIGVTALRM] = string("SIGVTALRM");
		signalNameMap[SIGPROF] = string("SIGPROF");
		signalNameMap[SIGWINCH] = string("SIGWINCH");
		signalNameMap[SIGPOLL] = string("SIGPOLL");
		signalNameMap[SIGIO] = string("SIGIO");
		signalNameMap[SIGPWR] = string("SIGPWR");
		signalNameMap[SIGSYS] = string("SIGSYS");
		signalNameMap[SIGUNUSED] = string("SIGUNUSED");
		
		return signalNameMap;
	}

	static map<int, string> SIGNAL_NAMES = getSignalNameMap();

	static const int NO_SIGNAL = 0;
	static const int ERROR_SIGNALS[] = { SIGFPE, SIGILL, SIGINT, SIGSEGV, NO_SIGNAL };
	static const int EXIT_SIGNALS[] = { SIGABRT, SIGQUIT, SIGTERM, NO_SIGNAL };
	
	static const int NO_ERROR = 0;
	
	class ErrorUtils
	{
	public:
		template<class E>
		__inline__ static void throwIfError(string message = string())
		{
			throwIfError<E>(errno, message);
		}
		
		template<class E>
		__inline__ static void throwIfError(FILE* stream, string message = string())
		{
			throwIfError<E>(ferror(stream), message);
		}
		
		template<class E>
		__inline__ static void throwIfError(int errorNum, string message = string())
		{
			if (ErrorUtils::isError(errorNum))
			{
				throw E(getError(errorNum, message));
			}
		}
		
		static void exceptionSignalHandler(int signalNum);
		static void printingSignalHandler(int signalNum);
		
		static void handleErrorSignals(sighandler_t signalHandler);
		static void handleExitSignals(sighandler_t signalHandler);
		static void handleSignals(const int* signalNums, sighandler_t signalHandler);
		static sighandler_t handleSignal(int signalNum, sighandler_t signalHandler);
		
		static string getSignalMessage(int signalNum, string message = string());
		static string getSignalName(int signalNum);
		
		static bool isSignalIgnored(int signalNum);
		static bool isSignalHandled(int signalNum);
		static bool isSignal(int signalNum);
		
		static sighandler_t getSignalHandler(int signalNum);
		
		static string getError(string message = string());
		static string getError(FILE* stream, string message = string());
		static string getError(int errorNum, string message = string());
		
		static bool hasError();
		static bool hasError(FILE* stream);
		
		static bool isError(int errorNum);
		
	protected:
		static void tempSignalHandler(int signalNum)
		{
		}
		
	private:
		ErrorUtils()
		{
		}
	};
}

#endif
