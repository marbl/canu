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

static const char* rcsid_ERRORUTILS_H = "$Id: ErrorUtils.h,v 1.4 2011-09-06 17:14:15 mkotelbajcvi Exp $";

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

#ifdef __APPLE__
typedef void (*sighandler_t)(int)
#endif

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
		//signalNameMap[SIGSTKFLT] = string("SIGSTKFLT");
		//signalNameMap[SIGCLD] = string("SIGCLD");
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
		//signalNameMap[SIGPOLL] = string("SIGPOLL");
		signalNameMap[SIGIO] = string("SIGIO");
		//signalNameMap[SIGPWR] = string("SIGPWR");
		signalNameMap[SIGSYS] = string("SIGSYS");
		//signalNameMap[SIGUNUSED] = string("SIGUNUSED");
		
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
		inline static void handleErrorSignals(sighandler_t signalHandler)
		{
			handleSignals(ERROR_SIGNALS, signalHandler);
		}

		inline static void handleExitSignals(sighandler_t signalHandler)
		{
			handleSignals(EXIT_SIGNALS, signalHandler);
		}

		inline static void handleSignals(const int* signalNums, sighandler_t signalHandler)
		{
			for (size_t a = 0; signalNums[a] != NO_SIGNAL; a++)
			{
				handleSignal(signalNums[a], signalHandler);
			}
		}

		inline static sighandler_t handleSignal(int signalNum, sighandler_t signalHandler)
		{
			return signal(signalNum, signalHandler);
		}

		inline static string& getSignalMessage(int signalNum, string& buffer)
		{
			if (isSignal(signalNum))
			{
				buffer += strsignal(signalNum);
			}
			
			return buffer;
		}

		inline static string getSignalName(int signalNum)
		{
			return isSignal(signalNum) ? SIGNAL_NAMES[signalNum] : string();
		}

		inline static bool isSignalIgnored(int signalNum)
		{
			return getSignalHandler(signalNum) == SIG_IGN;
		}

		inline static bool isSignalHandled(int signalNum)
		{
			__sighandler_t signalHandler = getSignalHandler(signalNum);
			
			return (signalHandler != SIG_ERR) && (signalHandler != SIG_IGN) && (signalHandler != SIG_DFL);
		}

		inline static bool isSignal(int signalNum)
		{
			return SIGNAL_NAMES.count(signalNum) != 0;
		}

		inline static __sighandler_t getSignalHandler(int signalNum)
		{
			__sighandler_t signalHandler = signal(signalNum, tempSignalHandler);
			signal(signalNum, signalHandler);
			
			return signalHandler;
		}

		inline static string& getError(string& buffer)
		{
			return getError(errno, buffer);
		}

		inline static string& getError(FILE* stream, string& buffer)
		{
			return getError(ferror(stream), buffer);
		}

		inline static string& getError(int errorNum, string& buffer)
		{
			if (isError(errorNum))
			{
				buffer += strerror(errorNum);
			}
			
			return buffer;
		}

		inline static bool hasError(FILE* stream)
		{
			return isError(ferror(stream));
		}

		inline static bool hasError()
		{
			return isError(errno);
		}

		inline static bool isError(int errorNum)
		{
			return errorNum > NO_ERROR;
		}
		
		static void exceptionSignalHandler(int signalNum);
		static void printingSignalHandler(int signalNum);
		
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
