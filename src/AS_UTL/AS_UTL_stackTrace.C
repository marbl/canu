
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2013, J. Craig Venter Institute
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

#ifndef AS_UTL_STACKTRACE_H
#define AS_UTL_STACKTRACE_H

static const char *rcsid_AS_UTL_STACKTRACE_H = "$Id: AS_UTL_stackTrace.C,v 1.1 2013-03-14 06:24:13 brianwalenz Exp $";

#include "AS_global.h"

//  Derived from
//    http://sourceware.org/git/?p=glibc.git;a=blob;f=debug/segfault.c
//
//  Linux    g++ -rdynamic -g3 -o st1 st1.C
//  FreeBSD  CC -o unwind unwind.C -I/usr/local/include -L/usr/local/lib -lunwind

#include <alloca.h>
#include <execinfo.h>  //  backtrace
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cxxabi.h>

#define WRITE_STRING(S) write(2, S, strlen(S))


#ifdef LIBUNWIND

#include <libunwind.h>

void
AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx) {

  WRITE_STRING("\nFailed with '");

#if 0
  //  GNUisms below _NSIG And _sys_siglist
  if ((sig_num < 0) || (sig_num >= _NSIG) || (_sys_siglist[sig_num] == NULL))
    WRITE_STRING("Unknown sig_num");
  else
    WRITE_STRING (_sys_siglist[sig_num]);
#endif

  WRITE_STRING("(unknown)'\n\n");

  unw_cursor_t   cursor;
  unw_context_t  uc;
  unw_word_t     ip, sp;
  int            depth = 0;

  unw_getcontext(&uc);           //  Get state
  unw_init_local(&cursor, &uc);  //  Initialize state cursor

  while (unw_step(&cursor) > 0) {
    unw_get_reg(&cursor, UNW_REG_IP, &ip);
    unw_get_reg(&cursor, UNW_REG_SP, &sp);

    unw_word_t off;
    char       name[256];

    if (unw_get_proc_name (&cursor, name, sizeof (name), &off) == 0) {
      if (off)
        fprintf(stderr, "%02d <%s + 0x%lx>  ip=%lx sp=%lx\n", depth, name, (long)off, ip, sp);
      else
        fprintf(stderr, "%02d <%s>  ip=%lx sp=%lx\n", depth, name, ip, sp);
    } else {
      fprintf(stderr, "%02d <???>  ip=%lx sp=%lx\n", depth, name, ip, sp);
    }

    depth++;
  }
}

#elif defined(_GNU_SOURCE)

void
AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx) {
  void  *arr[256];
  int32  cnt = backtrace(arr, 256);

  //  Report the signal we failed on, be careful to not allocate memory.

  WRITE_STRING("\nFailed with '");

  if ((sig_num < 0) || (sig_num >= _NSIG) || (_sys_siglist[sig_num] == NULL))
    WRITE_STRING("Unknown sig_num");
  else
    WRITE_STRING(_sys_siglist[sig_num]);

  WRITE_STRING("'\n");

  //  Dump a full backtrace, even including the signal handler frames, instead of not dumping anything.

  //  The _fd version doesn't allocate memory, the standard backtrace_symbols does.  Demangling
  //  names also allocates memory.  We'll do the safe one first, just to get something, then try the
  //  others.

  WRITE_STRING("\nBacktrace (mangled):\n\n");

  backtrace_symbols_fd(arr, cnt, 2);

  //  Now that we have something, try to generate a pretty backtrace.

  WRITE_STRING("\nBacktrace (demangled):\n\n");

  char **messages = backtrace_symbols(arr, cnt);    

  for (int32 i=0; (i < cnt) && (messages != NULL); ++i) {
    char *mang = NULL;
    char *obgn = NULL;
    char *oend = NULL;

    //  find parantheses and +address offset surrounding mangled name

    for (char *p=messages[i]; (*p) && (oend == NULL); ++p) {
      if      (*p == '(') 
        mang = p; 

      else if (*p == '+') 
        obgn = p;

      else if (*p == ')')
        oend = p;
    }

    //  If all three are valid, attempt to demangle the symbol

    if ((mang != NULL) &&
        (obgn != NULL)  &&
        (oend != NULL)  &&
        (mang < obgn)) {

      *mang++ = 0;
      *obgn++ = 0;
      *oend++ = 0;

      int32 status = 0;
      char *name   = abi::__cxa_demangle(mang, 0, 0, &status);

      if (status > 0)
        //  Failed.
        name = mang;

      fprintf(stderr, "[%d] %s::%s + %s %s\n", i, messages[i], name, obgn, oend);

      //free(name);
    } else {
      fprintf(stderr, "[%d] %s\n", i, messages[i]);
    }
  }

  WRITE_STRING("\n");

  //  Pass the signal through, only so a core file can get generated.

  struct sigaction sa;

  sa.sa_handler = SIG_DFL;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;

  sigaction(sig_num, &sa, NULL);

  raise(sig_num);
}

#else

void
AS_UTL_catchCrash(int sig_num, siginfo_t *info, void *ctx) {
}

#endif


void
AS_UTL_installCrashCatcher(void) {
  struct sigaction sigact;

  sigact.sa_sigaction = AS_UTL_catchCrash;
  sigact.sa_flags     = SA_RESTART | SA_SIGINFO;

  //  Don't especially care if these fail or not.

  //sigaction(SIGINT,  &sigact, NULL);  //  02 - Interrupt
  sigaction(SIGILL,  &sigact, NULL);  //  04 - Illegal instruction
  sigaction(SIGFPE,  &sigact, NULL);  //  06 - Floating point exception
  sigaction(SIGBUS,  &sigact, NULL);  //  10 - Bus error
  sigaction(SIGSEGV, &sigact, NULL);  //  11 - Segmentation fault
}


#endif  //  AS_UTL_STACKTRACE_H

