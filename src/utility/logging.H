
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2019-OCT-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef LOGGING_H
#define LOGGING_H

#include "AS_global.H"
#include "files.H"

const uint32 logFileLevelNameLen = 64;

class logFileInstance;
class logFileLevel;




class logFileHandle {
public:
  logFileHandle()              { _index = UINT32_MAX; };
private:
  logFileHandle(uint32 index)  { _index = index;      };

private:
  uint32  _index;

  friend class logFile;
};




class logFile {
public:
  logFile(char const *prefix=NULL, uint64 maxSize=512 * 1024 * 1024);
  ~logFile();

  //  All logging is written to files with the same 'prefix', set at
  //  logFile() construction time (or via setPrefix().  getPrefix()
  //  will return that name.
  //
  //  getLogName() will return the current log file name being written
  //  to, something like 'prefix.###.name'.  If called in threaded
  //  code, it does NOT append the thread number; it returns the same
  //  string regardless of thread status.
  //
  void        setPrefix(char const *prefix);
  char const *getPrefix(void);
  char const *getLogName(void);

  //  Each log file is restricted to some maximum size in bytes.  Files will
  //  roll over to the next numbered file.  This is to keep files to a
  //  reasonable size so they can be loaded into an editor.
  //
  void        setMaxSize(uint64 size);

  //  Each call to setName() will:
  //    close any existing log files
  //    increment the index number of the logging
  //    open new log files 'prefix.###.name'.
  //
  void        setName(char const *name);


  //  CONFIGURE LOGGING LEVELS
  //
  //  Add a logging level with name 'levelName' at verbosity 0 (highest).
  //
  //  Add a logging level with name 'levelName' at the specified verbosity.
  //  Vebosity ranges from 0 (highest, always output) to MAX_INT.
  //
  logFileHandle  addLevel(char const *levelName,                   bool enabled=false);
  logFileHandle  addLevel(char const *levelName, uint32 verbosity, bool enabled=false);


  //  ENABLE/DISABLE LOGGING LEVELS
  //
  //  Enable (disable) logging to the given levelName or verbosity.
  //
  //  The first enable() is a convenience function for command line
  //  processing.  It is expecting the option string as the first argument,
  //  and a level name as the second argument.  If a NULL level name is
  //  supplied, the global verbosity level is increased.  In both cases, the
  //  length of the option string (excluding any leading dashes) is used as
  //  the verbosity increment.  It returns a suitable increment to arg.
  //
  //    arg += enable("-vvvv", NULL);        //  Adds zero to arg.
  //    arg += enable("-DDDD", "levelName);  //  Adds one to arg.
  //
  int32       enable   (char const *optionString, char const *levelName);

  void        enable   (char const *levelName, uint32 verbosity=0);
  void        disable  (char const *levelName);
  void        increment(char const *levelName);

  void        enable   (logFileHandle levelName, uint32 verbosity=0);
  void        disable  (logFileHandle levelName);
  void        increment(logFileHandle levelName);

  void        enable   (uint32 verbosity);
  void        increment(void);


  //  OUTPUT LOGGING AND STATUS MESSAGES
  //
  //  Status messages go to stderr, log messages go to the currently active
  //  log file.
  //
  //  A message is written if:
  //    the specified levelName is enabled.            (e.g., -D debug)
  //    the specified verbosity is no more than set.   (e.g., -V or -VV)
  //    both conditions are met                        (e.g., -D debug -VV)
  //
  //  Log messages are buffered.  flush() will write the buffer to disk.
  //
private:
  void        writeStatus(char const *fmt, va_list ap);
  void        writeLog   (char const *fmt, va_list ap);

  bool        verbosityEnabled(uint32 verbosity);
  bool        levelEnabled(logFileHandle level, uint32 verbosity=0);

public:
  void        writeStatus(char const *fmt, ...);
  void        writeLog   (char const *fmt, ...);

  void        writeStatus(logFileHandle levelName, char const *fmt, ...);
  void        writeLog   (logFileHandle levelName, char const *fmt, ...);

  void        writeStatus(uint32 verbosity, char const *fmt, ...);
  void        writeLog   (uint32 verbosity, char const *fmt, ...);

  void        writeStatus(logFileHandle levelName, uint32 verbosity, char const *fmt, ...);
  void        writeLog   (logFileHandle levelName, uint32 verbosity, char const *fmt, ...);

  void        flush(void);

private:
  logFileInstance             *_mainI;
  logFileInstance            **_threadI;

  uint32                        findLevelIndex(char const *levelName);

  uint32                       _levelsLen;
  uint32                       _levelsMax;
  logFileLevel               **_levels;

  uint32                       _verbosity;
};


#endif  //  LOGGING_H

