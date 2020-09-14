
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_Logging.H"

#include <stdarg.h>


class logFileInstance {
public:
  logFileInstance()   { clear(); };
  ~logFileInstance()  { close(); };

  void  clear(void);

  void  set(char const *prefix_, int32 order_, char const *label_, int32 thread_, bool stderr_);

  void  write(char const *fmt, va_list ap);

  void  rotate(uint64 limit);
  void  open(void);
  void  close(void);

  void  flush(void);

  const
  char *prefix(void)  { return(_prefix); };

private:
  FILE   *_file;
  char    _prefix[FILENAME_MAX+1];   //  "%s.%03u.%s"
  char    _name  [FILENAME_MAX+1];   //  "%s.%03u.%s.thr%03d"
  char    _path  [FILENAME_MAX+1];   //  "%s.%03u.%s.thr%03d.num%03d.log"
  uint32  _part;
  uint64  _length;
};



//  Set a new "label" for the log files.
//
//  If the label is nullptr, just close the files and make it look like we
//  were just constructed.
//
//  If stderr_ is false, _file will remain unset, and we'll open a file for
//  it in open().
//
void
logFileInstance::set(char const *prefix_, int32 order_, char const *label_, int32 thread_, bool stderr_) {

  assert(prefix_ != nullptr);

  if (label_ == nullptr) {
    close();
    return;
  }

  if (stderr_ == true)
    _file = stderr;

  snprintf(_prefix, FILENAME_MAX, "%s.%03u.%s",         prefix_, order_, label_);
  snprintf(_name,   FILENAME_MAX, "%s.%03u.%s.thr%03d", prefix_, order_, label_, thread_);
}



//  If the current file size is more than the limit, close the current file
//  and increment the part number.  The next file will be opened on the next
//  write.
void
logFileInstance::rotate(uint64 limit) {

  if ((_file == nullptr) ||   //  Is nullptr on the first call to rotate(), before open() is called.
      (_file == stderr) ||
      (_length < limit))
    return;

  fprintf(_file, "logFile()--  size " F_U64 " exceeds limit of " F_U64 "; rotate to new file.\n",
          _length, limit);

  AS_UTL_closeFile(_file, _name);

  _file    = nullptr;
  _part   += 1;
  _length  = 0;
}



//  Open a log file.  If there is already a file open, do nothing - this is
//  actually the usual case.
//  
void
logFileInstance::open(void) {

  if (_file != nullptr)
    return;

  snprintf(_path, FILENAME_MAX, "%s.num%03d.log", _name, _part);

  _file = fopen(_path, "w");

  if (_file == nullptr) {
    writeStatus("setLogFile()-- Failed to open logFile '%s': %s.\n", _path, strerror(errno));
    writeStatus("setLogFile()-- Will now log to stderr instead.\n");
    _file = stderr;
  }
}



void
logFileInstance::clear(void) {
  _file      = nullptr;
  _prefix[0] = 0;
  _name[0]   = 0;
  _path[0]   = 0;
  _part      = 0;
  _length    = 0;
}



void
logFileInstance::close(void) {
  AS_UTL_closeFile(_file, _path);
  clear();
}



void
logFileInstance::flush(void) {
  if (_file != nullptr)
    fflush(_file);
}



void
logFileInstance::write(const char *fmt, va_list ap) {
  assert(_file != nullptr);
  _length += vfprintf(_file, fmt, ap);
}





//  NONE of the logFileMain/logFileThread is implemented


logFileInstance    logFileMain;           //  For writes during non-threaded portions
logFileInstance   *logFileThread = nullptr;  //  For writes during threaded portions.
uint32             logFileOrder  = 0;
uint64             logFileFlags  = 0;

uint64 LOG_OVERLAP_SCORING             = 0x0000000000000001;  //  Debug, scoring of overlaps
uint64 LOG_BEST_EDGES                  = 0x0000000000000002;
uint64 LOG_BEST_OVERLAPS               = 0x0000000000000004;
uint64 LOG_ERROR_PROFILES              = 0x0000000000000008;
uint64 LOG_OPTIMIZE_POSITIONS          = 0x0000000000000010;
uint64 LOG_CHUNK_GRAPH                 = 0x0000000000000020;  //  Report the chunk graph as we build it
uint64 LOG_BUILD_UNITIG                = 0x0000000000000040;  //  Report building of initial tigs (both unitig creation and read placement)
uint64 LOG_PLACE_UNPLACED              = 0x0000000000000080;  //  Report placing of unplaced reads
uint64 LOG_ORPHAN_DETAIL               = 0x0000000000000100;
uint64 LOG_SPLIT_DISCONTINUOUS         = 0x0000000000000200;  //
uint64 LOG_INTERMEDIATE_TIGS           = 0x0000000000000400;  //  At various spots, dump the current tigs
uint64 LOG_SET_PARENT_AND_HANG         = 0x0000000000000800;  //
uint64 LOG_STDERR                      = 0x0000000000001000;  //  Write ALL logging to stderr, not the files.

uint64 LOG_PLACE_READ                  = 0x8000000000000000;  //  Internal use only.

char const *logFileFlagNames[64] = { "overlapScoring",
                                     "bestEdges",
                                     "bestOverlaps",
                                     "errorProfiles",
                                     "optimizePositions",
                                     "chunkGraph",
                                     "buildUnitig",
                                     "placeUnplaced",
                                     "orphans",
                                     "splitDiscontinuous",   //  Update made it to here, need repeats
                                     "intermediateTigs",
                                     "setParentAndHang",
                                     "stderr",
                                     nullptr
};



//  Closes the current logFile, opens a new one called
//  'prefix.logFileOrder.label'.
//
//  If 'label' is nullptr, the logFile is reset to stderr.
//
void
setLogFile(char const *prefix, char const *label) {

  //  Allocate space.  Unfortunately, this leaks.

  if (logFileThread == nullptr)
    logFileThread = new logFileInstance [omp_get_max_threads()];

  //  Close out the old.

  logFileMain.close();

  for (int32 tn=0; tn<omp_get_max_threads(); tn++)
    logFileThread[tn].close();

  //  Move to the next iteration.

  logFileOrder++;

  //  Set up for that iteration.

  logFileMain.set(prefix, logFileOrder, label, 0, logFileFlagSet(LOG_STDERR));

  for (int32 tn=0; tn<omp_get_max_threads(); tn++)
    logFileThread[tn].set(prefix, logFileOrder, label, tn+1, logFileFlagSet(LOG_STDERR));

  //  File open is delayed until it is used.

}



char const *
getLogFilePrefix(void) {
  return(logFileMain.prefix());
}



void
writeStatus(char const *fmt, ...) {
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}



void
writeLog(char const *fmt, ...) {
  va_list ap;
  int32   nt = omp_get_num_threads();
  int32   tn = omp_get_thread_num();

  logFileInstance  *lf = (nt == 1) ? (&logFileMain) : (&logFileThread[tn]);

  //  Close big files and make sure we have a file opened.

  lf->rotate(512 * 1024 * 1024);
  lf->open();

  //  Write the log.

  va_start(ap, fmt);
  lf->write(fmt, ap);
  va_end(ap);
}



void
flushLog(void) {
  int32   nt = omp_get_num_threads();
  int32   tn = omp_get_thread_num();

  if (nt == 1)
    logFileMain.flush();
  else
    logFileThread[tn].flush();
}
