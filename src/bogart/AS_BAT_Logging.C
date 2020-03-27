
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
  logFileInstance() {
    file      = stderr;
    prefix[0] = 0;
    name[0]   = 0;
    part      = 0;
    length    = 0;
  };
  ~logFileInstance() {
    if ((name[0] != 0) && (file)) {
      fprintf(stderr, "WARNING: open file '%s'\n", name);
      AS_UTL_closeFile(file, name);
    }
  };

  void  set(char const *prefix_, int32 order_, char const *label_, int32 tn_) {
    if (label_ == NULL) {
      file      = stderr;
      prefix[0] = 0;
      name[0]   = 0;
      part      = 0;
      length    = 0;
      return;
    }

    snprintf(prefix, FILENAME_MAX, "%s.%03u.%s",         prefix_, order_, label_);
    snprintf(name, FILENAME_MAX,   "%s.%03u.%s.thr%03d", prefix_, order_, label_, tn_);
  };

  void  rotate(void) {

    assert(name[0] != 0);

    AS_UTL_closeFile(file, name);

    file   = NULL;
    length = 0;

    part++;
  }

  void  open(void) {
    char    path[FILENAME_MAX];

    assert(file == NULL);
    assert(name[0] != 0);

    snprintf(path, FILENAME_MAX, "%s.num%03d.log", name, part);

    errno = 0;
    file = fopen(path, "w");
    if (errno) {
      writeStatus("setLogFile()-- Failed to open logFile '%s': %s.\n", path, strerror(errno));
      writeStatus("setLogFile()-- Will now log to stderr instead.\n");
      file = stderr;
    }
  };

  void  close(void) {
    AS_UTL_closeFile(file, name);

    file      = NULL;
    prefix[0] = 0;
    name[0]   = 0;
    part      = 0;
    length    = 0;
  };

  FILE   *file;
  char    prefix[FILENAME_MAX];
  char    name[FILENAME_MAX];
  uint32  part;
  uint64  length;
};


//  NONE of the logFileMain/logFileThread is implemented


logFileInstance    logFileMain;           //  For writes during non-threaded portions
logFileInstance   *logFileThread = NULL;  //  For writes during threaded portions.
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
                                     NULL
};

//  Closes the current logFile, opens a new one called 'prefix.logFileOrder.label'.  If 'label' is
//  NULL, the logFile is reset to stderr.
void
setLogFile(char const *prefix, char const *label) {

  assert(prefix != NULL);

  //  Allocate space.

  if (logFileThread == NULL)
    logFileThread = new logFileInstance [omp_get_max_threads()];

  //  If writing to stderr, that's all we needed to do.

  if (logFileFlagSet(LOG_STDERR))
    return;

  //  Close out the old.

  logFileMain.close();

  for (int32 tn=0; tn<omp_get_max_threads(); tn++)
    logFileThread[tn].close();

  //  Move to the next iteration.

  logFileOrder++;

  //  Set up for that iteration.

  logFileMain.set(prefix, logFileOrder, label, 0);

  for (int32 tn=0; tn<omp_get_max_threads(); tn++)
    logFileThread[tn].set(prefix, logFileOrder, label, tn+1);

  //  File open is delayed until it is used.

}



char *
getLogFilePrefix(void) {
  return(logFileMain.prefix);
}



void
writeStatus(char const *fmt, ...) {
  va_list           ap;

  va_start(ap, fmt);

  vfprintf(stderr, fmt, ap);

  va_end(ap);
}



void
writeLog(char const *fmt, ...) {
  va_list           ap;
  int32             nt = omp_get_num_threads();
  int32             tn = omp_get_thread_num();

  logFileInstance  *lf = (nt == 1) ? (&logFileMain) : (&logFileThread[tn]);

  //  Rotate the log file please, HAL.

  uint64  maxLength = 512 * 1024 * 1024;

  if ((lf->name[0] != 0) &&
      (lf->length  > maxLength)) {
    fprintf(lf->file, "logFile()--  size " F_U64 " exceeds limit of " F_U64 "; rotate to new file.\n",
            lf->length, maxLength);
    lf->rotate();
  }

  //  Open the file if needed.

  if (lf->file == NULL)
    lf->open();

  //  Write the log.

  va_start(ap, fmt);

  lf->length += vfprintf(lf->file, fmt, ap);

  va_end(ap);
}



void
flushLog(void) {
  int32             nt = omp_get_num_threads();
  int32             tn = omp_get_thread_num();

  logFileInstance  *lf = (nt == 1) ? (&logFileMain) : (&logFileThread[tn]);

  if (lf->file != NULL)
    fflush(lf->file);
}
