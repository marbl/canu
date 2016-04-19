
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_Logging.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-JUL-29 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-19
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_Logging.H"

class logFileInstance {
public:
  logFileInstance() {
    file    = stderr;
    name[0] = 0;
    part    = 0;
    length  = 0;
  };
  ~logFileInstance() {
    if ((name[0] != 0) && (file)) {
      fprintf(stderr, "WARNING: open file '%s'\n", name);
      fclose(file);
    }
  };

  void  set(char const *prefix, int32 order, char const *label, int32 tn) {
    if (label == NULL) {
      file    = stderr;
      name[0] = 0;
      part    = 0;
      length  = 0;
      return;
    }

    sprintf(name, "%s.%03u.%s.thr%03d", prefix, order, label, tn);
  };

  void  rotate(void) {

    assert(name[0] != 0);

    fclose(file);

    file   = NULL;
    length = 0;

    part++;
  }

  void  open(void) {
    char    path[FILENAME_MAX];

    assert(file == NULL);
    assert(name[0] != 0);

    sprintf(path, "%s.num%03d.log", name, part);

    errno = 0;
    file = fopen(path, "w");
    if (errno) {
      fprintf(stderr, "setLogFile()-- Failed to open logFile '%s': %s.\n", path, strerror(errno));
      fprintf(stderr, "setLogFile()-- Will now log to stderr instead.\n");
      file = stderr;
    }
  };

  void  close(void) {
    if ((file != NULL) && (file != stderr))
      fclose(file);

    file    = NULL;
    name[0] = 0;
    part    = 0;
    length  = 0;
  };

  FILE   *file;
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
uint64 LOG_ALL_BEST_EDGES              = 0x0000000000000002;
uint64 LOG_CHUNK_GRAPH                 = 0x0000000000000004;  //  Report the chunk graph as we build it
uint64 LOG_BUILD_UNITIG                = 0x0000000000000008;  //  Report building of initial unitigs (both unitig creation and fragment placement)
uint64 LOG_PLACE_UNPLACED              = 0x0000000000000010;  //  Report placing of unplaced reads
uint64 LOG_BUBBLE_DETAIL               = 0x0000000000000020;
uint64 LOG_SPLIT_DISCONTINUOUS         = 0x0000000000000040;  //
uint64 LOG_INTERMEDIATE_UNITIGS        = 0x0000000000000080;  //  At various spots, dump the current unitigs
uint64 LOG_SET_PARENT_AND_HANG         = 0x0000000000000100;  //
uint64 LOG_STDERR                      = 0x0000000000000200;  //  Write ALL logging to stderr, not the files.

uint64 LOG_PLACE_FRAG                  = 0x8000000000000000;  //  Internal use only.

char const *logFileFlagNames[64] = { "overlapScoring",
                                     "allBestEdges",
                                     "chunkGraph",
                                     "buildUnitig",
                                     "placeUnplaced",
                                     "bubbles",
                                     "splitDiscontinuous",   //  Update made it to here, need repeats
                                     "intermediateUnitigs",
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

  if (label != NULL)
    fprintf(stderr,  "setLogFile()-- Now logging to '%s.%03d.%s'\n", prefix, logFileOrder, label);
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
    fprintf(lf->file, "logFile()--  size "F_U64" exceeds limit of "F_U64"; rotate to new file.\n",
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
