
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
 *    kmer/meryl/build-threads.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2004-APR-13 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2006-MAY-14 to 2014-APR-11
 *      are Copyright 2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"

void
runSegment(merylArgs *args, uint64 segment);

pthread_mutex_t        segmentMutex;
uint64                 segmentNext;
uint64                 segmentMax;
uint32                *segmentDone;


void*
buildThread(void *U) {
  uint64      segment = uint32ZERO;
  merylArgs  *args = (merylArgs *)U;

  while (segment < segmentMax) {
    pthread_mutex_lock(&segmentMutex);
    segment = segmentNext++;
    pthread_mutex_unlock(&segmentMutex);

    if (segment < segmentMax) {
      runSegment(args, segment);
      segmentDone[segment]++;
    }
  }

  if (args->beVerbose)
    fprintf(stderr, "Thread exits.\n");

  return(0L);
}


void
runThreaded(merylArgs *args) {

  //  Clear stuff
  //
  segmentNext = uint64ZERO;
  segmentMax  = args->segmentLimit;
  segmentDone = new uint32 [segmentMax];
  for (uint64 s=0; s<segmentMax; s++)
    segmentDone[s] = uint32ZERO;

  //  Initialize threads
  //
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_mutex_init(&segmentMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  //  Start the threads
  //
  for (uint64 i=0; i<args->numThreads; i++)
    pthread_create(&threadID, &threadAttr, buildThread, (void *)args);

  //  Wait for the threads to complete
  //
  struct timespec  shortNap;
  shortNap.tv_sec  = 1;
  shortNap.tv_nsec = 0;

  uint64 s=0;
  while (s < segmentMax) {
    if (segmentDone[s] == 0)
      nanosleep(&shortNap, 0L);
    else
      s++;
  }

  if (args->beVerbose)
    fprintf(stderr, "Threads all done, cleaning up.\n");

  //  Cleanup
  //
  pthread_attr_destroy(&threadAttr);
  pthread_mutex_destroy(&segmentMutex);

  delete [] segmentDone;
}
