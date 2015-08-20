
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
 *    Brian P. Walenz from 2006-JUL-24 to 2010-NOV-02
 *      are Copyright 2006,2010 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#include <pthread.h>
#include <semaphore.h>

//  Tests if add is atomic.

long int     count     = 0;
long int     counts[8] = { 0 };
pthread_t    threadID[8];

double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}

void*
workerThread(void *idx) {
  long int   &c = (*(long int *)idx);

  double  start = getTime();
  while (getTime() - start < 5) {
    c++;  count++;
    c++;  count++;
    c++;  count++;
    c++;  count++;
    c++;  count++;
    c++;  count++;
  }
}

int
main(int argc, char **argv) {
  pthread_attr_t   threadAttr;

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  int  numThreads = 5;
  int  sum        = 0;

  for (int i=0; i<numThreads; i++)
    pthread_create(threadID+i, &threadAttr, workerThread, (void *)(counts + i));

  for (int i=0; i<numThreads; i++) {
    pthread_join(threadID[i], 0L);
    fprintf(stderr, "thread %2ld  %ld\n", i, counts[i]);
    sum += counts[i];
  }

  fprintf(stderr, "total sum: %ld  total global: %ld\n", sum, count);
}


