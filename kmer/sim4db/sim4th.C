#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim4.H"
#include "libbri.H"
#include "fasta-cache.H"
#include "time.H"

#include <umalloc.h>

//  Simple threaded implementation.  For testing thread safety of sim4 itself.

sim4parameters  sim4params;
FastACache     *ESTs        = 0L;
FastACache     *GENs        = 0L;
u32bit          finished    = 0;

void*
theThread(void *U) {
  u32bit    workDone  = 0;

  FastASequenceInCore  *GEN = 0L;
  FastASequenceInCore  *EST = 0L;

  char path[80];
  sprintf(path, "out-%lu", (u64bit)U);
  int  fileout = fileno(fopen(path, "w"));

  for (u32bit genIdx=0; genIdx < GENs->fasta()->getNumberOfSequences(); genIdx++) {
    GEN = GENs->getSequence(genIdx);

    for (u32bit estIdx=0; estIdx < ESTs->fasta()->getNumberOfSequences(); estIdx++) {

      fprintf(stderr, "%lu At %u load\n", (u64bit)U, workDone);
      fflush(stderr);

      EST = ESTs->getSequence(estIdx);

      fprintf(stderr, "%lu At %u loaded\n", (u64bit)U, workDone); fflush(stderr);

      sim4command     *P4 = new sim4command(EST, GEN, 0, GEN->sequenceLength(), true, true);
      fprintf(stderr, "%lu At %u created command\n", (u64bit)U, workDone); fflush(stderr);

      Sim4            *S4 = new Sim4(&sim4params);
      fprintf(stderr, "%lu At %u created sim4\n", (u64bit)U, workDone); fflush(stderr);

      sim4polishList  *l4 = S4->run(P4);
      fprintf(stderr, "%lu At %u ran\n", (u64bit)U, workDone); fflush(stderr);

      delete S4;
      delete P4;

      workDone++;

#if 0
      for (u32bit i=0; L4[i]; i++) {
        char *o = s4p_polishToString(L4[i]);
        errno = 0;
        write(fileno(stdout), o, strlen(o) * sizeof(char));
        if (errno) {
          fprintf(stderr, "Couldn't write to stdout?\n");
          exit(1);
        }
        free(o);
      }
#endif

      delete l4;

      //  These are private to the cache!
      //delete EST;

      _dump_allocated_delta(fileout);
    }

    //  These are private to the cache!
    //delete GEN;
  }
  finished++;
  fprintf(stderr, "%lu At %u ALL DONE!\n", (u64bit)U, workDone); fflush(stderr);
  return(0L);
}


int
main(int argc, char **argv) {

  if (argc != 3) {
    fprintf(stderr, "usage: %s est-file gen-file\n", argv[0]);
    exit(1);
  }

  ESTs = new FastACache(argv[1], 0, true);
  GENs = new FastACache(argv[2], 0, true);

  fprintf(stderr, "After open cache\n");
  _dump_allocated_delta(fileno(stderr));
  fprintf(stderr, "----------------------------------------\n");

  pthread_attr_t   threadAttr;
  u32bit           threadIDX = 0;
  pthread_t        threadID[64];

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  fprintf(stderr, "After thread create\n");
  _dump_allocated_delta(fileno(stderr));
  fprintf(stderr, "----------------------------------------\n");

  for (u64bit i=0; i<16; i++)
    pthread_create(threadID + threadIDX++, &threadAttr, theThread, (void *)i);

  while (finished != 16)
    sleep(3);

  fprintf(stderr, "Main exits.\n");

  delete ESTs;
  delete GENs;

  return(0);
}
