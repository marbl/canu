// This file is part of sim4db.
// Copyright (c) 2005 Brian Walenz
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "sim4db.H"

//  A threaded sim4db implementation
//
//  Unthreaded
//  1101.581u 198.865s 28:57.02 74.8%       281+5357k 294664+10116io 28pf+0w
//
//  Threaded, two threads, busy wait on workers
//  1574.594u 216.661s 16:59.35 175.7%      286+-3135k 503909+10116io 19pf+0w
//
//  8 threads, possibly busy wait on workers
//  1649.177u 205.446s 17:28.67 176.8%      285+-649k 509999+10116io 25pf+0w
//
//  3 threads, workers sleep, 64k cache, 16k queue
//  1232.354u 237.168s 16:22.09 149.6%      288+-13241k 509893+10116io 24pf+0w



//  Define this to enable the loader to submit things to the queue in
//  batches, hopefully, will reduce the contention on the lock.
//
#define LOAD_BATCHES


//  We keep a count of the number of times the workers try to access
//  the mutex but it is locked.
//
//  This is not an issue, even disabling LOAD_BATCHES (waste of time
//  that was!)
//
//#define COUNT_LOCKED



//  Run options set from the command line.
//
char             *cdnaFileName     = 0L;
char             *scriptFileName   = 0L;
char             *databaseFileName = 0L;
char             *outputFileName   = 0L;
char             *statsFileName    = 0L;
char             *touchFileName    = 0L;

bool              beVerbose        = false;       //  Print progress
bool              beYesNo          = false;       //  Print each script line as we process, with answer

sim4parameters    sim4params;


//  We used to use five arrays holding the state for each compute.  all-vs-all
//  usually blew out the memory on just those arrays, so we switched over
//  to a linked list of structs.
//
struct state_s {
public:
  sim4command          *input;
  char                 *script;
  sim4polishList       *output;
  FastASequenceInCore  *gendelete;
  FastASequenceInCore  *estdelete;
  state_s              *next;

  state_s() {
    input = 0L;
    script = 0L;
    output = 0L;
    gendelete = 0L;
    estdelete = 0L;
    next = 0L;
  };
};

pthread_mutex_t        stateMutex;
state_s               *rootState = 0L;  //  Where output takes stuff from, the tail
state_s               *currState = 0L;  //  Where computes happen, the middle
state_s               *headState = 0L;  //  Where input is put, the head

u32bit                 numberLoaded   = 0;
u32bit                 numberComputed = 0;
u32bit                 numberOutput   = 0;

#ifdef COUNT_LOCKED
u32bit            mutexLocked = 0;
#endif

state_s*
newState(void) {
  return(new state_s);
}

void
freeState(state_s *s) {
  delete s;
}

FastAWrapper     *GENs             = 0L;
FastACache       *ESTs             = 0L;

//  I see 500/sec on a dual P4 Xeon 2.8 with local disk, so
//  loaderQueueSize gives us 16 seconds of buffer.  ESTs are ~500 bp
//  long, keeping 256,000 in the cache will eat up 128MB.  I can't
//  characterize what types of sequences will see benefits from
//  increasing or decreasing this.
//
//  Increased loaderQueueSize to make sure that we don't hit a bad
//  region in the middle of a run and empty the queue.
//
u32bit            loaderQueueSize  =  64 * 1024;
u32bit            loaderCacheSize  = 128 * 1024;


//  Add a new state to the head of the queue.  This is just gross.
//
void
loaderAdd(state_s *thisState, sim4command *cmd) {
  pthread_mutex_lock(&stateMutex);
  thisState->input = cmd;
  thisState->next  = 0L;
  if (headState == 0L) {
    rootState      = thisState;
    currState      = thisState;
    headState      = thisState;
  } else {
    headState->next = thisState;
  }
  headState        = thisState;
  pthread_mutex_unlock(&stateMutex);
  numberLoaded++;
}

#ifdef LOAD_BATCHES

//  Build a list of states to add in one swoop
//
void
loaderSave(state_s *&tail, state_s *&head, state_s *thisState, sim4command *cmd) {

  thisState->input = cmd;
  thisState->next  = 0L;

  if (tail) {
    tail->next = thisState;
    head       = thisState;
  } else {
    tail = head = thisState;
  }
  numberLoaded++;
}

//  Add a bunch of new states to the queue.
//
void
loaderAppend(state_s *&tail, state_s *&head) {
  pthread_mutex_lock(&stateMutex);

  if (headState == 0L) {
    rootState      = tail;
    currState      = tail;
    headState      = head;
  } else {
    headState->next = tail;
  }
  headState        = head;

  pthread_mutex_unlock(&stateMutex);

  tail = 0L;
  head = 0L;
}

#endif




void*
loader(void *) {
  u32bit                lastGENiid = ~u32bitZERO;
  FastASequenceInCore  *lastGENseq = 0L;

  //  We can batch several loads together before we push them onto the
  //  queue, this should reduce the number of times the loader needs to
  //  lock the queue.
  //
#ifdef LOAD_BATCHES
  state_s               *tail = 0L;  //  The first thing loaded
  state_s               *head = 0L;  //  The last thing loaded
  u32bit                 batchSize  = 0;
  u32bit                 batchLimit = 128;
#endif

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 10000000;

  //  Open our script file
  //
  readBuffer  *scriptFile = new readBuffer(scriptFileName);

  bool  moreToLoad = true;

  while (moreToLoad) {

    //  Zzzzzzz....
    while (numberLoaded - numberComputed > loaderQueueSize)
      nanosleep(&naptime, 0L);

    bool                  doForward = true;
    bool                  doReverse = true;
    u32bit                ESTiid = 0;
    u32bit                GENiid = 0;
    u32bit                GENlo  = 0;
    u32bit                GENhi  = 0;

    state_s  *thisState = newState();
    thisState->script = getNextScript(ESTiid, GENiid, GENlo, GENhi, doForward, doReverse, scriptFile);

    if (thisState->script) {
      FastASequenceInCore  *ESTseq = 0L;
      FastASequenceInCore  *GENseq = 0L;

      //  If we already have the GENseq, use that, otherwise, register it for deletion.
      //
      if (lastGENiid == GENiid) {
        GENseq = lastGENseq;
      } else {
        //  Technically, we're deleting this on the state AFTER it's
        //  used, but we can't guarantee that that state is still
        //  around.
        //
        thisState->gendelete = lastGENseq;

        GENs->find(GENiid);
        GENseq = GENs->getSequence();

        lastGENiid = GENiid;
        lastGENseq = GENseq;
      }

      //  For now, we just copy the EST from the cache.
      //
      ESTseq               = ESTs->getSequence(ESTiid)->copy();
      thisState->estdelete = ESTseq;

#ifdef LOAD_BATCHES
      loaderSave(tail, head, thisState, new sim4command(ESTseq, GENseq, GENlo, GENhi, doForward, doReverse));
      batchSize++;
      if (batchSize >= batchLimit)
        loaderAppend(tail, head);
#else
      loaderAdd(thisState, new sim4command(ESTseq, GENseq, GENlo, GENhi, doForward, doReverse));
#endif
      thisState = newState();
    } else {
      //  Didn't read, must be all done!  Push on the end-of-input marker state.
      //
#ifdef LOAD_BATCHES
#endif
      loaderAdd(newState(), 0L);
      moreToLoad = false;
    }
  }

  delete scriptFile;

  return(0L);
}



void*
loaderAll(void *) {
  FastASequenceInCore  *ESTseq = 0L;
  FastASequenceInCore  *GENseq = 0L;

  bool                  doForward = true;
  bool                  doReverse = true;

  u32bit                numESTs = ESTs->fasta()->getNumberOfSequences();

  //  It's super easy for us to fill up the queue, so the delay here
  //  is pretty large.
  //
  struct timespec   naptime;
  naptime.tv_sec      = 1;
  naptime.tv_nsec     = 0;

  //  We ping-pong through the ESTs.  At least we can use the cache a little
  //  bit on the ends.
  //
  //  Of course, this really REALLY bites us for the reverse when
  //  we're not in the cache.

  for (u32bit GENiid=0; GENiid < GENs->getNumberOfSequences(); GENiid++) {
    GENs->find(GENiid);
    GENseq = GENs->getSequence();

    for (u32bit ESTiid=0; ESTiid < numESTs; ESTiid++) {
      while (numberLoaded - numberComputed > loaderQueueSize)
        nanosleep(&naptime, 0L);

      state_s  *thisState = newState();

      ESTseq           = ESTs->getSequence(ESTiid)->copy();
      thisState->estdelete = ESTseq;

      if (ESTiid == numESTs-1)
        thisState->gendelete = GENseq;

      loaderAdd(thisState, new sim4command(ESTseq,
                                           GENseq, 0, GENseq->sequenceLength(),
                                           doForward, doReverse));
    }

    //  Ping pong!

    GENiid++;
    if (GENiid == GENs->getNumberOfSequences())
      return(0L);

    GENs->find(GENiid);
    GENseq = GENs->getSequence();

    for (u32bit ESTiid=numESTs; ESTiid-- > 0;) {
      while (numberLoaded - numberComputed > loaderQueueSize)
        nanosleep(&naptime, 0L);

      state_s  *thisState = newState();

      ESTseq           = ESTs->getSequence(ESTiid)->copy();
      thisState->estdelete = ESTseq;

      if (ESTiid == 0)
        thisState->gendelete = GENseq;

      loaderAdd(thisState, new sim4command(ESTseq,
                                           GENseq, 0, GENseq->sequenceLength(),
                                           doForward, doReverse));
    }
  } 

  return(0L);
}





void*
worker(void *) {
  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 500000000;

  while (currState && currState->input) {

    //  Grab the next state.  We don't grab it if it's the last in the
    //  queue (else we would fall off the end) UNLESS it really is the
    //  last one.
    //
    state_s  *thisState = 0L;

#ifdef COUNT_LOCKED
    if (pthread_mutex_trylock(&stateMutex)) {
      //  mutex already locked, wait for it.
      mutexLocked++;
      pthread_mutex_lock(&stateMutex);
    }
#else
    pthread_mutex_lock(&stateMutex);
#endif

    if ((currState) && ((currState->next != 0L) || (currState->input == 0L))) {
      thisState = currState;
      currState = currState->next;
    }
    pthread_mutex_unlock(&stateMutex);

    //  Execute
    //
    if (thisState && thisState->input) {
      Sim4            *sim = new Sim4(&sim4params);
      thisState->output    = sim->run(thisState->input);
      delete sim;
      numberComputed++;
    } else {
      nanosleep(&naptime, 0L);
    }
  }

  return(0L);
}





int
openOutputFile(char *outputFileName) {
  int  fOutput = 0;

  if (strcmp(outputFileName, "-") == 0) {
    fOutput = fileno(stdout);
  } else {
    errno = 0;
    fOutput = open(outputFileName,
                   O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "Couldn't open the output file '%s': %s\n", outputFileName, strerror(errno)), exit(1);
  }
  return(fOutput);
}


void*
stats(void *) {
  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 250000000ULL;

  double  startTime = getTime() - 0.001;

  while (rootState) {
#ifdef COUNT_LOCKED
    fprintf(stderr, " (%6.1f/s) (out="u32bitFMTW(8)" + "u32bitFMTW(8)" = cpu = "u32bitFMTW(8)" + "u32bitFMTW(8)" ("u32bitFMTW(8)" locked) = in = "u32bitFMTW(8)"\r",
            numberOutput / (getTime() - startTime),
            numberOutput,
            numberComputed - numberOutput,
            numberComputed,
            numberLoaded - numberComputed,
            mutexLocked,
            numberLoaded);
#else
    fprintf(stderr, " (%6.1f/s) (out="u32bitFMTW(8)" + "u32bitFMTW(8)" = cpu = "u32bitFMTW(8)" + "u32bitFMTW(8)" = in = "u32bitFMTW(8)"\r",
            numberOutput / (getTime() - startTime),
            numberOutput,
            numberComputed - numberOutput,
            numberComputed,
            numberLoaded - numberComputed,
            numberLoaded);
#endif
    fflush(stderr);
    nanosleep(&naptime, 0L);
  }
  return(0L);
}



int
main(int argc, char **argv) {

  double  mainStartTime = getTime();

  parseCommandLine(argc, argv);

  //  Open input files
  //
  GENs = new FastAWrapper(databaseFileName);
  ESTs = new FastACache(cdnaFileName,     loaderCacheSize, false);

  GENs->openIndex();

  //  Launch some threads
  //
  pthread_attr_t   threadAttr;
  pthread_t        threadID;

  pthread_mutex_init(&stateMutex, NULL);

  pthread_attr_init(&threadAttr);
  pthread_attr_setscope(&threadAttr, PTHREAD_SCOPE_SYSTEM);
  pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_DETACHED);
  pthread_attr_setschedpolicy(&threadAttr, SCHED_OTHER);

  //  Fire off the loader, or the all-vs-all loader
  pthread_create(&threadID, &threadAttr, (scriptFileName) ? loader : loaderAll, 0L);

  //  Wait for it to actually load something (otherwise all the threads exit)
  while (!rootState && !currState && !headState)
    sleep(1);

  //  And the stats
  pthread_create(&threadID, &threadAttr, stats, 0L);

  //  Give it some lead time to load stuff
  sleep(32);

  //  Fire off some workers
  pthread_create(&threadID, &threadAttr, worker, 0L);
  pthread_create(&threadID, &threadAttr, worker, 0L);
  pthread_create(&threadID, &threadAttr, worker, 0L);

  //  Open the output file
  int fOutput = openOutputFile(outputFileName);

  struct timespec   naptime;
  naptime.tv_sec      = 0;
  naptime.tv_nsec     = 100000000;

  //  Wait for output to appear.
  //
  while (rootState && rootState->input) {

    //  Wait a bit if there is no output on this node (we caught up to
    //  the computes), or there is no next node (we caught up to the
    //  input, yikes!)
    //
    if ((rootState->output == 0L) || (rootState->next == 0L)) {
      nanosleep(&naptime, 0L);
    } else {
      sim4polishList  &L4 = *(rootState->output);

      for (u32bit i=0; L4[i]; i++) {
        char *o = s4p_polishToString(L4[i]);

        errno = 0;
        //write(fOutput, o, strlen(o) * sizeof(char));
        if (errno)
          fprintf(stderr, "Couldn't write the output file '%s': %s\n", outputFileName, strerror(errno)), exit(1);

        free(o);
      }

      if (beYesNo) {
        if (L4[0])
          fprintf(stdout, "%s -Y "u32bitFMT" "u32bitFMT"\n",
                  rootState->script,
                  L4[0]->percentIdentity,
                  L4[0]->querySeqIdentity);
        else
          fprintf(stdout, "%s -N 0 0\n",
                  rootState->script);
      }

      state_s  *saveState = rootState->next;

      delete    rootState->input;
      delete [] rootState->script;
      delete    rootState->output;
      delete    rootState->gendelete;
      delete    rootState->estdelete;
      delete    rootState;

      rootState = saveState;

      numberOutput++;
    }
  }


  //  Only close the file if it isn't stdout
  //
  if (strcmp(outputFileName, "-") != 0)
    close(fOutput);

  if (statsFileName) {
    FILE  *statsFile = fopen(statsFileName, "w");
    if (statsFile) {
      write_rusage(statsFile);
      fprintf(statsFile, "clockTime:      %f\n", getTime() - mainStartTime);
      fclose(statsFile);
    }
  }

  if (touchFileName) {
    FILE  *touchFile = fopen(touchFileName, "w");
    fclose(touchFile);
  }

  return(0);
}
