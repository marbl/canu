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
#include "sweatShop.H"

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





//  XXX  Both loader and loaderAll leave the last gen sequence undeleted!







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

u32bit            loaderCacheSize  = 1024;        //  Size of the EST cache

//  Things that the various threads need -- we should put these into a
//  sim4thGlobal class, but we already have a bunch of globals
//  anyway....
//
readBuffer           *scriptFile       = 0L;

sim4parameters        sim4params;

FastAWrapper         *GENs             = 0L;
FastACache           *ESTs             = 0L;

u32bit                lastGENiid       = ~u32bitZERO;
u32bit                lastESTiid       = ~u32bitZERO;
FastASequenceInCore  *lastGENseq       = 0L;

int                   fOutput          = 0;



class sim4thWork {
public:
  sim4command          *input;
  char                 *script;
  sim4polishList       *output;
  FastASequenceInCore  *gendelete;
  FastASequenceInCore  *estdelete;

  sim4thWork() {
    input = 0L;
    script = 0L;
    output = 0L;
    gendelete = 0L;
    estdelete = 0L;
  };
};




void*
loader(void *U) {
  bool                  doForward = true;
  bool                  doReverse = true;
  u32bit                ESTiid = 0;
  u32bit                GENiid = 0;
  u32bit                GENlo  = 0;
  u32bit                GENhi  = 0;

  sim4thWork *p = new sim4thWork();
  
  p->script = getNextScript(ESTiid, GENiid, GENlo, GENhi, doForward, doReverse, scriptFile);

  if (p->script) {
    FastASequenceInCore  *ESTseq = 0L;
    FastASequenceInCore  *GENseq = 0L;

    //  If we already have the GENseq, use that, otherwise, register it for deletion.
    //
    if (lastGENiid == GENiid) {
      GENseq = lastGENseq;
    } else {

      //  Register it for deletion.  Technically, we're deleting this
      //  on the state AFTER it's used, but we can't guarantee that
      //  that state is still around.  The writer is deleting this, so
      //  by the time it gets here, it already wrote everyone that
      //  used this, which kind of implies that everyone that needs
      //  this is already computed.
      //
      p->gendelete = lastGENseq;

      GENs->find(GENiid);
      GENseq = GENs->getSequence();

      lastGENiid = GENiid;
      lastGENseq = GENseq;
    }

    //  The cache can, and does, overwrite the EST sequence we care
    //  about.  For now, we just copy the EST from the cache.
    //
    ESTseq         = ESTs->getSequence(ESTiid)->copy();
    p->estdelete   = ESTseq;

    p->input       = new sim4command(ESTseq, GENseq, GENlo, GENhi, doForward, doReverse);
  } else {
    delete p;
    p = 0L;
  }

  return(p);
}



void*
loaderAll(void *) {

  sim4thWork  *p = new sim4thWork();

  //  Previous implementations "Ping-pong'd" through the ESTs.  The
  //  idea being we would use the cache on the ends.  We can't easily
  //  do that here, so we always go forward.

  //  Flip around the end, if needed.
  if (lastESTiid >= ESTs->fasta()->getNumberOfSequences()) {
    lastESTiid   = 0;
    p->gendelete = lastGENseq;
    lastGENseq   = 0L;

    if (lastGENiid == ~u32bitZERO)  //  happens on the first time through
      lastGENiid = 0;
    else
      lastGENiid++;
  }

  //  If we've run out of sequences, we're done!
  if (lastGENiid >= GENs->getNumberOfSequences()) {
    delete p;
    return(0L);
  }

  //  Update the genomic sequence?
  if (lastGENseq == 0L) {
    GENs->find(lastGENiid);
    lastGENseq = GENs->getSequence();
  }

  //  Grab the EST sequence
  p->estdelete = ESTs->getSequence(lastESTiid++)->copy();

  //  build the command
  p->input     = new sim4command(p->estdelete,
                                 lastGENseq, 0, lastGENseq->sequenceLength(),
                                 true, true);

  return(p);
}





void
worker(void *U, void *S) {
  sim4thWork  *p = (sim4thWork *)S;

  Sim4       *sim = new Sim4(&sim4params);
  p->output       = sim->run(p->input);
  delete sim;
}


void
writer(void *U, void *S) {
  sim4thWork  *p = (sim4thWork *)S;

  sim4polishList  &L4 = *(p->output);

  for (u32bit i=0; L4[i]; i++) {
    char *o = s4p_polishToString(L4[i]);

    errno = 0;
    write(fOutput, o, strlen(o) * sizeof(char));
    if (errno)
      fprintf(stderr, "Couldn't write the output file '%s': %s\n", outputFileName, strerror(errno)), exit(1);

    free(o);
  }

  if (beYesNo) {
    if (L4[0])
      fprintf(stdout, "%s -Y "u32bitFMT" "u32bitFMT"\n",
              p->script,
              L4[0]->percentIdentity,
              L4[0]->querySeqIdentity);
    else
      fprintf(stdout, "%s -N 0 0\n",
              p->script);
  }

  //  Release this compute

  delete    p->input;
  delete [] p->script;
  delete    p->output;
  delete    p->gendelete;
  delete    p->estdelete;
  delete    p;
}








int
openOutputFile(char *outputFileName) {
  int  f = 0;

  if (strcmp(outputFileName, "-") == 0) {
    f = fileno(stdout);
  } else {
    errno = 0;
    f = open(outputFileName,
             O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
             S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "Couldn't open the output file '%s': %s\n", outputFileName, strerror(errno)), exit(1);
  }
  return(f);
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

  //  Open the output file
  fOutput = openOutputFile(outputFileName);

  sweatShop  *ss = 0L;

  //  If we have a script, read work from there, otherwise,
  //  do an all-vs-all.
  //
  if (scriptFileName) {
    scriptFile = new readBuffer(scriptFileName);
    ss = new sweatShop(loader,
                       worker,
                       writer);
  } else {
    //ss = new sweatShop(loader,
    //                   workerAll,
    //                   writer);
  }

  ss->run();

  //  Only close the file if it isn't stdout
  //
  if (strcmp(outputFileName, "-") != 0)
    close(fOutput);

  delete scriptFile;

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

  exit(0);
}
