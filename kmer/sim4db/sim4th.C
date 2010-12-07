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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <signal.h>
#include <math.h>

#include <pthread.h>
#include <semaphore.h>

#include "bio++.H"
#include "sim4.H"
#include "sweatShop.H"

//  XXX  Both loader and loaderAll leave the last gen sequence undeleted!

readBuffer       *scriptFile       = 0L;

seqCache         *GENs             = 0L;
seqCache         *ESTs             = 0L;

u32bit            lastGENiid       = ~u32bitZERO;
u32bit            lastESTiid       = ~u32bitZERO;
seqInCore        *lastGENseq       = 0L;

int               fOutput          = 0;
int               fYesNo           = 0;

char             *cdnaFileName     = 0L;
char             *scriptFileName   = 0L;
char             *databaseFileName = 0L;
char             *outputFileName   = 0L;
char             *yesnoFileName    = 0L;
char             *touchFileName    = 0L;

bool              pairwise         = false;

bool              beVerbose        = false;
bool              beYesNo          = false;

u32bit            numThreads       = 2;
u32bit            loaderCacheSize  = 1024;

sim4parameters    sim4params;

//  Parse the command line to create a sim4command object
//
//  [-f|-r] -e ESTid -D GENid GENlo GENhi
//
//    -f  Forward only
//    -r  Reverse only
//    -D  genSeqIID genLo genHi
//    -e  estSeqIID
//
//
char*
getNextScript(u32bit     &ESTiid,
              u32bit     &GENiid, u32bit &GENlo, u32bit &GENhi,
              bool       &doForward,
              bool       &doReverse) {

  char x = scriptFile->read();

  //  Skip any white space in the file
  //
  while ((scriptFile->eof() == false) && (whitespaceSymbol[x]))
    x = scriptFile->read();

  //  Exit if we're all done.
  //
  if (scriptFile->eof())
    return(0L);

  u32bit  linePos = 0;
  u32bit  lineMax = 128;
  char   *line    = new char [lineMax];

  //  Copy the line from the readBuffer into our storage
  //
  while ((scriptFile->eof() == false) && (x != '\n')) {
    line[linePos++] = x;
    x = scriptFile->read();
  }
  line[linePos] = 0;

  //  Decode the line
  //
  u32bit         argWords = 0;
  splitToWords   words(line);

  while (words.getWord(argWords)) {
    switch (words.getWord(argWords)[1]) {
      case 'f':
        doForward = true;
        doReverse = false;
        break;
      case 'r':
        doForward = false;
        doReverse = true;
        break;
      case 'D':
        GENiid = strtou32bit(words.getWord(++argWords), 0L);
        GENlo  = strtou32bit(words.getWord(++argWords), 0L);
        GENhi  = strtou32bit(words.getWord(++argWords), 0L);
        break;
      case 'e':
        ESTiid = strtou32bit(words.getWord(++argWords), 0L);
        break;
      default:
        //fprintf(stderr, "Unknown option '%s'\n", words.getWord(argWords));
        break;
    }

    argWords++;
  }

  return(line);
}




class sim4thWork {
public:
  sim4command          *input;
  char                 *script;
  sim4polishList       *output;
  seqInCore            *gendelete;
  seqInCore            *estdelete;

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
  
  p->script = getNextScript(ESTiid, GENiid, GENlo, GENhi, doForward, doReverse);

  if (p->script) {
    seqInCore  *ESTseq = 0L;
    seqInCore  *GENseq = 0L;

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

      GENseq = GENs->getSequenceInCore(GENiid);

      lastGENiid = GENiid;
      lastGENseq = GENseq;
    }

    //  The cache can, and does, overwrite the EST sequence we care
    //  about.  For now, we just copy the EST from the cache.
    //
    ESTseq         = ESTs->getSequenceInCore(ESTiid)->copy();
    p->estdelete   = ESTseq;

    p->input       = new sim4command(ESTseq, GENseq, GENlo, GENhi, doForward, doReverse);
  } else {
    delete p;
    p = 0L;
  }

  return(p);
}


void*
loaderPairwise(void *) {

  //  Align cDNA i to genomic i.

  if (lastGENiid == ~u32bitZERO)  //  happens on the first time through
    lastGENiid = 0;
  if (lastESTiid == ~u32bitZERO)  //  happens on the first time through
    lastESTiid = 0;

  //  If we've run out of sequences, we're done!
  if ((lastGENiid >= GENs->getNumberOfSequences()) ||
      (lastESTiid >= ESTs->getNumberOfSequences()))
    return(0L);

  sim4thWork  *p = new sim4thWork();

  //  Grab the GEN sequence
  p->gendelete = GENs->getSequenceInCore(lastGENiid++);

  //  Grab the EST sequence
  p->estdelete = ESTs->getSequenceInCore(lastESTiid++)->copy();

  //  build the command
  p->input     = new sim4command(p->estdelete,
                                 p->gendelete, 0, p->gendelete->sequenceLength(),
                                 true, true);

  return(p);
}


void*
loaderAll(void *) {

  sim4thWork  *p = new sim4thWork();

  //  Previous implementations "Ping-pong'd" through the ESTs.  The
  //  idea being we would use the cache on the ends.  We can't easily
  //  do that here, so we always go forward.

  //  Flip around the end, if needed.
  if (lastESTiid >= ESTs->getNumberOfSequences()) {
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
    lastGENseq = GENs->getSequenceInCore(lastGENiid);
  }

  //  Grab the EST sequence
  p->estdelete = ESTs->getSequenceInCore(lastESTiid++)->copy();

  //  build the command
  p->input     = new sim4command(p->estdelete,
                                 lastGENseq, 0, lastGENseq->sequenceLength(),
                                 true, true);

  return(p);
}




void
worker(void *U, void *T, void *S) {
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
    char *o = L4[i]->s4p_polishToString(sim4params.getOutputFormat());

    errno = 0;
    write(fOutput, o, strlen(o) * sizeof(char));
    if (errno)
      fprintf(stderr, "Couldn't write the output file '%s': %s\n", outputFileName, strerror(errno)), exit(1);

    delete [] o;
  }

  if (yesnoFileName) {
    char  str[128];

    if (L4[0])
      sprintf(str, "%s -Y "u32bitFMT" "u32bitFMT"\n",
              p->script, L4[0]->_percentIdentity, L4[0]->_querySeqIdentity);
    else
      sprintf(str, "%s -N 0 0\n", p->script);

    write(fYesNo, str, strlen(str) * sizeof(char));
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
openOutputFile(char *name) {
  int  f = 0;

  if (name == 0L)
    return(0);

  if (strcmp(name, "-") == 0) {
    f = fileno(stdout);
  } else {
    errno = 0;
    f = open(name,
             O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
             S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "Couldn't open the output file '%s': %s\n", name, strerror(errno)), exit(1);
  }
  return(f);
}



int
main(int argc, char **argv) {

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-alignments", 4) == 0) {
      sim4params.setPrintAlignments(true);

    } else if (strncmp(argv[arg], "-alwaysprint", 4) == 0) {
      sim4params.setFindAllExons(true);
      sim4params.setAlwaysReport(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-cdna", 3) == 0) {
      cdnaFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-cut", 3) == 0) {
      double x = atof(argv[++arg]);
      if (x < 0.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 0.0 (you gave %f)!\n", x);
        x = 0.0;
      }
      if (x > 1.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 1.0 (you gave %f)!\n", x);
        x = 1.0;
      }
      sim4params.setPolyTailPercent(x);

    } else if (strncmp(argv[arg], "-genomic", 4) == 0) {
      databaseFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-minc", 5) == 0) {
      sim4params.setFindAllExons(true);
      sim4params.setMinCoverage(atoi(argv[++arg]) / 100.0);

    } else if (strncmp(argv[arg], "-mini", 5) == 0) {
      sim4params.setFindAllExons(true);
      sim4params.setMinPercentExonIdentity(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-minl", 5) == 0) {
      sim4params.setFindAllExons(true);
      sim4params.setMinCoverageLength(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-nod", 4) == 0) {
      sim4params.setIncludeDefLine(false);

    } else if (strncmp(argv[arg], "-non", 4) == 0) {
      sim4params.setDontForceCanonicalSplicing(true);

    } else if (strncmp(argv[arg], "-f", 2) == 0) {
      sim4params.setForceStrandPrediction(true);

    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      outputFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-po", 3) == 0) {
      sim4params.setIgnorePolyTails(false);

    } else if (strncmp(argv[arg], "-sc", 3) == 0) {
      scriptFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-sp", 3) == 0) {
      sim4params.setSpliceModel(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-pa", 3) == 0) {
      pairwise = true;

    } else if (strncmp(argv[arg], "-to", 3) == 0) {
      touchFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;

    } else if (strncmp(argv[arg], "-YN", 3) == 0) {
      yesnoFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-threads", 3) == 0) {
      numThreads = strtou32bit(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-H", 2) == 0) {
      sim4params.setRelinkWeight(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-K", 2) == 0) {
      sim4params.setMSPThreshold1(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-C", 2) == 0) {
      sim4params.setMSPThreshold2(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-Z", 2) == 0) {
      sim4params.setSpacedSeed(argv[++arg]);

    } else if (strncmp(argv[arg], "-Ma", 3) == 0) {
      sim4params.setMSPLimitAbsolute(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-Mp", 3) == 0) {
      sim4params.setMSPLimitPercent(atof(argv[++arg]));

    } else if (strncmp(argv[arg], "-interspecies", 2) == 0) {
      sim4params.setInterspecies(true);

    } else if (strcmp(argv[arg], "-gff3") == 0) {
      sim4params.setOutputFormat(S4P_POLISH_GFF3);

    } else {
      fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if ((err) ||
      (cdnaFileName == 0L) ||
      (databaseFileName == 0L) ||
      (outputFileName == 0L)) {
    fprintf(stderr, "usage: %s -genomic g.fasta -cdna c.fasta -output o.sim4db [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "       -v            print status to stderr while running\n");
    fprintf(stderr, "       -V            print script lines (stderr) as they are processed\n");
    fprintf(stderr, "       -YN           print script lines (to given file) as they are processed, annotated with yes/no\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -cdna         use these cDNA sequences\n");
    fprintf(stderr, "       -genomic      use these genomic sequences\n");
    fprintf(stderr, "       -script       use this script file\n");
    fprintf(stderr, "       -pairwise     do pairs of sequences\n");
    fprintf(stderr, "       -output       write output to this file\n");
    fprintf(stderr, "       -touch        create this file when the program finishes execution\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -threads      Use n threads.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -mincoverage  iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum PERCENT COVERAGE\n");
    fprintf(stderr, "       -minidentity  iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum PERCENT EXON IDENTITY\n");
    fprintf(stderr, "       -minlength    iteratively find all exon models with the specified\n");
    fprintf(stderr, "                     minimum ABSOLUTE COVERAGE (number of bp matched)\n");
    fprintf(stderr, "       -alwaysreport always report <number> exon models, even if they\n");
    fprintf(stderr, "                     are below the quality thresholds\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         If no mincoverage or minidentity or minlength is given, only\n");
    fprintf(stderr, "         the best exon model is returned.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         You will probably want to specify ALL THREE of mincoverage,\n");
    fprintf(stderr, "         minidentity and minlength!  Don't assume the default values\n");
    fprintf(stderr, "         are what you want!\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "         You will DEFINITELY want to specify at least one of mincoverage,\n");
    fprintf(stderr, "         minidentity and minlength with alwaysreport!  If you don't, mincoverage\n");
    fprintf(stderr, "         will be set to 90 and minidentity to 95 -- to reduce the number of\n");
    fprintf(stderr, "         spurious matches when a good match is found.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -nodeflines   don't include the defline in the output\n");
    fprintf(stderr, "       -alignments   print alignments\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -polytails    DON'T mask poly-A and poly-T tails.\n");
    fprintf(stderr, "       -cut          Trim marginal exons if A/T %% > x (poly-AT tails)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -noncanonical Don't force canonical splice sites\n");
    fprintf(stderr, "       -splicemodel  Use the following splice model: 0 - original sim4;\n");
    fprintf(stderr, "                     1 - GeneSplicer; 2 - Glimmer (default: 0)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -forcestrand  Force the strand prediction to always be\n");
    fprintf(stderr, "                     'forward' or 'reverse'\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "       -interspecies Use sim4cc for inter-species alignments\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  The following are for use only by immortals.\n");
    fprintf(stderr, "       -Z            set the (spaced) seed pattern\n");
    fprintf(stderr, "       -H            set the relink weight factor\n");
    fprintf(stderr, "       -K            set the first MSP threshold\n");
    fprintf(stderr, "       -C            set the second MSP threshold\n");
    fprintf(stderr, "       -Ma           set the limit of the number of MSPs allowed\n");
    fprintf(stderr, "       -Mp           same, as percentage of bases in cDNA\n");
    fprintf(stderr, "                     NOTE:  If used, both -Ma and -Mp must be specified!\n");
    exit(1);
  }

  //  Open input files
  //
  GENs = new seqCache(databaseFileName);
  ESTs = new seqCache(cdnaFileName, loaderCacheSize, false);

  //  Open the output file
  fOutput = openOutputFile(outputFileName);
  fYesNo  = openOutputFile(yesnoFileName);

  sweatShop  *ss = 0L;

  err = sim4params.setSpliceMutex();
  if (err) {
    fprintf(stderr, "sim4th::main()--  Failed to initialize splice mutex: %s.\n", strerror(err));
    exit(1);
  }


  //  If we have a script, read work from there, otherwise,
  //  do an all-vs-all.
  //
  if (scriptFileName) {
    scriptFile = new readBuffer(scriptFileName);
    ss = new sweatShop(loader,
                       worker,
                       writer);
  } else if (pairwise) {
    ss = new sweatShop(loaderPairwise,
                       worker,
                       writer);
  } else {
    ss = new sweatShop(loaderAll,
                       worker,
                       writer);
  }

  ss->setNumberOfWorkers(numThreads);
  ss->run(0L, beVerbose);

  //  Only close the file if it isn't stdout
  //
  if (strcmp(outputFileName, "-") != 0)
    close(fOutput);

  if (yesnoFileName)
    close(fYesNo);

  delete scriptFile;


  exit(0);
}
