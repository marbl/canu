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

#include "sim4.H"
#include "bri++.H"
#include "fasta.H"
#include "buildinfo-sim4db.h"
#include "buildinfo-libbri.h"
#include "buildinfo-libsim4.h"

#ifdef USE_CACHE
//  Using the FastACache for EST sequences has no effect on the enormous
//  system time we see in this code
#include "fasta-cache.H"
#endif


//  Run options set from the command line.
//
char           *cdnaFileName     = 0L;
char           *scriptFileName   = 0L;
char           *databaseFileName = 0L;
char           *outputFileName   = 0L;
char           *statsFileName    = 0L;
char           *touchFileName    = 0L;

bool            beVerbose        = false;       //  Print progress
bool            beExplicit       = false;       //  Print each script line as we process
bool            beYesNo          = false;       //  Print each script line as we process, with answer

sim4parameters  sim4params;

const char *usage = 
"usage: %s [options]\n"
"\n"
"       -v            print status to stderr while running\n"
"       -V            print script lines (stderr) as they are processed\n"
"       -YN           print script lines (stdout) as they are processed, annotated with yes/no\n"
"\n"
"       -cdna         use these cDNA sequences\n"
"       -genomic      use these genomic sequences\n"
"       -script       use this script file\n"
"       -output       write output to this file\n"
"       -stats        write execution statistics to this file\n"
"       -touch        create this file when the program finishes execution\n"
"\n"
"       -mincoverage  iteratively find all exon models with the specified\n"
"                     minimum PERCENT COVERAGE\n"
"       -minidentity  iteratively find all exon models with the specified\n"
"                     minimum PERCENT EXON IDENTITY\n"
"       -minlength    iteratively find all exon models with the specified\n"
"                     minimum ABSOLUTE COVERAGE (number of bp matched)\n"
"       -alwaysreport always report <number> exon models, even if they\n"
"                     are below the quality thresholds\n"
"\n"
"         If no mincoverage or minidentity or minlength is given, only\n"
"         the best exon model is returned.\n"
"\n"
"         You will probably want to specify ALL THREE of mincoverage,\n"
"         minidentity and minlength!  Don't assume the default values\n"
"         are what you want!\n"
"\n"
"         You will DEFINITELY want to specify at least one of mincoverage,\n"
"         minidentity and minlength with alwaysreport!  If you don't, mincoverage\n"
"         will be set to 90 and minidentity to 95 -- to reduce the number of\n"
"         spurious matches when a good match is found.\n"
"\n"
"       -nodeflines   don't include the defline in the output\n"
"       -alignments   print alignments\n"
"\n"
"       -polytails    DON'T mask poly-A and poly-T tails.\n"
"       -cut          Trim marginal exons if A/T %% > x (poly-AT tails)\n"
"\n"
"       -noncanonical Don't force canonical splice sites\n"
"\n"
"       -forcestrand  Force the strand prediction to always be\n"
"                     'forward' or 'reverse'\n"
"\n"
"       -interspecies Configure sim4 for better inter-species alignments\n"
"\n"
"  The following are for use only by immortals.\n"
"       -H            set the relink weight factor\n"
"       -K            set the first MSP threshold\n"
"       -C            set the second MSP threshold\n"
"       -Ma           set the limit of the number of MSPs allowed\n"
"       -Mp           same, as percentage of bases in cDNA\n"
"                     NOTE:  If used, both -Ma and -Mp must be specified!\n"
;




void
parseCommandLine(int argc, char **argv) {
  int arg = 1;

  while (arg < argc) {
    if        (strncmp(argv[arg], "--buildinfo", 3) == 0) {
        buildinfo_sim4db(stderr);
        buildinfo_libbri(stderr);
        buildinfo_libsim4(stderr);
        exit(1);
    } else if (strncmp(argv[arg], "-alignments", 4) == 0) {
      sim4params.setPrintAlignments(true);
    } else if (strncmp(argv[arg], "-alwaysprint", 4) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setAlwaysReport(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-cdna", 3) == 0) {
      arg++;
      cdnaFileName = argv[arg];
    } else if (strncmp(argv[arg], "-cut", 3) == 0) {
      arg++;
      double x = atof(argv[arg]);
      if (x < 0.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 0.0 (you gave %f)!\n", x);
        x = 0.0;
      }
      if (x > 1.0) {
        fprintf(stderr, "WARNING:  -cut adjusted to 1.0 (you gave %f)!\n", x);
        x = 1.0;
      }
      sim4params.setPolyTailPercent(x);
    } else if (strncmp(argv[arg], "-genomic", 2) == 0) {
      arg++;
      databaseFileName = argv[arg];
    } else if (strncmp(argv[arg], "-minc", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinCoverage(atoi(argv[arg]) / 100.0);
    } else if (strncmp(argv[arg], "-mini", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinPercentExonIdentity(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-minl", 5) == 0) {
      sim4params.setFindAllExons(true);
      arg++;
      sim4params.setMinCoverageLength(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-nod", 4) == 0) {
      sim4params.setIncludeDefLine(false);
    } else if (strncmp(argv[arg], "-non", 4) == 0) {
      sim4params.setDontForceCanonicalSplicing(true);
    } else if (strncmp(argv[arg], "-f", 2) == 0) {
      sim4params.setForceStrandPrediction(true);
    } else if (strncmp(argv[arg], "-o", 2) == 0) {
      arg++;
      outputFileName = argv[arg];
    } else if (strncmp(argv[arg], "-p", 2) == 0) {
      sim4params.setIgnorePolyTails(false);
    } else if (strncmp(argv[arg], "-sc", 3) == 0) {
      arg++;
      scriptFileName = argv[arg];
    } else if (strncmp(argv[arg], "-st", 3) == 0) {
      arg++;
      statsFileName = argv[arg];
    } else if (strncmp(argv[arg], "-to", 3) == 0) {
      arg++;
      touchFileName = argv[arg];
    } else if (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;
    } else if (strncmp(argv[arg], "-V", 2) == 0) {
      beExplicit = true;
    } else if (strncmp(argv[arg], "-YN", 3) == 0) {
      beYesNo = true;
    } else if (strncmp(argv[arg], "-H", 2) == 0) {
      arg++;
      sim4params.setRelinkWeight(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-K", 2) == 0) {
      arg++;
      sim4params.setMSPThreshold1(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-C", 2) == 0) {
      arg++;
      sim4params.setMSPThreshold2(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-Ma", 3) == 0) {
      arg++;
      sim4params.setMSPLimitAbsolute(atoi(argv[arg]));
    } else if (strncmp(argv[arg], "-Mp", 3) == 0) {
      arg++;
      sim4params.setMSPLimitPercent(atof(argv[arg]));
    } else if (strncmp(argv[arg], "-interspecies", 2) == 0) {
      sim4params.setInterspecies(true);
    } else {
      fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
    }

    arg++;
  }

#if 0
  //
  //  XXX:  I don't think this is ever needed, but not sure.
  //
  if ((_findAllExons) &&
      (_minCoverage == 0.0) &&
      (_minPercentExonIdentity == 0) &&
      (_minCoverageLength == 0)) {
    _minCoverage            = 0.90;
    _minCoverageLength      = 0;
    _minPercentExonIdentity = 90;
  }
#endif
}





void
sim4db(char          **scriptLines,
       u32bit          scriptLinesNum,
       FastAWrapper   *GENs,
#ifdef USE_CACHE
       FastACache     *ESTs,
#else
       FastAWrapper   *ESTs,
#endif
       int             fOutput) {
  double   startTime = getTime() - 1e-5;

  FastASequenceInCore  *ESTseq = 0L;
  u32bit                ESTiid = 0;
  FastASequenceInCore  *GENseq = 0L;
  u32bit                GENiid = 0;
  u32bit                GENlo  = 0;
  u32bit                GENhi  = 0;


  for (u32bit workDone=0; workDone<scriptLinesNum; workDone++) {
    if (beExplicit) {
      fprintf(stderr, "At %u (%6.3f per second) -- '%s'\n",
              workDone, workDone / (getTime() - startTime), scriptLines[workDone]);
      fflush(stderr);

#ifdef MEMORY_DEBUG
    _dump_allocated_delta(0);
#endif

    } else {
      if ((beVerbose) && ((workDone & 0xff) == 0xff)) {
        fprintf(stderr, " At %u (%6.3f per second)\r",
                workDone, workDone / (getTime() - startTime));
        fflush(stderr);
      }
    }


    //  Parse the command line to create a sim4command object
    //
    //  [-f|-r] -e <one-or-more-ESTid> -D GENid GENlo GENhi
    //
    //    -f  Forward only
    //    -r  Reverse only
    //    -s  strand indicator string
    //    -D  genSeqIID genLo genHi
    //    -d  dbfile (REMOVED)
    //    -e  estSeqIID
    //
    //  XXX:  est-list currently must be exactly one EST; giving
    //  multiple ESTs offers no speedup and is broken.
    //
    bool           doForward = true;
    bool           doReverse = true;

    {
      u32bit         argWords = 0;
      splitToWords   words(scriptLines[workDone]);

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
            GENiid = atoi(words.getWord(argWords + 1));
            GENlo  = atoi(words.getWord(argWords + 2));
            GENhi  = atoi(words.getWord(argWords + 3));
            argWords += 3;

            if ((GENseq == 0L) || (GENseq->getIID() != GENiid)) {
              delete GENseq;
              GENs->find(GENiid);
              GENseq = GENs->getSequence();
            }

            if ((GENlo == 0) && (GENhi == 0))
              GENhi = GENseq->sequenceLength();
            break;
          case 'e':
            ESTiid = atoi(words.getWord(argWords + 1));
            argWords++;

            if ((ESTseq == 0L) || (ESTseq->getIID() != ESTiid)) {
#ifdef USE_CACHE
              ESTseq = ESTs->getSequence(ESTiid);
#else
              delete ESTseq;
              ESTs->find(ESTiid);
              ESTseq = ESTs->getSequence();
#endif
            }
            break;
          default:
            //fprintf(stderr, "Unknown option '%s'\n", words.getWord(argWords));
            break;
        }

        argWords++;
      }
    }

    sim4command     *P4 = new sim4command(ESTseq,
                                          GENseq,
                                          GENlo,
                                          GENhi,
                                          doForward,
                                          doReverse);
    Sim4            *S4 = new Sim4(&sim4params);
    sim4polishList  *l4 = S4->run(P4);
    sim4polishList  &L4 = *l4;

    delete S4;
    delete P4;

    for (u32bit i=0; L4[i]; i++) {
      char *o = s4p_polishToString(L4[i]);

      errno = 0;
      write(fOutput, o, strlen(o) * sizeof(char));
      if (errno) {
        fprintf(stderr, "Couldn't write the output file '%s'.\n%s\n",
                outputFileName, strerror(errno));
        exit(1);
      }

      free(o);
    }

    if (beYesNo) {
      if (L4[0])
        fprintf(stdout, "%s -Y %d %d\n",
                scriptLines[workDone],
                L4[0]->percentIdentity,
                L4[0]->querySeqIdentity);
      else
        fprintf(stdout, "%s -N 0 0\n",
                scriptLines[workDone]);
    }

    delete l4;

#ifdef MEMORY_DEBUG
    if (beExplicit)
      _dump_allocated_delta(0);
#endif

  }
}




void
sim4dball(FastAWrapper  *GENs,
          FastAWrapper  *ESTs,
          int            fOutput) {
  double    startTime = getTime() - 1e-10;
  u32bit    workDone  = 0;

  FastASequenceInCore  *GEN = 0L;
  FastASequenceInCore  *EST = 0L;

  for (u32bit genIdx=0; genIdx < GENs->getNumberOfSequences(); genIdx++) {
    GENs->find(genIdx);
    GEN = GENs->getSequence();

    for (u32bit estIdx=0; estIdx < ESTs->getNumberOfSequences(); estIdx++) {
      ESTs->find(estIdx);
      EST = ESTs->getSequence();

      if (beExplicit) {
        fprintf(stderr, "At %u (%6.3f per second) -- '-e %u -D %u'\n",
                workDone, workDone / (getTime() - startTime), estIdx, genIdx);
        fflush(stderr);
      } else {
        if ((beVerbose) && ((workDone & 0x1f) == 0x1f)) {
          fprintf(stderr, " At %u (%6.3f per second)\r",
                  workDone, workDone / (getTime() - startTime));
          fflush(stderr);
        }
      }
      workDone++;

      sim4command     *P4 = new sim4command(EST,
                                            GEN,
                                            0,
                                            GEN->sequenceLength(),
                                            true,
                                            true);

      Sim4            *S4 = new Sim4(&sim4params);
      sim4polishList  *l4 = S4->run(P4);
      sim4polishList  &L4 = *l4;

      delete S4;
      delete P4;

      for (u32bit i=0; L4[i]; i++) {
        char *o = s4p_polishToString(L4[i]);

        errno = 0;
        write(fOutput, o, strlen(o) * sizeof(char));
        if (errno) {
          fprintf(stderr, "Couldn't write the output file '%s'.\n%s\n",
                  outputFileName, strerror(errno));
          exit(1);
        }

        free(o);
      }

      delete l4;
      delete EST;
    }
    delete GEN;
  }
}







int
main(int argc, char **argv) {
  u32bit       scriptLinesNum  = 0;
  char        *scriptLinesData = 0L;
  char       **scriptLines     = 0L;

  double       mainStartTime = getTime();

  parseCommandLine(argc, argv);

  if (cdnaFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No EST file?\n");
    exit(1);
  }

  if (databaseFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No database file?\n");
    exit(1);
  }

  if (outputFileName == 0L) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "--No output file?\n");
    exit(1);
  }


  //  Read (actually, just open) the sequences.
  //
#ifdef USE_CACHE
  FastAWrapper  *GENs = new FastAWrapper(databaseFileName);
  FastACache    *ESTs = new FastACache(cdnaFileName, 0, true);

  GENs->openIndex();
#else
  FastAWrapper  *GENs = new FastAWrapper(databaseFileName);
  FastAWrapper  *ESTs = new FastAWrapper(cdnaFileName, 2048);

  GENs->openIndex();
  ESTs->openIndex();
#endif


  //
  //  Read the script lines.
  //

  if (scriptFileName) {
    FILE           *F;
    struct stat     Fstat;

    errno = 0;
    if (stat(scriptFileName, &Fstat)) {
      fprintf(stderr, "sim4db: Couldn't stat the script file '%s'\n", scriptFileName);
      fprintf(stderr, "sim4db: %s\n", strerror(errno));
      exit(1);
    }

    size_t          scriptFileLen = Fstat.st_size;

    errno = 0;
    F = fopen(scriptFileName, "r");
    if (F == 0L) {
      fprintf(stderr, "sim4db: Couldn't open the script file '%s'\n", scriptFileName);
      fprintf(stderr, "sim4db: %s\n", strerror(errno));
      exit(1);
    }

    //  Allocate space for the script lines.
    //
    scriptLinesData = new char [scriptFileLen];

    //  Suck in the whole file.
    //
    errno = 0;
    if (fread(scriptLinesData, sizeof(char), scriptFileLen, F) == 0) {
#ifdef TRUE64BIT
      fprintf(stderr, "sim4db: Couldn't read %lu bytes from '%s'\nsim4db: %s\n",
              scriptFileLen, scriptFileName, strerror(errno));
#else
      fprintf(stderr, "sim4db: Couldn't read %d bytes from '%s'\nsim4db: %s\n",
              scriptFileLen, scriptFileName, strerror(errno));
#endif
      exit(1);
    }

    //  Done with the file.  Close it.
    //
    fclose(F);

    //  Count the number of lines in the file
    //
    for (size_t i=0; i<scriptFileLen; i++) {
      if (scriptLinesData[i] == '\n') {
        scriptLinesData[i] = 0;
        scriptLinesNum++;
      }
    }

    //  Allocate space for the script line pointers
    //
    scriptLines = new char * [scriptLinesNum + 1];

    //  Set the pointers to the start of each line
    //
    scriptLinesNum = 0;
    scriptLines[scriptLinesNum++] = scriptLinesData;

    for (size_t i=0; i<scriptFileLen; i++) {
      if (scriptLinesData[i] == 0)
        scriptLines[scriptLinesNum++] = scriptLinesData + i + 1;
    }

    //  We overcounted the lines by one.  Remove the last.
    //
    scriptLinesNum--;
  }



  //  Open the output file.  If the filename is "-", use stdout.
  //
  int fOutput;

  if (strcmp(outputFileName, "-") == 0) {
    fOutput = fileno(stdout);
  } else {
    errno = 0;
    fOutput = open(outputFileName,
                   O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC,
                   S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno) {
      fprintf(stderr, "Couldn't open the output file '%s'.\n%s\n",
              outputFileName, strerror(errno));
      exit(1);
    }
  }


#ifdef MEMORY_DEBUG
  if (beExplicit)
    _dump_allocated_delta(0);
#endif


  if (scriptFileName) {
    sim4db(scriptLines, scriptLinesNum, GENs, ESTs, fOutput);
  } else {
#ifndef USE_CACHE
    sim4dball(GENs, ESTs, fOutput);
#endif
  }

  //  Only close the file if it isn't stdout
  //
  if (strcmp(outputFileName, "-") != 0) {
    errno = 0;
    close(fOutput);
    if (errno) {
      fprintf(stderr, "Couldn't close the output file '%s'.\n%s\n",
              outputFileName, strerror(errno));
      exit(1);
    }
  }

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

  delete ESTs;
  delete GENs;

  delete [] scriptLinesData;
  delete [] scriptLines;

  return(0);
}
