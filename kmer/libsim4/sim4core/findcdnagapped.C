//                    Confidential -- Do Not Distribute
//   Copyright (c) 2002 Applera Corporation through the Celera Genomics Group
//                           All Rights Reserved.
//
//  brian.walenz@celera.com, Fri Sep 27 15:47:27 EDT 2002

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include "libbri.H"
#include "sim4reader.h"
#include "sim4.H"
#include "sim4db.H"


//  -a -i /dev5/walenz/dbest-20020331/1-polishes-good-with-cdna-gaps -g /dev5/walenz/GENOMES/hum_scaffold_axis.fasta -c /dev5/walenz/ESTs/dbEST_human_31mar02.fasta


#define  MIN_EXON_LENGTH       20
#define  MIN_PERCENT_IDENTITY  90

#define  LOOP_SCALE       1.0
#define  LOOP_CONSTANT   10


sim4dbParameters   dbParams;


//  Filter by cDNA gaps.
//
//  -i inputpolishes -g gappedpolishes -u ungappedpolishes
//

FILE*
openForRead(char *name) {
  errno = 0;
  int file = open(name, O_RDONLY | O_LARGEFILE);
  if (errno) {
    fprintf(stderr, "%s: Couldn't open() '%s' for reading.\n%s\n", name, strerror(errno));
    exit(1);
  }

  errno = 0;
  FILE *F = fdopen(file, "r");
  if (errno) {
    fprintf(stderr, "%s: Couldn't fdopen() '%s'.\n%s\n", name, strerror(errno));
    exit(1);
  }

  return(F);
}


FILE*
openForWrite(char *name) {
  errno = 0;
  int file = open(name,
                  O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE,
                  S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  if (errno) {
    fprintf(stderr, "%s: Couldn't open() '%s' for writing.\n%s\n", name, strerror(errno));
    exit(1);
  }

  errno = 0;
  FILE *F = fdopen(file, "w");
  if (errno) {
    fprintf(stderr, "%s: Couldn't fdopen() '%s'.\n%s\n", name, strerror(errno));
    exit(1);
  }

  return(F);
}



bool
lowComplexityExon(char *s, int minexonlength, double qualitythreshold) {
  int    cnt[5][5] = { {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};
  int    map[256]  = {0};
  int    i, j, l = 0;
  int    a=0, b=0, c=0;
  double qual = 0.0;

  if (s == 0L)
    return(false);

  map['A'] = map['a'] = 1;
  map['C'] = map['c'] = 2;
  map['G'] = map['g'] = 3;
  map['T'] = map['t'] = 4;

  for (i=0, j=1, l=0; s[j]; i++, j++, l++)
    cnt[map[s[i]]][map[s[j]]]++;

  if (l > minexonlength)
    return(false);

  for (i=0; i<5; i++) {
    for (j=0; j<5; j++) {
      if (a < cnt[i][j]) {
        c = b;
        b = a;
        a = cnt[i][j];
      } else if (b < cnt[i][j]) {
        c = b;
        b = cnt[i][j];
      } else if (c < cnt[i][j]) {
        c = cnt[i][j];
      }
    }
  }

  qual = (double)(a+b+c) / (double)(l);

  return(qual > qualitythreshold);
}


bool
junkExon(sim4polishExon *e, u32bit minExonLength, u32bit minPercentIdentity) {

  if ((e->estTo - e->estFrom + 1 >= minExonLength) &&
      (e->percentIdentity >= minPercentIdentity) &&
      (lowComplexityExon(e->estAlignment, minExonLength, 0.75) == false))
    return(false);
  else
    return(true);
}


bool
trimJunk(sim4polish *p) {

  //fprintf(stdout, "TRIMJUNK------------------------------------------------------\n");

  //  Find the first and last cDNA gaps.
  //
  u32bit   firstGapIdx = ~u32bitZERO;
  u32bit   lastGapIdx  = ~u32bitZERO;

  for (u32bit i=1; i<p->numExons; i++) {
    if (p->exons[i-1].estTo+1 != p->exons[i].estFrom) {
      if (i < firstGapIdx)
        firstGapIdx = i;

      lastGapIdx = i;
    }
  }

  //  If we don't find a gap, return false
  //
  if (lastGapIdx == ~u32bitZERO)
    return(false);

  //  See if we can trim off junk exons before the first gap
  //
  bool  doTrimFirst = true;
  bool  doTrimLast  = true;

  for (u32bit i=0; i<firstGapIdx; i++) {
    if (junkExon(p->exons+i, MIN_EXON_LENGTH, MIN_PERCENT_IDENTITY) == false)
      doTrimFirst = false;
  }

  for (u32bit i=lastGapIdx; i<p->numExons; i++) {
    if (junkExon(p->exons+i, MIN_EXON_LENGTH, MIN_PERCENT_IDENTITY) == false)
      doTrimLast = false;
  }

  if (doTrimLast) {
    for (u32bit i=p->numExons-1; i>=lastGapIdx; i--) {
      s4p_deleteExon(p, i);
    }
  }

  if (doTrimFirst) {
    for (u32bit i=0; i<firstGapIdx; i++) {
      s4p_deleteExon(p, 0);
    }
  }

  return(doTrimFirst | doTrimLast);
}


bool
gapIsALoop(sim4polish *p, u32bit gap) {
  u32bit c = p->exons[gap].estFrom - p->exons[gap-1].estTo;
  u32bit g = p->exons[gap].genFrom - p->exons[gap-1].genTo;

  return((LOOP_SCALE * (double)c + LOOP_CONSTANT) > (double)g);
}


bool
gapIsLowQuality(SequenceManager *C, sim4polish *p, u32bit gap) {
  unsigned char  *seq = C->getSequence(p->estID);
  int             cnt = 0;
  double          lim = 0.1 * (double)(p->exons[gap].estFrom - p->exons[gap-1].estTo - 1);

  for (int j=p->exons[gap-1].estTo+1; j<p->exons[gap].estFrom; j++)
    if (seq[j] == 'N')
      cnt++;

  return((double)cnt >= lim);
}










//
//  Run sim4 on two sequences, using only a portion of each.
//
sim4polish *
sim4db(sim4polish *p,
       unsigned char *est, int estlow, int esthigh,
       unsigned char *gen, int genlow, int genhigh) {

  //  We trick the sim4parameters class into running by creating a
  //  new sequencemanger with the EST / genomic region we want to
  //  run.
  //
  //  SequenceManager checks that the limits make sense.
  //
  SequenceManager *E = new SequenceManager(est, estlow, esthigh);
  SequenceManager *G = new SequenceManager(gen, genlow, genhigh);

  //  Not the greatest method, but we need to make a command string
  //
  char cmd[32];

  if (p->matchOrientation == MATCH_FORWARD) {
    sprintf(cmd, "-f -e 0 -D 0 0 0\n");
  } else {
    sprintf(cmd, "-r -e 0 -D 0 0 0\n");
  }

  //sprintf(cmd, "-e 0 -D 0 0 0\n");

  sim4parameters  *P4 = new sim4parameters(cmd, E, G);
  Sim4            *S4 = new Sim4(dbParams._mspThresh1, dbParams._mspThresh2);
  char            *O4 = S4->run(P4);

  delete S4;
  delete P4;

  sim4polish *x = s4p_stringToPolish(O4);

  //  Massage the polish the make it correct
  //
  if (x) {
    x->estID = p->estID;
    x->genID = p->genID;
    x->genLo = genlow;

    if (x->estDefLine)
      free(x->estDefLine);

    if (x->genDefLine)
      free(x->genDefLine);

    x->estDefLine = strdup(p->estDefLine);
    x->genDefLine = strdup(p->genDefLine);
  }


  delete [] O4;

  delete E;
  delete G;

  return(x);
}









int
main(int argc, char **argv) {
  bool        analyze    = false;
  char       *inputName  = 0L;
  char       *genomeName = 0L;
  char       *cdnaName   = 0L;
  int         processed  = 0;

  if (argc == 1) {
    fprintf(stderr, "usage: %s -i inputpolishes -g gappedpolishes -u ungappedpolishes\n", argv[0]);
    exit(1);
  }

  int arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-analyze", 2) == 0) {
      analyze = true;
    } else if (strncmp(argv[arg], "-input", 2) == 0) {
      inputName = argv[++arg];
    } else if (strncmp(argv[arg], "-genome", 2) == 0) {
      genomeName = argv[++arg];
    } else if (strncmp(argv[arg], "-cdna", 2) == 0) {
      cdnaName = argv[++arg];
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
    }

    arg++;
  }

  if (inputName == 0L) {
    fprintf(stderr, "Hmmm.  But what -input shall I use?\n");
    exit(1);
  }

  if (genomeName == 0L) {
    fprintf(stderr, "Hmmm.  But what -genome shall I use?\n");
    exit(1);
  }

  if (cdnaName == 0L) {
    fprintf(stderr, "Hmmm.  But what -cdna shall I use?\n");
    exit(1);
  }

  FILE             *I = openForRead(inputName);
  SequenceManager  *G = new SequenceManager(genomeName, 1);
  SequenceManager  *C = new SequenceManager(cdnaName,   256);


  //  Create some output files
  //
  char   *outputName           = new char [strlen(inputName) + 256];
  char   *analysisName         = new char [strlen(inputName) + 256];
#ifdef WRITE_MANY_FILES
  char   *resUhohName          = new char [strlen(inputName) + 256];
  char   *resOKName            = new char [strlen(inputName) + 256];
  char   *resGoodName          = new char [strlen(inputName) + 256];
  char   *resStrangeName       = new char [strlen(inputName) + 256];
  char   *resBadName           = new char [strlen(inputName) + 256];
  char   *resCheckName         = new char [strlen(inputName) + 256];
#endif

  sprintf(outputName,           "%s.fixedGaps",    inputName);
  sprintf(analysisName,         "%s.log.analysis", inputName);
#ifdef WRITE_MANY_FILES
  sprintf(resUhohName,          "%s.res.uhoh",    inputName);
  sprintf(resOKName,            "%s.res.ok",      inputName);
  sprintf(resGoodName,          "%s.res.good",    inputName);
  sprintf(resStrangeName,       "%s.res.strange", inputName);
  sprintf(resBadName,           "%s.res.bad",     inputName);
  sprintf(resCheckName,         "%s.res.check",   inputName);
#endif

  FILE   *outputFile           = openForWrite(outputName);
  FILE   *analysisFile         = 0L;
#ifdef WRITE_MANY_FILES
  FILE   *resUhohFile          = openForWrite(resUhohName);
  FILE   *resOKFile            = openForWrite(resOKName);
  FILE   *resGoodFile          = openForWrite(resGoodName);
  FILE   *resStrangeFile       = openForWrite(resStrangeName);
  FILE   *resBadFile           = openForWrite(resBadName);
  FILE   *resCheckFile         = openForWrite(resCheckName);
#endif

  if (analyze)
    analysisFile               = openForWrite(analysisName);


  //  Configure sim4 to generate alignments, but not deflines
  //
  dbParams._printAlignments = true;
  dbParams._includeDefLine  = false;


  //  Stats
  //
  int    gapMapForward  = 0;
  int    gapMapReverse  = 0;

  int    resultTrim     = 0;
  int    resultLoop     = 0;
  int    resultLowQ     = 0;
  int    resultLoopLowQ = 0;

  int    resultUhoh     = 0;
  int    resultOK       = 0;
  int    resultGood     = 0;
  int    resultStrange  = 0;
  int    resultBad      = 0;
  int    resultCheck    = 0;




  while (!feof(I)) {
    sim4polish *p = readPolish(I);

    if (p) {
      sim4polish *o = copyPolish(p);

      //fprintf(stdout, "--------------------------------------------------------------\n");
      //printPolish(stdout, p);


      ///////////////////////////////////////
      //
      //  Trim junk "recursively".  Doesn't affect many polishes, but
      //  the few it does make it worth the effort.
      //
      bool    doTrimming = true;
      u32bit  numTrims   = 0;

      while (doTrimming) {
        numTrims++;
        doTrimming = trimJunk(p);
      }

      bool   hasagap         = false;
      bool   moreprocessing  = false;
      bool   looped          = false;
      bool   lowquality      = false;


      //  For each cDNA gap, we test for loops and low quality cDNA
      //  sequence.  If any cDNA gap is explained by neither, we must
      //  continue processing.  If all are explained by at least one,
      //  we can stop.
      //
      for (u32bit i=1; i<p->numExons; i++) {
        if (p->exons[i-1].estTo+1 != p->exons[i].estFrom) {
          hasagap = true;

          bool lp = gapIsALoop(p, i);
          bool lq = gapIsLowQuality(C, p, i);

          looped     |= lp;
          lowquality |= lq;

          //  if the gap is not looped, and not low quality, then
          //  we need to explore the entire match more
          //
          if (!lp && !lq)
            moreprocessing = true;
        }
      }



      if (hasagap == false) {
        resultTrim++;

        //
        //  Trimming junk removed the gap
        //
        //  hasOnlyLowQExonsFile

        if (analyze) {
          fprintf(analysisFile, "---ONLY LOWQ EXONS BEGIN\n");
          printPolish(analysisFile, o);
          printPolish(analysisFile, p);
          fprintf(analysisFile, "---ONLY LOWQ EXONS END\n");
        }

        printPolish(outputFile, p);
      } else {

        //
        //  Still have a gap.
        //

        if (moreprocessing == false) {

          //
          //  All gaps explained.
          //  Just print it.
          //

          printPolish(outputFile, p);

          if (looped && lowquality) {
            resultLoopLowQ++;
            if (analyze) {
              fprintf(analysisFile, "---BOTH LOOPED AND LOWQUALITY BEGIN\n");
              printPolish(analysisFile, o);
              printPolish(analysisFile, p);
              fprintf(analysisFile, "---BOTH LOOPED AND LOWQUALITY END\n");
            }
          } else if (looped) {
            resultLoop++;
            if (analyze) {
              fprintf(analysisFile, "---ONLY LOOPED GAP BEGIN\n");
              printPolish(analysisFile, o);
              printPolish(analysisFile, p);
              fprintf(analysisFile, "---ONLY LOOPED GAP END\n");
            }
          } else if (lowquality) {
            resultLowQ++;
            if (analyze) {
              fprintf(analysisFile, "---ONLY LOWQUALITY GAP BEGIN\n");
              printPolish(analysisFile, o);
              printPolish(analysisFile, p);
              fprintf(analysisFile, "---ONLY LOWQUALITY GAP END\n");
            }
          }
        } else {

          //
          //  At least one gap has no explanation.  Examine all the gaps
          //  in detail.
          //
          //  analysisFile has a log of our analysis:
          //
          //  1) sim4 of the cDNA gap vs the corresponding genomic
          //  gap, with C=12, K=15, W=10, filtering at 80% coverage
          //  and 90% identity
          //
          //  2) sim4 of the cDNA vs the 50kb extended genomic region,
          //  with H=2000

          bool foundAnExon = false;


          //  #1
          //
          u32bit  numGaps = 0;

          fprintf(analysisFile, "---ANALYSIS BEGIN\n");
          printPolish(analysisFile, o);

          for (u32bit i=1; i<p->numExons; i++) {
            if (p->exons[i-1].estTo+1 != p->exons[i].estFrom) {
              numGaps++;
              fprintf(analysisFile, "gap polish %d[%d-%d] vs %d[%d-%d]\n",
                      p->estID, p->exons[i-1].estTo+1, p->exons[i].estFrom,
                      p->genID, p->genLo + p->exons[i-1].genTo+1, p->genLo + p->exons[i].genFrom);

              //  Can't set 'W'

              dbParams._mspThresh1   = 15;  //  K
              dbParams._mspThresh2   = 12;  //  C
              dbParams._relinkWeight = 500; //  H

              //  Beware of the 1-base to 0-space conversion!
              sim4polish *x;
              if (p->matchOrientation == MATCH_FORWARD)
                x = sim4db(p,
                           C->getSequence(p->estID), p->exons[i-1].estTo, p->exons[i].estFrom,
                           G->getSequence(p->genID), p->genLo + p->exons[i-1].genTo, p->genLo + p->exons[i].genFrom);
              else
                x = sim4db(p,
                           C->getSequence(p->estID), p->estLen - p->exons[i-1].estTo + 1, p->estLen - p->exons[i].estFrom + 1,
                           G->getSequence(p->genID), p->genLo + p->exons[i-1].genTo, p->genLo + p->exons[i].genFrom);

              if (x) {
                if ((x->percentIdentity >= 90) && (x->querySeqIdentity >= 70)) {
                  foundAnExon = true;

                  if (p->matchOrientation == MATCH_FORWARD)
                    gapMapForward++;
                  else
                    gapMapReverse++;

                  fprintf(analysisFile, "Found something at intron %d!\n", i);
                  printPolish(analysisFile, x);
                } else {
                  fprintf(analysisFile, "Found something USELESS at intron %d!\n", i);
                  printPolish(analysisFile, x);
                }

                destroyPolish(x);
              }
            }
          }



          //  #2
          //
          fprintf(analysisFile, "repolish   %d[%d-%d] vs %d[%d-%d]\n",
                  p->estID, 0, p->estLen,
                  p->genID, p->genLo + p->exons[0].genFrom - 50000, p->genLo + p->exons[p->numExons-1].genTo + 50000);

          dbParams._mspThresh1   = 0;    //  K (use default)
          dbParams._mspThresh2   = 0;    //  C (use default)
          dbParams._relinkWeight = 2000; //  H

          sim4polish *x = sim4db(p,
                                 C->getSequence(p->estID), 0, p->estLen,
                                 G->getSequence(p->genID), p->genLo + p->exons[0].genFrom - 50000, p->genLo + p->exons[p->numExons-1].genTo + 50000);

          if (x) {
            fprintf(analysisFile, "repolished from %d coverage and %d identity to %d coverage and %d identity\n",
                    p->querySeqIdentity, p->percentIdentity,
                    x->querySeqIdentity, x->percentIdentity);

            printPolish(analysisFile, x);

            //  XXX:  Which to put into the output file??

            printPolish(outputFile, x);
          } else {
            fprintf(stderr, "WARNING: didn't find anything by repolishing?\n");
            fprintf(analysisFile, "WARNING: didn't find anything by repolishing?\n");
          }


          //
          //  finally, we compare #1 with #2 to make sure that both
          //  methods find (or both don't find) any missing exons
          //


          int totalGapSize = 0;

          for (u32bit i=1; i<p->numExons; i++)
            totalGapSize += p->exons[i].estFrom - p->exons[i-1].estTo - 1;

          bool   gotBetter     = false;
          bool   stayedTheSame = false;
          bool   gotWorse      = false;

#if 0
          gotBetter     = (p->numMatches  < x->numMatches);
          stayedTheSame = (p->numMatches == x->numMatches);
          gotWorse      = !gotBetter && !stayedTheSame;
#else
          gotBetter     = (p->numMatches  < x->numMatches - totalGapSize/2) && (p->percentIdentity-2 <= x->percentIdentity);
          stayedTheSame = (p->numMatches == x->numMatches)                  && (p->percentIdentity-2 <= x->percentIdentity);
          gotWorse      = !gotBetter && !stayedTheSame;
#endif

          if (gotBetter) {
            if (foundAnExon) {
              fprintf(analysisFile, "RESULT=GOOD: it got better AND we found an exon separately\n");
              resultGood++;
            } else {
              fprintf(analysisFile, "RESULT=STRANGE: it got better without finding an exon separately\n");
              resultStrange++;
            }
          }

          if (stayedTheSame) {
            if (foundAnExon) {
              fprintf(analysisFile, "RESULT=UHOH: it stayed the same AND we found an exon separately\n");
              resultUhoh++;
            } else {
              fprintf(analysisFile, "RESULT=OK: no change\n");
              resultOK++;
            }
          }

          if (gotWorse) {
            if (foundAnExon) {
              fprintf(analysisFile, "RESULT=BAD: it got worse AND we found an exon separately\n");
              resultBad++;
            } else {
              fprintf(analysisFile, "RESULT=CHECK: it got worse without finding an exon separately\n");
              resultCheck++;
            }
          }

          fprintf(analysisFile, "---ANALYSIS END\n");
        }
      }

      destroyPolish(p);
    }

    processed++;

    fprintf(stderr, "%8d] | gap: %5d/%5d | jk:%6d lp:%6d lq:%6d pq:%6d | >+:%5d >-:%5d | =+:%5d =-:%5d | <+:%5d <-:%5d\r",
            processed,
            gapMapForward, gapMapReverse,
            resultTrim, resultLoop, resultLowQ, resultLoopLowQ,
            resultGood, resultStrange, resultUhoh, resultOK, resultBad, resultCheck);
    fflush(stderr);
  }

  fprintf(stderr, "%8d] | gap: %5d/%5d | jk:%6d lp:%6d lq:%6d pq:%6d | >+:%5d >-:%5d | =+:%5d =-:%5d | <+:%5d <-:%5d\n",
          processed,
          gapMapForward, gapMapReverse,
          resultTrim, resultLoop, resultLowQ, resultLoopLowQ,
          resultGood, resultStrange, resultUhoh, resultOK, resultBad, resultCheck);
}
