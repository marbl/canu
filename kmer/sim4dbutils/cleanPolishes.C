#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "sim4.H"

char const *usage =
"usage: %s [-threshold t] [-savejunk] [-quiet] [-debug]\n"
"  -threshold    Introns bigger than this are candidates for trimming (default = 100000).\n"
"  -quiet        Don't print unmodified matches\n"
"  -beforeafter  Save (in separate files) the before/after of each modified match\n"
"  -segregate    Save (in separate files) the after of each modified match\n"
"  -savejunk     Also print the trimmed pieces (as separate matches)\n"
"\n";

//#define  MIN_EXON_LENGTH       50
//#define  MIN_PERCENT_IDENTITY  88

#define  MIN_EXON_LENGTH       20
#define  MIN_PERCENT_IDENTITY  90


FILE *openDebugFile(char const *name) {
  FILE *f;

  errno=0;
  f = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "Can't debug file '%s' for writing!\n%s\n", name, strerror(errno));
    exit(1);
  }

  return(f);
}


bool
lowComplexityExon(char *s) {
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

  if (l > MIN_EXON_LENGTH)
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

  return(qual > 0.75);
}


//  Delete exons before/after a specific intron.
//
void
trimExonsBefore(int intronSplit, sim4polish *p) {
  for (int i=0; i<intronSplit; i++)
    p->s4p_deleteExon(0);
}

void
trimExonsAfter(int intronSplit, sim4polish *p) {
  for (int i=p->_numExons-1; i>=intronSplit; i--)
    p->s4p_deleteExon(i);
}



int
main(int argc, char ** argv) {
  int  totMatches = 0;
  int  oneExon    = 0;
  int  smaIntron  = 0;
  int  junkFirst  = 0;
  int  junkLast   = 0;
  int  junkBoth   = 0;
  int  splitOnGap = 0;
  int  goodQual   = 0;
  int  probGood   = 0;

  bool    filter        = true;
  bool    saveJunk      = false;
  u32bit  intronLimit   = 100000;

  //  Before / after files
  //
  bool  beforeafter   = false;
#if 0
  FILE *splGood       = 0L;
  FILE *splProbGood   = 0L;
#endif
  FILE *splJunkLeft   = 0L;
  FILE *splJunkRight  = 0L;
  FILE *splJunkBoth   = 0L;
  FILE *splIntronGap  = 0L;

  //  Segregation files
  //
  bool  segregate     = false;
#if 0
  FILE *filtOne       = 0L;
  FILE *filtAllSmall  = 0L;
#endif
  FILE *filtGood      = 0L;
  FILE *filtProbGood  = 0L;
  FILE *filtJunkLeft  = 0L;
  FILE *filtJunkRight = 0L;
  FILE *filtJunkBoth  = 0L;
  FILE *filtIntronGap = 0L;

  bool  hasBeenWarned = false;

  bool  beVerbose     = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-threshold", 2) == 0) {
      intronLimit = atoi(argv[++arg]);
    } else if (strncmp(argv[arg], "-quiet", 2) == 0) {
      fprintf(stderr, "QUIET MODE ENABLED -- non-modified matches not output!\n");
      filter = false;
    } else if (strncmp(argv[arg], "-beforeafter", 2) == 0) {
      fprintf(stderr, "DEBUG MODE ENABLED -- many 'spl.*' files created!\n");
      beforeafter  = true;
#if 0
      splGood      = openDebugFile("spl.good");
      splProbGood  = openDebugFile("spl.probGood");
#endif
      splJunkLeft  = openDebugFile("spl.junkLeft");
      splJunkRight = openDebugFile("spl.junkRight");
      splJunkBoth  = openDebugFile("spl.junkBoth");
      splIntronGap = openDebugFile("spl.intronGap");
    } else if (strncmp(argv[arg], "-segregate", 3) == 0) {
      fprintf(stderr, "SEGREGATION MODE ENABLED -- many 'filt.*' files created!\n");
      segregate     = true;
#if 0
      filtOne       = openDebugFile("filt.filtOne");
      filtAllSmall  = openDebugFile("filt.allSmall");
#endif
      filtGood      = openDebugFile("filt.good");
      filtProbGood  = openDebugFile("filt.probGood");
      filtJunkLeft  = openDebugFile("filt.junkLeft");
      filtJunkRight = openDebugFile("filt.junkRight");
      filtJunkBoth  = openDebugFile("filt.junkBoth");
      filtIntronGap = openDebugFile("filt.intronGap");
    } else if (strncmp(argv[arg], "-savejunk", 3) == 0) {
      saveJunk = true;
    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;
    }

    arg++;
  }

  if (isatty(fileno(stdin))) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "error: I cannot read polishes from the terminal!\n");
    exit(1);
  }

  if (isatty(fileno(stdout)) && filter) {
    fprintf(stderr, usage, argv[0]);
    fprintf(stderr, "error: Please redirect the polishes (stdout) to a file.\n");
    exit(1);
  }

  if (beVerbose)
    fprintf(stderr, "A big intron is one that is at least "u32bitFMT"bp long.\n", intronLimit);

  sim4polish *p = new sim4polish(stdin);
  while (p->_numExons > 0) {
    u32bit exA;
    u32bit exB;

    if (p->_numExons == 1) {
      oneExon++;
      if (filter)
        p->s4p_printPolish(stdout);
#if 0
      if (segregate)
        p->s4p_printPolish(filtOneExon);
#endif
    } else {

      //  Find the big intron.  We assume there is only one big intron.
      //
      u32bit biggestIntron = 0;
      u32bit intronSplit   = 0;
      u32bit intronOri     = 0;

      for (exA=0, exB=1; exB < p->_numExons; exA++, exB++) {
        u32bit dist = p->_exons[exB]._genFrom - p->_exons[exA]._genTo + 1;
        if (dist > biggestIntron) {
          biggestIntron = dist;
          intronSplit   = exB;
          intronOri     = p->_exons[exA]._intronOrientation;
        }
      }

      if (intronOri == 0) {
        fprintf(stderr, "didn't find the largest intron? (got zero)?\n");
        exit(1);
      }

      if (intronOri == SIM4_INTRON_NONE) {
        fprintf(stderr, "biggest intron isn't an intron? (got none)?\n");
        exit(1);
      }

      if (biggestIntron < intronLimit) {
        smaIntron++;
        if (filter)
          p->s4p_printPolish(stdout);
#if 0
        if (segregate)
          p->s4p_printPolish(filtAllSmall);
#endif
      } else {

        //  Declare the split obvious if all exons on either side are
        //  below MIN_EXON_LENGTH, difficult otherwise.
        //
        bool  killFirst = true;
        bool  killLast  = true;

        for (u32bit i=0; i<intronSplit; i++)
          if ((p->_exons[i]._estTo - p->_exons[i]._estFrom + 1 >= MIN_EXON_LENGTH) &&
              (p->_exons[i]._percentIdentity >= MIN_PERCENT_IDENTITY) &&
              (lowComplexityExon(p->_exons[i]._estAlignment) == false))
            killFirst = false;

        for (u32bit i=intronSplit; i<p->_numExons; i++)
          if ((p->_exons[i]._estTo - p->_exons[i]._estFrom + 1 >= MIN_EXON_LENGTH) &&
              (p->_exons[i]._percentIdentity >= MIN_PERCENT_IDENTITY) &&
              (lowComplexityExon(p->_exons[i]._estAlignment) == false))
            killLast = false;


        //  Sometimes, all exons look crappy.  If they have a large
        //  intron too, just kill the match.
        //
        if ((killFirst == true) && (killLast == true)) {
          junkBoth++;

          if ((hasBeenWarned == false) &&
              ((p->_exons[0]._estAlignment == 0L) || (p->_exons[0]._genAlignment == 0L))) {
            hasBeenWarned = true;
            fprintf(stderr, "cleanPolishes: Need alignments to recompute scores correctly!\n");
          }

          sim4polish *a = new sim4polish(p);
          sim4polish *b = new sim4polish(p);
          trimExonsAfter(intronSplit, a);
          trimExonsBefore(intronSplit, b);

          if (filter && saveJunk) {
            a->s4p_printPolish(stdout);
            b->s4p_printPolish(stdout);
          }

          if (beforeafter) {
            fprintf(splJunkBoth, "====================\n");
            p->s4p_printPolish(splJunkBoth);
            a->s4p_printPolish(splJunkBoth);
            b->s4p_printPolish(splJunkBoth);
          }

          if (segregate) {
            a->s4p_printPolish(filtJunkBoth);
            b->s4p_printPolish(filtJunkBoth);
          }

          delete a;
          delete b;
        }

        //  If the first half (before the big intron) is crappy, delete
        //  those exons.
        //
        if ((killFirst == true) && (killLast == false)) {
          junkFirst++;

          sim4polish *a = new sim4polish(p);
          sim4polish *b = new sim4polish(p);
          trimExonsAfter(intronSplit, a);
          trimExonsBefore(intronSplit, b);

          if (filter) {
            if (saveJunk)
              a->s4p_printPolish(stdout);
            b->s4p_printPolish(stdout);
          }

          if (beforeafter) {
            fprintf(splJunkLeft, "====================\n");
            p->s4p_printPolish(splJunkLeft);
            a->s4p_printPolish(splJunkLeft);
            b->s4p_printPolish(splJunkLeft);
          }

          if (segregate) {
            a->s4p_printPolish(filtJunkLeft);
            b->s4p_printPolish(filtJunkLeft);
          }

          delete a;
          delete b;
        }

        if ((killFirst == false) && (killLast == true)) {
          junkLast++;

          sim4polish *a = new sim4polish(p);
          sim4polish *b = new sim4polish(p);
          trimExonsAfter(intronSplit, a);
          trimExonsBefore(intronSplit, b);

          if (filter) {
            a->s4p_printPolish(stdout);
            if (saveJunk)
              b->s4p_printPolish(stdout);
          }

          if (beforeafter) {
            fprintf(splJunkRight, "====================\n");
            p->s4p_printPolish(splJunkRight);
            a->s4p_printPolish(splJunkRight);
            b->s4p_printPolish(splJunkRight);
          }

          if (segregate) {
            a->s4p_printPolish(filtJunkRight);
            b->s4p_printPolish(filtJunkRight);
          }

          delete a;
          delete b;
        }

        if ((killFirst == false) && (killLast == false)) {
          if (intronOri == SIM4_INTRON_GAP) {
            splitOnGap++;

            //  Break the polish into two pieces, one before and one
            //  after the large intron.  This is done by copying the
            //  entire polish, then deleting one half from each.
            //
            //  XXX If we want to update the strand prediction of the
            //  split pieces, we should
            //
            //  a) make sure that all the intron signals agree
            //  b) make sure that the percent identites of each exon are > 90%
            //
            //  For now, we don't.

            sim4polish *a = new sim4polish(p);
            sim4polish *b = new sim4polish(p);
            trimExonsBefore(intronSplit, a);
            trimExonsAfter(intronSplit, b);

            if (filter) {
              a->s4p_printPolish(stdout);
              b->s4p_printPolish(stdout);
            }

            if (beforeafter) {
              fprintf(splIntronGap, "====================\n");
              p->s4p_printPolish(splIntronGap);
              a->s4p_printPolish(splIntronGap);
              b->s4p_printPolish(splIntronGap);
            }

            if (segregate) {
              a->s4p_printPolish(filtIntronGap);
              b->s4p_printPolish(filtIntronGap);
            }

            delete a;
            delete b;
          } else {

            //  If there is a valid strand prediction and
            //      a) all exons >= 90%
            //      b) all exons >= 95%
            //      c) all exons >= 95%, except first and last, which can be >= 90%
            //  save the match as is.
            //

            bool  qualIsC = ((p->_exons[0]._percentIdentity >= 90) &&
                             (p->_exons[p->_numExons-1]._percentIdentity >= 90));

            for (exA=1; exA < p->_numExons-1; exA++)
              if (p->_exons[exA]._percentIdentity < 95)
                qualIsC = false;

            //  If the match looks good, but just has a large intron, keep it.
            //
            if (qualIsC &&
                ((p->_strandOrientation == SIM4_STRAND_POSITIVE) ||
                 (p->_strandOrientation == SIM4_STRAND_NEGATIVE))) {
              goodQual++;
              if (filter)
                p->s4p_printPolish(stdout);
              if (segregate)
                p->s4p_printPolish(filtGood);
            } else {
              probGood++;
              if (filter)
                p->s4p_printPolish(stdout);
              if (segregate)
                p->s4p_printPolish(filtProbGood);
            }
          }
        }

      }  //  Has a big intron
    }  //  More than one exon

    totMatches++;

    delete p;
    p = new sim4polish(stdin);
  }

  if (beforeafter) {
#if 0
    fclose(splGood);
    fclose(splProbGood);
#endif
    fclose(splJunkLeft);
    fclose(splJunkRight);
    fclose(splJunkBoth);
    fclose(splIntronGap);
  }

  if (segregate) {
#if 0
    fclose(filtOne);
    fclose(filtAllSmall);
#endif
    fclose(filtGood);
    fclose(filtProbGood);
    fclose(filtJunkLeft);
    fclose(filtJunkRight);
    fclose(filtJunkBoth);
    fclose(filtIntronGap);
  }

  if (beVerbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "oneExon:         %7d\n", oneExon);
    fprintf(stderr, "allSmallIntrons: %7d\n", smaIntron);
    fprintf(stderr, "good:            %7d\n", goodQual);
    fprintf(stderr, "probably good:   %7d\n", probGood);
    fprintf(stderr, "junkExonsLeft:   %7d\n", junkFirst);
    fprintf(stderr, "junkExonsRight:  %7d\n", junkLast);
    fprintf(stderr, "junkExonsBoth:   %7d\n", junkBoth);
    fprintf(stderr, "intronOnGap:     %7d\n", splitOnGap);
    fprintf(stderr, "total:           %7d\n", totMatches);
  }

  return(0);
}
