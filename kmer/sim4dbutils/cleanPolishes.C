#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "sim4.H"

//#define  MIN_EXON_LENGTH       50
//#define  MIN_PERCENT_IDENTITY  88

#define  MIN_EXON_LENGTH       20
#define  MIN_PERCENT_IDENTITY  90


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
  sim4polishWriter *splGood       = 0L;
  sim4polishWriter *splProbGood   = 0L;
#endif
  sim4polishWriter *splJunkLeft   = 0L;
  sim4polishWriter *splJunkRight  = 0L;
  sim4polishWriter *splJunkBoth   = 0L;
  sim4polishWriter *splIntronGap  = 0L;

  //  Segregation files
  //
  bool  segregate     = false;
#if 0
  sim4polishWriter *filtOne       = 0L;
  sim4polishWriter *filtAllSmall  = 0L;
#endif
  sim4polishWriter *filtGood      = 0L;
  sim4polishWriter *filtProbGood  = 0L;
  sim4polishWriter *filtJunkLeft  = 0L;
  sim4polishWriter *filtJunkRight = 0L;
  sim4polishWriter *filtJunkBoth  = 0L;
  sim4polishWriter *filtIntronGap = 0L;

  sim4polishStyle style = sim4polishStyleDefault;

  bool  hasBeenWarned = false;

  bool  beVerbose     = false;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-threshold", 2) == 0) {
      intronLimit = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-quiet", 2) == 0) {
      fprintf(stderr, "QUIET MODE ENABLED -- non-modified matches not output!\n");
      filter = false;

    } else if (strncmp(argv[arg], "-beforeafter", 2) == 0) {
      fprintf(stderr, "DEBUG MODE ENABLED -- many 'spl.*' files created!\n");
      beforeafter  = true;

    } else if (strncmp(argv[arg], "-segregate", 3) == 0) {
      fprintf(stderr, "SEGREGATION MODE ENABLED -- many 'filt.*' files created!\n");
      segregate     = true;

    } else if (strncmp(argv[arg], "-gff3", 5) == 0) {
      style = sim4polishGFF3;

    } else if (strncmp(argv[arg], "-savejunk", 3) == 0) {
      saveJunk = true;

    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      beVerbose = true;

    } else {
      err++;
    }

    arg++;
  }
  if ((err) ||
      (isatty(fileno(stdin))) ||
      (isatty(fileno(stdout)) && filter)) {
    fprintf(stderr, "usage: %s [-threshold t] [-savejunk] [-quiet] [-debug]\n", argv[0]);
    fprintf(stderr, "  -threshold    Introns bigger than this are candidates for trimming (default = 100000).\n");
    fprintf(stderr, "  -quiet        Don't print unmodified matches\n");
    fprintf(stderr, "  -beforeafter  Save (in separate files) the before/after of each modified match\n");
    fprintf(stderr, "  -segregate    Save (in separate files) the after of each modified match\n");
    fprintf(stderr, "  -savejunk     Also print the trimmed pieces (as separate matches)\n");

    if (isatty(fileno(stdin)))
      fprintf(stderr, "error: I cannot read polishes from the terminal!\n");

    if (isatty(fileno(stdout)) && filter)
      fprintf(stderr, "error: Please redirect the polishes (stdout) to a file.\n");

    exit(1);
  }

  if (beVerbose)
    fprintf(stderr, "A big intron is one that is at least "u32bitFMT"bp long.\n", intronLimit);

  if (beforeafter) {
#if 0
    splGood      = new sim4polishWriter("spl.good",      style);
    splProbGood  = new sim4polishWriter("spl.probGood",  style);
#endif
    splJunkLeft  = new sim4polishWriter("spl.junkLeft",  style);
    splJunkRight = new sim4polishWriter("spl.junkRight", style);
    splJunkBoth  = new sim4polishWriter("spl.junkBoth",  style);
    splIntronGap = new sim4polishWriter("spl.intronGap", style);
  }

  if (segregate) {
#if 0
    filtOne       = new sim4polishWriter("filt.filtOne",   style);
    filtAllSmall  = new sim4polishWriter("filt.allSmall",  style);
#endif
    filtGood      = new sim4polishWriter("filt.good",      style);
    filtProbGood  = new sim4polishWriter("filt.probGood",  style);
    filtJunkLeft  = new sim4polishWriter("filt.junkLeft",  style);
    filtJunkRight = new sim4polishWriter("filt.junkRight", style);
    filtJunkBoth  = new sim4polishWriter("filt.junkBoth",  style);
    filtIntronGap = new sim4polishWriter("filt.intronGap", style);
  }

  sim4polishWriter *W = new sim4polishWriter("-", style);
  sim4polishReader *R = new sim4polishReader("-");
  sim4polish       *p = 0L;

  if (R->getsim4polishStyle() != style)
    fprintf(stderr, "warning: input format and output format differ.\n");

  while (R->nextAlignment(p)) {
    u32bit exA;
    u32bit exB;

    if (p->_numExons == 1) {
      oneExon++;
      if (filter)
        W->writeAlignment(p);
#if 0
      if (segregate)
        filtOneExon->writeAlignment(p);
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
          W->writeAlignment(p);
#if 0
        if (segregate)
          filtAllSmall->writeAlignment(p);
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
            W->writeAlignment(a);
            W->writeAlignment(b);
          }

          if (beforeafter) {
            //fprintf(splJunkBoth, "====================\n");
            splJunkBoth->writeAlignment(p);
            splJunkBoth->writeAlignment(a);
            splJunkBoth->writeAlignment(b);
          }

          if (segregate) {
            filtJunkBoth->writeAlignment(a);
            filtJunkBoth->writeAlignment(b);
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
              W->writeAlignment(a);
            W->writeAlignment(b);
          }

          if (beforeafter) {
            //fprintf(splJunkLeft, "====================\n");
            splJunkLeft->writeAlignment(p);
            splJunkLeft->writeAlignment(a);
            splJunkLeft->writeAlignment(b);
          }

          if (segregate) {
            filtJunkLeft->writeAlignment(a);
            filtJunkLeft->writeAlignment(b);
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
            W->writeAlignment(a);
            if (saveJunk)
              W->writeAlignment(b);
          }

          if (beforeafter) {
            //fprintf(splJunkRight, "====================\n");
            splJunkRight->writeAlignment(p);
            splJunkRight->writeAlignment(a);
            splJunkRight->writeAlignment(b);
          }

          if (segregate) {
            filtJunkRight->writeAlignment(a);
            filtJunkRight->writeAlignment(b);
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
              W->writeAlignment(a);
              W->writeAlignment(b);
            }

            if (beforeafter) {
              //fprintf(splIntronGap, "====================\n");
              splIntronGap->writeAlignment(p);
              splIntronGap->writeAlignment(a);
              splIntronGap->writeAlignment(b);
            }

            if (segregate) {
              filtIntronGap->writeAlignment(a);
              filtIntronGap->writeAlignment(b);
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
                W->writeAlignment(p);
              if (segregate)
                filtGood->writeAlignment(p);
            } else {
              probGood++;
              if (filter)
                W->writeAlignment(p);
              if (segregate)
                filtProbGood->writeAlignment(p);
            }
          }
        }

      }  //  Has a big intron
    }  //  More than one exon

    totMatches++;
  }

  delete R;
  delete W;

  if (beforeafter) {
#if 0
    delete splGood;
    delete splProbGood;
#endif
    delete splJunkLeft;
    delete splJunkRight;
    delete splJunkBoth;
    delete splIntronGap;
  }

  if (segregate) {
#if 0
    delete filtOne;
    delete filtAllSmall;
#endif
    delete filtGood;
    delete filtProbGood;
    delete filtJunkLeft;
    delete filtJunkRight;
    delete filtJunkBoth;
    delete filtIntronGap;
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
