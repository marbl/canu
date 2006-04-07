#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "sim4.H"

//  Writes polished from stdin as atac-format matches.  Splits polishes on any indel to generate gapless
//  atac matches (type 'u').
//
//  Does no cleanup.


void
indelRedo(char *a, char *b) {
  u32bit  orig = 0;
  u32bit  copy = 0;

  while (a[orig] && b[orig]) {
    if ((a[orig] != '-') ||
        (b[orig] != '-')) {
      if (orig != copy) {
        a[copy] = a[orig];
        b[copy] = b[orig];
      }
      copy++;
    }
    orig++;
  }

  a[copy] = 0;
  b[copy] = 0;
}


u32bit
indelFixAlignment(char *a, char *b) {
  bool    redo  = false;
  u32bit  len   = strlen(a) - 1;
  u32bit  fixed = 0;

  //fprintf(stdout, "fixIndel\n");
  //fprintf(stdout, "%s\n%s\n", a, b);

  for (u32bit i=2; i<len; i++) {

    //  -Ac
    //  cA-  two gaps -> two mismatches
    if ((a[i-2] == '-') && (b[i] == '-')) {
      a[i-2] = toUpper[a[i-1]];   a[i-1] = toUpper[a[i]];      a[i] = '-';
      b[i-2] = toUpper[b[i-2]];   b[i-1] = toUpper[b[i-1]];    b[i] = '-';
      fixed++;
      redo = true;
    }

    if ((a[i] == '-') && (b[i-2] == '-')) {
      a[i-2] = toUpper[a[i-2]];   a[i-1] = toUpper[a[i-1]];    a[i] = '-';
      b[i-2] = toUpper[b[i-1]];   b[i-1] = toUpper[b[i]];      b[i] = '-';
      fixed++;
      redo = true;
    }
  }

  if (redo) {
    //fprintf(stdout, "%s\n%s\n", a, b);
    //fprintf(stdout, "Fixed "u32bitFMT" 1 base wide indel\n", fixed);
    indelRedo(a, b);
  }

  redo  = false;
  len   = strlen(a) - 1;

  for (u32bit i=3; i<len; i++) {

    //  cAt-
    //  -Agg  two gaps, one mismatch -> three mismatches
    //        we also would do two gaps -> three mismatches
    if ((a[i] == '-') && (b[i-3] == '-')) {
      a[i-3] = toUpper[a[i-3]];   a[i-2] = toUpper[a[i-2]];   a[i-1] = toUpper[a[i-1]];    a[i] = '-';
      b[i-3] = toUpper[b[i-2]];   b[i-2] = toUpper[b[i-1]];   b[i-1] = toUpper[b[i]];      b[i] = '-';
      fixed++;
      redo = true;
    }

    if ((a[i-3] == '-') && (b[i] == '-')) {
      a[i-3] = toUpper[a[i-2]];   a[i-2] = toUpper[a[i-1]];   a[i-1] = toUpper[a[i]];      a[i] = '-';
      b[i-3] = toUpper[b[i-3]];   b[i-2] = toUpper[b[i-2]];   b[i-1] = toUpper[b[i-1]];    b[i] = '-';
      fixed++;
      redo = true;
    }
  }

  if (redo) {
    //fprintf(stdout, "%s\n%s\n", a, b);
    //fprintf(stdout, "Fixed "u32bitFMT" 2 base wide indel\n", fixed);
    indelRedo(a, b);
  }

  return(fixed);
}



int
main(int argc, char **argv) {
  char  *nickname1 = 0L;
  char  *nickname2 = 0L;
  bool   flip      = false;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-1", 2) == 0) {
      nickname1 = argv[++arg];
    } else if (strncmp(argv[arg], "-2", 2) == 0) {
      nickname2 = argv[++arg];
    } else if (strncmp(argv[arg], "-f", 2) == 0) {
      flip = true;
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((nickname1 == 0L) || (nickname2 == 0L)) {
    fprintf(stderr, "usage: %s [-f] -1 nickname1 -2 nickname2 < matches.sim4db > matches.atac\n", argv[0]);
    exit(1);
  }

  u32bit   dupRecordIID = 0;
  u32bit   dupParentIID = 0;

  u32bit   totalFixed   = 0;

  while (!feof(stdin)) {
    sim4polish *p = s4p_readPolish(stdin);

    if (p != 0L) {

      //  Parse the defline to find the genomic region our 'est'
      //  (unfortunate sim4db term) is from.  Search for our
      //  information in the defline
      //
      //    extracted from iid (\d+) pos (\d+) (\d+) 

      splitToWords  W(p->estDefLine);

      u32bit   i=0;
      while ((i < W.numWords()) && (strcmp(W[i], "iid") != 0))
        i++;
      if ((i == 0) || (i == W.numWords()))
        fprintf(stderr, "Failed to match est defline '%s'\n", p->estDefLine), exit(1);

      u32bit  qSeqIID = strtou32bit(W[i+1], 0L);
      u32bit  qSeqBeg = strtou32bit(W[i+3], 0L);
      u32bit  qSeqEnd = strtou32bit(W[i+4], 0L);  //  Not used


      W.split(p->genDefLine);

      i=0;
      while ((i<W.numWords()) && (strcmp(W[i], "iid") != 0))
        i++;
      if ((i == 0) || (i == W.numWords()))
        fprintf(stderr, "Failed to match gen defline '%s'\n", p->genDefLine), exit(1);

      u32bit  gSeqIID = strtou32bit(W[i+1], 0L);
      u32bit  gSeqBeg = strtou32bit(W[i+3], 0L);
      u32bit  gSeqEnd = strtou32bit(W[i+4], 0L);  //  Not used

      bool    fwd  = (p->matchOrientation == SIM4_MATCH_FORWARD);


      //  Fix the coords
      //
      if (fwd) {
        //  Forward is easy!  Just add.

        for (u32bit exon=0; exon<p->numExons; exon++) {
          sim4polishExon *e = p->exons + exon;

          e->estFrom += qSeqBeg;
          e->estTo   += qSeqBeg;
          e->genFrom += gSeqBeg + p->genLo;
          e->genTo   += gSeqBeg + p->genLo;
        }

        p->genLo = 0;
      } else {
        //  Reverse is not easy.  Need to reverse complement the query positions.

        for (u32bit exon=0; exon<p->numExons; exon++) {
          sim4polishExon *e = p->exons + exon;

          //  First, reverse the query relative to our extracted piece
          //
          u32bit f = (qSeqEnd - qSeqBeg) - e->estTo   + 2;  //  Extra +1 to offset -1 when we set qBeg
          u32bit t = (qSeqEnd - qSeqBeg) - e->estFrom + 2;

          //  Now we can just offset stuff.
          e->estFrom  = qSeqBeg + t;  //  Really the end!
          e->estTo    = qSeqBeg + f;  //  Really the begin!
          e->genFrom += gSeqBeg + p->genLo;
          e->genTo   += gSeqBeg + p->genLo;
        }

        p->genLo = 0;
      }



      for (u32bit exon=0; exon<p->numExons; exon++) {
        sim4polishExon *e = p->exons + exon;

        //  Parse the alignment to find ungapped blocks

        u32bit  aPos = 0;

        u32bit  qBeg = e->estFrom - 1;
        u32bit  gBeg = e->genFrom - 1;

        u32bit  mLen = 0;

        totalFixed += indelFixAlignment(e->estAlignment, e->genAlignment);

        bool  notDone = true;  //  There should be a way to get rid of this stupid variable....
        while (notDone) {
          notDone = ((e->estAlignment[aPos] != 0) &&
                     (e->genAlignment[aPos] != 0));

          //  If we find the end of a gapless block, emit a match

          if ((e->estAlignment[aPos] == '-') || (e->estAlignment[aPos] == 0) ||
              (e->genAlignment[aPos] == '-') || (e->genAlignment[aPos] == 0)) {

            //  If there is an indel at the start (which probably
            //  shouldn't happen anyway!), or possibly at the end,
            //  then our length is zero, and we should not emit
            //  anything.
            //
            if (mLen > 0) {
              if (flip == false) {
                fprintf(stdout, "M u dupr"u32bitFMT" dupp"u32bitFMT" %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %s\n",
                        dupRecordIID,
                        dupParentIID,
                        //nickname1, qSeqIID, qBeg, mLen,
                        nickname1, qSeqIID, (fwd) ? qBeg : qBeg - mLen, mLen,
                        nickname2, gSeqIID, gBeg, mLen,
                        (fwd) ? "1" : "-1");
              } else {
                fprintf(stdout, "M u dupr"u32bitFMT" dupp"u32bitFMT" %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %s\n",
                        dupRecordIID,
                        dupParentIID,
                        nickname2, gSeqIID, gBeg, mLen,
                        //nickname1, qSeqIID, qBeg, mLen,
                        nickname1, qSeqIID, (fwd) ? qBeg : qBeg - mLen, mLen,
                        (fwd) ? "1" : "-1");
              }
              dupRecordIID++;

              //  Adjust our begin and end positions to the end of this record
              if (fwd)  qBeg += mLen;
              else      qBeg -= mLen;
              gBeg += mLen;

              mLen  = 0;
            }

            //  Skip whatever caused us to emit a gapless block
            while ((e->estAlignment[aPos] == '-') ||
                   (e->genAlignment[aPos] == '-')) {
              if (e->estAlignment[aPos] == '-')
                gBeg++;
              if (e->genAlignment[aPos] == '-')
                if (fwd) qBeg++;
                else     qBeg--;
              aPos++;
            }
          } else {
            //  Not the end of a gapless block, extend this match by one
            mLen++;
            aPos++;
          }

        }  //  over all positions in the alignemnt
      }  //  over all exons

      dupParentIID++;
      s4p_destroyPolish(p);
    }
  }

  fprintf(stderr, "Fixed "u32bitFMT" indel/mismatches.\n", totalFixed);

  return(0);
}

