// This file is part of A2Amapper.
// Copyright (c) 2006 J. Craig Venter Institute
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
#include <unistd.h>

#include "bio++.H"
#include "atac.H"


void
usage(char *name) {
  fprintf(stderr, "usage: %s [] -m matches -l log\n", name);
  fprintf(stderr, "   When it works, fill this in...\n");
}


int
main(int argc, char *argv[]) {
  char         *matchesFile = 0L;
  FILE         *logFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-l") == 0) {
      errno = 0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open logfile '%s': %s\n", strerror(errno), argv[arg]), exit(1);
    } else {
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
    usage(argv[0]), exit(1);
  if (logFile == 0L)
    usage(argv[0]), exit(1);

  atacFile       AF(matchesFile);
  atacMatchList &ML = *AF.matches();
  atacMatchOrder MO(ML);

  //  Sort by either axis.
  MO.sortA();

  //  We need to compute the identity of the gap; our metric (thanks to Nelson) is
  //  if ("long" and "not low identity") or ("short"), close the gap

  //  We could use the seqCache, but with only a handful of gaps, we
  //  just let the OS cache stuff.

  seqCache           *C1 = new seqCache(AF.assemblyFileA(),    2, false);
  seqCache           *C2 = new seqCache(AF.assemblyFileB(), 1024, false);

  seqInCore  *S1 = 0L;
  seqInCore  *S2 = 0L;

  for (u32bit iter=0; iter<10; iter++) {
    u32bit  gapsize = 1000;
    u32bit  fgaps = 0;
    u32bit  rgaps = 0;

    int mergeuid = 1;

    for (u32bit i=1; i<MO.numMatches(); i++) {
      atacMatch  *l = MO[i-1];
      atacMatch  *r = MO[i];

      bool  joinMatches = false;

      u32bit  gap1 = 0;
      u32bit  gap2 = 0;

      if ((l->iid1 == r->iid1) &&  //  Matches are between the same sequences
          (l->iid2 == r->iid2) &&
          (l->fwd2 == r->fwd2)) {  //  Matches are the same orientation
        

        if (l->fwd2 == true) {
          if ((l->pos1 + l->len1 <= r->pos1) &&  //  Matches are ordered correctly (should be, from the sort)
              (l->pos2 + l->len2 <= r->pos2)) {

            gap1 = r->pos1 - (l->pos1 + l->len1);
            gap2 = r->pos2 - (l->pos2 + l->len2);

            if ((gap1 == gap2) &&
                (gap1 <= gapsize)) {
              S1 = C1->getSequenceInCore(l->iid1);
              S2 = C2->getSequenceInCore(l->iid2);

              char *s1 = S1->sequence() + l->pos1 + l->len1;
              char *s2 = S2->sequence() + l->pos2 + l->len1;

              u32bit  identities = 0;
              u32bit  n1 = 0;
              u32bit  n2 = 0;
              for (u32bit p=0; p<gap1; p++) {
                if (toUpper[(int)s1[p]] == toUpper[(int)s2[p]])
                  identities++;
                if (toUpper[(int)s1[p]] == 'N')
                  n1++;
                if (toUpper[(int)s2[p]] == 'N')
                  n2++;
              }

              if ((100*n1 < 20*gap1) &&            //  gap is not N and
                  (100*n2 < 20*gap1) &&            //  gap is not N and
                  ((100*identities < 80*gap1) ||   //  (gap is high identity
                   ((gap1 < 11) && (2*gap1 <= l->len1) && (2*gap1 <= r->len1)) ||   //   (gap is short, and the flanks are big
                   ((gap1 < 11) && (100*identities < 90*gap1)))) {                  //   (gap is short and high quality

                //  ALSO need to check that the gap is not actually
                //  mapped in sequence 2.  Not really, just make sure
                //  these two matches are in the same run.
                //
                if (strcmp(l->parentuid, r->parentuid) != 0) {
                  fprintf(logFile, "HEY!  F gap of size "u32bitFMT" not in a run?\n", gap1);
                  l->print(logFile, AF.labelA(), AF.labelB());
                  r->print(logFile, AF.labelA(), AF.labelB());
                } else {
                  fgaps++;

                  joinMatches = true;
                  
                  //fprintf(logFile, "potential f gap of size L "u32bitFMTW(4)" (n1="u32bitFMTW(4)" n2="u32bitFMTW(4)" ident="u32bitFMTW(4)"/"u32bitFMTW(4)")!\n",
                  //        gap1, n1, n2, identities, gap1);
                  //l->print(logFile, AF.labelA(), AF.labelB());
                  //r->print(logFile, AF.labelA(), AF.labelB());
                }
              }
            }
          }
        }  //  was a forward match



        if (l->fwd2 == false) {
          if ((l->pos1 + l->len1 <= r->pos1) &&  //  Matches are ordered correctly (should be, from the sort)
              (r->pos2 + r->len2 <= l->pos2)) {

            gap1 = r->pos1 - (l->pos1 + l->len1);
            gap2 = l->pos2 - (r->pos2 + r->len2);

            if ((gap1 == gap2) &&
                (gap1 <= gapsize)) {

              S1 = C1->getSequenceInCore(l->iid1);
              S2 = C2->getSequenceInCore(l->iid2);

              char *s1 = S1->sequence() + l->pos1 + l->len1;
              char *s2 = S2->sequence() + r->pos2 + r->len2;

              u32bit  identities = 0;
              u32bit  n1 = 0;
              u32bit  n2 = 0;
              for (u32bit p=0, q=gap1-1; p<gap1; p++, q--) {
                if (toUpper[(int)s1[p]] == toUpper[complementSymbol[(int)s2[q]]])
                  identities++;
                if (toUpper[(int)s1[p]] == 'N')
                  n1++;
                if (toUpper[(int)s2[q]] == 'N')
                  n2++;
              }

              //  Gap is short, flanks are big
              //  Gap is short, flanks are short and gap is good

              if ((100*n1 < 20*gap1) &&                                             //  gap is not N and
                  (100*n2 < 20*gap1) &&                                             //  gap is not N and
                  ((100*identities < 80*gap1) ||                                    //  (gap is high identity
                   ((gap1 < 11) && (2*gap1 <= l->len1) && (2*gap1 <= r->len1)) ||   //   (gap is short, and the flanks are big
                   ((gap1 < 11) && (100*identities < 90*gap1)))) {                  //   (gap is short and high quality

                //  ALSO need to check that the gap is not actually
                //  mapped in sequence 2.  Not really, just make sure
                //  these two matches are in the same run.
                //
                if (strcmp(l->parentuid, r->parentuid) != 0) {
                  fprintf(logFile, "HEY!  R gap of size "u32bitFMT" not in a run?\n", gap1);
                  l->print(logFile, AF.labelA(), AF.labelB());
                  r->print(logFile, AF.labelA(), AF.labelB());
                } else {
                  rgaps++;

                  joinMatches = true;

                  //fprintf(logFile, "potential r gap of size L "u32bitFMTW(4)" (n1="u32bitFMTW(4)" n2="u32bitFMTW(4)" ident="u32bitFMTW(4)"/"u32bitFMTW(4)")!\n",
                  //        gap1, n1, n2, identities, gap1);
                  //l->print(logFile, AF.labelA(), AF.labelB());
                  //r->print(logFile, AF.labelA(), AF.labelB());
                }
              }


            }
          }
        }
      }

      if (joinMatches) {
        fprintf(logFile, "CLOSE "u32bitFMT"----------------------------------------\n", gap1);
        l->print(logFile, AF.labelA(), AF.labelB());
        r->print(logFile, AF.labelA(), AF.labelB());

        MO.mergeMatches(l, r, mergeuid);

        l->print(logFile, AF.labelA(), AF.labelB());

        mergeuid++;
        i--;
      }
    }

    fprintf(logFile, "At gapSize="u32bitFMT" closed "u32bitFMT" f-gaps and "u32bitFMT" r-gaps.\n", gapsize, fgaps, rgaps);

    if (fgaps + rgaps == 0)
      iter = 10;
  }


#if 0
  //  This analyzes an atac mapping, looking for a signature that indicates a bad
  //  alignment.  If we have an alignment of:
  //      XXXXXXaC-YYYYYY
  //      XXXXXX-CtYYYYYY
  //  this will generate three matches, instead of one match with mismatches in it.
  //  We scan the FORWARD matches for this pattern, and report any we find.
  //
  //  We only found 3 on huref4 vs b35.  Further development here was stopped.

  for (u32bit i=2; i<ML.numMatches(); i++) {
    atacMatch  *l = ML[i-2];
    atacMatch  *m = ML[i-1];
    atacMatch  *r = ML[i];

    if (m->len1 < 3) {  //  The match in the middle is small

      if ((l->iid1 == r->iid1) &&  //  Matches are between the same sequences
          (l->iid2 == r->iid2) &&
          (l->fwd2 == r->fwd2)) {  //  Matches are the same orientation

        if (l->fwd2 == true) {
          if ((l->pos1 + l->len1 <= r->pos1) &&  //  Matches are ordered correctly (should be, from the sort)
              (l->pos2 + l->len2 <= r->pos2)) {

            u32bit  gapl1 = m->pos1 - (l->pos1 + l->len1);
            u32bit  gapl2 = m->pos2 - (l->pos2 + l->len2);
            u32bit  gapr1 = r->pos1 - (m->pos1 + m->len1);
            u32bit  gapr2 = r->pos2 - (m->pos2 + m->len2);

            if ((gapl1 + gapr1 == gapl2 + gapr2) && (gapl1 + gapr1 < 5)) {
              fprintf(logFile, "potential f fix of size L "u32bitFMT" "u32bitFMT" and R "u32bitFMT" "u32bitFMT"!\n",
                      gapl1, gapl2, gapr1, gapr2);
              l->print(logFile, "A", "B");
              m->print(logFile, "A", "B");
              r->print(logFile, "A", "B");
            }
          } else {
            fprintf(logFile, "sort is forward broken.\n");
          }
        }  //  was a forward match
      }
    }
  }
#endif


  //  Write the new output to stdout -- we preserve runs here, but
  //  discard everything else.
  //
  AF.writeHeader(stdout);

  for (u32bit i=0; i<MO.numMatches(); i++)
    MO[i]->print(stdout, AF.labelA(), AF.labelB());

  for (u32bit i=0; i<AF.runs()->numberOfMatches(); i++)
    AF.runs()->getMatch(i)->print(stdout, AF.labelA(), AF.labelB());

  return(0);
}
