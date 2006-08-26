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

//  Reads an atac mapping, and a list of features.  Features on one
//  axis are projected to the other axis using the atac map.

int
main(int argc, char **argv) {
  char         *matchesFile = 0L;
  char         *featureFile = 0L;
  FILE         *logFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-f") == 0) {
      featureFile = argv[++arg];
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
  if (featureFile == 0L)
    usage(argv[0]), exit(1);
  if (logFile == 0L)
    usage(argv[0]), exit(1);

  atacMatchList    ML(matchesFile, 'm');
  atacMatchOrder   MO(ML);
  atacFeatureList  FL(featureFile);

  //  Project features from A to B.
  MO.sortA();
  FL.sort();

  u32bit  mid = 0;
  u32bit  fid = 0;
  u32bit  pid = 0;

  while ((mid < MO.numberOfMatches()) &&
         (fid < FL.numberOfFeatures())) {
    atacMatch    *m = MO[mid];
    atacFeature  *f = FL[fid];

    if (m->iid1 < f->iid) {
      mid++;
      continue;
    }
    if (f->iid < m->iid1) {
      fid++;
      continue;
    }

    //  Same sequences now!

    if (m->pos1 + m->len1 < f->pos) {
      //  match ends before the feature
      mid++;
      continue;
    }

    if (f->pos + f->len < m->pos1) {
      //  Feature begins before match
      fid++;
      continue;
    }

    //  Feature and match now overlap!


    //
    //  This does A -> B -- ONLY.
    //



    //  If feature is completely in match, this is easy.
    //
    if ((m->pos1 <= f->pos) && ((f->pos + f->len) <= (m->pos1 + m->len1))) {
      u32bit beg;

      if (m->fwd2 == true) {
        beg = m->pos2 + f->pos - m->pos1;
      } else {
        beg = m->pos2 + m->len2 - (f->pos - m->pos1) - f->len;
      }

      if (f->len > 0)
        fprintf(stdout, "M u Aprojected"u32bitFMT" %s.%s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
                pid,
                f->featureuid, m->matchuid,
                ML.labelA(), f->iid, f->pos, f->len,
                ML.labelB(), m->iid2, beg, f->len, (m->fwd2) ? 1 : -1);
      pid++;
      fid++;
      continue;
    }

    //  If match is completely within feature, super easy!
    //
    if ((f->pos < m->pos1) && (m->pos1 + m->len1) < (f->pos + f->len)) {
      if (m->len1 > 0)
        fprintf(stdout, "M u Bprojected"u32bitFMT" %s.%s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
                pid,
                f->featureuid, m->matchuid,
                ML.labelA(), m->iid1, m->pos1, m->len1,
                ML.labelB(), m->iid2, m->pos2, m->len2, (m->fwd2) ? 1 : -1);
      pid++;
      fid++;
      continue;
    }


    //  Dang, feature isn't completely in match.  Guess where feature
    //  could be ending?  Or just project as much as possible?

    if (f->pos < m->pos1) {
      u32bit len = f->len - (m->pos1 - f->pos);
      u32bit beg;

      if (m->fwd2 == true) {
        beg = m->pos2;
      } else {
        beg = m->pos2 + m->len2 - len;
      }

      if (len > 0)
        fprintf(stdout, "M u Cprojected"u32bitFMT" %s.%s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
                pid,
                f->featureuid, m->matchuid,
                ML.labelA(), f->iid, m->pos1, len,
                ML.labelB(), m->iid2, beg, len, (m->fwd2) ? 1 : -1);
      pid++;
      fid++;
      continue;
    }

    if (m->pos1 + m->len1 < f->pos + f->len) {
      u32bit len = m->pos1 + m->len1 - f->pos;
      u32bit beg;

      if (m->fwd2 == true) {
        beg = m->pos2 + m->len2 - len;
      } else {
        beg = m->pos2;
      }

      if (len > 0)
        fprintf(stdout, "M u Dprojected"u32bitFMT" %s.%s %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" 1 %s:"u32bitFMT" "u32bitFMT" "u32bitFMT" %d\n",
                pid,
                f->featureuid, m->matchuid,
                ML.labelA(), f->iid, f->pos, len,
                ML.labelB(), m->iid2, beg, len, (m->fwd2) ? 1 : -1);
      pid++;
      fid++;
      continue;
    }

    fprintf(stderr, "projectFeatures:  Unhandled case?\n");
    m->print(stdout, "A", "B");
    f->print(stdout, "A");

    assert(0);
  }
}
