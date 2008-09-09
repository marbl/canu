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

//  Reads a set of atac matches, trims off ends that are mismatch.
//  Computes the percent identity of the resulting match.
//  Outputs the trimmed match if it is above some percent identity.

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-d identity] [-i identity] -m matches\n", name);
  fprintf(stderr, "  -d    discard the match if it is below this percent identity\n");
}

int
main(int argc, char *argv[]) {
  char         *matchesFile      = 0L;
  double        discardThreshold = 0.0;
  u32bit        discardLength    = 0;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      matchesFile = argv[++arg];
    } else if (strcmp(argv[arg], "-d") == 0) {
      discardThreshold = atof(argv[++arg]);
      if (discardThreshold > 1.0)
        discardThreshold /= 100;
    } else if (strcmp(argv[arg], "-l") == 0) {
      discardLength = atoi(argv[++arg]);
    } else {
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (matchesFile == 0L)
    usage(argv[0]), exit(1);

  atacFile       AF(matchesFile);
  atacMatchList &ML = *AF.matches();
  seqCache     Acache(AF.assemblyFileA(), 32, false);
  seqCache     Bcache(AF.assemblyFileB(), 32, false);

  for (u32bit i=0; i<ML.numMatches(); i++) {
    atacMatch            *m = ML.getMatch(i);

    u32bit identities = 0;

    //char   *a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
    //char   *b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;
    //u32bit  p, q;


    //  Trim the match
    //
    if (m->fwd2) {
      char   *a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
      char   *b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;
      u32bit  p = 0;

      while ((m->len1 > 0) && (toUpper[(int)a[p]] != toUpper[(int)b[p]])) {
        m->pos1++;
        m->pos2++;
        m->len1--;
        m->len2--;
        p++;
      }

      a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
      b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;
      p = m->len1-1;
      while ((m->len1 > 0) && (toUpper[(int)a[p]] != toUpper[(int)b[p]])) {
        m->len1--;
        m->len2--;
        p--;
      }

    } else {
      char   *a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
      char   *b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;
      u32bit  p = 0;
      u32bit  q = m->len2 - 1;

      while ((m->len1 > 0) && (toUpper[(int)a[p]] != complementSymbol[toUpper[(int)b[q]]])) {
        m->pos1++;
        m->len1--;
        m->len2--;
        p++;
        q--;
      }

      a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
      b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;
      p = m->len1 - 1;
      q = 0;
      while ((m->len1 > 0) && (toUpper[(int)a[p]] != complementSymbol[toUpper[(int)b[q]]])) {
        m->len1--;
        m->pos2++;
        m->len2--;
        p--;
        q++;
      }
    }

    if (m->len1 > 0) {
      char *a = Acache.getSequenceInCore(m->iid1)->sequence() + m->pos1;
      char *b = Bcache.getSequenceInCore(m->iid2)->sequence() + m->pos2;

      if (m->fwd2) {
        for (u32bit p=0; p<m->len1; p++) {
          if (toUpper[(int)a[p]] == toUpper[(int)b[p]])
            identities++;
        }
      } else {
        for (u32bit p=0, q=m->len2-1; p<m->len1; p++, q--) {
          if (toUpper[(int)a[p]] == toUpper[complementSymbol[(int)b[q]]])
            identities++;
        }
      }

      double   myIdentity = (double)identities / m->len1;

      if ((myIdentity > discardThreshold) && (m->len1 > discardLength)) {
        m->print(stdout, AF.labelA(), AF.labelB());
      } else {
        fprintf(stderr, "match "u32bitFMT" is only %6.2f%% identity and "u32bitFMT" long:  ",
                i, 100.0 * identities / m->len1, m->len1);
        m->print(stderr, AF.labelA(), AF.labelB());
        if (m->len1 < 200) {
          char   tmp[1000];

          strncpy(tmp, a, m->len1);
          tmp[m->len1] = 0;
          fprintf(stderr, "  %s\n", tmp);

          strncpy(tmp, b, m->len1);
          tmp[m->len1] = 0;
          fprintf(stderr, "  %s\n", tmp);
        }
      }
    }
  }

  return(0);
}
