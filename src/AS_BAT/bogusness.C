
/**************************************************************************
 * Copyright (C) 2010, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: bogusness.C,v 1.3 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_bogusUtil.H"

//  Reads snapper/nucmer output, figures out what the minimal set of alignments that cover each
//  aligned sequence, then compares those alignments against the ideal unitigs.
//
//  It reports:
//    bubbles that should have been popped.
//    unitigs that span a unique-repeat junction (one end in unique, one end in repeat).
//    unitigs that include a repeat.
//    unitigs that include a repeat, and are chimeric at the repeat.
//    unitigs that should be joined.
//
//

#define IDEAL_REPT      0  //  A repeat unitig
#define IDEAL_UNIQ      1  //  A unique unitig
#define IDEAL_REPTUNIQ  2  //  A portion of this repeat unitig might be separable
#define IDEAL_UNIQWEAK  3  //  A portion of this unique unitig might be unjoinable

char *types[4] = { "REPT", "UNIQ", "SEPR", "WEAK" };

class idealUnitig {
public:
  idealUnitig() {
    bgn  = 0;
    end  = 0;
    type = 0;
  };
  idealUnitig(int32 b, int32 e, char t) {
    bgn  = b;
    end  = e;
    type = t;
  };
  ~idealUnitig() {
  };
  int32  bgn;
  int32  end;
  char   type;
};

void
loadIdealUnitigs(char *idealName,
                 vector<idealUnitig>  &ideal) {

  errno = 0;
  FILE *F = fopen(idealName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", idealName, strerror(errno)), exit(1);

  char          L[1024];
  splitToWords  S;
  char          t = 0;

  fgets(L, 1024, F);

  while (!feof(F)) {
    chomp(L);

    S.split(L);

    if        (strcmp(S[2], "REPT") == 0)
      t = (S.numWords() == 4) ? IDEAL_REPT : IDEAL_UNIQWEAK;

    else if (strcmp(S[2], "UNIQ") == 0)
      t = (S.numWords() == 4) ? IDEAL_UNIQ : IDEAL_REPTUNIQ;

    else
      fprintf(stderr, "Unknown type in '%s'\n", L), exit(1);

    ideal.push_back(idealUnitig(S(0), S(1), t));

    fgets(L, 1024, F);
  }

  fclose(F);

  fprintf(stderr, "Loaded %lu ideal unitigs.\n", ideal.size());
}





bool
isUnitigContained(int32  frgIID,
                  idealUnitig  &ideal,
                  genomeAlignment  &cover,
                  int32  &alen, int32 &ulen, int32 &ilen,
                  double &ufrac,
                  double &ifrac) {

  if ((ideal.bgn    <= cover.genBgn) &&
      (cover.genEnd <= ideal.end)) {
    alen  = cover.genEnd - cover.genBgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigContaining(int32  frgIID,
                   idealUnitig  &ideal,
                   genomeAlignment  &cover,
                   int32  &alen, int32 &ulen, int32 &ilen,
                   double &ufrac,
                   double &ifrac) {

  if ((cover.genBgn <= ideal.bgn) &&
      (ideal.end    <= cover.genEnd)) {
    alen  = ideal.end - ideal.bgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigEnding(int32  frgIID,
               idealUnitig  &ideal,
               genomeAlignment  &cover,
               int32  &alen, int32 &ulen, int32 &ilen,
               double &ufrac,
               double &ifrac) {

  if ((cover.genBgn <= ideal.bgn) &&
      (ideal.bgn    <= cover.genEnd)) {
    alen  = cover.genEnd - ideal.bgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}


bool
isUnitigBeginning(int32  frgIID,
                  idealUnitig  &ideal,
                  genomeAlignment  &cover,
                  int32  &alen, int32 &ulen, int32 &ilen,
                  double &ufrac,
                  double &ifrac) {

  if ((cover.genBgn <= ideal.end) &&
      (ideal.end    <= cover.genEnd)) {
    alen  = ideal.end - cover.genBgn;
    ulen  = cover.genEnd - cover.genBgn;
    ilen  = ideal.end - ideal.bgn;
    ufrac = 100.0 * alen / ulen;
    ifrac = 100.0 * alen / ilen;
    return(true);
  }

  return(false);
}





int
main(int argc, char **argv) {
  char  *nucmerName   = 0L;
  char  *snapperName  = 0L;
  char  *idealName = 0L;

  vector<idealUnitig>    ideal;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-nucmer") == 0) {
      nucmerName = argv[++arg];

    } else if (strcmp(argv[arg], "-snapper") == 0) {
      snapperName = argv[++arg];

    } else if (strcmp(argv[arg], "-ideal") == 0) {
      idealName = argv[++arg];

    } else {
      err++;
    }
    arg++;
  }
  if ((nucmerName == 0L) && (snapperName == 0L))
    fprintf(stderr, "ERROR: No input matches supplied (either -nucmer or -snapper).\n"), err++;
  if (idealName == 0L)
    fprintf(stderr, "ERROR: No ideal unitigs supplied (-ideal).\n"), err++;
  if (err) {
    exit(1);
  }


  loadIdealUnitigs(idealName, ideal);


  vector<longestAlignment>   longest;
  vector<genomeAlignment>    genome;
  map<string, int32>         IIDmap;
  vector<string>             IIDname;
  uint32                     IIDnext = 0;

  if (nucmerName)
    loadNucmer(nucmerName, longest, genome, IIDmap, IIDname, IIDnext, NULL);

  if (snapperName)
    loadSnapper(snapperName, longest, genome, IIDmap, IIDname, IIDnext, NULL);

  //  Process the matches fragment by fragment.  Find the longest, or the minimal set that spans the
  //  unitig.  Cases:
  //
  //    - one alignment covers whole unitig
  //
  //    - one alignment coverss whole unitig, some repeat matches contained
  //      - keep only that longest alignment
  //
  //    - multiple alignments cover the uniitg
  //      - keep all spanning alignments
  //
  //    - multiple alignments cover the unitig, some repeat matches contained
  //
  //  But basically, discard any alignment that is contained in some other alignment.  We're not
  //  trying to explain the genome using the unitigs, we're trying to explain (evaluate) the unitigs
  //  using the genome.


  sort(genome.begin(), genome.end(), byFragmentID);

  longest.resize(IIDnext);


  for (uint32 bgn=0, lim=genome.size(); bgn<lim; ) {
    uint32 frgIID = genome[bgn].frgIID;

    vector<genomeAlignment>  cover;

    uint32 end = bgn + 1;
    while ((end < genome.size()) &&
           (genome[bgn].frgIID == genome[end].frgIID))
      end++;

    //  Delete any alignments that are contained in some other alignment.

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  &I = genome[i];

      if ((I.frgBgn == 0) && (I.frgEnd == 0))
        //  [i] already deleted
        continue;

      for (uint32 j=bgn; j<end; j++) {
        genomeAlignment  &J = genome[j];

        if (i == j)
          //  We are not contained in ourself
          continue;

        if ((J.frgBgn == 0) && (J.frgEnd == 0))
          //  [j] already deleted
          continue;

        if ((I.frgBgn <= J.frgBgn) &&
            (J.frgEnd <= I.frgEnd))
          //  J contained in I.
          J.frgBgn = J.frgEnd = 0;
      }
    }

    //  Copy whatever is left over to a list of 'covering' alignments

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  &I = genome[i];

      if ((I.frgBgn == 0) && (I.frgEnd == 0))
        //  [i] contained
        continue;
    
      cover.push_back(I);
    }

    //  Attempt to classify

#if 0
    int32  utgBgn = cover[0].frgBgn;
    int32  utgEnd = cover[0].frgEnd;

    for (uint32 c=0; c<cover.size(); c++) {
      utgBgn = MIN(utgBgn, cover[c].frgBgn);
      utgEnd = MIN(utgEnd, cover[c].frgEnd);
    }

    if (utgEnd < 1000) {
      bgn = end;
      continue;
    }
#endif

    const char *fmt = "| %s || %d of %d || %d-%d || %d-%d || %s || %s || %05d || %d-%d || %dbp || %.2f%% || %.2f%% || %s\n";

    for (uint32 c=0; c<cover.size(); c++) {
      int32  alen  = 0;
      int32  ulen  = 0;
      int32  ilen  = 0;
      double ufrac = 0.0;
      double ifrac = 0.0;

      //  Pass one, is our unitig (region) completely 100% explained by any single ideal unitig?
      //  What usually happens is that our unitig is 100% contained in a unique ideal unitig, but
      //  there is (less than a fragment) of overlap with the adjacent repeat unitig.

      bool explained = false;

      for (uint32 i=0; i<ideal.size(); i++) {
        if (isUnitigContained(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          fprintf(stdout, fmt,
                  IIDname[frgIID].c_str(),
                  c+1, cover.size(),
                  cover[c].frgBgn, cover[c].frgEnd, cover[c].genBgn, cover[c].genEnd,
                  "CONTAINED",
                  types[ideal[i].type], i,
                  ideal[i].bgn, ideal[i].end,
                  alen, ufrac, ifrac,
                  (ideal.size() == 1) ? " (CORRECT)" : " (EXPLAINED)");
          explained = true;
          break;
        }
      }

      if (explained)
        continue;

      for (uint32 i=0; i<ideal.size(); i++) {
        if (isUnitigContained(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          fprintf(stdout, fmt,
                  IIDname[frgIID].c_str(),
                  c+1, cover.size(),
                  cover[c].frgBgn, cover[c].frgEnd, cover[c].genBgn, cover[c].genEnd,
                  "CONTAINED",
                  types[ideal[i].type], i,
                  ideal[i].bgn, ideal[i].end,
                  alen, ufrac, ifrac,
                  "");
          continue;
        }

        if (isUnitigContaining(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          fprintf(stdout, fmt,
                  IIDname[frgIID].c_str(),
                  c+1, cover.size(),
                  cover[c].frgBgn, cover[c].frgEnd, cover[c].genBgn, cover[c].genEnd,
                  "CONTAINS ",
                  types[ideal[i].type], i,
                  ideal[i].bgn, ideal[i].end,
                  alen, ufrac, ifrac,
                  "");
          continue;
        }

        if (isUnitigBeginning(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          fprintf(stdout, fmt,
                  IIDname[frgIID].c_str(),
                  c+1, cover.size(),
                  cover[c].frgBgn, cover[c].frgEnd, cover[c].genBgn, cover[c].genEnd,
                  "BEGINSin ",
                  types[ideal[i].type], i,
                  ideal[i].bgn, ideal[i].end,
                  alen, ufrac, ifrac,
                  "");
        }

        if (isUnitigEnding(frgIID, ideal[i], cover[c], alen, ulen, ilen, ufrac, ifrac)) {
          fprintf(stdout, fmt,
                  IIDname[frgIID].c_str(),
                  c+1, cover.size(),
                  cover[c].frgBgn, cover[c].frgEnd, cover[c].genBgn, cover[c].genEnd,
                  "ENDSin   ",
                  types[ideal[i].type], i,
                  ideal[i].bgn, ideal[i].end,
                  alen, ufrac, ifrac,
                  "");
        }
      }
    }

    bgn = end;
  }

  return(0);
}
