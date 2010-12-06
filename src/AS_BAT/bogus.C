
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

const char *mainid = "$Id: bogus.C,v 1.2 2010-12-06 08:03:48 brianwalenz Exp $";

#include "AS_BAT_bogusUtil.H"

int
main(int argc, char **argv) {
  char  *nucmerName   = 0L;
  char  *snapperName  = 0L;
  char  *refName = 0L;
  char  *outputPrefix   = 0L;
  char   outputName[FILENAME_MAX];
  FILE  *outputFile = 0L;

  //  When comparing coords to if two alignments span the same piece of the fragment,
  //  allow (twice) this much slop in the comparison.  E.g., Abgn +- 5 == Bbgn
  //
  int32  alignWobble  = 5;

  //  When constructing the unique coverage map, trim each read by this amount, on each end, to
  //  simulate the minimum overlap length needed by unitigger.  This amount is automagically added
  //  back in on output.
  int32  minOverlap   = 40;

  //  When testing if a repeat alignment is contained in a non-repeat alignment, the repeat must
  //  start at least this many bases from the end of the non-repeat alignment.  In other words, a
  //  fragment with a repeat in the middle can be uniquely placed (by overlaps) with only 20 bases of
  //  unique sequence on the end.
  int32  uniqEnd      = 10;


  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-nucmer") == 0) {
      nucmerName = argv[++arg];

    } else if (strcmp(argv[arg], "-snapper") == 0) {
      snapperName = argv[++arg];

    } else if (strcmp(argv[arg], "-reference") == 0) {
      refName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-wobble") == 0) {
      alignWobble = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-overlap") == 0) {
      minOverlap = atoi(argv[++arg]);

    } else {
      err++;
    }
    arg++;
  }
  if ((nucmerName == 0L) && (snapperName == 0L))
    fprintf(stderr, "ERROR: No input matches supplied (either -nucmer or -snapper).\n"), err++;
  if (refName == 0L)
    fprintf(stderr, "ERROR: No reference sequence supplied (-reference).\n"), err++;
  if (outputPrefix == 0L)
    fprintf(stderr, "ERROR: No output prefix supplied (-output).\n"), err++;
  if (err) {
    exit(1);
  }

  vector<longestAlignment>   longest;
  vector<genomeAlignment>    genome;
  map<string, int32>         IIDmap;
  vector<string>             IIDname;
  uint32                     IIDnext = 0;

  //  Load all the matches into genomeAlignment.  Generate longestAlignment for each fragment.
  //  genomeAlignment::isLognest and genomeAlignment::isRepeat are computed later.

  sprintf(outputName, "%s.inputs", outputPrefix);
  outputFile = fopen(outputName, "w");

  if (nucmerName)
    loadNucmer(nucmerName, longest, genome, IIDmap, IIDname, IIDnext, outputFile);

  if (snapperName)
    loadSnapper(snapperName, longest, genome, IIDmap, IIDname, IIDnext, outputFile);

  fclose(outputFile);


  //  Process the matches fragment by fragment.  Find the longest, count the number
  //  of duplicates of the longest match, label as repeat/unique.

  sprintf(outputName, "%s.processing", outputPrefix);
  outputFile = fopen(outputName, "w");

  sort(genome.begin(), genome.end(), byFragmentID);

  longest.resize(IIDnext);

  for (uint32 bgn=0, lim=genome.size(); bgn<lim; ) {
    uint32 frgIID = genome[bgn].frgIID;

    uint32 end = bgn + 1;
    while ((end < genome.size()) &&
           (genome[bgn].frgIID == genome[end].frgIID))
      end++;

    //fprintf(stderr, "%d bgn %d end %d\n", frgIID, bgn, end);

    //  Find the longest alignment.  Set 'num' to zero, we'll fill it out correctly in a minute.
    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  *A = &genome[i];

      if (longest[frgIID].len < A->frgEnd - A->frgBgn) {
        longest[frgIID].bgn = A->frgBgn;
        longest[frgIID].end = A->frgEnd;
        longest[frgIID].len = A->frgEnd - A->frgBgn;
        longest[frgIID].num = 0;
      }

      if (longest[frgIID].frgLen < A->frgEnd) {
        longest[frgIID].rptBgn = 0;
        longest[frgIID].rptEnd = A->frgEnd;
        longest[frgIID].frgLen = A->frgEnd;
      }
    }

    //  Now that we know the longest, count the number of matches that have approximately the same
    //  span.  If that count is one, mark all but that span as repeat, otherwise, mark everthing as
    //  repeat.

    assert(longest[frgIID].num == 0);

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  *A = &genome[i];

      A->isRepeat = true;

      if (((longest[frgIID].bgn - alignWobble <= A->frgBgn) && (A->frgBgn <= longest[frgIID].bgn + alignWobble)) &&
          ((longest[frgIID].end - alignWobble <= A->frgEnd) && (A->frgEnd <= longest[frgIID].end + alignWobble))) {
        longest[frgIID].num++;
        A->isRepeat = false;
      }
    }

    assert(longest[frgIID].num > 0);

    //  More than one longest?  Then all alignments are repeats.
    if (longest[frgIID].num > 1)
      for (uint32 i=bgn; i<end; i++)
        genome[i].isRepeat = true;

    //  Exactly one longest?  Transfer any repeats at the end of the fragment, then add new matches
    //  for those repeats.  If those new matches are spanned, they'll get thrown out, and we'll
    //  maintain a UNIQ unitig.  If they are not spanned, we'll have a REPT-in-UNIQ situation, and
    //  we should break the UNIQ unitig.

    if (longest[frgIID].num == 1) {
      genomeAlignment *LONG = 0L;

      for (uint32 i=bgn; i<end; i++) {
        genomeAlignment  *A = &genome[i];

        if (A->isRepeat == false) {
          assert(LONG == 0L);
          LONG = A;
          continue;
        }

        if ((A->frgBgn < alignWobble) &&
            (A->frgEnd > longest[frgIID].rptBgn))
          //  First N bases of the read is part of a repeat
          longest[frgIID].rptBgn = A->frgEnd;

        if ((A->frgEnd > longest[frgIID].frgLen - alignWobble) &&
            (A->frgBgn < longest[frgIID].rptEnd))
          //  Last N bases of the read is part of a repeat
          longest[frgIID].rptEnd = A->frgBgn;
      }

      assert(LONG != 0L);

      fprintf(outputFile, "TRIM longest %8d to %8d %8d\n",
              frgIID,
              longest[frgIID].rptBgn, longest[frgIID].rptEnd);

      if (longest[frgIID].rptBgn > 0)
        addAlignment(genome,
                     LONG->frgIID,
                     0, longest[frgIID].rptBgn, false,
                     LONG->genBgn,
                     LONG->genBgn + longest[frgIID].rptBgn);

      if (longest[frgIID].rptEnd < LONG->frgEnd)
        addAlignment(genome,
                     LONG->frgIID,
                     longest[frgIID].rptEnd, LONG->frgEnd, false,
                     LONG->genBgn + longest[frgIID].rptEnd,
                     LONG->genBgn + LONG->frgEnd);

    }

    bgn = end;
  }  //  Over all matches, by fragment

  fclose(outputFile);


  sort(genome.begin(), genome.end(), byGenomePosition);


  //  Find the spanned repeats.  For eacn non-repeat fragment, search forward for any repeats that
  //  it spans.  Actually, we just do any fragment that it spans, so a fragment could be non-repeat
  //  and spanned.
  //
  //  We assume that if a read is spanned by a repeat, that read is also a repeat.
  //
  for (uint32 uniq=0; uniq<genome.size(); uniq++) {
    if (genome[uniq].isRepeat == true)
      continue;
    if (genome[uniq].isSpanned == true)
      continue;

    for (uint32 test = uniq + 1; ((test < genome.size()) &&
                                  (genome[test].genBgn <= genome[uniq].genEnd)); test++) {
      assert(genome[uniq].genBgn <= genome[test].genBgn);
      assert(genome[uniq].genBgn <= genome[test].genEnd);

      if ((genome[uniq].genBgn + uniqEnd <= genome[test].genBgn) &&
          (genome[test].genEnd           <= genome[uniq].genEnd - uniqEnd))
        genome[test].isSpanned = true;
    }
  }


  //  Now, throw all the non-spanned repeats into an intervalList, squash them to get intervals,
  //  and report the repeat regions.

  sprintf(outputName, "%s.intervalLog", outputPrefix);
  outputFile = fopen(outputName, "w");

  intervalList  REPT;
  intervalList  UNIQ;

  for (uint32 i=0; i<genome.size(); i++) {
    int32  frgIID = genome[i].frgIID;
    int32  bgnOff = longest[frgIID].rptBgn;
    int32  endOff = longest[frgIID].frgLen - longest[frgIID].rptEnd;
    int32  frgOff = (genome[i].isReverse) ? endOff : bgnOff;
    int32  len = genome[i].genEnd - genome[i].genBgn;

    if (genome[i].isSpanned == true) {
      if (outputFile)
        fprintf(outputFile, "SPAN frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
                frgIID,
                genome[i].frgBgn, genome[i].frgEnd,
                genome[i].isReverse ? 'r' : 'f',
                genome[i].genBgn, genome[i].genEnd);
      continue;
    }

    if (genome[i].isRepeat == true) {
      REPT.add(genome[i].genBgn, len);
      if (outputFile)
        fprintf(outputFile, "REPT frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
                frgIID,
                genome[i].frgBgn, genome[i].frgEnd,
                genome[i].isReverse ? 'r' : 'f',
                genome[i].genBgn, genome[i].genEnd);
      continue;
    }

    //  Allow the fragment to extend into the repeat region.  We've added that region
    //  as a separate repeat alignment.  If something else spans it, the alignment will be removed,
    //  and we'll get a good overlap.
    bgnOff = 0;
    endOff = 0;

    len -= bgnOff;
    len -= endOff;
    len -= 2 * minOverlap;

    if (len <= 0)
      continue;

    UNIQ.add(genome[i].genBgn + frgOff + minOverlap,
             len);

    if (outputFile)
      fprintf(outputFile, "UNIQ frgIID %8d frg %8d,%8d rev %c gen %8d,%8d mod %8d,%8d\n",
              frgIID,
              genome[i].frgBgn, genome[i].frgEnd,
              genome[i].isReverse ? 'r' : 'f',
              genome[i].genBgn, genome[i].genEnd,
              genome[i].genBgn + frgOff + minOverlap,
              genome[i].genBgn + frgOff + minOverlap + len);
  }

  fclose(outputFile);

  REPT.merge();
  UNIQ.merge();

  bool         *REPTvalid = new bool [REPT.numberOfIntervals()];
  bool         *UNIQvalid = new bool [UNIQ.numberOfIntervals()];

  for (uint32 i=0; i<REPT.numberOfIntervals(); i++)
    REPTvalid[i] = true;
  for (uint32 i=0; i<UNIQ.numberOfIntervals(); i++)
    UNIQvalid[i] = true;

  //  Search for exceptions -- one completely contained in the other

  int32   UNIQexceptions=0;
  int32   REPTexceptions=0;

  for (uint32 ir=0; ir<REPT.numberOfIntervals(); ir++) {
    for (uint32 iu=0; iu<UNIQ.numberOfIntervals(); iu++) {
      if ((REPT.lo(ir) <= UNIQ.lo(iu)) &&
          (UNIQ.hi(iu) <= REPT.hi(ir))) {
        fprintf(stderr, "EXCEPTION:  UNIQ %ld,%ld len=%ld ct=%d contained in REPT %ld,%ld len=%ld ct=%d\n",
                UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.ct(iu),
                REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.ct(ir));
        UNIQvalid[iu] = false;
        UNIQexceptions++;
      }
      if ((UNIQ.lo(iu) <= REPT.lo(ir)) &&
          (REPT.hi(ir) <= UNIQ.hi(iu))) {
        fprintf(stderr, "EXCEPTION:  REPT %ld,%ld len=%ld  ct=%d contained in UNIQ %ld,%ld len=%ld ct=%d\n",
                REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.ct(ir),
                UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.ct(iu));
        REPTvalid[ir] = false;
        REPTexceptions++;
      }
    }
  }

  fprintf(stderr, "Found %d REPT intervals, and %d REPT exceptions.\n", REPT.numberOfIntervals() - REPTexceptions, REPTexceptions);
  fprintf(stderr, "Found %d UNIQ intervals, and %d UNIQ exceptions.\n", UNIQ.numberOfIntervals() - UNIQexceptions, UNIQexceptions);

  //  DEBUG (I guess) report the intervals.

  sprintf(outputName, "%s.intervals", outputPrefix);
  outputFile = fopen(outputName, "w");

  for (uint32 ir=0, iu=0; ((ir < REPT.numberOfIntervals()) ||
                           (iu < UNIQ.numberOfIntervals())); ) {
    int64  lor = (ir < REPT.numberOfIntervals()) ? REPT.lo(ir) : 999999999;
    int64  lou = (iu < UNIQ.numberOfIntervals()) ? UNIQ.lo(iu) : 999999999;

    if (lor < lou) {
      fprintf(outputFile, "%8ld %8ld REPT %d%s\n",
              REPT.lo(ir),
              REPT.hi(ir),
              REPT.ct(ir),
              (REPTvalid[ir]) ? "" : " weak");
      ir++;
    } else {
      fprintf(outputFile, "%8ld %8ld UNIQ %d%s\n",
              UNIQ.lo(iu) - minOverlap,
              UNIQ.hi(iu) + minOverlap,
              UNIQ.ct(iu),
              (UNIQvalid[iu]) ? "" : " separation");
      iu++;
    }
  }

  fclose(outputFile);

  //  Load the reference sequence.  ASSUMES it is all on one line!

  FILE    *F      = fopen(refName, "r");
  char    *refhdr = new char [1024];
  char    *refseq = new char [16 * 1024 * 1024];
  char    *refano = new char [16 * 1024 * 1024];
  int32    reflen = 0;

  fgets(refhdr,             1024, F);
  fgets(refseq, 16 * 1026 * 1024, F);

  fclose(F);

  reflen = strlen(refseq);

  //  Label the reference with REPEAT or UNIQUE, output fasta.  (yes, this is pointless)

  for (int32 p=0; p<reflen; p++)
    refano[p] = 'G';

  for (uint32 i=0; i<UNIQ.numberOfIntervals(); i++) {
    assert(UNIQ.lo(i) >= 0);
    assert(UNIQ.hi(i) <= reflen);
    assert(UNIQ.lo(i) < UNIQ.hi(i));

    if (UNIQvalid[i] == false)
      continue;

    for (int32 p=UNIQ.lo(i); p<UNIQ.hi(i); p++)
      refano[p] = 'U';
  }

  for (uint32 i=0; i<REPT.numberOfIntervals(); i++) {
    assert(REPT.lo(i) >= 0);
    assert(REPT.hi(i) <= reflen);
    assert(REPT.lo(i) < REPT.hi(i));

    if (REPTvalid[i] == false)
      continue;

    for (int32 p=REPT.lo(i); p<REPT.hi(i); p++)
      refano[p] = 'R';
  }

  //  Output a multifasta of ideal unitigs.

  sprintf(outputName, "%s.fasta", outputPrefix);
  outputFile = fopen(outputName, "w");

  for (uint32 ir=0, iu=0; ((ir < REPT.numberOfIntervals()) ||
                           (iu < UNIQ.numberOfIntervals())); ) {
    int64  lor = (ir < REPT.numberOfIntervals()) ? REPT.lo(ir) : 999999999;
    int64  lou = (iu < UNIQ.numberOfIntervals()) ? UNIQ.lo(iu) : 999999999;

    if (lor < lou) {
      if (REPTvalid[ir] == true) {
        char save = refseq[REPT.hi(ir)];

        refseq[REPT.hi(ir)] = 0;
        fprintf(outputFile, ">utg%08dr\n%s\n",
                ir + iu,
                refseq + REPT.lo(ir));
        refseq[REPT.hi(ir)] = save;
      }

      ir++;
    }

    if (lou < lor) {
      if (UNIQvalid[iu] == true) {
        char save = refseq[UNIQ.hi(iu)];

        refseq[UNIQ.hi(iu)] = 0;
        fprintf(outputFile, ">utg%08du\n%s\n",
                ir + iu,
                refseq + UNIQ.lo(iu));
        refseq[UNIQ.hi(iu)] = save;
      }

      iu++;
    }
  }

  fclose(outputFile);
}
