
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

const char *mainid = "$Id: bogus.C,v 1.6 2011-02-11 03:25:37 brianwalenz Exp $";

#include "AS_BAT_bogusUtil.H"

class longestAlignment {
public:
  longestAlignment() {
    bgn = end = len = num = 0;
    rptBgn = rptEnd = frgLen = 0;
  };

  int32   bgn;  //  Begin coord on the fragment for this region
  int32   end;  //  End coord on the fragment for this region
  int32   len;  //  Length of the match
  int32   num;  //  Number of matches on this region

  int32   rptBgn;
  int32   rptEnd;

  int32   frgLen;
};


//  MUCH, much easier to modularize the code if all this stuff is globally available.

vector<referenceSequence>  refList;
map<string,uint32>         refMap;

vector<longestAlignment>   longest;
vector<genomeAlignment>    genome;
map<string, int32>         IIDmap;       //  Maps an ID string to an IID.
vector<string>             IIDname;      //  Maps an IID to an ID string.
vector<uint32>             IIDcount;     //  Maps an IID to the number of alignments

intervalList  REPT;
intervalList  UNIQ;

bool         *REPTvalid = NULL;
bool         *UNIQvalid = NULL;

int32        *REPTvalidParent = NULL;
int32        *UNIQvalidParent = NULL;

FILE         *gffOutput      = NULL;
FILE         *intervalOutput = NULL;


static
void
writeInputsAsGFF3(char *outputPrefix) {

  for (uint32 i=0; i<genome.size(); i++) {
    fprintf(gffOutput, "%s\t.\tbogus_raw_input\t%d\t%d\t.\t%c\t.\tID=align%08d-frag%08d\n",  //  ;Name=%s;Note=%d-%d
            refList[genome[i].genIID].rsrefName,        //  reference name
            genome[i].genBgn,                           //  begin position in base-based
            genome[i].genEnd,                           //  end position
            (genome[i].isReverse) ? '-' : '+',          //  strand
            i,                                          //  ID - match id
            genome[i].frgIID);                          //  ID - frag id
  }
}



static
void
processMatches(int32 alignWobble, int32 uniqEnd) {

  longest.resize(IIDname.size());

  for (uint32 bgn=0, lim=genome.size(); bgn<lim; ) {
    uint32 frgIID = genome[bgn].frgIID;

    uint32 end = bgn + 1;
    while ((end < genome.size()) &&
           (genome[bgn].frgIID == genome[end].frgIID))
      end++;

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

      if (longest[frgIID].rptBgn > 0)
        addAlignment(genome,
                     LONG->frgIID,
                     0, longest[frgIID].rptBgn, false,
                     LONG->chnBgn,
                     LONG->chnBgn + longest[frgIID].rptBgn,
                     LONG->genIID,
                     LONG->genBgn,
                     LONG->genBgn + longest[frgIID].rptBgn);

      if (longest[frgIID].rptEnd < LONG->frgEnd)
        addAlignment(genome,
                     LONG->frgIID,
                     longest[frgIID].rptEnd, LONG->frgEnd, false,
                     LONG->chnBgn + longest[frgIID].rptEnd,
                     LONG->chnBgn + LONG->frgEnd,
                     LONG->genIID,
                     LONG->genBgn + longest[frgIID].rptEnd,
                     LONG->genBgn + LONG->frgEnd);

    }

    bgn = end;
  }  //  Over all matches, by fragment
}



//  Find the spanned repeats.  For eacn non-repeat fragment, search forward for any repeats that
//  it spans.  Actually, we just do any fragment that it spans, so a fragment could be non-repeat
//  and spanned.
//
//  We assume that if a read is spanned by a repeat, that read is also a repeat.
//
static
void
findSpannedMatches(int32 uniqEnd) {

  for (uint32 uniq=0; uniq<genome.size(); uniq++) {
    if (genome[uniq].isRepeat == true)
      continue;
    if (genome[uniq].isSpanned == true)
      continue;

    for (uint32 test = uniq + 1; ((test < genome.size()) &&
                                  (genome[test].chnBgn <= genome[uniq].chnEnd)); test++) {
      assert(genome[uniq].chnBgn <= genome[test].chnBgn);
      assert(genome[uniq].chnBgn <= genome[test].chnEnd);

      if ((genome[uniq].chnBgn + uniqEnd <= genome[test].chnBgn) &&
          (genome[test].chnEnd           <= genome[uniq].chnEnd - uniqEnd))
        genome[test].isSpanned = true;
    }
  }
}



static
void
buildIntervals(char *outputPrefix, int32 fragTrim, bool includeRaw) {

  for (uint32 i=0; i<genome.size(); i++) {
    int32  frgIID = genome[i].frgIID;
    int32  bgnOff = longest[frgIID].rptBgn;
    int32  endOff = longest[frgIID].frgLen - longest[frgIID].rptEnd;
    int32  frgOff = (genome[i].isReverse) ? endOff : bgnOff;
    int32  len    = genome[i].chnEnd - genome[i].chnBgn;
    char  *refhdr = refList[genome[i].genIID].rsrefName;

    if (genome[i].isSpanned == true) {
#if 0
      fprintf(gffOutput, "\n");
      fprintf(gffOutput, "#SPAN frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
              frgIID,
              genome[i].frgBgn, genome[i].frgEnd,
              genome[i].isReverse ? 'r' : 'f',
              genome[i].genBgn, genome[i].genEnd);
#endif
      if (includeRaw)
        fprintf(gffOutput, "%s\t.\tbogus_span_input\t%d\t%d\t.\t%c\t.\tID=SPAN%08d-frag%08d\n", //;Name=%s;Note=%d-%d\n",
                refhdr,                                //  reference name
                genome[i].genBgn, genome[i].genEnd,    //  reference position
                (genome[i].isReverse) ? '-' : '+',     //  strand
                i,                                     //  ID - match id
                genome[i].frgIID);                     //  ID - frag id
      //IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
      //genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
      continue;
    }

    if (genome[i].isRepeat == true) {
      REPT.add(genome[i].chnBgn, len);
#if 0
      fprintf(gffOutput, "\n");
      fprintf(gffOutput, "#REPT frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
              frgIID,
              genome[i].frgBgn, genome[i].frgEnd,
              genome[i].isReverse ? 'r' : 'f',
              genome[i].genBgn, genome[i].genEnd);
#endif
      if (includeRaw)
        fprintf(gffOutput, "%s\t.\tbogus_rept_input\t%d\t%d\t.\t%c\t.\tID=REPT%08d-frag%08d\n", //;Name=%s;Note=%d-%d\n",
                refhdr,                                //  reference name
                genome[i].genBgn, genome[i].genEnd,    //  reference position
                (genome[i].isReverse) ? '-' : '+',     //  strand
                i,                                     //  ID - match id
                genome[i].frgIID);                     //  ID - frag id
      //IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
      //genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
      continue;
    }

    //  Allow the fragment to extend into the repeat region.  We've added that region
    //  as a separate repeat alignment.  If something else spans it, the alignment will be removed,
    //  and we'll get a good overlap.
    bgnOff = 0;
    endOff = 0;

    len -= bgnOff;
    len -= endOff;
    len -= fragTrim * 2;

    if (len <= 0)
      continue;

    UNIQ.add(genome[i].chnBgn + frgOff + fragTrim,
             len);

#if 0
    fprintf(gffOutput, "\n");
    fprintf(gffOutput, "#UNIQ frgIID %8d frg %8d,%8d rev %c gen %8d,%8d mod %8d,%8d\n",
            frgIID,
            genome[i].frgBgn, genome[i].frgEnd,
            genome[i].isReverse ? 'r' : 'f',
            genome[i].genBgn, genome[i].genEnd,
            genome[i].genBgn + frgOff + fragTrim,
            genome[i].genBgn + frgOff + fragTrim + len);
#endif
    if (includeRaw)
      fprintf(gffOutput, "%s\t.\tbogus_uniq_input\t%d\t%d\t.\t%c\t.\tID=UNIQ%08d-frag%08d\n", //;Name=%s;Note=%d-%d\n",
              refhdr,                                //  reference name
              genome[i].genBgn, genome[i].genEnd,    //  reference position
              (genome[i].isReverse) ? '-' : '+',     //  strand
              i,                                     //  ID - match id
              genome[i].frgIID);                     //  ID - frag id
    //IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
    //genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
  }
}



void
markWeak(void) {
  REPTvalid = new bool [REPT.numberOfIntervals()];
  UNIQvalid = new bool [UNIQ.numberOfIntervals()];

  REPTvalidParent = new int32 [REPT.numberOfIntervals()];
  UNIQvalidParent = new int32 [UNIQ.numberOfIntervals()];

  for (uint32 i=0; i<REPT.numberOfIntervals(); i++) {
    REPTvalid[i] = true;
    REPTvalidParent[i] = 0;
  }
  for (uint32 i=0; i<UNIQ.numberOfIntervals(); i++) {
    UNIQvalid[i] = true;
    UNIQvalidParent[i] = 0;
  }

  int32   UNIQexceptions=0;
  int32   REPTexceptions=0;

  for (uint32 ir=0; ir<REPT.numberOfIntervals(); ir++) {
    for (uint32 iu=0; iu<UNIQ.numberOfIntervals(); iu++) {
      if ((REPT.lo(ir) <= UNIQ.lo(iu)) &&
          (UNIQ.hi(iu) <= REPT.hi(ir))) {
        //fprintf(stderr, "EXCEPTION:  UNIQ %ld,%ld len=%ld ct=%d contained in REPT %ld,%ld len=%ld ct=%d\n",
        //        UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.ct(iu),
        //        REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.ct(ir));
        UNIQvalid[iu] = false;
        UNIQvalidParent[iu] = ir;
        UNIQexceptions++;
      }
      if ((UNIQ.lo(iu) <= REPT.lo(ir)) &&
          (REPT.hi(ir) <= UNIQ.hi(iu))) {
        //fprintf(stderr, "EXCEPTION:  REPT %ld,%ld len=%ld  ct=%d contained in UNIQ %ld,%ld len=%ld ct=%d\n",
        //        REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.ct(ir),
        //        UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.ct(iu));
        REPTvalid[ir] = false;
        REPTvalidParent[ir] = iu;
        REPTexceptions++;
      }
    }
  }

  fprintf(stderr, "Found %d REPT intervals, and %d REPT weak intervals.\n", REPT.numberOfIntervals() - REPTexceptions, REPTexceptions);
  fprintf(stderr, "Found %d UNIQ intervals, and %d UNIQ weak intervals.\n", UNIQ.numberOfIntervals() - UNIQexceptions, UNIQexceptions);
}




int
main(int argc, char **argv) {
  uint32   nucmerNamesLen  = 0;
  uint32   snapperNamesLen = 0;
  char    *nucmerNames[1024];
  char    *snapperNames[1024];
  char    *refName             = 0L;
  char    *outputPrefix        = 0L;

  //  When comparing coords to if two alignments span the same piece of the fragment,
  //  allow (twice) this much slop in the comparison.  E.g., Abgn +- 5 == Bbgn
  //
  int32  alignWobble  = 5;

  //  When constructing the unique coverage map, trim each read by this amount, on each end, to
  //  simulate the minimum overlap length needed by unitigger.  This amount is automagically added
  //  back in on output.
  int32  fragTrim     = 40 / 2;

  //  When testing if a repeat alignment is contained in a non-repeat alignment, the repeat must
  //  start at least this many bases from the end of the non-repeat alignment.  In other words, a
  //  fragment with a repeat in the middle can be uniquely placed (by overlaps) with only 20 bases of
  //  unique sequence on the end.
  int32  uniqEnd      = 40;

  //  Loading jbrowse with all the raw input reads on large(r) genomes either
  //  fails, or takes forever.  This disables the raw read output.
  bool   includeRaw = true;


  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-nucmer") == 0) {
      nucmerNames[nucmerNamesLen++] = argv[++arg];

    } else if (strcmp(argv[arg], "-snapper") == 0) {
      snapperNames[snapperNamesLen++] = argv[++arg];

    } else if (strcmp(argv[arg], "-reference") == 0) {
      refName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-wobble") == 0) {
      alignWobble = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-overlap") == 0) {
      fragTrim = atoi(argv[++arg]) / 2;

    } else if (strcmp(argv[arg], "-noraw") == 0) {
      includeRaw = false;

    } else {
      err++;
    }
    arg++;
  }
  if ((nucmerNamesLen == 0) && (snapperNamesLen == 0))
    fprintf(stderr, "ERROR: No input matches supplied (either -nucmer or -snapper).\n"), err++;
  if (refName == 0L)
    fprintf(stderr, "ERROR: No reference sequence supplied (-reference).\n"), err++;
  if (outputPrefix == 0L)
    fprintf(stderr, "ERROR: No output prefix supplied (-output).\n"), err++;
  if (err) {
    exit(1);
  }

  {
    char   outputName[FILENAME_MAX];

    errno = 0;

    sprintf(outputName, "%s.intervals", outputPrefix);
    intervalOutput = fopen(outputName, "w");

    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n",
              outputName, strerror(errno)), exit(1);

    sprintf(outputName, "%s.gff3", outputPrefix);
    gffOutput = fopen(outputName, "w");
  
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n",
              outputName, strerror(errno)), exit(1);

    fprintf(gffOutput, "##gff-version 3\n");
  }

  //  Load the reference sequence.  ASSUMES it is all on one line!

  loadReferenceSequence(refName, refList, refMap);

  //  Load all the matches into genomeAlignment.  Generate longestAlignment for each fragment.
  //  genomeAlignment::isLognest and genomeAlignment::isRepeat are computed later.

  for (uint32 nn=0; nn<nucmerNamesLen; nn++)
    loadNucmer(nucmerNames[nn], genome, IIDmap, IIDname, refList, refMap);

  for (uint32 nn=0; nn<snapperNamesLen; nn++)
    loadSnapper(snapperNames[nn], genome, IIDmap, IIDname, refList, refMap);

  //if (includeRaw)
  //  writeInputsAsGFF3(outputPrefix);

  //  Process the matches fragment by fragment.  Find the longest, count the number
  //  of duplicates of the longest match, label as repeat/unique.

  sort(genome.begin(), genome.end(), byFragmentID);

  processMatches(alignWobble, uniqEnd);

  sort(genome.begin(), genome.end(), byGenomePosition);

  findSpannedMatches(uniqEnd);

  //  Now, throw all the non-spanned repeats into an intervalList, squash them to get intervals,
  //  and report the repeat regions.

  buildIntervals(outputPrefix, fragTrim, includeRaw);

  REPT.merge();
  UNIQ.merge();


  //  Search for exceptions -- one completely contained in the other

  markWeak();

  //
  //  Write the output
  //

  for (uint32 ir=0, iu=0; ((ir < REPT.numberOfIntervals()) ||
                           (iu < UNIQ.numberOfIntervals())); ) {

    //  UNIQ regions are offset by fragTrim on each side

    int64  lor = (ir < REPT.numberOfIntervals()) ? REPT.lo(ir)            : 999999999;
    int64  lou = (iu < UNIQ.numberOfIntervals()) ? UNIQ.lo(iu) - fragTrim : 999999999;


    //  Search the refList for the reference sequence we are in.  We should never span reference
    //  sequences (which isn't tested, as we only know the low coordinate at this time).

    char  *refhdr = NULL;
    int64  refbgn = 0;
    int64  refend = 0;
    int64  refcnt = 0;

    if (lor < lou) {
      for (uint32 rr=0; rr<refList.size(); rr++)
        if ((refList[rr].rschnBgn <= lor) && (lor <= refList[rr].rschnEnd)) {
          refhdr = refList[rr].rsrefName;
          refbgn = REPT.lo(ir) - refList[rr].rschnBgn;
          refend = REPT.hi(ir) - refList[rr].rschnBgn;
          refcnt = REPT.ct(ir);
        }

      assert(refcnt != 0);

      fprintf(intervalOutput, "%s\t%8ld\t%8ld\tREPT\t%ld%s\n",
              refhdr, refbgn, refend, refcnt, (REPTvalid[ir]) ? "" : " weak");

      if (REPTvalid[ir])
        fprintf(gffOutput, "%s\t.\tbogus_rept_interval\t%ld\t%ld\t.\t.\t.\tID=REPT%04d\n",
                refhdr, refbgn, refend, ir);
      else
        fprintf(gffOutput, "%s\t.\tbogus_weak_interval\t%ld\t%ld\t.\t.\t.\tParent=UNIQ%04d\n",
                refhdr, refbgn, refend, REPTvalidParent[ir]);

      ir++;

    } else {
      for (uint32 rr=0; rr<refList.size(); rr++)
        if ((refList[rr].rschnBgn <= lou) && (lou <= refList[rr].rschnEnd)) {
          refhdr = refList[rr].rsrefName;
          refbgn = UNIQ.lo(iu) - fragTrim - refList[rr].rschnBgn;
          refend = UNIQ.hi(iu) + fragTrim - refList[rr].rschnBgn;
          refcnt = UNIQ.ct(iu);
        }

      assert(refcnt != 0);

      fprintf(intervalOutput, "%s\t%8ld\t%8ld\tUNIQ\t%ld%s\n",
              refhdr, refbgn, refend, refcnt, (UNIQvalid[iu]) ? "" : " separation");

      if (UNIQvalid[iu])
        fprintf(gffOutput, "%s\t.\tbogus_uniq_interval\t%ld\t%ld\t.\t.\t.\tID=UNIQ%04d\n",
                refhdr, refbgn, refend, iu);
      else
        fprintf(gffOutput, "%s\t.\tbogus_sepr_interval\t%ld\t%ld\t.\t.\t.\tParent=REPT%04d\n",
                refhdr, refbgn, refend, UNIQvalidParent[iu]);

      iu++;
    }
  }

  fclose(gffOutput);
  fclose(intervalOutput);

  //  See CVS version 1.3 for writing rept/uniq fasta

  return(0);
}
