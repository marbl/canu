
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/bogus.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2014-DEC-23
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "bogusUtil.H"

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

intervalList<int32>  REPT;
intervalList<int32>  UNIQ;

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
    fprintf(gffOutput, "%s\t.\tbogus_raw_input\t%d\t%d\t.\t%c\t.\tID=align%08d-frag%08d;Name=%s;Note=%d-%d\n",
            refList[genome[i].genIID].rsrefName,        //  reference name
            genome[i].genBgn,                           //  begin position in base-based
            genome[i].genEnd,                           //  end position
            (genome[i].isReverse) ? '-' : '+',          //  strand
            i,                                          //  ID - match id
            genome[i].frgIID,                           //  ID - frag id
            IIDname[genome[i].frgIID].c_str(),          //  Name - actual sequence name
            genome[i].frgBgn, genome[i].frgEnd);        //  Note - position on frag
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

#ifdef DEBUG
    for (uint32 ii=bgn; ii<end; ii++)
      fprintf(stderr, "processMatches()-- INIT %u-%u frg %u at chn %u-%u gen %u %u-%u\n",
              bgn, end, genome[ii].frgIID, genome[ii].chnBgn, genome[ii].chnEnd, genome[ii].genIID, genome[ii].genBgn, genome[ii].genEnd);
#endif

    //  Find the longest alignment.  Set 'num' to zero, we'll fill it out correctly in a minute.
    //  Also, based on the longest alignment, compute an indel rate.

    double indelRate = 1.0;

    for (uint32 i=bgn; i<end; i++) {
      genomeAlignment  *A = &genome[i];

      if (longest[frgIID].len < A->frgEnd - A->frgBgn) {
        longest[frgIID].bgn = A->frgBgn;
        longest[frgIID].end = A->frgEnd;
        longest[frgIID].len = A->frgEnd - A->frgBgn;
        longest[frgIID].num = 0;

        indelRate = (A->genEnd - A->genBgn) / (double)(A->frgEnd - A->frgBgn);
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

      //  Indel screws up this compute.  We have the repeat region marked on the read, and need to transfer it to the genome.
      //  We estimate, globally, the indel rate, and scale the repeat region length.

      if (longest[frgIID].rptBgn > 0) {
        uint32 rptLen = (longest[frgIID].rptBgn - 0) * indelRate;

        if (rptLen > 0) {
#ifdef DEBUG
          fprintf(stderr, "addAlignment()-- RPTBGN longest %u-%u frg %u at frg %u-%u chn %u-%u gen %u %u-%u rptlen %u indelrate %f\n",
                  longest[frgIID].rptBgn, longest[frgIID].rptEnd, LONG->frgIID, LONG->frgBgn, LONG->frgEnd, LONG->chnBgn, LONG->chnEnd, LONG->genIID, LONG->genBgn, LONG->genEnd, rptLen, indelRate);
#endif
          addAlignment(genome,
                       LONG->frgIID,
                       0, longest[frgIID].rptBgn, false,
                       LONG->chnBgn,
                       LONG->chnBgn + rptLen,
                       LONG->identity,
                       LONG->genIID,
                       LONG->genBgn,
                       LONG->genBgn + rptLen);
#ifdef DEBUG
          uint32 ii = genome.size() - 1;
          fprintf(stderr, "addAlignment()-- RPTBGN FINI frg %u at frg %u-%u chn %u-%u gen %u %u-%u\n",
                  genome[ii].frgIID, genome[ii].frgBgn, genome[ii].frgEnd, genome[ii].chnBgn, genome[ii].chnEnd, genome[ii].genIID, genome[ii].genBgn, genome[ii].genEnd);
#endif
        }
      }

      if (longest[frgIID].rptEnd < LONG->frgEnd) {
        uint32 rptLen = (LONG->frgEnd - longest[frgIID].rptEnd) * indelRate;

        if (rptLen > 0) {
#ifdef DEBUG
          fprintf(stderr, "addAlignment()-- RPTEND longest %u-%u frg %u at frg %u-%u chn %u-%u gen %u %u-%u rptlen %u indelrate %f\n",
                  longest[frgIID].rptBgn, longest[frgIID].rptEnd, LONG->frgIID, LONG->frgBgn, LONG->frgEnd, LONG->chnBgn, LONG->chnEnd, LONG->genIID, LONG->genBgn, LONG->genEnd, rptLen, indelRate);
#endif
          addAlignment(genome,
                       LONG->frgIID,
                       longest[frgIID].rptEnd, LONG->frgEnd, false,
                       LONG->chnEnd - rptLen,
                       LONG->chnEnd,
                       LONG->identity,
                       LONG->genIID,
                       LONG->genEnd - rptLen,
                       LONG->genEnd);
#ifdef DEBUG
          uint32 ii = genome.size() - 1;
          fprintf(stderr, "addAlignment()-- RPTEND FINI frg %u at frg %u-%u chn %u-%u gen %u %u-%u\n",
                  genome[ii].frgIID, genome[ii].frgBgn, genome[ii].frgEnd, genome[ii].chnBgn, genome[ii].chnEnd, genome[ii].genIID, genome[ii].genBgn, genome[ii].genEnd);
#endif
        }
      }
    }

    bgn = end;
  }  //  Over all matches, by fragment

#ifdef DEBUG
  for (uint32 ii=0; ii<genome.size(); ii++)
    fprintf(stderr, "processMatches()-- FINI frg %u at chn %u-%u gen %u %u-%u\n",
            genome[ii].frgIID, genome[ii].chnBgn, genome[ii].chnEnd, genome[ii].genIID, genome[ii].genBgn, genome[ii].genEnd);
#endif
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
#ifdef DEBUG
    fprintf(stderr, "findSpannedMatched()-- frg %u at chn %u-%u gen %u %u-%u\n",
            genome[uniq].frgIID, genome[uniq].chnBgn, genome[uniq].chnEnd, genome[uniq].genIID, genome[uniq].genBgn, genome[uniq].genEnd);
#endif

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


#undef REPORT_INTERVALS_IN_GFF

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
#ifdef REPORT_INTERVALS_IN_GFF
      fprintf(gffOutput, "\n");
      fprintf(gffOutput, "#SPAN frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
              frgIID,
              genome[i].frgBgn, genome[i].frgEnd,
              genome[i].isReverse ? 'r' : 'f',
              genome[i].genBgn, genome[i].genEnd);
#endif
      if (includeRaw)
        fprintf(gffOutput, "%s\t.\tbogus_span_input\t%d\t%d\t.\t%c\t.\tID=SPAN%08d-frag%08d;Name=%s;Note=%d-%d\n",
                refhdr,                                //  reference name
                genome[i].genBgn, genome[i].genEnd,    //  reference position
                (genome[i].isReverse) ? '-' : '+',     //  strand
                i,                                     //  ID - match id
                genome[i].frgIID,                      //  ID - frag id
                IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
                genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
      continue;
    }

    if (genome[i].isRepeat == true) {
      REPT.add(genome[i].chnBgn, len);
#ifdef REPORT_INTERVALS_IN_GFF
      fprintf(gffOutput, "\n");
      fprintf(gffOutput, "#REPT frgIID %8d frg %8d,%8d rev %c gen %8d,%8d\n",
              frgIID,
              genome[i].frgBgn, genome[i].frgEnd,
              genome[i].isReverse ? 'r' : 'f',
              genome[i].genBgn, genome[i].genEnd);
#endif
      if (includeRaw)
        fprintf(gffOutput, "%s\t.\tbogus_rept_input\t%d\t%d\t.\t%c\t.\tID=REPT%08d-frag%08d;Name=%s;Note=%d-%d\n",
                refhdr,                                //  reference name
                genome[i].genBgn, genome[i].genEnd,    //  reference position
                (genome[i].isReverse) ? '-' : '+',     //  strand
                i,                                     //  ID - match id
                genome[i].frgIID,                      //  ID - frag id
                IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
                genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
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

#ifdef REPORT_INTERVALS_IN_GFF
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
      fprintf(gffOutput, "%s\t.\tbogus_uniq_input\t%d\t%d\t.\t%c\t.\tID=UNIQ%08d-frag%08d;Name=%s;Note=%d-%d\n",
              refhdr,                                //  reference name
              genome[i].genBgn, genome[i].genEnd,    //  reference position
              (genome[i].isReverse) ? '-' : '+',     //  strand
              i,                                     //  ID - match id
              genome[i].frgIID,                      //  ID - frag id
              IIDname[genome[i].frgIID].c_str(),     //  Name - actual sequence name
              genome[i].frgBgn, genome[i].frgEnd);   //  Note - position on frag
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
        //        UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.count(iu),
        //        REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.count(ir));
        UNIQvalid[iu] = false;
        UNIQvalidParent[iu] = ir;
        UNIQexceptions++;
      }
      if ((UNIQ.lo(iu) <= REPT.lo(ir)) &&
          (REPT.hi(ir) <= UNIQ.hi(iu))) {
        //fprintf(stderr, "EXCEPTION:  REPT %ld,%ld len=%ld  ct=%d contained in UNIQ %ld,%ld len=%ld ct=%d\n",
        //        REPT.lo(ir), REPT.hi(ir), REPT.hi(ir) - REPT.lo(ir), REPT.count(ir),
        //        UNIQ.lo(iu), UNIQ.hi(iu), UNIQ.hi(iu) - UNIQ.lo(iu), UNIQ.count(iu));
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

  //  Output intervals must be at least minLength bases long, and be formed from at least minFrags
  //  mappings.
  uint32   minLength  = 0;
  uint32   minFrags   = 0;

  //  Input matches must be at least minIdentity.
  double   minIdentity = 0;

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

    } else if (strcmp(argv[arg], "-minlength") == 0) {
      minLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-minfrags") == 0) {
      minFrags = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-minidentity") == 0) {
      minIdentity = atof(argv[++arg]);

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
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
    loadNucmer(nucmerNames[nn], genome, IIDmap, IIDname, refList, refMap, minIdentity);

  for (uint32 nn=0; nn<snapperNamesLen; nn++)
    loadSnapper(snapperNames[nn], genome, IIDmap, IIDname, refList, refMap, minIdentity);

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

  //  Extend UNIQ to cover gaps in alignments.
  //  Extend REPT to cover gaps in alignments.
  //  Merge adjacent/overlapping UNIQ-UNIQ or REPT-REPT intervals.



  //
  //  Write the output
  //

  for (uint32 ir=0, iu=0; ((ir < REPT.numberOfIntervals()) ||
                           (iu < UNIQ.numberOfIntervals())); ) {

    //  UNIQ regions are offset by fragTrim on each side

    int64  lor = (ir < REPT.numberOfIntervals()) ? REPT.lo(ir)            : 999999999;
    int64  hir = (ir < REPT.numberOfIntervals()) ? REPT.hi(ir)            : 999999999;
    int64  lou = (iu < UNIQ.numberOfIntervals()) ? UNIQ.lo(iu) - fragTrim : 999999999;
    int64  hiu = (iu < UNIQ.numberOfIntervals()) ? UNIQ.hi(iu) + fragTrim : 999999999;

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
          refcnt = REPT.count(ir);
        }

      if (refcnt == 0) {
        fprintf(stderr, "DIDN'T FIND REGION.\n");
        for (uint32 rr=0; rr<refList.size(); rr++) {
          fprintf(stderr, "  %3u rschnBgn %6d REPT %5d %ld-%ld UNIQ %5d %ld-%ld rschnEnd %6d REPT\n",
                  rr,
                  refList[rr].rschnBgn,
                  ir, lor, hir,
                  iu, lou, hiu,
                  refList[rr].rschnEnd);
        }
      }
      //assert(refcnt != 0);

      if ((refcnt > 0) && (minFrags <= refcnt) && (minLength <= refend - refbgn)) {
        fprintf(intervalOutput, "%s\t%8"F_S64P"\t%8"F_S64P"\tREPT\t"F_S64"%s\n",
                refhdr, refbgn, refend, refcnt, (REPTvalid[ir]) ? "" : " weak");

        if (REPTvalid[ir])
          fprintf(gffOutput, "%s\t.\tbogus_rept_interval\t"F_S64"\t"F_S64"\t.\t.\t.\tID=REPT%04d;fragCount="F_S64"\n",
                  refhdr, refbgn, refend, ir, refcnt);
        else
          fprintf(gffOutput, "%s\t.\tbogus_weak_interval\t"F_S64"\t"F_S64"\t.\t.\t.\tParent=UNIQ%04d;fragCount="F_S64"\n",
                  refhdr, refbgn, refend, REPTvalidParent[ir], refcnt);
      }

      ir++;

    } else {
      for (uint32 rr=0; rr<refList.size(); rr++)
        if ((refList[rr].rschnBgn <= lou) && (lou <= refList[rr].rschnEnd)) {
          refhdr = refList[rr].rsrefName;
          refbgn = UNIQ.lo(iu) - fragTrim - refList[rr].rschnBgn;
          refend = UNIQ.hi(iu) + fragTrim - refList[rr].rschnBgn;
          refcnt = UNIQ.count(iu);
        }

      //  Not sure why some data sets (long pacbio for example) trigger this.
#if 0
      if (refcnt == 0) {
        fprintf(stderr, "DIDN'T FIND REGION.\n");
        for (uint32 rr=0; rr<refList.size(); rr++) {
          fprintf(stderr, "  %3u rschnBgn %6d REPT %5d %ld-%ld UNIQ %5d %ld-%ld rschnEnd %6d UNIQ\n",
                  rr,
                  refList[rr].rschnBgn,
                  ir, lor, hir,
                  iu, lou, hiu,
                  refList[rr].rschnEnd);
        }
      }
#endif
      //assert(refcnt != 0);

      if ((refcnt > 0) && (minFrags <= refcnt) && (minLength <= refend - refbgn)) {
        fprintf(intervalOutput, "%s\t%8"F_S64P"\t%8"F_S64P"\tUNIQ\t"F_S64"%s\n",
                refhdr, refbgn, refend, refcnt, (UNIQvalid[iu]) ? "" : " separation");

        if (UNIQvalid[iu])
          fprintf(gffOutput, "%s\t.\tbogus_uniq_interval\t"F_S64"\t"F_S64"\t.\t.\t.\tID=UNIQ%04d;fragCount="F_S64"\n",
                  refhdr, refbgn, refend, iu, refcnt);
        else
          fprintf(gffOutput, "%s\t.\tbogus_sepr_interval\t"F_S64"\t"F_S64"\t.\t.\t.\tParent=REPT%04d;fragCount="F_S64"\n",
                  refhdr, refbgn, refend, UNIQvalidParent[iu], refcnt);
      }

      iu++;
    }
  }

  fclose(gffOutput);
  fclose(intervalOutput);

  //  See CVS version 1.3 for writing rept/uniq fasta

  return(0);
}
