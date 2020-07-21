
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "runtime.H"
#include "sqStore.H"
#include "tgStore.H"

#include "edlib.H"

#include "strings.H"
#include "sequence.H"

#include "gfa.H"
#include "bed.H"

#define IS_GFA   1
#define IS_BED   2

#define STEP_BP        500

class sequence {
public:
  sequence() {
    seq = NULL;
    len = 0;
  };
  ~sequence() {
    delete [] seq;
  };

  void  set(tgTig *tig) {
    len = tig->length();
    seq = new char [len + 1];

    if (tig->consensusExists())
       memcpy(seq, tig->bases(), len);
    else
       memset(seq, 0, len);

    seq[len] = 0;
  };

  char   *seq;
  uint32  len;
};



class sequences {
public:
  sequences(char *tigName, uint32 tigVers) {
    if (tigVers <= 0)
       fprintf(stderr, "Error: invalid version of tigstore: %d, not supported\n", tigVers);
    assert(tigVers > 0);

    tgStore *tigStore = new tgStore(tigName, tigVers);

    b    = 0;
    e    = tigStore->numTigs();
    seqs = new sequence [e+1];
    used = new uint32   [e+1];

    for (uint32 ti=b; ti < e; ti++) {
      tgTig *tig = tigStore->loadTig(ti);

      used[ti] = 0;

      if (tig == NULL)
        continue;

      seqs[ti].set(tig);

      tigStore->unloadTig(ti);
    }

    delete tigStore;
  };

  ~sequences() {
    delete [] seqs;
    delete [] used;
  };

  sequence &operator[](uint32 xx) {
    if (xx < e)
      return(seqs[xx]);

    fprintf(stderr, "ERROR: sequence id %u out of range b=%u e=%u\n", xx, b, e);

    assert(xx < e);
    return(seqs[0]);
  };

  uint32    b;
  uint32    e;
  sequence *seqs;
  uint32   *used;
};



void
dotplot(uint32 Aid, bool Afwd, char *Aseq,
        uint32 Bid, bool Bfwd, char *Bseq) {
  char   Aname[16];
  char   Bname[16];
  char   Pname[64];
  FILE  *F;

  sprintf(Aname, "tig%08u%c",  Aid, (Afwd) ? '+' : '-');
  sprintf(Bname, "tig%08u%c",  Bid, (Bfwd) ? '+' : '-');
  sprintf(Pname, "plot-%s-%s", Aname, Bname);

  F = AS_UTL_openOutputFile(Pname, '.', "sh");
  fprintf(F, "#!/bin/sh\n");
  fprintf(F, "\n");
  fprintf(F, "nucmer --maxmatch --nosimplify -p %s %s.fasta %s.fasta\n", Pname, Aname, Bname);
  fprintf(F, "show-coords -l -o -r -T %s.delta | expand -t 8 > %s.coords\n", Pname, Pname);
  fprintf(F, "mummerplot --fat -t png -p %s %s.delta\n", Pname, Pname);
  fprintf(F, "echo mummerplot --fat -p %s %s.delta\n", Pname, Pname);
  AS_UTL_closeFile(F, Pname, '.', "sh");

  F = AS_UTL_openOutputFile(Aname, '.', "fasta");
  fprintf(F, ">%s\n%s\n", Aname, Aseq);
  AS_UTL_closeFile(F, Aname, '.', "fasta");

  F = AS_UTL_openOutputFile(Bname, '.', "fasta");
  fprintf(F, ">%s\n%s\n", Bname, Bseq);
  AS_UTL_closeFile(F, Bname, '.', "fasta");

  sprintf(Pname, "sh plot-%s-%s.sh", Aname, Bname);
  system(Pname);
}



bool
checkLink(gfaLink   *link,
          sequences &seqs,
          sequences &seqs_orig,
          double     erate,
          bool       beVerbose,
          bool       doPlot) {

  char   *Aseq = seqs[link->_Aid].seq, *Arev = NULL;
  char   *Bseq = seqs[link->_Bid].seq, *Brev = NULL;

  int32  Abgn, Aend, Alen = seqs[link->_Aid].len;
  int32  Bbgn, Bend, Blen = seqs[link->_Bid].len;

  double ratio = max((double)Alen/seqs_orig[link->_Aid].len, (double)Blen/seqs_orig[link->_Bid].len);

  EdlibAlignResult  result  = { 0, NULL, NULL, 0, NULL, 0, 0 };

  int32  AalignLen = 0;
  int32  BalignLen = 0;
  int32  editDist  = 0;
  int32  alignLen  = 0;
  int32  maxEdit   = 0;

  //  NOTE!  edlibAlign calls the 'A' sequence the 'query' and
  //  the 'B' sequence the 'target' (aka, 'reference').

  link->alignmentLength(AalignLen, BalignLen, alignLen);

  //  Regardless of if we find an new alignment or not, remove the old one.
  //  If we don't find a new one, we'll discard the link.

  delete [] link->_cigar;
  link->_cigar = NULL;

  if (link->_Afwd == false)
    Aseq = Arev = reverseComplementCopy(Aseq, Alen);
  if (link->_Bfwd == false)
    Bseq = Brev = reverseComplementCopy(Bseq, Blen);

  //  Ty to find the end coordinate on B.  Align the last bits of A to B.
  //
  //   -------(---------]     v--??
  //             [------------)------
  //
  // if we fail to find the specified link, we will shrink A to see if we can find it that way and B will be trimmed to the right point
  // no need to retry on the second time because we'll have the right location of B and will trim A to it

  Abgn = max(Alen - AalignLen, 0);
  Aend =     Alen;

  Bbgn = 0;
  Bend = min(Blen, (int32)(1.10 * ratio * BalignLen));  //  Expand by whatever factor the contig grew during consensus plus 10%

  maxEdit = (int32)ceil(alignLen * erate);

  while (result.numLocations == 0) {
    if (beVerbose)
      fprintf(stderr, "LINK tig%08u %c %17s    tig%08u %c %17s   Aalign %6u Balign %6u align %6u\n",
              link->_Aid, (link->_Afwd) ? '+' : '-', "",
              link->_Bid, (link->_Bfwd) ? '+' : '-', "",
              AalignLen, BalignLen, alignLen);

    if (beVerbose)
      fprintf(stderr, "TEST tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  maxEdit=%6d  (extend B)",
              link->_Aid, (link->_Afwd) ? '+' : '-', Abgn, Aend,
              link->_Bid, (link->_Bfwd) ? '+' : '-', Bbgn, Bend,
              maxEdit);

    result = edlibAlign(Aseq + Abgn, Aend-Abgn,  //  The 'query'
                        Bseq + Bbgn, Bend-Bbgn,  //  The 'target'
                        edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

    if (result.numLocations > 0) {
      if (beVerbose)
        fprintf(stderr, "\n");
      Bend = Bbgn + result.endLocations[0] + 1;  // 0-based to space-based
      edlibFreeAlignResult(result);
      break; // found it, stop
    } else {
      if (beVerbose)
        fprintf(stderr, " - FAILED\n");
      if (Abgn + STEP_BP >= Alen) break; // we ran out of sequence, stop
      Abgn += STEP_BP;
    }
  }

  //  Do the same for A.  Aend and Bbgn never change; Bend was set above.
  //
  //   ------(--------------]
  //         ^--??  [-------]-----------
  //

  Abgn = max(Alen - (int32)(1.10 * ratio * AalignLen), 0);  //  Expand by whatever factor the contig grew during consensus plus 10%
  // update max edit by new length
  maxEdit = (int32)ceil((Bend-Bbgn+1) * erate);

  if (beVerbose)
    fprintf(stderr, "     tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  maxEdit=%6d  (extend A)",
            link->_Aid, (link->_Afwd) ? '+' : '-', Abgn, Aend,
            link->_Bid, (link->_Bfwd) ? '+' : '-', Bbgn, Bend,
            maxEdit);

  //  NEEDS to be MODE_HW because we need to find the suffix alignment.

  result = edlibAlign(Bseq + Bbgn, Bend-Bbgn,  //  The 'query'
                      Aseq + Abgn, Aend-Abgn,  //  The 'target'
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  if (result.numLocations > 0) {
    if (beVerbose)
      fprintf(stderr, "\n");
    Abgn = Abgn + result.startLocations[0];
    edlibFreeAlignResult(result);
  } else {
    if (beVerbose)
      fprintf(stderr, " - FAILED\n");
  }

  //  One more alignment, this time, with feeling - notice EDLIB_MODE_MW and EDLIB_TASK_PATH.

  if (beVerbose)
    fprintf(stderr, "     tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  maxEdit=%6d  (final)",
            link->_Aid, (link->_Afwd) ? '+' : '-', Abgn, Aend,
            link->_Bid, (link->_Bfwd) ? '+' : '-', Bbgn, Bend,
            maxEdit);

  result = edlibAlign(Aseq + Abgn, Aend-Abgn,
                      Bseq + Bbgn, Bend-Bbgn,
                      edlibNewAlignConfig(2 * maxEdit, EDLIB_MODE_NW, EDLIB_TASK_PATH));


  bool   success = false;

  if (result.numLocations > 0) {
    if (beVerbose)
      fprintf(stderr, "\n");

    editDist = result.editDistance;
    alignLen = ((Aend - Abgn) + (Bend - Bbgn) + (editDist)) / 2;
    alignLen = result.alignmentLength;     //  Edlib 'alignmentLength' is populated only for TASK_PATH

    link->_cigar = edlibAlignmentToCigar(result.alignment,
                                         result.alignmentLength, EDLIB_CIGAR_STANDARD);

    edlibFreeAlignResult(result);

    success = true;
  } else {
    if (beVerbose)
      fprintf(stderr, " - FAILED\n");
  }

  if (beVerbose)
    fprintf(stderr, "     tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d   %.4f\n",
            link->_Aid, (link->_Afwd) ? '+' : '-', Abgn, Aend,
            link->_Bid, (link->_Bfwd) ? '+' : '-', Bbgn, Bend,
            (double)editDist / alignLen);

  //  Make a plot.

  if ((success == false) && (doPlot == true))
    dotplot(link->_Aid, link->_Afwd, Aseq,
            link->_Bid, link->_Bfwd, Bseq);

  //  Cleanup for the next link.

  delete [] Arev;
  delete [] Brev;

  if (beVerbose)
    fprintf(stderr, "\n");

  return(success);
}



//   Align all of B into A.  Extend A as needed to make the whole thing fit.
//   Abgn, Aend and score are updated with the alignment.
//
bool
checkRecord_align(char const *label,
                  char const *Aname, char const *Aseq, int32 Alen, int32 &Abgn, int32 &Aend,
                  char const *Bname, char const *Bseq, int32 Blen,
                  int32 &score,
                  bool   beVerbose) {

  EdlibAlignResult  result  = { 0, NULL, NULL, 0, NULL, 0, 0 };

  int32  editDist    = 0;
  int32  alignLen    = 0;
  int32  alignScore  = 0;
  int32  maxEdit     = (int32)ceil(Blen * 0.03);  //  Should be the same sequence, but allow for a little difference.
  int32  step        = (int32)ceil(Blen * 0.15);

  Aend = min(Aend + 2 * step, Alen);        //  Limit Aend to the actual length of the contig (consensus can shrink)
  Abgn = max(Aend - Blen - 2 * step, 0);    //  Then push Abgn back to make space for the unitig.

 tryAgain:
  if (beVerbose)
    fprintf(stderr, "ALIGN %5s utg %s len=%7d to ctg %s %9d-%9d len=%9d",
            label,
            Bname, Blen,
            Aname, Abgn, Aend, Alen);

#if 0
  char N[FILENAME_MAX+1];
  FILE *F;

  char  ach = Aseq[Aend];   Aseq[Aend] = 0;
  char  bch = Bseq[Bend];   Bseq[Bend] = 0;

  sprintf(N, "compare%04d-%04d-ctg%04d.fasta", record->_Aid, record->_Bid, record->_Aid);
  F = AS_UTL_openOutputFile(N);
  fprintf(F, ">ctg%04d\n%s\n", record->_Aid, Aseq + Abgn);
  AS_UTL_closeFile(F, N);

  sprintf(N, "compare%04d-%04d-utg%04d.fasta", record->_Aid, record->_Bid, record->_Bid);
  F = AS_UTL_openOutputFile(N);
  fprintf(F, ">utg%04d\n%s\n", record->_Bid, Bseq + Bbgn);
  AS_UTL_closeFile(F, N);

  Aseq[Aend] = ach;
  Bseq[Bend] = bch;
#endif

  result = edlibAlign(Bseq,        Blen,       //  The 'query'   (unitig)
                      Aseq + Abgn, Aend-Abgn,  //  The 'target'  (contig)
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  //  Got an alignment?  Process and report, and maybe try again.

  if (result.numLocations > 0) {
    int32  nAbgn = Abgn + result.startLocations[0];
    int32  nAend = Abgn + result.endLocations[0] + 1;  // 0-based to space-based
    char  *cigar = NULL;

    editDist   = result.editDistance;
    alignLen   = ((nAend - nAbgn) + (Blen) + (editDist)) / 2;
    alignScore = 1000 - (int32)(1000.0 * editDist / alignLen);

    //  If there's an alignment, we can get a cigar string and better alignment length.
    if ((result.alignment != NULL) && (result.alignmentLength > 0)) {
      cigar      = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
      alignLen   = result.alignmentLength;
    }

    edlibFreeAlignResult(result);

    if (beVerbose)
      fprintf(stderr, " - POSITION from %9d-%-9d to %9d-%-9d score %5d/%9d = %4d%s%s\n",
              Abgn, Aend,
              nAbgn, nAend,
              editDist, alignLen, alignScore,
              (cigar != NULL) ? " align " : "",
              (cigar != NULL) ? cigar     : "");

    delete [] cigar;

    //  If it's a full alignment -- if the A region was big enough to have unaligned bases -- then
    //  we're done.  Update the result and get out of here.

    if (((Abgn  < nAbgn) || (Abgn == 0)) &&
        ((nAend < Aend)  || (Aend == Alen))) {

      Abgn  = nAbgn;
      Aend  = nAend;
      score = alignScore;

      return(true);
    }

    //  Otherwise, we ran out of A sequence to align to before we ran out of stuff to align.  Extend
    //  the A region and try again.

    if (Abgn == nAbgn)
      Abgn = max(Abgn - step, 0);

    if (Aend == nAend)
      Aend = min(Aend + step, Alen);

    goto tryAgain;
  }

  //  Didn't get a good alignment.

  //  We fail for one of two reasons - either not enough bases in the reference, or too high of
  //  error.  Unitigs are supposed to be from the same sequence, but they might be lower coverage
  //  and therefore higher error.  It's more likely they are misplaced.


  if ((Aend - Abgn < 10 * Blen) &&
      (maxEdit < Blen * 0.35)) {
    if (beVerbose)
      fprintf(stderr, " - FAILED, RELAX\n");

    Abgn = max(Abgn - step, 0);
    Aend = min(Aend + step, Alen);
    maxEdit *= 1.2;

    goto tryAgain;
  }

  if (beVerbose)
    fprintf(stderr, " - ABORT, ABORT, ABORT!\n");

  return(false);
}



bool
checkRecord(bedRecord   *record,
            sequences   &ctgs,
            sequences   &ctgs_orig,
            sequences   &utgs,
            bool         beVerbose,
            bool         UNUSED(doPlot)) {

  char   *Aseq = ctgs[record->_Aid].seq;
  char   *Bseq = utgs[record->_Bid].seq, *Brev = NULL;

  int32  Alen = ctgs[record->_Aid].len;
  int32  Blen = utgs[record->_Bid].len;

  double ratio = (double)Alen/ctgs_orig[record->_Aid].len;

  int32  Abgn  = max(0,    (int32)(floor(record->_bgn*ratio)));
  int32  Aend  = min(Alen, (int32)(ceil (record->_end*ratio)));

  bool   success    = true;
  int32  alignScore = 0;

  if (record->_Bfwd == false)
    Bseq = Brev = reverseComplementCopy(Bseq, Blen);

  //  If Bseq (the unitig) is small, just align the full thing.

  if (Blen < 50000) {
    success &= checkRecord_align("ALL",
                                 record->_Aname, Aseq, Alen, Abgn, Aend,
                                 record->_Bname, Bseq, Blen,
                                 alignScore,
                                 beVerbose);
  }

  //  Otherwise, we need to try to align only the ends of the unitig.
  //
  //        -----------------------[AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA]------------------
  //                                     BBBBB..............BBBBB
  //

  else {
    int32   AbgnL = Abgn,         AendL = Abgn + 50000;
    int32   AbgnR = Aend - 50000, AendR = Aend;

    char   *BseqL = Bseq;
    char   *BseqR = Bseq + Blen - 50000;

#if 0
    success &= checkRecord_align("ALL",
                                 record->_Aname, Aseq, Alen, Abgn, Aend,
                                 record->_Bname, Bseq, Blen,
                                 alignScore,
                                 beVerbose);
#endif

    success &= checkRecord_align("LEFT",
                                 record->_Aname, Aseq,  Alen, AbgnL, AendL,
                                 record->_Bname, BseqL, 50000,
                                 alignScore,
                                 beVerbose);

    success &= checkRecord_align("RIGHT",
                                 record->_Aname, Aseq,  Alen, AbgnR, AendR,
                                 record->_Bname, BseqR, 50000,
                                 alignScore,
                                 beVerbose);

    Abgn = AbgnL;
    Aend = AendR;
  }

  delete [] Brev;

  //  If successful, save the coordinates.  Because we're usually not aligning the whole
  //  unitig to the contig, we can't save the score.

  if (success) {
    record->_bgn   = Abgn;
    record->_end   = Aend;
    record->_score = 0;   //alignScore;
  }

  return(success);
}



//
//  Try to find an alignment for each link in the GFA file.  If found, output a new link
//  with correct CIGAR string.  If not found, discard the link.
//
void
processGFA(char     *tigName,
           uint32    tigVers,
           char     *inGFA,
           char     *otGFA,
           double    erate,
           uint32    verbosity) {

  //  Load the GFA file.

  fprintf(stderr, "-- Reading GFA '%s'.\n", inGFA);

  gfaFile   *gfa  = new gfaFile(inGFA);

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", tigName, tigVers-1);

  sequences *seqs_origp = new sequences(tigName, tigVers-1);
  sequences &seqs_orig  = *seqs_origp;

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", tigName, tigVers);

  sequences *seqsp = new sequences(tigName, tigVers);
  sequences &seqs  = *seqsp;

  //  Set GFA lengths based on the sequences we loaded.

  fprintf(stderr, "-- Resetting sequence lengths.\n");

  for (uint32 ii=0; ii<gfa->_sequences.size(); ii++)
    gfa->_sequences[ii]->_length = seqs[gfa->_sequences[ii]->_id].len;

  //  Align!

  uint32  passCircular = 0;
  uint32  failCircular = 0;

  uint32  passNormal = 0;
  uint32  failNormal = 0;

  uint32  iiLimit      = gfa->_links.size();
  uint32  iiNumThreads = omp_get_max_threads();
  uint32  iiBlockSize  = (iiLimit < 1000 * iiNumThreads) ? iiNumThreads : iiLimit / 999;

  fprintf(stderr, "-- Aligning " F_U32 " links using " F_U32 " threads and %.2f error rate.\n", iiLimit, iiNumThreads, erate*100);

#pragma omp parallel for schedule(dynamic, iiBlockSize)
  for (uint32 ii=0; ii<iiLimit; ii++) {
    gfaLink *link = gfa->_links[ii];

    if (link->_Aid == link->_Bid) {
      if (verbosity > 0)
        fprintf(stderr, "Processing circular link for tig %u\n", link->_Aid);

      if (link->_Afwd != link->_Bfwd)
        fprintf(stderr, "WARNING: %s %c %s %c -- circular to the same end!?\n",
                link->_Aname, link->_Afwd ? '+' : '-',
                link->_Bname, link->_Bfwd ? '+' : '-');

      bool  pN = checkLink(link, seqs, seqs_orig, erate, (verbosity > 0), false);

      if (pN == true)
        passCircular++;
      else
        failCircular++;
    }

    //  Now the usual case.

    else {
      if (verbosity > 0)
        fprintf(stderr, "Processing link between tig %u %s and tig %u %s\n",
                link->_Aid, link->_Afwd ? "-->" : "<--",
                link->_Bid, link->_Bfwd ? "-->" : "<--");

      bool  pN = checkLink(link, seqs, seqs_orig, erate, (verbosity > 0), false);

      if (pN == true)
        passNormal++;
      else
        failNormal++;
    }

    //  If the cigar exists, we found an alignment.  If not, delete the link.

    if (link->_cigar == NULL) {
      if (verbosity > 0)
        fprintf(stderr, "  Failed to find alignment.\n");
      delete gfa->_links[ii];
      gfa->_links[ii] = NULL;
    }
  }

  fprintf(stderr, "-- Writing GFA '%s'.\n", otGFA);

  gfa->saveFile(otGFA);

  fprintf(stderr, "-- Cleaning up.\n");

  delete seqsp;
  delete gfa;

  fprintf(stderr, "-- Aligned %6u ciruclar tigs, failed %6u\n", passCircular, failCircular);
  fprintf(stderr, "-- Aligned %6u   linear tigs, failed %6u\n", passNormal,   failNormal);
}



//
//  Find an alignment between the unitig (the feature) and the contig (the 'chromosome').
//  Output updated coordiates.
//
void
processBED(char   *tigName,
           uint32  tigVers,
           char   *seqName,
           uint32  seqVers,
           char   *inBED,
           char   *otBED,
           uint32  verbosity) {

  //  Load the BED file.

  fprintf(stderr, "-- Reading BED '%s'.\n", inBED);

  bedFile   *bed  = new bedFile(inBED);

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", tigName, tigVers);

  sequences *utgsp = new sequences(tigName, tigVers);
  sequences &utgs  = *utgsp;

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", seqName, seqVers-1);

  sequences *ctgs_origp = new sequences(seqName, seqVers-1);
  sequences &ctgs_orig  = *ctgs_origp;

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", seqName, seqVers);

  sequences *ctgsp = new sequences(seqName, seqVers);
  sequences &ctgs  = *ctgsp;

  //  Align!

  uint32  pass = 0;
  uint32  fail = 0;

  uint32  iiLimit      = bed->_records.size();
  uint32  iiNumThreads = omp_get_max_threads();
  uint32  iiBlockSize  = (iiLimit < 1000 * iiNumThreads) ? iiNumThreads : iiLimit / 999;

  fprintf(stderr, "-- Aligning " F_U32 " records using " F_U32 " threads.\n", iiLimit, iiNumThreads);

#pragma omp parallel for schedule(dynamic, iiBlockSize)
  for (uint32 ii=0; ii<iiLimit; ii++) {
    bedRecord *record = bed->_records[ii];

    if (checkRecord(record, ctgs, ctgs_orig, utgs, (verbosity > 0), false)) {
      pass++;
    } else {
      delete bed->_records[ii];
      bed->_records[ii] = NULL;
      fail++;
    }
  }

  fprintf(stderr, "-- Writing BED '%s'.\n", otBED);

  bed->saveFile(otBED);

  fprintf(stderr, "-- Cleaning up.\n");

  delete utgsp;
  delete ctgsp;
  delete bed;

  fprintf(stderr, "-- Aligned %6u unitigs to contigs, failed %6u\n", pass, fail);
}



//
//  Infer a graph from the positions of unitigs (features) in contigs (chromosomes).  Generate a GFA
//  input and toss that up to processGFA.
//
void
processBEDtoGFA(char   *tigName,
                uint32  tigVers,
                char   *inBED,
                char   *otGFA,
                double  erate,
                uint32  verbosity) {

  int32  minOlap = 100;

  //  We only really need the sequence lengths here, but eventually, we'll want to generate
  //  alignments for all the overlaps, and so we'll need the sequences too.
  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", tigName, tigVers-1);

  sequences *seqs_origp = new sequences(tigName, tigVers-1);
  sequences &seqs_orig  = *seqs_origp;

  fprintf(stderr, "-- Loading sequences from tigStore '%s' version %u.\n", tigName, tigVers);

  sequences *seqsp = new sequences(tigName, tigVers);
  sequences &seqs  = *seqsp;

  //  Load the BED file and allocate an output GFA.

  fprintf(stderr, "-- Reading BED '%s'.\n", inBED);

  bedFile   *bed  = new bedFile(inBED);
  gfaFile   *gfa  = new gfaFile("H\tVN:Z:1.0");

  //  Iterate over sequences, looking for overlaps in contigs.  Stupid, O(n^2) but seems fast enough.

  uint32  iiLimit      = bed->_records.size();
  uint32  iiNumThreads = omp_get_max_threads();
  uint32  iiBlockSize  = (iiLimit < 1000 * iiNumThreads) ? iiNumThreads : iiLimit / 999;

  fprintf(stderr, "-- Aligning " F_U32 " records using " F_U32 " threads.\n", iiLimit, iiNumThreads);

#pragma omp parallel for schedule(dynamic, iiBlockSize)
  for (uint64 ii=0; ii<bed->_records.size(); ii++) {
    for (uint64 jj=ii+1; jj<bed->_records.size(); jj++) {

      if (bed->_records[ii]->_Aid != bed->_records[jj]->_Aid)                 //  Different contigs?
        continue;                                                             //  No overlap.

      if ((bed->_records[ii]->_end < bed->_records[jj]->_bgn + minOlap) ||    //  No (thick) intersection?
          (bed->_records[jj]->_end < bed->_records[ii]->_bgn + minOlap))      //
        continue;                                                             //  No overlap.

      //  Overlap!

      //fprintf(stderr, "OVERLAP %s %d-%d - %s %d-%d\n",
      //        bed->_records[ii]->_Bname, bed->_records[ii]->_bgn, bed->_records[ii]->_end,
      //        bed->_records[jj]->_Bname, bed->_records[jj]->_bgn, bed->_records[jj]->_end);

      int32  olapLen = 0;

      if (bed->_records[ii]->_bgn < bed->_records[jj]->_end)
        olapLen = bed->_records[ii]->_end - bed->_records[jj]->_bgn;

      if (bed->_records[jj]->_bgn < bed->_records[ii]->_end)
        olapLen = bed->_records[jj]->_end - bed->_records[ii]->_bgn;

      assert(olapLen > 0);

      char   cigar[81];

      sprintf(cigar, "%dM", olapLen);

      gfaLink *link = new gfaLink(bed->_records[ii]->_Bname, bed->_records[ii]->_Bid, true,
                                  bed->_records[jj]->_Bname, bed->_records[jj]->_Bid, true,
                                  cigar);

      bool  pN = checkLink(link, seqs, seqs_orig, erate, (verbosity > 0), false);

#pragma omp critical
      {
        if (pN)
          gfa->_links.push_back(link);
        else
          gfa->_links.push_back(link);

        //  Remember sequences we've hit.

        seqs.used[bed->_records[ii]->_Bid]++;
        seqs.used[bed->_records[jj]->_Bid]++;
      }
    }
  }

  //  Add sequences.  We could have done this as we're running through making edges, but we then
  //  need to figure out if we've seen a sequence already.

  char   seqName[80];

  for (uint32 ii=0; ii<seqs.e; ii++)
    if (seqs.used[ii] > 0) {
      sprintf(seqName, "utg%08u", ii);
      gfa->_sequences.push_back(new gfaSequence(seqName, ii, seqs[ii].len));
    }

  //  Write the file, cleanup, done!

  gfa->saveFile(otGFA);

  delete gfa;
  delete bed;

  delete seqsp;
}



int
main (int argc, char **argv) {
  char     *tigName         = NULL;         //  For GFA and BED, the source of the tigs
  uint32    tigVers         = UINT32_MAX;

  char     *seqName         = NULL;         //  For BED, the source of the 'chromosomes'
  uint32    seqVers         = UINT32_MAX;   //  The -C option (either chromosome or container)

  char     *inGraph         = NULL;
  char     *otGraph         = NULL;

  uint32    graphType       = IS_GFA;

  uint32    verbosity      = 0;

  // not used for bed alignments, that is fixed at a lower value
  double    erate          = 0.10;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-C") == 0) {
      seqName = argv[++arg];
      seqVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-gfa") == 0) {
      graphType = IS_GFA;
    } else if (strcmp(argv[arg], "-bed") == 0) {
      graphType = IS_BED;

    } else if (strcmp(argv[arg], "-i") == 0) {
      inGraph = argv[++arg];
    } else if (strcmp(argv[arg], "-o") == 0) {
      otGraph = argv[++arg];

    } else if (strcmp(argv[arg], "-V") == 0) {
      verbosity++;

    } else if (strcmp(argv[arg], "-t") == 0) {
      omp_set_num_threads(atoi(argv[++arg]));

    } else if (strcmp(argv[arg], "-e") == 0) {
      erate = atof(argv[++arg]);

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (tigName == NULL)
    err++;
  if (inGraph == NULL)
    err++;
  if (otGraph == NULL)
    err++;

  if ((tigName) && (tigVers == 0))
    err++;
  if ((seqName) && (seqVers == 0))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  Validates a GFA by generating alignments.\n");
    fprintf(stderr, "  Optionally writes new GFA with updated CIGAR string (NOT IMPLEMENTED).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -T t v         Load tigs from tgStore 't', version 'v'.\n");
    fprintf(stderr, "    -C t v         For BED format, the source of the 'chromosomes'.  Similar to -T.\n");
    fprintf(stderr, "                     Consensus sequence must exist for -T and -C (usually in v=2)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -i input       Input graph.\n");
    fprintf(stderr, "    -o output      Output graph.\n");
    fprintf(stderr, "                     Graph are either GFA (v1) or BED format.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -gfa           The input and output graphs are in GFA (v1) format.\n");
    fprintf(stderr, "    -bed           The input graph is in BED format.  If -C is supplied, the\n");
    fprintf(stderr, "                   output will also be BED, and will have updated positions.\n");
    fprintf(stderr, "                   If -C is not supplied, the output will be GFA (v1) of the\n");
    fprintf(stderr, "                   overlaps inferred from the BED positions.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -V             Increase chatter.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t threads     Use 'threads' computational threads.\n");
    fprintf(stderr, "\n");

    if (tigName == NULL)
      fprintf(stderr, "ERROR: no tigStore (-T) supplied.\n");
    if (inGraph == NULL)
      fprintf(stderr, "ERROR: no input GFA (-i) supplied.\n");
    if (otGraph == NULL)
      fprintf(stderr, "ERROR: no output GFA (-o) supplied.\n");

    if ((tigName) && (tigVers == 0))
      fprintf(stderr, "ERROR: invalid tigStore version (-T) supplied.\n");
    if ((seqName) && (seqVers == 0))
      fprintf(stderr, "ERROR: invalid tigStore version (-C) supplied.\n");

    exit(1);
  }

  if (graphType == IS_GFA)
    processGFA(tigName, tigVers, inGraph, otGraph, erate, verbosity);

  if ((graphType == IS_BED) && (seqName != NULL))
    processBED(tigName, tigVers, seqName, seqVers, inGraph, otGraph, verbosity);

  if ((graphType == IS_BED) && (seqName == NULL))
    processBEDtoGFA(tigName, tigVers, inGraph, otGraph, erate, verbosity);

  fprintf(stderr, "Bye.\n");

  exit(0);
}
