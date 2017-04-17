
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2017-APR-04
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"

#include "edlib.H"

#include "splitToWords.H"
#include "AS_UTL_reverseComplement.H"

#include "gfa.H"


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
    len = tig->length(false);
    seq = new char [len + 1];

    memcpy(seq, tig->bases(false), len);

    seq[len] = 0;
  };

  char   *seq;
  uint32  len;
};



void
dotplot(uint32 Aid, bool Afwd, char *Aseq,
        uint32 Bid, bool Bfwd, char *Bseq) {
  char   Aname[128], Afile[128];
  char   Bname[128], Bfile[128];
  char   Pname[128], Pfile[128];
  FILE  *F;

  sprintf(Aname, "tig%08u%c",       Aid, (Afwd) ? '+' : '-');
  sprintf(Afile, "tig%08u%c.fasta", Aid, (Afwd) ? '+' : '-');
  sprintf(Bname, "tig%08u%c",       Bid, (Bfwd) ? '+' : '-');
  sprintf(Bfile, "tig%08u%c.fasta", Bid, (Bfwd) ? '+' : '-');
  sprintf(Pname, "plot-%s-%s",    Aname, Bname);
  sprintf(Pfile, "plot-%s-%s.sh", Aname, Bname);

  F = fopen(Pfile, "w");
  fprintf(F, "#!/bin/sh\n");
  fprintf(F, "\n");
  fprintf(F, "nucmer --maxmatch --nosimplify -p %s %s.fasta %s.fasta\n", Pname, Aname, Bname);
  fprintf(F, "show-coords -l -o -r -T %s.delta | expand -t 8 > %s.coords\n", Pname, Pname);
  fprintf(F, "mummerplot --fat -t png -p %s %s.delta\n", Pname, Pname);
  fprintf(F, "echo mummerplot --fat -p %s %s.delta\n", Pname, Pname);
  fclose(F);

  F = fopen(Afile, "w");
  fprintf(F, ">%s\n%s\n", Aname, Aseq);
  fclose(F);

  F = fopen(Bfile, "w");
  fprintf(F, ">%s\n%s\n", Bname, Bseq);
  fclose(F);

  sprintf(Pfile, "sh plot-%s-%s.sh", Aname, Bname);

  system(Pfile);
}




bool
checkLink(gfaLink  *link,
          sequence *seqs,
          bool      beVerbose,
          bool      doPlot) {

  char   *Aseq = seqs[link->_Aid].seq;
  char   *Bseq = seqs[link->_Bid].seq;

  int32  Abgn, Aend, Alen = seqs[link->_Aid].len;
  int32  Bbgn, Bend, Blen = seqs[link->_Bid].len;

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
    reverseComplementSequence(Aseq, Alen);
  if (link->_Bfwd == false)
    reverseComplementSequence(Bseq, Blen);


  //  Ty to find the end coordinate on B.  Align the last bits of A to B.
  //
  //   -------(---------]     v--??
  //             [------------)------
  //

  Abgn = max(Alen - AalignLen, 0);
  Aend =     Alen;

  Bbgn = 0;
  Bend = min(Blen, (int32)(1.10 * BalignLen));  //  Allow 25% gaps over what the GFA said?

  maxEdit = (int32)ceil(alignLen * 0.12);

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
  } else {
    if (beVerbose)
      fprintf(stderr, " - FAILED\n");
  }

  //  Do the same for A.  Aend and Bbgn never change; Bend was set above.
  //
  //   ------(--------------]
  //         ^--??  [-------]-----------
  //

  Abgn = max(Alen - (int32)(1.10 * AalignLen), 0);  //  Allow 25% gaps over what the GFA said?

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

  //  Cleanup for the next link->

  if (link->_Afwd == false)
    reverseComplementSequence(Aseq, Alen);
  if (link->_Bfwd == false)
    reverseComplementSequence(Bseq, Blen);

  if (beVerbose)
    fprintf(stderr, "\n");

  return(success);
}




int
main (int argc, char **argv) {
  char    *tigName         = NULL;
  uint32   tigVers         = UINT32_MAX;

  char    *inGFA           = NULL;
  char    *otGFA           = NULL;

  uint32    verbosity      = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

      if (tigVers == 0)
        fprintf(stderr, "invalid tigStore version (-T store version partition) '-t %s %s %s'.\n", argv[arg-2], argv[arg-1], argv[arg]), exit(1);

    } else if (strcmp(argv[arg], "-i") == 0) {
      inGFA = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      otGFA = argv[++arg];

    } else if (strcmp(argv[arg], "-V") == 0) {
      verbosity++;

    } else if (strcmp(argv[arg], "-t") == 0) {
      omp_set_num_threads(atoi(argv[++arg]));

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (tigName == NULL)
    err++;
  if (inGFA == NULL)
    err++;
  if (otGFA == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  Validates a GFA by generating alignments.\n");
    fprintf(stderr, "  Optionally writes new GFA with updated CIGAR string (NOT IMPLEMENTED).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -G g           Load reads from gkStore 'g'\n");
    fprintf(stderr, "    -T t v         Load tigs from tgStore 't', version 'v'.\n");
    fprintf(stderr, "                     Consensus sequence must exist (usually in v=2)\n");
    fprintf(stderr, "    -i input.gfa\n");
    fprintf(stderr, "    -o output.gfa\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -V             Increase chatter\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -t threads     Use 'threads' computational threads.\n");
    fprintf(stderr, "\n");

    if (tigName == NULL)
      fprintf(stderr, "ERROR: no tigStore (-T) supplied.\n");
    if (inGFA == NULL)
      fprintf(stderr, "ERROR: no input GFA (-i) supplied.\n");
    if (otGFA == NULL)
      fprintf(stderr, "ERROR: no output GFA (-o) supplied.\n");

    exit(1);
  }

  fprintf(stderr, "-- Opening tigStore '%s' version %u.\n", tigName, tigVers);
  tgStore *tigStore = new tgStore(tigName, tigVers);

  //  Load the GFA file.

  fprintf(stderr, "-- Reading GFA '%s'.\n", inGFA);
  gfaFile  *gfa = new gfaFile(inGFA);

  //  Load all consensus sequences

  uint32  b = 0;
  uint32  e = tigStore->numTigs();

  sequence  *seqs = new sequence [e+1];

  fprintf(stderr, "-- Loading tigs %u to %u.\n", b, e);

  for (uint32 ti=b; ti < e; ti++) {
    tgTig *tig = tigStore->loadTig(ti);

    if (tig == NULL)
      continue;

    seqs[ti].set(tig);

    tigStore->unloadTig(ti);
  }

  //  Set GFA lengths based on the sequences we loaded.

  fprintf(stderr, "-- Resetting sequence lengths.\n", inGFA);

  for (uint32 ii=0; ii<gfa->_sequences.size(); ii++)
    gfa->_sequences[ii]->_length = seqs[gfa->_sequences[ii]->_id].len;

  //  Done with the stores.

  fprintf(stderr, "-- Closing tigStore '%s'.\n", tigName);

  delete tigStore;

  //  Align!

  uint32  passCircular = 0;
  uint32  failCircular = 0;

  uint32  passNormal = 0;
  uint32  failNormal = 0;

  uint32  iiLimit      = gfa->_links.size();
  uint32  iiNumThreads = omp_get_max_threads();
  uint32  iiBlockSize  = (iiLimit < 1000 * iiNumThreads) ? iiNumThreads : iiLimit / 999;

  fprintf(stderr, "-- Aligning " F_U32 " links using " F_U32 " threads.\n", iiLimit, iiNumThreads);

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

      bool  pN = checkLink(link, seqs, (verbosity > 0), false);

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

      bool  pN = checkLink(link, seqs, (verbosity > 0), false);

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

  delete [] seqs;
  delete    gfa;

  fprintf(stderr, "-- Aligned %6u ciruclar tigs, failed %6u\n", passCircular, failCircular);
  fprintf(stderr, "-- Aligned %6u   linear tigs, failed %6u\n", passNormal,   failNormal);
  fprintf(stderr, "-- Bye.\n");

  exit(0);
}
