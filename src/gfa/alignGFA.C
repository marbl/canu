 
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
 *    Brian P. Walenz beginning on 2017-MAR-23
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




class gfaEdge {
public:
  gfaEdge(char *inLine) {
    splitToWords  W(inLine);

    strcpy(Aname, W[1]);
    strcpy(Bname, W[3]);

    Aid   = atoi(Aname + 3);
    Bid   = atoi(Bname + 3);

    Afwd  = (W[2][0] == '+');
    Bfwd  = (W[4][0] == '+');

    cigar    = NULL;

    olapLen = atoi(W[5]);
  };

  ~gfaEdge() {
    delete [] cigar;
    cigar = NULL;
  };


  char    Aname[128];
  uint32  Aid;
  bool    Afwd;

  char    Bname[128];
  uint32  Bid;
  bool    Bfwd;

  char   *cigar;

  int32   olapLen;
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
checkEdge(gfaEdge  &edge,
          sequence *seqs,
          bool      beVerbose,
          bool      doPlot) {

  char   *Aseq = seqs[edge.Aid].seq;
  char   *Bseq = seqs[edge.Bid].seq;

  int32  Abgn, Aend, Alen = seqs[edge.Aid].len;
  int32  Bbgn, Bend, Blen = seqs[edge.Bid].len;

  EdlibAlignResult  result  = { 0, NULL, NULL, 0, NULL, 0, 0 };

  int32             maxEdit  = (int32)ceil(edge.olapLen * 0.12);

  if (edge.Afwd == false)
    reverseComplementSequence(Aseq, Alen);
  if (edge.Bfwd == false)
    reverseComplementSequence(Bseq, Blen);


  //  Ty to find the end coordinate on B.  Align the last bits of A to B.
  //
  //   -------(---------]     v--??
  //             [------------)------
  //

  Abgn = max(Alen - edge.olapLen, 0);
  Aend =     Alen;

  Bbgn = 0;
  Bend = min(Blen, (int32)(1.25 * edge.olapLen));  //  Allow 25% gaps?

  if (beVerbose)
    fprintf(stdout, "TEST tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  %6dM  (extend B)\n",
            edge.Aid, (edge.Afwd) ? '+' : '-', Abgn, Aend,
            edge.Bid, (edge.Bfwd) ? '+' : '-', Bbgn, Bend,
            edge.olapLen);

  result = edlibAlign(Aseq + Abgn, Aend-Abgn,
                      Bseq + Bbgn, Bend-Bbgn,
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  if (result.numLocations > 0) {
    Bend = Bbgn + result.endLocations[0] + 1;  // 0-based to space-based
    edlibFreeAlignResult(result);
  }

  //  Do the same for A.  Aend and Bbgn never change; Bend was set above.
  //
  //   ------(--------------]
  //         ^--??  [-------]-----------
  //

  Abgn = max(Alen - (int32)(1.25 * edge.olapLen), 0);

  if (beVerbose)
    fprintf(stdout, "     tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  %6dM  (extend A)\n",
            edge.Aid, (edge.Afwd) ? '+' : '-', Abgn, Aend,
            edge.Bid, (edge.Bfwd) ? '+' : '-', Bbgn, Bend,
            edge.olapLen);

  result = edlibAlign(Bseq + Bbgn, Bend-Bbgn,
                      Aseq + Abgn, Aend-Abgn,
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_HW, EDLIB_TASK_LOC));

  if (result.numLocations > 0) {
    Abgn = Abgn + result.startLocations[0];
    edlibFreeAlignResult(result);
  }

  //  One more alignment, this time, with feeling - notice EDLIB_MODE_MW and EDLIB_TASK_PATH.

  result = edlibAlign(Aseq + Abgn, Aend-Abgn,
                      Bseq + Bbgn, Bend-Bbgn,
                      edlibNewAlignConfig(maxEdit, EDLIB_MODE_NW, EDLIB_TASK_PATH));


  int32  editDist = 0;
  int32  alignLen = 0;
  bool   success = false;

  if (result.numLocations > 0) {
    editDist = result.editDistance;
    alignLen = ((Aend - Abgn) + (Bend - Bbgn) + (editDist)) / 2;
    alignLen = result.alignmentLength;     //  Edlib 'alignmentLength' is populated only for TASK_PATH

    edge.cigar = edlibAlignmentToCigar(result.alignment,
                                       result.alignmentLength, EDLIB_CIGAR_STANDARD);

    edlibFreeAlignResult(result);

    success = true;
  }

  if (beVerbose)
    fprintf(stdout, "     tig%08u %c %8d-%-8d    tig%08u %c %8d-%-8d  %6dM  (final)  %6dM %.4f\n",
            edge.Aid, (edge.Afwd) ? '+' : '-', Abgn, Aend,
            edge.Bid, (edge.Bfwd) ? '+' : '-', Bbgn, Bend,
            edge.olapLen,
            alignLen, (double)editDist / alignLen);

  //  Make a plot.

  if ((success == false) && (doPlot == true))
    dotplot(edge.Aid, edge.Afwd, Aseq,
            edge.Bid, edge.Bfwd, Bseq);

  //  Cleanup for the next edge.

  if (edge.Afwd == false)
    reverseComplementSequence(Aseq, Alen);
  if (edge.Bfwd == false)
    reverseComplementSequence(Bseq, Blen);

  if (beVerbose)
    fprintf(stdout, "\n");

  return(success);
}




int
main (int argc, char **argv) {
  char    *gkpName         = NULL;

  char    *tigName         = NULL;
  uint32   tigVers         = UINT32_MAX;

  char    *inGFA           = NULL;
  char    *otGFA           = NULL;

  uint32    verbosity      = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
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

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "  Validates a GFA by generating alignments.\n");
    fprintf(stderr, "  Optionally writes new GFA with updated CIGAR string (NOT IMPLEMENTED).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -G g           Load reads from gkStore 'g'\n");
    fprintf(stderr, "    -T t v         Load tigs from tgStore 't', version 'v'.\n");
    fprintf(stderr, "                     Consensus sequence must exist (usually in v=2)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -V             Increase chatter\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -i input.gfa\n");
    fprintf(stderr, "    -o output.gfa\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  gkStore  *gkpStore          = NULL;
  tgStore  *tigStore          = NULL;

  if (gkpName) {
    fprintf(stderr, "-- Opening gkpStore '%s'.\n", gkpName);
    gkpStore = gkStore::gkStore_open(gkpName, gkStore_readOnly);
  }

  if (tigName) {
    fprintf(stderr, "-- Opening tigStore '%s' version %u.\n", tigName, tigVers);
    tigStore = new tgStore(tigName, tigVers);
  }

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

  //  Done with the stores.

  gkpStore->gkStore_close();

  delete tigStore;

  //  Open the GFAs.

  errno = 0;
  FILE *in = fopen(inGFA, "r");
  if (errno)
    fprintf(stderr, "Failed to open input GFA '%s': %s\n", inGFA, strerror(errno)), exit(1);

  errno = 0;
  FILE *ot = (otGFA != NULL) ? fopen(otGFA, "w") : NULL;
  if (errno)
    fprintf(stderr, "Failed to open input GFA '%s': %s\n", otGFA, strerror(errno)), exit(1);

  //  Read the garbage.

  char L[1048576];

  uint32  passCircular = 0;
  uint32  failCircular = 0;

  uint32  passNormal = 0;
  uint32  failNormal = 0;

  for (fgets(L, 1048576, in);
       feof(in) == false;
       fgets(L, 1048576, in)) {

    //  Bail on non links.
    if (L[0] != 'L') {
      if (ot)
        fputs(L, ot);
      continue;
    }

    //  Find the data.

    gfaEdge   edge(L);

    //  Handle circular stuff.  Some care -- not done yet -- needs to be taken so that we
    //  don't just align self-to-self:
    //    L tig00001121 + tig00001121 + 31504M
    //    L tig00001121 - tig00001121 - 31504M

    if (edge.Aid == edge.Bid) {
      fprintf(stderr, "Processing circular edge for tig %u\n", edge.Aid);

      if (edge.Afwd != edge.Bfwd)
        fprintf(stderr, "WARNING: %s %c %s %c -- circular to the same end!?\n",
                edge.Aname, edge.Afwd ? '+' : '-',
                edge.Bname, edge.Bfwd ? '+' : '-');

      bool  pN = checkEdge(edge, seqs, false, false);

      if (pN == true)
        passCircular++;
      else
        failCircular++;
    }

    //  Now the usual case.

    else {
      fprintf(stderr, "Processing edge between tig %u %s and tig %u %s\n",
              edge.Aid, edge.Afwd ? "-->" : "<--",
              edge.Bid, edge.Bfwd ? "-->" : "<--");

      bool  pN = checkEdge(edge, seqs, false, false);

      if (pN == true)
        passNormal++;
      else
        failNormal++;
    }

    //  If the cigar exists, we found an alignment.

    if ((ot) && (edge.cigar))
      fprintf(ot, "L\t%s\t%c\t%s\t%c\t%s\n",
              edge.Aname, edge.Afwd ? '+' : '-',
              edge.Bname, edge.Bfwd ? '+' : '-',
              edge.cigar);

    delete [] edge.cigar;
    edge.cigar = NULL;
  }

  delete [] seqs;

  if (in)
    fclose(in);
  if (ot)
    fclose(ot);

  fprintf(stderr, "passCircular  %u\n", passCircular);
  fprintf(stderr, "failCircular  %u\n", failCircular);
  fprintf(stderr, "\n");
  fprintf(stderr, "passNormal    %u\n", passNormal);
  fprintf(stderr, "failNormal    %u\n", failNormal);

  exit(0);
}














