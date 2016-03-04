
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
 *    src/alignment/AS_ALN_forcns.C
 *    src/alignment/alignment-drivers.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-DEC-23 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "aligners.H"
#include "gkStore.H"  //  For AS_MAX_READLEN
#include "AS_UTL_reverseComplement.H"




class Optimal_Overlap_Data {
public:
  char         h_alignA[AS_MAX_READLEN + AS_MAX_READLEN + 2];
  char         h_alignB[AS_MAX_READLEN + AS_MAX_READLEN + 2];
  int          h_trace[AS_MAX_READLEN + AS_MAX_READLEN + 2];
  ALNoverlap   o;
};


ALNoverlap *
Optimal_Overlap_AS_forCNS(char *a, char *b,
                          int begUNUSED, int endUNUSED,
                          int ahang, int bhang,
                          int opposite,
                          double erate, double thresh, int minlen,
                          CompareOptions what) {
  static Optimal_Overlap_Data  *ood = NULL;

  if (ood == NULL)
    ood = new Optimal_Overlap_Data;

  memset(ood->h_alignA, 0, sizeof(char) * (AS_MAX_READLEN + AS_MAX_READLEN + 2));
  memset(ood->h_alignB, 0, sizeof(char) * (AS_MAX_READLEN + AS_MAX_READLEN + 2));
  memset(ood->h_trace,  0, sizeof(int)  * (AS_MAX_READLEN + AS_MAX_READLEN + 2));

  alignLinker_s   al;

  //if (VERBOSE_MULTIALIGN_OUTPUT >= 3)
  //  fprintf(stderr, "Optimal_Overlap_AS_forCNS()--  Begins\n");

#if 0
  if (erate > AS_MAX_ERROR_RATE) {
    //fprintf(stderr, "Optimal_Overlap_AS_forCNS()--  erate=%f >= AS_MAX_ERROR_RATE=%f, reset to max\n", erate, (double)AS_MAX_ERROR_RATE);
    erate = AS_MAX_ERROR_RATE;
  }
  assert((0.0 <= erate) && (erate <= AS_MAX_ERROR_RATE));
#endif

  if (opposite)
    reverseComplementSequence(b, strlen(b));

  //if (VERBOSE_MULTIALIGN_OUTPUT >= 3) {
  //  fprintf(stderr, "ALIGN %s\n", a);
  //  fprintf(stderr, "ALIGN %s\n", b);
  //}

  alignLinker(ood->h_alignA,
              ood->h_alignB,
              a,
              b,
              &al,
              true,   //  Looking for global end-to-end alignments
              false,  //  Count matches to N as matches
              ahang, bhang);
   if (al.alignLen == 0) {
      return NULL;
   }

   //if (VERBOSE_MULTIALIGN_OUTPUT >= 3) {
   //  fprintf(stderr, "ALIGN %d %d-%d %d-%d opposite=%d\n", al.alignLen, al.begI, al.endI, al.begJ, al.endJ, opposite);
   //  fprintf(stderr, "ALIGN '%s'\n", ood->h_alignA);
   //  fprintf(stderr, "ALIGN '%s'\n", ood->h_alignB);
   //}

  if (opposite) {
    reverseComplementSequence(b, strlen(b));

    reverseComplementSequence(ood->h_alignA, al.alignLen);
    reverseComplementSequence(ood->h_alignB, al.alignLen);

    int x = al.begJ;
    al.begJ = al.lenB - al.endJ;
    al.endJ = al.lenB - x;
  }

  //  We don't expect partial overlaps here.  At least one fragment
  //  must have an alignment to the very start.
  //
  //  ECR depends on this return value; it is allowed to fail
  //  when building a new unitig multialign.  For example:
  //
  //  <-----------------------
  //        ------>
  //
  //  When ECR tries to extend the second fragment, it checks that
  //  the extended fragment overlaps the next contig.  It does not
  //  check that the extended bits agree with the first fragment,
  //  leaving that up to "does the unitig rebuild".
  //
  if ((al.begJ != 0) && (al.begI != 0))
    return(NULL);

  ood->o.begpos  = (al.begI           > 0) ? (al.begI)           : -(al.begJ);
  ood->o.endpos  = (al.lenB - al.endJ > 0) ? (al.lenB - al.endJ) : -(al.lenA - al.endI);
  ood->o.length  = al.alignLen;
  ood->o.diffs   = 0;
  ood->o.comp    = opposite;
  ood->o.trace   = ood->h_trace;

  {
    int x=0;

    int tp = 0;
    int ap = al.begI;
    int bp = al.begJ;

    for (x=0; x<al.alignLen; x++) {
      if (ood->h_alignA[x] == '-') {
        ood->h_trace[tp++] = -(ap + 1);
        ap--;
      }
      if (ood->h_alignB[x] == '-') {
        ood->h_trace[tp++] = bp + 1;
        bp--;
      }

      //  Count the differences.
      //
      //  But allow N's and lowercase as matches.  If either letter is N, then the other letter is
      //  NOT N (if both letters were N, both would be lowercase n, representing a match).  This
      //  just subtracts out the diff we added in above.
      //
      bool  diff   = false;
      bool  ignore = false;

      if (toupper(ood->h_alignA[x]) != toupper(ood->h_alignB[x]))
        diff = true;

      if ((ood->h_alignA[x] == 'N') || (ood->h_alignA[x] == 'n') ||
          (ood->h_alignB[x] == 'N') || (ood->h_alignB[x] == 'n'))
        ignore = true;

      if (islower(ood->h_alignA[x]) && (ood->h_alignB[x] == '-'))
        ignore = true;

      if ((diff == true) && (ignore == false))
        ood->o.diffs++;

      bp++;
      ap++;
    }

    ood->h_trace[tp] = 0;

    //if (VERBOSE_MULTIALIGN_OUTPUT >= 4) {
    //  fprintf(stderr, "trace");
    //  for (x=0; x<tp; x++)
    //    fprintf(stderr, " %d", ood->h_trace[x]);
    //  fprintf(stderr, "\n");
    //  fprintf(stderr, "A: %4d-%4d %4d %s\n", al.begI, al.endI, al.lenA, ood->h_alignA);
    //  fprintf(stderr, "B: %4d-%4d %4d %s\n", al.begJ, al.endJ, al.lenB, ood->h_alignB);
    //}
  }

  //if (VERBOSE_MULTIALIGN_OUTPUT >= 3) {
  //  fprintf(stderr, "ERATE:   diffs=%d / length=%d = %f\n", ood->o.diffs, ood->o.length, (double)ood->o.diffs / ood->o.length);
  //  fprintf(stderr, "Optimal_Overlap_AS_forCNS()--  Ends\n");
  //}

  if ((double)ood->o.diffs / ood->o.length <= erate)
    return(&ood->o);

  return(NULL);
}
