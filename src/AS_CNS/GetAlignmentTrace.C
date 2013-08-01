
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

static char *rcsid = "$Id$";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"
#include "MicroHetREZ.H"
#include "AS_UTL_reverseComplement.H"

#define SHOW_ATTEMPT   2
#define SHOW_ACCEPTED  3

int32    numScores = 0;
double lScoreAve = 0.0;
double aScoreAve = 0.0;
double bScoreAve = 0.0;

double acceptThreshold = 0.1;  //1.0 / 3.0;

// init value is 200; this could be set to the amount you extend the clear
// range of seq b, plus 10 for good measure
extern int32 MaxBegGap;

// init value is 200; this could be set to the amount you extend the
// clear range of seq a, plus 10 for good measure
extern int32 MaxEndGap;


typedef struct CNS_AlignParams {
  int32 bandBgn;      //  A sequence band begin/end - the align must
  int32 bandEnd;      //  fall in this range
  int32 maxBegGap;
  int32 maxEndGap;
  int32 ahang;        //  Actual hangs, used by Optimal_Overlap_AS_forCNS
  int32 bhang;
  int32 opposite;
  double erate;
  double thresh;
  int32 minlen;
  CompareOptions what;
} CNS_AlignParams;

static
void
CNS_AlignParams_init(CNS_AlignParams *ap) {
  //CNS_AlignParams LOCAL_DEFAULT_PARAMS = {0,0,0,0,0,AS_CNS_ERROR_RATE,CNS_DP_THRESH,CNS_DP_MINLEN,AS_FIND_ALIGN};

  ap->bandBgn   = 0;
  ap->bandEnd   = 0;
  ap->maxBegGap = 0;
  ap->maxEndGap = 0;
  ap->ahang     = 0;
  ap->bhang     = 0;
  ap->opposite  = 0;
  ap->erate     = AS_CNS_ERROR_RATE;
  ap->thresh    = CNS_DP_THRESH;
  ap->minlen    = CNS_DP_MINLEN;
  ap->what      = AS_FIND_ALIGN;
}



static
int
InvertTrace(int32 alen, int32 blen, ALNoverlap *O) {
  int32 aend=alen+2;
  int32 bend=blen+2;
  int32 n_dels=0;
  int32 *otrace=O->trace;
  int32 *t=otrace;
  int32 *s;
  int32 tmp;
  while ( *t != 0 ) {
    n_dels++; t++;
  }
  t=otrace;
  s=t+n_dels-1;
  while (  s - t > 0 ) {
    tmp = *t;
    if ( *s < 0 ) {
      *t = - (aend + *s);
    } else {
      *t = (bend - *s);
    }
    if ( tmp < 0 ) {
      *s = - (aend + tmp);
    } else {
      *s = (bend - tmp);
    }
    t++;s--;
  }
  if ( s == t ) {
    if ( *s < 0 ) {
      *s = - (aend + *s);
    } else {
      *s = (bend - *s);
    }
  }
  tmp =O->begpos;
  O->begpos = - O->endpos;
  O->endpos = - tmp;
  return 1;
}



static
ALNoverlap *
Compare(char *a, int32 alen,char *b, int32 blen, AS_ALN_Aligner *alignFunction, CNS_AlignParams *params) {
  ALNoverlap *O;

  int32 maxbegdef=MaxBegGap;
  int32 maxenddef=MaxEndGap;

  if ( params->bandBgn > alen)
    params->bandBgn = alen;

  if ( params->bandEnd > alen)
    params->bandEnd = alen;

  if ( params->bandEnd <-blen)
    params->bandEnd = -blen;

  if ( params->bandBgn <-blen)
    params->bandBgn = -blen;

  if ( params->erate > AS_MAX_ERROR_RATE)
    params->erate = AS_MAX_ERROR_RATE;

  MaxBegGap = params->maxBegGap;
  MaxEndGap = params->maxEndGap;

  O = (*alignFunction)(a, b,
                       params->bandBgn, params->bandEnd,
                       params->ahang,   params->bhang,
                       params->opposite,
                       params->erate,
                       params->thresh,
                       params->minlen,
                       params->what);

  MaxBegGap = maxbegdef;
  MaxEndGap = maxenddef;

#ifdef DEBUG_GET_ALIGNMENT_TRACE
  if (O)
    fprintf(stderr, "Compare()--  O: beg:%d end:%d len:%d diffs:%d comp:%d\n",
            O->begpos, O->endpos, O->length, O->diffs, O->comp);
  else
    fprintf(stderr, "Compare()--  No O.\n");
#endif

  return O;
}

static
void
ReportOverlap(FILE *fp, AS_ALN_Aligner *alignFunction, CNS_AlignParams params,
                   int32 aiid,char atype,int32 biid,char btype,ALNoverlap *O,int32 expected_hang) {

  if ((O == NULL) || (fp == NULL))
    return;

  if ( alignFunction == DP_Compare ) {
    fprintf(fp,"DP_Compare ");
  } else if (alignFunction == Local_Overlap_AS_forCNS ) {
    fprintf(fp,"Local_Overlap_AS_forCNS ");
  } else if (alignFunction == Optimal_Overlap_AS_forCNS ) {
    fprintf(fp,"Optimal_Overlap_AS_forCNS ");
  } else {
    fprintf(fp,"An alternate aligner ");
  }

  fprintf(fp,"found overlap between %d (%c) and %d (%c) ahang: %d, bhang: %d (expected hang was %d)\n",
          aiid,atype,biid,btype,O->begpos,O->endpos,expected_hang);
  fprintf(fp,"Alignment params: %d %d %d %d %d %5.2f %g %d %d\n", params.bandBgn, params.bandEnd,params.maxBegGap,params.maxEndGap,params.opposite,
          params.erate,params.thresh,params.minlen, params.what);
  if (O->begpos < 0 ) fprintf(fp,"Beware, encountered unexpected negative ahang!\n");
}

static
int
ScoreOverlap(ALNoverlap *O,
             int32         expected_length,
             int32         ahang_input,
             int32         bhang_input,
             double      maxerate,
             int32         alignment_context,
             double     *lScore_out,
             double     *aScore_out,
             double     *bScore_out) {

  if (O == NULL)
    return(0);

  if (((double)O->diffs / (double)O->length) > maxerate) {
    if (VERBOSE_MULTIALIGN_OUTPUT)
      fprintf(stderr,"GetAlignmentTrace()-- Overlap rejected.  erate %f > max allowed %f\n",
              (double)O->diffs / (double)O->length, maxerate);
    return(0);
  }

  double lScore = (O->length - expected_length) / (double)expected_length;
  double aScore = (O->begpos - ahang_input)     / (double)100.0;
  double bScore = (O->endpos - bhang_input)     / (double)100.0;

  if (lScore < 0)  lScore = -lScore;
  if (aScore < 0)  aScore = -aScore;
  if (bScore < 0)  bScore = -bScore;

  //  CGW merge doesn't compute the length very accurately.
  if (alignment_context == GETALIGNTRACE_MERGE) {
    if ((lScore < 0.01) &&
        ((aScore - bScore < 0.01) &&
         (bScore - aScore < 0.01))) {
      //  In merging, we occasionally find the correct alignment, but it seems to be shifted left or
      //  right.  This hack notices that quirk.
      lScore = 0.0;
      aScore = 0.0;
      bScore = 0.0;
    }

    //  But, in general, we want to be less strict for merges, trusting pretty much everything.
    lScore /= 2;
    aScore /= 5;
    bScore /= 5;
  }

  //  By design, we don't expect to have negative ahangs, so penalize the score if so (and if this
  //  is unexpected).
  //
  //  This same test is used above, at the end of a fragment overlap.
  //
  if ((O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF)) {
    fprintf(stderr,"GetAlignmentTrace()-- NEGATIVE HANG.  lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).\n",
            lScore, O->length, expected_length,
            aScore, O->begpos, ahang_input,
            bScore, O->endpos, bhang_input);
    aScore *= 10;
  }

  if (lScore_out) {
    *lScore_out = lScore;
    *aScore_out = aScore;
    *bScore_out = bScore;
  }

  //  Decide if these scores are good enough to accept the overlap.

  int32  isGood  = 0;
  int32  isGreat = 0;

  if (aScore < acceptThreshold)    isGood++;
  if (bScore < acceptThreshold)    isGood++;
  if (lScore < acceptThreshold)    isGood++;

  if (aScore < acceptThreshold/2)  isGreat++;
  if (bScore < acceptThreshold/2)  isGreat++;
  if (lScore < acceptThreshold/2)  isGreat++;


  switch (alignment_context) {
    case GETALIGNTRACE_UNITIG:
      if ((isGood >= 3) || (isGreat >= 2)) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr,"GetAlignmentTrace()-- Overlap ACCEPTED!  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).  (UNITIG)\n",
                  acceptThreshold,
                  lScore, O->length, expected_length,
                  aScore, O->begpos, ahang_input,
                  bScore, O->endpos, bhang_input);
        return(1);
      }
      break;

    case GETALIGNTRACE_CONTIGF:
      if (isGood >= 2) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr,"GetAlignmentTrace()-- Overlap ACCEPTED!  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).  (CONTIGF)\n",
                  acceptThreshold,
                  lScore, O->length, expected_length,
                  aScore, O->begpos, ahang_input,
                  bScore, O->endpos, bhang_input);
        return(1);
      }
      break;

    case GETALIGNTRACE_CONTIGU:
      if (isGood >= 2) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr,"GetAlignmentTrace()-- Overlap ACCEPTED!  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).  (CONTIGU)\n",
                  acceptThreshold,
                  lScore, O->length, expected_length,
                  aScore, O->begpos, ahang_input,
                  bScore, O->endpos, bhang_input);
        return(1);
      }
      break;

    case GETALIGNTRACE_MERGE:
      if ((isGood >= 2) || (isGreat >= 1)) {
        if (VERBOSE_MULTIALIGN_OUTPUT)
          fprintf(stderr,"GetAlignmentTrace()-- Overlap ACCEPTED!  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).  (MERGE)\n",
                  acceptThreshold,
                  lScore, O->length, expected_length,
                  aScore, O->begpos, ahang_input,
                  bScore, O->endpos, bhang_input);
        return(1);
      }
      break;
  }


  //  BAD!
  //ReportOverlap(stderr,alignFunction,params,afrag->iid,afrag->type,bfrag->iid,bfrag->type,O,ahang_input);
  //PrintALNoverlap(stderr, aseq, bseq, O);
  if (VERBOSE_MULTIALIGN_OUTPUT)
    fprintf(stderr,"GetAlignmentTrace()-- Overlap rejected.  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).\n",
            acceptThreshold,
            lScore, O->length, expected_length,
            aScore, O->begpos, ahang_input,
            bScore, O->endpos, bhang_input);

  return(0);
}

int
ScoreOverlap(ALNoverlap *O,
             int32         expected_length,
             int32         ahang_input,
             int32         bhang_input,
             double      maxerate,
             double     *lScore_out,
             double     *aScore_out,
             double     *bScore_out) {
   return ScoreOverlap(O, expected_length, ahang_input, bhang_input, maxerate, GETALIGNTRACE_MERGE, lScore_out, aScore_out, bScore_out);
}

//*********************************************************************************
// Look for the required overlap between two fragments, and return trace
//*********************************************************************************

int
GetAlignmentTrace(int32                        afid,
                  char                        *aseq_input,
                  int32                        bfid,
                  int32                       *ahang,
                  int32                       *bhang,
                  int32                        expected_length,
                  VA_TYPE(int32)              *trace,
                  OverlapType                 *otype,
                  AS_ALN_Aligner              *alignFunction,
                  int32                          show_olap,
                  int32                          allow_big_endgaps,
                  GetAlignmentTraceContext     alignment_context,
                  double                       input_erate) {
  int32   i;

  Fragment *afrag = NULL,          *bfrag = NULL;
  uint32    aiid  = 0,              biid  = 0;
  char     *aseq  = NULL,          *bseq  = NULL;
  char     *arev  = NULL,          *brev  = NULL;
  int32       alen  = 0,              blen  = 0;
  int32       ahang_input = *ahang,   bhang_input = *bhang;

  ALNoverlap        *O = NULL;

  CNS_AlignParams params;
  CNS_AlignParams paramsDefault;

  CNS_AlignParams_init(&paramsDefault);

  if (VERBOSE_MULTIALIGN_OUTPUT)
    show_olap = SHOW_OLAP;

  int32     CNS_TIGHTSEMIBANDWIDTH = 18;  //  Switch back to 6
  int32     CNS_LOOSESEMIBANDWIDTH = 100;
  int32     CNS_DP_THIN_MINLEN     = 10;

  //  Clear the 'output' early.  This is needed so we can force
  //  fragments in MultiAlignContig.
  //
  ResetVA_int32(trace);

  if (AS_CNS_ERROR_RATE > 0.06)
    CNS_TIGHTSEMIBANDWIDTH = 100;

  if (aseq_input) {
    assert(afid == -1);

    afrag = NULL;
    aseq  = aseq_input;
    aiid  = 0;
  } else {
    assert(afid >= 0);      //  If this triggers, see version 1.232 for what was here.

    afrag = GetFragment(fragmentStore,afid);
    aseq  = Getchar(sequenceStore,afrag->sequence);
    aiid  = afrag->iid;
  }

  bfrag = GetFragment(fragmentStore,bfid);
  bseq  = Getchar(sequenceStore,bfrag->sequence);
  biid  = bfrag->iid;

  alen  = strlen(aseq);
  blen  = strlen(bseq);

  //
  //  Set our default master parameters
  //
  //  The old logic was to:
  //
  //  atype is unitig or contig  ->  default_erate = 2 * input
  //  btype is unitig            ->  paramsDefault erate = 2 * CNS_ERATE
  //
  //  Then confusingly use either default_erate or paramsDefault.erate

  paramsDefault.maxBegGap = (allow_big_endgaps > 0) ? allow_big_endgaps : MaxBegGap;
  paramsDefault.maxEndGap = (allow_big_endgaps > 0) ? allow_big_endgaps : MaxEndGap;

  paramsDefault.bandBgn   = ahang_input - CNS_TIGHTSEMIBANDWIDTH;
  paramsDefault.bandEnd   = ahang_input + CNS_TIGHTSEMIBANDWIDTH;

  paramsDefault.ahang     = ahang_input;
  paramsDefault.bhang     = bhang_input;

  paramsDefault.erate     = input_erate;

  if (((afrag) && (afrag->type == AS_UNITIG)) ||
      ((bfrag) && (bfrag->type == AS_UNITIG)))
    paramsDefault.erate  *= 2;


  ////////////////////////////////////////
  //
  //  Compare with the default parameters.  Optimial_Overlap doesn't benefit from any
  //  of the tricks, so if that is the aligner, we're done after this attempt.
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 1\n");
#endif
  params = paramsDefault;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if ((O) || (alignFunction == Optimal_Overlap_AS_forCNS))
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // look for potentially narrower overlap
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 2\n");
#endif
  params        = paramsDefault;
  params.minlen = CNS_DP_THIN_MINLEN;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  //  If there are N's in the sequence, or they are consensus sequences, loosen the erate.
  //
  if ((bfrag->type==AS_UNITIG) || (strchr(aseq,'N') != NULL) || (strchr(bseq,'N') != NULL)) {
#ifdef DEBUG_GET_ALIGNMENT_TRACE
    fprintf(stderr, "GetAlignmentTrace()--  Compare 3\n");
#endif
    params       = paramsDefault;
    params.erate = paramsDefault.erate * 2;
    O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
    if ((O) && (allow_neg_hang == 0) && (O->begpos < 0))
      O = NULL;
    if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
      O = NULL;
    if (O)
      goto GetAlignmentTrace_ScoreOverlap;
  }


  ////////////////////////////////////////
  //
  //  The rest are for CONTIG and MERGE only.  No fragment-to-fragment alignments here.
  //
  ////////////////////////////////////////

  if (aseq_input)
    goto GetAlignmentTrace_ScoreOverlap;

  if ((alignment_context == GETALIGNTRACE_UNITIG) ||
      (alignment_context == GETALIGNTRACE_CONTIGF))
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // broaden scope out, and look for potentially narrower overlap
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 4\n");
#endif
  params         = paramsDefault;
  params.bandBgn = ahang_input-2*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = ahang_input+2*CNS_LOOSESEMIBANDWIDTH;
  params.erate   = paramsDefault.erate * 2;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // broaden even more, loosen erate
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 5\n");
#endif
  params         = paramsDefault;
  params.bandBgn = ahang_input-3*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = ahang_input+3*CNS_LOOSESEMIBANDWIDTH;
  params.erate   = paramsDefault.erate * 2;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // broaden even more, loosen erate
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 6\n");
#endif
  params         = paramsDefault;
  params.bandBgn = ahang_input-5*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = ahang_input+5*CNS_LOOSESEMIBANDWIDTH;
  params.erate   = paramsDefault.erate * 2;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // broaden even more, loosen erate
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 7\n");
#endif
  params         = paramsDefault;
  params.bandBgn = ahang_input-2*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = ahang_input+2*CNS_LOOSESEMIBANDWIDTH;
  params.erate   = paramsDefault.erate * 2;
  params.minlen  = CNS_DP_THIN_MINLEN;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // broaden even more, loosen erate
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 8\n");
#endif
  params         = paramsDefault;
  params.bandBgn = -blen;
  params.bandEnd = alen;
  params.erate   = paramsDefault.erate * 2;
  params.minlen  = CNS_DP_THIN_MINLEN;
  O = Compare(aseq,alen,bseq,blen,alignFunction,&params);
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  // perhaps a poor prefix is terminating the search, or causing an alternate
  // overlap to be found
  // try from other end
  //

  arev = (char *)safe_malloc(sizeof(char) * (alen + 1));
  brev = (char *)safe_malloc(sizeof(char) * (blen + 1));

  strcpy(arev, aseq);
  strcpy(brev, bseq);

  reverseComplementSequence(arev, alen);
  reverseComplementSequence(brev, blen);


  // calculate the hang if coming from the right instead
  //
  // note: the calc may be problematic: we really would like to have
  // the bhang, and the equation above gives exactly the bhang only if
  // the number of gaps in A and B is the same -- otherwise, we are
  // off by the number of gaps
  //


  ////////////////////////////////////////
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 9\n");
#endif
  params         = paramsDefault;
  params.bandBgn = alen - ahang_input - blen - CNS_TIGHTSEMIBANDWIDTH;
  params.bandEnd = alen - ahang_input - blen + CNS_TIGHTSEMIBANDWIDTH;
  O = Compare(arev,alen,brev,blen,alignFunction,&params);
  if (O) {
    InvertTrace(alen, blen, O);
  }
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 10\n");
#endif
  params         = paramsDefault;
  params.bandBgn = alen - ahang_input - blen - 2*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = alen - ahang_input - blen + 2*CNS_LOOSESEMIBANDWIDTH;
  params.erate   = paramsDefault.erate * 2;
  O = Compare(arev,alen,brev,blen,alignFunction,&params);
  if (O) {
    InvertTrace(alen, blen, O);
  }
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;


  ////////////////////////////////////////
  //
  //  try full length of fragments, due to troubles estimating the original bhang
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 11\n");
#endif
  params         = paramsDefault;
  params.bandBgn = -blen;
  params.bandEnd = alen;
  params.erate   = paramsDefault.erate * 2;
  O = Compare(arev,alen,brev,blen,alignFunction,&params);
  if (O) {
    InvertTrace(alen, blen, O);
  }
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;
    

  ////////////////////////////////////////
  //
  // here, we'll try to swap the fragments too
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 12\n");
#endif
  params         = paramsDefault;
  params.bandBgn = -alen - ahang_input - blen - 2*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = -alen - ahang_input - blen + 2*CNS_LOOSESEMIBANDWIDTH;
  params.ahang   = -params.ahang;
  params.bhang   = -params.bhang;
  O = Compare(brev,blen,arev,alen,alignFunction,&params);
  if (O) {
    for (i=0; O->trace[i] != 0; i++)
      O->trace[i] = -O->trace[i];
    O->begpos = -O->begpos;
    O->endpos = -O->endpos;
    InvertTrace(alen,blen,O);
  }
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;



  ////////////////////////////////////////
  //
  // this still isn't a good overlap
  // try to see whether just swapping the fragments to see if that locates the overlap
  //
#ifdef DEBUG_GET_ALIGNMENT_TRACE
  fprintf(stderr, "GetAlignmentTrace()--  Compare 13\n");
#endif

  params         = paramsDefault;
  params.bandBgn = ahang_input-3*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = ahang_input+3*CNS_LOOSESEMIBANDWIDTH;
  params.bandEnd = alen-CNS_DP_MINLEN;

  i              = params.bandBgn;
  params.bandBgn = -params.bandEnd;
  params.bandEnd = -i;

  params.ahang   = -params.ahang;
  params.bhang   = -params.bhang;

  O = Compare(bseq,blen,aseq,alen,alignFunction,&params);
  if (O) {
    for (i=0; O->trace[i] != 0; i++)
      O->trace[i] = -O->trace[i];
    O->begpos = -O->begpos;
    O->endpos = -O->endpos;
  }
  if ((O) && (allow_neg_hang == 0) && (O->begpos < 0) && (0 < ahang_input + CNS_NEG_AHANG_CUTOFF))
    O = NULL;
  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, NULL, NULL, NULL) == 0)
    O = NULL;
  if (O)
    goto GetAlignmentTrace_ScoreOverlap;

  //
  //  End of MERGE or UNITIG tricks
  //

 GetAlignmentTrace_ScoreOverlap:

  safe_free(arev);
  safe_free(brev);

  {
    const char *aligner = "(something else)";
    if (alignFunction == DP_Compare)                aligner = "DP_Compare";
    if (alignFunction == Local_Overlap_AS_forCNS)   aligner = "Local_Overlap";
    if (alignFunction == Optimal_Overlap_AS_forCNS) aligner = "Optimal_Overlap";

    if ( O == NULL )
      return(FALSE);

    if (show_olap)
      fprintf(stderr,"GetAlignmentTrace()-- Overlap found between %d (%c) and %d (%c) expected hangs: a=%d b=%d erate=%f aligner=%s\n",
              aiid, (afrag) ? afrag->type : AS_CONTIG, biid, bfrag->type, ahang_input, bhang_input, input_erate, aligner);
  }

  //
  //  From this point32 on, we have an Overlap.  Score it and decide if
  //  it is the one we are looking for.
  //

  double  lScore = 0.0;
  double  aScore = 0.0;
  double  bScore = 0.0;

  if (ScoreOverlap(O, expected_length, ahang_input, bhang_input, params.erate, alignment_context, &lScore, &aScore, &bScore) == 0) {
    //  Bad.
    //
    //  Should never occur as we throw out bad overlaps as soon as we generate them.
    //
    assert(0);

    ReportOverlap(stderr, alignFunction, params, aiid, (afrag) ? afrag->type : AS_CONTIG, biid, bfrag->type, O, ahang_input);
    PrintALNoverlap("", aseq, bseq, O);
    fprintf(stderr,"GetAlignmentTrace()-- Overlap rejected.  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).\n",
            acceptThreshold,
            lScore, O->length, expected_length,
            aScore, O->begpos, ahang_input,
            bScore, O->endpos, bhang_input);

    return(FALSE);
  }

  if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ACCEPTED)
    if ((lScore > 0.1) || (aScore > 0.1) || (bScore > 0.1))
      fprintf(stderr,"GetAlignmentTrace()-- Overlap accepted.  accept=%f lScore=%f (%d vs %d) aScore=%f (%d vs %d) bScore=%f (%d vs %d).\n",
              acceptThreshold,
              lScore, O->length, expected_length,
              aScore, O->begpos, ahang_input,
              bScore, O->endpos, bhang_input);

  //
  //  From this point32 on, we have a Good Overlap.
  //

  lScoreAve += lScore;
  aScoreAve += aScore;
  bScoreAve += bScore;
  numScores++;

  if (show_olap) {
    ReportOverlap(stderr, alignFunction, params, aiid, (afrag) ? afrag->type : AS_CONTIG, biid, bfrag->type, O, ahang_input);
    PrintALNoverlap("", NULL, NULL, O);  //  Replace with a,b to print the bases in the align
  }

  *otype = (O->endpos<0)?AS_CONTAINMENT:AS_DOVETAIL;
  *ahang = O->begpos;
  *bhang = O->endpos;

  {
    int32 *tmp = O->trace;

    ResetVA_int32(trace);

    while ( *tmp != 0) {
      AppendVA_int32(trace,tmp);
      tmp++;
    }
  }

  return(TRUE);
}


int
GetAlignmentTraceDriver(Fragment                     *afrag,
                        char                         *aseq,
                        Fragment                     *bfrag,
                        int32                        *ahang,
                        int32                        *bhang,
                        int32                         expected_length,
                        VA_TYPE(int32)               *trace,
                        OverlapType                  *otype,
                        GetAlignmentTraceContext      alignment_context,
                        int32                           max_gap) {
  double AS_CNS_ERROR_RATE_SAVE = AS_CNS_ERROR_RATE;
  int32  aiid = (afrag) ? afrag->iid  : 0;
  int32  alid = (afrag) ? afrag->lid  : -1;
  char   atyp = (afrag) ? afrag->type : AS_CONTIG;
  int32  biid = (bfrag) ? bfrag->iid  : 0;
  int32  blid = (bfrag) ? bfrag->lid  : -1;
  char   btyp = (bfrag) ? bfrag->type : AS_CONTIG;

  for (; (AS_CNS_ERROR_RATE < AS_CNS_ERROR_RATE_SAVE + 0.03); AS_CNS_ERROR_RATE += 0.0025) {

    //  try with the default aligner

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
      fprintf(stderr, "GetAlignmentTraceDriver()-- Attempting alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d and erate %1.4f (DP_Compare)\n",
              aiid, atyp, biid, btyp,
              *ahang,
              AS_CNS_ERROR_RATE);
    if (GetAlignmentTrace(alid, aseq, blid, ahang, bhang, expected_length, trace, otype, DP_Compare, DONT_SHOW_OLAP, 0, alignment_context, AS_CNS_ERROR_RATE)) {
      AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE_SAVE;
      return(TRUE);
    }

    // try again, perhaps with alternate overlapper

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
      fprintf(stderr, "GetAlignmentTraceDriver()-- Attempting alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d erate %1.4f (Local_Overlap)\n",
              aiid, atyp, biid, btyp,
              *ahang,
              AS_CNS_ERROR_RATE);
    if (GetAlignmentTrace(alid, aseq, blid, ahang, bhang, expected_length, trace, otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, alignment_context, AS_CNS_ERROR_RATE)) {
      AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE_SAVE;
      return(TRUE);
    }

    // try again, perhaps with dynamic programming, but only for unitigs (extend once we fix up MultiAlignContig)

    if (alignment_context == GETALIGNTRACE_UNITIG) {
      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
        fprintf(stderr, "GetAlignmentTraceDriver()-- Attempting alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d erate %1.4f (Optimal_Overlap)\n",
                aiid, atyp, biid, btyp,
                *ahang,
                AS_CNS_ERROR_RATE);
      if (GetAlignmentTrace(alid, aseq, blid, ahang, bhang, expected_length, trace, otype, Optimal_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, alignment_context, AS_CNS_ERROR_RATE)) {
        AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE_SAVE;
        return(TRUE);
      }
    }

    //  try again, perhaps allowing larger end gaps (contigs only)

    if ((alignment_context == GETALIGNTRACE_CONTIGU) && (max_gap > 0)) {
      if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
        fprintf(stderr, "GetAlignmentTraceDriver()-- Attempting alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d erate %1.4f max_gap %d (Local_Overlap)\n",
                aiid, atyp, biid, btyp,
                *ahang,
                AS_CNS_ERROR_RATE,
                max_gap);
      if (GetAlignmentTrace(alid, aseq,blid, ahang, bhang, expected_length, trace, otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, max_gap, alignment_context, AS_CNS_ERROR_RATE)) {
        AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE_SAVE;
        return(TRUE);
      }
    }
  }

  AS_CNS_ERROR_RATE = AS_CNS_ERROR_RATE_SAVE;

  //  If all those failed, and we are not consensus, fail.

  if (thisIsConsensus == 0)
    return(FALSE);

  //  Otherwise, try again, but greatly relax the length, ahang and bhang correctness criteria.
  //  This is a last ditch effort to get an alignment, used only in consensus (not cgw).
  //
  //  The attepts are always reported, since this is an exceptional event.
  //
  if ((alignment_context == GETALIGNTRACE_CONTIGU) ||
      (alignment_context == GETALIGNTRACE_CONTIGF)) {
    int32     oldVB = VERBOSE_MULTIALIGN_OUTPUT;
    double  oldAT = acceptThreshold;

    VERBOSE_MULTIALIGN_OUTPUT = MAX(oldVB, 1);
    acceptThreshold           = 1000.0;

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
      fprintf(stderr, "GetAlignmentTraceDriver()-- Attemping alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d erate %1.4f max_gap %d (DP_Compare) LAST DITCH\n",
              aiid, atyp, biid, btyp,
              *ahang,
              AS_CNS_ERROR_RATE,
              max_gap);
    if (GetAlignmentTrace(alid, aseq, blid, ahang, bhang, expected_length, trace, otype, DP_Compare, DONT_SHOW_OLAP, 0, alignment_context, AS_CNS_ERROR_RATE)) {
      VERBOSE_MULTIALIGN_OUTPUT = oldVB;
      acceptThreshold           = oldAT;
      return(TRUE);
    }

    if (VERBOSE_MULTIALIGN_OUTPUT >= SHOW_ATTEMPT)
      fprintf(stderr, "GetAlignmentTraceDriver()-- Attemping alignment of afrag %d (%c) and bfrag %d (%c) with ahang %d erate %1.4f max_gap %d (Local_Overlap) LAST DITCH\n",
              aiid, atyp, biid, btyp,
              *ahang,
              AS_CNS_ERROR_RATE,
              max_gap);
    if (GetAlignmentTrace(alid, aseq, blid, ahang, bhang, expected_length, trace, otype, Local_Overlap_AS_forCNS, DONT_SHOW_OLAP, 0, alignment_context, AS_CNS_ERROR_RATE)) {
      VERBOSE_MULTIALIGN_OUTPUT = oldVB;
      acceptThreshold           = oldAT;
      return(TRUE);
    }

    VERBOSE_MULTIALIGN_OUTPUT = oldVB;
    acceptThreshold           = oldAT;
  }

  return(FALSE);
}
