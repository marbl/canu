
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
 *    Brian P. Walenz from 2015-JUN-18 to 2015-JUL-01
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "analyzeAlignment.H"



//  Return the substitution vote corresponding to  Ch .
static
Vote_Value_t
Matching_Vote(char ch) {

  switch  (ch) {
    case 'A':  return(A_SUBST);  break;
    case 'C':  return(C_SUBST);  break;
    case 'G':  return(G_SUBST);  break;
    case 'T':  return(T_SUBST);  break;
  }

  fprintf(stderr, "Matching_Vote()-- invalid letter '%c'\n", ch);

  return(NO_VOTE);
}


static
char
Matching_Char(Vote_Value_t vv) {

  switch  (vv) {
    case A_SUBST:  return('A');  break;
    case C_SUBST:  return('C');  break;
    case G_SUBST:  return('G');  break;
    case T_SUBST:  return('T');  break;
    default:       return('?');  break;
  }

  return('?');
}




//  This is expecting:
//    aSeq and bSeq to be pointers to the start of the sequence that was aligned with Prefix_Edit_Distance
//    aLen and bLen to be the strlen of those strings
//    aOffset to be the offset from base zero of the read
//    deltaLen and delta as passed back from the aligner
//
void
analyzeAlignment::analyze(char  *aSeq, int32 aLen,  int32 aOffset,
                          char  *bSeq, int32 bLen,
                          int32  deltaLen,
                          int32 *delta) {

  assert(aLen >= 0);
  assert(bLen >= 0);

  int32  ct = 0;

  _readSub[ct]   = -1;
  _algnSub[ct]   = -1;
  _voteValue[ct] = A_SUBST;   // Dummy value
  ct++;

  int32  i = 0;
  int32  j = 0;
  int32  p = 0;

  uint32  nMatch    = 0;
  uint32  nMismatch = 0;
  uint32  nInsert   = 0;
  uint32  nDelete   = 0;

  for (int32 k=0; k<deltaLen; k++) {
#ifdef DEBUG
    fprintf(stderr, "#0\n");
    fprintf(stderr, "#1  k=%d deltalen=%d delta=%d  i=%d out of %d   j=%d out of %d\n",
            k, deltaLen, delta[k], i, aLen, j, bLen);
#endif

    //  Add delta[k] matches or mismatches

    for (int32 m=1; m<abs(delta[k]); m++) {
      if (aSeq[i] != bSeq[j])
        nMismatch++;
      else
        nMatch++;

      if (aSeq[i] != bSeq[j]) {
        _readSub[ct] = i;
        _algnSub[ct] = p;

        switch (bSeq[j]) {
          case 'A':  _voteValue[ct] = A_SUBST;  break;
          case 'C':  _voteValue[ct] = C_SUBST;  break;
          case 'G':  _voteValue[ct] = G_SUBST;  break;
          case 'T':  _voteValue[ct] = T_SUBST;  break;
          default :
            fprintf(stderr, "ERROR:[1] Bad sequence '%c' 0x%02x)\n", bSeq[j], bSeq[j]);
            assert(0);
        }

        ct++;
      }

      i++;  assert(i <= aLen);
      j++;  assert(j <= bLen);
      p++;
    }

#ifdef DEBUG
    fprintf(stderr, "#2  match %u mismatch %u insert %u delete %u\n", nMatch, nMismatch, nInsert, nDelete);
#endif

    //  If a negative delta, insert a base.

    if (delta[k] < 0) {
      _readSub[ct] = i - 1;
      _algnSub[ct] = p;

#ifdef DEBUG
      fprintf(stderr, "INSERT %c at %d #%d\n", bSeq[j], i-1, p);
#endif
      nInsert++;

      switch (bSeq[j]) {
        case 'A':  _voteValue[ct] = A_INSERT;  break;
        case 'C':  _voteValue[ct] = C_INSERT;  break;
        case 'G':  _voteValue[ct] = G_INSERT;  break;
        case 'T':  _voteValue[ct] = T_INSERT;  break;
        default :
          fprintf(stderr, "ERROR:[2] Bad sequence '%c' 0x%02x)\n", bSeq[j], bSeq[j]);
          assert(0);
      }

      ct++;

      j++;  assert(j <= bLen);
      p++;
    }

#ifdef DEBUG
    fprintf(stderr, "#3  match %u mismatch %u insert %u delete %u\n", nMatch, nMismatch, nInsert, nDelete);
#endif

    //  If a positive deta, delete the base.

    if (delta[k] > 0) {
      _readSub[ct] = i;
      _algnSub[ct] = p;
      _voteValue[ct]  = DELETE;

#ifdef DEBUG
      fprintf(stderr, "DELETE %c at %d #%d\n", aSeq[i], i, p);
#endif
      nDelete++;

      ct++;

      i++;  assert(i <= aLen);
      p++;
    }

#ifdef DEBUG
    fprintf(stderr, "#4  match %u mismatch %u insert %u delete %u\n", nMatch, nMismatch, nInsert, nDelete);
#endif
  }

  // No more deltas.  While there is still sequence, add matches or mismatches.

#ifdef DEBUG
  fprintf(stderr, "#5  match %u mismatch %u insert %u delete %u\n", nMatch, nMismatch, nInsert, nDelete);
  fprintf(stderr, "#5  k=DONE   i=%d out of %d   j=%d out of %d\n", i, aLen, j, bLen);
#endif

  while (i < aLen) {
    //fprintf(stderr, "k=DONE   i=%d out of %d   j=%d out of %d\n", i, aLen, j, bLen);

    if (aSeq[i] != bSeq[j])
      nMismatch++;
    else
      nMatch++;

    if (aSeq[i] != bSeq[j]) {
      _readSub[ct] = i;
      _algnSub[ct] = p;

      switch (bSeq[j]) {
        case 'A':  _voteValue[ct] = A_SUBST;  break;
        case 'C':  _voteValue[ct] = C_SUBST;  break;
        case 'G':  _voteValue[ct] = G_SUBST;  break;
        case 'T':  _voteValue[ct] = T_SUBST;  break;
        default :
          fprintf(stderr, "ERROR:[3] Bad sequence '%c' 0x%02x)\n", bSeq[j], bSeq[j]);
          assert(0);
      }

      ct++;
    }

    i++;  assert(i <= aLen);  //  Guaranteed, we're looping on this
    j++;  assert(j <= bLen);
    p++;
  }

#ifdef DEBUG
  fprintf(stderr, "#6  match %u mismatch %u insert %u delete %u\n", nMatch, nMismatch, nInsert, nDelete);
#endif

  _readSub[ct] = i;
  _algnSub[ct] = p;


  //  For each identified change, add votes for some region around the change.

  //fprintf(stderr, "Found %u changes.\n", ct);

  for (int32 i=1; i<=ct; i++) {
    int32  prev_match = _algnSub[i] - _algnSub[i-1] - 1;
    int32  p_lo = (i == 1 ? 0 : End_Exclude_Len);
    int32  p_hi = (i == ct ? prev_match : prev_match - End_Exclude_Len);

    //  If distance to previous match is bigger than 'kmer' size, make a new vote.
    //  This operates one ahead of where votes are added - we add votes for _readSub[i-1] when at [i].

    if (prev_match >= Kmer_Len) {
      fprintf(stderr, "adjust ct %d pos %d - lo %d hi %d\n", i, _readSub[i-1], p_lo, p_hi);

      fprintf(stderr, "  match vote %u to %u\n", aOffset + _readSub[i-1] + 1, aOffset + _readSub[i-1] + p_lo + 1);
      for (int32 p=0;  p<p_lo;  p++) {
        castVote(Matching_Vote(aSeq[_readSub[i-1] + p + 1]), aOffset + _readSub[i-1] + p + 1);
      }

      fprintf(stderr, "  no insert %u to %u\n", aOffset + _readSub[i-1] + p_lo + 1, aOffset + _readSub[i-1] + p_hi + 1);
      for (int32 p=p_lo;  p<p_hi;  p++) {
        int32 k = aOffset + _readSub[i-1] + p + 1;

        if (_vote[k].confirmed < VOTETALLY_MAX)
          _vote[k].confirmed++;

        if ((p < p_hi - 1) &&
            (_vote[k].no_insert < VOTETALLY_MAX))
          _vote[k].no_insert++;
      }

      fprintf(stderr, "  match vote %u to %u\n", aOffset + _readSub[i-1] + p_hi + 1, aOffset + _readSub[i-1] + prev_match + 1);
      for (int32 p=p_hi; p<prev_match; p++) {
        castVote(Matching_Vote(aSeq[_readSub[i-1] + p + 1]), aOffset + _readSub[i-1] + p + 1);
      }
    }

    //  Don't allow consecutive inserts.  If we aren't the last change, and there is non-adjacent
    //  previous (or this and the previous votes are not insertions), do another vote.

    if ((i < ct) &&
        ((prev_match > 0) ||
         (_voteValue[i-1] <= T_SUBST) ||
         (_voteValue[i]   <= T_SUBST))) {
      int32 next_match = _algnSub[i+1] - _algnSub[i] - 1;

      if (prev_match + next_match >= Vote_Qualify_Len)
        castVote(_voteValue[i], aOffset + _readSub[i]);
    }
  }
}




void
analyzeAlignment::outputDetails(uint32 j) {
  fprintf(stderr, "%3"F_U32P": %c  conf %3"F_U64P"  deletes %3"F_U64P" | subst %3"F_U64P" %3"F_U64P" %3"F_U64P" %3"F_U64P" | no_insert %3"F_U64P" insert %3"F_U64P" %3"F_U64P" %3"F_U64P" %3"F_U64P"\n",
          j,
          _seq[j],
          _vote[j].confirmed,
          _vote[j].deletes,
          _vote[j].a_subst,
          _vote[j].c_subst,
          _vote[j].g_subst,
          _vote[j].t_subst,
          _vote[j].no_insert,
          _vote[j].a_insert,
          _vote[j].c_insert,
          _vote[j].g_insert,
          _vote[j].t_insert);
};


void
analyzeAlignment::outputDetails(void) {
  fprintf(stderr, ">%d\n", _readID);
  for (uint32 j=0; _seq[j] != '\0'; j++)
    outputDetails(j);
};




void
analyzeAlignment::generateCorrections(FILE *corFile) {

  _corLen = 0;

  _cor[_corLen].keep_left   = (_lDegree < Degree_Threshold);
  _cor[_corLen].keep_right  = (_rDegree < Degree_Threshold);
  _cor[_corLen].type        = IDENT;
  _cor[_corLen].pos         = 0;
  _cor[_corLen].readID      = _readID;

  _corLen++;
  resizeArray(_cor, _corLen, _corMax, _corLen+1);

  uint32   passedLowConfirmed = 0;
  uint32   substitutions      = 0;
  uint32   skippedTooFew      = 0;   //  0 or 1 votes
  uint32   skippedTooWeak     = 0;   //  No vote more than 50%
  uint32   skippedNoChange    = 0;   //  Is a substitution vote, but it's the same as the base that is there
  uint32   skippedHaplo       = 0;   //  More than one significant vote, and we're not correcting haplotypes
  uint32   skippedConfirmed   = 0;   //  ??


  uint32   passedInsert       = 0;
  uint32   insertions         = 0;
  uint32   skippedInsTotal    = 0;
  uint32   skippedInsMax      = 0;
  uint32   skippedInsHaplo    = 0;
  uint32   skippedInsTooMany  = 0;

  for (uint32 j=0; j<_seqLen; j++) {

    outputDetails(j);

    if  (_vote[j].confirmed < 2) {
      Vote_Value_t  vval      = DELETE;
      uint64        max       = _vote[j].deletes;
      bool          is_change = true;

      passedLowConfirmed++;

      if  (_vote[j].a_subst > max) {
        vval      = A_SUBST;
        max       = _vote[j].a_subst;
        is_change = (_seq[j] != 'A');
      }

      if  (_vote[j].c_subst > max) {
        vval      = C_SUBST;
        max       = _vote[j].c_subst;
        is_change = (_seq[j] != 'C');
      }

      if  (_vote[j].g_subst > max) {
        vval      = G_SUBST;
        max       = _vote[j].g_subst;
        is_change = (_seq[j] != 'G');
      }

      if  (_vote[j].t_subst > max) {
        vval      = T_SUBST;
        max       = _vote[j].t_subst;
        is_change = (_seq[j] != 'T');
      }

      uint64 haplo_ct  =  ((_vote[j].deletes >= Min_Haplo_Occurs) +
                           (_vote[j].a_subst >= Min_Haplo_Occurs) +
                           (_vote[j].c_subst >= Min_Haplo_Occurs) +
                           (_vote[j].g_subst >= Min_Haplo_Occurs) +
                           (_vote[j].t_subst >= Min_Haplo_Occurs));

      uint64 total  = (_vote[j].deletes +
                       _vote[j].a_subst +
                       _vote[j].c_subst +
                       _vote[j].g_subst +
                       _vote[j].t_subst);

      //  The original had a gargantuan if test (five clauses, all had to be true) to decide if a
      //  record should be output.  It was negated into many small tests if we should skip the
      //  output.  A side effect is that we can abort a little earlier (skipping the two clauses
      //  above)....but we don't bother.

      //  (total > 1)
      if (total <= 1) {
        fprintf(stderr, "FEW   total = "F_U64" <= 1\n", total);
        skippedTooFew++;
        continue;
      }

      //  (2 * max > total)
      if (2 * max <= total) {
        fprintf(stderr, "WEAK  2*max = "F_U64" <= total = "F_U64"\n", 2*max, total);
        skippedTooWeak++;
        continue;
      }

      //  (is_change == true)
      if (is_change == false) {
        fprintf(stderr, "SAME  is_change = %s\n", (is_change) ? "true" : "false");
        skippedNoChange++;
        continue;
      }

      //  ((haplo_ct < 2) || (Use_Haplo_Ct == false))
      if ((haplo_ct >= 2) && (Use_Haplo_Ct == true)) {
        fprintf(stderr, "HAPLO haplo_ct="F_U64" >= 2 AND Use_Haplo_Ct = %s\n", haplo_ct, (Use_Haplo_Ct) ? "true" : "false");
        skippedHaplo++;
        continue;
      }

      //  ((_vote[j].confirmed == 0) ||
      //   ((_vote[j].confirmed == 1) && (max > 6)))
      if ((_vote[j].confirmed > 0) &&
          ((_vote[j].confirmed != 1) || (max <= 6))) {
        fprintf(stderr, "INDET confirmed = "F_U64" max = "F_U64"\n", _vote[j].confirmed, max);
        skippedConfirmed++;
        continue;
      }

      //  Otherwise, output.

      substitutions++;

      fprintf(stderr, "SUBSTITUTE position "F_U32" to %c\n", j, Matching_Char(vval));

      _cor[_corLen].type       = vval;
      _cor[_corLen].pos        = j;
      _cor[_corLen].readID     = _readID;

      _corLen++;
      resizeArray(_cor, _corLen, _corMax, _corLen+1);
    }  //  confirmed < 2


    if  (_vote[j].no_insert < 2) {
      Vote_Value_t  ins_vote = A_INSERT;
      uint64        ins_max  = _vote[j].a_insert;

      passedInsert++;

      if  (ins_max < _vote[j].c_insert) {
        ins_vote = C_INSERT;
        ins_max  = _vote[j].c_insert;
      }

      if  (ins_max < _vote[j].g_insert) {
        ins_vote = G_INSERT;
        ins_max  = _vote[j].g_insert;
      }

      if  (ins_max < _vote[j].t_insert) {
        ins_vote = T_INSERT;
        ins_max  = _vote[j].t_insert;
      }

      uint64 ins_haplo_ct = ((_vote[j].a_insert >= Min_Haplo_Occurs) +
                             (_vote[j].c_insert >= Min_Haplo_Occurs) +
                             (_vote[j].g_insert >= Min_Haplo_Occurs) +
                             (_vote[j].t_insert >= Min_Haplo_Occurs));

      uint64 ins_total = (_vote[j].a_insert +
                          _vote[j].c_insert +
                          _vote[j].g_insert +
                          _vote[j].t_insert);

      if (ins_total <= 1) {
        fprintf(stderr, "FEW   ins_total = "F_U64" <= 1\n", ins_total);
        skippedInsTotal++;
        continue;
      }

      if (2 * ins_max >= ins_total) {
        fprintf(stderr, "WEAK  2*ins_max = "F_U64" <= ins_total = "F_U64"\n", 2*ins_max, ins_total);
        skippedInsMax++;
        continue;
      }

      if ((ins_haplo_ct >= 2) && (Use_Haplo_Ct == true)) {
        fprintf(stderr, "HAPLO ins_haplo_ct="F_U64" >= 2 AND Use_Haplo_Ct = %s\n", ins_haplo_ct, (Use_Haplo_Ct) ? "true" : "false");
        skippedInsHaplo++;
        continue;
      }

      if ((_vote[j].no_insert > 0) &&
          ((_vote[j].no_insert != 1) || (ins_max <= 6))) {
        fprintf(stderr, "INDET no_insert = "F_U64" ins_max = "F_U64"\n", _vote[j].no_insert, ins_max);
        skippedInsTooMany++;
        continue;
      }

      //  Otherwise, output.

      insertions++;

      fprintf(stderr, "INSERT position "F_U32" to %c\n", j, Matching_Char(ins_vote));

      _cor[_corLen].type       = ins_vote;
      _cor[_corLen].pos        = j;
      _cor[_corLen].readID     = _readID;

      _corLen++;
      resizeArray(_cor, _corLen, _corMax, _corLen+1);
    }  //  insert < 2
  }

  fprintf(stderr, "Processed corrections: made %6u subs and %6u inserts - possible %6u (few %6u weak %6u same %6u haplo %6u confirmed %6u) inserts %6u (total %6u max %6u haplo %6u confirmed %6u)\n",
          substitutions,
          insertions,
          passedLowConfirmed,
          skippedTooFew,
          skippedTooWeak,
          skippedNoChange,
          skippedHaplo,
          skippedConfirmed,
          passedInsert,
          skippedInsTotal,
          skippedInsMax,
          skippedInsHaplo,
          skippedInsTooMany);

  if (corFile)
    AS_UTL_safeWrite(corFile, _cor, "corrections", sizeof(Correction_Output_t), _corLen);
}




#if 0
  fprintf(stderr, "Corrected "F_U64" bases with "F_U64" substitutions, "F_U64" deletions and "F_U64" insertions.\n",
          G->basesLen,
          changes[A_SUBST] + changes[C_SUBST] + changes[G_SUBST] + changes[T_SUBST],
          changes[DELETE],
          changes[A_INSERT] + changes[C_INSERT] + changes[G_INSERT] + changes[T_INSERT]);
#endif


void
analyzeAlignment::generateCorrectedRead(Adjust_t *fadj, uint32 *fadjLen,
                                        uint64   *changes) {

  //  oseq = original sequence
  //  fseq - fixed sequence

  _corSeqLen = 0;

  uint32  corPos = 1;        //  First entry is the ident block.
  uint32  adjVal = 0;        //  Corrected reads start at position zero, really!

  for (uint32 i=0; i<_seqLen; i++) {

    //  No more corrections, or no more corrections for this read -- just copy bases till the end.

    if ((corPos == _corLen) || (_cor[corPos].readID != _readID)) {
      //fprintf(stderr, "no more corrections at i=%u, copy rest of read as is\n", i);
      while (i < _seqLen)
        _corSeq[_corSeqLen++] = _filter[_seq[i++]];
      break;
    }

    //  Not at a correction -- copy the base.
    if (i < _cor[corPos].pos) {
      _corSeq[_corSeqLen++] = _filter[_seq[i]];
      continue;
    }

    if ((i != _cor[corPos].pos) &&
        (i != _cor[corPos].pos + 1))
      fprintf(stderr, "i=%d corPos=%d _cor[corPos].pos=%d\n", i, corPos, _cor[corPos].pos);
    assert((i == _cor[corPos].pos) ||
           (i == _cor[corPos].pos + 1));

    if (changes)
      changes[_cor[corPos].type]++;

    switch (_cor[corPos].type) {
      case DELETE:  //  Delete base
        //fprintf(stderr, "DELETE %u pos %u adjust %d\n", (*fadjLen), i+1, adjVal-1);
        if (fadj) {
          fadj[(*fadjLen)].adjpos = i + 1;
          fadj[(*fadjLen)].adjust = --adjVal;
          (*fadjLen)++;
        }
        break;

      case A_SUBST:  _corSeq[_corSeqLen++] = 'A';  break;
      case C_SUBST:  _corSeq[_corSeqLen++] = 'C';  break;
      case G_SUBST:  _corSeq[_corSeqLen++] = 'G';  break;
      case T_SUBST:  _corSeq[_corSeqLen++] = 'T';  break;

      case A_INSERT:
        if (i != _cor[corPos].pos + 1) {                // Insert not immediately after subst
          //fprintf(stderr, "A i=%d != _cor[%d].pos+1=%d\n", i, corPos, _cor[corPos].pos+1);
          _corSeq[_corSeqLen++] = _filter[_seq[i++]];
        }
        _corSeq[_corSeqLen++] = 'A';

        if (fadj) {
          fadj[(*fadjLen)].adjpos = i + 1;
          fadj[(*fadjLen)].adjust = ++adjVal;
          (*fadjLen)++;
        }
        i--;  //  Undo the automagic loop increment
        break;

      case C_INSERT:
        if (i != _cor[corPos].pos + 1) {
          //fprintf(stderr, "C i=%d != _cor[%d].pos+1=%d\n", i, corPos, _cor[corPos].pos+1);
          _corSeq[_corSeqLen++] = _filter[_seq[i++]];
        }
        _corSeq[_corSeqLen++] = 'C';

        if (fadj) {
          fadj[(*fadjLen)].adjpos = i + 1;
          fadj[(*fadjLen)].adjust = ++adjVal;
          (*fadjLen)++;
        }
        i--;
        break;

      case G_INSERT:
        if (i != _cor[corPos].pos + 1) {
          //fprintf(stderr, "G i=%d != _cor[%d].pos+1=%d\n", i, corPos, _cor[corPos].pos+1);
          _corSeq[_corSeqLen++] = _filter[_seq[i++]];
        }
        _corSeq[_corSeqLen++] = 'G';

        if (fadj) {
          fadj[(*fadjLen)].adjpos = i + 1;
          fadj[(*fadjLen)].adjust = ++adjVal;
          (*fadjLen)++;
        }
        i--;
        break;

      case T_INSERT:
        if (i != _cor[corPos].pos + 1) {
          //fprintf(stderr, "T i=%d != _cor[%d].pos+1=%d\n", i, corPos, _cor[corPos].pos+1);
          _corSeq[_corSeqLen++] = _filter[_seq[i++]];
        }
        _corSeq[_corSeqLen++] = 'T';

        if (fadj) {
          fadj[(*fadjLen)].adjpos = i + 1;
          fadj[(*fadjLen)].adjust = ++adjVal;
          (*fadjLen)++;
        }
        i--;
        break;

      default:
        fprintf (stderr, "ERROR:  Illegal vote type\n");
        break;
    }

    corPos++;
  }

  //  Terminate the sequence.

  _corSeq[_corSeqLen] = 0;
}

