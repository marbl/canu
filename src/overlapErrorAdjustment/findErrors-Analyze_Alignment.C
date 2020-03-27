
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

#include "findErrors.H"
#include <string>
#include <vector>

//  Add vote val to G.reads[sub] at sequence position  p
static void
Cast_Vote(feParameters *G,
          Vote_Value_t val,
          int32        pos,
          int32        sub) {
  Vote_Tally_t &vote = G->reads[sub].vote[pos];
  //fprintf(stderr, "Casting vote val %d at pos %d\n", val, pos);
  switch (val) {
    case DELETE:
      //fprintf(stderr, "Casting deletion\n");
      if (vote.deletes < MAX_VOTE)
        vote.deletes++;
      break;
    case A_SUBST:
      //fprintf(stderr, "Casting A_SUBST\n");
      if (vote.a_subst < MAX_VOTE)
        vote.a_subst++;
      break;
    case C_SUBST:
      //fprintf(stderr, "Casting C_SUBST\n");
      if (vote.c_subst < MAX_VOTE)
        vote.c_subst++;
      break;
    case G_SUBST:
      //fprintf(stderr, "Casting G_SUBST\n");
      if (vote.g_subst < MAX_VOTE)
        vote.g_subst++;
      break;
    case T_SUBST:
      //fprintf(stderr, "Casting T_SUBST\n");
      if (vote.t_subst < MAX_VOTE)
        vote.t_subst++;
      break;
    case A_INSERT: //fallthrough
    case C_INSERT: //fallthrough
    case G_INSERT: //fallthrough
    case T_INSERT: //fallthrough
      //fprintf(stderr, "Casting insertion of char %c\n", VoteChar(val));
      vote.insertions += VoteChar(val);
      break;
    default :
      fprintf(stderr, "ERROR:  Illegal vote type\n");
      assert(false);
  }
}



//  Return the substitution vote corresponding to  Ch .
static
Vote_Value_t
Matching_Vote(char ch) {

  switch  (ch) {
    case 'a':  return(A_SUBST);  break;
    case 'c':  return(C_SUBST);  break;
    case 'g':  return(G_SUBST);  break;
    case 't':  return(T_SUBST);  break;
  }

  fprintf(stderr, "Matching_Vote()-- invalid letter '%c'\n", ch);

  return(NO_VOTE);
}




//  Analyze the delta-encoded alignment in  delta[0 .. (deltaLen - 1)]
//  between  a_part  and  b_part  and store the resulting votes
//  about the a sequence in  G->reads[sub]. The alignment starts
//   a_offset  bytes in from the start of the a sequence in  G->reads[sub] .
//   a_len  and  b_len  are the lengths of the prefixes of  a_part  and
//   b_part , resp., that align.

//b read should be the "primary" overlaps for which are being analyzed
//BUT the votes are cast for the "a" read, with a "shifted" id == sub
//TODO make globalvote a local variable
void
Analyze_Alignment(Thread_Work_Area_t *wa,
                  char   *a_part, int32 a_len, int32 a_offset,
                  char   *b_part, //int32 b_len,
                  int32   sub) {

  assert(a_len >= 0);
  //assert(b_len >= 0);

  // ===== COLLECTING EVENTS =====
  //  Event counter. Each individual (1bp) mismatch/insertion/deletion is an event
  int32  ct = 0;

  wa->globalvote[ct].frag_sub  = -1;
  wa->globalvote[ct].align_sub = -1;
  wa->globalvote[ct].vote_val  = A_SUBST;   // Dummy value
  ct++;

  //position in a_part
  int32  i = 0;
  //position in b_part
  int32  j = 0;
  //position in "alignment" of a_part and b_part
  int32  p = 0;

  //fprintf(stderr, "Deltas:");
  //for (int32 k=0; k<wa->ped.deltaLen; k++) {
  //  fprintf(stderr, " %d ", wa->ped.delta[k]);
  //}
  //fprintf(stderr, "\n");
  for (int32 k = 0; k < wa->ped.deltaLen; k++) {
    //fprintf(stderr, "k=%d deltalen=%d  i=%d our of %d   j=%d out of %d\n", k, wa->ped.deltaLen, i, a_len, j, b_len);

    //  Add delta[k] - 1 matches or mismatches; +-1 encodes the 'continuation' of the insertion/deletion
    for (int32 m=1; m<abs(wa->ped.delta[k]); m++) {
      if (a_part[i] != b_part[j]) {
        wa->globalvote[ct].frag_sub  = i;
        wa->globalvote[ct].align_sub = p;
        wa->globalvote[ct].vote_val = SubstVote(b_part[j]);
        //fprintf(stderr, "Vote subst %c -> %c at %d\n", a_part[i], b_part[j], i);
        //fprintf(stderr, "ct %d\n", ct);
        ct++;
      }

      i++;  //assert(i <= a_len);
      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a negative delta, insert a base.

    if (wa->ped.delta[k] < 0) {
      //fprintf(stderr, "INSERT %c at %d #%d\n", b_part[j], i-1, p);
      wa->globalvote[ct].frag_sub  = i;
      wa->globalvote[ct].align_sub = p;
      wa->globalvote[ct].vote_val = InsVote(b_part[j]);
      //fprintf(stderr, "Vote ins %c at %d\n", b_part[j], i);
      //fprintf(stderr, "vote_val %d\n", InsVote(b_part[j]));
      //fprintf(stderr, "ct %d\n", ct);
      ct++;

      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a positive delta, delete the base.

    if (wa->ped.delta[k] > 0) {
      //fprintf(stderr, "DELETE %c at %d #%d\n", a_part[i], i, p);
      wa->globalvote[ct].frag_sub  = i;
      wa->globalvote[ct].align_sub = p;
      wa->globalvote[ct].vote_val  = DELETE;
      //fprintf(stderr, "Vote del at %d\n", i);
      //fprintf(stderr, "ct %d\n", ct);
      ct++;

      i++;  assert(i <= a_len);
      p++;
    }
  }

  // No more deltas.  While there is still sequence, add matches or mismatches.

  //fprintf(stderr, "k=done   i=%d our of %d   j=%d out of %d\n", i, a_len, j, b_len);

  while (i < a_len) {
    //fprintf(stderr, "k=done   i=%d our of %d   j=%d out of %d\n", i, a_len, j, b_len);

    if (a_part[i] != b_part[j]) {
      wa->globalvote[ct].frag_sub  = i;
      wa->globalvote[ct].align_sub = p;
      //fprintf(stderr, "Vote subst %c -> %c at %d\n", a_part[i], b_part[j], i);
      wa->globalvote[ct].vote_val = SubstVote(b_part[j]);
      //fprintf(stderr, "ct %d\n", ct);
      ct++;
    }

    i++;  //assert(i <= a_len);  //  Guaranteed, we're looping on this
    j++;  //assert(j <= b_len);
    p++;
  }

  assert(i <= a_len);
  wa->globalvote[ct].frag_sub  = i;
  wa->globalvote[ct].align_sub = p;

  //  For each identified change, add votes for some region around the change.
  //
  //  This is adding extra votes if the distance between two errors is larger than a kmer.
  //  Not sure why there are no 'matching bases' in this region.
  //
  //  X == changes, mismatch or indel
  //
  //                          ------- <- confirmed count added
  //                          -----   <- no_insert count added
  //  matching-bases} X 1 2 3 1 2 3 4 3 2 1 X {matching-bases
  //                    -----         -----
  //                    match         match
  //                    votes         votes
  //

  // ===== PROCESSING COLLECTED EVENTS =====
  assert(ct >= 1);
  //fprintf(stdout, "wa->G->Kmer_Len %d\n", wa->G->Kmer_Len);

  for (int32 event_idx = 1; event_idx <= ct; event_idx++) {
    // ===== CASTING MATCH/CONFIRMED/NO_INSERT VOTES BETWEEN EVENTS event_idx AND event_idx-1 =====
    const auto &prev_event = wa->globalvote[event_idx - 1];
    const int32 prev_event_end = prev_event.frag_sub + (prev_event.vote_val < A_INSERT ? 1 : 0);
    const int32 prev_event_dist = wa->globalvote[event_idx].frag_sub - prev_event_end;
    const int32 p_lo = (event_idx == 1) ? 0 : wa->G->End_Exclude_Len;
    const int32 p_hi = (event_idx == ct) ? prev_event_dist : prev_event_dist - wa->G->End_Exclude_Len;

    //  If distance to previous match is bigger than 'kmer' size cast flanking matching votes & mark confirmed positions
    if (prev_event_dist >= wa->G->Kmer_Len) {
      for (int32 p = 0; p < prev_event_dist; ++p) {
        const int32 part_pos = prev_event_end + p;
        const int32 a_pos = a_offset + part_pos;

        if (p < p_lo) {
          Cast_Vote(wa->G,
                    Matching_Vote(a_part[part_pos]),
                    a_pos,
                    sub);
        } else if (p < p_hi) {
          //p_lo <= p < p_hi
          if (wa->G->reads[sub].vote[a_pos].confirmed < MAX_VOTE)
            wa->G->reads[sub].vote[a_pos].confirmed++;

          if ((p < p_hi - 1) &&
              (wa->G->reads[sub].vote[a_pos].no_insert < MAX_VOTE))
            wa->G->reads[sub].vote[a_pos].no_insert++;
        } else {
          //p_hi <= p < prev_event_dist
          Cast_Vote(wa->G,
                    Matching_Vote(a_part[part_pos]),
                    a_pos,
                    sub);
        }
      }
    }

    // ===== ENDED CASTING MATCH/CONFIRMED/NO_INSERT VOTES =====

    //  Don't allow consecutive inserts.  If we aren't the last change, and there is non-adjacent
    //  previous (or this and the previous votes are not insertions), do another vote.

    // ===== CASTING EVENT event_idx =====
    if ((event_idx < ct)) { // && ((prev_match > 0) || (wa->globalvote[event_idx-1].vote_val <= T_SUBST) || (wa->globalvote[event_idx].vote_val <= T_SUBST)))
      //fprintf(stderr, "Try casting event #%d\n", event_idx);
      int32 next_match = wa->globalvote[event_idx + 1].align_sub - wa->globalvote[event_idx].align_sub - 1;

      // if our vote is outside of the bounds (meaning we have gaps at the start or end of the alignment), skip the vote
      assert(a_offset + wa->globalvote[event_idx].frag_sub >= 0);
      assert(wa->globalvote[event_idx].frag_sub <= a_len);
      //if (a_offset + wa->globalvote[event_idx].frag_sub < 0 || a_offset + wa->globalvote[event_idx].frag_sub >= a_len) {
      //  fprintf(stderr, "Corrected coord %d, a_len %d", a_offset + wa->globalvote[event_idx].frag_sub, a_len);
      //  fprintf(stderr, "Fail\n");
      //  continue;
      //}

      //TODO re-enable in some form?
      //Checking that sum of distances to the previous/next event is >= 9
      //if (prev_match + next_match >= wa->G->Vote_Qualify_Len)
      Cast_Vote(wa->G, wa->globalvote[event_idx].vote_val, a_offset + wa->globalvote[event_idx].frag_sub, sub);
    }
  }

  // ===== Finalizing cast insertions =====
  for (int32 a_pos = a_offset; a_pos < a_offset + a_len ; ++a_pos) {
    auto &insertions_str = wa->G->reads[sub].vote[a_pos].insertions;
    if (insertions_str.size() > 0 && insertions_str.back() != Vote_Tally_t::INSERTIONS_DELIM) {
      insertions_str += Vote_Tally_t::INSERTIONS_DELIM;
      wa->G->reads[sub].vote[a_pos].insertion_cnt++;
      //fprintf(stderr, "Increasing insertion count at position %d\n", a_pos);
    }
  }
}
