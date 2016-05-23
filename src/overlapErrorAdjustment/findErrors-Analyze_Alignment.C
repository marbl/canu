
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
 *    Brian P. Walenz on 2015-JUN-18
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Sergey Koren beginning on 2016-MAR-22
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Brian P. Walenz beginning on 2016-MAY-18
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"






//  Add vote val to G.reads[sub] at sequence position  p
static
void
Cast_Vote(feParameters *G,
          Vote_Value_t val,
          int32        pos,
          int32        sub) {
  int32  v=0;

  switch (val) {
    case DELETE:    if (G->reads[sub].vote[pos].deletes  < MAX_VOTE)  v = ++G->reads[sub].vote[pos].deletes;   break;
    case A_SUBST:   if (G->reads[sub].vote[pos].a_subst  < MAX_VOTE)  v = ++G->reads[sub].vote[pos].a_subst;   break;
    case C_SUBST:   if (G->reads[sub].vote[pos].c_subst  < MAX_VOTE)  v = ++G->reads[sub].vote[pos].c_subst;   break;
    case G_SUBST:   if (G->reads[sub].vote[pos].g_subst  < MAX_VOTE)  v = ++G->reads[sub].vote[pos].g_subst;   break;
    case T_SUBST:   if (G->reads[sub].vote[pos].t_subst  < MAX_VOTE)  v = ++G->reads[sub].vote[pos].t_subst;   break;
    case A_INSERT:  if (G->reads[sub].vote[pos].a_insert < MAX_VOTE)  v = ++G->reads[sub].vote[pos].a_insert;  break;
    case C_INSERT:  if (G->reads[sub].vote[pos].c_insert < MAX_VOTE)  v = ++G->reads[sub].vote[pos].c_insert;  break;
    case G_INSERT:  if (G->reads[sub].vote[pos].g_insert < MAX_VOTE)  v = ++G->reads[sub].vote[pos].g_insert;  break;
    case T_INSERT:  if (G->reads[sub].vote[pos].t_insert < MAX_VOTE)  v = ++G->reads[sub].vote[pos].t_insert;  break;
    case NO_VOTE:
      break;
    default :
      fprintf(stderr, "ERROR:  Illegal vote type\n");
      break;
  }

  //  Largely useless, just too much output.
  //fprintf(stderr, "Cast_Vote()-- sub %d at %d vote %d\n", sub, p, val);
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

void
Analyze_Alignment(Thread_Work_Area_t *wa,
                  char   *a_part, int32 a_len, int32 a_offset,
                  char   *b_part, int32 b_len,
                  int32   sub) {

  assert(a_len >= 0);
  assert(b_len >= 0);

  int32  ct = 0;

  //  Necessary??
  //memset(wa->globalvote, 0, sizeof(Vote_t) * AS_MAX_READLEN);

  wa->globalvote[ct].frag_sub  = -1;
  wa->globalvote[ct].align_sub = -1;
  wa->globalvote[ct].vote_val  = A_SUBST;   // Dummy value
  ct++;

  int32  i = 0;
  int32  j = 0;
  int32  p = 0;

  for (int32 k=0; k<wa->ped.deltaLen; k++) {
    //fprintf(stderr, "k=%d deltalen=%d  i=%d our of %d   j=%d out of %d\n", k, wa->ped.deltaLen, i, a_len, j, b_len);

    //  Add delta[k] matches or mismatches

    for (int32 m=1; m<abs(wa->ped.delta[k]); m++) {
      if (a_part[i] != b_part[j]) {
        wa->globalvote[ct].frag_sub  = i;
        wa->globalvote[ct].align_sub = p;

        switch (b_part[j]) {
          case 'a':  wa->globalvote[ct].vote_val = A_SUBST;  break;
          case 'c':  wa->globalvote[ct].vote_val = C_SUBST;  break;
          case 'g':  wa->globalvote[ct].vote_val = G_SUBST;  break;
          case 't':  wa->globalvote[ct].vote_val = T_SUBST;  break;
          default :
            fprintf(stderr, "ERROR:[1] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
            assert(0);
        }

        ct++;
      }

      i++;  //assert(i <= a_len);
      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a negative delta, insert a base.

    if (wa->ped.delta[k] < 0) {
      wa->globalvote[ct].frag_sub  = i - 1;
      wa->globalvote[ct].align_sub = p;

      //fprintf(stderr, "INSERT %c at %d #%d\n", b_part[j], i-1, p);

      switch (b_part[j]) {
        case 'a':  wa->globalvote[ct].vote_val = A_INSERT;  break;
        case 'c':  wa->globalvote[ct].vote_val = C_INSERT;  break;
        case 'g':  wa->globalvote[ct].vote_val = G_INSERT;  break;
        case 't':  wa->globalvote[ct].vote_val = T_INSERT;  break;
        default :
          fprintf(stderr, "ERROR:[2] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
          assert(0);
      }

      ct++;

      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a positive deta, delete the base.

    if (wa->ped.delta[k] > 0) {
      wa->globalvote[ct].frag_sub  = i;
      wa->globalvote[ct].align_sub = p;
      wa->globalvote[ct].vote_val  = DELETE;

      //fprintf(stderr, "DELETE %c at %d #%d\n", a_part[i], i, p);

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

      switch (b_part[j]) {
        case 'a':  wa->globalvote[ct].vote_val = A_SUBST;  break;
        case 'c':  wa->globalvote[ct].vote_val = C_SUBST;  break;
        case 'g':  wa->globalvote[ct].vote_val = G_SUBST;  break;
        case 't':  wa->globalvote[ct].vote_val = T_SUBST;  break;
        default :
          fprintf(stderr, "ERROR:[3] Bad sequence '%c' 0x%02x)\n", b_part[j], b_part[j]);
          assert(0);
      }

      ct++;
    }

    i++;  //assert(i <= a_len);  //  Guaranteed, we're looping on this
    j++;  //assert(j <= b_len);
    p++;
  }

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

  for (int32 i=1; i<=ct; i++) {
    int32  prev_match = wa->globalvote[i].align_sub - wa->globalvote[i - 1].align_sub - 1;
    int32  p_lo = (i == 1 ? 0 : wa->G->End_Exclude_Len);
    int32  p_hi = (i == ct ? prev_match : prev_match - wa->G->End_Exclude_Len);

    //  If distance to previous match is bigger than 'kmer' size, make a new vote.

    if (prev_match >= wa->G->Kmer_Len) {
      for (int32 p=0;  p<p_lo;  p++)
        Cast_Vote(wa->G,
                  Matching_Vote(a_part[wa->globalvote[i-1].frag_sub + p + 1]),
                            a_offset + wa->globalvote[i-1].frag_sub + p + 1,
                  sub);


      for (int32 p=p_lo;  p<p_hi;  p++) {
        int32 k = a_offset + wa->globalvote[i-1].frag_sub + p + 1;

        if (wa->G->reads[sub].vote[k].confirmed < MAX_VOTE)
          wa->G->reads[sub].vote[k].confirmed++;

        if ((p < p_hi - 1) &&
            (wa->G->reads[sub].vote[k].no_insert < MAX_VOTE))
          wa->G->reads[sub].vote[k].no_insert++;
      }

      for (int32 p=p_hi; p<prev_match; p++)
        Cast_Vote(wa->G,
                  Matching_Vote(a_part[wa->globalvote[i-1].frag_sub + p + 1]),
                            a_offset + wa->globalvote[i-1].frag_sub + p + 1,
                  sub);
    }

    //  Don't allow consecutive inserts.  If we aren't the last change, and there is non-adjacent
    //  previous (or this and the previous votes are not insertions), do another vote.

    if ((i < ct) &&
        ((prev_match > 0) ||
         (wa->globalvote[i-1].vote_val <= T_SUBST) ||
         (wa->globalvote[i  ].vote_val <= T_SUBST))) {
      int32 next_match = wa->globalvote[i + 1].align_sub - wa->globalvote[i].align_sub - 1;

      // if our vote is outside of the bounds (meaning we have gaps at the start or end of the alignment), skip the vote
      if (a_offset + wa->globalvote[i].frag_sub < 0 || a_offset + wa->globalvote[i].frag_sub >= a_len) {
         continue;
      }

      if (prev_match + next_match >= wa->G->Vote_Qualify_Len)
        Cast_Vote(wa->G,
                             wa->globalvote[i].vote_val,
                  a_offset + wa->globalvote[i].frag_sub,
                  sub);
    }
  }
}



