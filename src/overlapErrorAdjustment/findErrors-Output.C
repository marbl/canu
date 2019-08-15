
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
 *    Brian P. Walenz beginning on 2018-JAN-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "findErrors.H"

void
Output_Details(feParameters *G, uint32 i) {

  fprintf(stderr, ">%d\n", G->bgnID + i);

  for  (uint32 j=0;  G->reads[i].sequence[j] != '\0';  j++)
    fprintf(stderr, "%3d: %c  conf %3d  deletes %3d | subst %3d %3d %3d %3d | no_insert %3d insert %3d %3d %3d %3d\n",
            j,
            j >= G->reads[i].clear_len ? toupper (G->reads[i].sequence[j]) : G->reads[i].sequence[j],
            G->reads[i].vote[j].confirmed,
            G->reads[i].vote[j].deletes,
            G->reads[i].vote[j].a_subst,
            G->reads[i].vote[j].c_subst,
            G->reads[i].vote[j].g_subst,
            G->reads[i].vote[j].t_subst,
            G->reads[i].vote[j].no_insert,
            G->reads[i].vote[j].a_insert,
            G->reads[i].vote[j].c_insert,
            G->reads[i].vote[j].g_insert,
            G->reads[i].vote[j].t_insert);
}

enum class VoteCheckRes {
  NO_VARIANT, UNPROCESSED, FILTERED, CORRECTED_SD, CORRECTED_I
};

void
FPrint_Vote(FILE *fp, char base, const Vote_Tally_t &vote) {
  if (vote.all_but(base) == 0)
    fprintf(fp, "%c", base);
  else
    fprintf(fp, "[%c conf:no_ins %d:%d | del %d | subst %d:%d:%d:%d | ins %d %d %d %d]",
            base,
            vote.confirmed,
            vote.no_insert,
            vote.deletes,
            vote.a_subst, vote.c_subst, vote.g_subst, vote.t_subst,
            vote.a_insert, vote.c_insert, vote.g_insert, vote.t_insert);
}

void
FPrint_Votes(FILE *fp, const Frag_Info_t &read, uint32 j, uint32 loc_r) {
  assert(j < read.clear_len);
  uint32 s = j;
  uint32 e = j + 1;
  uint32 gathered_r = 0;
  while (gathered_r < loc_r) {
    if (s == 0)
      break;
    --s;
    if (read.vote[s].all_but(read.sequence[s]) == 0)
      ++gathered_r;
    else
      gathered_r = 0;
  }
  
  gathered_r = 0;
  while (gathered_r < loc_r) {
    if (e == read.clear_len)
      break;
    if (read.vote[e].all_but(read.sequence[e]) == 0)
      ++gathered_r;
    else
      gathered_r = 0;
    ++e;
  }

  if (e - s > 20) {
    //fprintf(fp, "Not printing locality for position %3d\n", j);
    fprintf(fp, "Warning: locality for position %3d\n was too big. Taking 10 on either side:", j);
    s = (j > 10) ? j - 10 : 0;
    e = (j + 12 < read.clear_len) ? j + 12 : read.clear_len;
  } 
 
  fprintf(fp, "Printing votes for positions [%3d, %3d) (around position %3d): \n", s, e, j);
  for (uint32 i = s; i < e; ++i) {
    if (i == j)
      fprintf(fp, "*");
    FPrint_Vote(fp, read.sequence[i], read.vote[i]);
    if (i == j)
      fprintf(fp, "*");
  }
  fprintf(fp, "\n");
}

Vote_Value_t 
Check_Del_Subst(const Vote_Tally_t &vote, char base, bool use_haplo_cnt) {
  Vote_Value_t  vote_t      = DELETE;
  int32         max       = vote.deletes;
  bool          is_change = true;

  if  (vote.a_subst > max) {
    vote_t    = A_SUBST;
    max       = vote.a_subst;
    is_change = (base != 'a');
  }

  if  (vote.c_subst > max) {
    vote_t    = C_SUBST;
    max       = vote.c_subst;
    is_change = (base != 'c');
  }

  if  (vote.g_subst > max) {
    vote_t    = G_SUBST;
    max       = vote.g_subst;
    is_change = (base != 'g');
  }

  if  (vote.t_subst > max) {
    vote_t    = T_SUBST;
    max       = vote.t_subst;
    is_change = (base != 't');
  }

  int32 haplo_ct  =  ((vote.deletes >= MIN_HAPLO_OCCURS) +
      (vote.a_subst >= MIN_HAPLO_OCCURS) +
      (vote.c_subst >= MIN_HAPLO_OCCURS) +
      (vote.g_subst >= MIN_HAPLO_OCCURS) +
      (vote.t_subst >= MIN_HAPLO_OCCURS));

  //  The original had a gargantuajn if test (five clauses, all had to be true) to decide if a record should be output.
  //  It was negated into many small tests if we should skip the output.
  //  A side effect is that we can abort a little earlier in two cases (and we don't even bother).


  //fprintf(stderr, "TEST   read %d position %d type %d -- ", i, j, vote);

  //  (total > 1)
  if (vote.total() <= 1) {
    //fprintf(stderr, "FEW   total = %d <= 1\n", total);
    fprintf(stderr, "HERE1\n");
    return NO_VOTE;
  }

  //  (2 * max > total)
  if (2 * max <= vote.total()) {
    //fprintf(stderr, "WEAK  2*max = %d <= total = %d\n", 2*max, total);
    fprintf(stderr, "HERE2\n");
    return NO_VOTE;
  }

  //  (is_change == true)
  if (is_change == false) {
    //fprintf(stderr, "SAME  is_change = %d\n", is_change);
    fprintf(stderr, "HERE3\n");
    return NO_VOTE;
  }

  //  ((haplo_ct < 2) || (G->Use_Haplo_Ct == false))
  if ((haplo_ct >= 2) && use_haplo_cnt) {
    //fprintf(stderr, "HAPLO haplo_ct=%d >= 2 AND Use_Haplo_Ct = %d\n", haplo_ct, G->Use_Haplo_Ct);
    fprintf(stderr, "HERE4\n");
    return NO_VOTE;
  }

  //  ((vote.confirmed == 0) ||
  //   ((vote.confirmed == 1) && (max > 6)))
  if (vote.confirmed > 2) {
    //fprintf(stderr, "INDEL confirmed = %d", vote.confirmed);
    fprintf(stderr, "HERE5\n");
    return NO_VOTE;
  }

  if (vote.confirmed == 1 && max <= 6) {
    //fprintf(stderr, "INDEL confirmed = %d max = %d\n", vote.confirmed, max);
    fprintf(stderr, "HERE6\n");
    return NO_VOTE;
  }

  return vote_t;
}

Vote_Value_t 
Check_Insert(const Vote_Tally_t &vote, char base, bool use_haplo_cnt) {
    Vote_Value_t  ins_vote_t = A_INSERT;
    int32         ins_max  = vote.a_insert;

    if  (ins_max < vote.c_insert) {
      ins_vote_t = C_INSERT;
      ins_max  = vote.c_insert;
    }

    if  (ins_max < vote.g_insert) {
      ins_vote_t = G_INSERT;
      ins_max  = vote.g_insert;
    }

    if  (ins_max < vote.t_insert) {
      ins_vote_t = T_INSERT;
      ins_max  = vote.t_insert;
    }

    int32 ins_haplo_ct = ((vote.a_insert >= MIN_HAPLO_OCCURS) +
                          (vote.c_insert >= MIN_HAPLO_OCCURS) +
                          (vote.g_insert >= MIN_HAPLO_OCCURS) +
                          (vote.t_insert >= MIN_HAPLO_OCCURS));

    //fprintf(stderr, "TEST   read %d position %d type %d (insert) -- ", i, j, ins_vote);

    if (vote.ins_total() <= 1) {
      //fprintf(stderr, "FEW   ins_total = %d <= 1\n", ins_total);
      fprintf(stderr, "OPPA1\n");
      return NO_VOTE;
    }

    if (2 * ins_max <= vote.ins_total()) {
      //fprintf(stderr, "WEAK  2*ins_max = %d <= ins_total = %d\n", 2*ins_max, ins_total);
      fprintf(stderr, "OPPA2\n");
      return NO_VOTE;
    }

    if ((ins_haplo_ct >= 2) && use_haplo_cnt) {
      //fprintf(stderr, "HAPLO ins_haplo_ct=%d >= 2 AND Use_Haplo_Ct = %d\n", ins_haplo_ct, G->Use_Haplo_Ct);
      fprintf(stderr, "OPPA3\n");
      return NO_VOTE;
    }

    if (vote.no_insert >= 2) {
      //fprintf(stderr, "INDEL no_insert = %d\n", vote.no_insert);
      fprintf(stderr, "OPPA4\n");
      return NO_VOTE;
    }

    if (vote.no_insert == 1 && ins_max <= 6) {
      //fprintf(stderr, "INDEL no_insert = %d ins_max = %d\n", vote.no_insert, ins_max);
      fprintf(stderr, "OPPA5\n");
      return NO_VOTE;
    }

    return ins_vote_t;
}

//TODO special case of two reads voting for different bases
VoteCheckRes 
Check_Position(const feParameters *G, const Frag_Info_t &read, uint32 j, Correction_Output_t *out) {
  Vote_Tally_t vote = read.vote[j];
  char base = read.sequence[j];

  static const uint32 STRONG_CONFIRMATION_READ_CNT = 2;

  //FIXME understand why it can fail on a fully confirmed base
  //assert(vote.no_insert >= vote.confirmed);

  if (vote.all() == 0)
    return VoteCheckRes::NO_VARIANT;

  if (vote.no_insert < vote.confirmed) {
    fprintf(stderr, "WARN: no_insert %d ; confirmed %d \n", vote.no_insert, vote.confirmed);
    FPrint_Vote(stderr, read.sequence[j], vote);
  }

  FPrint_Votes(stderr, read, j, /*locality radius*/5);

  VoteCheckRes code = VoteCheckRes::UNPROCESSED;

  if (vote.confirmed < STRONG_CONFIRMATION_READ_CNT) {
    Vote_Value_t vote_t = Check_Del_Subst(vote, base, G->Use_Haplo_Ct);
    if (vote_t == NO_VOTE) {
      code = VoteCheckRes::FILTERED;
    } else {
      //fprintf(stderr, "CORRECT!\n");
      out->type       = vote_t;
      out->pos        = j;
      return VoteCheckRes::CORRECTED_SD;
    }
  }  

  if  (vote.no_insert < STRONG_CONFIRMATION_READ_CNT) {
    Vote_Value_t vote_t = Check_Insert(vote, base, G->Use_Haplo_Ct);
    if (vote_t == NO_VOTE) {
      code = VoteCheckRes::FILTERED;
    } else {
    //fprintf(stderr, "INSERT!\n");
      //fprintf(stderr, "CORRECT!\n");
      out->type       = vote_t;
      out->pos        = j;
      return VoteCheckRes::CORRECTED_I;
    }
  }

  return code;
}


void
Output_Corrections(feParameters *G) {
  Correction_Output_t  out;

  FILE *fp = AS_UTL_openOutputFile(G->outputFileName);
  fprintf(stderr, "Output file: %s\n", G->outputFileName);

  for (uint32 i=0; i<G->readsLen; i++) {
    //More debug ouptput
    //if (i == 0)
    //Output_Details(G, i);
    const Frag_Info_t &read = G->reads[i];

    out.keep_left   = (read.left_degree  < G->Degree_Threshold);
    out.keep_right  = (read.right_degree < G->Degree_Threshold);
    out.type        = IDENT;
    out.pos         = 0;
    out.readID      = G->bgnID + i;

    fprintf(stderr, "read %d clear_len %lu\n", out.readID, G->reads[i].clear_len);
    writeToFile(out, "correction1", fp);
    fprintf(stderr, "written");

    if (G->reads[i].sequence == NULL) {
      fprintf(stderr, "Deleted fragment");
      // Deleted fragment
      continue;
    } else {
      fprintf(stderr, "Checking positions\n");
    }

    for (uint32 j=0; j<read.clear_len; j++) {
      VoteCheckRes res = Check_Position(G, read, j, &out);
      switch (res) { 
        case VoteCheckRes::NO_VARIANT:
          break;
        case VoteCheckRes::CORRECTED_SD:
          writeToFile(out, "correction2", fp);
          fprintf(stderr, "Read:pos %d:%d -- corrected substitution/deletion\n", out.readID, j);
          break;
        case VoteCheckRes::CORRECTED_I:
          fprintf(stderr, "Read:pos %d:%d -- corrected insertion\n", out.readID, j);
          writeToFile(out, "correction3", fp);
          break;
        case VoteCheckRes::FILTERED:
          fprintf(stderr, "Read:pos %d:%d -- filtered out\n", out.readID, j);
          //writeToFile(out, "correction4", fp);
          break;
        case VoteCheckRes::UNPROCESSED:
          fprintf(stderr, "Read:pos %d:%d -- unprocessed\n", out.readID, j);
          //writeToFile(out, "correction5", fp);
          break;
        default:
          assert(false);
      }
    }
  }

  AS_UTL_closeFile(fp, G->outputFileName);
}
