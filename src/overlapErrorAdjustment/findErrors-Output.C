
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
#include <map>

//void
//Output_Details(feParameters *G, uint32 i) {
//
//  fprintf(stderr, ">%d\n", G->bgnID + i);
//
//  for  (uint32 j=0;  G->reads[i].sequence[j] != '\0';  j++) {
//    const Vote_Tally_t &vote = G->reads[i].vote[j];
//    fprintf(stderr, "%3d: %c  conf %3d  deletes %3d | subst %3d %3d %3d %3d | no_insert %3d insert %3d sequences %s\n",
//            j,
//            j >= G->reads[i].clear_len ? toupper (G->reads[i].sequence[j]) : G->reads[i].sequence[j],
//            vote.confirmed,
//            vote.deletes,
//            vote.a_subst,
//            vote.c_subst,
//            vote.g_subst,
//            vote.t_subst,
//            vote.no_insert,
//            vote.insertion_cnt,
//            vote.insertions.empty() ? "" : vote.insertions.c_str());
//  }
//}

void
FPrint_Vote(FILE *fp, char base, const Vote_Tally_t &vote) {
  if (vote.all_but(base) == 0)
    fprintf(fp, "%c", base);
  else
    fprintf(fp, "[%c conf:no_ins %d:%d | del %d | subst %d:%d:%d:%d | ins %d sequences '%s']",
            base,
            vote.confirmed,
            vote.no_insert,
            vote.deletes,
            vote.a_subst, vote.c_subst, vote.g_subst, vote.t_subst,
            vote.insertion_cnt,
            vote.insertions.empty() ? "" : vote.insertions.c_str());
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
  }

  if  (vote.c_subst > max) {
    vote_t    = C_SUBST;
    max       = vote.c_subst;
  }

  if  (vote.g_subst > max) {
    vote_t    = G_SUBST;
    max       = vote.g_subst;
  }

  if  (vote.t_subst > max) {
    vote_t    = T_SUBST;
    max       = vote.t_subst;
  }

  int32 haplo_ct  =  ((vote.deletes >= MIN_HAPLO_OCCURS) +
      (vote.a_subst >= MIN_HAPLO_OCCURS) +
      (vote.c_subst >= MIN_HAPLO_OCCURS) +
      (vote.g_subst >= MIN_HAPLO_OCCURS) +
      (vote.t_subst >= MIN_HAPLO_OCCURS));

  if (vote_t != DELETE && base == VoteChar(vote_t)) {
    //fprintf(stderr, "SAME  base = %c, vote = %c\n", base, VoteChar(vote_t));
    return NO_VOTE;
  }

  if (vote.total() <= 1) {
    //fprintf(stderr, "FEW   total = %d <= 1\n", vote.total());
    return NO_VOTE;
  }

  if (2 * max <= vote.total()) {
    //fprintf(stderr, "WEAK  2*max = %d <= total = %d\n", 2*max, vote.total());
    return NO_VOTE;
  }

  if ((haplo_ct >= 2) && use_haplo_cnt) {
    //fprintf(stderr, "HAPLO haplo_ct=%d >= 2\n", haplo_ct);
    return NO_VOTE;
  }

  if (vote.confirmed > 2) {
    //fprintf(stderr, "Support for no correction: confirmed = %d\n", vote.confirmed);
    //Can not be triggered because confirmed < STRONG_CONFIRMATION_READ_CNT = 2
    assert(false);
    return NO_VOTE;
  }

  if (vote.confirmed == 1 && max <= 6) {
    //fprintf(stderr, "No correction was supported & small weight of vote: confirmed = %d ins_max = %d\n", vote.confirmed, max);
    return NO_VOTE;
  }

  return vote_t;
}

std::string
Check_Insert(const Vote_Tally_t &vote, char base, bool use_haplo_cnt) {

  std::map<std::string, uint32> insert_cnts;
  for (const auto &ins : vote.insertions_list()) {
    assert(!ins.empty());
    insert_cnts[ins] += 1;
  }

  int32 ins_haplo_ct = 0;

  int32 ins_max = 0;
  std::string ins_vote;
  for (const auto &ins_cnt : insert_cnts) {
    if (ins_cnt.second >= MIN_HAPLO_OCCURS) {
      ins_haplo_ct++;
    }
    if (ins_cnt.second > ins_max) {
      ins_max = ins_cnt.second;
      ins_vote = ins_cnt.first;
    }
  }

  //fprintf(stderr, "TEST   read %d position %d type %d (insert) -- ", i, j, ins_vote);

  if (vote.ins_total() <= 1) {
    //fprintf(stderr, "FEW   ins_total = %d <= 1\n", vote.ins_total());
    return "";
  }

  if (2 * ins_max <= vote.ins_total()) {
    //fprintf(stderr, "WEAK  2*ins_max = %d <= ins_total = %d\n", 2*ins_max, vote.ins_total());
    return "";
  }

  if ((ins_haplo_ct >= 2) && use_haplo_cnt) {
    //fprintf(stderr, "HAPLO ins_haplo_ct=%d >= 2\n", ins_haplo_ct);
    return "";
  }

  if (vote.no_insert >= 2) {
    //fprintf(stderr, "Support for no insert: no_insert = %d\n", vote.no_insert);
    //Can not be triggered because no_insert < STRONG_CONFIRMATION_READ_CNT = 2
    assert(false);
    return "";
  }

  if (vote.no_insert == 1 && ins_max <= 6) {
    //fprintf(stderr, "No insert was supported & small weight of vote: no_insert = %d ins_max = %d\n", vote.no_insert, ins_max);
    return "";
  }

  return ins_vote;
}

//TODO consider special case of two reads voting for different bases
// return false if nothing happened on the position and true otherwise
bool
Report_Position(const feParameters *G, const Frag_Info_t &read, uint32 pos,
    //Correction_Output_t out, std::ostream &os) {
    Correction_Output_t out, FILE *fp) {
  Vote_Tally_t vote = read.vote[pos];
  char base = read.sequence[pos];

  static const uint32 STRONG_CONFIRMATION_READ_CNT = 2;

  if (vote.all_but(base) == 0)
    return false;

  //Printing votes around position
  //FPrint_Votes(stderr, read, pos, /*locality radius*/5);

  bool corrected = false;

  if (vote.no_insert < STRONG_CONFIRMATION_READ_CNT) {
    //fprintf(stderr, "Checking read:pos %d:%d for insertion\n", out.readID, pos);
    std::string ins_str = Check_Insert(vote, base, G->Use_Haplo_Ct);
    if (ins_str.empty()) {
      //fprintf(stderr, "Read:pos %d:%d -- filtered out\n", out.readID, pos);
    } else {
      //fprintf(stderr, "Read:pos %d:%d -- corrected insertion. Insertion string ='%s'\n", out.readID, pos, ins_str.c_str());
      for (char c : ins_str) {
        out.type       = InsVote(c);
        out.pos        = pos;
        writeToFile(out, "correction2", fp);
      }
      corrected = true;
    }
  }

  if (vote.confirmed < STRONG_CONFIRMATION_READ_CNT) {
    //fprintf(stderr, "Checking read:pos %d:%d for del/subst\n", out.readID, pos);
    Vote_Value_t vote_t = Check_Del_Subst(vote, base, G->Use_Haplo_Ct);
    if (vote_t == NO_VOTE) {
      //fprintf(stderr, "Read:pos %d:%d -- filtered out\n", out.readID, pos);
    } else {
      out.type       = vote_t;
      out.pos        = pos;
      //fprintf(stderr, "Read:pos %d:%d -- corrected substitution/deletion\n", out.readID, pos);
      writeToFile(out, "correction3", fp);
      corrected = true;
    }
  }

  return corrected;
}

void
Output_Corrections(feParameters *G) {
  FILE *fp = AS_UTL_openOutputFile(G->outputFileName);
  //std::ofstream os(G->outputFileName);
  fprintf(stderr, "Output file: %s\n", G->outputFileName);

  for (uint32 read_idx = 0; read_idx < G->readsLen; ++read_idx) {
    //More debug ouptput
    //if (read_idx == 0)
    //Output_Details(G, read_idx);
    const Frag_Info_t &read = G->reads[read_idx];

    Correction_Output_t  out;
    out.keep_left   = (read.left_degree  < G->Degree_Threshold);
    out.keep_right  = (read.right_degree < G->Degree_Threshold);
    out.type        = IDENT;
    out.pos         = 0;
    out.readID      = G->bgnID + read_idx;

    //fprintf(stderr, "read %d clear_len %lu\n", out.readID, G->reads[read_idx].clear_len);
    //Writing stub event for the read
    writeToFile(out, "correction1", fp);

    if (read.sequence == NULL) {
      //fprintf(stderr, "Deleted fragment");
      // Deleted fragment
      continue;
    }

    //fprintf(stderr, "Checking positions\n");

    for (uint32 pos = 0; pos < read.clear_len; pos++) {
      Report_Position(G, read, pos, out, fp);
    }
  }

  AS_UTL_closeFile(fp, G->outputFileName);
}
