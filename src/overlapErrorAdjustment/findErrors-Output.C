
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

//FIXME update?
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
    fprintf(fp, "[%c conf:conf_no_ins %d:%d | del %d | subst %d:%d:%d:%d | no_ins:ins %d:%d sequences '%s']",
            base,
            vote.confirmed,
            vote.conf_no_insert,
            vote.deletes,
            vote.a_subst, vote.c_subst, vote.g_subst, vote.t_subst,
            vote.no_insert,
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


//  Return the highest count DELETE or SUBST vote.
//  Returns NO_VOTE if
//    The SUBST is what is there already
//    There is only one vote
//    The highest vote is less than 50% of votes
//    There are >1 votes (DEL/SUBST) with weight >= Haplo_Confirm
//    The vote is 'confirmed == 1' and weight <= 6
//
//  'confirmed' means that a read has this base and it is at least 3 bases
//  away from an alignment difference.
//
Vote_Value_t
Check_Del_Subst(const Vote_Tally_t &vote, char base, int32 Haplo_Confirm) {
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

  int32 haplo_ct  =  ((vote.deletes >= Haplo_Confirm) +
                      (vote.a_subst >= Haplo_Confirm) +
                      (vote.c_subst >= Haplo_Confirm) +
                      (vote.g_subst >= Haplo_Confirm) +
                      (vote.t_subst >= Haplo_Confirm));

  if (haplo_ct > 1) {
    //fprintf(stderr, "HAPLO haplo_ct=%d >= 2\n", haplo_ct);
    return NO_VOTE;
  }

  if (max <= 6 * vote.confirmed) {
    //fprintf(stderr, "No correction was supported & small weight of vote: confirmed = %d max = %d\n", vote.confirmed, max);
    return NO_VOTE;
  }

  return vote_t;
}

//  Returns a string to insert.
//    Empty string if only one read supports the insertion
//    Empty string if the insertion is not dominant (including reads that vote for no insertion)
//    Empty string if more than 2 votes for an insertion (or one insert and one no insertion)
//    Empty string if EXACTLY one read confirms no insertion and 6 or fewer vote for an insertion.
//
std::string
Check_Insert(const Vote_Tally_t &vote, char base, int32 Haplo_Confirm) {

  std::map<std::string, uint32> insert_cnts;
  for (const auto &ins : vote.insertions_list()) {
    assert(!ins.empty());
    insert_cnts[ins]++;
  }

  int32 ins_haplo_ct = 0;

  int32 ins_max = 0;
  std::string ins_vote;
  for (const auto &ins_cnt : insert_cnts) {
    if (ins_cnt.second >= Haplo_Confirm) {
      ins_haplo_ct++;
    }
    if (ins_cnt.second > ins_max) {
      ins_max = ins_cnt.second;
      ins_vote = ins_cnt.first;
    }
  }

  //considering empty insertion as a valid vote
  if (vote.no_insert >= Haplo_Confirm)
    ins_haplo_ct++;

  //fprintf(stderr, "TEST   read %d position %d type %d (insert) -- ", i, j, ins_vote);

  if (vote.insertion_cnt <= 1) {
    //fprintf(stderr, "FEW   ins_total = %d <= 1\n", vote.ins_total());
    return "";
  }

  //considering empty insertion as a valid vote
  if (2 * ins_max <= vote.insertion_cnt + vote.no_insert) {
    //fprintf(stderr, "WEAK  2*ins_max = %d <= total = %d\n", 2*ins_max, vote.insertion_cnt + vote.no_insert);
    return "";
  }

  //no need to check non-insertions here, since we are in no_insert < STRONG_CONFIRMATION_READ_CNT = 2 case
  if (ins_haplo_ct > 1) {
    //fprintf(stderr, "HAPLO ins_haplo_ct=%d >= 2\n", ins_haplo_ct);
    return "";
  }

  if (ins_max <= 6 * vote.conf_no_insert) {
    //fprintf(stderr, "No insert was supported & small weight of vote: no_insert = %d ins_max = %d\n", vote.no_insert, ins_max);
    return "";
  }

  return ins_vote;
}


//  Returns true if a Haplo_Confirm delete is present AND 
//

bool Is_Het_Del(const Vote_Tally_t &vote, int32 Haplo_Confirm) {
  //TODO consider using STRONG_CONFIRMATION_READ_CNT
  //if (vote.no_insert >= Haplo_Confirm && vote.insertion_cnt >= Haplo_Confirm)
  //  return true;

  //int32 haplo_ct  = ((vote.deletes >= Haplo_Confirm) +
  //    (vote.a_subst >= Haplo_Confirm) +
  //    (vote.c_subst >= Haplo_Confirm) +
  //    (vote.g_subst >= Haplo_Confirm) +
  //    (vote.t_subst >= Haplo_Confirm));

  //return haplo_ct >= 2;
  return((vote.deletes                >= Haplo_Confirm) &&
         (vote.total() - vote.deletes >= Haplo_Confirm));
}

std::vector<uint32>
Find_Het_Del_Positions(const feParameters *G, const Frag_Info_t &read, int32 Haplo_Confirm) {
  std::vector<uint32> answer;
  for (uint32 pos = 0; pos < read.clear_len; pos++) {
    if (Is_Het_Del(read.vote[pos], Haplo_Confirm)) {
      answer.push_back(pos);
    }
  }
  return answer;
}

//TODO consider special case of two reads voting for different bases
// return false if nothing happened on the position and true otherwise
bool
Report_Position(const feParameters *G, const Frag_Info_t &read, uint32 pos, bool block_deletion,
    //Correction_Output_t out, std::ostream &os) {
    Correction_Output_t out, FILE *fp) {
  Vote_Tally_t vote = read.vote[pos];
  char base = read.sequence[pos];

  //static const uint32 STRONG_CONFIRMATION_READ_CNT = 2;

  if (vote.all_but(base) == 0)
    return false;

  //Printing votes around position
  //FPrint_Votes(stderr, read, pos, /*locality radius*/5);

  bool corrected = false;

  if (vote.conf_no_insert < G->Base_Confirm || vote.conf_no_insert_rc < G->Opposite_Confirm) {
    //fprintf(stderr, "Checking read:pos %d:%d for insertion\n", out.readID, pos);
    std::string ins_str = Check_Insert(vote, base, G->Haplo_Confirm);
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

  if (vote.confirmed < G->Base_Confirm || vote.confirmed_rc < G->Opposite_Confirm) {
    //fprintf(stderr, "Checking read:pos %d:%d for del/subst\n", out.readID, pos);
    Vote_Value_t vote_t = Check_Del_Subst(vote, base, G->Haplo_Confirm);
    if (vote_t == NO_VOTE) {
      //fprintf(stderr, "Read:pos %d:%d -- filtered out\n", out.readID, pos);
    } else if (vote_t == DELETE && block_deletion) {
      //fprintf(stderr, "Read:pos %d:%d -- deletion blocked\n", out.readID, pos);
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
  FILE *fp = merylutil::openOutputFile(G->outputFileName);
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

    std::vector<uint32> het_del_pos;

    if (G->Het_Del_Freeze > 0) {
      //no reason to use different value;
      //freezing one position seems enough to not falsely correct ambiguous deletion
      assert(G->Het_Del_Freeze == 1);
      het_del_pos = Find_Het_Del_Positions(G, read, G->Haplo_Confirm);
    }

    auto het_it = het_del_pos.begin();
    for (uint32 pos = 0; pos < read.clear_len; pos++) {
      bool block_deletion = false;
      while (het_it != het_del_pos.end() && pos > *het_it + G->Het_Del_Freeze)
        ++het_it;
      //here het_it points to a value >= pos - freeze
      if (het_it != het_del_pos.end() && *het_it <= pos + G->Het_Del_Freeze) {
        //fprintf(stderr, "Ignoring position %d too close to a het at position %d\n");
        //blocking deletion near position with heterozygous deletion -- too ambiguous
        block_deletion = true;
      }
      Report_Position(G, read, pos, block_deletion, out, fp);
    }
  }

  merylutil::closeFile(fp, G->outputFileName);
}
