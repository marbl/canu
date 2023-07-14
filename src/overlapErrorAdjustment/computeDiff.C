
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

#include "computeDiff.H"
#include <algorithm>

static
int32
Nucl2Int(char nucl) {
  switch (nucl) {
    case 'a':
      return 0;
    case 'c':
      return 1;
    case 'g':
      return 2;
    case 't':
      return 3;
    default:
      assert(false);
  }
}

//If prefix_len is negative -- go back in the string
static
int32
Convert2Int(const char* seq, int32 prefix_len) {
  assert(prefix_len != 0 && std::abs(prefix_len) < 16);
  int32 ans = 0;
  if (prefix_len > 0) {
    for (int32 i = 0; i < prefix_len; ++i) {
      assert(seq[i] != '\0');
      ans = (ans << 2) | Nucl2Int(seq[i]);
    }
  } else {
    for (int32 i = 0; i < -prefix_len; ++i) {
      ans |= Nucl2Int(*(seq - i)) << (2 * i);
    }
  }
  assert(ans < (1 << (2 * std::abs(prefix_len))));
  return ans;
}

//Collects kmer stats for the region of length |reg_len|
//If reg_len is negative -- go back in the string
static
void
CollectKmerStat(const char* seq, int32 reg_len, int32 kmer_len, int32 *stats) {
  assert(kmer_len > 0);
  assert(reg_len != 0);
  memset(stats, 0, sizeof(int32) * (1 << (2 * kmer_len)));
  if (reg_len > 0) {
    for (int32 i = 0; (i + kmer_len) <= reg_len; i = i + kmer_len) {
      stats[Convert2Int(seq + i, kmer_len)]++;
    }
  } else {
    for (int32 i = 0; (i + kmer_len) <= -reg_len; i = i + kmer_len) {
      stats[Convert2Int(seq - i, -kmer_len)]++;
    }
  }
}

static
bool
CheckTrivialDNA(const char* seq, int32 remaining, int32 offset, int32 size_factor, int32 repeat_cnt) {
  //TODO configure trivial DNA analysis
  //static const int32 SIZE_FACTOR = 6;
  //static const int32 REPEAT_NUM = 5;
  //static const int32 SIZE_FACTOR = 4;
  //static const int32 REPEAT_NUM = 3;
  static const int32 MIN_K = 2;
  static const int32 MAX_K = 5;

  int32 stats_buff[1 << (2 * MAX_K)];
  for (int32 k = MIN_K; k <= MAX_K; ++k) {
    const int32 possible_kmer_cnt = 1 << (2 * k);
    int32 reg_len = k * size_factor;

    //exploring sequence to the right
    //fprintf(stderr, "checking upstream k=%d, init_shift=%d\n", k, std::max(-k, -offset));
    for (int32 shift = std::max(-k, -offset); shift < k; ++shift) {

      if (reg_len + shift > remaining)
        break;
      CollectKmerStat(seq + shift, reg_len, k, stats_buff);
      if (*std::max_element(stats_buff, stats_buff + possible_kmer_cnt) >= repeat_cnt) {
        //comment out!
        //char subbuff[reg_len + 1];
        //memcpy(subbuff, seq + shift, reg_len);
        //subbuff[reg_len] = '\0';
        //fprintf(stderr, "Trivial DNA (k=%d) upstream\n", k);
        //fprintf(stderr, "%s\n", subbuff);
        return true;
      }
    }

    //fprintf(stderr, "checking downstream k=%d, init_shift=%d\n", k, std::max(-k, -remaining));
    //exploring sequence to the left
    for (int32 shift = std::max(-k, -remaining); shift < k; ++shift) {
      if (reg_len + shift > offset)
        break;
      CollectKmerStat(seq - shift - 1, -reg_len, k, stats_buff);
      if (*std::max_element(stats_buff, stats_buff + possible_kmer_cnt) >= repeat_cnt) {
        //comment out!
        //char subbuff[reg_len + 1];
        //memcpy(subbuff, seq - shift - reg_len, reg_len);
        //subbuff[reg_len] = '\0';
        //fprintf(stderr, "Trivial DNA (k=%d) downstream\n", k);
        //fprintf(stderr, "%s\n", subbuff);
        return true;
      }
    }
  }
  return false;
}

static
bool
CheckTrivialDNA(const char* a_part, const char* b_part,
                int32 a_len, int32 b_len,
                int32 a_pos, int32 b_pos,
                int32 size_factor, int32 repeat_cnt) {
  //fprintf(stderr, "Checking for trivial DNA in around positions %d, %d\n", a_pos, b_pos);

  if (CheckTrivialDNA(a_part + a_pos, /*remaining*/a_len - a_pos, /*offset*/a_pos, size_factor, repeat_cnt)) {
    //fprintf(stderr, "Trivial DNA detected in A around position %d\n", a_pos);
    return true;
  }
  if (CheckTrivialDNA(b_part + b_pos, /*remaining*/b_len - b_pos, /*offset*/b_pos, size_factor, repeat_cnt)) {
    //fprintf(stderr, "Trivial DNA detected in B around position %d\n", b_pos);
    return true;
  }
  //fprintf(stderr, "NON-Trivial DNA!!!\n", a_pos);
  return false;
}

std::pair<size_t, size_t>
ComputeErrors(const char* const a_part, const char* const b_part,
    int32 delta_len, int32 *deltas,
    int32 a_len, int32 b_len,
    bool check_trivial_dna,
    uint32 ignore_flank) {

  static const int32 MM_SIZE_FACTOR = 6;
  static const int32 MM_REPEAT_NUM = 5;

  //  The following two IND constants were tested against v2.1.1 and v2.2.
  //    In v2.1.1, IND_SIZE_FACTOR=6 and IND_REPEAT_NUM=5.
  //    They were changed to 4 and 3, respectively in commit
  //      cb94432adf36cbb77c7ec22517a50e3051abd066
  //
  //  The changed values resulted in many more misassemblies on
  //  human-chm13-hifi-20k, from 9 in v2.1.1 to 14 with the new values.  (The
  //  effect is also observable in dmel f1, but I don't have exact numbers.)
  //
  //  Switching back to the v2.1.1 values (6 and 5) actually drops
  //  misassembles down to 6, due to other changes made around commit
  //  cb94432.
  //
  static const int32 IND_SIZE_FACTOR = 6;
  static const int32 IND_REPEAT_NUM = 5;

  //  Event counter. Each individual (1bp) mismatch/insertion/deletion is an event
  int32  all_ct = 0;
  //  Processed event counter
  int32  ct = 0;
  //position in a_part
  int32  i = 0;
  //position in b_part
  int32  j = 0;
  //position in "alignment" of a_part and b_part
  int32  p = 0;

  auto cnt_event_f = [&](int32 size_factor, int32 repeat_cnt) {
    if (i < ignore_flank ||
        j < ignore_flank ||
        i + ignore_flank >= a_len ||
        j + ignore_flank >= b_len) {
      return false;
    }
    if (check_trivial_dna &&
        CheckTrivialDNA(a_part, b_part, a_len, b_len, i, j,
                         size_factor, repeat_cnt)) {
      return false;
    }
    return true;
  };

  for (int32 k=0; k < delta_len; k++) {
    //fprintf(stderr, "k=%d (deltalen=%d) deltas[k]=%d i=%d out of %d  j=%d out of %d\n", k, delta_len, deltas[k], i, a_len, j, b_len);

    //  Add delta[k] - 1 matches or mismatches; +-1 encodes the 'continuation' of the insertion/deletion
    for (int32 m=1; m<abs(deltas[k]); m++) {
      if (a_part[i] != b_part[j]) {
        //Substitution at i in a_part (p in "alignment")
        //fprintf(stderr, "SUBST %c -> %c at %d #%d\n", a_part[i], b_part[j], i, p);

        all_ct++;
        ct += (int32) cnt_event_f(MM_SIZE_FACTOR, MM_REPEAT_NUM);
      }

      i++;  //assert(i <= a_len);
      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a negative delta, insert a base.

    if (deltas[k] < 0) {
      //Insertion at i - 1 in a_part (p in "alignment")
      //fprintf(stderr, "INSERT %c at %d #%d\n", b_part[j], i-1, p);
      all_ct++;
      ct += (int32) cnt_event_f(IND_SIZE_FACTOR, IND_REPEAT_NUM);

      j++;  //assert(j <= b_len);
      p++;
    }

    //  If a positive delta, delete the base.

    if (deltas[k] > 0) {
      //Deletion at i in a_part (p in "alignment")
      //fprintf(stderr, "DELETE %c at %d #%d\n", a_part[i], i, p);
      all_ct++;
      ct += (int32) cnt_event_f(IND_SIZE_FACTOR, IND_REPEAT_NUM);

      i++;  //assert(i <= a_len);
      p++;
    }
  }

  // No more deltas.  While there is still sequence, add matches or mismatches.
  while (i < a_len) {
    if (a_part[i] != b_part[j]) {
      //fprintf(stderr, "SUBST %c -> %c at %d #%d\n", a_part[i], b_part[j], i, p);
      //Substitution at i in a_part (p in "alignment")
      all_ct++;
      ct += (int32) cnt_event_f(MM_SIZE_FACTOR, MM_REPEAT_NUM);
    }

    i++;  //assert(i <= a_len);  //  Guaranteed, we're looping on this
    j++;  //assert(j <= b_len);
    p++;
  }

  //if (all_ct > 0) {
  //  fprintf(stderr, "Reported %d out of %d\n", ct, all_ct);
  //}

  assert(i == a_len);
  assert(j == b_len);

  return std::make_pair(ct, p);
}
