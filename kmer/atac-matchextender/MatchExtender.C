// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Dan Fasulo
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <cassert>
#include <csignal>
#include <iostream>
#include <vector>
#include <algorithm>

#include "MEMatch.H"
#include "atacreader.H"
#include "bio++.H"

//
// CONSTANTS
//

int   DEF_MIN_END_RUN_LEN   = 10;
int   DEF_MAX_MM_BLOCK      = 3;
int   DEF_MIN_BLOCK_SEP     = 20;
float DEF_MIN_IDENTITY      = 0.95;
int   DEF_MAX_NBR_SEP       = 100;
int   DEF_MAX_NBR_PATH_MM   = 5;
int   DEF_VERBOSE           = 0;

//
// GLOBAL PARAMETERS
//

MEMatch *curmatch_G         = NULL;   // The match currently being operated on


//
// SIMPLE COMPARISON UTILITY TYPE TO SORT MATCHES
//

class MEMatchPtrCmp {
public:
  int operator()(const MEMatch *m1, const MEMatch *m2) {
    return (int) (*m1 < *m2);
  }
};



static void
print_usage(void) {
  cerr << "Usage: MatchExtender [options] <infile> <outfile>\n"
       << "   <infile>    - Input ATAC file with parentless exact matches\n"
       << "   <outfile>   - Output ATAC file\n"
       << "Options:\n"
       << "   -x <string> - Prefix of axis1 (read from infile otherwise)\n"
       << "   -y <string> - Prefix of axis2 (read from infile otherwise)\n"
       << "   -v          - Verbose (debugging) mode.  Lots of output!\n"
       << "   -h          - Print this help message\n"
       << "Parameters:\n"
       << "Param    Meaning                                       Def. value\n"
       << "-----------------------------------------------------------------\n"
       << " -E n    Minimum exact match run at ends               " << DEF_MIN_END_RUN_LEN << '\n'
       << " -B n    Maximum of n mm's without separating block    " << DEF_MAX_MM_BLOCK << '\n'
       << " -S n    Minimum block size to separate mm groups      " << DEF_MIN_BLOCK_SEP << '\n'
       << " -I n    Minimum % identity of match, 0 <= n <= 1.0    " << DEF_MIN_IDENTITY << '\n'
       << " -P n    Maximum separation of 'nearby' matches        " << DEF_MAX_NBR_SEP << '\n'
       << " -D n    Max mm's allowed to connect nearby matches    " << DEF_MAX_NBR_PATH_MM << '\n';
    }


static void
log_parameters(void) {
  cerr << " * Parameter values:\n"
       << " *     Minimum exact match run at ends of matches            = " << DEF_MIN_END_RUN_LEN << '\n'
       << " *     Maximum number of mismatches without separating block = " << DEF_MAX_MM_BLOCK << '\n'
       << " *     Minimum exact run separating regions of mismatches    = " << DEF_MIN_BLOCK_SEP << '\n'
       << " *     Minimum % identity of match                           = " << DEF_MIN_IDENTITY << '\n'
       << " *     Maximum separation of 'nearby' matches                = " << DEF_MAX_NBR_SEP << '\n'
       << " *     Maximum mismatches when connecting nearby matches     = " << DEF_MAX_NBR_PATH_MM << '\n';
}


static bool
trim_to_pct(vector<MEMatch *>& matches, u32bit midx, float pct) {

  if (DEF_VERBOSE >= 3)
    fprintf(stderr, "trim_to_pct()\n");

  bool     result       = false;        // Return true iff trimming was performed
  SeqOff_t start        = 0;
  SeqOff_t len          = 0;
  SeqOff_t best_start   = 0;
  SeqOff_t best_len     = 0;

  MEMatch *m = matches[midx];

  FastAAccessor   A(m->s1(), false);
  FastAAccessor   B(m->s2(), m->isReversed());

  B.setReverseComplementRange(m->pos2(), m->len());

  for (start = 0; start < m->len(); ++start) {
    SeqOff_t best_run_len = 0;
    SeqOff_t sum          = 0;

    A.setPosition(m->pos1() + start);
    B.setPosition(m->pos2() + start);

    for (len = 1; len <= (m->len() - start); ++len) {
      char c1 = *A;
      char c2 = *B;

      if (validSymbol[c1] &&
          validSymbol[c2] &&
          IUPACidentity[c1][c2])
	sum++;

      if (((float) sum / (float) len) >= pct)
	best_run_len = len;

      ++A;
      ++B;
    }

    // Special case : if the whole string is okay, don't check any
    // subranges
    //
    if ((start == 0) && (best_run_len == m->len()))
      return false;

    if (best_run_len > best_len) {
      best_len = best_run_len;
      best_start = start;   
    }
  }

  if (best_len < m->len()) {
    if (DEF_VERBOSE >= 2)
      cout << "Trimming to substring with start = " << best_start 
	   << " and len = " << best_len << " for % identity\n";
    result = true;
    if (start > 0)
      m->trimFront(best_start);

    if (best_len < m->len())
      m->trimEnd(m->len() - best_len);
  }

  return result;
}


static void
extend_match_backward(vector<MEMatch *>& matches,
                      u32bit midx,
		      SeqOff_t min_start_pos) {

  if (DEF_VERBOSE >= 3)
    fprintf(stderr, "extend_match_backward()\n");

  // Assumes when traveling backwards that we will never run into
  // another match (otherwise, that match would have been forward
  // extended previously).

  int num_recent_mismatches = 0;
  MEMatch *m = matches[midx];
  int good_run_len = (int) m->len();
  SeqOff_t num_pending = 0;

  FastAAccessor   A(m->s1(), false);
  FastAAccessor   B(m->s2(), m->isReversed());

  B.setReverseComplementRange(m->pos2(), m->len());

  A.setPosition(m->pos1());
  B.setPosition(m->pos2());

  while ((A.getPosition() > min_start_pos) && (B.isValid())) {
    --A;
    --B;

    char c1 = *A;
    char c2 = *B;

    /* If we've got a match (and N's don't count) ... */
    if (validSymbol[c1] &&
        validSymbol[c2] &&
        IUPACidentity[c1][c2]) {
      good_run_len++; 

      /* If we've gone long enough, erase our mismatch record */
      if (good_run_len == DEF_MIN_BLOCK_SEP) {
	num_recent_mismatches = 0;
      }

      //  If we're in the middle of a long good run, add the character to
      //  the match
      //
      if (good_run_len > DEF_MIN_END_RUN_LEN) {
	m->extendBackward(1);
      } else if (good_run_len == DEF_MIN_END_RUN_LEN) {
        //  else if we just made the minimum extension length, add
        //  all of the pending characters.
        //
	m->extendBackward(num_pending + 1);
	num_pending = 0;
      } else {
        //  Otherwise, this character is pending.  However, still do output
        //  if we're run out of sequence.
        //
	num_pending++;
	if ((A.getPosition() == min_start_pos) || (B.getPosition() == 0)) {
	  m->extendBackward(num_pending);
	}
      }
    } else {
      //  Else this character is a mismatch.
      //
      num_pending++;
      good_run_len = 0;
      num_recent_mismatches++;
      if (num_recent_mismatches > DEF_MAX_MM_BLOCK)
	break;
    }
  }
}


static bool
can_reach_nearby_match(MEMatch *src, MEMatch *dest) {

  if (DEF_VERBOSE >= 3)
    fprintf(stderr, "can_reach_nearby_match()\n");

  int num_mismatch = 0;

  if (dest->pos1() - (src->pos1() + src->len()) > (SeqOff_t) DEF_MAX_NBR_SEP)
    return false;

  FastAAccessor   A(src->s1(), false);
  FastAAccessor   B(src->s2(), src->isReversed());

  B.setReverseComplementRange(src->pos2(), src->len());

  A.setPosition(src->pos1() + src->len());
  B.setPosition(src->pos2() + src->len());

  while ((num_mismatch <= DEF_MAX_NBR_PATH_MM) &&
         (A.getPosition() < dest->pos1())) {
    char c1 = *A;
    char c2 = *B;

    if (validSymbol[c1] &&
        validSymbol[c2] &&
        !IUPACidentity[c1][c2])
      num_mismatch++;

    ++A;
    ++B;
  }

  return (num_mismatch <= DEF_MAX_NBR_PATH_MM);
}


/* Stops and returns true if we hit the next match */
static bool
extend_match_forward(vector<MEMatch *>& matches, u32bit midx, MEMatch *target) {

  if (DEF_VERBOSE >= 3)
    fprintf(stderr, "extend_match_forward()\n");

  MEMatch *m = matches[midx];
  int num_recent_mismatches = 0;
  SeqOff_t num_pending = 0;

  int             good_run_len = (int) m->len();

  FastAAccessor   A(m->s1(), false);
  FastAAccessor   B(m->s2(), m->isReversed());

  B.setReverseComplementRange(m->pos2(), m->len());

  A.setPosition(m->pos1() + m->len());
  B.setPosition(m->pos2() + m->len());

  while (A.isValid() && B.isValid()) {
    char c1 = *A;
    char c2 = *B;

    // If there's a match ...
    //
    if (validSymbol[c1] &&
        validSymbol[c2] &&
        IUPACidentity[c1][c2]) {
      good_run_len++;

      // Pass Go and collect $200
      //
      if (good_run_len == DEF_MIN_BLOCK_SEP)
	num_recent_mismatches = 0;

      
      //  Not enough good characters yet, so output is pending ...
      //
      if (good_run_len < DEF_MIN_END_RUN_LEN) {
	num_pending++;	

	//  If we've got a short good run but have hit the end of
	//  a sequence, do extension.
        //
	if (((A.getPosition() + 1) == A.getLength()) || ((B.getPosition() + 1) == B.getLength())) 
	  m->extendForward(num_pending);	

      } else if (good_run_len == DEF_MIN_END_RUN_LEN) {
        //  Else if we've just made the minimum good end run length,
        //  commit to extending.
        //
	m->extendForward(num_pending + 1);
	num_pending = 0;
      } else if (good_run_len > DEF_MIN_END_RUN_LEN) {
        //  Otherwise, if we're above the minimum good run length,
	//  extend by one character.
	m->extendForward(1);
      }

      //  If we've run into (and possibly over) another seed match,
      //  consume it and return so the main loop can restart
      //
      if (target && (m->canMergeWith(*target))) {
	m->consume(*target);
	return true;
      }
    } else {
      // Otherwise, process mismatch ...
      num_pending++;
      good_run_len = 0;
      num_recent_mismatches++;
      if (num_recent_mismatches > DEF_MAX_MM_BLOCK)
	break;
    }

    ++A;
    ++B;
  }
  
  return false;
}


static u32bit
extend_matches_on_diagonal(vector<MEMatch *>& matches, u32bit diag_start) {

  if (DEF_VERBOSE >= 3)
    fprintf(stderr, "extend_matches_on_diagonal()\n");

  u64bit diag_id = matches[diag_start]->diagID();
  u32bit idx, next_idx;
  MEMatch *m, *prev_m = NULL, *next_m = NULL;

  //  Back extend each match as far as possible (but never over the
  //  preceding match
  //
  for (idx = diag_start; 
       (idx < matches.size()) && (matches[idx]->diagID() == diag_id);
       ++idx) {
    curmatch_G = m = matches[idx];
    if (DEF_VERBOSE >= 2) {
      cout << "Before back extension:\n";
      m->textDump(cout, true, 5);    
    }

    if (prev_m)
      extend_match_backward(matches, idx, prev_m->pos1() + prev_m->len());
    else
      extend_match_backward(matches, idx, 0);

    if (DEF_VERBOSE >= 2) {
      cout << "After back extension:\n";
      m->textDump(cout, true, 5);    
    }

    if (DEF_VERBOSE >= 2)
      cout << "\n==============\n\n";
    
    prev_m = m;
  }

  /* Now forward extend each match ... */


  idx = diag_start; 
  while ((idx < matches.size()) && 
	 (matches[idx]->diagID() == diag_id)) {
    /* Set current and next match pointers */
    if (matches[idx]->isDeleted()) {
      idx++;
      continue;
    }
    
    curmatch_G = m = matches[idx];
    next_m = NULL;
    next_idx = idx + 1;
    while ((next_idx < matches.size()) && 
	   (matches[next_idx]->diagID() == diag_id))
      if (!matches[next_idx]->isDeleted()) {
	next_m = matches[next_idx];
	break;
      }
      else
	next_idx++;

    /* First, try to reach the next match with the simple "maximum of
       k mismatches" rule. If we made it, consume the next match and
       start the loop again with the same match (now extended) */
    if (next_m && can_reach_nearby_match(m, next_m)) {
      m->consume(*next_m);
      next_m->setDeleted();
      if (DEF_VERBOSE >= 2) {
	cout << "Extended through next match via neighbor search:\n";
	m->textDump(cout, true, 5);    
      }
      continue;
    }

    /* Otherwise, try to make it to the next match with the character-at-
       a-time extension rules.  If we make it, restart the loop with the
       same match (now extended).  Otherwise, trim the extended match
       as necessary and move on to the next match. */
    if (extend_match_forward(matches, idx, next_m)) {
      next_m->setDeleted();
      if (DEF_VERBOSE >= 2) {
	cout << "Extended through next match via forward extension:\n";
	m->textDump(cout, true, 5);
      }
      continue;
    } else if (DEF_VERBOSE >= 2) {
      cout << "Failed to make next match.  Final extended version:\n";
      m->textDump(cout, true, 5);
    }
    
    //  Didn't make it, so trim and move on
    //
    if (trim_to_pct(matches, idx, DEF_MIN_IDENTITY)) {
      if (DEF_VERBOSE >= 2) {
	cout << "After trimming:\n";
	m->textDump(cout, true, 5);
      }
    } else {
      if (DEF_VERBOSE >= 2) {
        cout << "No trimming done.\n";
      }
    }

    if (DEF_VERBOSE >= 2)
      cout << "\n==============\n\n";

    ++idx;
  }
  
  curmatch_G = NULL;
  return idx;
}




int
main(int argc, char *argv[]) {
  int cmd_opt;
  char *args = "vE:B:S:I:P:D:x:y:h";
  bool cmd_err_flag = false;
  string prefix1, prefix2;
  string infile, outfile;
  time_t start_t, end_t;
  u32bit dp;

  //
  // Process arguments
  //

  while (!cmd_err_flag && ((cmd_opt = getopt(argc, argv, args)) != -1)) {
    switch (cmd_opt) {

    case 'v':
      DEF_VERBOSE++;
      break;

    case 'E':
      DEF_MIN_END_RUN_LEN = atoi(optarg);
      break;

    case 'B':
      DEF_MAX_MM_BLOCK = atoi(optarg);
      break;

    case 'S':
      DEF_MIN_BLOCK_SEP = atoi(optarg);
      break;

    case 'I':
      DEF_MIN_IDENTITY = (float) atof(optarg);
      break;

    case 'P':
      DEF_MAX_NBR_SEP = atoi(optarg);
      break;

    case 'D':
      DEF_MAX_NBR_PATH_MM = atoi(optarg);
      break;
      
    case 'h':
      print_usage();
      exit(0);

    default:
      cerr << "Error - unrecongized option: " << cmd_opt << endl;
      cmd_err_flag = true;
    }
  }

  if (cmd_err_flag) {
    print_usage();
    exit(1);
  }

  /* Input file */
  if (optind == argc) {
    cerr << "Error: input ATAC file required.\n";
    print_usage();
    exit(1);
  }

  infile = argv[optind++];
  
  /* Output file */
  if (optind == argc) {
    cerr << "Error: output ATAC file required.\n";
    print_usage();
    exit(1);
  }

  outfile = argv[optind++];

  // Log parameter values
  log_parameters();
  cerr << " * Input file = " << infile << '\n';
  cerr << " * Output file = " << outfile << '\n';

  //  XXX:  This could probably be moved into the option parsing, removing
  //  the "temporary" globals
  //
  MatchExtenderParameters p;
  p.minEndRunLen  = DEF_MIN_END_RUN_LEN;
  p.maxMMBlock    = DEF_MAX_MM_BLOCK;
  p.minBlockSep   = DEF_MIN_BLOCK_SEP;
  p.minIdentity   = DEF_MIN_IDENTITY;
  p.maxNbrSep     = DEF_MAX_NBR_SEP;
  p.maxNbrPathMM  = DEF_MAX_NBR_PATH_MM;
  p.beVerbose     = DEF_VERBOSE;

  //
  // Set up the ATAC dataset
  //

  cerr << " * Make new ATACreader ..." << endl;
  ATACreader atac_me(&p, infile, outfile);
  cerr << " * Made." << endl;

  cerr << " * First parsing pass to extract global variables ..." << endl;
  if (atac_me.processPreamble() == false) {
    fprintf(stderr, "Failed to process the preamble.\n");
    exit(1);
  }
  cerr << " * Parsed." << endl;

  prefix1 = atac_me.filePrefix1();
  prefix2 = atac_me.filePrefix2();

  atac_me.setFastaFiles(prefix1, prefix2);

  //
  // Main loop: process each set of hits sharing the same two sequences
  // a group at a time
  //

  string axis1, axis2;
  vector<MEMatch *> fwd_matches(1000000); // Guess a good initial capacity
  vector<MEMatch *> rev_matches(1000000);
  FastASequenceInCore *s1 = 0L;
  FastASequenceInCore *s2 = 0L;
  u32bit diag_start;

  // XXX removed to reduce clutter
  //install_sig_handler();

  start_t = time(NULL);
  while (!atac_me.endOfMatches()) {

    if (DEF_VERBOSE > 0)
      cerr << " * Reading next batch of matches and sequence information ..."
	   << endl;
    assert(atac_me.getNextMatchBatches(fwd_matches, rev_matches, s1, s2, 
				       axis1, axis2));

    fprintf(stderr, " * Processing matches between %s and %s\r", axis1.c_str(), axis2.c_str());
    fflush(stderr);

    if (DEF_VERBOSE > 0)
      cerr << " * Sorting matches" << endl;

    if (fwd_matches.size() > 0)
      sort(fwd_matches.begin(), fwd_matches.end(), MEMatchPtrCmp());
    if (rev_matches.size() > 0)
      sort(rev_matches.begin(), rev_matches.end(), MEMatchPtrCmp());

    if (DEF_VERBOSE > 0)
      cerr << " * Extending matches" << endl;

    diag_start = 0;
    while (diag_start < fwd_matches.size()) {
      diag_start = extend_matches_on_diagonal(fwd_matches, diag_start);
    }

    diag_start = 0;
    while (diag_start < rev_matches.size()) {
      diag_start = extend_matches_on_diagonal(rev_matches, diag_start);
    }

    if (DEF_VERBOSE > 0)
      cerr << " * Writing extended matches" << endl;

    atac_me.writeMatchBatch(fwd_matches, axis1, axis2);    
    atac_me.writeMatchBatch(rev_matches, axis1, axis2);

    for (dp = 0; dp < fwd_matches.size(); ++dp)
      MEMatch::freeMatch(fwd_matches[dp]);
    for (dp = 0; dp < rev_matches.size(); ++dp)
      MEMatch::freeMatch(rev_matches[dp]);

    //  If we're using the FastACache, we don't own the sequences, and
    //  we don't need to destroy them.
    //
    s1 = 0L;
    s2 = 0L;
  }

  fprintf(stderr, "\n");
  
  end_t = time(NULL);
  cerr << " * Run took " << ((end_t - start_t) / 60) << "m "
       << ((end_t - start_t) % 60) << "s\n" 
       << " * Read " << atac_me.numRead() << " matches\n"
       << " * Wrote " << atac_me.numWritten() << " matches\n"
       << " *" << endl;
  
  cerr << " * Normal exit.  Bye!\n";

  return 0;
}

