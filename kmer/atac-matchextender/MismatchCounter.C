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
#include <fstream>
#include <vector>
#include <algorithm>

#include "MEMatch.H"
#include "atacreader.H"

//
// GLOBAL PARAMETERS
//

MEMatch *curmatch_G = NULL;   // The match currently being operated on

int   DEF_VERBOSE           = 0;


static void
grow_histogram(vector<unsigned long>& hist, int max_val) {
  unsigned int cur_size = hist.size();
  unsigned int new_size = (unsigned int) max_val + 1;
  unsigned int i;

  hist.resize(new_size);
  for (i = cur_size; i < new_size; ++i)
    hist[i] = 0;
}


static void
write_histogram(vector<unsigned long>& hist, char *file) {
  unsigned int i, max_nz;
  ofstream os;

  os.open(file);
  if (!os) {
    cerr << "Error: Could not open " << file << " for output" << endl;
    return;
  }

  max_nz = hist.size() - 1;
  while ((max_nz != 0) && (hist[max_nz] == 0))
    max_nz--;

  for (i = 1; i <= max_nz; ++i)
    os << i << '\t' << hist[i] << '\n';

  os.close();
}


static void
print_usage(void) {
  cerr << "Usage: MismatchCounter [options] <infile>\n"
       << "   <infile>    - Input ATAC file with parentless exact matches\n"
       << "Options:\n"
       << "   -x <string> - Prefix of axis1 (read from infile otherwise)\n"
       << "   -y <string> - Prefix of axis2 (read from infile otherwise)\n"
       << "   -v          - Verbose (debugging) mode.  Lots of output!\n"
       << "   -h          - Print this help message\n"
       << "Parameters:\n"
       << "Param    Meaning                                       Def. value\n"
       << "-----------------------------------------------------------------\n"
       << " -H <file> Set output histogram file to <file>           (null)"
       << endl;
}


static void
log_parameters(void) {
  cerr << " * Parameter values: (none)\n";
}


static unsigned long
count_mismatch_blocks(MEMatch *m, vector<unsigned long>& hist) {
  char c1 = 0;
  char c2 = 0;
  SeqOff_t l = 0;
  unsigned long num_mismatch = 0;
  unsigned int cur_blk_size = 0;

  char    *s1  = m->s1()->sequence();
  SeqOff_t s1p = m->pos1();

  char    *s2  = m->s2()->sequence();
  SeqOff_t s2p = m->pos2();
  SeqOff_t s2r = m->rcPos2() + s2p;

  //
  //  XXX: Just do one case for forward and one case for complement.
  //  Complement the base in real time.  This saves us from storing a
  //  reverse-complelemt copy of the chromosome
  //

  //if (m->isReversed()) {
  //  s2  = m->s2rc()->sequence();
  //  s2p = m->rcPos2();
  //}

  //  XXX: do we need to reversecomplement the coordinates or just
  //  decrement the positions?

  for (l=0; l<m->len(); ++l, s1p++, s2p++) {
    c1 = s1[s1p];

    if (m->isReversed())
      c2 = complementSymbol[s2[s2r - s2p]];
    else
      c2 = s2[s2p];

    //  validsymbol[] is true if ACGTacgt, then we can use
    //  IUPACidentity to test if A==a.
    //
    if (validSymbol[c1] &&
        validSymbol[c2] &&
        IUPACidentity[c1][c2]) {
      if (cur_blk_size == 0)
	num_mismatch++;
      cur_blk_size++;
    } else if (cur_blk_size > 0) {
      if (cur_blk_size >= hist.size())
	grow_histogram(hist, cur_blk_size);
      hist[cur_blk_size]++;
      cur_blk_size = 0;
    }
  }

  return num_mismatch;
}


int
main(int argc, char *argv[]) {
  int                     cmd_opt;
  char                   *args = "H:vx:y:h";
  bool                    cmd_err_flag = false;
  string                  prefix1;
  string                  prefix2;
  string                  infile;
  bool                    override_pre1 = false;
  bool                    override_pre2 = false;
  time_t                  start_t;
  time_t                  end_t;
  u32bit                  dp;
  unsigned long           result = 0;
  vector<unsigned long>   hist(20);
  char                   *hist_file = NULL;

  //
  // Process arguments
  //

  while (!cmd_err_flag && ((cmd_opt = getopt(argc, argv, args)) != -1)) {
    switch (cmd_opt) {
      
    case 'H':
      hist_file = optarg;
      break;

    case 'v':
      DEF_VERBOSE++;
      break;

    case 'x':
      prefix1 = optarg;
      override_pre1 = true;
      cerr << " * Prefix 1 = " << prefix1 << endl;
      break;

    case 'y':
      prefix2 = optarg;
      override_pre2 = true;
      cerr << " * Prefix 2 = " << prefix2 << endl;
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
  
  // Log parameter values
  log_parameters();
  cerr << " * Input file = " << infile << '\n';

  MatchExtenderParameters p;
  
  //
  // Set up the ATAC dataset
  //

  for (unsigned int i = 0; i < hist.size(); ++i)
    hist[i] = 0;

  ATACreader    atac_me(&p, infile);

  cerr << " * First parsing pass to extract global variables ..." << endl;
  atac_me.processPreamble();

  if (!override_pre1) {
    prefix1 = atac_me.filePrefix1();
    cerr << " * Prefix 1 = " << prefix1 << endl;
 }

  if (!override_pre2) {
    prefix2 = atac_me.filePrefix2();
    cerr << " * Prefix 2 = " << prefix2 << endl;
  }

  atac_me.setFastaFiles(prefix1 + ".fasta", prefix2 + ".fasta");


  //
  // Main loop: process each set of hits sharing the same two sequences
  // a group at a time
  //

  string axis1, axis2;
  vector<MEMatch *> fwd_matches(1000000); // Guess a good initial capacity
  vector<MEMatch *> rev_matches(1000000);
  FastASequenceInCore *s1 = 0L;
  FastASequenceInCore *s2 = 0L;
  unsigned int mi;

  //install_sig_handler();

  start_t = time(NULL);
  while (!atac_me.endOfMatches()) {

    if (DEF_VERBOSE > 0)
      cerr << " * Reading next batch of matches and sequence information ..."
	   << endl;
    assert(atac_me.getNextMatchBatches(fwd_matches, rev_matches, s1, s2, 
				       axis1, axis2));

    cerr << " * Processing matches between " << axis1 << " and " 
	 << axis2 << endl;

    for (mi = 0; mi < fwd_matches.size(); ++mi)
      result += count_mismatch_blocks(fwd_matches[mi], hist);
    
    for (mi = 0; mi < rev_matches.size(); ++mi)
      result += count_mismatch_blocks(rev_matches[mi], hist);

    for (dp = 0; dp < fwd_matches.size(); ++dp)
      MEMatch::freeMatch(fwd_matches[dp]);
    for (dp = 0; dp < rev_matches.size(); ++dp)
      MEMatch::freeMatch(rev_matches[dp]);

    s1 = 0L;
    s2 = 0L;
  }

  if (hist_file)
    write_histogram(hist, hist_file);
  
  end_t = time(NULL);
  cerr << " * Run took " << ((end_t - start_t) / 60) << "m "
       << ((end_t - start_t) % 60) << "s\n" 
       << " * Read " << atac_me.numRead() << " matches\n"
       << " *" << endl;
  
  cout << result << endl;

  cerr << " * Normal exit.  Bye!\n";

  return 0;
}

