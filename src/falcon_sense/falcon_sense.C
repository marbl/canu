
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
 *    Sergey Koren beginning on 2016-FEB-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "gkStore.H"
#include "splitToWords.H"
#include "AS_UTL_fasta.H"

#include "falcon.H"

#include <omp.h>
#include <vector>
#include <string>

using namespace std;

int
main (int argc, char **argv) {
  uint32 threads = 0;
  uint32 min_cov = 4;
  uint32 min_len = 500;
  uint32 min_ovl_len = 500;
  double min_idy = 0.5;
  uint32 K = 8;
  uint32 max_read_len = AS_MAX_READLEN;

  argc = AS_configure(argc, argv);


  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "--n_core") == 0) {
      threads = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "--min_cov") == 0) {
      min_cov = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "--min_idt") == 0) {
      min_idy = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "--min_len") == 0) {
       min_len = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "--min_ovl_len") == 0) {
       min_ovl_len = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "--max_read_len") == 0) {
       max_read_len = atoi(argv[++arg]);
       if (max_read_len <= 0 || max_read_len > AS_MAX_READLEN) {
          max_read_len = AS_MAX_READLEN;
       }

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }

  if (err) {
     fprintf(stderr, "Invalid usage");
     exit(1);
  }

  if (threads > 0) {
    omp_set_num_threads(threads);
  } else {
    omp_set_num_threads(omp_get_max_threads());
  }

  // read in a loop and get consensus of each read
  vector<string> seqs;

  string seed;
  char *A = new char[AS_MAX_READLEN * 2];

  while (true) {
    fgets(A, AS_MAX_READLEN * 2, stdin);
    splitToWords W(A);

    if (W[0][0] == '+') {
       uint32 splitSeqID = 0;
       FConsensus::consensus_data *consensus_data_ptr = FConsensus::generate_consensus( seqs, min_cov, K, min_idy, min_ovl_len, max_read_len );
       char * split = strtok(consensus_data_ptr->sequence, "acgt");
       while (split != NULL) {
          if (strlen(split) > min_len) {
             AS_UTL_writeFastA(stdout, split, strlen(split), 60, ">%s_%d\n", seed.c_str(), splitSeqID);
             splitSeqID++;
          }
          split = strtok(NULL, "acgt");
       }
       FConsensus::free_consensus_data( consensus_data_ptr );
       seqs.clear();
       seed.clear();
    } else if (W[0][0] == '-') {
       break;
    } else {
       if (seed.length() == 0) {
          seed = W[0];
          seqs.push_back(string(W[1]));
       }
       else if (strlen(W[1]) > min_ovl_len) {
          seqs.push_back(string(W[1]));
       }
    }
  }

  delete[] A;
}
