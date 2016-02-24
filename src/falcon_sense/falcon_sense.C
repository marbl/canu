#include "AS_global.H"
#include "gkStore.H"
#include "splitToWords.H"
#include "AS_UTL_fasta.H"

#include "falcon.H"

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
  bool doSort = false;
  uint32 K = 8;

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

    } else if (strcmp(argv[arg], "--sort") == 0) {
       doSort = true;

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

  // since we're align-level parallelism we won't guarantee alignment order unless there's only one thread so enforce that if a user really wants sort maintained
  if (doSort) 
     threads = 1;

  if (threads > 0) {
    omp_set_num_threads(threads);
  } else {
    omp_set_num_threads(omp_get_max_threads());
  }

  // read in a loop and get consensus of each read
  // when we're sorting we store results in a multimap so they're sorted by length, otherwise in a vector by order of appearance
  multimap<uint32, string> sortedSeqs;
  vector<string> seqs;

  string seed;
  string seedSeq;
  char *A = new char[AS_MAX_READLEN * 2];

  while (true) {
    fgets(A, AS_MAX_READLEN * 2, stdin);
    splitToWords W(A);

    if (W[0][0] == '+') {
       if (doSort) {
          // take longest in order except for seed
          seqs.push_back(seedSeq);
          for (multimap<uint32, string>::reverse_iterator it = sortedSeqs.rbegin(); it != sortedSeqs.rend(); it++)
             seqs.push_back(it->second);
       }
       uint32 splitSeqID = 0;
       FConsensus::consensus_data *consensus_data_ptr = FConsensus::generate_consensus( seqs, min_cov, K, min_idy, min_ovl_len );
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
       sortedSeqs.clear();
       seed.clear();
    } else if (W[0][0] == '-') {
       break;
    } else {
       if (seed.length() == 0) {
          seed = W[0];
          if (doSort)
             seedSeq = string(W[1]);
          else
             seqs.push_back(string(W[1]));
       }
       else if (strlen(W[1]) > min_ovl_len) {
          if (doSort)
             sortedSeqs.insert(pair<uint32, string>(strlen(W[1]), string(W[1])));
          else
             seqs.push_back(string(W[1]));
       }
    } 
  }

  delete[] A;
}
