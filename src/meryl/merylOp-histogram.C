
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "meryl.H"

void
merylOperation::reportHistogram(void) {

  if (_inputs.size() > 1)
    fprintf(stderr, "ERROR: told to dump a histogram for more than one input!\n"), exit(1);

  if (_inputs[0]->_operation)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  if (_inputs[0]->_sequence)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  //  Tell the stream to load and return the histogram.

  kmerCountStatistics *stats = _inputs[0]->_stream->stats();

  //  Now just dump it.

  for (uint32 ii=0; ii<stats->histogramLength(); ii++)
    fprintf(stdout, F_U64 "\t" F_U64 "\n",
            stats->histogramValue(ii),
            stats->histogramOccurrences(ii));
}



void
merylOperation::reportStatistics(void) {

  if (_inputs.size() > 1)
    fprintf(stderr, "ERROR: told to dump a histogram for more than one input!\n"), exit(1);

  if (_inputs[0]->_operation)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  if (_inputs[0]->_sequence)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  //  Tell the stream to load and return the histogram.

  kmerCountStatistics *stats = _inputs[0]->_stream->stats();

  //  Now just dump it.

  uint64  nUniverse = uint64MASK(kmer::merSize() * 2) + 1;
  uint64  sDistinct = 0;
  uint64  sTotal    = 0;

  fprintf(stdout, "Number of %u-mers that are:\n", kmer::merSize());
  fprintf(stdout, "  unique   %20" F_U64P "  (exactly one instance of the kmer is in the input)\n", stats->numUnique());
  fprintf(stdout, "  distinct %20" F_U64P "  (non-redundant kmer sequences in the input)\n", stats->numDistinct());
  fprintf(stdout, "  present  %20" F_U64P "  (...)\n", stats->numTotal());
  fprintf(stdout, "  missing  %20" F_U64P "  (non-redundant kmer sequences not in the input)\n", nUniverse - stats->numDistinct());
  fprintf(stdout, "\n");
  fprintf(stdout, "             number of   cumulative   cumulative     presence\n");
  fprintf(stdout, "              distinct     fraction     fraction   in dataset\n");
  fprintf(stdout, "frequency        kmers     distinct        total       (1e-6)\n");
  fprintf(stdout, "--------- ------------ ------------ ------------ ------------\n");

  for (uint32 ii=0; ii<stats->histogramLength(); ii++) {
    uint64  value = stats->histogramValue(ii);
    uint64  occur = stats->histogramOccurrences(ii);

    sDistinct  += occur;
    sTotal     += occur * value;

    fprintf(stdout, "%9" F_U64P " %12" F_U64P " %12.4f %12.4f %12.6f\n",
            value,
            occur,
            (double)sDistinct / stats->numDistinct(),
            (double)sTotal    / stats->numTotal(),
            (double)value     / stats->numTotal() * 1000000.0);
  }

  assert(sDistinct == stats->numDistinct());
  assert(sTotal    == stats->numTotal());
}
