
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

#ifndef TRIM_STAT_H
#define TRIM_STAT_H

#include "runtime.H"

class trimStat {
public:
  trimStat() {
    nReads = 0;
    nBases = 0;
  };

  trimStat &operator+=(uint32 bases) {
    nReads += 1;
    nBases += bases;

    histo.push_back(bases);

    return(*this);
  };

  void       generatePlots(char const *outputPrefix, char const *outputName, uint32 binwidth) {
    char  N[FILENAME_MAX];
    FILE *F;

    snprintf(N, FILENAME_MAX, "%s.%s.dat", outputPrefix, outputName);
    F = AS_UTL_openOutputFile(N);
    for (uint64 ii=0; ii<histo.size(); ii++)
      fprintf(F, F_U32"\n", histo[ii]);
    AS_UTL_closeFile(F, N);

    snprintf(N, FILENAME_MAX, "%s.%s.gp", outputPrefix, outputName);
    F = AS_UTL_openOutputFile(N);
    fprintf(F, "set title '%s'\n", outputName);
    fprintf(F, "set xlabel 'length, bin width = %u'\n", binwidth);
    fprintf(F, "set ylabel 'number'\n");
    fprintf(F, "\n");
    fprintf(F, "binwidth=%u\n", binwidth);
    fprintf(F, "set boxwidth binwidth\n");
    fprintf(F, "bin(x,width) = width*floor(x/width) + binwidth/2.0\n");
    fprintf(F, "\n");
    fprintf(F, "set terminal png size 1024,1024\n");
    fprintf(F, "set output '%s.%s.png'\n", outputPrefix, outputName);
    fprintf(F, "plot [] [0:] '%s.%s.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes title ''\n", outputPrefix, outputName);
    AS_UTL_closeFile(F, N);

    snprintf(N, FILENAME_MAX, "gnuplot %s.%s.gp > /dev/null 2>&1", outputPrefix, outputName);

    system(N);
  };

  uint32          nReads;
  uint64          nBases;

  vector<uint32>  histo;
};

#endif  //  TRIM_STAT_H
