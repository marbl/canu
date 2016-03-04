
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
 *    Brian P. Walenz beginning on 2016-JAN-19
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef TRIM_STAT_H
#define TRIM_STAT_H

#include "AS_global.H"

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

  void       generatePlots(char *outputPrefix, char *outputName, uint32 binwidth) {
    char  N[FILENAME_MAX];
    FILE *F;

    sprintf(N, "%s.%s.dat", outputPrefix, outputName);
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    for (uint64 ii=0; ii<histo.size(); ii++)
      fprintf(F, F_U32"\n", histo[ii]);
    fclose(F);

    sprintf(N, "%s.%s.gp", outputPrefix, outputName);
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", N, strerror(errno)), exit(1);
    fprintf(F, "set title '%s'\n", outputName);
    fprintf(F, "set xlabel 'length, bin width = %u'\n", binwidth);
    fprintf(F, "set ylabel 'number'\n");
    fprintf(F, "\n");
    fprintf(F, "binwidth=%u\n", binwidth);
    fprintf(F, "set boxwidth binwidth\n");
    fprintf(F, "bin(x,width) = width*floor(x/width) + binwidth/2.0\n");
    fprintf(F, "\n");
    fprintf(F, "set terminal png size 1024,1024\n");
    fprintf(F, "set output '%s.%s.lg.png'\n", outputPrefix, outputName);
    fprintf(F, "plot [] [0:] '%s.%s.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes title ''\n", outputPrefix, outputName);
    fprintf(F, "\n");
    fprintf(F, "set terminal png size 256,256\n");
    fprintf(F, "set output '%s.%s.sm.png'\n", outputPrefix, outputName);
    fprintf(F, "replot\n");
    fclose(F);

    sprintf(N, "gnuplot %s.%s.gp > /dev/null 2>&1", outputPrefix, outputName);

    system(N);
  };

  uint32          nReads;
  uint64          nBases;

  vector<uint32>  histo;
};

#endif  //  TRIM_STAT_H
