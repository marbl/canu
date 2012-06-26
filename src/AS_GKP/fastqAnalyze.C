
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2010, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: fastqAnalyze.C,v 1.4 2012-06-26 00:09:19 brianwalenz Exp $";

#include "AS_global.h"


int
main(int argc, char **argv) {
  char   *inName = NULL;
  char   *otName = NULL;

  bool    originalIsSolexa   = false;
  bool    originalIsIllumina = false;
  bool    originalIsSanger   = false;

  bool    correctToSanger    = false;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      otName          = argv[++arg];
      correctToSanger = true;

    } else if (strcmp(argv[arg], "-solexa") == 0) {
      originalIsSolexa = true;

    } else if (strcmp(argv[arg], "-illumina") == 0) {
      originalIsIllumina = true;

    } else if (strcmp(argv[arg], "-sanger") == 0) {
      originalIsSanger = true;

    } else if (inName == NULL) {
      inName = argv[arg];

    } else {
      err++;
    }

    arg++;
  }
  if ((err) || (inName == NULL)) {
    fprintf(stderr, "usage: %s [options] [-o output.fastq] input.fastq\n", argv[0]);
    fprintf(stderr, "  If no options are given, input.fastq is analyzed and a best guess for the\n");
    fprintf(stderr, "  QC encoding is output.  Otherwise, the QV encoding is converted to Sanger-style\n");
    fprintf(stderr, "  using this guess.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  In some cases, the encoding cannot be determined.  When this occurs, no guess is\n");
    fprintf(stderr, "  output.  For conversion, you can force the input QV type with:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -solexa     input QV is solexa\n");
    fprintf(stderr, "  -illumina   input QV is illumina\n");
    fprintf(stderr, "  -sanger     input QV is sanger\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o          sanger-style-output.fastq\n");
    exit(1);
  }

  uint32 numValid  = 5;
  uint32 numTrials = 100000;

  char   A[1024];
  char   B[1024];
  char   C[1024];
  char   D[1024];

  if ((originalIsSolexa   == false) &&
      (originalIsIllumina == false) &&
      (originalIsSanger   == false)) {
    errno = 0;
    FILE *F = fopen(inName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName, strerror(errno)), exit(1);

    //  Initially, it could be any of these.
    //
    bool   isNotSanger    = false;
    bool   isNotSolexa    = false;
    bool   isNotIllumina3 = false;
    bool   isNotIllumina5 = false;
    bool   isNotIllumina8 = false;

    uint32 qvCounts[256]  = {0};

    while ((!feof(F)) &&
           (numValid != 1) &&
           (numTrials > 0)) {
      fgets(A, 1024, F);
      fgets(B, 1024, F);
      fgets(C, 1024, F);
      fgets(D, 1024, F);  chomp(D);

      for (uint32 x=0; D[x] != 0; x++) {
        if (D[x] < '!')  isNotSanger    = true;
        if (D[x] < ';')  isNotSolexa    = true;
        if (D[x] < '@')  isNotIllumina3 = true;  //  Illumina 1.3
        if (D[x] < 'B')  isNotIllumina5 = true;  //  Illumina 1.5
        if (D[x] < '!')  isNotIllumina8 = true;  //  Illumina 1.5

        if ('I' < D[x])  isNotSanger    = true;
        if ('h' < D[x])  isNotSolexa    = true;
        if ('h' < D[x])  isNotIllumina3 = true;  //  Illumina 1.3
        if ('h' < D[x])  isNotIllumina5 = true;  //  Illumina 1.5
        if ('J' < D[x])  isNotIllumina8 = true;  //  Illumina 1.5

        qvCounts[D[x]]++;

        //fprintf(stderr, "%d%d%d%d%d %c %d\n",
        //        isNotSanger, isNotSolexa, isNotIllumina3, isNotIllumina5, isNotIllumina8, D[x], x);
      }

      //fprintf(stderr, "%d%d%d%d%d '%s'\n",
      //        isNotSanger, isNotSolexa, isNotIllumina3, isNotIllumina5, isNotIllumina8, D);

      numValid = 0;

      if (isNotSanger    == false)  numValid++;
      if (isNotSolexa    == false)  numValid++;
      if (isNotIllumina3 == false)  numValid++;
      if (isNotIllumina5 == false)  numValid++;
      if (isNotIllumina8 == false)  numValid++;

      numTrials--;
    }

    fclose(F);

    fprintf(stdout, "%s --", inName);

    if (isNotSanger    == false)  fprintf(stdout, " SANGER");
    if (isNotSolexa    == false)  fprintf(stdout, " SOLEXA");
    if (isNotIllumina3 == false)  fprintf(stdout, " ILLUMINA_1.3+");
    if (isNotIllumina5 == false)  fprintf(stdout, " ILLUMINA_1.5+");
    if (isNotIllumina8 == false)  fprintf(stdout, " ILLUMINA_1.8+");
    if (numValid       == 0)      fprintf(stdout, " NO_VALID_ENCODING");

    fprintf(stdout, "\n");

    if (numValid == 0) {
      fprintf(stdout, "QV histogram:\n");

      for (uint32 c=0, i=0; i<12; i++) {
        fprintf(stdout, "%3d: ", i * 10);

        for (uint32 j=0; j<10; j++, c++)
          fprintf(stdout, " %8d/%c", qvCounts[c], isprint(c) ? c : '.');

        fprintf(stdout, "\n");
      }
    }

    if (isNotSanger    == false)  originalIsSanger   = true;
    if (isNotSolexa    == false)  originalIsSolexa   = true;
    if (isNotIllumina3 == false)  originalIsIllumina = true;
    if (isNotIllumina5 == false)  originalIsIllumina = true;
    if (isNotIllumina8 == false)  originalIsSanger   = true;
  }

  if (correctToSanger == false)
    exit(0);

  numValid = 0;

  if (originalIsSanger    == true)  numValid++;
  if (originalIsSolexa    == true)  numValid++;
  if (originalIsIllumina  == true)  numValid++;

  if (numValid == 0)
    fprintf(stderr, "No QV decision made.  No valid encoding found.  Specify a QV encoding to convert from.\n"), exit(0);

  if (numValid > 1)
    fprintf(stderr, "No QV decision made.  Multiple valid encodings found.  Specify a QV encoding to convert from.\n"), exit(0);


  if ((correctToSanger) &&
      (originalIsSanger == true))
    fprintf(stderr, "No QV changes needed; original is in sanger format already.\n"), exit(0);


  if ((correctToSanger) &&
      (originalIsSanger == false)) {
    errno = 0;
    FILE *F = fopen(inName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for reading: %s\n", inName, strerror(errno)), exit(1);

    errno = 0;
    FILE *O = fopen(otName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", otName, strerror(errno)), exit(1);

    while (!feof(F)) {
      fgets(A, 1024, F);
      fgets(B, 1024, F);
      fgets(C, 1024, F);
      fgets(D, 1024, F);  chomp(D);

      if (feof(F))
        break;

      for (uint32 x=0; D[x] != 0; x++) {
        if (originalIsSolexa) {
          double qs  = D[x] - '@';
          qs /= 10.0;
          qs  = 10.0 * log10(pow(10.0, qs) + 1);
          D[x] = lround(qs) + '0';
        }

        if (originalIsIllumina) {
          D[x] -= '@';
          D[x] += '!';
        }
      }

      fprintf(O, "%s%s%s%s\n", A, B, C, D);
    }

    fclose(F);
    fclose(O);
  }

  exit(0);
}
