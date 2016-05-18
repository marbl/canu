
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
 *    Sergey Koren beginning on 2016-MAY-16
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "ovStore.H"
#include "splitToWords.H"

#include "AS_UTL_decodeRange.H"
#include "stddev.H"

#include <vector>
#include <algorithm>
#include <map>

using namespace std;

int
main(int argc, char **argv) {
  char           *scoreFileName    = NULL;
  uint32         deviations = 6;
  float          mass=0.98;
  bool           isOvl=false;

  argc = AS_configure(argc, argv);

  int32     arg = 1;
  int32     err = 0;
  while (arg < argc) {
    if (strcmp(argv[arg], "-S") == 0) {
      scoreFileName = argv[++arg];

    } else if (strcmp(argv[arg], "-d") == 0) {
       deviations = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-m") == 0) {
       mass = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
       isOvl=true;

    } else {
      fprintf(stderr, "ERROR:  invalid arg '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }

  if (scoreFileName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\n");

    exit(1);
  }

  errno = 0;
  FILE     *scoreFile   = (scoreFileName == NULL) ? NULL : (scoreFileName[0] == '-' ? stdin : fopen(scoreFileName, "r"));
  if (errno)
    fprintf(stderr, "ERROR: failed to open '%s' for reading: %s\n", scoreFileName, strerror(errno)), exit(1);

  // read the file and store best hits
  char  ovStr[1024];
  ovOverlap   ov(NULL);
  map<uint32, uint32> readToLength;
  map<uint32, double> readToIdy;
  double mean, median, stddev, mad;
  mean = median = stddev = mad = 0.0;

  while (fgets(ovStr, 1024, scoreFile) != NULL) {
      splitToWords  W(ovStr);

      if (isOvl) {
         ov.a_iid = W(0);
         ov.b_iid = W(1);
         if (ov.a_iid == ov.b_iid)
            continue;
         ov.dat.ovl.ahg5 = W(4);
         ov.dat.ovl.ahg3 = W(6);
         ov.dat.ovl.bhg5 = W(6);
         ov.dat.ovl.bhg3 = W(7);
         ov.span(W(3));
         ov.erate(atof(W[8]));
         ov.flipped(W[3][0] == 'I' ? true : false);

      } else {
         ov.a_iid = W(0);
         ov.b_iid = W(1);

         if (ov.a_iid == ov.b_iid)
            continue;

         assert(W[4][0] == '0');

         ov.dat.ovl.ahg5 = W(5);
         ov.dat.ovl.ahg3 = W(7) - W(6);

         if (W[8][0] == '0') {
            ov.dat.ovl.bhg5 = W(9);
            ov.dat.ovl.bhg3 = W(11) - W(10);
            ov.flipped(false);
         } else {
            ov.dat.ovl.bhg3 = W(9);
            ov.dat.ovl.bhg5 = W(11) - W(10);
            ov.flipped(true);
         }
         ov.erate(atof(W[2]));
         ov.span(W(10)-W(9));
      }

      if (ov.erate() == 0.0)
         ov.erate(0.01); // round up when we can't estimate accurately

      if (readToLength.find(ov.b_iid) == readToLength.end() || readToLength[ov.b_iid] < ov.span()) {
         readToLength[ov.b_iid] = ov.span();
         readToIdy[ov.b_iid] = ov.erate();
      }
  }
  fclose(scoreFile);

  stdDev<double>  edgeStats;

  //  Find the overlap for every best edge.

  double  *absdev    = new double [readToLength.size() + 1];
  double  *erates    = new double [readToLength.size() + 1];
  uint32   eratesLen = 0;


  for (map<uint32, double>::iterator it=readToIdy.begin(); it != readToIdy.end(); ++it) {
     edgeStats.insert(erates[eratesLen++] = it->second);
  }

  mean   = edgeStats.mean();
  stddev = edgeStats.stddev();

  fprintf(stderr, "with %u points - mean %f stddev %f - would use overlaps below %f fraction error\n", edgeStats.size(), mean, stddev, mean + deviations * stddev);

  //  Find the median and absolute deviations.

  sort(erates, erates+eratesLen);

  median = erates[ eratesLen / 2 ];

  double massCutoff = 0;
  uint32 totalBelow = 0;
  for (uint32 ii=0; ii<eratesLen/2; ii++) {
    absdev[ii] = median - erates[ii];
    if ((double)totalBelow / eratesLen < mass) {
       massCutoff = erates[ii];
       totalBelow++;
    }
  }

  for (uint32 ii=eratesLen/2; ii<eratesLen; ii++) {
    absdev[ii] = erates[ii] - median;
    if ((double)totalBelow / eratesLen < mass) {
       massCutoff = erates[ii];
       totalBelow++;
    }
  }

  sort(absdev, absdev+eratesLen);

  assert(absdev[0] >= 0.0);

  mad    = absdev[eratesLen/2];

  delete [] absdev;
  delete [] erates;

  fprintf(stderr, "with %u points - median %f mad %f - would use overlaps below %f fraction error\n",
           edgeStats.size(), median, mad, median + deviations * 1.4826 * mad);

   fprintf(stderr, "with %u points - mass of %d is below %f\n", edgeStats.size(), totalBelow, massCutoff);

  if (scoreFile)
    fclose(scoreFile);

  fprintf(stdout, "%.3f\n",  massCutoff /* median + deviations * 1.4826 * mad*/);
  exit(0);
}
