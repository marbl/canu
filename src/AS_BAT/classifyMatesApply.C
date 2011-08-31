
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: classifyMatesApply.C,v 1.2 2011-08-31 17:42:18 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include "classifyMates.H"

//  One per fragment, stores the results from multiple classification runs.
class classifyMatesSummary {
public:
  classifyMatesSummary() {
    isTested = 0;
    doDelete = 0;
    numNO    = 0;
    numPE    = 0;
    numMP    = 0;
  };

  void      setTested(void) { isTested = 1; };
  bool      getTested(void) { return(isTested > 0); };

  void      addNO(void) { if (numNO < 255)  numNO++; };
  void      addPE(void) { if (numPE < 255)  numPE++; };
  void      addMP(void) { if (numMP < 255)  numMP++; };

  uint8     getNO(void) { return(numNO); };
  uint8     getPE(void) { return(numPE); };
  uint8     getMP(void) { return(numMP); };

  void      setDelete(void)  { doDelete = 1; };
  bool      getDelete(void)  { return(doDelete == 1); };

  void      setUnMate(void)  { doUnMate = 1; };
  bool      getUnMate(void)  { return(doUnMate == 1); };

private:
  uint32    isTested   : 1;  //  We have attempted classification on this read
  uint32    doDelete   : 1;  //  Delete the read, it's junk
  uint32    doUnMate   : 1;  //  Remove the mate edge, it's mate read is junk
  uint32    numNO      : 8;
  uint32    numPE      : 8;
  uint32    numMP      : 8;
};



int
main(int argc, char **argv) {
  char      *gkpStoreName      = NULL;

  bool       justPrint         = false;

  char     **resultsNames      = new char * [argc];
  uint32     resultsNamesLen   = 0;

  char      *outputName        = NULL;
  FILE      *outputFile        = NULL;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-p") == 0) {
      justPrint = true;

    } else if (strcmp(argv[arg], "-r") == 0) {
      resultsNames[resultsNamesLen++] = argv[++arg];

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputName = argv[++arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }

    arg++;
  }
  if (gkpStoreName == 0L)
    err++;
  if (resultsNamesLen == 0)
    err++;
  if (outputName == 0L)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -G gkpStore -r resultsFile [-r resultsFile]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G gkpStore      Read fragments from here\n");
    fprintf(stderr, "  -r results       Read results from here; any number of -r options can be supplied\n");
    fprintf(stderr, "  -p               Also dump the results to stdout\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o editFile      Output gatekeeper edit file\n");
    fprintf(stderr, "\n");
    if (gkpStoreName == 0L)
      fprintf(stderr, "No gatekeeper store (-G) supplied.\n");
    if (resultsNamesLen == 0)
      fprintf(stderr, "No results files (-r) supplied.\n");
    if (outputName == 0L)
      fprintf(stderr, "No output file name (-o) supplied.\n");
    exit(1);
  }

  //  Open the output

  errno = 0;
  outputFile = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open output file '%s': %s\n", outputName, strerror(errno)), exit(1);

  //  Load gatekepeer

  gkStore               *gkpStore  = new gkStore(gkpStoreName, FALSE, FALSE);
  gkFragment             frg;

  uint32                 numFrags = gkpStore->gkStore_getNumFragments();
  uint32                 numLibs  = gkpStore->gkStore_getNumLibraries();

  //  We gloss over the distinction between PE and MP reads here.  It is writen assuming an MP
  //  library, with MP-as-innie and PE-as-outtie.  It works for any classification though.

  classifyMatesSummary  *results  = new classifyMatesSummary [numFrags + 1];
  classifyMatesResult    result;

  //  For each input file, load the results.

  for (uint32 i=0; i<resultsNamesLen; i++) {
    classifyMatesResultFile  *RF = new classifyMatesResultFile(resultsNames[i]);

    while (RF->read(result)) {
      AS_IID   f = result.readIID;
      AS_IID   m = result.mateIID;

      if (justPrint)
        result.print(stdout);

      results[f].setTested();
      results[m].setTested();

      //  If any single result says the read (or mate) is junk, delete it.

      if ((result.fragSpur     == true) ||
          (result.fragChimer   == true) ||
          (result.fragJunction == true)) {
        results[f].setDelete();
        results[f].setUnMate();
        results[m].setUnMate();
      }

      if ((result.mateSpur     == true) ||
          (result.mateChimer   == true) ||
          (result.mateJunction == true)) {
        results[m].setDelete();
        results[f].setUnMate();
        results[m].setUnMate();
      }

      //  Rembmer if any result failed to classify.  This is for the final accounting summary only.

      if (result.classified == false) {
        results[f].addNO();
        results[m].addNO();
        continue;
      }

      //  Otherwise, it must have been classified as either PE or MP.

      if (result.innie) {
        results[f].addMP();
        results[m].addMP();
      } else {
        results[f].addPE();
        results[m].addPE();
      }
    }

    delete RF;
  }

  //  Process the classifications.  If any are labeled both PE and MP, remove the mate.

  uint32  numDL = 0;
  uint32  numUM = 0;
  uint32  numNO = 0;
  uint32  numPE = 0;
  uint32  numMP = 0;
  uint32  numAB = 0;
  uint32  numIG = 0;

  for (uint32 i=1; i<numFrags; i++) {
    uint32  no = results[i].getNO();
    uint32  mp = results[i].getMP();
    uint32  pe = results[i].getPE();

    if     (results[i].getDelete()) {
      fprintf(outputFile, "frg iid "F_U32" isdeleted 1\n", i);
      fprintf(outputFile, "frg iid "F_U32" mateiid 0\n", i);
      numUM++;
      numDL++;

    } else if (results[i].getUnMate()) {
      fprintf(outputFile, "frg iid "F_U32" mateiid 0\n", i);
      numUM++;

    } else if ((mp > 0) && (pe > 0)) {
      //fprintf(outputFile, "frg iid "F_U32" mateiid 0\n", i);
      numAB++;

    } else if (mp > 0) {
      numMP++;

    } else if (pe > 0) {
      fprintf(outputFile, "frg iid "F_U32" mateiid 0\n", i);
      numPE++;

    } else if (no > 0) {
      numNO++;

    } else {
      numIG++;
    }
  }

  fprintf(stderr, "NO "F_U32" PE "F_U32" MP "F_U32" AB "F_U32" DL "F_U32" UM "F_U32" IG "F_U32"\n",
          numNO, numPE, numMP, numAB, numDL, numUM, numIG);

  delete [] results;
  delete    gkpStore;

  return(0);
}
