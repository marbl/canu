
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

const char *mainid = "$Id: classifyMatesApply.C,v 1.4 2012-05-18 21:28:22 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_OVS_overlapStore.h"
#include "AS_PER_gkpStore.h"

#include "classifyMates.H"

//  One per fragment, stores the results from multiple classification runs.
class classifyMatesSummary {
public:
  classifyMatesSummary() {
    isTested   = 0;
    doDelete   = 0;
    doUnMate   = 0;

    isSpur     = 0;
    isChimer   = 0;
    isJunction = 0;

    isClassified = 0;
    isLimited    = 0;
    isExhausted  = 0;

    numNO      = 0;
    numPE      = 0;
    numMP      = 0;
  };

  void      setTested(void) { isTested = 1; };
  bool      getTested(void) { return(isTested > 0); };

  void      setDelete(void)  { doDelete = 1; };
  bool      getDelete(void)  { return(doDelete == 1); };

  void      setUnMate(void)  { doUnMate = 1; };
  bool      getUnMate(void)  { return(doUnMate == 1); };


  void      setIsSpur(void)     { isSpur = true; };
  void      setIsChimer(void)   { isChimer = true; };
  void      setIsJunction(void) { isJunction = true; };

  bool      getIsSpur(void)     { return(isSpur); };
  bool      getIsChimer(void)   { return(isChimer); };
  bool      getIsJunction(void) { return(isJunction); };


  void      setIsClassified(void)   { isClassified = true; };
  void      setIsLimited(void)      { isLimited = true; };
  void      setIsExhausted(void)    { isExhausted = true; };

  bool      getIsClassified(void)   { return(isClassified); };
  bool      getIsLimited(void)      { return(isLimited); };
  bool      getIsExhausted(void)    { return(isExhausted); };


  void      addNO(void) { if (numNO < 255)  numNO++; };
  void      addPE(void) { if (numPE < 255)  numPE++; };
  void      addMP(void) { if (numMP < 255)  numMP++; };

  uint8     getNO(void) { return(numNO); };
  uint8     getPE(void) { return(numPE); };
  uint8     getMP(void) { return(numMP); };


private:
  uint32    isTested     : 1;  //  We have attempted classification on this read
  uint32    doDelete     : 1;  //  Delete the read, it's junk
  uint32    doUnMate     : 1;  //  Remove the mate edge, it's mate read is junk

  //  Any classification called this read...
  uint32    isSpur       : 1;
  uint32    isChimer     : 1;
  uint32    isJunction   : 1;

  uint32    isClassified : 1;
  uint32    isLimited    : 1;
  uint32    isExhausted  : 1;

  //  Number of votes for....
  uint32    numNO        : 8;  //  ...no classification
  uint32    numPE        : 8;  //  ...PE pair
  uint32    numMP        : 8;  //  ...MP pair
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

  uint32  numAspur     = 0;  //  These stats are only correct if each pair is processed once
  uint32  numAchimer   = 0;
  uint32  numAjunction = 0;

  uint32  numBspur     = 0;
  uint32  numBchimer   = 0;
  uint32  numBjunction = 0;

  uint32  numABspurChimer[4][4];

  memset(numABspurChimer, 0, sizeof(uint32) * 16);

  for (uint32 i=0; i<resultsNamesLen; i++) {
    classifyMatesResultFile  *RF = new classifyMatesResultFile(resultsNames[i]);

    while (RF->read(result)) {
      AS_IID   f = result.readIID;
      AS_IID   m = result.mateIID;

      if (justPrint)
        result.print(stdout);

      results[f].setTested();
      results[m].setTested();

      if (result.fragSpur     == true)  { numAspur++;       results[f].setIsSpur(); }
      if (result.fragChimer   == true)  { numAchimer++;     results[f].setIsChimer(); }
      if (result.fragJunction == true)  { numAjunction++;   results[f].setIsJunction(); }

      if (result.mateSpur     == true)  { numBspur++;       results[m].setIsSpur(); }
      if (result.mateChimer   == true)  { numBchimer++;     results[m].setIsChimer(); }
      if (result.mateJunction == true)  { numBjunction++;   results[m].setIsJunction(); }

      assert(result.fragJunction == false);  //  unused
      assert(result.mateJunction == false);  //  unused

      assert(result.fragSpur + result.fragChimer != 2);
      assert(result.mateSpur + result.mateChimer != 2);

      uint32 a = result.fragSpur + 2 * result.fragChimer;
      uint32 b = result.mateSpur + 2 * result.mateChimer;

      numABspurChimer[a][b]++;

      if (result.classified   == true) {
        results[f].setIsClassified();
        results[m].setIsClassified();
      }

      if (result.limited      == true) {
        results[f].setIsLimited();
        results[m].setIsLimited();
      }

      if (result.exhausted    == true) {
        results[f].setIsExhausted();
        results[m].setIsExhausted();
      }

      //  If any single result says the read (or mate) is junk, delete it.

      if ((result.fragSpur     == true) ||
          (result.fragChimer   == true) ||
          (result.fragJunction == true)) {
        results[f].setDelete();  //  Delete the frag
        results[f].setUnMate();  //  and unmate
        results[m].setUnMate();
      }

      if ((result.mateSpur     == true) ||
          (result.mateChimer   == true) ||
          (result.mateJunction == true)) {
        results[m].setDelete();  //  Delete the mate
        results[f].setUnMate();  //  and unmate
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

  uint32  numSpur     = 0;
  uint32  numChimer   = 0;
  uint32  numJunction = 0;

  uint32  numClassified = 0;
  uint32  numLimited    = 0;
  uint32  numExhausted  = 0;

  //  Process the classifications.  If any are labeled both PE and MP, remove the mate.

  uint32  numDL = 0;
  uint32  numUM = 0;
  uint32  numNO = 0;
  uint32  numPE = 0;
  uint32  numMP = 0;
  uint32  numAB = 0;
  uint32  numIG = 0;

  for (uint32 i=1; i<=numFrags; i++) {
    uint32  no = results[i].getNO();
    uint32  mp = results[i].getMP();
    uint32  pe = results[i].getPE();

    if     (results[i].getDelete()) {
      fprintf(outputFile, "frg iid "F_U32" isdeleted 1\n", i);
      fprintf(outputFile, "frg iid "F_U32" mateiid 0\n", i);
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

    if (results[i].getIsSpur())        numSpur++;
    if (results[i].getIsChimer())      numChimer++;
    if (results[i].getIsJunction())    numJunction++;

    if (results[i].getIsClassified())  numClassified++;
    if (results[i].getIsLimited())     numLimited++;
    if (results[i].getIsExhausted())   numExhausted++;
  }

  fprintf(stderr, "Read Fate:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   unclassified (remains mated)\n", numNO);
  fprintf(stderr, "%10"F_U32P"   paired-end (break the mate)\n", numPE);
  fprintf(stderr, "%10"F_U32P"   mate-pair (remains mated)\n", numMP);
  fprintf(stderr, "%10"F_U32P"   paired-end + mate-pair (remains mated)\n", numAB);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   delete (read is spur/chimer/junction)\n", numDL);
  fprintf(stderr, "%10"F_U32P"   unmate (pair is spur/chimer/junction)\n", numUM);
  fprintf(stderr, "%10"F_U32P"   ignored (not searched)\n", numIG);
  fprintf(stderr, "\n");

  //fprintf(stderr, "NO "F_U32" PE "F_U32" MP "F_U32" AB "F_U32" DL "F_U32" UM "F_U32" IG "F_U32"\n",
  //        numNO, numPE, numMP, numAB, numDL, numUM, numIG);

  fprintf(stderr, "Read Classification:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   A spur reads\n",     numAspur);
  fprintf(stderr, "%10"F_U32P"   A chimer reads\n",   numAchimer);
  fprintf(stderr, "%10"F_U32P"   A junction reads\n", numAjunction);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   B spur reads\n",     numBspur);
  fprintf(stderr, "%10"F_U32P"   B chimer reads\n",   numBchimer);
  fprintf(stderr, "%10"F_U32P"   B junction reads\n", numBjunction);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   spur reads\n",       numSpur);
  fprintf(stderr, "%10"F_U32P"   chimer reads\n",     numChimer);
  fprintf(stderr, "%10"F_U32P"   junction reads\n",   numJunction);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   good/good\n",     numABspurChimer[0][0]);
  fprintf(stderr, "%10"F_U32P"   good/spur\n",     numABspurChimer[0][1]);
  fprintf(stderr, "%10"F_U32P"   good/chimer\n",   numABspurChimer[0][2]);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   spur/good\n",     numABspurChimer[1][0]);
  fprintf(stderr, "%10"F_U32P"   spur/spur\n",     numABspurChimer[1][1]);
  fprintf(stderr, "%10"F_U32P"   spur/chimer\n",   numABspurChimer[1][2]);
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   chimer/good\n",   numABspurChimer[2][0]);
  fprintf(stderr, "%10"F_U32P"   chimer/spur\n",   numABspurChimer[2][1]);
  fprintf(stderr, "%10"F_U32P"   chimer/chimer\n", numABspurChimer[2][2]);
  fprintf(stderr, "\n");
  fprintf(stderr, "Read Termination:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "%10"F_U32P"   classified reads\n", numClassified);
  fprintf(stderr, "%10"F_U32P"   limited reads\n",    numLimited);
  fprintf(stderr, "%10"F_U32P"   exhausted reads\n",  numExhausted);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");

  delete [] results;
  delete    gkpStore;

  return(0);
}
