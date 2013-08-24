
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2013, J. Craig Venter Institute.
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

const char *mainid = "$Id:  $";

#include "AS_global.H"

#include "AS_UTL_splitToWords.H"

#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

#include <vector>

using namespace std;

class readMap {
public:
  readMap() {
    rIID = UINT32_MAX;
    rFWD = false;
    rCNT = 0;

    tIID = UINT32_MAX;
    tBGN = 0;
    tEND = 0;

    good = false;
    proc = false;
  };

  uint32  rIID;
  bool    rFWD;
  uint32  rCNT;

  uint32  tIID;
  uint32  tBGN;
  uint32  tEND;

  bool    good;
  bool    proc;
};


class ungapToGap {
public:
  uint32 *gapToUngap;
};



int
main(int argc, char **argv) {
  char  *gkpName = NULL;

  char  *tigName = NULL;
  int32  tigVers = -1;

  vector<char *>  alignMapNames;

  bool   doConsensus = false;
  bool   doModify    = true;

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

  bool   ignoreContains = false;
  double ignoreContainT = 0.75;

  bool   inplace = false;
  bool   loadall = false;

  CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                          CNS_OPTIONS_MIN_ANCHOR_DEFAULT,
                          CNS_OPTIONS_DO_PHASING_DEFAULT };

  vector<readMap>     RM;
  vector<ungapToGap>  UG;



  //  Comminucate to MultiAlignment_CNS.c that we are doing consensus and not cgw.
  thisIsConsensus = 1;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

      if (tigVers <= 0)
        fprintf(stderr, "invalid tigStore version (-t store version) '-t %s %s'.\n", argv[arg-1], argv[arg]), exit(1);

    } else if (strcmp(argv[arg], "-m") == 0) {
      alignMapNames.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-r") == 0) {
      doConsensus = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else if (strcmp(argv[arg], "-n") == 0) {
      doModify = false;

    } else if (strcmp(argv[arg], "-nocontains") == 0) {
      ignoreContains = true;
      ignoreContainT = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-inplace") == 0) {
      inplace = true;

    } else if (strcmp(argv[arg], "-loadall") == 0) {
      loadall = true;

    } else {
      err++;
    }

    arg++;
  }
  if (gkpName == NULL)
    err++;
  if (tigName == NULL)
    err++;
  if (alignMapNames.size() == 0)
    err++;
  if (err) {
    exit(1);
  }

  //
  //  Load the alignment map.  The challenge here is to parse the unitig and read names
  //  into correct IIDs.  We assume that:
  //    Reads were dumped with -dumpfasta and have names ">UID,IID"
  //    Unitigs were dumped with tigStore -d consensus and have names "utgIID"
  //

  uint32  totAligns = 0;
  uint32  totUnique = 0;
  uint32  totDups   = 0;

  for (uint32 an=0; an<alignMapNames.size(); an++) {
    uint32  numAligns = 0;
    uint32  numUnique = 0;
    uint32  numDups   = 0;

    errno = 0;
    FILE *M = fopen(alignMapNames[an], "r");
    if (errno)
      fprintf(stderr, "failed to open '%s': %s\n", alignMapNames[an], strerror(errno)), exit(1);

    char  L[1024];

    fgets(L, 1024, M);  //  Header line
    fgets(L, 1024, M);  //  First data

    while (!feof(M)) {
      numAligns++;

      for (uint32 xx=0; !isspace(L[xx]); xx++)
        if (L[xx] == ',')
          L[xx] = ' ';

      splitToWords S(L);
      readMap      rm;

      rm.rIID = S(1);
      rm.rFWD = (S(5) < S(6));
      rm.rCNT = 1;

      rm.tIID = atoi(S[7] + 3);
      rm.tBGN = S(9);
      rm.tEND = S(10);

      if (RM.size() < rm.rIID)
        RM.resize(2 * rm.rIID);

      rm.rCNT = RM[rm.rIID].rCNT + 1;

      RM[rm.rIID] = rm;

      fgets(L, 1024, M);
    }

    fclose(M);

    fprintf(stderr, "%s - %u aligns\n", alignMapNames[an], numAligns);

    totAligns += numAligns;
  }

  for (uint32 xx=0; xx<RM.size(); xx++) {
    if (RM[xx].rCNT == 1)
      totUnique++;
    if (RM[xx].rCNT > 1)
      totDups++;
  }

  fprintf(stderr, "Loaded %u aligns, %u reads had a unique alignment and %u reads had multiple aligns.\n",
          totAligns, totUnique, totDups);

  //
  //  Update deletion status in gkpStore.  This processes every read, regardless.
  //

  fprintf(stderr, "Processing mate pairs, updating gkpStore.\n");

  uint32   pairsToSame = 0;
  uint32   pairsToDiff = 0;

  gkpStore = new gkStore(gkpName, FALSE, TRUE);  //  last arg - TRUE - writable

  for (uint32 ff=0, mm=0; ff<RM.size(); ff++) {
    gkFragment  read;
    gkFragment  mate;

    if (RM[ff].rCNT != 1)
      //  Not mapped, mapped too much
      continue;

    gkpStore->gkStore_getFragment(ff, &read, GKFRAGMENT_INF);

    mm = read.gkFragment_getMateIID();

    if (mm == 0)
      //  No mate, pacbio read?
      continue;

    gkpStore->gkStore_getFragment(mm, &mate, GKFRAGMENT_INF);

    if ((RM[ff].rCNT != 1) ||
        (RM[mm].rCNT != 1)) {
      //  One or both reads has too few or too many mappings.  Don't use the pair.
      RM[ff].good = false;
      RM[mm].good = false;
      continue;
    }

    RM[ff].good = true;
    RM[mm].good = true;

    read.gkFragment_setIsDeleted(false);
    mate.gkFragment_setIsDeleted(false);

    gkpStore->gkStore_setFragment(&read);
    gkpStore->gkStore_setFragment(&mate);

    if (RM[ff].tIID == RM[mm].tIID)
      pairsToSame++;
    else
      pairsToDiff++;
  }

  delete gkpStore;

  fprintf(stderr, "Will add %u pairs in the same unitig\n", pairsToSame);
  fprintf(stderr, "Will add %u pairs in different unitigs\n", pairsToDiff);

  //
  //  Open stores.  gkpStore cannot be opened fr writing, because then we can't loadall.
  //

  gkpStore = new gkStore(gkpName, FALSE, FALSE);  //  last arg - FALSE - not writable
  tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);  //  Write back to the same version

  if (loadall) {
    fprintf(stderr, "Loading all reads into memory.\n");
    gkpStore->gkStore_load(0, 0, GKFRAGMENT_QLT);  //  fails if gkStore is writable
  }

  //
  //  Rebuild unitigs, stuff them back into the same version.
  //

  //  Argh, really should convert this to a vector right now....
  VA_TYPE(int32) *unGappedOffsets = CreateVA_int32(1024 * 1024);

  for (uint32 bb=0; bb<RM.size(); bb++) {
    if ((RM[bb].tIID == UINT32_MAX) ||   //  not mapped into a unitig
        (RM[bb].good == false)      ||   //  not useful mate pair
        (RM[bb].proc == true))           //  already added
      continue;

    MultiAlignT *ma = tigStore->loadMultiAlign(RM[bb].tIID, TRUE);

    uint32   ungapLength = GetMultiAlignUngappedLength(ma);
    uint32   gapLength   = GetMultiAlignLength(ma);

    vector<int32>  ungapToGap;

    GetMultiAlignUngapToGap(ma, ungapToGap);

    fprintf(stderr, "Loaded UTG %u offset size %lu\n", ma->maID, ungapToGap.size());

    //for (uint32 xx=0; xx<ungapToGap.size(); xx++)
    //  fprintf(stderr, "ungap %u -> gap %u\n", xx, ungapToGap[xx]);

    uint32  readsAdded = 0;

    for (uint32 ee=bb; ee<RM.size(); ee++) {
      if ((RM[ee].tIID != ma->maID) ||
          (RM[ee].good == false))
        continue;

      uint32  bgn = ungapToGap[RM[ee].tBGN];
      uint32  end = ungapToGap[RM[ee].tEND];

      readsAdded++;

      fprintf(stderr, "bb=%u ee=%u ADD frag %u to unitig %u at %u,%u (from %u,%u)\n",
              bb, ee,
              RM[ee].rIID, RM[ee].tIID, bgn, end, RM[ee].tBGN, RM[ee].tEND);

      //  Add a read to the unitig.

      IntMultiPos  frg;

      frg.type         = AS_READ;
      frg.ident        = RM[ee].rIID;
      frg.contained    = 0;
      frg.parent       = 0;
      frg.ahang        = 0;
      frg.bhang        = 0;
      frg.position.bgn = (RM[ee].rFWD) ? bgn : end,
      frg.position.end = (RM[ee].rFWD) ? end : bgn,
      frg.delta_length = 0;
      frg.delta        = NULL;

      AppendVA_IntMultiPos(ma->f_list, &frg);

      //  Mark that we've processed this read.
      assert(RM[ee].proc == false);
      RM[ee].proc = true;
    }

    fprintf(stderr, "Added %u reads to unitig %u (previously %lu reads)\n",
            readsAdded,
            ma->maID,
            GetNumIntMultiPoss(ma->f_list) - readsAdded);

    if (doConsensus) {
      fprintf(stderr, "Regenerating consensus.\n");

      if (MultiAlignUnitig(ma, gkpStore, &options, NULL)) {
        if (showResult)
          PrintMultiAlignT(stdout, ma, gkpStore, false, false, AS_READ_CLEAR_LATEST);
      } else {
        fprintf(stderr, "MultiAlignUnitig()-- unitig %d failed.\n", ma->maID);
        numFailures++;
      }
    }

    if (doModify) {
      fprintf(stderr, "Updating unitig %u\n", ma->maID);
      tigStore->insertMultiAlign(ma, TRUE, FALSE);
    }
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}

