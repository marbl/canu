
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
 *  This file is derived from:
 *
 *    src/AS_CNS/addReadsToUnitigs.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2013-AUG-24 to 2014-APR-22
 *      are Copyright 2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-APR-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "splitToWords.H"

#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

#include <vector>
#include <string>
#include <map>

using namespace std;

class readMap {
public:
  readMap() {
    good = false;
    proc = false;

    rFWD = false;
    rIID = UINT32_MAX;
    rCNT = 0;

    tIID = UINT32_MAX;
    tBGN = 0;
    tEND = 0;
  };

  bool    good;
  bool    proc;

  bool    rFWD;
  uint32  rIID;
  uint32  rCNT;

  uint32  tIID;
  uint32  tBGN;
  uint32  tEND;
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

#ifdef UNFINISHED_ADD_TO_SINGLETON
  bool   doPlaceUnmapped  = true;
  bool   doDeleteUnmapped = false;
#endif

  int32  numFailures = 0;
  int32  numSkipped  = 0;

  bool   showResult = false;

  char               *lookupFile = NULL;
  map<string,uint32>  lookupIID;

#ifdef UNFINISHED_ADD_TO_SINGLETON
  vector<bool>        iidInTig;       //  true if the read is already in a unitig
  vector<uint32>      iidInTigByLib;  //  count of reads in tigs, by library
  vector<uint32>      iidInLib;       //  count of reads, by library
#endif

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
      while ((argv[arg+1] != NULL) && (AS_UTL_fileExists(argv[arg+1], false, false) == true))
        alignMapNames.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-lookup") == 0) {
      lookupFile = argv[++arg];

    } else if (strcmp(argv[arg], "-r") == 0) {
      doConsensus = true;

    } else if (strcmp(argv[arg], "-v") == 0) {
      showResult = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      VERBOSE_MULTIALIGN_OUTPUT++;

    } else if (strcmp(argv[arg], "-loadall") == 0) {
      loadall = true;

    } else if (strcmp(argv[arg], "-n") == 0) {
      doModify = false;

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
  if (lookupFile == NULL)
    err++;
  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version -m coords\n", argv[0]);
    fprintf(stderr, "  -g gkpStore           gatekeeper store\n");
    fprintf(stderr, "  -t tigStore version   tigStore and version to modify\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m map-file           input map coords\n");
    fprintf(stderr, "  -M fastqUIDmap        gatekeeper output fastqUIDmap for read name to IID translation\n");
    fprintf(stderr, "\n");
#if 0
    fprintf(stderr, "unmapped reads:  default is to promote to singleton unitigs\n");
    fprintf(stderr, "  -U                    leave unmapped reads alone (will crash CGW)\n");
    fprintf(stderr, "  -D                    delete unmapped reads from gkpStore\n");
#else
    fprintf(stderr, "unmapped reads: all reads that are mapped and eligible for addition must be\n");
    fprintf(stderr, "marked as deleted before running this program.  reads that are added will be\n");
    fprintf(stderr, "undeleted.  reads that are not added will remain deleted.\n");
#endif
    fprintf(stderr, "\n");
    fprintf(stderr, "consensus:  default is to not rebuild consensus\n");
    fprintf(stderr, "  -r                    rebuild consensus including the new reads\n");
    fprintf(stderr, "  -v                      show result\n");
    fprintf(stderr, "  -V                      verbose\n");
    fprintf(stderr, "  -loadall                load all reads in gkpStore into memory (faster consensus)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -n                    do all the work, but discard the result\n");

    if (gkpName == NULL)
      fprintf(stderr, "ERROR: no gkpStore (-g) supplied.\n");
    if (tigName == NULL)
      fprintf(stderr, "ERROR: no tigStore (-t) supplied.\n");
    if (alignMapNames.size() == 0)
      fprintf(stderr, "ERROR: no map-file (-m) inputs supplied.\n");
    if (lookupFile == NULL)
      fprintf(stderr, "ERROR: no fasqUIDmap (-M) supplied.\n");

    exit(1);
  }


  {
    fprintf(stderr, "Loading Name to IID map from '%s'\n", lookupFile);

    errno = 0;
    FILE *LF = fopen(lookupFile, "r");
    if (errno)
      fprintf(stderr, "Failed to open fastqUIDmap '%s'\n", lookupFile);

    char  LL[1024];
    fgets(LL, 1024, LF);

    while (!feof(LF)) {
      chomp(LL);  //  Shouldn't be necessary, but splitToWords isn't doing it!

      splitToWords  SW(LL);

      if        (SW.numWords() == 3) {
        lookupIID[string(SW[2])] = SW(1);

      } else if (SW.numWords() == 6) {
        lookupIID[string(SW[2])] = SW(1);
        lookupIID[string(SW[5])] = SW(4);
        //fprintf(stderr, "'%s' - %u -- '%s' - %u\n",
        //        SW[2], SW(1), SW[5], SW(4));

      } else {
      }

      fgets(LL, 1024, LF);
    }

    fprintf(stderr, "Loaded "F_SIZE_T" name to IIDs\n", lookupIID.size());
  }


#ifdef UNFINISHED_ADD_TO_SINGLETON
  {
    fprintf(stderr, "Loading tig placement from tigStore '%s'\n", tigName);

    tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, false, false, false);  //  Read only

    for (uint32 ti=0; ti<tigStore->numUnitigs(); ti++) {
      if (tigStore->isDeleted(ti, true))
        continue;

      MultiAlignT *ma = tigStore->loadMultiAlign(ti, true);

      if (ma == NULL)
        continue;

      uint32   fiMax = GetNumIntMultiPoss(ma->f_list);

      for (uint32 fi=0; fi<fiMax; fi++) {
        IntMultiPos  *imp = GetIntMultiPos(ma->f_list, fi);

        iidInTig[imp->ident] = true;
      }

      tigStore->unloadMultiAlign(ti, true);
    }

    delete tigStore;
  }
#endif


  //
  //  Load the alignment map.  The challenge here is to parse the unitig and read names
  //  into correct IIDs.  We assume that:
  //    Reads were dumped with -dumpfasta and have names ">UID,IID"
  //    Unitigs were dumped with tigStore -d consensus and have names "utgIID"
  //    Alignments are in the convertToExtent -extended format
  //


  uint32  totAligns = 0;
  uint32  totUnmap  = 0;
  uint32  totUnique = 0;
  uint32  totDups   = 0;

  uint32  lowID     = UINT32_MAX;
  uint32  maxID     = 0;

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

      //  Some gatekeeper dumps contain NAME,IID.  Strip out the ,IID, BEFORE converting to words.
      bool  cF = false;
      for (uint32 xx=0; !isspace(L[xx]); xx++) {
        if (L[xx] == ',')
          cF = true;
        if (cF == true)
          L[xx] = ' ';
      }

      splitToWords S(L);
      readMap      rm;

      rm.rIID = S(1);
      rm.rFWD = (S(4) < S(5));
      rm.rCNT = 1;

      assert(S[6][0] == 'u');  //  Unitig must be from a tigStore dump, and will be
      assert(S[6][1] == 't');  //  named 'utg#######'.
      assert(S[6][2] == 'g');

      rm.tIID = atoi(S[6] + 3);
      rm.tBGN = S(8);
      rm.tEND = S(9);

      //  We no longer use the encoded IID, but always lookup the read name in the map.
      rm.rIID = lookupIID[string(S[0])];

      if (rm.rIID == 0)
        fprintf(stderr, "FAILED to find IID for read '%s'\n", S[0]);

      //fprintf(stderr, "%d %d %d -- %d %d %d\n", rm.rIID, rm.rFWD, rm.rCNT, rm.tIID, rm.tBGN, rm.tEND);

      if (RM.size() <= rm.rIID)
        RM.resize(2 * rm.rIID);

      rm.rCNT = RM[rm.rIID].rCNT + 1;

      RM[rm.rIID] = rm;

      lowID = (lowID < rm.rIID) ? lowID : rm.rIID;
      maxID = (maxID < rm.rIID) ? rm.rIID : maxID;

      fgets(L, 1024, M);
    }

    fclose(M);

    fprintf(stderr, "%s - %u aligns\n", alignMapNames[an], numAligns);

    totAligns += numAligns;
  }

  for (uint32 xx=lowID; xx<=maxID; xx++) {
    if (RM[xx].rCNT == 0)
      totUnmap++;
    if (RM[xx].rCNT == 1)
      totUnique++;
    if (RM[xx].rCNT > 1)
      totDups++;
  }

  fprintf(stderr, "Loaded %u aligns from ID %u to %u, %u reads unmapped, %u reads had a unique alignment and %u reads had multiple aligns.\n",
          totAligns, lowID, maxID, totUnmap, totUnique, totDups);

  //
  //  Update deletion status in gkpStore.  This processes every read, regardless.
  //

  fprintf(stderr, "Processing mate pairs, updating gkpStore.\n");
  gkpStore = gkStore::gkStore_open(gkpName, false, true);  //  last arg - TRUE - writable

  uint32   unpaired    = 0;
  uint32   multiple    = 0;
  uint32   pairsToSame = 0;
  uint32   pairsToDiff = 0;

#ifdef UNFINISHED_ADD_TO_SINGLETON
  //  Need to process all reads, since we don't know where the first/last unmapped read is!
  //  We could instead process from the first to last deleted read in gkpStore, or ask which
  //  libraries were being added.
  lowID = 1;
  maxID = gkpStore->gkStore_getNumFragments();
#endif

  for (uint32 ff=lowID, mm=lowID; ff<=maxID; ff++) {
    gkFragment  read;
    gkFragment  mate;

    //  I think this is just a short circuit of the two checks of 'one or both reads has too
    //  few/many mappings' below.
    //if (RM[ff].rCNT != 1)
    //  //  Not mapped, mapped too much
    //  continue;

    gkpStore->gkStore_getFragment(ff, &read, GKFRAGMENT_INF);

    mm = read.gkFragment_getMateIID();

    if (mm == 0)
      //  No mate, pacbio read?
      continue;

    if (mm < ff)
      //  Already processed.
      continue;

    gkpStore->gkStore_getFragment(mm, &mate, GKFRAGMENT_INF);

    if ((RM[ff].rCNT == 0) ||
        (RM[mm].rCNT == 0)) {
      //  One or both reads has too few mappings.
      unpaired++;
      RM[ff].good = false;
      RM[mm].good = false;
      continue;
    }

    if ((RM[ff].rCNT > 1) ||
        (RM[mm].rCNT > 1)) {
      //  One or both reads has too many mappings.
      multiple++;
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

  gkpStore->gkStore_close();

  fprintf(stderr, "Will NOT add %u pairs - one read failed to map.\n", unpaired);
  fprintf(stderr, "Will NOT add %u pairs - multiple mappings.\n", multiple);
  fprintf(stderr, "Will add %u pairs in the same unitig\n", pairsToSame);
  fprintf(stderr, "Will add %u pairs in different unitigs\n", pairsToDiff);

  //
  //  Open stores.  gkpStore cannot be opened for writing, because then we can't loadall.
  //

  gkpStore = gkStore::gkStore_open(gkpName, false, false);                              //  last arg - false - not writable
  tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, true, true, false);  //  Write back to the same version

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

    MultiAlignT *ma = tigStore->loadMultiAlign(RM[bb].tIID, true);

    uint32   ungapLength = GetMultiAlignUngappedLength(ma);
    uint32   gapLength   = GetMultiAlignLength(ma);

    vector<int32>  ungapToGap;

    GetMultiAlignUngapToGap(ma, ungapToGap);

    //fprintf(stderr, "Loaded UTG %u offset size %lu\n", ma->maID, ungapToGap.size());

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

      fprintf(stdout, "bb=%u ee=%u ADD frag %u to unitig %u at %u,%u (from ungapped %u,%u)\n",
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
      tigStore->insertMultiAlign(ma, true, false);
    }
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}

