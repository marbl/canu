
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

const char *mainid = "$Id$";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"

#include  <map>
#include  <vector>
#include  <algorithm>

using namespace std;

#define F10U         "%10"F_U64P
#define F00U         "%"F_U64P
#define UNDEFINED  (uint64)0


map<AS_UID,AS_IID>  uid2iid;

vector<uint64>     allScaffolds_contigs;
vector<uint64>     allScaffolds_bases;
vector<uint64>     allScaffolds_span;
vector<uint64>     allScaffolds_gaps;

vector<uint64>     smallScaffolds_contigs;
vector<uint64>     smallScaffolds_bases;
vector<uint64>     smallScaffolds_span;
vector<uint64>     smallScaffolds_gaps;

vector<uint64>     largeScaffolds_contigs;
vector<uint64>     largeScaffolds_bases;
vector<uint64>     largeScaffolds_span;
vector<uint64>     largeScaffolds_gaps;

vector<uint64>     allContigs;
vector<uint64>     scaffContigs;
vector<uint64>     smallContigs;
vector<uint64>     largeContigs;
vector<uint64>     degenContigs;

vector<uint64>     allContigs_reads;
vector<uint64>     scaffContigs_reads;
vector<uint64>     smallContigs_reads;
vector<uint64>     largeContigs_reads;
vector<uint64>     degenContigs_reads;

vector<uint64>     all_reads;
vector<uint64>     surrogate_reads;
vector<uint64>     singleton_reads;

vector<uint64>     allUnitigs;
vector<uint64>     uniqueUnitigs;
vector<uint64>     notRezUnitigs;
vector<uint64>     surrogateUnitigs;
vector<uint64>     otherUnitigs;

vector<uint32>     contigLength;
vector<uint32>     contigGCBases;

//  frgLength (from the asm file) should be the same as frgClearLen[0] (from the store)
vector<uint32>     frgLength;

uint64             surrogateInstances;

uint64             varsInScaffolds;
uint64             varsInDegenerates;

uint64             readsWithNoMate;
uint64             readsWithGoodMate;
uint64             readsWithBadShortMate;
uint64             readsWithBadLongMate;
uint64             readsWithSameOrientMate;
uint64             readsWithOuttieMate;
uint64             readsWithBothChaffMate;
uint64             readsWithChaffMate;
uint64             readsWithBothDegenMate;
uint64             readsWithDegenMate;
uint64             readsWithBothSurrogateMate;
uint64             readsWithSurrogateMate;
uint64             readsWithDiffScafMate;
uint64             readsWithUnassignedMate;

uint64             totalReadsFromAFG;

vector<uint64>     scaffoldLinkWeight;
vector<uint64>     contigLinkWeight;

uint64             readsInContig;
uint64             readsInSmallContig;
uint64             readsInLargeContig;
uint64             readsInDegenerate;
uint64             readsInSurrogatePlaced;
uint64             readsInSurrogateUnplaced;
uint64             readsInChaff;





void
processMDI(SnapMateDistMesg *mdi) {
  uint32  samples = 0;

  for (int i=0; i<mdi->num_buckets; i++)
    samples += mdi->histogram[i];

  fprintf(stdout, "MDI\t%s\t%.6f\t%.6f\t%d\n",
          AS_UID_toString(mdi->erefines),
          mdi->mean, mdi->stddev, samples);
}


void
processAFG(AugFragMesg *afg) {
  uint32  iid = afg->iaccession;
  uint32  len = afg->clear_rng.end - afg->clear_rng.bgn;

  if (frgLength.capacity() <= iid)
    frgLength.reserve(iid + 10 * 1024 * 1024);

  frgLength[iid] = len;

  uid2iid[afg->eaccession] = afg->iaccession;

  totalReadsFromAFG++;

  if (afg->chaff) {
    readsInChaff++;
    singleton_reads.push_back(len);
  }

  all_reads.push_back(len);

  switch (afg->mate_status) {
    case 'Z':
      readsWithUnassignedMate++;
      break;
    case 'G':
      readsWithGoodMate++;
      break;
    case 'C':
      readsWithBadShortMate++;
      break;
    case 'L':
      readsWithBadLongMate++;
      break;
    case 'S':
      readsWithSameOrientMate++;
      break;
    case 'O':
      readsWithOuttieMate++;
      break;
    case 'N':
      readsWithNoMate++;
      break;
    case 'H':
      readsWithBothChaffMate++;
      break;
    case 'A':
      readsWithChaffMate++;
      break;
    case 'D':
      readsWithBothDegenMate++;
      break;
    case 'E':
      readsWithDegenMate++;
      break;
    case 'U':
      readsWithBothSurrogateMate++;
      break;
    case 'R':
      readsWithSurrogateMate++;
      break;
    case 'F':
      readsWithDiffScafMate++;
      break;
    default:
      readsWithUnassignedMate++;
      break;
  }
}


void
processAMP(AugMatePairMesg *amp) {
  char   mateStatus[256] = {0};

  //fprintf(mates, "%s\t%s\t%s\n",
  //        AS_UID_toString(amp->fragment1),
  //        AS_UID_toString(amp->fragment2),
  //        decodeMateStatus(amp->mate_status, mateStatus));
}


void
processUTG(SnapUnitigMesg *utg) {
  uint32    unitigLength       = 0;
  uint32    unitigLengthGapped = strlen(utg->consensus);

  for (uint32 i=0; i<unitigLengthGapped; i++) {
    if (utg->consensus[i] != '-')
      unitigLength++;
  }

  allUnitigs.push_back(unitigLength);

  switch (utg->status) {
    case AS_UNIQUE:
      uniqueUnitigs.push_back(unitigLength);
      break;
    case AS_NOTREZ:
      notRezUnitigs.push_back(unitigLength);
      break;
    case AS_SEP:
      surrogateUnitigs.push_back(unitigLength);
      break;
    case AS_UNASSIGNED:
      otherUnitigs.push_back(unitigLength);
      break;
    default:
      assert(0);
      break;
  }
}



void
processULK(SnapUnitigLinkMesg *ulk) {
}



void
processCCOfrags(SnapConConMesg *cco, vector<uint64> &reads) {

  for (int32 i=0; i<cco->num_pieces; i++) {
    uint32  iid = uid2iid[cco->pieces[i].eident];

    reads.push_back(frgLength[iid]);
  }
}



void
processCCO(SnapConConMesg *cco) {
  uint32    lenGapped    = strlen(cco->consensus);
  uint32    len          = 0;
  uint32    gc           = 0;
  uint32    isDegenerate = 0;
  uint32    iid          = cco->iaccession;

  for (uint32 i=0; i<lenGapped; i++) {
    if (cco->consensus[i] != '-')
      len++;
    if ((cco->consensus[i] == 'g') ||
        (cco->consensus[i] == 'G') ||
        (cco->consensus[i] == 'c') ||
        (cco->consensus[i] == 'C'))
      gc++;
  }

  allContigs.push_back(len);
  processCCOfrags(cco, allContigs_reads);

  if (contigLength.capacity() <= iid) {
    contigLength.reserve(iid + 1 * 1024 * 1024);
    contigGCBases.reserve(iid + 1 * 1024 * 1024);
  }

  contigLength[iid]  = len;
  contigGCBases[iid] = gc;

  uid2iid[cco->eaccession] = cco->iaccession;

  //  By definition, a degenerate contig has one unitig and is unplaced.
  //  In reality, those two conditions always occur together.

  if ((cco->placed == AS_UNPLACED) && (cco->num_unitigs == 1)) {
    degenContigs.push_back(len);
    processCCOfrags(cco, degenContigs_reads);
    varsInDegenerates += cco->num_vars;

  } else {
    scaffContigs.push_back(len);
    processCCOfrags(cco, scaffContigs_reads);

    if (len < 10000) {
      smallContigs.push_back(len);
      processCCOfrags(cco, smallContigs_reads);
    } else {
      largeContigs.push_back(len);
      processCCOfrags(cco, largeContigs_reads);
    }

    varsInScaffolds += cco->num_vars;
  }
}



void
processCLK(SnapContigLinkMesg *clk) {
}



#warning should be a library function
int32
computeGapSize(double gapsize) {
  if (gapsize <= 20.0)
    return(20);
  return((int32)gapsize);
}



void
processSCF(SnapScaffoldMesg *scf) {
  uint32           numContigs = 0;
  uint32           baseLength = 0;
  uint32           spanLength = 0;
  vector<uint32>   gapLengths;

  //  Tricky.  If scf->num_contig_pairs==0, there is a single contig
  //  in this scaffold.

  int  singleContig = 0;
  if (scf->num_contig_pairs == 0) {
    singleContig = 1;
    scf->num_contig_pairs = 1;
  }

  //  CTP/contig pairs
  for (int32 i=0; i<scf->num_contig_pairs; i++) {

    //  First contig-pair, print the first contig.
    //
    if (i == 0) {
      AS_IID ctgIID = uid2iid[scf->contig_pairs[i].econtig1];

      spanLength += contigLength[ctgIID];
      baseLength += contigLength[ctgIID];

      numContigs++;
    }

    //  Not first contig-pair, or there is more than one pair, print
    //  the seocnd contig.
    //
    if ((i > 0) || (singleContig == 0)) {
      AS_IID ctgIID = uid2iid[scf->contig_pairs[i].econtig2];
      uint32 gapLen = computeGapSize(scf->contig_pairs[i].mean);

      gapLengths.push_back(gapLen);

      spanLength += gapLen;
      spanLength += contigLength[ctgIID];
      baseLength += contigLength[ctgIID];

      numContigs++;
    }
  }

  allScaffolds_contigs.push_back(numContigs);
  allScaffolds_bases.push_back(baseLength);
  allScaffolds_span.push_back(spanLength);
  allScaffolds_gaps.insert(allScaffolds_gaps.end(), gapLengths.begin(), gapLengths.end());

  if (spanLength < 2000) {
    smallScaffolds_contigs.push_back(numContigs);
    smallScaffolds_bases.push_back(baseLength);
    smallScaffolds_span.push_back(spanLength);
    smallScaffolds_gaps.insert(smallScaffolds_gaps.end(), gapLengths.begin(), gapLengths.end());
  } else {
    largeScaffolds_contigs.push_back(numContigs);
    largeScaffolds_bases.push_back(baseLength);
    largeScaffolds_span.push_back(spanLength);
    largeScaffolds_gaps.insert(largeScaffolds_gaps.end(), gapLengths.begin(), gapLengths.end());
  }
}



void
processSLK(SnapScaffoldLinkMesg *slk) {
}




template <class T>
uint64
getNum(vector<T> v) {
  return((uint64)v.size());
}

template <class T>
uint64
getSum(vector<T> v) {
  if (v.size() == 0)
    return(0);
  uint64 r = 0;
  for (uint64 i=0; i<v.size(); i++)
    r += v[i];
  return(r);
}

template <class T>
uint64
getAve(vector<T> v) {
  if (v.size() == 0)
    return(0);
  uint64 r = 0;
  for (uint64 i=0; i<v.size(); i++)
    r += v[i];
  return(r / v.size());
}

template <class T>
uint64
getMin(vector<T> v) {
  if (v.size() == 0)
    return(0);
  uint64 r = v[0];
  for (uint64 i=0; i<v.size(); i++)
    r = MIN(r, v[i]);
  return(r);
}

template <class T>
uint64
getMax(vector<T> v) {
  if (v.size() == 0)
    return(0);
  uint64 r = v[0];
  for (uint64 i=0; i<v.size(); i++)
    r = MAX(r, v[i]);
  return(r);
}

template <class T>
uint64
getN(uint32 n, vector<T> v) {
  if (v.size() == 0)
    return(0);
  sort(v.rbegin(), v.rend());
  uint64 m = getSum(v) * n / 100;
  uint64 s = 0;
  uint64 i = 0;
  while (s < m) {
    s += v[i];
    i++;
  }
  return(v[i-1]);
}

template <class T>
double
getNfrac(uint32 n, vector<T> v) {
  if (v.size() == 0)
    return(0);
  sort(v.rbegin(), v.rend());
  uint64 M = getSum(v);
  uint64 m = M * n / 100;
  uint64 s = 0;
  uint64 i = 0;
  while (s < m) {
    s += v[i];
    i++;
  }
  return(100.0 * s / M);
}

template <class T>
uint64
getNidx(uint32 n, vector<T> v) {
  if (v.size() == 0)
    return(0);
  sort(v.rbegin(), v.rend());
  uint64 m = getSum(v) * n / 100;
  uint64 s = 0;
  uint64 i = 0;
  while (s < m) {
    s += v[i];
    i++;
  }
  return(i);
}






void
generateStatisticsFromASMFile(char *progName, char *asmName) {
  GenericMesg *pmesg       = NULL;

  errno = 0;
  FILE *asmFile = fopen(asmName, "r");
  if (errno)
    fprintf(stderr, "%s: failed to open '%s' for reading: %s\n",
            progName, asmName, strerror(errno));

  while(ReadProtoMesg_AS(asmFile, &pmesg) != EOF){
    switch(pmesg->t){
      case MESG_MDI:
        processMDI((SnapMateDistMesg *)pmesg->m);
        break;

      case MESG_AFG:
        processAFG((AugFragMesg *)pmesg->m);
        break;
      case MESG_AMP:
        processAMP((AugMatePairMesg *)pmesg->m);
        break;

      case MESG_UTG:
        processUTG((SnapUnitigMesg *)pmesg->m);
        break;
      case MESG_ULK:
        processULK((SnapUnitigLinkMesg *)pmesg->m);
        break;

      case MESG_CCO:
        processCCO((SnapConConMesg *)pmesg->m);
        break;
      case MESG_CLK:
        processCLK((SnapContigLinkMesg *)pmesg->m);
        break;

      case MESG_SCF:
        processSCF((SnapScaffoldMesg *)pmesg->m);
        break;
      case MESG_SLK:
        processSLK((SnapScaffoldLinkMesg *)pmesg->m);
        break;

      default:
        break;
    }
  }
}



void
generateStatisticsFromGKPStore(char *progName, char *gkpName) {
  //  Open the gkStore, but discard UID info for now.
  gkStore    *gkp = new gkStore(gkpName, FALSE, TRUE);
  gkStream   *gks = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);
  gkFragment  gkf;

  //for (uint32 c=0; c<AS_READ_CLEAR_NUM; c++)
  //  frgClearLength[c].reserve(gkp->gkStore_getNumFragments());

  class frgInfoPerLibrary {
  public:
    frgInfoPerLibrary() {
      numActiveFrag  = 0;
      numDeletedFrag = 0;
      numMatedFrag   = 0;
    };
    ~frgInfoPerLibrary() {
    };

    uint32             numActiveFrag;
    uint32             numDeletedFrag;
    uint32             numMatedFrag;

    vector<uint32>     readLen;

    vector<uint32>     clearBgn[AS_READ_CLEAR_NUM];
    vector<uint32>     clearEnd[AS_READ_CLEAR_NUM];
    vector<uint32>     clearLen[AS_READ_CLEAR_NUM];
  };

  frgInfoPerLibrary    global;
  frgInfoPerLibrary   *lib    = new frgInfoPerLibrary [gkp->gkStore_getNumLibraries() + 1];

  //  Scan the fragments.  Remember:
  //    untrimmed read length
  //    clear lengths for all clear types
  //    number deleted
  //    number alive
  //    number mated
  //  Complicated by also remembering all this per library.
  //
  while (gks->next(&gkf)) {
    AS_IID  l = gkf.gkFragment_getLibraryIID();

    if (gkf.gkFragment_getIsDeleted()) {
      global.numDeletedFrag++;
      lib[l].numDeletedFrag++;
      continue;
    }

    global.numActiveFrag++;
    lib[l].numActiveFrag++;

    if (gkf.gkFragment_getMateIID() > 0) {
      global.numMatedFrag++;
      lib[l].numMatedFrag++;
    }

    global.readLen.push_back(gkf.gkFragment_getSequenceLength());
    lib[l].readLen.push_back(gkf.gkFragment_getSequenceLength());

    for (uint32 c=0; c<AS_READ_CLEAR_NUM; c++) {
      uint32  b = gkf.gkFragment_getClearRegionBegin(c);
      uint32  e = gkf.gkFragment_getClearRegionEnd(c);

      if (b > e)
        continue;

      global.clearBgn[c].push_back(b);
      lib[l].clearBgn[c].push_back(b);

      global.clearEnd[c].push_back(e);
      lib[l].clearEnd[c].push_back(e);

      global.clearLen[c].push_back(e - b);
      lib[l].clearLen[c].push_back(e - b);
    }
  }

  uint64   basesInReads = getSum(global.readLen);

  fprintf(stdout, "<< all reads, untrimmed >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num reads                  = "F10U"\n", getNum(global.readLen));
  fprintf(stdout, "bases in reads             = "F10U" %7.3f%%\n", getSum(global.readLen), 100.0 * getSum(global.readLen) / basesInReads);
  fprintf(stdout, "bases in reads ave         = "F10U"\n", getAve(global.readLen));
  fprintf(stdout, "bases in reads min         = "F10U"\n", getMin(global.readLen));
  fprintf(stdout, "bases in reads max         = "F10U"\n", getMax(global.readLen));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in reads N25         = "F10U", N=%06.2f I="F00U"\n", getN(25, global.readLen), getNfrac(25, global.readLen), getNidx(25, global.readLen));
  fprintf(stdout, "bases in reads N50         = "F10U", N=%06.2f I="F00U"\n", getN(50, global.readLen), getNfrac(50, global.readLen), getNidx(50, global.readLen));
  fprintf(stdout, "bases in reads N75         = "F10U", N=%06.2f I="F00U"\n", getN(75, global.readLen), getNfrac(75, global.readLen), getNidx(75, global.readLen));
  fprintf(stdout, "\n");

  for (uint32 c=1; c<AS_READ_CLEAR_NUM; c++) {
    if (getNum(global.clearLen[c]) == 0)
      continue;
    fprintf(stdout, "<< all reads, trimmed %s >>\n", AS_READ_CLEAR_NAMES[c]);
    fprintf(stdout, "\n");
    fprintf(stdout, "num reads                  = "F10U"\n", getNum(global.clearLen[c]));
    fprintf(stdout, "bases in reads             = "F10U" %7.3f%%\n", getSum(global.clearLen[c]), 100.0 * getSum(global.clearLen[c]) / basesInReads);
    fprintf(stdout, "bases in reads ave         = "F10U"\n", getAve(global.clearLen[c]));
    fprintf(stdout, "bases in reads min         = "F10U"\n", getMin(global.clearLen[c]));
    fprintf(stdout, "bases in reads max         = "F10U"\n", getMax(global.clearLen[c]));
    fprintf(stdout, "\n");
    fprintf(stdout, "bases in reads N25         = "F10U", N=%06.2f I="F00U"\n", getN(25, global.clearLen[c]), getNfrac(25, global.clearLen[c]), getNidx(25, global.clearLen[c]));
    fprintf(stdout, "bases in reads N50         = "F10U", N=%06.2f I="F00U"\n", getN(50, global.clearLen[c]), getNfrac(50, global.clearLen[c]), getNidx(50, global.clearLen[c]));
    fprintf(stdout, "bases in reads N75         = "F10U", N=%06.2f I="F00U"\n", getN(75, global.clearLen[c]), getNfrac(75, global.clearLen[c]), getNidx(75, global.clearLen[c]));
    fprintf(stdout, "\n");
  }



  for (uint32 l=0; l<gkp->gkStore_getNumLibraries() + 1; l++) {
    fprintf(stdout, "<< reads in library %d, untrimmed >>\n", l);
    fprintf(stdout, "\n");
    fprintf(stdout, "num reads                  = "F10U"\n", getNum(lib[l].readLen));
    fprintf(stdout, "bases in reads             = "F10U" %7.3f%%\n", getSum(lib[l].readLen), 100.0 * getSum(lib[l].readLen) / basesInReads);
    fprintf(stdout, "bases in reads ave         = "F10U"\n", getAve(lib[l].readLen));
    fprintf(stdout, "bases in reads min         = "F10U"\n", getMin(lib[l].readLen));
    fprintf(stdout, "bases in reads max         = "F10U"\n", getMax(lib[l].readLen));
    fprintf(stdout, "\n");
    fprintf(stdout, "bases in reads N25         = "F10U", N=%06.2f I="F00U"\n", getN(25, lib[l].readLen), getNfrac(25, lib[l].readLen), getNidx(25, lib[l].readLen));
    fprintf(stdout, "bases in reads N50         = "F10U", N=%06.2f I="F00U"\n", getN(50, lib[l].readLen), getNfrac(50, lib[l].readLen), getNidx(50, lib[l].readLen));
    fprintf(stdout, "bases in reads N75         = "F10U", N=%06.2f I="F00U"\n", getN(75, lib[l].readLen), getNfrac(75, lib[l].readLen), getNidx(75, lib[l].readLen));
    fprintf(stdout, "\n");

    for (uint32 c=1; c<AS_READ_CLEAR_NUM; c++) {
      if (getNum(lib[l].clearLen[c]) == 0)
        continue;
      fprintf(stdout, "<< reads in library %d, trimmed %s >>\n", l, AS_READ_CLEAR_NAMES[c]);
      fprintf(stdout, "\n");
      fprintf(stdout, "num reads                  = "F10U"\n", getNum(lib[l].clearLen[c]));
      fprintf(stdout, "bases in reads             = "F10U" %7.3f%%\n", getSum(lib[l].clearLen[c]), 100.0 * getSum(lib[l].clearLen[c]) / basesInReads);
      fprintf(stdout, "bases in reads ave         = "F10U"\n", getAve(lib[l].clearLen[c]));
      fprintf(stdout, "bases in reads min         = "F10U"\n", getMin(lib[l].clearLen[c]));
      fprintf(stdout, "bases in reads max         = "F10U"\n", getMax(lib[l].clearLen[c]));
      fprintf(stdout, "\n");
      fprintf(stdout, "bases in reads N25         = "F10U", N=%06.2f I="F00U"\n", getN(25, lib[l].clearLen[c]), getNfrac(25, lib[l].clearLen[c]), getNidx(25, lib[l].clearLen[c]));
      fprintf(stdout, "bases in reads N50         = "F10U", N=%06.2f I="F00U"\n", getN(50, lib[l].clearLen[c]), getNfrac(50, lib[l].clearLen[c]), getNidx(50, lib[l].clearLen[c]));
      fprintf(stdout, "bases in reads N75         = "F10U", N=%06.2f I="F00U"\n", getN(75, lib[l].clearLen[c]), getNfrac(75, lib[l].clearLen[c]), getNidx(75, lib[l].clearLen[c]));
      fprintf(stdout, "\n");
    }
  }

  delete [] lib;
}



void
generateStatisticsFromStores(char *progName,
                             char *gkpName,
                             char *ovlName,
                             char *tigName, uint32 tigVersion,
                             char *ckpName, uint32 ckpVersion) {

}




int main (int argc, char *argv[]) {
  char   *outputPrefix       = NULL;
  char   *gkpName            = NULL;
  char   *ovlName            = NULL;
  char   *tigName            = NULL;
  uint32  tigVers            = 0;
  char   *ckpName            = NULL;
  uint32  ckpVers            = 0;
  char   *asmName            = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];


    } else if (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-O") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-C") == 0) {
      ckpName = argv[++arg];
      ckpVers = atoi(argv[++arg]);


    } else if (strcmp(argv[arg], "-a") == 0) {
      asmName = argv[++arg];


    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((asmName == NULL) || (outputPrefix == NULL) || (err)) {
    fprintf(stderr, "usage: %s -o prefix -a asmFile [-h]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Assembly Statistics from a .asm file:\n");
    fprintf(stderr, "  -a asmFile       read the assembly from here\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Assembly Statistics from binary stores:\n");
    fprintf(stderr, "  -G gkpStore      gatekeeper store\n");
    fprintf(stderr, "  -O ovlStore      overlap store\n");
    fprintf(stderr, "  -T tigStore v    tigStore and version\n");
    fprintf(stderr, "  -C ckpname v     scaffolder checkpoint file name and version\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Outputs:\n");
    fprintf(stderr, "  -o prefix        write the output here\n");
    fprintf(stderr, "  -h               print help\n");

    fprintf(stderr, "\n");
    exit(1);
  }


  if (gkpName) {
    generateStatisticsFromGKPStore(argv[0], gkpName);
    exit(0);
  }

  if (ovlName) {
    //generateStatisticsFromOVLStore(argv[0], ovlName);
    exit(0);
  }

  if (tigName && ckpName) {
    //generateStatisticsFromStores(gkpName, ovlName, tigName, tigVersion, ckpName, ckpVersion);
    exit(0);
  }

  if (tigName) {
    //generateStatisticsFromStores(gkpName, ovlName, tigName, tigVersion, ckpName, ckpVersion);
    exit(0);
  }

  if (asmName) {
    generateStatisticsFromASMFile(argv[0], asmName);
    exit(0);
  }


  //writeStatistics();

  //  NOTE!  The "N=%06.2f" format is IMPORTANT.  This will print N=100.00 and N=059.32, both of
  //  which have no whitespace between the N= and the number.

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Scaffolds]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "total scaffolds            = "F10U"\n", getNum(allScaffolds_contigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< contigs in scaffolds >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in scaffolds         = "F10U"\n", getSum(allScaffolds_bases));
  fprintf(stdout, "ave bases in scaffold      = "F10U"\n", getAve(allScaffolds_bases));
  fprintf(stdout, "min bases in scaffold      = "F10U"\n", getMin(allScaffolds_bases));
  fprintf(stdout, "max bases in scaffold      = "F10U"\n", getMax(allScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(25, allScaffolds_bases), getNfrac(25, allScaffolds_bases), getNidx(25, allScaffolds_bases));
  fprintf(stdout, "N50 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(50, allScaffolds_bases), getNfrac(50, allScaffolds_bases), getNidx(50, allScaffolds_bases));
  fprintf(stdout, "N75 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(75, allScaffolds_bases), getNfrac(75, allScaffolds_bases), getNidx(75, allScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "span of scaffolds          = "F10U"\n", getSum(allScaffolds_span));
  fprintf(stdout, "ave span of scaffold       = "F10U"\n", getAve(allScaffolds_span));
  fprintf(stdout, "min span of scaffold       = "F10U"\n", getMin(allScaffolds_span));
  fprintf(stdout, "max span of scaffold       = "F10U"\n", getMax(allScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(25, allScaffolds_span), getNfrac(25, allScaffolds_span), getNidx(25, allScaffolds_span));
  fprintf(stdout, "N50 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(50, allScaffolds_span), getNfrac(50, allScaffolds_span), getNidx(50, allScaffolds_span));
  fprintf(stdout, "N75 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(75, allScaffolds_span), getNfrac(75, allScaffolds_span), getNidx(75, allScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "intra-scaffold gaps        = "F10U"\n", getNum(allScaffolds_gaps));
  fprintf(stdout, "ave intra-scaffold gap     = "F10U"\n", getAve(allScaffolds_gaps));
  fprintf(stdout, "min intra-scaffold gap     = "F10U"\n", getMin(allScaffolds_gaps));
  fprintf(stdout, "max intra-scaffold gap     = "F10U"\n", getMax(allScaffolds_gaps));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< small scaffolds >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in scaffolds         = "F10U"\n", getSum(smallScaffolds_bases));
  fprintf(stdout, "ave bases in scaffold      = "F10U"\n", getAve(smallScaffolds_bases));
  fprintf(stdout, "min bases in scaffold      = "F10U"\n", getMin(smallScaffolds_bases));
  fprintf(stdout, "max bases in scaffold      = "F10U"\n", getMax(smallScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(25, smallScaffolds_bases), getNfrac(25, smallScaffolds_bases), getNidx(25, smallScaffolds_bases));
  fprintf(stdout, "N50 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(50, smallScaffolds_bases), getNfrac(50, smallScaffolds_bases), getNidx(50, smallScaffolds_bases));
  fprintf(stdout, "N75 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(75, smallScaffolds_bases), getNfrac(75, smallScaffolds_bases), getNidx(75, smallScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "span of scaffolds          = "F10U"\n", getSum(smallScaffolds_span));
  fprintf(stdout, "ave span of scaffold       = "F10U"\n", getAve(smallScaffolds_span));
  fprintf(stdout, "min span of scaffold       = "F10U"\n", getMin(smallScaffolds_span));
  fprintf(stdout, "max span of scaffold       = "F10U"\n", getMax(smallScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(25, smallScaffolds_span), getNfrac(25, smallScaffolds_span), getNidx(25, smallScaffolds_span));
  fprintf(stdout, "N50 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(50, smallScaffolds_span), getNfrac(50, smallScaffolds_span), getNidx(50, smallScaffolds_span));
  fprintf(stdout, "N75 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(75, smallScaffolds_span), getNfrac(75, smallScaffolds_span), getNidx(75, smallScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "intra-scaffold gaps        = "F10U"\n", getNum(smallScaffolds_gaps));
  fprintf(stdout, "ave intra-scaffold gap     = "F10U"\n", getAve(smallScaffolds_gaps));
  fprintf(stdout, "min intra-scaffold gap     = "F10U"\n", getMin(smallScaffolds_gaps));
  fprintf(stdout, "max intra-scaffold gap     = "F10U"\n", getMax(smallScaffolds_gaps));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< large scaffolds >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in scaffolds         = "F10U"\n", getSum(largeScaffolds_bases));
  fprintf(stdout, "ave bases in scaffold      = "F10U"\n", getAve(largeScaffolds_bases));
  fprintf(stdout, "min bases in scaffold      = "F10U"\n", getMin(largeScaffolds_bases));
  fprintf(stdout, "max bases in scaffold      = "F10U"\n", getMax(largeScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(25, largeScaffolds_bases), getNfrac(25, largeScaffolds_bases), getNidx(25, largeScaffolds_bases));
  fprintf(stdout, "N50 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(50, largeScaffolds_bases), getNfrac(50, largeScaffolds_bases), getNidx(50, largeScaffolds_bases));
  fprintf(stdout, "N75 scaffold bases         = "F10U", N=%06.2f I="F00U"\n", getN(75, largeScaffolds_bases), getNfrac(75, largeScaffolds_bases), getNidx(75, largeScaffolds_bases));
  fprintf(stdout, "\n");
  fprintf(stdout, "span of scaffolds          = "F10U"\n", getSum(largeScaffolds_span));
  fprintf(stdout, "ave span of scaffold       = "F10U"\n", getAve(largeScaffolds_span));
  fprintf(stdout, "min span of scaffold       = "F10U"\n", getMin(largeScaffolds_span));
  fprintf(stdout, "max span of scaffold       = "F10U"\n", getMax(largeScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "N25 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(25, largeScaffolds_span), getNfrac(25, largeScaffolds_span), getNidx(25, largeScaffolds_span));
  fprintf(stdout, "N50 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(50, largeScaffolds_span), getNfrac(50, largeScaffolds_span), getNidx(50, largeScaffolds_span));
  fprintf(stdout, "N75 scaffold span          = "F10U", N=%06.2f I="F00U"\n", getN(75, largeScaffolds_span), getNfrac(75, largeScaffolds_span), getNidx(75, largeScaffolds_span));
  fprintf(stdout, "\n");
  fprintf(stdout, "intra-scaffold gaps        = "F10U"\n", getNum(largeScaffolds_gaps));
  fprintf(stdout, "ave intra-scaffold gap     = "F10U"\n", getAve(largeScaffolds_gaps));
  fprintf(stdout, "min intra-scaffold gap     = "F10U"\n", getMin(largeScaffolds_gaps));
  fprintf(stdout, "max intra-scaffold gap     = "F10U"\n", getMax(largeScaffolds_gaps));
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Contigs]\n");
  fprintf(stdout, "\n");

  uint64 basesInContigs = getSum(allContigs);

  fprintf(stdout, "<< all contigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num contigs                = "F10U"\n", getNum(allContigs));
  fprintf(stdout, "bases in contigs           = "F10U" %7.3f%%\n", getSum(allContigs), 100.0 * getSum(allContigs) / basesInContigs);
  fprintf(stdout, "bases in contigs ave       = "F10U"\n", getAve(allContigs));
  fprintf(stdout, "bases in contigs min       = "F10U"\n", getMin(allContigs));
  fprintf(stdout, "bases in contigs max       = "F10U"\n", getMax(allContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in contigs N25       = "F10U", N=%06.2f I="F00U"\n", getN(25, allContigs), getNfrac(25, allContigs), getNidx(25, allContigs));
  fprintf(stdout, "bases in contigs N50       = "F10U", N=%06.2f I="F00U"\n", getN(50, allContigs), getNfrac(50, allContigs), getNidx(50, allContigs));
  fprintf(stdout, "bases in contigs N75       = "F10U", N=%06.2f I="F00U"\n", getN(75, allContigs), getNfrac(75, allContigs), getNidx(75, allContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< contigs in scaffolds >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num contigs                = "F10U"\n", getNum(scaffContigs));
  fprintf(stdout, "bases in contigs           = "F10U" %7.3f%%\n", getSum(scaffContigs), 100.0 * getSum(scaffContigs) / basesInContigs);
  fprintf(stdout, "bases in contigs ave       = "F10U"\n", getAve(scaffContigs));
  fprintf(stdout, "bases in contigs min       = "F10U"\n", getMin(scaffContigs));
  fprintf(stdout, "bases in contigs max       = "F10U"\n", getMax(scaffContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in contigs N25       = "F10U", N=%06.2f I="F00U"\n", getN(25, scaffContigs), getNfrac(25, scaffContigs), getNidx(25, scaffContigs));
  fprintf(stdout, "bases in contigs N50       = "F10U", N=%06.2f I="F00U"\n", getN(50, scaffContigs), getNfrac(50, scaffContigs), getNidx(50, scaffContigs));
  fprintf(stdout, "bases in contigs N75       = "F10U", N=%06.2f I="F00U"\n", getN(75, scaffContigs), getNfrac(75, scaffContigs), getNidx(75, scaffContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< small contigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num contigs                = "F10U"\n", getNum(smallContigs));
  fprintf(stdout, "bases in contigs           = "F10U" %7.3f%%\n", getSum(smallContigs), 100.0 * getSum(smallContigs) / basesInContigs);
  fprintf(stdout, "bases in contigs ave       = "F10U"\n", getAve(smallContigs));
  fprintf(stdout, "bases in contigs min       = "F10U"\n", getMin(smallContigs));
  fprintf(stdout, "bases in contigs max       = "F10U"\n", getMax(smallContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in contigs N25       = "F10U", N=%06.2f I="F00U"\n", getN(25, smallContigs), getNfrac(25, smallContigs), getNidx(25, smallContigs));
  fprintf(stdout, "bases in contigs N50       = "F10U", N=%06.2f I="F00U"\n", getN(50, smallContigs), getNfrac(50, smallContigs), getNidx(50, smallContigs));
  fprintf(stdout, "bases in contigs N75       = "F10U", N=%06.2f I="F00U"\n", getN(75, smallContigs), getNfrac(75, smallContigs), getNidx(75, smallContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< large contigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num contigs                = "F10U"\n", getNum(largeContigs));
  fprintf(stdout, "bases in contigs           = "F10U" %7.3f%%\n", getSum(largeContigs), 100.0 * getSum(largeContigs) / basesInContigs);
  fprintf(stdout, "bases in contigs ave       = "F10U"\n", getAve(largeContigs));
  fprintf(stdout, "bases in contigs min       = "F10U"\n", getMin(largeContigs));
  fprintf(stdout, "bases in contigs max       = "F10U"\n", getMax(largeContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in contigs N25       = "F10U", N=%06.2f I="F00U"\n", getN(25, largeContigs), getNfrac(25, largeContigs), getNidx(25, largeContigs));
  fprintf(stdout, "bases in contigs N50       = "F10U", N=%06.2f I="F00U"\n", getN(50, largeContigs), getNfrac(50, largeContigs), getNidx(50, largeContigs));
  fprintf(stdout, "bases in contigs N75       = "F10U", N=%06.2f I="F00U"\n", getN(75, largeContigs), getNfrac(75, largeContigs), getNidx(75, largeContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< degenerate contigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num contigs                = "F10U"\n", getNum(degenContigs));
  fprintf(stdout, "bases in contigs           = "F10U" %7.3f%%\n", getSum(degenContigs), 100.0 * getSum(degenContigs) / basesInContigs);
  fprintf(stdout, "bases in contigs ave       = "F10U"\n", getAve(degenContigs));
  fprintf(stdout, "bases in contigs min       = "F10U"\n", getMin(degenContigs));
  fprintf(stdout, "bases in contigs max       = "F10U"\n", getMax(degenContigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "bases in contigs N25       = "F10U", N=%06.2f I="F00U"\n", getN(25, degenContigs), getNfrac(25, degenContigs), getNidx(25, degenContigs));
  fprintf(stdout, "bases in contigs N50       = "F10U", N=%06.2f I="F00U"\n", getN(50, degenContigs), getNfrac(50, degenContigs), getNidx(50, degenContigs));
  fprintf(stdout, "bases in contigs N75       = "F10U", N=%06.2f I="F00U"\n", getN(75, degenContigs), getNfrac(75, degenContigs), getNidx(75, degenContigs));
  fprintf(stdout, "\n");

  uint64  basesInUnitigs = getSum(allUnitigs);

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Unitigs]\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "<< all unitigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num unitigs                = "F10U"\n", getNum(allUnitigs));
  fprintf(stdout, "bases in unitigs           = "F10U" %7.3f%%\n", getSum(allUnitigs), 100.0 * getSum(allUnitigs) / basesInUnitigs);
  fprintf(stdout, "bases in unitigs ave       = "F10U"\n", getAve(allUnitigs));
  fprintf(stdout, "bases in unitigs min       = "F10U"\n", getMin(allUnitigs));
  fprintf(stdout, "bases in unitigs max       = "F10U"\n", getMax(allUnitigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< unique unitigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num unitigs                = "F10U"\n", getNum(uniqueUnitigs));
  fprintf(stdout, "bases in unitigs           = "F10U" %7.3f%%\n", getSum(uniqueUnitigs), 100.0 * getSum(uniqueUnitigs) / basesInUnitigs);
  fprintf(stdout, "bases in unitigs ave       = "F10U"\n", getAve(uniqueUnitigs));
  fprintf(stdout, "bases in unitigs min       = "F10U"\n", getMin(uniqueUnitigs));
  fprintf(stdout, "bases in unitigs max       = "F10U"\n", getMax(uniqueUnitigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< not rez unitigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num unitigs                = "F10U"\n", getNum(notRezUnitigs));
  fprintf(stdout, "bases in unitigs           = "F10U" %7.3f%%\n", getSum(notRezUnitigs), 100.0 * getSum(notRezUnitigs) / basesInUnitigs);
  fprintf(stdout, "bases in unitigs ave       = "F10U"\n", getAve(notRezUnitigs));
  fprintf(stdout, "bases in unitigs min       = "F10U"\n", getMin(notRezUnitigs));
  fprintf(stdout, "bases in unitigs max       = "F10U"\n", getMax(notRezUnitigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< surrogate unitigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num unitigs                = "F10U"\n", getNum(surrogateUnitigs));
  fprintf(stdout, "bases in unitigs           = "F10U" %7.3f%%\n", getSum(surrogateUnitigs), 100.0 * getSum(surrogateUnitigs) / basesInUnitigs);
  fprintf(stdout, "bases in unitigs ave       = "F10U"\n", getAve(surrogateUnitigs));
  fprintf(stdout, "bases in unitigs min       = "F10U"\n", getMin(surrogateUnitigs));
  fprintf(stdout, "bases in unitigs max       = "F10U"\n", getMax(surrogateUnitigs));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< other unitigs >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "num unitigs                = "F10U"\n", getNum(otherUnitigs));
  fprintf(stdout, "bases in unitigs           = "F10U" %7.3f%%\n", getSum(otherUnitigs), 100.0 * getSum(otherUnitigs) / basesInUnitigs);
  fprintf(stdout, "bases in unitigs ave       = "F10U"\n", getAve(otherUnitigs));
  fprintf(stdout, "bases in unitigs min       = "F10U"\n", getMin(otherUnitigs));
  fprintf(stdout, "bases in unitigs max       = "F10U"\n", getMax(otherUnitigs));
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Surrogates]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Polymorphism]\n");
  fprintf(stdout, "\n");

  uint64  readsWithMate = totalReadsFromAFG - readsWithNoMate;

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Mate Pairs]\n");
  fprintf(stdout, "reads with good            = "F10U" %7.3f%%\n", readsWithGoodMate,          100.0 * readsWithGoodMate          / readsWithMate);
  fprintf(stdout, "reads with bad short       = "F10U" %7.3f%%\n", readsWithBadShortMate,      100.0 * readsWithBadShortMate      / readsWithMate);
  fprintf(stdout, "reads with bad long        = "F10U" %7.3f%%\n", readsWithBadLongMate,       100.0 * readsWithBadLongMate       / readsWithMate);
  fprintf(stdout, "reads with same orient     = "F10U" %7.3f%%\n", readsWithSameOrientMate,    100.0 * readsWithSameOrientMate    / readsWithMate);
  fprintf(stdout, "reads with outtie          = "F10U" %7.3f%%\n", readsWithOuttieMate,        100.0 * readsWithOuttieMate        / readsWithMate);
  fprintf(stdout, "reads with both chaff      = "F10U" %7.3f%%\n", readsWithBothChaffMate,     100.0 * readsWithBothChaffMate     / readsWithMate);
  fprintf(stdout, "reads with chaff           = "F10U" %7.3f%%\n", readsWithChaffMate,         100.0 * readsWithChaffMate         / readsWithMate);
  fprintf(stdout, "reads with both degen      = "F10U" %7.3f%%\n", readsWithBothDegenMate,     100.0 * readsWithBothDegenMate     / readsWithMate);
  fprintf(stdout, "reads with degen           = "F10U" %7.3f%%\n", readsWithDegenMate,         100.0 * readsWithDegenMate         / readsWithMate);
  fprintf(stdout, "reads with surrogate       = "F10U" %7.3f%%\n", readsWithBothSurrogateMate, 100.0 * readsWithBothSurrogateMate / readsWithMate);
  fprintf(stdout, "reads with both surrogate  = "F10U" %7.3f%%\n", readsWithSurrogateMate,     100.0 * readsWithSurrogateMate     / readsWithMate);
  fprintf(stdout, "reads with diff scaff      = "F10U" %7.3f%%\n", readsWithDiffScafMate,      100.0 * readsWithDiffScafMate      / readsWithMate);
  fprintf(stdout, "reads with unassigned      = "F10U" %7.3f%%\n", readsWithUnassignedMate,    100.0 * readsWithUnassignedMate    / readsWithMate);
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Reads]\n");
  fprintf(stdout, "unmated reads              = "F10U" %7.3f%%\n", readsWithNoMate, 100.0 * readsWithNoMate / totalReadsFromAFG);
  fprintf(stdout, "mated reads                = "F10U" %7.3f%%\n", readsWithMate,   100.0 * readsWithMate   / totalReadsFromAFG);
  fprintf(stdout, "\n");
  fprintf(stdout, "<< all reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(all_reads), 100.0 * getNum(all_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(all_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(all_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< all contig reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(allContigs_reads), 100.0 * getNum(allContigs_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(allContigs_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(allContigs_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< large contig reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(largeContigs_reads), 100.0 * getNum(largeContigs_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(largeContigs_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(largeContigs_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< small contig reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(smallContigs_reads), 100.0 * getNum(smallContigs_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(smallContigs_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(smallContigs_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< surrogate reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(surrogate_reads), 100.0 * getNum(surrogate_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(surrogate_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(surrogate_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< degenerate reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(degenContigs_reads), 100.0 * getNum(degenContigs_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(degenContigs_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(degenContigs_reads));
  fprintf(stdout, "\n");
  fprintf(stdout, "<< singleton reads >>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "number of reads            = "F10U" %7.3f%%\n", getNum(singleton_reads), 100.0 * getNum(singleton_reads) / totalReadsFromAFG);
  fprintf(stdout, "sum clear range            = "F10U"\n",         getSum(singleton_reads));
  fprintf(stdout, "ave clear range            = "F10U"\n",         getAve(singleton_reads));
  fprintf(stdout, "\n");


  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Libraries]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[GC Content]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Read Base Counts]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "========================================\n");
  fprintf(stdout, "[Depth]\n");
  fprintf(stdout, "\n");

  fprintf(stdout, "\n");

  exit(0);
}
