
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#include "AS_UTL_decodeRange.H"
#include "AS_UTL_intervalList.H"

#include "MultiAlign.H"
#include "MultiAlignStore.H"
#include "MultiAlignMatePairAnalysis.H"
#include "MultiAlignSizeAnalysis.H"
#include "MultiAlignment_CNS.H"
#include "MultiAlignment_CNS_private.H"

#define DUMP_PROPERTIES       1
#define DUMP_FRAGS            2
#define DUMP_UNITIGS          3
#define DUMP_CONSENSUS        4
#define DUMP_CONSENSUSGAPPED  5
#define DUMP_LAYOUT           6
#define DUMP_MULTIALIGN       7
#define DUMP_MATEPAIR         8
#define DUMP_SIZES            9
#define DUMP_COVERAGE        10
#define DUMP_THINOVERLAP     11
#define DUMP_FMAP            12

#define OPERATION_UNITIGLIST  1
#define OPERATION_CONTIGLIST  2
#define OPERATION_PROPERTIES  3
#define OPERATION_TIG         4
#define OPERATION_EDIT        5
#define OPERATION_REPLACE     6
#define OPERATION_BUILD       7
#define OPERATION_COMPRESS    8

void
changeProperties(MultiAlignStore *tigStore,
                 char            *editName) {
  char  editLine[1024];

  errno = 0;
  FILE *editFile = fopen(editName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", editName, strerror(errno)), exit(1);

  fgets(editLine, 1024, editFile);
  while (!feof(editFile)) {
    chomp(editLine);

    //  Decode the string into three pieces, the operation, the tig id and the value.  Note that
    //  these are just pointers to the start of the piece, and that the pieces are NOT
    //  zero-terminated (so we can print an error of the whole line if needed).
    //
    char   *op = editLine;
    char   *tp = editLine;
    char   *vp = editLine;

    //  Skip whitespace before the operation
    while (isspace(*op))
      op++;

    //  Skip the operation, then skip whitespace before the tig.
    tp = op;
    while (!isspace(*tp))
      tp++;
    while (isspace(*tp))
      tp++;

    //  Skip the tig, then skip whitespace before the value.
    vp = tp;
    while (!isspace(*vp))
      vp++;
    while (isspace(*vp))
      vp++;

    int32   tid = atoi(tp);

    if        (strncmp(op, "unitig_coverage_stat", 20) == 0) {
      tigStore->setUnitigCoverageStat(tid, atof(vp));

    } else if (strncmp(op, "unitig_microhet_prob", 20) == 0) {
      tigStore->setUnitigMicroHetProb(tid, atof(vp));

    } else if (strncmp(op, "unitig_status", 13) == 0) {
      UnitigStatus  st = tigStore->getUnitigStatus(tid);
      switch (*vp) {
        case AS_UNIQUE:      st = AS_UNIQUE;      break;
        case AS_NOTREZ:      st = AS_NOTREZ;      break;
        case AS_SEP:         st = AS_SEP;         break;
        case AS_UNASSIGNED:  st = AS_UNASSIGNED;  break;
        default:
          fprintf(stderr, "unknown unitig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigStatus(tid, st);

    } else if (strncmp(op, "unitig_unique_rept", 18) == 0) {
      UnitigFUR ur = tigStore->getUnitigFUR(tid);
      switch (*vp) {
        case AS_FORCED_NONE:    ur = AS_FORCED_NONE;    break;
        case AS_FORCED_UNIQUE:  ur = AS_FORCED_UNIQUE;  break;
        case AS_FORCED_REPEAT:  ur = AS_FORCED_REPEAT;  break;
        default:
          fprintf(stderr, "unknown unitig_unique_rept in '%s'\n", editLine);
          break;
      }
      tigStore->setUnitigFUR(tid, ur);

    } else if (strncmp(op, "contig_status", 13) == 0) {
      ContigStatus  st = tigStore->getContigStatus(tid);
      switch (*vp) {
        case AS_PLACED:      st = AS_PLACED;    break;
        case AS_UNPLACED:    st = AS_UNPLACED;  break;
        default:
          fprintf(stderr, "unknown contig_status in '%s'\n", editLine);
          break;
      }
      tigStore->setContigStatus(tid, st);

    } else {
      fprintf(stderr, "unknown edit '%s'\n", editLine);
    }

    fgets(editLine, 1024, editFile);
  }
}



void
dumpProperties(MultiAlignStore *tigStore,
               int32 tigID,
               int32 tigIsUnitig,
               MultiAlignT *ma) {

  fprintf(stdout, "maID                "F_S32"\n", ma->maID);
  fprintf(stdout, "unitigCoverageStat  %f\n",      ma->data.unitig_coverage_stat);
  fprintf(stdout, "unitigMicrohetProb  %f\n",      ma->data.unitig_microhet_prob);
  fprintf(stdout, "unitigStatus        %c/%d\n",   ma->data.unitig_status, ma->data.unitig_status);
  fprintf(stdout, "unitigFUR           %c/%d\n",   ma->data.unitig_unique_rept, ma->data.unitig_unique_rept);
  fprintf(stdout, "contigStatus        %c/%d\n",   ma->data.contig_status, ma->data.contig_status);

#if GCCONTENT
  float gcContent = 0.0;
  int ulen = 0;
  int glen = 0;

  if (ma->consensus) {
    char *cns = Getchar(ma->consensus, 0);

    if (cns && *cns) {
      int gcs = 0;

      while (*cns) {
        glen++;

        if (*cns != '-') {
          if (*cns == 'G' || *cns == 'C') { gcs++; }
          ulen++;
        }

        cns++;
      }

      if (ulen) { gcContent = (float) gcs / (float) ulen; }
    }
  }

  fprintf(stdout, "gcContent           %f\n",      gcContent);
  fprintf(stdout, "uLen                %d\n",      ulen);
  fprintf(stdout, "gLen                %d\n",      glen);
#endif

  fprintf(stdout, "numFrags            "F_U32" (vs "F_U64")\n", ma->data.num_frags, (uint64)GetNumIntMultiPoss(ma->f_list));
  fprintf(stdout, "numUnitigs          "F_U32" (vs "F_U64")\n", ma->data.num_unitigs, (uint64)GetNumIntUnitigPoss(ma->u_list));

  tigStore->dumpMultiAlignR(tigID, tigIsUnitig);
}


void
dumpFrags(MultiAlignStore *tigStore,
          int32 tigID,
          int32 tigIsUnitig,
          MultiAlignT *ma) {

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    fprintf(stdout, "FRG %7d %5d,%5d\n",
            imp->ident, imp->position.bgn, imp->position.end);
  }
}


void
dumpUnitigs(MultiAlignStore *tigStore,
            int32 tigID,
            int32 tigIsUnitig,
            MultiAlignT *ma) {

  for (uint32 i=0; i<GetNumIntUnitigPoss(ma->u_list); i++) {
    IntUnitigPos *iup = GetIntUnitigPos(ma->u_list, i);

    fprintf(stdout, "UTG %7d %5d,%5d\n",
            iup->ident, iup->position.bgn, iup->position.end);
  }
}


void
dumpConsensus(MultiAlignStore *tigStore,
              int32            tigID,
              int32            tigIsUnitig,
              MultiAlignT     *ma,
              bool             withGaps,
              uint32           minCoverage) {

  if (ma->consensus == NULL)
    return;

  char   *cns    = Getchar(ma->consensus, 0);

  if ((cns == NULL) || (cns[0] == 0))
    return;

  //  If a minCoverage is specified, convert the low coverage bases to underscores, which will be
  //  filtered later.

  if (minCoverage > 0) {
    intervalList  allL;

    for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

      int32   bgn = MIN(imp->position.bgn, imp->position.end);
      int32   end = MAX(imp->position.bgn, imp->position.end);

      allL.add(bgn, end - bgn);
    }

    intervalDepth  ID(allL);

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
      if (ID.de(ii) >= minCoverage)
        continue;

      for (uint32 pp=ID.lo(ii); pp<ID.hi(ii); pp++)
        cns[pp] = '_';
    }
  }

  //  Now filter out gaps in the consensus.

  if (withGaps == false) {
    char *o = cns;
    char *n = cns;

    while (*n) {
      if (*n != '-')
        *o++ = *n;
      n++;
    }

    *o = 0;
  }

  //  If no min coverage, we can just dump the consensus and be done.  Plus we output a few bits of
  //  useful info with the sequence.

  if (minCoverage == 0) {
    if (tigIsUnitig)
      fprintf(stdout, ">utg%d len="F_U64" reads="F_U32" status=%c microHet=%.2f covStat=%.2f\n%s\n",
              ma->maID, GetNumchars(ma->consensus) - 1, ma->data.num_frags,
              ma->data.unitig_status,
              ma->data.unitig_microhet_prob,
              ma->data.unitig_coverage_stat,
              cns);
    else
      fprintf(stdout, ">ctg%d len="F_U64" reads="F_U32" unitigs="F_U32" status=%c\n%s\n",
              ma->maID, GetNumchars(ma->consensus) - 1, ma->data.num_frags,
              ma->data.num_unitigs,
              ma->data.contig_status,
              cns);
    return;
  }

  //  Otherwise, we need to find subsequences in the consensus.  The useful bits of info aren't
  //  valid anymore.

  uint32  cnsLen = strlen(cns);
  uint32  part   = 0;

  for (uint32 bgn=0; bgn<cnsLen; bgn++)
    if (cns[bgn] == '_')
      cns[bgn] = 0;

  for (uint32 bgn=0; bgn<cnsLen; bgn++) {
    while ((cns[bgn] == 0) && (bgn < cnsLen))
      bgn++;

    if (bgn >= cnsLen)
      break;

    uint32 end = bgn + 1;

    while ((cns[end] != 0) && (end < cnsLen))
      end++;

    fprintf(stdout, ">%s%d.%u bgn=%u end=%u len=%u\n%s\n",
            (tigIsUnitig) ? "utg" : "cns",
            ma->maID, part, bgn, end, end-bgn,
            cns + bgn);

    bgn = end + 1;

    part++;
  }
}



void
dumpCoverage(MultiAlignStore *tigStore,
             int32            tigID,
             int32            tigIsUnitig,
             MultiAlignT     *ma,
             uint32           minCoverage,
             uint32           maxCoverage,
             uint64          *coverageHistogram,
             uint32           coverageHistogramLen,
             char            *outPrefix) {
  intervalList  allL;

  uint32        maxPos = 0;

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    int32   bgn = MIN(imp->position.bgn, imp->position.end);
    int32   end = MAX(imp->position.bgn, imp->position.end);
    int32   len = end - bgn;

    if (maxPos < end)
      maxPos = end;

    allL.add(bgn, len);
  }

  maxPos++;  //  Now the C-style maxPos.

  intervalDepth  ID(allL);

  intervalList   minL;
  intervalList   maxL;

  uint32  maxDepth    = 0;
  double  aveDepth    = 0;
  double  sdeDepth    = 0;

#if 0
  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    if ((ID.de(ii) < minCoverage) && (ID.lo(ii) != 0) && (ID.hi(ii) != maxPos)) {
      fprintf(stderr, "%s %d low coverage interval %ld %ld max %u coverage %u\n",
              (tigIsUnitig) ? "unitig" : "contig", tigID, ID.lo(ii), ID.hi(ii), maxPos, ID.de(ii));
      minL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
    }

    if (maxCoverage <= ID.de(ii)) {
      fprintf(stderr, "%s %d high coverage interval %ld %ld max %u coverage %u\n",
              (tigIsUnitig) ? "unitig" : "contig", tigID, ID.lo(ii), ID.hi(ii), maxPos, ID.de(ii));
      maxL.add(ID.lo(ii), ID.hi(ii) - ID.lo(ii) + 1);
    }
  }
#endif

  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    if (ID.de(ii) > maxDepth)
      maxDepth = ID.de(ii);

    aveDepth += (ID.hi(ii) - ID.lo(ii) + 1) * ID.de(ii);

    if (ID.de(ii) < coverageHistogramLen)
      coverageHistogram[ID.de(ii)] += ID.hi(ii) - ID.lo(ii) + 1;
    else
      fprintf(stderr, "deep coverage %d\n", ID.de(ii));
  }

  aveDepth /= maxPos;

  for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++) {
    sdeDepth += (ID.hi(ii) - ID.lo(ii) + 1) * (ID.de(ii) - aveDepth) * (ID.de(ii) - aveDepth);
  }

  sdeDepth = sqrt(sdeDepth / maxPos);

  if (maxDepth > 1000)
    fprintf(stderr, "DEEP unitig %u of length %u with maxDepth %u\n",
            tigID, maxPos, maxDepth);

  allL.merge();
  minL.merge();
  maxL.merge();


#if 0
  if      ((minL.numberOfIntervals() > 0) && (maxL.numberOfIntervals() > 0))
    fprintf(stderr, "%s %d has %u intervals, %u regions below %u coverage and %u regions at or above %u coverage\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            allL.numberOfIntervals(),
            minL.numberOfIntervals(), minCoverage,
            maxL.numberOfIntervals(), maxCoverage);
  else if (minL.numberOfIntervals() > 0)
    fprintf(stderr, "%s %d has %u intervals, %u regions below %u coverage\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            allL.numberOfIntervals(),
            minL.numberOfIntervals(), minCoverage);
  else if (maxL.numberOfIntervals() > 0)
    fprintf(stderr, "%s %d has %u intervals, %u regions at or above %u coverage\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            allL.numberOfIntervals(),
            maxL.numberOfIntervals(), maxCoverage);
  else
    fprintf(stderr, "%s %d has %u intervals\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            allL.numberOfIntervals());
#endif

  if (outPrefix) {
    char  outName[FILENAME_MAX];

    sprintf(outName, "%s.%s%08u.depth", outPrefix, (tigIsUnitig) ? "utg" : "ctg", tigID);

    FILE *outFile = fopen(outName, "w");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

    for (uint32 ii=0; ii<ID.numberOfIntervals(); ii++)
      fprintf(outFile, "%lu\t%lu\t%u\n", ID.lo(ii), ID.hi(ii), ID.de(ii));

    fclose(outFile);

    FILE *gnuPlot = popen("gnuplot > /dev/null 2>&1", "w");

    if (gnuPlot) {
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "set output '%s.%s%08u.png'\n", outPrefix, (tigIsUnitig) ? "utg" : "ctg", tigID);
      fprintf(gnuPlot, "set xlabel 'position'\n");
      fprintf(gnuPlot, "set ylabel 'coverage'\n");
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "plot '%s.%s%08u.depth' using 1:2 with lines title '%s %u length %u', \\\n",
              outPrefix,
              (tigIsUnitig) ? "utg" : "ctg",
              tigID,
              (tigIsUnitig) ? "unitig" : "contig", tigID, maxPos);
      fprintf(gnuPlot, "     %f title 'mean %.2f +- %.2f', \\\n", aveDepth, aveDepth, sdeDepth);
      fprintf(gnuPlot, "     %f title '' lt 0 lc 2, \\\n", aveDepth - sdeDepth);
      fprintf(gnuPlot, "     %f title '' lt 0 lc 2\n",     aveDepth + sdeDepth);

      fclose(gnuPlot);
    }
  }
}


void
dumpThinOverlap(MultiAlignStore *tigStore,
                int32            tigID,
                int32            tigIsUnitig,
                MultiAlignT     *ma,
                uint32           minOverlap) {
  intervalList  allL;
  intervalList  ovlL;

  uint32        maxPos = 0;

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    int32   bgn = MIN(imp->position.bgn, imp->position.end);
    int32   end = MAX(imp->position.bgn, imp->position.end);
    int32   len = end - bgn + 1;

    if (maxPos < end)
      maxPos = end;

    allL.add(bgn, len);
    ovlL.add(bgn, len);
  }

  allL.merge();
  ovlL.merge(minOverlap);

  if (ovlL.numberOfIntervals() == 1)
    return;


  if (maxPos < 10000)
    return;

#if 1
  for (uint32 i=0; i<ovlL.numberOfIntervals(); i++)
    fprintf(stderr, "%s %u IL "F_U64" "F_U64"\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            ovlL.lo(i), ovlL.hi(i));
#endif

  intervalList badL;

  for (uint32 i=1; i<ovlL.numberOfIntervals(); i++) {
    assert(ovlL.lo(i) < ovlL.hi(i-1));

    badL.add(ovlL.lo(i), ovlL.hi(i-1) - ovlL.lo(i));
    //badL.add(ovlL.hi(i-1), ovlL.lo(i) - ovlL.hi(i-1));
  }

#if 1
  for (uint32 i=0; i<badL.numberOfIntervals(); i++)
    fprintf(stderr, "%s %u BAD "F_U64" "F_U64"\n",
            (tigIsUnitig) ? "unitig" : "contig", tigID,
            badL.lo(i), badL.hi(i));
#endif

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    int32   bgn = MIN(imp->position.bgn, imp->position.end);
    int32   end = MAX(imp->position.bgn, imp->position.end);
    int32   len = end - bgn + 1;

    bool    report = false;

    for (uint32 oo=0; oo<badL.numberOfIntervals(); oo++)
      if ((badL.lo(oo) <= end) &&
          (bgn         <= badL.hi(oo))) {
        report = true;
        break;
      }

    if (report)
      fprintf(stderr, "%s %d frag %u %u-%u\n",
              (tigIsUnitig) ? "unitig" : "contig", tigID,
              imp->ident, imp->position.bgn, imp->position.end);
  }

  fprintf(stderr, "%s %d max %u has %u intervals, %u enforcing minimum overlap of %u\n",
          (tigIsUnitig) ? "unitig" : "contig", tigID, maxPos,
          allL.numberOfIntervals(),
          ovlL.numberOfIntervals(), minOverlap);
}





void
operationBuild(char *buildName, char *tigName,  int tigVers) {
  uint32  utgID = 0;
  uint32  ctgID = 0;
  uint32  orgID = 0;

  MultiAlignStore  *tigStore = NULL;
  MultiAlignT       *ma       = CreateEmptyMultiAlignT();
  bool               isUnitig = false;

  errno = 0;
  FILE         *F = fopen(buildName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", buildName, strerror(errno)), exit(1);

  if (AS_UTL_fileExists(tigName, TRUE, TRUE)) {
    fprintf(stderr, "ERROR: '%s' exists, and I will not clobber an existing store.\n", tigName);
    exit(1);
    tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);
  } else {
    tigStore = new MultiAlignStore(tigName);

    for (int32 v=1; v<tigVers; v++)
      tigStore->nextVersion();
  }

  while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
    if (ma->data.num_frags + ma->data.num_unitigs == 0)
      continue;

    orgID = ma->maID;

    if (isUnitig)
      ma->maID = utgID;
    else
      ma->maID = ctgID;

#if 0
    fprintf(stderr, "INSERTING %s %d (%d frags %d unitigs) (originally ID %d)\n",
            (isUnitig) ? "unitig" : "contig",
            ma->maID,
            ma->data.num_frags, ma->data.num_unitigs,
            orgID);
#endif

    tigStore->insertMultiAlign(ma, isUnitig, FALSE);

    if (isUnitig)
      utgID++;
    else
      ctgID++;
  }

  fclose(F);

  delete tigStore;
}



void
operationCompress(char *tigName, int tigVers) {
  MultiAlignStore  *tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, FALSE, FALSE, FALSE);
  bool              isUnitig = TRUE;

  uint32            nUtgErrors = 0;
  uint32            nCtgErrors = 0;

  uint32            nUtgCompress = 0;
  uint32            nCtgCompress = 0;

  //  Pass 0:  Fail if this isn't the latest version.  If we try to compress something that isn't the
  //  latest version, versions after this still point to the uncompressed tigs.

  //if (tigStore->

  //  Pass 1:  Check that we aren't going to pull a tig out of the future and place it in the past.

  for (uint32 ti=0; ti<tigStore->numUnitigs(); ti++) {
    if (tigStore->isDeleted(ti, isUnitig))
      continue;

    if (tigStore->getUnitigVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future unitig "F_U32" from version "F_U32" to previous version %d.\n",
              ti, tigStore->getUnitigVersion(ti), tigVers);
      nUtgErrors++;
    } else if (tigStore->getUnitigVersion(ti) < tigVers) {
      nUtgCompress++;
    }
  }

  isUnitig = FALSE;

  for (uint32 ti=0; ti<tigStore->numContigs(); ti++) {
    if (tigStore->isDeleted(ti, isUnitig))
      continue;

    if (tigStore->getContigVersion(ti) > tigVers) {
      fprintf(stderr, "WARNING:  Attempt to move future contig "F_U32" from version "F_U32" to previous version %d.\n",
              ti, tigStore->getContigVersion(ti), tigVers);
      nCtgErrors++;
    } else if (tigStore->getContigVersion(ti) < tigVers) {
      nCtgCompress++;
    }
  }

  if (nUtgErrors + nCtgErrors > 0) {
    fprintf(stderr, "Store can't be compressed; probably trying to compress to something that isn't the latest version.\n");
    fprintf(stderr, "  "F_U32" unitigs failed; "F_U32" compressable\n", nUtgErrors, nUtgCompress);
    fprintf(stderr, "  "F_U32" contigs failed; "F_U32" compressable\n", nCtgErrors, nCtgCompress);
    delete tigStore;
    exit(1);
  }

  //  Pass 2:  Actually do the moves

  if (nUtgCompress + nCtgCompress > 0) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, TRUE, FALSE);
  }

  if (nUtgCompress > 0) {
    isUnitig = TRUE;

    fprintf(stderr, "Compressing "F_U32" unitigs into version %d\n", nUtgCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numUnitigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti, isUnitig)) {
        continue;
      }

      if (tigStore->getUnitigVersion(ti) == tigVers)
        continue;

      MultiAlignT *ma = tigStore->loadMultiAlign(ti, isUnitig);

      if (ma == NULL) {
        //fprintf(stderr, "WARNING: unitig "F_U32" is NULL.\n", ti);
        continue;
      }

      tigStore->insertMultiAlign(ma, isUnitig, TRUE);
      tigStore->unloadMultiAlign(ti, isUnitig);
    }
  }

  if (nCtgCompress > 0) {
    isUnitig = FALSE;

    fprintf(stderr, "Compressing "F_U32" contigs into version %d\n", nCtgCompress, tigVers);

    for (uint32 ti=0; ti<tigStore->numContigs(); ti++) {
      if ((ti % 1000000) == 0)
        fprintf(stderr, "tig %d\n", ti);

      if (tigStore->isDeleted(ti, isUnitig))
        continue;

      if (tigStore->getContigVersion(ti) == tigVers)
        continue;

      MultiAlignT *ma = tigStore->loadMultiAlign(ti, isUnitig);

      if (ma == NULL) {
        //fprintf(stderr, "WARNING: contig "F_U32" is NULL.\n", ti);
        continue;
      }

      tigStore->insertMultiAlign(ma, isUnitig, TRUE);
      tigStore->unloadMultiAlign(ti, isUnitig);
    }
  }

  //  Now clean up the older files

  if (nUtgCompress + nCtgCompress > 0) {
    for (uint32 version=1; version<tigVers; version++) {
      fprintf(stderr, "Purge version "F_U32".\n", version);
      tigStore->purgeVersion(version);
    }
  }

  //  And the newer files

  delete tigStore;
}


void
dumpFmap(FILE         *out,
         MultiAlignT  *ma,
         bool          tigIsUnitig) {
  uint32   fiMax = GetNumIntMultiPoss(ma->f_list);

  for (uint32 fi=0; fi<fiMax; fi++) {
    IntMultiPos  *imp = GetIntMultiPos(ma->f_list, fi);

    fprintf(stdout, F_U32"\t"F_U32"\t"F_S32"\t"F_S32"\n",
            imp->ident, ma->maID, imp->position.bgn, imp->position.end);
  }
}




int
main (int argc, char **argv) {
  char          tmpName[FILENAME_MAX] = {0};
  char         *gkpName        = NULL;
  char         *tigName        = NULL;
  int           tigVers        = -1;
  int           tigPartU       = 0;
  int           tigPartC       = 0;
  bool          tigIDset       = false;
  uint32        tigIDbgn       = 0;
  uint32        tigIDend       = UINT32_MAX;
  int32         tigIsUnitig    = TRUE;
  uint32        opType         = 0;
  uint32        dumpFlags      = 0;
  uint32        dumpAll        = 0;
  char         *editName       = NULL;
  char         *replaceName    = NULL;
  bool          sameVersion    = true;
  char         *buildName      = NULL;

  uint32        minNreads      = 0;
  uint32        maxNreads      = UINT32_MAX;

  uint32        minCoverage    = 0;

  MultiAlignT  *ma             = NULL;
  int           showQV         = 0;
  int           showDots       = 1;

  matePairAnalysis  *mpa       = NULL;

  sizeAnalysis      *siz       = NULL;
  uint64             sizSize   = 0;

  uint64            *cov       = NULL;
  uint64             covMax    = 0;

  char              *outPrefix = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-up") == 0) {
      tigPartU = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-cp") == 0) {
      tigPartC = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-u") == 0) {
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigIsUnitig = TRUE;
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-U") == 0) {
      dumpAll     = TRUE;
      tigIsUnitig = TRUE;
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-c") == 0) {
      AS_UTL_decodeRange(argv[++arg], tigIDbgn, tigIDend);
      tigIsUnitig = FALSE;
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-C") == 0) {
      dumpAll     = TRUE;
      tigIsUnitig = FALSE;
      tigIDset    = true;

    } else if (strcmp(argv[arg], "-d") == 0) {
      arg++;

      opType = OPERATION_TIG;

      if      (strcmp(argv[arg], "properties") == 0)
        dumpFlags = DUMP_PROPERTIES;

      else if (strcmp(argv[arg], "frags") == 0)
        dumpFlags = DUMP_FRAGS;

      else if (strcmp(argv[arg], "unitigs") == 0)
        dumpFlags = DUMP_UNITIGS;

      else if (strcmp(argv[arg], "consensus") == 0) {
        dumpFlags = DUMP_CONSENSUS;

        if ((arg+1 < argc) && (isdigit(argv[arg+1][0])))
          minCoverage = atoi(argv[++arg]);
      }

      else if (strcmp(argv[arg], "consensusgapped") == 0)
        dumpFlags = DUMP_CONSENSUSGAPPED;

      else if (strcmp(argv[arg], "layout") == 0)
        dumpFlags = DUMP_LAYOUT;

      else if (strcmp(argv[arg], "multialign") == 0)
        dumpFlags = DUMP_MULTIALIGN;

      else if (strcmp(argv[arg], "matepair") == 0)
        dumpFlags = DUMP_MATEPAIR;

      else if (strcmp(argv[arg], "sizes") == 0)
        dumpFlags = DUMP_SIZES;
      
      else if (strcmp(argv[arg], "coverage") == 0)
        dumpFlags = DUMP_COVERAGE;

      else if (strcmp(argv[arg], "overlap") == 0)
        dumpFlags = DUMP_THINOVERLAP;

      else if (strcmp(argv[arg], "fmap") == 0)
        dumpFlags = DUMP_FMAP;

      else
        fprintf(stderr, "%s: Unknown dump option '-d %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-D") == 0) {
      arg++;

      if      (strcmp(argv[arg], "unitiglist") == 0)
        opType = OPERATION_UNITIGLIST;

      else if (strcmp(argv[arg], "contiglist") == 0)
        opType = OPERATION_CONTIGLIST;

      else if (strcmp(argv[arg], "properties") == 0)
        opType = OPERATION_PROPERTIES;

      else
        fprintf(stderr, "%s: Unknown dump option '-D %s'\n", argv[0], argv[arg]);

    } else if (strcmp(argv[arg], "-E") == 0) {
      opType = OPERATION_EDIT;
      editName = argv[++arg];

    } else if (strcmp(argv[arg], "-B") == 0) {
      opType = OPERATION_BUILD;
      buildName = argv[++arg];


    } else if (strcmp(argv[arg], "-R") == 0) {
      opType = OPERATION_REPLACE;
      replaceName = argv[++arg];

    } else if (strcmp(argv[arg], "-N") == 0) {
      sameVersion = false;

    } else if (strcmp(argv[arg], "-compress") == 0) {
      opType = OPERATION_COMPRESS;

    } else if (strcmp(argv[arg], "-nreads") == 0) {
      minNreads = atoi(argv[++arg]);
      maxNreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-w") == 0) {
      MULTIALIGN_PRINT_WIDTH = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-s") == 0) {
      MULTIALIGN_PRINT_SPACING = atoi(argv[++arg]);
      sizSize                  = atoll(argv[arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else {
      fprintf(stderr, "%s: Unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }

    arg++;
  }
  if ((err) || (gkpName == NULL) || (tigName == NULL)) {
    fprintf(stderr, "usage: %s -g <gkpStore> -t <tigStore> <v> [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -g <gkpStore>           Path to the gatekeeper store\n");
    fprintf(stderr, "  -t <tigStore> <v>       Path to the tigStore, and version to use\n");
    fprintf(stderr, "  -up <p>                 ...limit to unitigs in partition <p>\n");
    fprintf(stderr, "  -cp <p>                 ...limit to contigs in partition <p>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -D <operation>        Dump something about the store\n");
    fprintf(stderr, "     unitiglist         ...a list of the unitigs in the store\n");
    fprintf(stderr, "     contiglist         ...a list of the contigs in the store\n");
    fprintf(stderr, "     properties         ...a list of properties for ALL multialigns in the store (for -E)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -u id[-id]            Unitig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "  -c id[-id]            Contig to dump (for -d option); if A-B, dump tigs from id A to id B, inclusive\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U                    Dump ALL unitigs (for -d option)\n");
    fprintf(stderr, "  -C                    Dump ALL contigs (for -d option)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nreads min max       Dump tigs with between min and max reads (inclusive)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -d <operation>        Dump something about a multialign (-c or -u) in the store\n");
    fprintf(stderr, "     properties         ...properties\n");
    fprintf(stderr, "     frags              ...a list of fragments\n");
    fprintf(stderr, "     unitigs            ...a list of unitigs\n");
    fprintf(stderr, "     consensus [C]      ...the consensus sequence\n");
    fprintf(stderr, "                             if C supplied, only consensus with coverage >= C is output\n");
    fprintf(stderr, "     consensusgapped    ...the consensus sequence, with gaps as indicated in the multialignment\n");
    fprintf(stderr, "     layout             ...the layout\n");
    fprintf(stderr, "     multialign         ...the full multialignment\n");
    fprintf(stderr, "     matepair           ...an analysis of the mate pairs\n");
    fprintf(stderr, "     sizes              ...an analysis of sizes of the tigs\n");
    fprintf(stderr, "     coverage           ...an analysis of read coverage of the tigs\n");
    fprintf(stderr, "     overlap            ...an analysis of read overlaps in the tigs\n");
    fprintf(stderr, "     fmap               ...a map from fragment IID to unitig IID\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -E <editFile>         Change properties of multialigns\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -B <layout-file>      Construct a new store with unitigs in 'layout-file'.  Store versions\n");
    fprintf(stderr, "                        before that specified on the '-t' option are created but are empty.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -R <layout>           Replace a multialign with this one (type and id are from the layout)\n");
    fprintf(stderr, "                        The multialign is replaced in version <v> from -t.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -N                    Replace a multialign in the next version of the store.  This option is\n");
    fprintf(stderr, "                        needed if the version of the store to add a multialign does not exist.\n");
    fprintf(stderr, "                        The multialign is replaced in version <v>+1 from -t.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -compress             Move tigs from earlier versions into the specified version.  This removes\n");
    fprintf(stderr, "                        historical versions of unitigs/contigs, and can save tremendous storage space,\n");
    fprintf(stderr, "                        but makes it impossible to back up the assembly past the specified versions\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d multialign':\n");
    fprintf(stderr, "  -w width              Width of the page.\n");
    fprintf(stderr, "  -s spacing            Spacing between reads on the same line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d matepair' and '-d coverage':\n");
    fprintf(stderr, "  -o prefix             Output files will be written to 'prefix.*' in the current directory.\n");
    fprintf(stderr, "                        (defaults to 'tigStore' (the -t option) if not set.)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  For '-d sizes':\n");
    fprintf(stderr, "  -s genomesize         Denominator to use for n50 computation\n");
    exit(1);
  }

  //  To add a new multialign: In the layout, assign an id of -1 to the multialign (e.g., "unitig
  //  -1" or "contig -1").  Use -R to 'replace' this unitig in the store.  The store will assign the
  //  next unitig/contig ID to this new multialign.  WARNING!  The new multialign MUST be added to
  //  the latest version.
  //
  //  To delete a multialign: Remove ALL FRG and UTG lines, and set data.num_frags and
  //  data.num_unitigs to zero.  Use -R to 'replae' this unitig in the store.
  //  EXCEPT the code below will ignore treat these as EOF.
  //
  //  One can change partitioning by deleting a multialign from one partition and adding it to
  //  another partition.  Doing so WILL cause consensus to fail, as consensus is expecting a
  //  specific set of fragments in each partition.
  //
  //  It is not possible to add a new partition:
  //  MultiAlignStore::MultiAlignStore()-- ERROR, didn't find any unitigs or contigs in the store.  Correct version?


  if ((opType == OPERATION_BUILD) && (buildName != NULL)) {
    operationBuild(buildName, tigName, tigVers);
    exit(0);
  }


  if (opType == OPERATION_COMPRESS) {
    operationCompress(tigName, tigVers);
    exit(0);
  }


  gkpStore = new gkStore(gkpName, FALSE, FALSE);
  tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, FALSE, FALSE, FALSE);

  if (outPrefix == NULL)
    outPrefix = tigName;


  if ((opType == OPERATION_EDIT) && (editName != NULL)) {
    delete tigStore;
    tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, TRUE, FALSE);
    changeProperties(tigStore, editName);
  }


  if ((opType == OPERATION_REPLACE) && (replaceName != NULL)) {
    if (tigIDset) {
      fprintf(stderr, "ERROR:  -R is incompatible with -c, -u, -C and -U.  Did you mean -cp or -up instead?\n");
      exit(1);
    }

    errno = 0;
    FILE         *F = fopen(replaceName, "r");
    if (errno)
      fprintf(stderr, "Failed to open '%s': %s\n", replaceName, strerror(errno)), exit(1);

    fprintf(stderr, "Reading layouts from '%s'.\n", replaceName);

    delete tigStore;

    if (sameVersion)
      tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, true, false);  //  default
    else
      tigStore = new MultiAlignStore(tigName, tigVers, tigPartU, tigPartC, TRUE, false, false);

    MultiAlignT  *ma       = CreateEmptyMultiAlignT();
    bool          isUnitig = false;

    while (LoadMultiAlignFromHuman(ma, isUnitig, F) == true) {
      if (ma->data.num_frags + ma->data.num_unitigs == 0) {
        if (tigStore->isDeleted(ma->maID, isUnitig) == true) {
          fprintf(stderr, "DELETING %s %d -- ALREADY DELETED\n", (isUnitig) ? "unitig" : "contig", ma->maID);
        } else {
          fprintf(stderr, "DELETING %s %d\n", (isUnitig) ? "unitig" : "contig", ma->maID);
          tigStore->deleteMultiAlign(ma->maID, isUnitig);
        }
      } else {
        tigStore->insertMultiAlign(ma, isUnitig, FALSE);
        fprintf(stderr, "INSERTING %s %d\n", (isUnitig) ? "unitig" : "contig", ma->maID);
      }
    }

    fprintf(stderr, "Reading layouts from '%s' completed.\n", replaceName);

    fclose(F);
  }


  if (opType == OPERATION_UNITIGLIST) {
    tigStore->dumpMultiAlignRTable(true);
  }


  if (opType == OPERATION_CONTIGLIST) {
    tigStore->dumpMultiAlignRTable(false);
  }


  if (opType == OPERATION_PROPERTIES) {
    for (uint32 i=0; i<tigStore->numUnitigs(); i++) {
      if (tigStore->isDeleted(i, TRUE) == false) {
        fprintf(stdout, "unitig_coverage_stat %8u %d\n", i, tigStore->getUnitigCoverageStat(i));
        fprintf(stdout, "unitig_microhet_prob %8u %f\n", i, tigStore->getUnitigMicroHetProb(i));
        fprintf(stdout, "unitig_status        %8u %c\n", i, tigStore->getUnitigStatus(i));
        fprintf(stdout, "unitig_unique_rept   %8u %c\n", i, tigStore->getUnitigFUR(i));
      }
    }

    for (uint32 i=0; i<tigStore->numContigs(); i++) {
      if (tigStore->isDeleted(i, FALSE) == false) {
        fprintf(stdout, "contig_status        %8u %c\n", i, tigStore->getContigStatus(i));
      }
    }
  }


  if (opType == OPERATION_TIG) {
    if (tigIDset == false) {
      fprintf(stderr, "ERROR: No tig range set with -u, -c, -U or -C.\n");
      exit(1);
    }

    uint32   nTigs = (tigIsUnitig) ? tigStore->numUnitigs() : tigStore->numContigs();

    if (dumpAll == TRUE) {
      tigIDbgn = 0;
      tigIDend = nTigs - 1;
    }

    if (nTigs <= tigIDbgn) {
      fprintf(stderr, "ERROR: only "F_U32" %s in the store (IDs 0-"F_U32" inclusive); can't dump requested range "F_U32"-"F_U32"\n",
              nTigs,
              (tigIsUnitig) ? "unitigs" : "contigs",
              nTigs-1,
              tigIDbgn, tigIDend);
    }

    if (nTigs <= tigIDend)
      tigIDend = nTigs - 1;

    if (dumpFlags == DUMP_MATEPAIR)
      mpa = new matePairAnalysis(gkpName);

    if (dumpFlags == DUMP_SIZES)
      siz = new sizeAnalysis(sizSize);

    if (dumpFlags == DUMP_COVERAGE) {
      covMax = 1048576;
      cov    = new uint64 [covMax];

      memset(cov, 0, sizeof(uint64) * covMax);
    }

    for (uint32 ti=tigIDbgn; ti<=tigIDend; ti++) {
      uint32  Nreads = tigStore->getNumFrags(ti, tigIsUnitig);
  
      if ((Nreads    < minNreads) ||
          (maxNreads < Nreads))
        continue;

      ma = tigStore->loadMultiAlign(ti, tigIsUnitig);

      if (ma == NULL)
        continue;

      if (dumpFlags == DUMP_PROPERTIES)
        dumpProperties(tigStore, ti, tigIsUnitig, ma);

      if (dumpFlags == DUMP_FRAGS)
        dumpFrags(tigStore, ti, tigIsUnitig, ma);

      if (dumpFlags == DUMP_UNITIGS)
        dumpUnitigs(tigStore, ti, tigIsUnitig, ma);

      if (dumpFlags == DUMP_CONSENSUS)
        dumpConsensus(tigStore, ti, tigIsUnitig, ma, false, minCoverage);

      if (dumpFlags == DUMP_CONSENSUSGAPPED)
        dumpConsensus(tigStore, ti, tigIsUnitig, ma, true, minCoverage);

      if (dumpFlags == DUMP_LAYOUT)
        DumpMultiAlignForHuman(stdout, ma, tigIsUnitig);

      if (dumpFlags == DUMP_MULTIALIGN)
        PrintMultiAlignT(stdout, ma, gkpStore, showQV, showDots, (tigIsUnitig) ? AS_READ_CLEAR_OBTCHIMERA : AS_READ_CLEAR_LATEST);

      if (dumpFlags == DUMP_MATEPAIR)
        mpa->evaluateTig(ma);

      if (dumpFlags == DUMP_SIZES)
        siz->evaluateTig(ma, tigIsUnitig);

      if (dumpFlags == DUMP_COVERAGE)
        dumpCoverage(tigStore, ti, tigIsUnitig, ma, 2, UINT32_MAX, cov, covMax, outPrefix);

      if (dumpFlags == DUMP_THINOVERLAP)
        dumpThinOverlap(tigStore, ti, tigIsUnitig, ma, sizSize);

      if (dumpFlags == DUMP_FMAP)
        dumpFmap(stdout, ma, tigIsUnitig);

      tigStore->unloadMultiAlign(ti, tigIsUnitig);
    }
  }

  if (mpa) {
    mpa->finalize();
    mpa->printSummary(stdout);
    mpa->writeUpdate(outPrefix);
    mpa->drawPlots(outPrefix);
    delete mpa;
  }

  if (siz) {
    siz->finalize();
    siz->printSummary(stdout);
    delete siz;
  }

  if (cov) {
    char  N[FILENAME_MAX];
    FILE *F;

    sprintf(N, "%s.depthHistogram", outPrefix);

    errno = 0;
    F = fopen(N, "w");

    uint32   hMax = covMax - 1;
    while ((hMax > 0) && (cov[hMax] == 0))
      hMax--;

    for (uint32 i=0; i<=hMax; i++)
      fprintf(F, F_U32"\t"F_U64"\n", i, cov[i]);

    fclose(F);

    delete [] cov;

    FILE *gnuPlot = popen("gnuplot > /dev/null 2>&1", "w");

    if (gnuPlot) {
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "set output '%s.depthHistogram.png'\n", outPrefix);
      fprintf(gnuPlot, "set xlabel 'depth'\n");
      fprintf(gnuPlot, "set ylabel 'number of bases'\n");
      fprintf(gnuPlot, "set terminal 'png'\n");
      fprintf(gnuPlot, "plot '%s.depthHistogram' using 1:2 with lines title '%s %s %u-%u depthHistogram', \\\n",
              outPrefix, outPrefix, (tigIsUnitig) ? "unitigs" : "contigs", tigIDbgn, tigIDend);

      fclose(gnuPlot);
    }
  }

  delete gkpStore;
  delete tigStore;

  exit(0);
}
