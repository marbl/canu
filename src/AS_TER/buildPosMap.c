
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

const char *mainid = "$Id: buildPosMap.c,v 1.23 2012-11-27 18:33:15 brianwalenz Exp $";

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "AS_UTL_splitToWords.H"

#include <map>

using namespace std;

map<AS_UID,AS_IID>  uid2iid;
map<AS_UID,AS_IID>  uid2lib;

map<AS_UID,double>  covStat;
map<AS_UID,double>  microHet;

#define ORIR 'r'
#define ORIF 'f'


typedef struct {
  int      len;
  AS_UID   scfUID;
  int      scfBgn;
  int      scfEnd;
  char     scfOri;
} ctgInfo_t;

uint32       ctgInfoMax = 32 * 1024 * 1024;
ctgInfo_t   *ctgInfo    = NULL;

AS_IID       lastAFG = 1;

FILE *frags     = NULL;
FILE *mates     = NULL;
FILE *libraries = NULL;

FILE *utginf    = NULL;
FILE *deginf    = NULL;
FILE *ctginf    = NULL;
FILE *scfinf    = NULL;

FILE *utglkg    = NULL;
FILE *ctglkg    = NULL;
FILE *scflkg    = NULL;

FILE *frgutg    = NULL;
FILE *sfgutg    = NULL;

FILE *frgdeg    = NULL;
FILE *utgdeg    = NULL;
FILE *vardeg    = NULL;

FILE *frgctg    = NULL;
FILE *sfgctg    = NULL;
FILE *utgctg    = NULL;
FILE *varctg    = NULL;

FILE *frgscf    = NULL;
FILE *sfgscf    = NULL;
FILE *utgscf    = NULL;
FILE *ctgscf    = NULL;
FILE *varscf    = NULL;


FILE *utglen    = NULL;
FILE *deglen    = NULL;
FILE *ctglen    = NULL;
FILE *scflen    = NULL;


FILE *
openFile(char *label, char *prefix, int write) {
  char  N[FILENAME_MAX];
  FILE *F = NULL;

  sprintf(N, "%s.posmap.%s", prefix, label);

  errno = 0;
  F = fopen(N, (write) ? "w" : "r");
  if (errno) {
    fprintf(stderr, "Failed to open '%s' for %s: %s\n", N, (write) ? "write" : "read", strerror(errno));
    exit(1);
  }

  return(F);
}



void
processMDI(SnapMateDistMesg *mdi, char *prefix) {
  int32  samples = 0;

  for (int32 i=0; i<mdi->num_buckets; i++)
    samples += mdi->histogram[i];

  fprintf(libraries, "%s\t%.0f\t%.0f\t%d\t%d\t%d\n",
          AS_UID_toString(mdi->erefines),
          mdi->mean, mdi->stddev,
          mdi->min, mdi->max, samples);

  double  bgn   = mdi->min;
  double  width = (double)(mdi->max - mdi->min) / mdi->num_buckets;

  char  N[FILENAME_MAX];

  sprintf(N, "%s.posmap.libraries.%s", prefix, AS_UID_toString(mdi->erefines));

  errno = 0;
  FILE *F = fopen(N, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", N, strerror(errno)), exit(1);

  for (int32 i=0; i<mdi->num_buckets; i++, bgn += width)
    fprintf(F, "%.0f\t%.1f\t%.0f\t%d\n", bgn, bgn + width/2, bgn + width, mdi->histogram[i]);

  fclose(F);
}



char const *mateStatusStrings[15] = { 
  "ERROR",
  "unassigned",
  "good",
  "badShort",
  "badLong",
  "badSame",
  "badOuttie",
  "notMated",
  "bothChaff",
  "oneChaff",
  "bothDegen",
  "oneDegen",
  "bothSurrogate",
  "oneSurrogate",
  "diffScaffold"
};

char const *unitigStatusStrings[6] = {
  "ERROR",
  "unique",
  "surrogate",
  "degenerate",
  "chimer",
  "unresolved"
};

char const *unitigTypeStrings[7] = {
  "ERROR",
  "unique",
  "rock",
  "stone",
  "pebble",
  "singleton",
  "other"
};


char const *
decodeMateStatus(char mateStatusCode) {
  int32   mss = 0;

  switch (mateStatusCode) {
    case 'Z':  mss =  1;  break;
    case 'G':  mss =  2;  break;
    case 'C':  mss =  3;  break;
    case 'L':  mss =  4;  break;
    case 'S':  mss =  5;  break;
    case 'O':  mss =  6;  break;
    case 'N':  mss =  7;  break;
    case 'H':  mss =  8;  break;
    case 'A':  mss =  9;  break;
    case 'D':  mss = 10;  break;
    case 'E':  mss = 11;  break;
    case 'U':  mss = 12;  break;
    case 'R':  mss = 13;  break;
    case 'F':  mss = 14;  break;
    default:
      fprintf(stderr, "WARNING:  Invalid mate status %c\n", mateStatusCode);
      break;
  }

  return(mateStatusStrings[mss]);
}

char const *
decodeUnitigStatus(char unitigStatusCode) {
  int32  uss = 0;

  switch (unitigStatusCode) {
    case 'U':  uss = 1;  break;
    case 'S':  uss = 2;  break;
    case 'N':  uss = 3;  break;
    case 'C':  uss = 4;  break;
    case 'X':  uss = 5;  break;
    default:
      fprintf(stderr, "WARNING:  Invalid unitig status %c\n", unitigStatusCode);
      break;
  }

  return(unitigStatusStrings[uss]);
}

char const *
decodeUnitigType(char unitigTypeCode) {
  int32  uss = 0;

  switch (unitigTypeCode) {
    case 'U':  uss = 1;  break;
    case 'R':  uss = 2;  break;
    case 'S':  uss = 3;  break;
    case 'P':  uss = 4;  break;
    case 's':  uss = 5;  break;
    case 'X':  uss = 6;  break;
    default:
      fprintf(stderr, "WARNING:  Invalid unitig status %c\n", unitigTypeCode);
      break;
  }

  return(unitigTypeStrings[uss]);
}





void
processAFG(AugFragMesg *afg, gkStore *gkp) {
  char *libraryName = "N/A";

  if (gkp) {
    gkFragment  frg;

    while (lastAFG < afg->iaccession) {
      gkp->gkStore_getFragment(lastAFG, &frg, GKFRAGMENT_INF);

      fprintf(frags, "%s\t%d\t%d\t%s\t%s\t%s\n",
              AS_UID_toString(frg.gkFragment_getReadUID()),
              frg.gkFragment_getClearRegionBegin(),
              frg.gkFragment_getClearRegionEnd(),
              "deleted",
              "notMated",
              gkp->gkStore_getLibrary(frg.gkFragment_getLibraryIID())->libraryName);

      lastAFG++;
    }

    assert(lastAFG == afg->iaccession);

    gkp->gkStore_getFragment(lastAFG, &frg, GKFRAGMENT_INF);

    if (frg.gkFragment_getMateIID() > 0)
      uid2lib[afg->eaccession] = frg.gkFragment_getLibraryIID();

    libraryName = gkp->gkStore_getLibrary(frg.gkFragment_getLibraryIID())->libraryName;
  }

  fprintf(frags, "%s\t%d\t%d\t%s\t%s\t%s\n",
          AS_UID_toString(afg->eaccession),
          afg->clear_rng.bgn,
          afg->clear_rng.end,
          afg->chaff ? "chaff" : "placed",
          decodeMateStatus(afg->mate_status),
          libraryName);

  lastAFG++;
}


void
processAMP(AugMatePairMesg *amp, gkStore *gkp) {
  char *libraryName = "N/A";

  if (gkp) {
    AS_IID   lib1  = uid2lib[amp->fragment1];
    AS_IID   lib2  = uid2lib[amp->fragment2];

    assert(lib1 != 0);
    assert(lib2 != 0);
    assert(lib1 == lib2);

    libraryName = gkp->gkStore_getLibrary(lib1)->libraryName;
  }

  fprintf(mates, "%s\t%s\t%s\t%s\n",
          AS_UID_toString(amp->fragment1),
          AS_UID_toString(amp->fragment2),
          decodeMateStatus(amp->mate_status),
          libraryName);
}


void
processUTG(SnapUnitigMesg *utg) {
  int   *unitigGapToUngap   = NULL;
  int    unitigLengthGapped = strlen(utg->consensus);
  int    unitigLength       = 0;
  int    i = 0, j = 0;

  unitigGapToUngap = (int *)safe_calloc(unitigLengthGapped + 1, sizeof(int));

  unitigGapToUngap[0] = 0;

  for (i=0; i<unitigLengthGapped; i++) {
    if (utg->consensus[i] != '-')
      unitigLength++;
    unitigGapToUngap[i+1] = unitigLength;
  }

  fprintf(utglen, "%s\t%d\n",
          AS_UID_toString(utg->eaccession),
          unitigLength);
  fprintf(utginf, "%s\t%d\t%d\t%f\t%f\t%s\n",
          AS_UID_toString(utg->eaccession),
          unitigLength,
          utg->num_frags,
          utg->coverage_stat,
          utg->microhet_prob,
          decodeUnitigStatus(utg->status));

  covStat[utg->eaccession]  = utg->coverage_stat;
  microHet[utg->eaccession] = utg->microhet_prob;

  //  MPS/fragments
  for (i=0; i<utg->num_frags; i++) {
    int  bgn = unitigGapToUngap[utg->f_list[i].position.bgn];
    int  end = unitigGapToUngap[utg->f_list[i].position.end];
    char ori = ORIF;

    //  If bgn == end, then the fragment fell entirely within a gap.

    if (utg->f_list[i].position.bgn > utg->f_list[i].position.end) {
      ori = ORIR;
      end = unitigGapToUngap[utg->f_list[i].position.bgn];
      bgn = unitigGapToUngap[utg->f_list[i].position.end];
    }

    fprintf(frgutg, "%s\t%s\t%d\t%d\t%c\n",
            AS_UID_toString(utg->f_list[i].eident),
            AS_UID_toString(utg->eaccession),
            bgn, end, ori);
  }

  //  Remember what fragments are in surrogates
#warning not remembering surrogate fragments
  if (utg->status == AS_SEP) {
    for (i=0; i<utg->num_frags; i++) {
    }
  }

  safe_free(unitigGapToUngap);
}




void
processULK(SnapUnitigLinkMesg *ulk) {

  fprintf(utglkg, "%s\t%s\t%c\t%c\t%s\t%f\t%f\t%d\t%c",
          AS_UID_toString(ulk->eunitig1),
          AS_UID_toString(ulk->eunitig2),
          ulk->orientation.toLetter(),
          ulk->overlap_type,
          ulk->is_possible_chimera ? "chimeric" : ".",
          ulk->mean_distance,
          ulk->std_deviation,
          ulk->num_contributing,
          ulk->status);
  //  If overlap type indicates an overlap was counted, subtract one
  for (int32 i=0; i < ulk->num_contributing - (ulk->overlap_type != AS_NO_OVERLAP); i++)
    fprintf(utglkg, "\t%s,%s,%c",
            AS_UID_toString(ulk->jump_list[i].in1),
            AS_UID_toString(ulk->jump_list[i].in2),
            ulk->jump_list[i].type.toLetter());
  fprintf(utglkg, "\n");
}



void
processCCO(SnapConConMesg *cco) {
  int   *contigGapToUngap   = NULL;
  int    contigLengthGapped = strlen(cco->consensus);
  int    contigLengthUngap  = 0;
  int    i = 0, j = 0;
  int    isDegenerate = 0;

  //  By definition, a degenerate contig has one unitig and is unplaced.
  //  In reality, those two conditions always occur together.
  //
  FILE  *len = ctglen;
  FILE  *frg = frgctg;
  FILE  *utg = utgctg;
  FILE  *var = varctg;

  if ((cco->placed == AS_UNPLACED) && (cco->num_unitigs == 1)) {
    isDegenerate = 1;
    len = deglen;
    frg = frgdeg;
    utg = utgdeg;
    var = vardeg;
  }


  contigGapToUngap = (int *)safe_calloc(contigLengthGapped + 1, sizeof(int));

  contigGapToUngap[0] = 0;

  for (i=0; i<contigLengthGapped; i++) {
    if (cco->consensus[i] != '-')
      contigLengthUngap++;
    contigGapToUngap[i+1] = contigLengthUngap;
  }

  fprintf(len, "%s\t%d\n",
          AS_UID_toString(cco->eaccession),
          contigLengthUngap);

  if (isDegenerate == 0)
    fprintf(ctginf, "%s\t%d\t%d\t%d\t%d\n",
            AS_UID_toString(cco->eaccession),
            contigLengthUngap,
            cco->num_unitigs,
            cco->num_pieces,
            cco->num_vars);
  else
    fprintf(deginf, "%s\t%d\t%d\t%d\t%d\t%s\t%f\t%f\n",
            AS_UID_toString(cco->eaccession),
            contigLengthUngap,
            cco->num_unitigs,
            cco->num_pieces,
            cco->num_vars,
            AS_UID_toString(cco->unitigs[0].eident),
            covStat[cco->unitigs[0].eident],
            microHet[cco->unitigs[0].eident]);

  uid2iid[cco->eaccession] = cco->iaccession;

  if (ctgInfoMax <= cco->iaccession) {
    ctgInfoMax *= 2;
    ctgInfo     = (ctgInfo_t *)safe_realloc(ctgInfo, ctgInfoMax * sizeof(ctgInfo_t));
  }
  ctgInfo[cco->iaccession].len = contigLengthUngap;

  //  VAR/variants
  for (i=0; i<cco->num_vars; i++) {
    IntMultiVar *v = cco->vars + i;

    int  bgn = contigGapToUngap[v->position.bgn];
    int  end = contigGapToUngap[v->position.end];

    fprintf(var, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n",
            v->enc_var_seq,
            AS_UID_toString(cco->eaccession),
            bgn, end,
            v->num_reads,
            v->num_alleles_confirmed,
            v->min_anchor_size,
            v->var_length,
            v->enc_num_reads,
            v->enc_weights,
            v->enc_read_ids);
  }

  //  MPS/fragments
  for (i=0; i<cco->num_pieces; i++) {
    int  bgn = contigGapToUngap[cco->pieces[i].position.bgn];
    int  end = contigGapToUngap[cco->pieces[i].position.end];
    char ori = ORIF;

    //  If bgn == end, then the fragment fell entirely within a gap.

    if (cco->pieces[i].position.bgn > cco->pieces[i].position.end) {
      ori = ORIR;
      end = contigGapToUngap[cco->pieces[i].position.bgn];
      bgn = contigGapToUngap[cco->pieces[i].position.end];
    }

    fprintf(frg, "%s\t%s\t%d\t%d\t%c\n",
            AS_UID_toString(cco->pieces[i].eident),
            AS_UID_toString(cco->eaccession),
            bgn, end, ori);
  }

  //  UPS/unitigs
  for (i=0; i<cco->num_unitigs; i++) {
    int  bgn = contigGapToUngap[cco->unitigs[i].position.bgn];
    int  end = contigGapToUngap[cco->unitigs[i].position.end];
    char ori = ORIF;

    //  If this is a SURROGATE UNITIG, report "fraguid unitiguid contiguid"

    if (cco->unitigs[i].position.bgn > cco->unitigs[i].position.end) {
      ori = ORIR;
      end = contigGapToUngap[cco->unitigs[i].position.bgn];
      bgn = contigGapToUngap[cco->unitigs[i].position.end];
    }

    if (isDegenerate == 0)
      fprintf(utg, "%s\t%s\t%d\t%d\t%c\t%s\n",
              AS_UID_toString(cco->unitigs[i].eident),
              AS_UID_toString(cco->eaccession),
              bgn, end, ori,
              decodeUnitigType(cco->unitigs[i].type));
    else
      fprintf(utg, "%s\t%s\t%d\t%d\t%c\n",
              AS_UID_toString(cco->unitigs[i].eident),
              AS_UID_toString(cco->eaccession),
              bgn, end, ori);
  }

  safe_free(contigGapToUngap);
}



void
processCLK(SnapContigLinkMesg *clk) {

  fprintf(ctglkg, "%s\t%s\t%c\t%c\t%s\t%f\t%f\t%d\t%c",
          AS_UID_toString(clk->econtig1),
          AS_UID_toString(clk->econtig2),
          clk->orientation.toLetter(),
          clk->overlap_type,
          clk->is_possible_chimera ? "chimeric" : ".",
          clk->mean_distance,
          clk->std_deviation,
          clk->num_contributing,
          clk->status);
  //  If overlap type indicates an overlap was counted, subtract one
  for (int32 i=0; i < clk->num_contributing - (clk->overlap_type != AS_NO_OVERLAP); i++)
    fprintf(ctglkg, "\t%s,%s,%c",
            AS_UID_toString(clk->jump_list[i].in1),
            AS_UID_toString(clk->jump_list[i].in2),
            clk->jump_list[i].type.toLetter());
  fprintf(ctglkg, "\n");
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
  int     i = 0;
  int     bgn = 0, end = 0;
  char    ori = 0;
  int     scfLen = 0;
  AS_IID  ctgIID = 0;

  //  Tricky.  If scf->num_contig_pairs==0, there is a single contig
  //  in this scaffold.

  int  singleContig = 0;
  if (scf->num_contig_pairs == 0) {
    singleContig = 1;
    scf->num_contig_pairs = 1;
  }

  //  CTP/contig pairs
  for (i=0; i<scf->num_contig_pairs; i++) {

    //  First contig-pair, print the first contig.
    //
    if (i == 0) {
      ctgIID = uid2iid[scf->contig_pairs[i].econtig1];

      bgn     = scfLen;
      scfLen += ctgInfo[ctgIID].len;
      end     = scfLen;
      ori     = ORIF;

      if ((scf->contig_pairs[i].orient.isAnti()) ||
          (scf->contig_pairs[i].orient.isOuttie()))
        ori = ORIR;

      ctgInfo[ctgIID].scfUID = scf->eaccession;
      ctgInfo[ctgIID].scfBgn = bgn;
      ctgInfo[ctgIID].scfEnd = end;
      ctgInfo[ctgIID].scfOri = ori;

      fprintf(ctgscf, "%s\t%s\t%d\t%d\t%c\n",
              AS_UID_toString(scf->contig_pairs[i].econtig1),
              AS_UID_toString(scf->eaccession),
              bgn, end, ori);
    }

    //  Not first contig-pair, or there is more than one pair, print
    //  the seocnd contig.
    //
    if ((i > 0) || (singleContig == 0)) {
      ctgIID = uid2iid[scf->contig_pairs[i].econtig2];

      scfLen += computeGapSize(scf->contig_pairs[i].mean);
      bgn     = scfLen;
      scfLen += ctgInfo[ctgIID].len;
      end     = scfLen;
      ori     = ORIF;

      if ((scf->contig_pairs[i].orient.isAnti()) ||
          (scf->contig_pairs[i].orient.isInnie()))
        ori = ORIR;

      ctgInfo[ctgIID].scfUID = scf->eaccession;
      ctgInfo[ctgIID].scfBgn = bgn;
      ctgInfo[ctgIID].scfEnd = end;
      ctgInfo[ctgIID].scfOri = ori;

      fprintf(ctgscf, "%s\t%s\t%d\t%d\t%c\n",
              AS_UID_toString(scf->contig_pairs[i].econtig2),
              AS_UID_toString(scf->eaccession),
              bgn, end, ori);
    }
  }

  fprintf(scflen, "%s\t%d\n", AS_UID_toString(scf->eaccession), scfLen);
}



void
processSLK(SnapScaffoldLinkMesg *slk) {

  fprintf(scflkg, "%s\t%s\t%c\t%f\t%f\t%d",
          AS_UID_toString(slk->escaffold1),
          AS_UID_toString(slk->escaffold2),
          slk->orientation.toLetter(),
          slk->mean_distance,
          slk->std_deviation,
          slk->num_contributing);
  //  Unlike ULK and CLK, there is no overlap_type here, so all num_contributing are edges
  for (int32 i=0; i<slk->num_contributing; i++)
    fprintf(scflkg, "\t%s,%s,%c",
            AS_UID_toString(slk->jump_list[i].in1),
            AS_UID_toString(slk->jump_list[i].in2),
            slk->jump_list[i].type.toLetter());
  fprintf(scflkg, "\n");
}










void
transferFrg(FILE *ctg, FILE *scf) {
  char           line[1024]   = {0};
  splitToWords   W;
  int32          bgn = 0;
  int32          end = 0;

  rewind(ctg);

  fgets(line, 1024, ctg);

  while (!feof(ctg)) {
    chomp(line);
    W.split(line);

    //  frgUID
    //  ctgUID
    //  bgn
    //  end
    //  ori
    //  type (optional)

    AS_UID uid = AS_UID_lookup(W[1], NULL);
    AS_IID iid = uid2iid[uid];

    if (iid > 0) {
      uid  = ctgInfo[iid].scfUID;

      assert(AS_UID_isDefined(uid));

      if (ctgInfo[iid].scfOri == ORIF) {
        bgn = ctgInfo[iid].scfBgn + W(2);
        end = ctgInfo[iid].scfBgn + W(3);
      } else {
        end = ctgInfo[iid].scfEnd - W(2);
        bgn = ctgInfo[iid].scfEnd - W(3);
      }

      char ori = W[4][0];

      if (ori == ctgInfo[iid].scfOri)
        ori = ORIF;
      else
        ori = ORIR;

      if (W[5] == NULL)
        fprintf(scf, "%s\t%s\t%d\t%d\t%c\n",
                W[0], AS_UID_toString(uid), bgn, end, ori);
      else
        fprintf(scf, "%s\t%s\t%d\t%d\t%c\t%s\n",
                W[0], AS_UID_toString(uid), bgn, end, ori, W[5]);
    }

    fgets(line, 1024, ctg);
  }
}




#warning this function will overflow on large (bogus) var sequences

void
transferVar(FILE *ctg, FILE *scf) {
  int32  lineMax = 1 * 1024 * 1024;
  char  *line    = new char [lineMax];

  memset(line, 0, sizeof(char) * lineMax);

  rewind(ctg);

  fgets(line, lineMax, ctg);
  chomp(line);

  while (!feof(ctg)) {
    splitToWords  W(line);

    assert(line[lineMax-1] == 0);

    if (11 == W.numWords()) {
      int32 bgn = atoi(W[2]);
      int32 end = atoi(W[3]);

      AS_UID uid = AS_UID_lookup(W[1], NULL);
      AS_IID iid = uid2iid[uid];

      if (iid > 0) {
        uid  = ctgInfo[iid].scfUID;

        assert(AS_UID_isDefined(uid));

        if (ctgInfo[iid].scfOri == ORIF) {
          bgn = ctgInfo[iid].scfBgn + bgn;
          end = ctgInfo[iid].scfBgn + end;
        } else {
          int tmp;
          tmp = ctgInfo[iid].scfEnd - bgn;
          bgn = ctgInfo[iid].scfEnd - end;
          end = tmp;
        }

        fprintf(scf, "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                W[0], AS_UID_toString(uid), bgn, end, W[4], W[5], W[6], W[7], W[8], W[9], W[10]);
      }
    }

    fgets(line, lineMax, ctg);
    chomp(line);
  }

  delete [] line;
}




int main (int argc, char *argv[]) {
  char    *outputPrefix       = NULL;
  char    *asmName            = NULL;
  FILE    *asmFile            = stdin;
  char    *gkpName            = NULL;
  gkStore *gkp                = NULL;

  GenericMesg *pmesg       = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-i") == 0) {
      asmName = argv[++arg];

    } else if (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((outputPrefix == NULL) || (err)) {
    fprintf(stderr, "usage: %s -o prefix [-h] [-i prefix.asm | < prefix.asm]\n", argv[0]);
    fprintf(stderr, "  -o prefix        write the output here\n");
    fprintf(stderr, "  -i prefix.asm    read the assembly from here; default is to read stdin\n");
    fprintf(stderr, "  -h               print help\n");

    fprintf(stderr, "\n");
    exit(1);
  }

  if (asmName) {
    errno = 0;
    asmFile = fopen(asmName, "r");
    if (errno)
      fprintf(stderr, "Failed to open assembly file '%s': %s\n", asmName, strerror(errno)), exit(1);
  }

  if (gkpName) {
    gkp = new gkStore(gkpName, FALSE, FALSE);
    gkp->gkStore_metadataCaching(true);
  }

  ctgInfo = (ctgInfo_t *)safe_calloc(ctgInfoMax, sizeof(ctgInfo_t));

  frags     = openFile("frags",     outputPrefix, 1);
  mates     = openFile("mates",     outputPrefix, 1);
  libraries = openFile("libraries", outputPrefix, 1);

  utginf    = openFile("utginf",    outputPrefix, 1);
  deginf    = openFile("deginf",    outputPrefix, 1);
  ctginf    = openFile("ctginf",    outputPrefix, 1);
  scfinf    = openFile("scfinf",    outputPrefix, 1);

  utglkg    = openFile("utglkg",    outputPrefix, 1);
  ctglkg    = openFile("ctglkg",    outputPrefix, 1);
  scflkg    = openFile("scflkg",    outputPrefix, 1);

  frgutg    = openFile("frgutg",    outputPrefix, 1);
  sfgutg    = openFile("sfgutg",    outputPrefix, 1);

  frgdeg    = openFile("frgdeg",    outputPrefix, 1);
  utgdeg    = openFile("utgdeg",    outputPrefix, 1);
  vardeg    = openFile("vardeg",    outputPrefix, 1);

  frgctg    = openFile("frgctg",    outputPrefix, 1);
  sfgctg    = openFile("sfgctg",    outputPrefix, 1);
  utgctg    = openFile("utgctg",    outputPrefix, 1);
  varctg    = openFile("varctg",    outputPrefix, 1);

  frgscf    = openFile("frgscf",    outputPrefix, 1);
  sfgscf    = openFile("sfgscf",    outputPrefix, 1);
  utgscf    = openFile("utgscf",    outputPrefix, 1);
  ctgscf    = openFile("ctgscf",    outputPrefix, 1);
  varscf    = openFile("varscf",    outputPrefix, 1);

  utglen    = openFile("utglen",    outputPrefix, 1);
  deglen    = openFile("deglen",    outputPrefix, 1);
  ctglen    = openFile("ctglen",    outputPrefix, 1);
  scflen    = openFile("scflen",    outputPrefix, 1);


  while(ReadProtoMesg_AS(asmFile, &pmesg) != EOF){
    switch(pmesg->t){
      case MESG_MDI:
        processMDI((SnapMateDistMesg *)pmesg->m, outputPrefix);
        break;

      case MESG_AFG:
        processAFG((AugFragMesg *)pmesg->m, gkp);
        break;
      case MESG_AMP:
        processAMP((AugMatePairMesg *)pmesg->m, gkp);
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

  delete gkp;

  fclose(frags);
  fclose(mates);

  fclose(utginf);
  fclose(deginf);
  fclose(ctginf);
  fclose(scfinf);

  fclose(frgutg);
  fclose(sfgutg);

  fclose(frgdeg);
  fclose(utgdeg);
  fclose(vardeg);

  fclose(frgctg);
  fclose(sfgctg);
  fclose(utgctg);
  fclose(varctg);

  fclose(utglen);
  fclose(ctglen);
  fclose(deglen);
  fclose(scflen);

  fclose(ctgscf);

  //  Transfer the contig based info onto scaffolds.

  frgctg = openFile("frgctg", outputPrefix, 0);
  utgctg = openFile("utgctg", outputPrefix, 0);
  varctg = openFile("varctg", outputPrefix, 0);

  transferFrg(frgctg, frgscf);
  transferFrg(utgctg, utgscf);
  transferVar(varctg, varscf);

  fclose(frgctg);
  fclose(utgctg);
  fclose(varctg);

  fclose(frgscf);
  fclose(utgscf);
  fclose(varscf);

  exit(0);
}
