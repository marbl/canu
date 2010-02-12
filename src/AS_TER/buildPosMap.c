
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

const char *mainid = "$Id: buildPosMap.c,v 1.13 2010-02-12 20:33:12 brianwalenz Exp $";

#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <unistd.h>
#include  <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

HashTable_AS   *uid2iid      = NULL;

#define ORIR 'r'
#define ORIF 'f'


typedef struct {
  int      len;
  AS_UID   scfUID;
  int      scfBgn;
  int      scfEnd;
  char     scfOri;
} ctgInfo_t;

int          ctgInfoMax = 32 * 1024 * 1024;
int          ctgInfoLen = 0;
ctgInfo_t   *ctgInfo    = NULL;


FILE *frags  = NULL;
FILE *mates  = NULL;

FILE *utglkg = NULL;
FILE *ctglkg = NULL;
FILE *scflkg = NULL;

FILE *frgutg = NULL;

FILE *frgdeg = NULL;
FILE *utgdeg = NULL;
FILE *vardeg = NULL;

FILE *frgctg = NULL;
FILE *utgctg = NULL;
FILE *varctg = NULL;

FILE *frgscf = NULL;
FILE *utgscf = NULL;
FILE *ctgscf = NULL;
FILE *varscf = NULL;

FILE *utglen = NULL;
FILE *deglen = NULL;
FILE *ctglen = NULL;
FILE *scflen = NULL;


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
processMDI(SnapMateDistMesg *mdi) {
  int  samples = 0;
  int  i;

  for (i=0; i<mdi->num_buckets; i++)
    samples += mdi->histogram[i];

  fprintf(stderr, "MDI\t%s\t%.6f\t%.6f\t%d\n",
          AS_UID_toString(mdi->erefines),
          mdi->mean, mdi->stddev, samples);
}


char *
decodeMateStatus(char mateStatusCode, char *mateStatus) {

  switch (mateStatusCode) {
    case 'Z':
      strcpy(mateStatus, "unassigned");
      break;
    case 'G':
      strcpy(mateStatus, "good");
      break;
    case 'C':
      strcpy(mateStatus, "badShort");
      break;
    case 'L':
      strcpy(mateStatus, "badLong");
      break;
    case 'S':
      strcpy(mateStatus, "badSame");
      break;
    case 'O':
      strcpy(mateStatus, "badOuttie");
      break;
    case 'N':
      strcpy(mateStatus, "notMated");
      break;
    case 'H':
      strcpy(mateStatus, "bothChaff");
      break;
    case 'A':
      strcpy(mateStatus, "oneChaff");
      break;
    case 'D':
      strcpy(mateStatus, "bothDegen");
      break;
    case 'E':
      strcpy(mateStatus, "oneDegen");
      break;
    case 'U':
      strcpy(mateStatus, "bothSurrogate");
      break;
    case 'R':
      strcpy(mateStatus, "oneSurrogate");
      break;
    case 'F':
      strcpy(mateStatus, "diffScaffold");
      break;
    default:
      strcpy(mateStatus, "ERROR");
      break;
  }

  return(mateStatus);
}


void
processAFG(AugFragMesg *afg) {
  char   mateStatus[256] = {0};

  fprintf(frags, "%s\t%d\t%d\t%s\t%s\n",
          AS_UID_toString(afg->eaccession),
          afg->clear_rng.bgn,
          afg->clear_rng.end,
          afg->chaff    ? "chaff"    : "placed",
          decodeMateStatus(afg->mate_status, mateStatus));
}


void
processAMP(AugMatePairMesg *amp) {
  char   mateStatus[256] = {0};

  fprintf(mates, "%s\t%s\t%s\n",
          AS_UID_toString(amp->fragment1),
          AS_UID_toString(amp->fragment2),
          decodeMateStatus(amp->mate_status, mateStatus));
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

  fprintf(utglen, "%s\t%d\n", AS_UID_toString(utg->eaccession), unitigLength);

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

  fprintf(utglkg, "%s\t%s\t%c\t%c\t%c\t%f\t%f\t%d\t%c\n",
          AS_UID_toString(ulk->eunitig1),
          AS_UID_toString(ulk->eunitig2),
          ulk->orientation,
          ulk->overlap_type,
          ulk->is_possible_chimera ? 'chimeric' : '.',
          ulk->mean_distance,
          ulk->std_deviation,
          ulk->num_contributing,
          ulk->status);
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
    //  Degenerate, use deg files.
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

  fprintf(len, "%s\t%d\n", AS_UID_toString(cco->eaccession), contigLengthUngap);

  InsertInHashTable_AS(uid2iid,
                       AS_UID_toInteger(cco->eaccession), 0,
                       cco->iaccession, 0);
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

    fprintf(utg, "%s\t%s\t%d\t%d\t%c\n",
            AS_UID_toString(cco->unitigs[i].eident),
            AS_UID_toString(cco->eaccession),
            bgn, end, ori);
  }

  safe_free(contigGapToUngap);
}



void
processCLK(SnapContigLinkMesg *clk) {

  fprintf(ctglkg, "%s\t%s\t%c\t%c\t%c\t%f\t%f\t%d\t%c\n",
          AS_UID_toString(clk->econtig1),
          AS_UID_toString(clk->econtig2),
          clk->orientation,
          clk->overlap_type,
          clk->is_possible_chimera ? 'chimeric' : '.',
          clk->mean_distance,
          clk->std_deviation,
          clk->num_contributing,
          clk->status);
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
      ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                                 AS_UID_toInteger(scf->contig_pairs[i].econtig1),
                                                 0);

      bgn     = scfLen;
      scfLen += ctgInfo[ctgIID].len;
      end     = scfLen;
      ori     = ORIF;

      if ((scf->contig_pairs[i].orient == AS_ANTI) ||
          (scf->contig_pairs[i].orient == AS_OUTTIE))
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
      ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                                 AS_UID_toInteger(scf->contig_pairs[i].econtig2),
                                                 0);

      scfLen += computeGapSize(scf->contig_pairs[i].mean);
      bgn     = scfLen;
      scfLen += ctgInfo[ctgIID].len;
      end     = scfLen;
      ori     = ORIF;

      if ((scf->contig_pairs[i].orient == AS_ANTI) ||
          (scf->contig_pairs[i].orient == AS_INNIE))
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

  fprintf(scflkg, "%s\t%s\t%c\t%f\t%f\t%d\n",
          AS_UID_toString(slk->escaffold1),
          AS_UID_toString(slk->escaffold2),
          slk->orientation,
          slk->mean_distance,
          slk->std_deviation,
          slk->num_contributing);
}










void
transferFrg(FILE *ctg, FILE *scf) {
  char   frgUID[1024] = {0};
  char   ctgUID[1024] = {0};
  int    bgn = 0;
  int    end = 0;
  char   ori = 0;

  rewind(ctg);

  while (!feof(ctg)) {
    if (5 == fscanf(ctg, "%s %s %d %d %c ", frgUID, ctgUID, &bgn, &end, &ori)) {
      AS_UID uid = AS_UID_lookup(ctgUID, NULL);
      AS_IID iid = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                                     AS_UID_toInteger(uid),
                                                     0);

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

        if (ori == ctgInfo[iid].scfOri)
          ori = ORIF;
        else
          ori = ORIR;

        fprintf(scf, "%s\t%s\t%d\t%d\t%c\n",
                frgUID,
                AS_UID_toString(uid),
                bgn, end, ori);
      }
    }
  }
}





void
transferVar(FILE *ctg, FILE *scf) {
  char   varseq[16384] = {0};
  char   ctgUID[1024]  = {0};
  int    bgn = 0;
  int    end = 0;
  int    num_reads = 0;
  int    num_conf_alleles = 0;
  int    min_anchor_size = 0;
  int    var_length = 0;
  char   nr_conf_alleles[16384] = {0};
  char   weights[16384] = {0};
  char   conf_read_iids[16384] = {0};

  rewind(ctg);

  while (!feof(ctg)) {
    if (11 == fscanf(ctg, "%s %s %d %d %d %d %d %d %s %s %s ",
                     varseq, ctgUID,
                     &bgn, &end, &num_reads, &num_conf_alleles, &min_anchor_size, &var_length,
                     nr_conf_alleles,
                     weights,
                     conf_read_iids)) {
      AS_UID uid = AS_UID_lookup(ctgUID, NULL);
      AS_IID iid = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                                     AS_UID_toInteger(uid),
                                                     0);

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

        fprintf(scf, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n",
                varseq, AS_UID_toString(uid),
                bgn, end, num_reads, num_conf_alleles, min_anchor_size, var_length,
                nr_conf_alleles,
                weights,
                conf_read_iids);
      }
    }
  }
}




int main (int argc, char *argv[]) {
  char *outputPrefix       = NULL;

  GenericMesg *pmesg       = NULL;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-h") == 0) {
      err++;

    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((outputPrefix == NULL) || (err)) {
    fprintf(stderr, "usage: %s -o prefix [-h] < prefix.asm\n", argv[0]);
    fprintf(stderr, "  -o prefix        write the output here\n");
    fprintf(stderr, "  -h               print help\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  uid2iid = CreateScalarHashTable_AS();

  ctgInfo = (ctgInfo_t *)safe_calloc(ctgInfoMax, sizeof(ctgInfo_t));


  frags  = openFile("frags",  outputPrefix, 1);
  mates  = openFile("mates",  outputPrefix, 1);

  utglkg = openFile("utglkg", outputPrefix, 1);
  ctglkg = openFile("ctglkg", outputPrefix, 1);
  scflkg = openFile("scflkg", outputPrefix, 1);

  frgutg = openFile("frgutg", outputPrefix, 1);

  frgdeg = openFile("frgdeg", outputPrefix, 1);
  utgdeg = openFile("utgdeg", outputPrefix, 1);
  vardeg = openFile("vardeg", outputPrefix, 1);

  frgctg = openFile("frgctg", outputPrefix, 1);
  utgctg = openFile("utgctg", outputPrefix, 1);
  varctg = openFile("varctg", outputPrefix, 1);

  frgscf = openFile("frgscf", outputPrefix, 1);
  utgscf = openFile("utgscf", outputPrefix, 1);
  ctgscf = openFile("ctgscf", outputPrefix, 1);
  varscf = openFile("varscf", outputPrefix, 1);

  utglen = openFile("utglen", outputPrefix, 1);
  deglen = openFile("deglen", outputPrefix, 1);
  ctglen = openFile("ctglen", outputPrefix, 1);
  scflen = openFile("scflen", outputPrefix, 1);


  while(ReadProtoMesg_AS(stdin,&pmesg) != EOF){
    switch(pmesg->t){
      case MESG_MDI:
        //processMDI(pmesg->m);
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

  fclose(frags);
  fclose(mates);

  fclose(frgutg);

  fclose(frgdeg);
  fclose(utgdeg);
  fclose(vardeg);

  fclose(frgctg);
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

  DeleteHashTable_AS(uid2iid);

  exit(0);
}
