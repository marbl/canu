
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * The Celera Assembler is free software; you can redistribute it and/or modify
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

const char *mainid = "$Id: asmOutputFasta.c,v 1.13 2008-12-05 19:06:12 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_fasta.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_reverseComplement.h"

#include "MultiAlign.h"
#include "MultiAlignStore_CNS.h"


//  Parameters
//
int    minNumFragsInUnitig  = 0;

FILE   *UTGseqout = NULL;  //  bases; .fasta
FILE   *UTGqltout = NULL;  //  CA encoded quality values; .qv
FILE   *UTGquaout = NULL;  //  NCBI encoded quality values; .qual

FILE   *DEGseqout = NULL;
FILE   *DEGqltout = NULL;
FILE   *DEGquaout = NULL;

FILE   *CCOseqout = NULL;
FILE   *CCOqltout = NULL;
FILE   *CCOquaout = NULL;

FILE   *SCFseqout = NULL;
FILE   *SCFqltout = NULL;
FILE   *SCFquaout = NULL;


//  Shared Data
//
typedef struct {
  char  *cns;
  char  *qlt;
  int    len;
  int    isDegenerate;
} ctgData_t;

int             ctgDataMax = 0;
ctgData_t     **ctgData    = NULL;

HashTable_AS   *uid2iid;




int32
computeGapSize(double gapsize) {
  if (gapsize <= 20.0)
    return(20);
  return((int32)gapsize);
}



void
getseq(char *seq, char *qlt, MultiAlignT *ma) {
  char *c = Getchar(ma->consensus, 0);
  char *q = Getchar(ma->quality,   0);
  int   p = 0;

  while (*c != 0) {
    if (*c != '-') {
      seq[p] = *c;
      qlt[p] = *q;
      p++;
    }
    c++;
    q++;
  }

  seq[p] = 0;
  qlt[p] = 0;
}


//  De-gap the sequence in-place
int
degap(char *seq, char *qlt) {
  int n, i;

  for (i=n=0; seq[i]; i++)
    if (seq[i] != '-') {
      seq[n] = seq[i];
      qlt[n] = qlt[i];
      n++;
    }

  seq[n] = 0;
  qlt[n] = 0;

  return(n);
}


void
processIUM(IntUnitigMesg *ium_mesg) {
  assert(ium_mesg->length == strlen(ium_mesg->consensus));
  assert(ium_mesg->consensus[ium_mesg->length] == '\0');

  int len = degap(ium_mesg->consensus, ium_mesg->quality);

  if ((UTGseqout) && (ium_mesg->num_frags > minNumFragsInUnitig)) {
    AS_UTL_writeFastA(UTGseqout,
                      ium_mesg->consensus, len,
                      ">ium"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n",
                      ium_mesg->iaccession,
                      strlen(ium_mesg->consensus),
                      ium_mesg->num_frags,
                      ium_mesg->coverage_stat);
    AS_UTL_writeFastA(UTGqltout,
                      ium_mesg->quality, len,
                      ">ium"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n",
                      ium_mesg->iaccession,
                      strlen(ium_mesg->consensus),
                      ium_mesg->num_frags,
                      ium_mesg->coverage_stat);
    AS_UTL_writeQVFastA(UTGquaout,
                        ium_mesg->quality, len,
                        ">ium"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n",
                        ium_mesg->iaccession,
                        strlen(ium_mesg->consensus),
                        ium_mesg->num_frags,
                        ium_mesg->coverage_stat);
  }
}



void
processUTG(SnapUnitigMesg  *utg_mesg) {
  assert(utg_mesg->length == strlen(utg_mesg->consensus));
  assert(utg_mesg->consensus[utg_mesg->length] == '\0');

  int len = degap(utg_mesg->consensus, utg_mesg->quality);

  if ((UTGseqout) && (utg_mesg->num_frags > minNumFragsInUnitig)) {
    AS_UTL_writeFastA(UTGseqout,
                      utg_mesg->consensus, len,
                      ">utg%s length=%d num_frags="F_IID" Astat=%.2f\n",
                      AS_UID_toString(utg_mesg->eaccession),
                      strlen(utg_mesg->consensus),
                      utg_mesg->num_frags,
                      utg_mesg->coverage_stat);
    AS_UTL_writeFastA(UTGqltout,
                      utg_mesg->quality, len,
                      ">utg%s length=%d num_frags="F_IID" Astat=%.2f\n",
                      AS_UID_toString(utg_mesg->eaccession),
                      strlen(utg_mesg->consensus),
                      utg_mesg->num_frags,
                      utg_mesg->coverage_stat);
    AS_UTL_writeQVFastA(UTGquaout,
                        utg_mesg->quality, len,
                        ">utg%s length=%d num_frags="F_IID" Astat=%.2f\n",
                        AS_UID_toString(utg_mesg->eaccession),
                        strlen(utg_mesg->consensus),
                        utg_mesg->num_frags,
                        utg_mesg->coverage_stat);
  }
}



void
processCCO(SnapConConMesg *cco_mesg) {

  assert(strlen(cco_mesg->consensus) == cco_mesg->length);

  InsertInHashTable_AS(uid2iid,
                       AS_UID_toInteger(cco_mesg->eaccession), 0,
                       cco_mesg->iaccession, 0);

  ctgData_t   *cd = (ctgData_t *)safe_malloc(sizeof(ctgData_t));

  cd->cns          = (char *)safe_malloc(sizeof(char) * (cco_mesg->length + 1));
  cd->qlt          = (char *)safe_malloc(sizeof(char) * (cco_mesg->length + 1));
  cd->len          = cco_mesg->length;
  cd->isDegenerate = FALSE;

  //  By definition, a degenerate contig has one unitig and is unplaced.
  //  In reality, those two conditions always occur together.
  //
  cd->isDegenerate = (cco_mesg->placed == AS_UNPLACED) && (cco_mesg->num_unitigs == 1);

  memcpy(cd->cns, cco_mesg->consensus, cco_mesg->length);
  memcpy(cd->qlt, cco_mesg->quality,   cco_mesg->length);

  cd->cns[cco_mesg->length] = 0;
  cd->qlt[cco_mesg->length] = 0;

  cd->len = degap(cd->cns, cd->qlt);

  //fprintf(stderr, "CCO: "F_IID" len: %d -> %d\n", cco_mesg->iaccession, cco_mesg->length, cd->len);

  if (ctgDataMax <= cco_mesg->iaccession) {
    ctgDataMax = cco_mesg->iaccession + 32 * 1024 * 1024;
    ctgData    = (ctgData_t **)safe_realloc(ctgData, sizeof(ctgData_t) * ctgDataMax);
  }

  ctgData[cco_mesg->iaccession] = cd;

  if ((cd->isDegenerate == 0) && (CCOseqout)) {
    AS_UTL_writeFastA(CCOseqout,
                      cd->cns, cd->len,
                      ">ctg%s\n",
                      AS_UID_toString(cco_mesg->eaccession));
    AS_UTL_writeFastA(CCOqltout,
                      cd->qlt, cd->len,
                      ">ctg%s\n",
                      AS_UID_toString(cco_mesg->eaccession));
    AS_UTL_writeQVFastA(CCOquaout,
                        cd->qlt, cd->len,
                        ">ctg%s\n",
                        AS_UID_toString(cco_mesg->eaccession));
  }
  if ((cd->isDegenerate == 1) && (DEGseqout)) {
    AS_UTL_writeFastA(DEGseqout,
                      cd->cns, cd->len,
                      ">deg%s\n",
                      AS_UID_toString(cco_mesg->eaccession));
    AS_UTL_writeFastA(DEGqltout,
                      cd->qlt, cd->len,
                      ">deg%s\n",
                      AS_UID_toString(cco_mesg->eaccession));
    AS_UTL_writeQVFastA(DEGquaout,
                        cd->qlt, cd->len,
                        ">deg%s\n",
                        AS_UID_toString(cco_mesg->eaccession));
  }
}



void
processSCF(SnapScaffoldMesg *scf_mesg) {
  AS_IID        ctgIID = 0;
  int           stfLen = 0;
  MultiAlignT  *ctg    = NULL;
  int           ctgLen = 0;
  int           scfLen = 0;
  int           scfPos = 0;
  int           i;

  //  Probably should process it, but I don't see any error checking below, so why?
  if (SCFseqout == NULL)
    return;

  //  Orientation of each contig

  int  *reversed = (int *)safe_malloc(sizeof(int) * (scf_mesg->num_contig_pairs + 1));

  if (scf_mesg->num_contig_pairs == 0) {
    reversed[0] = FALSE;
  } else {
    int i;
    reversed[0] = ((scf_mesg->contig_pairs[0].orient == BA_AB) || (scf_mesg->contig_pairs[0].orient == BA_BA)) ? TRUE : FALSE;
    for (i=0; i<scf_mesg->num_contig_pairs; i++)
      reversed[i+1] = ((scf_mesg->contig_pairs[i].orient == BA_AB) || (scf_mesg->contig_pairs[i].orient == AB_BA)) ? !reversed[i] : reversed[i];
  }


  //  Total length of the scaffold

  ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                             AS_UID_toInteger(scf_mesg->contig_pairs[0].econtig1),
                                             0);
  scfLen = ctgData[ctgIID]->len;

  for (i=0; i<scf_mesg->num_contig_pairs; i++) {
    ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                               AS_UID_toInteger(scf_mesg->contig_pairs[i].econtig2),
                                               0);
    scfLen += computeGapSize(scf_mesg->contig_pairs[i].mean) + ctgData[ctgIID]->len;
  }


  //  Build scaffold sequence.  Append first contig, then iterate, appending gap then contig.

  char  *scfcns = (char *)safe_malloc(sizeof(char) * (scfLen + 1));
  char  *scfqlt = (char *)safe_malloc(sizeof(char) * (scfLen + 1));

  ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                             AS_UID_toInteger(scf_mesg->contig_pairs[0].econtig1),
                                             0);

  memcpy(scfcns + scfPos, ctgData[ctgIID]->cns, ctgData[ctgIID]->len);
  memcpy(scfqlt + scfPos, ctgData[ctgIID]->qlt, ctgData[ctgIID]->len);

  if (reversed[0])
    reverseComplement(scfcns + scfPos, scfqlt + scfPos, ctgData[ctgIID]->len);

  scfPos += ctgData[ctgIID]->len;


  for (i=0; i<scf_mesg->num_contig_pairs; i++) {
    int  gapSize = computeGapSize(scf_mesg->contig_pairs[i].mean);

    memset(scfcns + scfPos, 'N', gapSize);
    memset(scfqlt + scfPos, '0', gapSize);
    scfPos += gapSize;

    ctgIID = (AS_IID)LookupValueInHashTable_AS(uid2iid,
                                               AS_UID_toInteger(scf_mesg->contig_pairs[i].econtig2),
                                               0);

    memcpy(scfcns + scfPos, ctgData[ctgIID]->cns, ctgData[ctgIID]->len);
    memcpy(scfqlt + scfPos, ctgData[ctgIID]->qlt, ctgData[ctgIID]->len);

    if (reversed[i+1])
      reverseComplement(scfcns + scfPos, scfqlt + scfPos, ctgData[ctgIID]->len);

    scfPos += ctgData[ctgIID]->len;
  }

  //  Output

  AS_UTL_writeFastA(SCFseqout,
                    scfcns, scfPos,
                    ">scf%s\n", AS_UID_toString(scf_mesg->eaccession));
  AS_UTL_writeFastA(SCFqltout,
                    scfqlt, scfPos,
                    ">scf%s\n", AS_UID_toString(scf_mesg->eaccession));
  AS_UTL_writeQVFastA(SCFquaout,
                      scfqlt, scfPos,
                      ">scf%s\n", AS_UID_toString(scf_mesg->eaccession));

  safe_free(reversed);
  safe_free(scfcns);
  safe_free(scfqlt);
}




FILE *openOutput(char *prefix, char *template) {
  char  N[FILENAME_MAX];
  FILE *F;

  sprintf(N, template, prefix);

  errno = 0;
  F = fopen(N, "w");
  if (errno) {
    fprintf(stderr, "Failed to open output '%s': %s\n", N, strerror(errno));
    exit(1);
  }

  return(F);
}


int
main(int argc, char **argv) {
  GenericMesg      *pmesg;
  char              name[FILENAME_MAX];
  char             *prefix = NULL;
  char             *infile = NULL;

  int               dumpUnitigs     = 1;
  int               dumpDegenerates = 1;
  int               dumpContigs     = 1;
  int               dumpScaffolds   = 1;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-U") == 0) {
      dumpUnitigs = 0;
    } else if (strcmp(argv[arg], "-D") == 0) {
      dumpDegenerates = 0;
    } else if (strcmp(argv[arg], "-C") == 0) {
      dumpContigs = 0;
    } else if (strcmp(argv[arg], "-S") == 0) {
      dumpScaffolds = 0;
    } else if (strcmp(argv[arg], "-p") == 0) {
      prefix = argv[++arg];
    } else if (strcmp(argv[arg], "-n") == 0) {
      minNumFragsInUnitig = atoi(argv[++arg]);
    } else {
      if (infile == NULL) {
        infile = argv[arg];
      } else {
        fprintf(stderr, "unknown option '%s'\n", argv[arg]);
        err = 1;
      }
    }
    arg++;
  }
  if (prefix == NULL)
    err++;
  if (isatty(fileno(stdin)) && (infile == NULL))
    err++;
  if (err > 0) {
    fprintf(stderr, "usage: %s [options] -p prefix < asmfile\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p         write files named 'prefix.XXX.TYPE', etc.\n");
    fprintf(stderr, "                 XXX =  type of object\n");
    fprintf(stderr, "                        utg - unitig\n");
    fprintf(stderr, "                        deg - degenerate\n");
    fprintf(stderr, "                        ctg - contig\n");
    fprintf(stderr, "                        scf - scaffold\n");
    fprintf(stderr, "                 TYPE = type of sequence\n");
    fprintf(stderr, "                        fasta - acgt bases\n");
    fprintf(stderr, "                        qv    - Celera encoded quality values\n");
    fprintf(stderr, "                        qual  - NCBI encoded quality values\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U         do NOT dump unitigs\n");
    fprintf(stderr, "  -D         do NOT dump degenerates\n");
    fprintf(stderr, "  -C         do NOT dump contigs\n");
    fprintf(stderr, "  -S         do NOT dump scaffolds\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "UNITIG OPTIONS\n");
    fprintf(stderr, "  -n nf      dump only unitigs with at least nf reads\n");
    fprintf(stderr, "             in them.  Default is 0 (dump all unitigs).\n");
    fprintf(stderr, "\n");
#if 0
    fprintf(stderr, "\n");
    fprintf(stderr, "CONTIG OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SCAFFOLD OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
#endif
    exit(1);
  }

  uid2iid = CreateScalarHashTable_AS();

  if (dumpUnitigs) {
    UTGseqout = openOutput(prefix, "%s.utg.fasta");
    UTGqltout = openOutput(prefix, "%s.utg.qv");
    UTGquaout = openOutput(prefix, "%s.utg.qual");
  }
  if (dumpContigs) {
    DEGseqout = openOutput(prefix, "%s.deg.fasta");
    DEGqltout = openOutput(prefix, "%s.deg.qv");
    DEGquaout = openOutput(prefix, "%s.deg.qual");
  }
  if (dumpContigs) {
    CCOseqout = openOutput(prefix, "%s.ctg.fasta");
    CCOqltout = openOutput(prefix, "%s.ctg.qv");
    CCOquaout = openOutput(prefix, "%s.ctg.qual");
  }
  if (dumpScaffolds) {
    SCFseqout = openOutput(prefix, "%s.scf.fasta");
    SCFqltout = openOutput(prefix, "%s.scf.qv");
    SCFquaout = openOutput(prefix, "%s.scf.qual");
  }

  FILE *F = stdin;
  if (infile)
    F = fopen(infile, "r");

  while (ReadProtoMesg_AS(F, &pmesg) != EOF) {
    switch (pmesg->t) {
      case MESG_IUM:
        processIUM((IntUnitigMesg *)pmesg->m);
        break;
      case MESG_UTG:
        processUTG((SnapUnitigMesg *)pmesg->m);
        break;
      case MESG_CCO:
        processCCO((SnapConConMesg *)pmesg->m);
        break;
      case MESG_SCF:
        processSCF((SnapScaffoldMesg *)pmesg->m);
        break;
    }
  }

  errno = 0;

  if (dumpUnitigs) {
    fclose(UTGseqout);
    fclose(UTGqltout);
    fclose(UTGquaout);
  }
  if (dumpDegenerates) {
    fclose(DEGseqout);
    fclose(DEGqltout);
    fclose(DEGquaout);
  }
  if (dumpContigs) {
    fclose(CCOseqout);
    fclose(CCOqltout);
    fclose(CCOquaout);
  }
  if (dumpScaffolds) {
    fclose(SCFseqout);
    fclose(SCFqltout);
    fclose(SCFquaout);
  }

  if (errno) {
    fprintf(stderr, "WARNING:  Some files didn't close properly -- check disk space!\n");
    return(1);
  }

  fprintf(stderr, "All output successfully written.\n");

  return(0);
}
