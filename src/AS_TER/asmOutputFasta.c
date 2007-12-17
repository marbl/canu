
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

// $ Id:  $

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_fasta.h"
#include "AS_UTL_Hash.h"

#include "MultiAlign.h"
#include "MultiAlignStore_CNS.h"


//  Parameters
//
int    minNumFragsInUnitig  = 0;
char  *prefix               = NULL;

FILE   *UTGseqout = NULL;
FILE   *UTGqltout = NULL;

FILE   *CCOseqout = NULL;
FILE   *CCOqltout = NULL;

FILE   *DSCseqout = NULL;
FILE   *DSCqltout = NULL;

FILE   *SCFseqout = NULL;
FILE   *SCFqltout = NULL;


//  Shared Data
//
typedef struct {
  char  *cns;
  char  *qlt;
  int    len;
  int    isPlaced;
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



//  Reverse complement the ACGT sequence in seq, and reverse the QV
//  sequence in qlt.
//
void
reverseComplement(char *seq, char *qlt, int len) {
  char   inv[256] = {0};
  char   c=0;
  char  *s=seq,  *S=seq+len-1;
  char  *q=qlt,  *Q=qlt+len-1;

  inv['a'] = 't';
  inv['c'] = 'g';
  inv['g'] = 'c';
  inv['t'] = 'a';
  inv['n'] = 'n';
  inv['A'] = 'T';
  inv['C'] = 'G';
  inv['G'] = 'C';
  inv['T'] = 'A';
  inv['N'] = 'N';
  inv['-'] = '-';

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }

  if (s == S)
    *s = inv[*s];
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

  if (ium_mesg->num_frags > minNumFragsInUnitig) {
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
  }
}



void
processUTG(SnapUnitigMesg  *utg_mesg) {
  assert(utg_mesg->length == strlen(utg_mesg->consensus));
  assert(utg_mesg->consensus[utg_mesg->length] == '\0');

  int len = degap(utg_mesg->consensus, utg_mesg->quality);

  if (utg_mesg->num_frags > minNumFragsInUnitig) {
    AS_UTL_writeFastA(UTGseqout,
                      utg_mesg->consensus, len,
                      ">utg"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n",
                      utg_mesg->iaccession,
                      strlen(utg_mesg->consensus),
                      utg_mesg->num_frags,
                      utg_mesg->coverage_stat);
    AS_UTL_writeFastA(UTGqltout,
                      utg_mesg->quality, len,
                      ">utg"F_IID" length=%d num_frags="F_IID" Astat=%.2f\n",
                      utg_mesg->iaccession,
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
  cd->isPlaced     = FALSE;
  cd->isDegenerate = FALSE;
    
  //  By definition, a degenerate contig has one unitig and is
  //  unplaced.
  //
  //  one version used (numUtg == 1) && (unitig is not placed)
  //  one version used (numUtg == 1) && (contig is not placed)
  //
  //AS_IID utgIID =  (AS_IID)LookupValueInHashTable_AS(uid2iid, AS_UID_toInteger(cco_mesg->unitigs[0].eident), 0);
  //if (isPlacedUnitig[utgIID])  //  true if (utg_mesg->status == AS_SEP)
  //  cd->notDegenerate = TRUE;

  if (cco_mesg->placed  == AS_UNPLACED)
    cd->isPlaced = TRUE;

  if ((cco_mesg->num_unitigs == 1) &&
      (cco_mesg->placed      == AS_UNPLACED))
    cd->isDegenerate = TRUE;

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

  AS_UTL_writeFastA(CCOseqout,
                    cd->cns, cd->len,
                    ">ctg%s,"F_IID" placed=%s degenerate=%s\n",
                    AS_UID_toString(cco_mesg->eaccession),
                    cco_mesg->iaccession,
                    (cd->isPlaced)     ? "true" : "false",
                    (cd->isDegenerate) ? "true" : "false");
  AS_UTL_writeFastA(CCOqltout,
                    cd->qlt, cd->len,
                    ">ctg%s,"F_IID" placed=%s degenerate=%s\n",
                    AS_UID_toString(cco_mesg->eaccession),
                    cco_mesg->iaccession,
                    (cd->isPlaced)     ? "true" : "false",
                    (cd->isDegenerate) ? "true" : "false");
}



void
processDSC(SnapDegenerateScaffoldMesg *dsc_mesg) {
  AS_IID        ctgIID = LookupValueInHashTable_AS(uid2iid, AS_UID_toInteger(dsc_mesg->econtig), 0);
  ctgData_t    *cd     = ctgData[ctgIID];

  if (cd->isDegenerate == FALSE)
    return;

  //fprintf(stderr, "DSC: "F_IID" %s len: %d\n", ctgIID, AS_UID_toString(dsc_mesg->econtig), cd->len);

  AS_UTL_writeFastA(DSCseqout,
                    cd->cns, cd->len,
                    ">dsc%s\n", AS_UID_toString(dsc_mesg->econtig));
  AS_UTL_writeFastA(DSCqltout,
                    cd->qlt, cd->len,
                    ">dsc%s\n", AS_UID_toString(dsc_mesg->econtig));
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

  //  Orientation of each contig

  int  *reversed = (int *)safe_malloc(sizeof(int) * (scf_mesg->num_contig_pairs + 1));

  if (scf_mesg->num_contig_pairs == 0) {
    reversed[0] = FALSE;
  } else {
    int i;
    reversed[0] = ((scf_mesg->contig_pairs[0].orient == BA_AB) || (scf_mesg->contig_pairs[0].orient == AB_BA)) ? TRUE : FALSE;
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

  safe_free(reversed);
  safe_free(scfcns);
  safe_free(scfqlt);
}




FILE *openOutput(char *prefix, char *label) {
  char  N[FILENAME_MAX];
  FILE *F;

  sprintf(N, "%s.%s.fasta", prefix, label);

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
  char             *infile = NULL;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-U") == 0) {
    } else if (strcmp(argv[arg], "-C") == 0) {
    } else if (strcmp(argv[arg], "-S") == 0) {
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
  if (isatty(fileno(stdin)) && (infile == NULL))
    err++;
  if (err > 0) {
    fprintf(stderr, "usage: %s [options] -p prefix < asmfile\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -U         dump unitigs\n");
    fprintf(stderr, "  -C         dump contigs\n");
    fprintf(stderr, "  -S         dump scaffolds\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -p         write files named 'prefix.contig.fasta', etc.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "UNITIG OPTIONS\n");
    fprintf(stderr, "  -n nf      dump only unitigs with at least nf reads\n");
    fprintf(stderr, "             in them.  Default is 0 (dump all unitigs).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CONTIG OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "SCAFFOLD OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  uid2iid = CreateScalarHashTable_AS(32 * 1024);
  
  UTGseqout = openOutput(prefix, "utgcns");
  UTGqltout = openOutput(prefix, "utgqlt");

  CCOseqout = openOutput(prefix, "ctgcns");
  CCOqltout = openOutput(prefix, "ctgqlt");

  DSCseqout = openOutput(prefix, "dsccns");
  DSCqltout = openOutput(prefix, "dscqlt");

  SCFseqout = openOutput(prefix, "scfcns");
  SCFqltout = openOutput(prefix, "scfqlt");

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
      case MESG_DSC:
        processDSC((SnapDegenerateScaffoldMesg *)pmesg->m);
        break;
      case MESG_SCF:
        processSCF((SnapScaffoldMesg *)pmesg->m);
        break;
    }
  }

  errno = 0;

  fclose(UTGseqout);
  fclose(UTGqltout);

  fclose(CCOseqout);
  fclose(CCOqltout);

  fclose(DSCseqout);
  fclose(DSCqltout);

  fclose(SCFseqout);
  fclose(SCFqltout);

  if (errno) {
    fprintf(stderr, "WARNING:  Some files didn't close properly -- check disk space!\n");
    return(1);
  }

  fprintf(stderr, "All output successfully written.\n");

  return(0);
}
