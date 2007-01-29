
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

/* based on AS_CGW/ProcessScaffolds_CGW.c */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <unistd.h>
#include <dirent.h>
#include "assert.h"
#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "AS_UTL_ID_store.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

static UIDHashTable_AS *uid2iid;

FILE *fastaFile = NULL;
FILE *fastaDregsFile = NULL;
FILE *fastaQualsFile = NULL;
char fastaFileName[FILENAME_MAX] = {0};
char fastaDregsFileName[FILENAME_MAX] = {0};
char fastaQualsFileName[FILENAME_MAX] = {0};
int fasta = 0;
int fastaDregs = 0;
int fastaQuals = 0;
char *fastaIdent = NULL;

VA_TYPE(int32) *reversedVA        = NULL;
VA_TYPE(char)  *ctmp              = NULL;
VA_TYPE(char)  *qtmp              = NULL;
VA_TYPE(char)  *scaffold_sequence = NULL;
VA_TYPE(char)  *quality_sequence  = NULL;

MultiAlignStoreT *cstore;

VA_DEF(MultiPos);
VA_DEF(UnitigPos);
VA_DEF(SnapContigLinkMesg);

void CleanExit(int rc) {
  if( fastaFile != NULL ){
    fclose(fastaFile);
    unlink(fastaFileName);
  }
  exit(rc);
}

int IsForward(SeqInterval s) {
  return (s.bgn < s.end);
}

static void Complement(char *in, int len);
static void Reverse(char *in, int len);


int32 compute_gap(double gapsize){
  if (gapsize <= 20.0)
    return(20);
  return (int32)gapsize;
}


int FastaScaffold(FILE *out, FILE *qual, SnapScaffoldMesg *scaff){
  int               i;
  int               num_pairs = scaff->num_contig_pairs; 
  int32             contigid;
  MultiAlignT      *contig;
  SnapContigPairs  *cp = scaff->contig_pairs;
  int32            *reversed;
  char             *sseq;  //  scaffold sequence
  char             *qseq;  //  quality sequence
  int               running_length  = 0;
  char              nchar           = 'n';
  int               scaffold_length;
  int               contig_length;
  int               flen;
  int               ngaps;

  ResetVA_int32(reversedVA);
  ResetVA_char(ctmp);
  ResetVA_char(qtmp);
  ResetVA_char(scaffold_sequence);
  ResetVA_char(quality_sequence);

  reversed = Getint32(reversedVA,0);
  contigid = *LookupInUID2IIDHashTable_AS(uid2iid,cp[0].econtig1);

  // output label line for this scaffold
  fprintf(out,">" F_S64 " /type=%s\n",scaff->eaccession,fastaIdent);

  if (qual)
    fprintf(qual,">" F_S64 " /type=%s\n",scaff->eaccession,fastaIdent);

  // calculate length of sequence:
  //
  scaffold_length = GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,contigid));
  for ( i=0;i<num_pairs;i++ ) {
    scaffold_length += compute_gap(cp[i].mean);
    scaffold_length += GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,
									*LookupInUID2IIDHashTable_AS(uid2iid,cp[i].econtig2)));
  }

  EnableRangeVA_char(scaffold_sequence, scaffold_length+1);
  EnableRangeVA_char(quality_sequence,  scaffold_length+1);

  EnableRangeVA_int32(reversedVA, num_pairs + 1);
  reversed = Getint32(reversedVA,0);

  // Compute the orientation of each contig in the scaffold
  if(num_pairs == 0){
    reversed[0] = FALSE;
  }else{
    switch(cp[0].orient){
      case AB_AB:
      case AB_BA:
        reversed[0] = FALSE;
        break;
      case BA_AB:
      case BA_BA:
        reversed[0] = TRUE;
        break;
    }
  
    for ( i=0;i<num_pairs;i++ ) {
      switch(cp[i].orient){
        case AB_AB:
        case BA_BA:
          reversed[i+1] = reversed[i];
          break;
        case BA_AB:
        case AB_BA:
          reversed[i+1] = !reversed[i];
          break;
      }
    }
  }

  // append first contig to scaffold_sequence
  contig = GetMultiAlignInStore(cstore,contigid);
  assert(contig!=NULL);
  contig_length = GetMultiAlignUngappedLength(contig);
  GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);

  if(reversed[0]) {
    Complement(Getchar(ctmp,0), contig_length);
    Reverse(Getchar(qtmp,0), contig_length);
  }

  running_length+=contig_length;
  if ((running_length > GetNumchars(scaffold_sequence)) || (running_length > GetNumchars(quality_sequence)))
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");

  sseq = Getchar(scaffold_sequence,0);
  qseq = Getchar(quality_sequence,0);

  memcpy(sseq, Getchar(ctmp,0), contig_length);
  sseq+=contig_length;

  memcpy(qseq, Getchar(qtmp,0), contig_length);
  qseq+=contig_length;

  for ( i=0;i<num_pairs;i++ ) {
    ngaps = compute_gap(cp[i].mean);

    // append ngaps N's (or n's, depending on alternate_gap_spec)  to string as intercontig space

    memset(sseq, nchar, ngaps);
    memset(qseq, '0', ngaps);
    sseq+=ngaps;
    qseq+=ngaps;

    contigid      = *LookupInUID2IIDHashTable_AS(uid2iid,cp[i].econtig2);
    contig        = GetMultiAlignInStore(cstore,contigid);
    contig_length = GetMultiAlignUngappedLength(contig);

    GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);

    if(reversed[i+1]) {
      Complement(Getchar(ctmp,0), contig_length);
      Reverse(Getchar(qtmp,0), contig_length);
    }

    running_length += ngaps + contig_length;

    if ((running_length > GetNumchars(scaffold_sequence)) ||
        (running_length > GetNumchars(quality_sequence))) {
      fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
    }

    memcpy(sseq, Getchar(ctmp,0), contig_length);
    memcpy(qseq, Getchar(qtmp,0), contig_length);

    sseq += contig_length;
    qseq += contig_length;
  }

  // now, output the scaffold
  sseq = Getchar(scaffold_sequence,0);
  qseq = Getchar(quality_sequence,0);

  flen = strlen(sseq);

  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,sseq);
    sseq += 70;

    if (qual) {
      fprintf(qual,"%.*s\n",70,qseq);
      qseq += 70;
    }
    
  }
  return 1;
}

int FastaDegenerateScaffold(FILE *out, SnapDegenerateScaffoldMesg *scaff,VA_TYPE(int32) *is_placed)  {
  int32 contigid =  *LookupInUID2IIDHashTable_AS(uid2iid,scaff->econtig);
  int i;
  MultiAlignT *contig;
  char *sseq;
  int running_length=0;
  int scaffold_length,contig_length,flen;

  ResetVA_int32(reversedVA);
  ResetVA_char(ctmp);
  ResetVA_char(qtmp);
  ResetVA_char(scaffold_sequence);
  ResetVA_char(quality_sequence);

  // calculate length of sequence:
  scaffold_length = GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,contigid));

  EnableRangeVA_char(scaffold_sequence,scaffold_length+1);

  // append first contig to scaffold_sequence
  contig = GetMultiAlignInStore(cstore,contigid);

  contig_length = GetMultiAlignUngappedLength(contig);

  if ( GetNumUnitigPoss(contig->u_list) == 1 ) {
    uint32 iid=  *LookupInUID2IIDHashTable_AS(uid2iid,GetUnitigPos(contig->u_list,0)->eident);
    if ( *(Getint32(is_placed,iid)) ) return 0;
  }
  GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);

  sseq = Getchar(scaffold_sequence,0);

  running_length+=contig_length;
  if (running_length > GetNumchars(scaffold_sequence)) {
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
  }

  // output label line for this scaffold
  fprintf(out,">" F_S64 " /type=%s_Degenerate\n",scaff->econtig,fastaIdent);

  memcpy(sseq, Getchar(ctmp,0),contig_length);

  flen = strlen(sseq);
  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,sseq);
    sseq += 70;
  }
  return 1;
}

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-h] [-d] [-q] -f fasta.output.filename\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, " -f:  produce a multi-fasta file, with one fasta record for each scaffold,\n");
  fprintf(stderr, "      and 'n's used as intercontig gap placeholders.  For gaps > 50 bp,\n");
  fprintf(stderr, "      the number of 'n's used represents the calculated mean gap size.\n");
  fprintf(stderr, "      For others, 50 'n's are used.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -d:  also produce a multi-fasta dregs file, with one fasta record for\n");
  fprintf(stderr, "      each degenerate scaffold.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " -q:  also produce a multi-fasta quality file, with one fasta record for\n");
  fprintf(stderr, "      each scaffold.\n");
  fprintf(stderr, "\n");
  CleanExit(1);
}

int main(int argc, char *argv[]) {
  GenericMesg *pmesg;
  SnapConConMesg *contig;
  SnapUnitigMesg *unitig;
  MultiAlignT *ma;
  char fastaIdentifier[FILENAME_MAX];
  char  *suffix = NULL;
  VA_TYPE(int32) *is_placed;
  int32  placed=1;
  int32  unplaced=0;
  char   ch;

  cstore       = CreateMultiAlignStoreT(0);
  is_placed    = CreateVA_int32(0);

  optarg       = NULL;
  while ( ((ch = getopt(argc, argv, "dqf:h?")) != EOF)) {
    switch(ch) {
      case 'd':
        fastaDregs = 1;
        break;
      case 'q':
        fastaQuals = 1;
        break;
      case 'f':
        strcpy(fastaFileName, optarg);
        fastaFile = fopen(fastaFileName,"w");
        if (fastaFile == NULL ) {
          fprintf(stderr,"Failure to create fasta file %s\n", fastaFileName);
          CleanExit(1);
        }

        strcpy(fastaIdentifier,fastaFileName);
        fastaIdent = strrchr(fastaIdentifier,'/');
        if ( fastaIdent == NULL ) { 
          fastaIdent = fastaIdentifier;
        } else {
          fastaIdent++;
        }

        suffix = strrchr(fastaIdentifier,(int)'.'); 
        if(suffix!=NULL)
          *suffix = '\0'; // this cuts off the ext, so filename root can be used

        while ( (suffix = strchr(fastaIdent,'.')) != NULL )
          *suffix = '_';

        break;
      case 'h':
        usage(argv[0]);
        break;
      default:
        fprintf(stderr,"Invalid option -%c, try -h for usage\n",ch);
        CleanExit(1);
        break;
    }
  }

  if (fastaFile == 0L)
    usage(argv[0]);

  uid2iid = CreateUIDHashTable_AS(200000);

  if (fastaDregs){
    sprintf(fastaDregsFileName,"%s.dregs",fastaFileName);
    fastaDregsFile = fopen(fastaDregsFileName,"w");
    if (fastaDregsFile == NULL ) {
      fprintf(stderr,"Failure to create fasta file %s\n", fastaDregsFileName);
      CleanExit(1);
    }
  }

  if (fastaQuals){
    sprintf(fastaQualsFileName,"%s.quals",fastaFileName);
    fastaQualsFile = fopen(fastaQualsFileName,"w");
    if (fastaQualsFile == NULL ) {
      fprintf(stderr,"Failure to create fasta file %s\n", fastaQualsFileName);
      CleanExit(1);
    }
  }

  //  Create some variable arrays
  //
  reversedVA        = CreateVA_int32(1000);
  ctmp              = CreateVA_char(200000);
  qtmp              = CreateVA_char(200000);
  scaffold_sequence = CreateVA_char(200000);
  quality_sequence  = CreateVA_char(200000);
 
  
  while (ReadProtoMesg_AS(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_UTG)  {
      unitig = pmesg->m;
      InsertInUID2IIDHashTable_AS(uid2iid,unitig->eaccession,unitig->iaccession);
      SetVA_int32(is_placed, unitig->iaccession,  
                  ( unitig->status == AS_SEP )? &placed : &unplaced );
    }
    if (pmesg->t ==MESG_CCO)  {
      contig = pmesg->m;
      InsertInUID2IIDHashTable_AS(uid2iid,contig->eaccession,contig->iaccession);
      assert ( strlen(contig->consensus) == contig->length);
      ma = CreateMultiAlignTFromCCO(contig, -1,  0);
      SetMultiAlignInStore(cstore,ma->id,ma);
    }
    if (pmesg->t ==MESG_SCF)  {
      FastaScaffold(fastaFile, fastaQualsFile, (SnapScaffoldMesg *)pmesg->m);
    }
    if (pmesg->t ==MESG_DSC)  {
      if (fastaDregs) {
        FastaDegenerateScaffold(fastaDregsFile, (SnapDegenerateScaffoldMesg *)pmesg->m,is_placed);
      }
    }
  }

  if(fastaDregs)
    fclose(fastaDregsFile);
  if(fastaQuals)
    fclose(fastaQualsFile);

  fclose(fastaFile);

  exit(0);
}


/*** UTILITY ROUTINES ***/

/* Complement the sequence in fragment message a.  This include also
   revsersing the order of the quality values.  The operation does the
   complementation/reversal in place.  Calling it a second time on a
   given fragment restores it to its original state.                */

// Stolen from AS_ALN
static void Complement(char *in, int len)
{ static char WCinvert[256];
 static Firstime = 1;

 if (Firstime)          /* Setup complementation array */
   { 
     int i;
     Firstime = 0;
     for(i = 0; i < 256;i++){
       WCinvert[i] = '?';
     }
     WCinvert['a'] = 't';
     WCinvert['c'] = 'g';
     WCinvert['g'] = 'c';
     WCinvert['t'] = 'a';
     WCinvert['n'] = 'n';
     WCinvert['A'] = 'T';
     WCinvert['C'] = 'G';
     WCinvert['G'] = 'C';
     WCinvert['T'] = 'A';
     WCinvert['N'] = 'N';
     WCinvert['-'] = '-'; // added this to enable alignment of gapped consensi
   }
      
 { /* Complement and reverse sequence */

   { register char *s, *t;
   int c;

   s = in;
   t = in + (len-1);
   while (s < t)
     { // Sanity Check!
       assert(WCinvert[*t] != '?' &&
              WCinvert[*s] != '?');

       c = *s;
       *s++ = WCinvert[*t];
       *t-- = WCinvert[c];
     }
   if (s == t)
     *s = WCinvert[*s];
   }

 }
}


//  Stolen from Complement() above.
static void Reverse(char *in, int len) {
  register char *s, *t;
  int c;

  s = in;
  t = in + (len-1);
  while (s < t) {
    c = *s;
    *s++ = *t;
    *t-- = c;
  }
  if (s == t)
    *s = *s;
}
