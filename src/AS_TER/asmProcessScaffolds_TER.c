
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
// changed from ifndef i386 since systems other than OSF do not need mode.h (MP)
#ifdef _OSF_SOURCE
#include <sys/mode.h>
#endif
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
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

#define INTERSCAFFDIST 50

static UIDHashTable_AS *uid2iid;

FILE *fastaFile;
FILE *fastaDregsFile;
char fastaFileName[FILENAME_MAX];
char fastaDregsFileName[FILENAME_MAX];
int fasta;
int fastaDregs;
char *fastaIdent;
mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
VA_TYPE(int32) *scaff_index;
VA_TYPE(int32) *bactig_unitigs;
VA_TYPE(char) *scaffold_sequence=NULL;

VA_TYPE(int32) *cids=NULL;
VA_TYPE(int32) *reversedVA=NULL;
VA_TYPE(char) *ctmp=NULL;
VA_TYPE(char) *qtmp=NULL;

VA_DEF(MultiPos)

VA_DEF(UnitigPos)

VA_DEF(SnapContigLinkMesg)

void CleanExit(int rc) {
  char command[100+FILENAME_MAX];
  int i;
  if( fastaFile != NULL ){
    fclose(fastaFile);
    sprintf(command,"rm -f %s",fastaFileName);
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  exit(rc);
}

int HandleDir(char *filePathAndName, char *fileName) {
   mode_t mode = S_IRWXU | S_IRWXG | S_IROTH;
   char *suffix;
   char *DirName;
   char *FileName;
   DIR *Dir;
   suffix = strrchr(filePathAndName,(int)'/');
   if ( suffix != NULL ) {
      *suffix = '\0';
      DirName = filePathAndName; 
      if ( DirName != NULL ) {
        Dir = opendir(DirName);
        if ( Dir == NULL ) {
          if(mkdir(DirName,mode)){
            fprintf(stderr,"Failure to create directory %s\n", DirName);
            CleanExit(1);
          }
        }
      }
      *suffix = '/';
      FileName = filePathAndName;
    } else {
      FileName = filePathAndName;
    }
    strcpy(fileName,FileName);
    return 1;
}

MultiAlignStoreT *cstore;
uint32 global_coord;

int IsForward(SeqInterval s) {
  return (s.bgn < s.end);
}

static void Complement(char *in, int len);

#define IR_GAP_SIZE (50)

int32 compute_gap(double gapsize){
  if(gapsize<=0)return IR_GAP_SIZE;
  return (int32)gapsize;
}


int FastaScaffold(FILE *out, SnapScaffoldMesg *scaff){
  int i,num_pairs=scaff->num_contig_pairs; 
  int32 contigid;
  int reverse;
  MultiAlignT *contig;
  SnapContigPairs *cp=scaff->contig_pairs;
  int32 *reversed;
  char *fseq;
  char *cseq;
  int running_length=0;
  char nchar;
  int scaffold_length,contig_length,flen,ngaps;

  nchar = 'n';

  if(cids == NULL){
    cids=CreateVA_int32(1);
    reversedVA=CreateVA_int32(num_pairs + 1);
    ctmp=CreateVA_char(200000);
    qtmp=CreateVA_char(200000);
    scaffold_sequence = CreateVA_char(200000);
  }

    ResetVA_int32(cids);
    ResetVA_int32(reversedVA);
    ResetVA_char(ctmp);
    ResetVA_char(qtmp);
    ResetVA_char(scaffold_sequence);

    reversed = Getint32(reversedVA,0);
  contigid = *LookupInUID2IIDHashTable_AS(uid2iid,cp[0].econtig1);
  // output label line for this scaffold
 fprintf(out,">" F_S64 " /type=%s\n",scaff->eaccession,fastaIdent);

  scaffold_length = 0;
  // calculate length of sequence:
  scaffold_length += GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,contigid));
  for ( i=0;i<num_pairs;i++ ) {
    scaffold_length += compute_gap(cp[i].mean);
    scaffold_length += GetMultiAlignUngappedLength(GetMultiAlignInStore(cstore,
									*LookupInUID2IIDHashTable_AS(uid2iid,cp[i].econtig2)));
  }
  EnableRangeVA_char(scaffold_sequence,scaffold_length+1);
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
  cseq = Getchar(ctmp,0);
  fseq = Getchar(scaffold_sequence,0);

  if(reversed[0]){
    Complement(cseq,contig_length);
  }
  running_length+=contig_length;
  if (running_length > GetNumchars(scaffold_sequence)) {
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
  }
  memcpy(fseq, Getchar(ctmp,0),contig_length);
  fseq+=contig_length;

  for ( i=0;i<num_pairs;i++ ) {
    ngaps = compute_gap(cp[i].mean);

     // append ngaps N's (or n's, depending on alternate_gap_spec)  to string as intercontig space
     ResetVA_char(ctmp);  
     EnableRangeVA_char(ctmp,ngaps+1);
     running_length+=ngaps;
     memset(Getchar(ctmp,0),nchar,ngaps);
     memset(fseq, nchar, ngaps);
     fseq+=ngaps;
     contigid =  *LookupInUID2IIDHashTable_AS(uid2iid,cp[i].econtig2);
     contig = GetMultiAlignInStore(cstore,contigid);
     contig_length = GetMultiAlignUngappedLength(contig);
     GetMultiAlignUngappedConsensus(contig,ctmp,qtmp);
     if(reversed[i+1]){
       cseq = Getchar(ctmp,0);
       Complement(cseq,contig_length);
     }
     running_length+=contig_length;
     if (running_length > GetNumchars(scaffold_sequence)) {
       fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
     }
     memcpy(fseq, Getchar(ctmp,0),contig_length);
     fseq +=contig_length;
  }
  // now, output the scaffold
  fseq = Getchar(scaffold_sequence,0);
  flen = strlen(fseq);
  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,fseq);
    fseq += 70;
  }
  return 1;
}

int FastaDegenerateScaffold(FILE *out, SnapDegenerateScaffoldMesg *scaff,VA_TYPE(int32) *is_placed)  {
  int32 contigid =  *LookupInUID2IIDHashTable_AS(uid2iid,scaff->econtig);
  int reverse;
  int i;
  MultiAlignT *contig;
  char *fseq;
  int running_length=0;
  int scaffold_length,contig_length,flen,ngaps;

  if(cids == NULL){
    cids=CreateVA_int32(1);
    reversedVA=CreateVA_int32(1000);
    ctmp=CreateVA_char(200000);
    qtmp=CreateVA_char(200000);
    scaffold_sequence = CreateVA_char(200000);
  }
  ResetVA_int32(cids);
  ResetVA_int32(reversedVA);
  ResetVA_char(ctmp);
  ResetVA_char(qtmp);
  ResetVA_char(scaffold_sequence);

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

  fseq = Getchar(scaffold_sequence,0);

  running_length+=contig_length;
  if (running_length > GetNumchars(scaffold_sequence)) {
    fprintf(stderr,"FastaScaffold warning: unexpectedly long string in scaffold\n");
  }

  // output label line for this scaffold
  fprintf(out,">" F_S64 " /type=%s_Degenerate\n",scaff->econtig,fastaIdent);

  memcpy(fseq, Getchar(ctmp,0),contig_length);

  flen = strlen(fseq);
  for (i = 0; i < flen; i += 70) { 
    fprintf(out,"%.*s\n",70,fseq);
    fseq += 70;
  }
  return 1;
}

int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  SnapScaffoldMesg *isf;
  MesgReader   reader;
  ScreenMatch *mat;
  SnapContigPairs *pairs;
  SnapConConMesg *contig;
  SnapUnitigMesg *unitig;
  SnapContigLinkMesg *link;
  MultiAlignT *ma;
  MultiAlignT *ma1;
  MultiAlignT *ma2;
  MultiPos *frag1;
  MultiPos *frag2;
  int num_frag1;
  int num_frag2;
  int ma1_len;
  int ma2_len;
  int scaffold=-1;
  int num_pairs;
  int i;
  int isplaced = 1;
  char *suffix;
  FILE *out = NULL;
  FILE *sublist = NULL;
  char buffer[256];
  char fastaNameBuffer[FILENAME_MAX];
  char fastaIdentifier[FILENAME_MAX];
  DIR *fastaDir;
  char *fastaDirName;
  char storeNameBuffer[FILENAME_MAX];
  VA_TYPE(char) *dummy_consensus;
  VA_TYPE(int32) *link_index;
  VA_TYPE(SnapContigLinkMesg) *clinks;
  VA_TYPE(int32) *is_placed;
  int32 placed=1;
  int32 unplaced=0;
  int32 linkid;
  int32 num_links;
  char ch;
  char *sublist_file;
  ID_Arrayp  tig_iids;
  ID_Arrayp  tig_iids_found;
  cds_int64  this_id;
  int do_all = 0;
  do_all = 1;
  cstore = CreateMultiAlignStoreT(0);
  bactig_unitigs = CreateVA_int32(0);
  link_index = CreateVA_int32(0);
  scaff_index = CreateVA_int32(0);
  is_placed = CreateVA_int32(0);
  clinks = CreateVA_SnapContigLinkMesg(0);
  global_coord = 0;
  optarg = NULL;
  fasta = 0;
  fastaDregs = 0;
  if ( argc < 2 ) {
     fprintf(stderr,"Try -h for usage\n");
     exit(1);
  }
  while ( ((ch = getopt(argc, argv, "h?df:")) != EOF)) {
        switch(ch) {
	case 'd':
	  fastaDregs = 1;
	  break;
        case 'f':
          fasta = 1;
          strcpy(fastaNameBuffer, optarg);
          HandleDir(fastaNameBuffer,fastaFileName);
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
          if(suffix!=NULL) *suffix = '\0'; // this cuts off the ext, so filename root can be used
          while ( (suffix = strchr(fastaIdent,'.')) != NULL ) {
             *suffix = '_';
          }
          break;
        case 'h':
        case '?':
          {
            fprintf(stderr,"\n\nUsage: asmProcessScaffolds_TER [-d] [-h] [-f fasta.output.filename [-s]]\n");
            fprintf(stderr,"\n The -f option produces a multi-fasta file, with one fasta record for each scaffold,\n");
            fprintf(stderr," and 'n's used as intercontig gap placeholders.  For gaps > 50 bp, the number of 'n's\n");
            fprintf(stderr," used represents the calculated mean gap size.  For others, 50 'n's are used.\n");
            fprintf(stderr,"\n The -d option produces a multi-fasta dregs file, with one fasta record for each degenerate scaffold,\n");
            fprintf(stderr,"\n The -d option is only meaningful when specified with -f\n");
            CleanExit(1);
          }
        default:
          {
            fprintf(stderr,"Invalid option -%c, try -h for usage\n",ch);
            CleanExit(1);
          }
        }
  }

  uid2iid = CreateUIDHashTable_AS(200000);

  if(fasta&&fastaDregs){
    sprintf(fastaDregsFileName,"%s.dregs",fastaFileName);
    fastaDregsFile = fopen(fastaDregsFileName,"w");
    if (fastaDregsFile == NULL ) {
      fprintf(stderr,"Failure to create fasta file %s\n", fastaDregsFileName);
      CleanExit(1);
    }
  }
  
  reader = InputFileType_AS( stdin );
  dummy_consensus = CreateVA_char(200000);

 while (reader(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_UTG)  {
      unitig = pmesg->m;
      InsertInUID2IIDHashTable_AS(uid2iid,unitig->eaccession,unitig->iaccession);
      /*
      printf("utg inserted %lld : %d, lookup gives %d\n",
	     unitig->eaccession, 
	     unitig->iaccession, 
	     *LookupInUID2IIDHashTable_AS(uid2iid,unitig->eaccession)); 
	     */
      SetVA_int32(is_placed, unitig->iaccession,  
       ( unitig->status == AS_SEP )? &placed : &unplaced );
    }
    if (pmesg->t ==MESG_CCO)  {
      contig = pmesg->m;
      InsertInUID2IIDHashTable_AS(uid2iid,contig->eaccession,contig->iaccession);
      //      assert(LookupInUID2IIDHashTable_AS(uid2iid,contig->eaccession)!=NULL); 
      /*
	printf("ctg inserted %lld : %d, lookup gives %d\n",
	contig->eaccession, 
	contig->iaccession, 
	*LookupInUID2IIDHashTable_AS(uid2iid,contig->eaccession)); 
	*/
      assert ( strlen(contig->consensus) == contig->length);
      ma = CreateMultiAlignTFromCCO(contig, -1,  0);
      SetMultiAlignInStore(cstore,ma->id,ma);
    }
    if (pmesg->t ==MESG_SCF)  {
      if ( fasta ) {
        FastaScaffold(fastaFile, (SnapScaffoldMesg *)pmesg->m);
        fflush(fastaFile);
      }
    }
    if (pmesg->t ==MESG_DSC)  {
      if ( fasta  && fastaDregs) {
        FastaDegenerateScaffold(fastaDregsFile, (SnapDegenerateScaffoldMesg *)pmesg->m,is_placed);
        fflush(fastaDregsFile);
      }
    }
 }
 if (fasta) {
   if(fastaDregs)
     fclose(fastaDregsFile);

   fclose(fastaFile);
 }
 exit (0);
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
