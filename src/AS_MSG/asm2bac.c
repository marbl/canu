
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
/* RCS info
 * $Id: asm2bac.c,v 1.1.1.1 2004-04-14 13:52:09 catmandew Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_TER_alloc.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"
#include "PrimitiveVA_MSG.h"
#include "MultiAlignStore_CNS.h"

static FILE *gbatFile;
static char *gbatName;
VA_DEF(FragMesg)
VA_DEF(BactigMesg)
VA_DEF(IntUnitigLinkMesg)
VA_DEF(IntContigLinkMesg)
MultiAlignStoreT *cstore;
MultiAlignStoreT *ustore;

void CleanExit(int rc) {
  char command[100+FILENAME_MAX];
  if( gbatFile != NULL ){
    fclose(gbatFile);
    sprintf(command,"rm -f %s",gbatName);
    fprintf(stderr,"%s\n",command);
    system(command);
  }
  exit(rc);
}

CDS_UID_t grab_uid(int use_real_ids)  {
  // First, grab some UIDS from the server.
  static CDS_UID_t interval_UID[4];
  static int uid_blockSize=300;
  static int uid_init=0;
  int32  uidStatus;
  CDS_UID_t NEW_ID;
  if ( ! uid_init ) {
    uidStatus = get_uids(uid_blockSize,interval_UID,use_real_ids); // keep fake for initial testing
    if ( uidStatus != UID_CODE_OK ) {
      fprintf(stderr,"Abnormal return from get_uids (SYS). Return code %d\n",uidStatus);
      CleanExit(1);
    }
    uid_init=1;
  }
  uidStatus = get_next_uid(&NEW_ID, use_real_ids);
  if ( uidStatus != UID_CODE_OK ) {
      uidStatus = get_uids(uid_blockSize,interval_UID,use_real_ids);
      uidStatus = get_next_uid(&NEW_ID,use_real_ids);
      if ( uidStatus != UID_CODE_OK ) {
        fprintf(stderr,"Abnormal return from get_next_uid (SYS). Return code %d\n",uidStatus);
        CleanExit(1);
      }
  }
  return NEW_ID;
}

/* These are globals, though for clarity are passed in read_map arglist */
VA_TYPE(CDS_UID_t) *FRAG_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *UNITIG_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *CONTIG_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *SCAFF_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *DSCAFF_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *DIST_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *BACTIG_ID_MAP=NULL;
VA_TYPE(CDS_UID_t) *BAC_ID_MAP=NULL;

int read_map(FILE *mapFile){
  char lbuf[1024];
  int which_map=0;
  CDS_IID_t iid;
  CDS_UID_t uid;
  VA_TYPE(CDS_UID_t) *mptr;
  FRAG_ID_MAP = CreateVA_CDS_UID_t(1); // 1
  UNITIG_ID_MAP = CreateVA_CDS_UID_t(1); // 2 
  CONTIG_ID_MAP = CreateVA_CDS_UID_t(1); // 3
  SCAFF_ID_MAP = CreateVA_CDS_UID_t(1); // 4
  BACTIG_ID_MAP = CreateVA_CDS_UID_t(1); // 5 (contains the iids for the new Celera-only bactigs)
  DIST_ID_MAP = CreateVA_CDS_UID_t(1);  // 6
  // Screen map is 7
  // Repeat item map is 8
  DSCAFF_ID_MAP = CreateVA_CDS_UID_t(1); // 9
  // Note that DSCAFF is an exception since iids are of the contained contig
  //   the IntDegenerateScaffold itself 
  BAC_ID_MAP = CreateVA_CDS_UID_t(1); // 10 (contains the iids for the new Celera-only BACs)
  if ( mapFile == NULL ) {
     return 0;
  }
  while ( fgets(lbuf,1024,mapFile) ) {
    if ( isalpha(lbuf[0]) ) {
      if ( strncmp(lbuf,"Fragment",8) == 0 ) {
        which_map=1; 
      } else if ( strncmp(lbuf,"Unitig",6) == 0 ) {
        which_map=2; 
      } else if ( strncmp(lbuf,"Contig",6) == 0 ) {
        which_map=3; 
      } else if ( strncmp(lbuf,"Scaffold",8) == 0 ) {
        which_map=4; 
      } else if ( strncmp(lbuf,"Bactig",6) == 0 ) {
        which_map=5; 
      } else if ( strncmp(lbuf,"Distrib",7) == 0 ) {
        which_map=6; 
      } else if ( strncmp(lbuf,"Screen",6) == 0 ) {
        which_map=7; 
      } else if ( strncmp(lbuf,"Repeat",6) == 0 ) {
        which_map=8; 
      } else if ( strncmp(lbuf,"Degen",5) == 0 ) {
        which_map=9; 
      }
    } else {
     // is a iid-pair line, read  and place into proper map array
      if(sscanf(lbuf,F_IID " " F_UID,&iid,&uid) != 2)
        assert(0);
     switch (which_map) {
     case 1:
       mptr = FRAG_ID_MAP;
       break;
     case 2:
       mptr = UNITIG_ID_MAP;
       break;
     case 3:
       mptr = CONTIG_ID_MAP;
       break;
     case 4:
       mptr = SCAFF_ID_MAP;
       break;
     case 5:
       mptr = BACTIG_ID_MAP;
       break;
     case 6:
       mptr = DIST_ID_MAP;
       break;
     case 9:
       mptr = DSCAFF_ID_MAP;
       break;
     case 0:
     case 7:
     case 8:
     default:
       mptr = NULL;
       break;
     }
     if (mptr != NULL ) SetCDS_UID_t(mptr,iid,&uid);
    }
  }
  return 1;
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


static void Complement(char *in, int len);


int DumpScaffoldBac(FILE *batFile, IntScaffoldMesg *scaff, int length_cutoff, int REAL_UIDS, char *comment)  {
/*  Generate a Bac message (with accompanying DST record, Frag records,
    and dump to the batch File  */

//  For this scaffold, we need to generate:
//    a DST record (for length)
//    a BAC record (containing num_contigs BTG records
//    num_contigs FRG records

  int i,num_pairs=scaff->num_contig_pairs; 
  int num_contigs=num_pairs+1;
  int32 contigid;
  DistanceMesg dist; // will need uid
  BacMesg bac;  //  will need bacuid (use scaffold uid) AND sequid (get new)
  BactigMesg btgtmp;
  FragMesg frgtmp;
  GenericMesg pmesg;
  static VA_TYPE(int32) *cids=NULL;
  static VA_TYPE(FragMesg) *frgs=NULL;
  static VA_TYPE(BactigMesg) *btgs=NULL;
  static VA_TYPE(char) *ungappedSeq=NULL;
  static VA_TYPE(char) *ungappedQlt=NULL;
  IntContigPairs *cp=scaff->contig_pairs;
  time_t etm=time(0);

  if ( cids == NULL ) {
    cids=CreateVA_int32(5000);
    frgs=CreateVA_FragMesg(5000);
    btgs=CreateVA_BactigMesg(5000);
    ungappedSeq=CreateVA_char(200000);
    ungappedQlt=CreateVA_char(200000);
  } else {
    ResetVA_int32(cids);
    ResetVA_FragMesg(frgs);
    ResetVA_BactigMesg(btgs);
    ResetVA_char(ungappedSeq);
    ResetVA_char(ungappedQlt);
  }
  dist.mean = 0;
  dist.stddev = 0;
  for (i=0;i<num_contigs;i++) {
      int contigid = (num_contigs==1 || i < num_contigs-1)?cp[i].contig1:cp[i-1].contig2;
      MultiAlignT *contig=GetMultiAlignInStore(cstore,contigid);
      dist.mean += GetMultiAlignUngappedLength(contig);
      if ( i < num_contigs-1 ) {
        dist.mean += cp[i].mean;
        dist.stddev += cp[i].stddev;
      }
  }

  // Don't bother if not long enough to meet cutoff
  if ( dist.mean < length_cutoff ) return 0;

  dist.action = AS_ADD;
  dist.eaccession = grab_uid(REAL_UIDS);
  bac.action = AS_ADD;
  bac.ebac_id = *GetCDS_UID_t(SCAFF_ID_MAP,scaff->iaccession);
  bac.type = AS_UNFINISHED;
  bac.source = comment;
  bac.eseq_id = grab_uid(REAL_UIDS);
  bac.entry_time = etm;
  bac.elength = dist.eaccession;
  bac.num_bactigs = num_contigs;
  contigid = cp[0].contig1;
  frgtmp.action = AS_ADD;
  frgtmp.type = AS_BACTIG;
  frgtmp.elocale = bac.ebac_id;
  frgtmp.eseq_id = bac.eseq_id;
  frgtmp.entry_time = etm;
  frgtmp.source = comment;
  frgtmp.clear_rng.bgn = 0;


  for (i=0;i<num_contigs;i++) {
      // get contig multialignment out of cstore
      int contigid = (num_contigs==1 || i < num_contigs-1)?cp[i].contig1:cp[i-1].contig2;
      MultiAlignT *contig=GetMultiAlignInStore(cstore,contigid);
      btgtmp.eaccession = grab_uid(REAL_UIDS);
      btgtmp.length = GetMultiAlignUngappedLength(contig);
      AppendVA_int32(cids,&contigid);
      ResetVA_char(ungappedSeq);
      ResetVA_char(ungappedQlt);
      GetMultiAlignUngappedConsensus(contig, ungappedSeq, ungappedQlt);
      frgtmp.sequence = (char *) malloc(GetNumchars(ungappedSeq));
      memcpy(frgtmp.sequence,Getchar(ungappedSeq,0),GetNumchars(ungappedSeq));
      frgtmp.quality = (char *) malloc(GetNumchars(ungappedQlt));
      memcpy(frgtmp.quality,Getchar(ungappedQlt,0),GetNumchars(ungappedQlt));
      frgtmp.eaccession = btgtmp.eaccession;
      frgtmp.ebactig_id = btgtmp.eaccession;
      frgtmp.clear_rng.end = btgtmp.length;
      AppendVA_FragMesg(frgs,&frgtmp);
      AppendVA_BactigMesg(btgs,&btgtmp);
  }

  bac.bactig_list = GetBactigMesg(btgs,0);
  // Now, flush Dist/Bac/Frags to file
  assert(batFile != NULL);
  pmesg.t = MESG_DST;
  pmesg.m = &dist;
  WriteProtoMesg_AS(batFile,&pmesg);
  pmesg.t = MESG_BAC;
  pmesg.m = &bac;
  WriteProtoMesg_AS(batFile,&pmesg);
  pmesg.t = MESG_FRG;
  for (i=0;i<num_contigs;i++) {
    pmesg.m = GetFragMesg(frgs,i);
    WriteProtoMesg_AS(batFile,&pmesg);
  }
  fflush(batFile);
  for (i=0;i<num_contigs;i++) {
     FragMesg *f = GetFragMesg(frgs,i);
     free(f->sequence);
     free(f->quality);
  }
  return 1;
}

int DumpScaffoldCNS(IntScaffoldMesg *scaff, int length_cutoff)  {
/*  Generate a consensus file representing only this scaffold */

//  For this scaffold, we need to generate:
//    an IUM record for each contig
//    an ICM record for each contig
//    an ISF record 

  static VA_TYPE(int32) *cids=NULL;
  static VA_TYPE(FragMesg) *frgs=NULL;
  static VA_TYPE(BactigMesg) *btgs=NULL;
  static VA_TYPE(char) *ungappedSeq=NULL;
  static VA_TYPE(char) *ungappedQlt=NULL;

  if ( cids == NULL ) {
    cids=CreateVA_int32(5000);
    frgs=CreateVA_FragMesg(5000);
    btgs=CreateVA_BactigMesg(5000);
    ungappedSeq=CreateVA_char(200000);
    ungappedQlt=CreateVA_char(200000);
  } else {
    ResetVA_int32(cids);
    ResetVA_FragMesg(frgs);
    ResetVA_BactigMesg(btgs);
    ResetVA_char(ungappedSeq);
    ResetVA_char(ungappedQlt);
  }
  return 1;
}

int DumpDScaffoldBac(FILE *batFile, IntDegenerateScaffoldMesg *scaff, int length_cutoff, int REAL_UIDS)  {
//  For this scaffold, we need to generate:
//    a DST record (for length)
//    a BAC record (containing 1 BTG record)
//    a single FRG 
  int num_contigs=1;
  int32 contigid;
  DistanceMesg dist; // will need uid
  BacMesg bac;  //  will need bacuid (use scaffold uid) AND sequid (get new)
  BactigMesg btgtmp;
  FragMesg frgtmp;
  GenericMesg pmesg;
  static VA_TYPE(char) *ungappedSeq=NULL;
  static VA_TYPE(char) *ungappedQlt=NULL;
  time_t etm=time(0);

  if ( ungappedSeq == NULL ) {
    ungappedSeq=CreateVA_char(200000);
    ungappedQlt=CreateVA_char(200000);
  } else {
    ResetVA_char(ungappedSeq);
    ResetVA_char(ungappedQlt);
  }


  dist.mean = 0;
  dist.stddev = 0;
  dist.action = AS_ADD;
  dist.eaccession = grab_uid(REAL_UIDS);
  bac.action = AS_ADD;
  bac.ebac_id = *GetCDS_UID_t(DSCAFF_ID_MAP,scaff->icontig); // this is scaffold id for dscaff containing icontig
  bac.type = AS_UNFINISHED;
  bac.eseq_id = grab_uid(REAL_UIDS);
  bac.entry_time = etm;
  bac.elength = dist.eaccession;
  bac.num_bactigs = num_contigs;
  bac.source = NULL;
  contigid = scaff->icontig;
  frgtmp.action = AS_ADD;
  frgtmp.type = AS_BACTIG;
  frgtmp.elocale = bac.ebac_id;
  frgtmp.eseq_id = bac.eseq_id;
  frgtmp.entry_time = etm;
  frgtmp.source = NULL;
  frgtmp.clear_rng.bgn = 0;

  {
   MultiAlignT *contig=GetMultiAlignInStore(cstore,scaff->icontig);
   btgtmp.eaccession = *GetCDS_UID_t(CONTIG_ID_MAP,contigid);
   btgtmp.length = GetMultiAlignUngappedLength(contig);
   dist.mean += btgtmp.length;
   GetMultiAlignUngappedConsensus(contig, ungappedSeq, ungappedQlt);
   frgtmp.sequence = Getchar(ungappedSeq,0);
   frgtmp.quality = Getchar(ungappedQlt,0);
   frgtmp.eaccession = grab_uid(REAL_UIDS);
   frgtmp.clear_rng.end = btgtmp.length;
   bac.bactig_list = &btgtmp;
  }
  // Now, flush Dist/Bac/Frags to file
  if ( dist.mean > length_cutoff ) {
    assert(batFile != NULL);
    pmesg.t = MESG_DST;
    pmesg.m = &dist;
    WriteProtoMesg_AS(batFile,&pmesg);
    pmesg.t = MESG_BAC;
    pmesg.m = &bac;
    WriteProtoMesg_AS(batFile,&pmesg);
    pmesg.t = MESG_FRG;
    pmesg.m = &frgtmp;
    WriteProtoMesg_AS(batFile,&pmesg);
    fflush(batFile);
  }
  DeleteVA_char(ungappedSeq);
  DeleteVA_char(ungappedQlt);
  return 1;
}


int main(int argc, char *argv[])
{ GenericMesg *pmesg;
  GenericMesg bmesg;
  BatchMesg batch;
  MesgReader   reader;
  MultiAlignT *ma;
  char *comment=NULL;
  char nameBuffer[1000];
  char batchBuffer[1000];
  char commentBuffer[1000];
  char batchComment[]="Created from interally generated Celera-only assembly file.";
  FILE *mapFile = NULL;
  FILE *batFile = NULL; // this file will be keyed to UID of Batch
  int REAL_UIDS=0;
  int output_cns=0;
  int batches = 0;
  int bacs_in_batch = 0;
  int length_cutoff=5000;
  CDS_UID_t first_faux_uid=1000000;
  int ch;
  VA_TYPE(IntUnitigLinkMesg) *ulinks;
  VA_TYPE(IntContigLinkMesg) *clinks;
  VA_TYPE(char) *dummy_consensus;
  ustore = CreateMultiAlignStoreT(0);
  cstore = CreateMultiAlignStoreT(0);
  ulinks = CreateVA_IntUnitigLinkMesg(0);
  clinks = CreateVA_IntContigLinkMesg(0);
  optarg = NULL;
  if ( argc < 2 ) {
     fprintf(stderr,"Try -h for usage\n");
     exit(1);
  }
  while ( ((ch = getopt(argc, argv, "h?m:l:ui:c:")) != EOF)) {
        switch(ch) {
        case 'm':
          mapFile=fopen(optarg,"r"); 
          break;
        case 'l':
          length_cutoff=atoi(optarg); 
          break;
        case 'u':
          REAL_UIDS=1;
          break;
        case 'c':
          comment=optarg;
          break;
        case 'i':
          if ( ! REAL_UIDS ) {
            first_faux_uid=atoi(optarg);
            break;
          } else {
            fprintf(stderr,"Options -u and -i are incompatible\n");
          }
        case 'h':
        case '?':
          {
            fprintf(stderr,"\n\nasm2bac splits the contents of an asm file by scaffolds into a collection of BACs and CGW files\n");
            fprintf(stderr,"\n\nUsage: asm2bac -m map_file [-h] [-l scaffold_length_cutoff] [-u] [-i uid]< input_file \n");
            fprintf(stderr,"\n The -m flag specifies the iid->uid mapping for the assembly file\n");
            fprintf(stderr,"\n The -l option sets the cutoff value for scaffold length (default: 5 KB )\n");
            fprintf(stderr,"\n The -i option sets the first faux uid to use if not using real uids (default: 1000000)\n");
            fprintf(stderr,"\n    that is, only scaffolds > scaffold_length_cutoff will generate a BAC message and cgw file\n");
            fprintf(stderr,"\n The -u option indicates that REAL Celera UIDs should be used (should be used except for debugging)\n");
            CleanExit(1);
          }
        default:
          {
            fprintf(stderr,"Invalid option -%c, try -h for usage\n",ch);
            CleanExit(1);
          }
        }
  }
  reader = InputFileType_AS( stdin );
  dummy_consensus = CreateVA_char(500000);

  if ( mapFile != NULL ) {
     // read in mapping
     read_map(mapFile);
  } // else, maps will be NULL and identity mapping will be assumed
  if ( REAL_UIDS ) {
      check_environment();
  } else {
      set_start_uid(first_faux_uid);
  }

  if ( ! output_cns ) {

     // start up by creating a BAT message for this assembly
     // and creating a directory to hold  the component cgw files.
     //    Note:  generate the cgw files with natural order iid, and replace with iids from
     //           the overlay assembler gatekeeper when they become available.
     sprintf(nameBuffer,"Batch from Celera-only assembly");
     if ( comment != NULL ) {
       sprintf(commentBuffer,"%s\n%s",batchComment,comment);
     } else {
       sprintf(commentBuffer,"%s",batchComment);
     }
     batch.name = nameBuffer;
     batch.created = time(0);
     batch.iaccession = batches++;
     batch.eaccession = grab_uid(REAL_UIDS);
     batch.comment = batchComment;
     sprintf(batchBuffer,"Celera_Only_" F_UID ".cai",batch.eaccession);
     batFile = fopen(batchBuffer,"w");
     gbatFile = batFile;
     gbatName = batchBuffer;
     bacs_in_batch=0;
     bmesg.t = MESG_BAT;
     bmesg.m = &batch;
     WriteProtoMesg_AS(batFile,&bmesg);
 
  } else {
    // this just ensures that the output directory is available
     sprintf(nameBuffer,"Celera_UBAC-cns/output.cns");
     HandleDir(nameBuffer, batchBuffer);
  }

  while (reader(stdin,&pmesg) != EOF){
    if (pmesg->t ==MESG_IUM)  {
     if (0) {
       IntUnitigMesg  *unitig = (IntUnitigMesg *) pmesg->m;
      if ( strlen(unitig->consensus) != unitig->length) {
         char *cptr;
         ResetVA_char(dummy_consensus);
         EnableRangeVA_char(dummy_consensus,unitig->length+1);
         cptr = Getchar(dummy_consensus,0); 
         memset(cptr,'N',unitig->length);
         unitig->consensus = cptr;
         unitig->quality = cptr;
      }
      ma = CreateMultiAlignTFromIUM(unitig, -1,  0);
      SetMultiAlignInStore(ustore,ma->id,ma);
     }
    }
    if (pmesg->t ==MESG_IUL)  {
      //ulink = pmesg->m;
      // cast link to IntLink
      //linkid = GetNumint32s(ulinks);
      //SetVA_int32(ulinks, ulink->unitig1, &linkid);
      //AppendVA_IntUnitigLinkMesg(ulinks,ulink);
    }
    if (pmesg->t ==MESG_ICM)  {
      IntConConMesg  *contig = (IntConConMesg *) pmesg->m;
      if (contig->placed == AS_PLACED) {
        ma = CreateMultiAlignTFromICM(contig, -1,  0);
        SetMultiAlignInStore(cstore,ma->id,ma);
      }
    }
    if (pmesg->t ==MESG_ICL) {
      //clink = pmesg->m;
      // cast link to IntLink
      //linkid = GetNumint32s(clinks);
      //SetVA_int32(clinks, clink->contig1, &linkid);
      //AppendVA_IntContigLinkMesg(clinks,clink);
    }
    if (pmesg->t ==MESG_ISF) {
      if ( bacs_in_batch > 5000 )  {
        // create a new batch 
        batch.created = time(0);
        batch.iaccession = batches++;
        batch.eaccession = grab_uid(REAL_UIDS);
        sprintf(nameBuffer,"Batch " F_UID " from Celera-only assembly",batch.eaccession);
        batch.name = nameBuffer;
        sprintf(batchBuffer,"Celera_Only_" F_UID ".cai",batch.eaccession);
        batFile = fopen(batchBuffer,"w");
        gbatFile = batFile;
        bacs_in_batch=0;
        bmesg.t = MESG_BAT;
        bmesg.m = &batch;
        WriteProtoMesg_AS(batFile,&bmesg);
      }
      bacs_in_batch+= DumpScaffoldBac(batFile, (IntScaffoldMesg *)pmesg->m, length_cutoff, REAL_UIDS,comment);
      fflush(batFile);
    }
    if (pmesg->t ==MESG_IDS) {
      //DumpDScaffoldBac(batFile, (IntDegenerateScaffoldMesg *)pmesg->m, length_cutoff, REAL_UIDS);
      //fflush(batFile);
    }
 }
 fclose(batFile);
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
  static int Firstime = 1;

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
	  assert(WCinvert[(int )*t] != '?' &&
		 WCinvert[(int) *s] != '?');

	  c = *s;
          *s++ = WCinvert[(int) *t];
          *t-- = WCinvert[c];
        }
      if (s == t)
        *s = WCinvert[(int) *s];
    }

  }
}
