
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. Craig Venter Institute
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

static char const *rcsid = "$Id: AS_GKP_sff.c,v 1.19 2008-06-19 05:02:43 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_ALN_bruteforcedp.h"

//  Problems with this code:
//
//  1) It's possibly too aggressive when looking for linker.  It is
//  happy with 35 bases of alignment out of a 44bp linker, with 3
//  errors.
//
//  2) It needs to be done when OBT is done.  Ideally, after
//  merge-trimming, but before chimera.  Or, heck, during chimera.
//  But that creates a large headache -- we'd need to add new
//  fragments
//
//  To do that, we could abuse the clv to denote linker detected here,
//  but not acted on.
//
//
//  As it is now, most of the linker is detected, either by this or by
//  OBT.  On p.ging, half1, I only see
//
//  81530[218-0-0] 0[0-44] <25-0-100-forward-unknown>
//  edef=>E8YURXS01C8F4B,91131 mate=0,0 lib=LIBSFFE8YURXS01,1 clr=OBT,0,218 deleted=0
//  ddef=>linker
//  193-217 (20-44) <25-0-100>
//  gaattcaaaccctttcggttccaac
//  gaattcaaaccctttcggttccaac
//
//  after trimming.  This one was not detected because trimming got
//  rid of 30bp at the end of the read -- when gatekeeper was looking
//  for linker, it saw 25bp match, and 30bp more stuff on the other
//  side.  Not enough to split.


#undef  USE_454_TRIMMING

int    linkersLen = 1;
char  *linkers[1] = { "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC" };



typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32   magic_number;
  char     version[4];
  uint64   index_offset;
  uint32   index_length;
  uint32   number_of_reads;
  uint16   header_length;
  uint16   key_length;
  uint16   number_of_flows_per_read;
  uint8    flowgram_format_code;

  char    *flow_chars;        //  h->number_of_flows_per_read
  char    *key_sequence;      //  h->key_length

  void    *data_block;
  uint32   data_block_len;

  uint32   swap_endianess;

  //  Mate-finding DP related storage.
  //
  //  These really don't belong in here, but we know
  //  this is allocated.
  //
  char     alignA[AS_READ_MAX_LEN + AS_READ_MAX_LEN + 2];
  char     alignB[AS_READ_MAX_LEN + AS_READ_MAX_LEN + 2];
  dpCell   matrix[AS_READ_MAX_LEN][AS_READ_MAX_LEN];
} sffHeader;


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint32    magic_number;
  char      version[4];
  uint32    manifest_length;
  uint32    nothing;

  char     *manifest;
} sffManifest;


typedef struct {
  //  The next block is read in one swoop from the sff file.  DO NOT MODIFY!
  uint16   read_header_length;
  uint16   name_length;
  uint32   number_of_bases;
  uint16   clip_quality_left;
  uint16   clip_quality_right;
  uint16   clip_adapter_left;
  uint16   clip_adapter_right;

  char    *name;                 //  r->name_length
  uint16  *flowgram_values;      //  h->number_of_flows_per_read
  uint8   *flow_index_per_base;  //  r->number_of_bases
  char    *bases;                //  r->number_of_bases
  uint8   *quality_scores;       //  r->number_of_bases
  char    *quality;              //  quality_scores converted to CA-format qv

  int      final_length;         //  trimmed, processed read, ready for
  char    *final_bases;          //  loading.  NOT zero terminated.
  char    *final_quality;        //  DO NOT ZERO TERMINATE.

  void    *data_block;
  uint32   data_block_len;
} sffRead;


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x000000ff) | ((x <<  8) & 0x0000ff00);
  return(x);
}



static
void
readsff_manifest(FILE *sff, sffHeader *h, sffManifest *m) {

  if (h->index_length == 0)
    //  No manifest.
    return;

  if (AS_UTL_ftell(sff) != h->index_offset)
    //  Not at the manifest.
    return;

  if (m->manifest)
    //  Already got it?!
    return;

  AS_UTL_safeRead(sff, m, "readsff_manifest", sizeof(char), 16);

  if (h->swap_endianess) {
    m->magic_number    = uint32Swap(m->magic_number);
    m->manifest_length = uint32Swap(m->manifest_length);
  }

  m->manifest = (char *)safe_malloc(sizeof(char) * m->manifest_length + 1);
  AS_UTL_safeRead(sff, m->manifest, "readsff_manifest_text", sizeof(char), m->manifest_length);

  m->manifest[m->manifest_length] = 0;

  //  We only read the manifest.  There is still an index in there.

  uint64  padding_length = h->index_length - 16 - m->manifest_length;
  if (padding_length > 0) {
    //fprintf(stderr, "manifest pad "F_U64"\n", padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_manifest_pad", sizeof(char), padding_length);
    safe_free(junk);
  }
}


static
void
readsff_header(FILE *sff, sffHeader *h, sffManifest *m) {

  AS_UTL_safeRead(sff, h, "readsff_header_1", 31, 1);

  if (h->magic_number != 0x2e736666) {
    h->swap_endianess           = 1;
    h->magic_number             = uint32Swap(h->magic_number);
    h->index_offset             = uint64Swap(h->index_offset);
    h->index_length             = uint32Swap(h->index_length);
    h->number_of_reads          = uint32Swap(h->number_of_reads);
    h->header_length            = uint16Swap(h->header_length);
    h->key_length               = uint16Swap(h->key_length);
    h->number_of_flows_per_read = uint16Swap(h->number_of_flows_per_read);
  }

  assert(h->magic_number == 0x2e736666);

  uint32 newlen = h->number_of_flows_per_read + h->key_length + 2;
  if (h->data_block_len < newlen) {
    h->data_block_len = newlen;
    h->data_block     = safe_realloc(h->data_block, h->data_block_len);
  }

  memset(h->data_block, 0, h->data_block_len);

  h->flow_chars   = h->data_block;
  h->key_sequence = h->data_block + (h->number_of_flows_per_read + 1) * sizeof(char);

  AS_UTL_safeRead(sff,  h->flow_chars,   "readsff_header_2", sizeof(char), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff,  h->key_sequence, "readsff_header_3", sizeof(char), h->key_length);

  uint64  padding_length = h->header_length - 31 - h->number_of_flows_per_read - h->key_length;
  if (padding_length > 0) {
    //fprintf(stderr, "header pad "F_U64"\n", padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_header_4", sizeof(char), padding_length);
    safe_free(junk);
  }

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.
  //
  readsff_manifest(sff, h, m);
}


static
void
readsff_read(FILE *sff, sffHeader *h, sffRead *r) {

  AS_UTL_safeRead(sff, r, "readsff_read_1", 16, 1);

  if (h->swap_endianess) {
    r->read_header_length = uint16Swap(r->read_header_length);
    r->name_length        = uint16Swap(r->name_length);
    r->number_of_bases    = uint32Swap(r->number_of_bases);
    r->clip_quality_left  = uint16Swap(r->clip_quality_left);
    r->clip_quality_right = uint16Swap(r->clip_quality_right);
    r->clip_adapter_left  = uint16Swap(r->clip_adapter_left);
    r->clip_adapter_right = uint16Swap(r->clip_adapter_right);
  }

  //  Can you say UGLY?  Hey, it's a lot better than what I originally came up with.

  uint32 ss[6];
  ss[0] = (r->name_length + 1)          * sizeof(char);
  ss[1] = (h->number_of_flows_per_read) * sizeof(uint16) + ss[0];
  ss[2] = (r->number_of_bases)          * sizeof(uint8)  + ss[1];
  ss[3] = (r->number_of_bases + 1)      * sizeof(char)   + ss[2];
  ss[4] = (r->number_of_bases + 1)      * sizeof(uint8)  + ss[3];
  ss[5] = (r->number_of_bases + 1)      * sizeof(char)   + ss[4];

  if (r->data_block_len < ss[5]) {
    r->data_block_len = ss[5];
    r->data_block     = safe_realloc(r->data_block, r->data_block_len);
  }

  memset(r->data_block, 0, r->data_block_len);

  r->name                 = r->data_block;
  r->flowgram_values      = r->data_block + ss[0];
  r->flow_index_per_base  = r->data_block + ss[1];
  r->bases                = r->data_block + ss[2];
  r->quality_scores       = r->data_block + ss[3];
  r->quality              = r->data_block + ss[4];

  AS_UTL_safeRead(sff, r->name, "readsff_read_2", sizeof(char), r->name_length);
  r->name[r->name_length] = 0;

  uint64  padding_length = r->read_header_length - 16 - r->name_length;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 1 "F_U64"\n", padding_length);
    uint64  junk;
    AS_UTL_safeRead(sff, &junk, "readsff_read_3", sizeof(char), padding_length);
  }

  AS_UTL_safeRead(sff, r->flowgram_values,     "readsff_read_4", sizeof(uint16), h->number_of_flows_per_read);
  AS_UTL_safeRead(sff, r->flow_index_per_base, "readsff_read_5", sizeof(uint8),  r->number_of_bases);
  AS_UTL_safeRead(sff, r->bases,               "readsff_read_6", sizeof(char),   r->number_of_bases);
  AS_UTL_safeRead(sff, r->quality_scores,      "readsff_read_7", sizeof(uint8),  r->number_of_bases);

  int i;
  for (i=0; i<r->number_of_bases; i++)
    r->quality[i] = r->quality_scores[i] + '0';

  r->bases[r->number_of_bases] = 0;
  r->quality[r->number_of_bases] = 0;

  //  The padding_length is the number of bytes to make the above four
  //  chunks of data be of size that is divisible by 8.  The
  //  padding_length we compute directly below is the number of bytes
  //  we read past the last multiple of 8, and if that is non-zero, we
  //  need to read 8-padding_length bytes.
  //
  padding_length = (h->number_of_flows_per_read * sizeof(uint16) +
                    r->number_of_bases * sizeof(uint8) +
                    r->number_of_bases * sizeof(char) +
                    r->number_of_bases * sizeof(uint8)) % 8;
  if (padding_length > 0) {
    //fprintf(stderr, "read pad 2 "F_U64"\n", 8-padding_length);
    char *junk = (char *)safe_malloc(sizeof(char) * padding_length);
    AS_UTL_safeRead(sff, junk, "readsff_read_8", sizeof(char), 8 - padding_length);
    safe_free(junk);
  }
}



static
AS_UID
readsff_constructUIDFromName(char *name, int constructReadUID) {
  char  libname[32] = {0};

  assert(strlen(name) == 14);

  //  The read UID is just the name given by 454.  Thanks, 454!  The
  //  library UID is derived from the read UID by yanking out the
  //  position of the read, and replacing it with something obnoxious.
  //  We could have just truncated it to the last 9 letters, but
  //  that's kind of short.  Not that 15 is better.

  if (constructReadUID) {
    libname[0]  = name[0];
    libname[1]  = name[1];
    libname[2]  = name[2];
    libname[3]  = name[3];
    libname[4]  = name[4];
    libname[5]  = name[5];
    libname[6]  = name[6];
    libname[7]  = name[7];
    libname[8]  = name[8];
    libname[9]  = name[9];
    libname[10] = name[10];
    libname[11] = name[11];
    libname[12] = name[12];
    libname[13] = name[13];
    libname[14] = name[14];
    libname[15] = 0;
  } else {
    libname[0] = 'L';
    libname[1] = 'I';
    libname[2] = 'B';
    libname[3] = 'S';
    libname[4] = 'F';
    libname[5] = 'F';
    libname[6]  = name[0];
    libname[7]  = name[1];
    libname[8]  = name[2];
    libname[9]  = name[3];
    libname[10] = name[4];
    libname[11] = name[5];
    libname[12] = name[6];
    libname[13] = name[7];
    libname[14] = name[8];
    libname[15] = 0;
    libname[16] = 0;
  }

  return(AS_UID_load(libname));
}


static
AS_IID
readsff_constructLibraryIIDFromName(char *name) {
  AS_UID  uid = readsff_constructUIDFromName(name, 0);
  AS_IID  iid = getGatekeeperUIDtoIID(gkpStore, uid, NULL);

  if (iid == 0) {
    GateKeeperLibraryRecord  gkl;

    clearGateKeeperLibraryRecord(&gkl);

    gkl.libraryUID = uid;

    //  This crud is documented in AS_PER/AS_PER_gkpStore.h
    //  Zero is the default, we set to make it explicit

    gkl.spare2                      = 0;
    gkl.spare1                      = 0;

    gkl.discardReadsWithNs          = 1;
    gkl.doNotQVTrim                 = 1;
    gkl.goodBadQVThreshold          = 1;  //  Effectively, don't QV trim, redundant

    gkl.deletePerfectPrefixes       = 1;
    gkl.doNotTrustHomopolymerRuns   = 1;
    gkl.doNotOverlapTrim            = 0;
    gkl.isNotRandom                 = 0;

    gkl.hpsIsSomethingElse          = 0;
    gkl.hpsIsFlowGram               = 1;
    gkl.hpsIsPeakSpacing            = 0;

    gkl.orientation                 = AS_READ_ORIENT_UNKNOWN;

    gkl.mean                        = 0.0;
    gkl.stddev                      = 0.0;

    appendIndexStore(gkpStore->lib, &gkl);
    setGatekeeperUIDtoIID(gkpStore, gkl.libraryUID, getLastElemStore(gkpStore->lib), AS_IID_LIB);

    iid = getLastElemStore(gkpStore->lib);

    gkpStore->gkp.sffLibCreated++;
  }

  return(iid);
}



static
void
processRead(sffHeader *h,
            sffRead   *r, GateKeeperFragmentRecord *gkf) {

  clearGateKeeperFragmentRecord(gkf);

  gkf->readUID = readsff_constructUIDFromName(r->name, 1);

  //  Read already loaded?  Can't load again.  Set UID;s and IID's
  //  to zero to indicate this -- we'll catch it at the end.
  //
  if (getGatekeeperUIDtoIID(gkpStore, gkf->readUID, NULL)) {
    AS_GKP_reportError(AS_GKP_SFF_ALREADY_EXISTS,
                       AS_UID_toString(gkf->readUID));
    gkpStore->gkp.sffErrors++;

    gkf->readUID = AS_UID_undefined();
    return;
  }


  gkf->libraryIID  = readsff_constructLibraryIIDFromName(r->name);
  gkf->orientation = AS_READ_ORIENT_UNKNOWN;


  //  Set clear ranges
  //

#ifdef USE_454_TRIMMING
  //  Attempt to make sense of 454 supplied clear ranges.
  //
  //  These are base-based.  If either value is 0, that means the
  //  value was not computed.
  //
  //  We have a policy decision here.  If only one of the ranges is
  //  set, we can either ignore both, or set the unset one to the
  //  maximum.  We set it to the maximum.
  
  int  clq = r->clip_quality_left;
  int  crq = r->clip_quality_right;

  assert((r->clip_quality_left == 0) || (h->key_length <= r->clip_quality_left));
  assert((r->clip_adapter_left == 0) || (h->key_length <= r->clip_adapter_left));

  if (clq == 0)  clq = h->key_length + 1;
  if (crq == 0)  crq = r->number_of_bases;

  int  which;
  for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
    gkf->clearBeg[which] = clq - h->key_length - 1;
    gkf->clearEnd[which] = crq - h->key_length;
  }

  if ((r->clip_quality_left > 0) && (r->clip_quality_right > 0)) {
    gkf->hasQualityClear = 1;
    gkf->clearBeg[AS_READ_CLEAR_QLT] = r->clip_quality_left  - h->key_length - 1;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = r->clip_quality_right - h->key_length;
  } else if (r->clip_quality_left > 0) {
    gkf->hasQualityClear = 1;
    gkf->clearBeg[AS_READ_CLEAR_QLT] = r->clip_quality_left  - h->key_length - 1;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = r->number_of_bases - h->key_length;
  } else if (r->clip_quality_right > 0) {
    gkf->hasQualityClear = 1;
    gkf->clearBeg[AS_READ_CLEAR_QLT] = 0;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = r->clip_quality_right - h->key_length;
  } else {
    gkf->hasQualityClear = 0;
    gkf->clearBeg[AS_READ_CLEAR_QLT] = 0;
    gkf->clearEnd[AS_READ_CLEAR_QLT] = 0;
  }

  if ((r->clip_adapter_left > 0) && (r->clip_adapter_right > 0)) {
    gkf->hasVectorClear  = 1;
    gkf->clearBeg[AS_READ_CLEAR_VEC] = r->clip_adapter_left  - h->key_length - 1;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = r->clip_adapter_right - h->key_length;
  } else if (r->clip_adapter_left > 0) {
    gkf->hasVectorClear  = 1;
    gkf->clearBeg[AS_READ_CLEAR_VEC] = r->clip_adapter_left  - h->key_length - 1;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = r->number_of_bases - h->key_length;
  } else if (r->clip_adapter_right > 0) {
    gkf->hasVectorClear  = 1;
    gkf->clearBeg[AS_READ_CLEAR_VEC] = 0;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = r->clip_adapter_right - h->key_length;
  } else {
    gkf->hasVectorClear = 0;
    gkf->clearBeg[AS_READ_CLEAR_VEC] = 0;
    gkf->clearEnd[AS_READ_CLEAR_VEC] = 0;
  }
#else
  //  Don't trim at all.

  int  which;
  for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
    gkf->clearBeg[which] = 0;
    gkf->clearEnd[which] = r->number_of_bases - h->key_length;
  }
#endif


  //  Look for n's in the sequence.  This is a signature of an
  //  instrument problem.  Were we general, this would be
  //  discardReadsWithNs (for your grepping pleasure)
  {
    int  x = 0;

    for (x=0; x < r->number_of_bases; x++) {
      if ((r->bases[x] == 'n') || (r->bases[x] == 'N')) {

        //  Trim out the N?
        gkf->hasVectorClear = 1;
        gkf->clearBeg[AS_READ_CLEAR_VEC] = 0;
        gkf->clearEnd[AS_READ_CLEAR_VEC] = x - h->key_length;

        //  Nah, just delete it.
        gkf->deleted = 1;

        AS_GKP_reportError(AS_GKP_SFF_N, AS_UID_toString(gkf->readUID));
        break;
      }
    }
  }


  //  Too short?  Mark as deleted.
  //
  if (r->number_of_bases - h->key_length < AS_READ_MIN_LEN) {
    AS_GKP_reportError(AS_GKP_SFF_TOO_SHORT, r->name, r->number_of_bases, AS_READ_MIN_LEN);
    gkpStore->gkp.sffWarnings++;

    gkf->deleted = 1;
  }


  //  Too long?  Trim it.
  //
  if (r->number_of_bases - h->key_length > AS_READ_MAX_LEN) {
    AS_GKP_reportError(AS_GKP_SFF_TOO_LONG, r->name, r->number_of_bases, AS_READ_MAX_LEN);
    gkpStore->gkp.sffWarnings++;
      
    r->number_of_bases = AS_READ_MAX_LEN;

    r->bases  [AS_READ_MAX_LEN + h->key_length] = 0;
    r->quality[AS_READ_MAX_LEN + h->key_length] = 0;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      if (gkf->clearBeg[which] > AS_READ_MAX_LEN)
        gkf->clearBeg[which] = AS_READ_MAX_LEN;
      if (gkf->clearEnd[which] > AS_READ_MAX_LEN)
        gkf->clearEnd[which] = AS_READ_MAX_LEN;
    }
  }

  r->final_bases    = r->bases   + h->key_length;
  r->final_quality  = r->quality + h->key_length;
  r->final_length   = strlen(r->final_bases);
}


static
int
processMate(sffHeader *h,
            sffRead   *r,   GateKeeperFragmentRecord *gkf,
            sffRead   *rm1, GateKeeperFragmentRecord *gkm1,
            sffRead   *rm2, GateKeeperFragmentRecord *gkm2) {
  alignLinker_s  al = {0};
  int  linkerLength  = strlen(linkers[0]);

  //  Abort if this isn't a valid read -- processRead() might have already trashed it.
  //
  if ((gkf->deleted) || (AS_UID_isDefined(gkf->readUID) == 0))
    return(0);

  //  Otherwise, look for the mate linker.
  //
  alignLinker(h->alignA,
              h->alignB,
              linkers[0],
              r->final_bases,
              h->matrix,
              &al);


#ifdef USE_454_TRIMMING
#error mates do not support trimming with 454 trim points
#endif

#if 0
  int  lSize = al.begJ;
  int  rSize = al.lenB - al.endJ;
  fprintf(stderr, "fragment %d root=%s\n", getLastElemStore(gkpStore->frg) + 1, AS_UID_toString(gkf->readUID));
  fprintf(stderr, "alignLen=%d matches=%d lSize=%d rSize=%d\n", al.alignLen, al.matches, lSize, rSize);
  fprintf(stderr, "%3d-%3d %s\n", al.begI, al.endI, h->alignA);
  fprintf(stderr, "%3d-%3d %s\n", al.begJ, al.endJ, h->alignB);
#endif

  //  Did we find enough of the linker to do something?  We just need
  //  to throw out the obviously bad stuff.  When we get shorter and
  //  shorter, it's hard to define reasonable cutoffs.
  //
  //  Things that are called good here, but are actually bad, will be
  //  examined in OBT's chimera.  If there are no overlaps spanning,
  //  they'll be trimmed out, usually by being called chimeric.
  //
  int  goodAlignment = 0;
  int  bestAlignment = 0;

  if ((al.alignLen >=  5) && (al.matches + 1 >= al.alignLen))
    goodAlignment = 1;
  if ((al.alignLen >= 15) && (al.matches + 2 >= al.alignLen))
    goodAlignment = 1;
  if ((al.alignLen >= 30) && (al.matches + 3 >= al.alignLen))
    goodAlignment = 1;
  if ((al.alignLen >= 40) && (al.matches + 4 >= al.alignLen))
    goodAlignment = 1;

  if (goodAlignment == 0)
    return(0);

  if ((al.alignLen >= 42) && (al.matches + 2 >= al.alignLen))
    bestAlignment = 1;


  int  lSize = al.begJ;
  int  rSize = al.lenB - al.endJ;

  int  which;


  //  Adapter found on the left, but not enough to make a read.  Trim it out.
  //
  if ((bestAlignment) && (lSize < 64)) {
    r->final_length = rSize;

    r->final_bases   += al.endJ;
    r->final_quality += al.endJ;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      gkf->clearBeg[which] = 0;
      gkf->clearEnd[which] = r->final_length;
    }

    //  Recursively search for another copy of the linker.
    processMate(h, r, gkf, NULL, NULL, NULL, NULL);

    return(1);
  }


  //  Adapter found on the right, but not enough to make a read.  Trim it out.
  //
  if ((bestAlignment) && (rSize < 64)) {
    r->final_length = lSize;

    r->final_bases  [r->final_length] = 0;
    r->final_quality[r->final_length] = 0;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      gkf->clearBeg[which] = 0;
      gkf->clearEnd[which] = r->final_length;
    }

    //  Recursively search for another copy of the linker.
    processMate(h, r, gkf, NULL, NULL, NULL, NULL);

    return(1);
  }


  //  Adapter found in the middle, and enough to make two mated reads.
  //
  if ((bestAlignment) && (lSize >= 64) && (rSize >= 64)) {

    //  If we get here, and the two mate reads are null, we have
    //  found a second complete linker -- the original read is
    //  "seq-link-seq-link-seq" and we don't know where to break.
    //  The whole read is deleted in this case.
    //
    //  The mate is deleted later.
    //
    if ((rm1 == NULL) || (rm2 == NULL) || (gkm1 == NULL) || (gkm2 == NULL)) {
      gkf->deleted = 1;
      return(0);
    }

    memcpy(rm1, r, sizeof(sffRead));
    memcpy(rm2, r, sizeof(sffRead));

    memcpy(gkm1, gkf, sizeof(GateKeeperFragmentRecord));
    memcpy(gkm2, gkf, sizeof(GateKeeperFragmentRecord));

    //  1.  Make new UIDs for the two mated reads.  Nuke the old
    //  read.  Make the mates.
    //
    //  WARNING!  See that getLastElemStore() below?  It is forcing us to
    //  load the gkm1 read before the gkm2 read.
    {
      char  uid[64];
      strcpy(uid, AS_UID_toString(gkf->readUID));
      strcat(uid, "a");
      gkm1->readUID = AS_UID_load(uid);
      gkm1->readIID = getLastElemStore(gkpStore->frg) + 1;

      strcpy(uid, AS_UID_toString(gkf->readUID));
      strcat(uid, "b");
      gkm2->readUID = AS_UID_load(uid);
      gkm2->readIID = getLastElemStore(gkpStore->frg) + 2;

      gkf->readUID = AS_UID_undefined();

      gkm1->mateIID = gkm2->readIID;
      gkm2->mateIID = gkm1->readIID;

      gkm1->orientation = AS_READ_ORIENT_INNIE;
      gkm2->orientation = AS_READ_ORIENT_INNIE;
    }

    //need to enable more strict checking on things loaded -- no embedded nulls for exampe (is that valid in encoded?)

    //  2.  Construct new rm1, rm2.  Nuke the linker.  Reverse
    //  complement -- inplace -- the rm2 read.
    {
      int j;

      reverseComplement(r->final_bases, r->final_quality, lSize);

      rm1->final_length = lSize;

      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        gkm1->clearBeg[which] = 0;
        gkm1->clearEnd[which] = rm1->final_length;
      }

      for (j=al.begJ; j<al.endJ; j++) {
        r->final_bases[j]   = 0;
        r->final_quality[j] = 0;
      }

      //reverseComplement(r->final_bases + endJ, r->final_quality + endJ, rSize);

      rm2->final_length   = rSize;
      rm2->final_bases   += al.endJ;
      rm2->final_quality += al.endJ;

      for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
        gkm2->clearBeg[which] = 0;
        gkm2->clearEnd[which] = rm2->final_length;
      }
    }

    //  Recursively search for another copy of the linker.
    processMate(h, rm1, gkm1, NULL, NULL, NULL, NULL);
    processMate(h, rm2, gkm2, NULL, NULL, NULL, NULL);

    //  If either is deleted now, then we found another linker in
    //  the split read.
    //
    if ((gkm1->deleted) || (gkm2->deleted)) {
#warning sff errors not reported here
      gkm1->deleted = 1;
      gkm2->deleted = 1;
    }

    return(1);
  }  //  if match is in middle

  //  Significant alignment, but we didn't do anything with it.  Why?
  //  Overload some of the later clear ranges to convey this
  //  information to downstream processes.
  //
  //  We have 64 bits split into 4 16 bit words.  We want to store the
  //  position of the match in both the linker and the read, as well
  //  as length and number of matches.  Six things.  If we assume the
  //  linker is 255bp or smaller, we can pack.
  //
  //  This could benefit from having 64-bits of generic unioned data
  //  in the fragment record.

  gkf->sffLinkerDetectedButNotTrimmed = 1;

  gkf->hasVectorClear  = 0;
  gkf->hasQualityClear = 0;

  gkf->clearBeg[AS_READ_CLEAR_QLT] = (al.begI     << 8) | al.endI;     //  linker coords
  gkf->clearEnd[AS_READ_CLEAR_QLT] = (al.alignLen << 8) | al.matches;  //  quality

  gkf->clearBeg[AS_READ_CLEAR_VEC] = al.begJ;  //  read coords
  gkf->clearEnd[AS_READ_CLEAR_VEC] = al.endJ;

  if (0) {
    int  lSize = al.begJ;
    int  rSize = al.lenB - al.endJ;
    fprintf(stderr, "fragment %d root=%s\n", getLastElemStore(gkpStore->frg) + 1, AS_UID_toString(gkf->readUID));
    fprintf(stderr, "alignLen=%d matches=%d lSize=%d rSize=%d\n", al.alignLen, al.matches, lSize, rSize);
    fprintf(stderr, "%3d-%3d %s\n", al.begI, al.endI, h->alignA);
    fprintf(stderr, "%3d-%3d %s\n", al.begJ, al.endJ, h->alignB);
    fprintf(stderr, "%s\n", r->final_bases);
  }

  return(0);
}



static
void
addReadToStore(GateKeeperStore *gkp,
               sffRead *r,
               GateKeeperFragmentRecord *gkf) {
  char       encodedsequence[AS_FRAG_MAX_LEN+1] = {0};

  //  WARNING!  Search above for getLastElemStore() if you muck with this.
  gkf->readIID = getLastElemStore(gkpStore->frg) + 1;

  gkf->seqLen = r->final_length;
  gkf->hpsLen = 0;
  gkf->srcLen = 0;

  gkf->seqOffset = getLastElemStore(gkpStore->seq) + 1;
  gkf->qltOffset = getLastElemStore(gkpStore->qlt) + 1;
  gkf->hpsOffset = getLastElemStore(gkpStore->hps) + 1;
  gkf->srcOffset = getLastElemStore(gkpStore->src) + 1;

  setGatekeeperUIDtoIID(gkpStore, gkf->readUID, gkf->readIID, AS_IID_FRG);
  appendIndexStore(gkpStore->frg, gkf);

  appendStringStore(gkpStore->seq, r->final_bases, gkf->seqLen);

  encodeSequenceQuality(encodedsequence,
                        r->final_bases,
                        r->final_quality);
  appendStringStore(gkpStore->qlt, encodedsequence, gkf->seqLen);

  appendStringStore(gkpStore->hps, NULL, 0);
  appendStringStore(gkpStore->src, NULL, 0);

  gkpStore->gkp.sffLoaded++;
}



int
Load_SFF(FILE *sff, int searchForLinker) {

  sffHeader   *h   = (sffHeader   *)safe_calloc(sizeof(sffHeader),   1);
  sffManifest *m   = (sffManifest *)safe_calloc(sizeof(sffManifest), 1);
  sffRead     *r   = (sffRead     *)safe_calloc(sizeof(sffRead),     1);
  sffRead     *rm1 = (sffRead     *)safe_calloc(sizeof(sffRead),     1);
  sffRead     *rm2 = (sffRead     *)safe_calloc(sizeof(sffRead),     1);

  GateKeeperFragmentRecord *gkf  = (GateKeeperFragmentRecord *)safe_calloc(sizeof(GateKeeperFragmentRecord), 1);
  GateKeeperFragmentRecord *gkm1 = (GateKeeperFragmentRecord *)safe_calloc(sizeof(GateKeeperFragmentRecord), 1);
  GateKeeperFragmentRecord *gkm2 = (GateKeeperFragmentRecord *)safe_calloc(sizeof(GateKeeperFragmentRecord), 1);

  int         libraryUpdated = 0;

  readsff_header(sff, h, m);

  int         rn      = 0;

  //  Construct a gkpLibraryRecord for this sff file.  Well, this is
  //  where we'd LIKE to do it, but since the sff doesn't give us any
  //  reasonable way to make a UID from the header, we defer until we
  //  get the first read.  Then, we use the read timestamp, hash and
  //  region to make a library.

  gkpStore->gkp.sffInput += h->number_of_reads;

  //h->number_of_reads = 100;

  for (rn=0; rn < h->number_of_reads; rn++) {
    readsff_read(sff, h, r);

    gkf->readUID  = AS_UID_undefined();
    gkm1->readUID = AS_UID_undefined();
    gkm2->readUID = AS_UID_undefined();

    processRead(h, r, gkf);

    if (searchForLinker)
      processMate(h, r, gkf, rm1, gkm1, rm2, gkm2);

    if (AS_UID_isDefined(gkf->readUID)) {
      assert(AS_UID_isDefined(gkm1->readUID) == 0);
      assert(AS_UID_isDefined(gkm2->readUID) == 0);

      addReadToStore(gkpStore, r, gkf);
    }

    if (AS_UID_isDefined(gkm1->readUID) &&
        AS_UID_isDefined(gkm2->readUID)) {
      assert(AS_UID_isDefined(gkf->readUID) == 0);

      //  WARNING!  These two reads MUST be added in this order.
      //
      addReadToStore(gkpStore, rm1, gkm1);
      addReadToStore(gkpStore, rm2, gkm2);

      //  Update the library to refect that there are mates present.
      //  We guess at the insert size.
      //
      if (libraryUpdated == 0) {
        GateKeeperLibraryRecord  gkl;
        getIndexStore(gkpStore->lib, gkm1->libraryIID, &gkl);

        gkl.orientation                 = AS_READ_ORIENT_INNIE;
        gkl.mean                        = 3000.0;
        gkl.stddev                      =  300.0;

        setIndexStore(gkpStore->lib, gkm1->libraryIID, &gkl);

        libraryUpdated = 1;
      }
    }
  }

  //  Read the manifest?
  //
  if (m->manifest_length == 0) {
    //  We haven't read it yet.
    readsff_manifest(sff, h, m);
  }

  if (m->manifest != NULL) {
    //  Make sure that the reads have been rescored.

    if (strstr(m->manifest, "<qualityScoreVersion>1.1.03</qualityScoreVersion>") == NULL) {
      fprintf(stderr, "WARNING:  Fragments not rescored!\n");
    }
  }

  fprintf(stderr, "Added %d 454 reads.\n", rn);

  safe_free(h->data_block);
  safe_free(h);
  safe_free(m->manifest);
  safe_free(m);
  safe_free(r->data_block);
  safe_free(r);

  return(0);
}

