
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007-2008, J. CRaig Venter Institute.
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

const char *mainid = "$Id: sffToCA.c,v 1.48 2010-06-17 22:33:03 jasonmiller9704 Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_reverseComplement.h"
#include "AS_UTL/AS_UTL_fasta.h"
#include "AS_GKP_include.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_ALN_bruteforcedp.h"

//  For the exact-prefix dedup to work, a fragment must be larger than
//  the DEDUP_SPAN (valid values are 48 and 64).  After the dedup, we
//  search for mates, and those reads must only be larger than the
//  assembler minimum (which could be 30bp).
//
#define DEDUP_SPAN            48
#define FRAG_MIN_LEN          MAX(AS_READ_MIN_LEN, DEDUP_SPAN)
#define MATE_MIN_LEN              AS_READ_MIN_LEN


#define CLEAR_ALL        0x00
#define CLEAR_454        0x01
#define CLEAR_N          0x02
#define CLEAR_PAIR_N     0x04
#define CLEAR_DISCARD_N  0x08
#define CLEAR_ERRR       0xff

#define TRIM_NONE        0
#define TRIM_SOFT        1
#define TRIM_HARD        2
#define TRIM_CHOP        3
#define TRIM_ERRR        9

char *TRIM_NAMES[4]  = { "none", "soft", "hard", "chop" };

#define AS_LINKER_MAX_SEQS       50
#define AS_LINKER_CUSTOM_OFFSET  24

gkStore         *gkpStore    = NULL;
FILE            *logFile     = NULL;
uint32           clearAction = CLEAR_454;
uint32           clearSet    = 0;
uint32           trimAction  = TRIM_HARD;


typedef struct {
  uint32   readsInSFF;

  //  Length status:  Should add to numReadsInSFF.
  //
  uint32   lenTooShort;
  uint32   lenOK;
  uint32   lenTrimmedByN;
  uint32   lenTooLong;

  //  Linker status:  If we search for linker, these should add to numReadsInSFF.
  //
  uint32   notExaminedForLinker;  //  Already deleted (dup, short or N)
  uint32   noLinker;              //  No linker detected
  uint32   badLinker;             //  Inconsistent linker detected
  uint32   partialLinker;         //  Some linker detected, passed to OBT
  uint32   fullLinker;            //  Good linker

  //  Final status:  Should add to numReadsInSFF.
  //
  uint32   fragmentsOutput;
  uint32   matesOutput;
  uint32   deletedDuplicates;
  uint32   deletedTooShort;
  uint32   deletedByN;
} statistics;

statistics st = {0};

static
void
writeStatistics(char **argv, int argc, int firstFileArg, char *fragName, int haveLinker, char **linker, int *search, char *stsName) {

  errno = 0;
  FILE *statOut = fopen(stsName, "w");
  if (errno) {
    fprintf(stderr, "ERROR: Failed to open the stats file '%s': %s\n", stsName, strerror(errno));
    statOut = stderr;
  }

  fprintf(statOut, "PARAMETERS\n");
  while (firstFileArg < argc)
    fprintf(statOut, "input sff               %s\n", argv[firstFileArg++]);

  fprintf(statOut, "output fragments        %s\n", fragName);

  if (clearAction == CLEAR_ALL)
    fprintf(statOut, "clear range             all\n");
  if (clearAction & CLEAR_454)
    fprintf(statOut, "clear range             454\n");
  if (clearAction & CLEAR_N)
    fprintf(statOut, "clear range             n\n");
  if (clearAction & CLEAR_PAIR_N)
    fprintf(statOut, "clear range             pair-of-n\n");
  if (clearAction & CLEAR_DISCARD_N)
    fprintf(statOut, "clear range             discard-n\n");

  fprintf(statOut, "trimming                %s\n", TRIM_NAMES[trimAction]);

  if (search[0] == TRUE)
    fprintf(statOut, "linker                  %s (FLX)\n", linker[0]);

  if (search[1] == TRUE)
    fprintf(statOut, "linker                  %s (Titanium)\n", linker[1]);

  for (int32 linkerID=3; linkerID < AS_LINKER_MAX_SEQS; linkerID++)
    if (search[linkerID] == TRUE)
      fprintf(statOut, "linker                  %s\n", linker[linkerID]);

  fprintf(statOut, "\n");

  fprintf(statOut, "INPUT\n");
  fprintf(statOut, "numReadsInSFF           "F_U32"\n", st.readsInSFF);
  fprintf(statOut, "\n");
  fprintf(statOut, "LENGTH\n");
  fprintf(statOut, "too short               "F_U32"\n", st.lenTooShort);
  fprintf(statOut, "ok                      "F_U32"\n", st.lenOK);
  fprintf(statOut, "trimmed by N            "F_U32"\n", st.lenTrimmedByN);
  fprintf(statOut, "too long                "F_U32"\n", st.lenTooLong);
  fprintf(statOut, "                        -------\n");
  fprintf(statOut, "                        "F_U32"\n", st.lenTooShort + st.lenOK + st.lenTrimmedByN + st.lenTooLong);
  fprintf(statOut, "\n");

  if (haveLinker) {
    fprintf(statOut, "LINKER\n");
    fprintf(statOut, "not examined            "F_U32"\n", st.notExaminedForLinker);
    fprintf(statOut, "none detected           "F_U32"\n", st.noLinker);
    fprintf(statOut, "inconsistent            "F_U32"\n", st.badLinker);
    fprintf(statOut, "partial                 "F_U32"\n", st.partialLinker);
    fprintf(statOut, "good                    "F_U32"\n", st.fullLinker);
    fprintf(statOut, "                        -------\n");
    fprintf(statOut, "                        "F_U32"\n", st.notExaminedForLinker + st.noLinker + st.badLinker + st.partialLinker + st.fullLinker);
    fprintf(statOut, "\n");
  }

  fprintf(statOut, "OUTCOME\n");
  fprintf(statOut, "fragment                "F_U32"\n", st.fragmentsOutput);
  fprintf(statOut, "mate pair               "F_U32"\n", st.matesOutput);
  fprintf(statOut, "deleted inconsistent    "F_U32"\n", st.badLinker);
  fprintf(statOut, "deleted duplicate       "F_U32"\n", st.deletedDuplicates);
  fprintf(statOut, "deleted too short       "F_U32"\n", st.deletedTooShort);
  fprintf(statOut, "deleted N not allowed   "F_U32"\n", st.deletedByN);
  fprintf(statOut, "                        -------\n");
  fprintf(statOut, "                        "F_U32"\n", st.fragmentsOutput + st.matesOutput + st.badLinker + st.deletedDuplicates + st.deletedTooShort + st.deletedByN);

  if (statOut != stderr)
    fclose(statOut);

  assert(st.readsInSFF == st.lenTooShort + st.lenOK + st.lenTrimmedByN + st.lenTooLong);

  if (haveLinker)
    assert(st.readsInSFF == st.notExaminedForLinker + st.noLinker + st.badLinker + st.partialLinker + st.fullLinker);

  assert(st.readsInSFF == st.fragmentsOutput + st.matesOutput + st.badLinker + st.deletedDuplicates + st.deletedTooShort + st.deletedByN);
}



////////////////////////////////////////////////////////////////////////////////
//
//  Reads an SFF file, inserts all the reads (that are of the proper length)
//  into a gkpStore.
//

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

  uint8   *data_block;
  uint32   data_block_len;

  uint32   swap_endianess;
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

  uint8   *data_block;
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
    h->data_block     = (uint8 *)safe_realloc(h->data_block, sizeof(uint8) * h->data_block_len);
  }

  memset(h->data_block, 0, h->data_block_len);

  h->flow_chars   = (char *)h->data_block;
  h->key_sequence = (char *)h->data_block + (h->number_of_flows_per_read + 1) * sizeof(char);

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
    r->data_block     = (uint8 *)safe_realloc(r->data_block, sizeof(uint8) * r->data_block_len);
  }

  memset(r->data_block, 0, r->data_block_len);

  r->name                 = (char   *)(r->data_block);
  r->flowgram_values      = (uint16 *)(r->data_block + ss[0]);
  r->flow_index_per_base  = (uint8  *)(r->data_block + ss[1]);
  r->bases                = (char   *)(r->data_block + ss[2]);
  r->quality_scores       = (uint8  *)(r->data_block + ss[3]);
  r->quality              = (char   *)(r->data_block + ss[4]);

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
    char *junk = (char *)safe_malloc(sizeof(char) * (8 - padding_length));
    AS_UTL_safeRead(sff, junk, "readsff_read_8", sizeof(char), 8 - padding_length);
    safe_free(junk);
  }
}

// Process Read.
//
// return:
// The ruturn value indicates whether read should be added to the store.
// It is ok to delete a read and retain it to the store:
// this routine would set the deleted flag to 1 and return true.
// Current policy is we do not add deleted reads to the store.
//
// parameters:
// Pass in pointers to the input file header (h),
// the populated input read record (r),
// and the gatekeeper fragment record to be populated (fr). 
static
int
processRead(sffHeader *h,
            sffRead   *r, gkFragment *fr) {
  AS_UID  readUID;

  st.readsInSFF++;

  readUID = AS_UID_load(r->name);

  //  Read already loaded?  Can't load again.  Set UID;s and IID's
  //  to zero to indicate this -- we'll catch it at the end.
  //
  if (gkpStore->gkStore_getUIDtoIID(readUID, NULL)) {
    fprintf(stderr, "Read '%s' already exists.  Duplicate deleted.\n",
            AS_UID_toString(readUID));
    return(false);
  }

  ////////////////////////////////////////
  //
  //  Chop off any N's at the end of the read.  Titanium likes to do
  //  this to us.
  //
  while ((r->bases[r->number_of_bases-1] == 'N') ||
         (r->bases[r->number_of_bases-1] == 'n')) {
    r->number_of_bases--;
    r->bases[r->number_of_bases] = 0;
    r->quality[r->number_of_bases] = 0;
  }

  if (r->clip_adapter_right > r->number_of_bases)
    r->clip_adapter_right = r->number_of_bases;

  if (r->clip_quality_right > r->number_of_bases)
    r->clip_quality_right = r->number_of_bases;

  ////////////////////////////////////////
  //
  //  Check that the read is length is OK; some corrupt files will truncate sequence or quality.
  //  This is just a quick test for that case.  Using strlen here is overkill; when the fragment is
  //  added to the store, lengths are checked again.
  //
  if (r->bases[r->number_of_bases-1] == 0)
    fprintf(stderr, "ERROR:  Read '%s' sequence is truncated.  Corrupt file?\n",
            AS_UID_toString(readUID));
  if (r->quality[r->number_of_bases-1] == 0)
    fprintf(stderr, "ERROR:  Read '%s' quality values are truncated.  Corrupt file?\n",
            AS_UID_toString(readUID));
  assert(r->bases[r->number_of_bases-1] != 0);
  assert(r->quality[r->number_of_bases-1] != 0);

  ////////////////////////////////////////
  //
  //  Attempt to make sense of 454 supplied clear ranges.
  //
  //  These are base-based.  If either value is 0, that means the
  //  value was not computed.  In that case, we set it to the extent
  //  (max or min).
  //
  int clq = h->key_length;
  int crq = r->number_of_bases;

  int cla = h->key_length;
  int cra = r->number_of_bases;

  if (clearAction & CLEAR_454) {
    //  Left point should be zero or after the key
    assert((r->clip_quality_left == 0) || (h->key_length <= r->clip_quality_left));
    assert((r->clip_adapter_left == 0) || (h->key_length <= r->clip_adapter_left));

    //  Right point should be zero or before the end
    assert((r->clip_quality_right == 0) || (r->clip_quality_left <= r->number_of_bases));
    assert((r->clip_adapter_right == 0) || (r->clip_adapter_left <= r->number_of_bases));

    clq = MAX(r->clip_quality_left,  h->key_length + 1) - 1;
    cla = MAX(r->clip_adapter_left,  h->key_length + 1) - 1;

    crq = (r->clip_quality_right > 0) ? r->clip_quality_right : r->number_of_bases;
    cra = (r->clip_adapter_right > 0) ? r->clip_adapter_right : r->number_of_bases;
  }

  ////////////////////////////////////////
  //
  //  Find the CLEAR_N and CLEAR_PAIR_N points.  If we're allowing the
  //  use of the 454 clear ranges, don't check for N's before that.
  //  qThis assumes that clq and cla are NOT set if the 454 clear is
  //  NOT used.
  //
  int  cln = h->key_length;
  int  crn = r->number_of_bases;
  int  frn = r->number_of_bases;  //  first-n

  int  isTrimN = 0;  //  Remember if we changed the clear range

  if ((clearAction & CLEAR_N) ||
      (clearAction & CLEAR_DISCARD_N)) {
    int  f  = 0;
    int  b  = MAX(clq, cla);
    int  e  = r->number_of_bases;
    char *s = r->bases;

    for (f=b; f<e; f++)
      if ((s[f] == 'n') || (s[f] == 'N')) {
        isTrimN = 1;
        break;
      }

    if (clearAction & CLEAR_N) {
      cln = b;
      crn = f;
    }

    frn = f;
  }

  int  clp = h->key_length;
  int  crp = r->number_of_bases;

  if (clearAction & CLEAR_PAIR_N) {
    int   f = 0;
    int   b = MAX(clq, cla);
    int   e = r->number_of_bases - 1;
    char *s = r->bases;

    for (f=b; f<e; f++)
      if (((s[f+0] == 'n') || (s[f+0] == 'N')) &&
          ((s[f+1] == 'n') || (s[f+1] == 'N'))) {
        isTrimN = 1;
        break;
      }

    clp = b;
    crp = f;
  }


  ////////////////////////////////////////
  //
  //  Make sense of all these by blindly intersecting them together.
  //
  int clf = MAX(MAX(clq, cla), MAX(cln, clp));
  int crf = MIN(MIN(crq, cra), MIN(crn, crp));


  ////////////////////////////////////////
  //
  //  Now, decide how to set the clear ranges.
  //

  fr->maxBgn = 1;  //  No max yet.
  fr->maxEnd = 0;

  fr->vecBgn = 1;  //  There is no vector clear defined for 454 reads.
  fr->vecEnd = 0;

  fr->tntBgn = 1;  //  Nothing contaminated.
  fr->tntEnd = 0;

  switch (trimAction) {
    case TRIM_NONE:
      //  Set clear range to the whole untrimmed read.
      fr->clrBgn = h->key_length;
      fr->clrEnd = r->number_of_bases;
      break;

    case TRIM_SOFT:
      //  Set clear ranges to whatever we discovered above, but do not
      //  limit OBT to anything.
      fr->clrBgn = clf;
      fr->clrEnd = crf;
      break;

    case TRIM_HARD:
      //  Set clear ranges to whatever we discovered above, and limit
      //  OBT to those ranges.
      fr->clrBgn = fr->maxBgn = clf;
      fr->clrEnd = fr->maxEnd = crf;
      break;

    case TRIM_CHOP:
      //  Rewrite the read to remove the non-clear sequence.  We keep
      //  in the usually four base long key at the start (with some
      //  amount of pain since clf,clr include those four bases).
      memmove(r->bases   + h->key_length, r->bases   + clf, sizeof(char) * (crf - clf));
      memmove(r->quality + h->key_length, r->quality + clf, sizeof(char) * (crf - clf));

      r->number_of_bases = h->key_length + crf - clf;

      r->bases  [h->key_length + crf - clf] = 0;
      r->quality[h->key_length + crf - clf] = 0;

      fr->clrBgn = h->key_length;
      fr->clrEnd = h->key_length + crf - clf;
      break;
    default:
      break;
  }


  if ((r->clip_quality_left > r->clip_quality_right) ||
      (r->number_of_bases - h->key_length < FRAG_MIN_LEN)) {
    //  Reads too short will never be of any use, and they're not loaded.
    //
    //  The first test catches reads that 454 decided are completely trash.
    //
    //  The second makes sure that the bases are long enough, leaving it up to OBT to decide if
    //  they're any good.
    //
    st.lenTooShort++;
    st.deletedTooShort++;
    st.notExaminedForLinker++;  //  because this SFF read isn't even added to the store

    fprintf(logFile, "Read '%s' of length %d clear %d,%d is too short.  Read deleted.\n",
            r->name,
            r->number_of_bases - h->key_length,
            fr->clrBgn - h->key_length,
            fr->clrEnd - h->key_length);

    return(false);

  } else if (r->number_of_bases - h->key_length <= AS_READ_MAX_NORMAL_LEN) {
    //  Read is just right.
    if (isTrimN)
      st.lenTrimmedByN++;
    else
      st.lenOK++;

  } else {
    //  Reads too long can be loaded into the store, but until we fix overlaps,
    //  we cannot use them.  Truncate.
    //
    st.lenTooLong++;

    fprintf(logFile, "Read '%s' of length %d is too long.  Truncating to %d bases.\n",
            r->name, r->number_of_bases - h->key_length, AS_READ_MAX_NORMAL_LEN);

    r->number_of_bases = AS_READ_MAX_NORMAL_LEN + h->key_length;

    r->bases  [AS_READ_MAX_NORMAL_LEN + h->key_length] = 0;
    r->quality[AS_READ_MAX_NORMAL_LEN + h->key_length] = 0;

    if (fr->clrBgn > AS_READ_MAX_NORMAL_LEN - h->key_length)   fr->clrBgn = AS_READ_MAX_NORMAL_LEN + h->key_length;
    if (fr->clrEnd > AS_READ_MAX_NORMAL_LEN - h->key_length)   fr->clrEnd = AS_READ_MAX_NORMAL_LEN + h->key_length;
  }


  ////////////////////////////////////////
  //
  //  If told to, and there is still an N in the sequence, trash the
  //  whole thing.
  //
  if ((clearAction & CLEAR_DISCARD_N) && (frn < crf)) {
    st.deletedByN++;
    st.notExaminedForLinker++;  //  because this SFF read isn't even added to the store

    fprintf(logFile, "Read '%s' contains an N at position %d.  Read deleted.\n",
            AS_UID_toString(readUID), crn);

    return(false);
  }


  //  Finally, adjust everything to remove the key_length bases from the start.
  //
  fr->clrBgn -= h->key_length;
  fr->clrEnd -= h->key_length;

  if (fr->maxBgn < fr->maxEnd) {
    fr->maxBgn -= h->key_length;
    fr->maxEnd -= h->key_length;
  }

  r->final_bases    = r->bases   + h->key_length;
  r->final_quality  = r->quality + h->key_length;
  r->final_length   = r->number_of_bases - h->key_length;

  fr->gkFragment_setType(GKFRAGMENT_NORMAL);

  //  Construct a UID from the 454 read name
  fr->gkFragment_setReadUID(readUID);
  fr->gkFragment_setIsDeleted(0);

  fr->gkFragment_setLibraryIID(1);
  fr->gkFragment_setOrientation(AS_READ_ORIENT_UNKNOWN);

  //  Copy sequence to the gkFragment
  //
  memcpy(fr->gkFragment_getSequence(), r->final_bases,   sizeof(char) * (r->final_length + 1));
  memcpy(fr->gkFragment_getQuality(),  r->final_quality, sizeof(char) * (r->final_length + 1));

  fr->gkFragment_setLength(r->final_length);

  fr->gkFragment_getSequence()[r->final_length + 1] = 0;
  fr->gkFragment_getQuality() [r->final_length + 1] = 0;

  //  Check clear ranges.  Why allow equality?  The two N trimming options can set the clear range
  //  to 0,0 if there is an N in the first position.
  assert(fr->clrBgn <= fr->clrEnd);
  assert(fr->vecBgn >= fr->vecEnd);
  assert(fr->tntBgn >= fr->tntEnd);

  return(true);
} // processRead



int
loadSFF(char *sffName) {
  FILE                      *sff  = NULL;
  int                        fic  = 0;
  sffHeader                  h    = {0};
  sffManifest                m    = {0};
  sffRead                    r    = {0};
  gkFragment                 fr;
  int                        rn   = 0;

  fr.gkFragment_enableGatekeeperMode(gkpStore);

  errno = 0;

  fprintf(stderr, "loadSFF()-- Loading '%s'.\n", sffName);

  if        (strcasecmp(sffName + strlen(sffName) - 3, ".gz") == 0) {
    char  cmd[1024];
    sprintf(cmd, "gzip -dc %s", sffName);
    sff   = popen(cmd, "r");
    fic   = 1;
  } else if (strcasecmp(sffName + strlen(sffName) - 4, ".bz2") == 0) {
    char  cmd[1024];
    sprintf(cmd, "bzip2 -dc %s", sffName);
    sff   = popen(cmd, "r");
    fic   = 1;
  } else {
    sff   = fopen(sffName, "r");
    fic   = 0;
  }
  if (errno)
    fprintf(stderr, "ERROR!  Failed to open '%s': %s\n", sffName, strerror(errno)), exit(1);


  readsff_header(sff, &h, &m);

  for (rn=0; rn < h.number_of_reads; rn++) {
    readsff_read(sff, &h, &r);
    if (processRead(&h, &r, &fr))
      gkpStore->gkStore_addFragment(&fr);
  }

  //  Read the manifest if we haven't already done so.
  if (m.manifest_length == 0)
    readsff_manifest(sff, &h, &m);

  //  Make sure that the reads have been rescored.
  if (m.manifest != NULL)
    if (strstr(m.manifest, "<qualityScoreVersion>1.1.03</qualityScoreVersion>") == NULL)
      fprintf(stderr, "WARNING:  Fragments not rescored!\n");

  errno = 0;

  if (fic)
    pclose(sff);
  else
    fclose(sff);
  if (errno)
    fprintf(stderr, "WARNING!  Failed to close '%s': %s\n", sffName, strerror(errno));

  safe_free(h.data_block);
  safe_free(m.manifest);
  safe_free(r.data_block);

  return(0);
} 





////////////////////////////////////////////////////////////////////////////////
//
//  Removes all reads that are a perfect prefix of some other read.
//
//  The algorithm builds a 64-bit value from the first N bases, sorts
//  the hashes, then examines any clique of hash collisions for
//  perfect prefixes.

typedef struct {
  uint64    hash;
  uint32    iid;
} fragHash;

static
int
fragHashCompare(const void *a, const void *b) {
  fragHash const *A = (fragHash const *)a;
  fragHash const *B = (fragHash const *)b;

  if (A->hash < B->hash) return(-1);
  if (A->hash > B->hash) return( 1);
  return(0);
}

void
removeDuplicateReads(void) {
  uint32        fragsLen = 0;
  uint32        fragsMax = 0;
  gkFragment   *frags    = NULL;
  gkFragment    fr;

  fragHash   *fh    = new fragHash [gkpStore->gkStore_getNumFragments() + 1];
  uint32      fhLen = 0;

  uint64 map[256] = { 0 };
  uint32 s, h, n;

  for (s=0; s<256; s++)
    map[s] = 0;
  map['A'] = map['a'] = 0x00;
  map['C'] = map['c'] = 0x01;
  map['G'] = map['g'] = 0x02;
  map['T'] = map['t'] = 0x03;

  fprintf(stderr, "removeDuplicateReads()-- from %d to %d\n", 1, gkpStore->gkStore_getNumFragments() + 1);

  fr.gkFragment_enableGatekeeperMode(gkpStore);

  for (int32 thisElem=1; thisElem<=gkpStore->gkStore_getNumFragments(); thisElem++) {
    gkpStore->gkStore_getFragment(thisElem, &fr, GKFRAGMENT_SEQ);

    char *seq1      = fr.gkFragment_getSequence();

    uint32 seqLen   = fr.gkFragment_getSequenceLength();
    uint64 hash     = 0;

    assert(seqLen >= DEDUP_SPAN);

#if DEDUP_SPAN == 48
    //  Our "hash" is just the spaced seed "101" (repeating).  It
    //  covers the first 48 bases, picking out 32.
    //
    for (s=0, n=0; n<16; n++) {
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
    }
#elif DEDUP_SPAN == 64
    //  Our "hash" is just the spaced seed "1010" (repeating).  It
    //  covers the first 64 bases, picking out 32.
    //
    for (s=0, n=0; n<16; n++) {
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
      hash <<= 2;
      hash  |= map[seq1[s]];
      s++;
      s++;
    }
#else
#error invalid DEDUP_SPAN must be 48 or 64
#endif

    fh[fhLen].hash     = hash;
    fh[fhLen].iid      = thisElem;

    fhLen++;
  }

  qsort(fh, fhLen, sizeof(fragHash), fragHashCompare);

  uint32  beg = 0;
  uint32  end = 0;

  while (beg < fhLen) {

    //  We DO need to examine the whole clique (pairwise).  We cannot
    //  simply sort by size, because if we get three frags of the same
    //  size, it could be that #1 is a prefix of #3, and #2 is just of
    //  the same size.  Even there, we'd need to examine all pairs.

    end = beg + 1;

    //  First, find a pair of adjacent matches
    //
    while ((end < fhLen) &&
           (fh[beg].hash != fh[end].hash)) {
      beg++;
      end++;
    }

    //  Got a match?
    //
    if (end < fhLen) {
      uint32 b, e;

      //  Advance end to the end of the matches
      //
      while ((fh[beg].hash == fh[end].hash) && (end < fhLen))
        end++;

      //  Yeah, we could extend scope of this test to include the for
      //  loops, but those will stop quick enough.

      if ((beg + 1 < end) && (end-beg > 1000))
        fprintf(stderr, "Large potential duplicate set from "F_U32" to "F_U32" ("F_U32" things)\n", beg, end, end - beg);

      //  Load the fragments
      //
      if (end - beg > fragsMax) {
        delete [] frags;
        fragsMax = end - beg + 512;
        frags    = new gkFragment [fragsMax];
      }

      for (b=beg; b<end; b++)
        gkpStore->gkStore_getFragment(fh[b].iid, &frags[b-beg], GKFRAGMENT_SEQ);

      //  Compare all-vs-all in the range
      //
      for (b=beg; b<end; b++) {
        for (e=b+1; e<end; e++) {

          AS_IID     iid1 = fh[b].iid;
          AS_IID     iid2 = fh[e].iid;

          gkFragment  *fr1 = frags + b - beg;
          gkFragment  *fr2 = frags + e - beg;

          assert(iid1 == fr1->gkFragment_getReadIID());
          assert(iid2 == fr2->gkFragment_getReadIID());

          uint32 del1 = fr1->gkFragment_getIsDeleted();
          uint32 del2 = fr2->gkFragment_getIsDeleted();

          uint32 len1 = fr1->gkFragment_getSequenceLength();
          uint32 len2 = fr2->gkFragment_getSequenceLength();

          if ((del1) && (len1 < len2))
            continue;
          if ((del2) && (len2 < len1))
            continue;

          if (len1 == len2) {
            if ((del1) && (iid1 < iid2))
              continue;
            if ((del2) && (iid2 < iid1))
              continue;
          }

          if (del1 && del2)
            continue;

          char *seq1 = fr1->gkFragment_getSequence();
          char *seq2 = fr2->gkFragment_getSequence();

          uint32 len = MIN(len1, len2);

          if (strncmp(seq1, seq2, len) == 0) {

            //  A real collision.  Delete smaller of the two (either
            //  smaller sequence length or smaller iid).  We can skip
            //  the delete if it's already deleted.

            AS_UID     deletedUID = AS_UID_undefined();
            AS_IID     deletedIID = 0;
            uint32     deleted    = 0;

            if ((len == fr1->gkFragment_getSequenceLength()) &&
                (len == fr2->gkFragment_getSequenceLength())) {
              deletedIID = (iid1 < iid2) ? iid1 : iid2;
              deletedUID = (iid1 < iid2) ? fr1->gkFragment_getReadUID() : fr2->gkFragment_getReadUID();
              deleted    = (iid1 < iid2) ? del1 : del2;
            } else if (len == fr1->gkFragment_getSequenceLength()) {
              deletedIID = iid1;
              deletedUID = fr1->gkFragment_getReadUID();
              deleted    = del1;
            } else {
              deletedIID = iid2;
              deletedUID = fr2->gkFragment_getReadUID();
              deleted    = del2;
            }

            //  If we need to delete something, delete it, then update
            //  our cached copy.  We still need the sequence, as an
            //  even shorter fragment can be deleted by the one we
            //  just deleted.

            if (deleted == 0) {
              st.deletedDuplicates++;

              fprintf(logFile, "Delete read %s,%d a prefix of %s,%d\n",
                      AS_UID_toString(deletedUID), deletedIID,
                      (deletedIID == iid1) ? AS_UID_toString(fr2->gkFragment_getReadUID()) : AS_UID_toString(fr1->gkFragment_getReadUID()),
                      (deletedIID == iid1) ? iid2 : iid1);

              gkpStore->gkStore_delFragment(deletedIID);
              gkpStore->gkStore_getFragment(deletedIID, (deletedIID == iid1) ? fr1 : fr2, GKFRAGMENT_SEQ);
            }
          }
        }
      }
    }

    beg = end;
  }

  delete [] fh;
  delete [] frags;

  fprintf(stderr, "removeDuplicateReads()-- finished\n");
}


////////////////////////////////////////////////////////////////////////////////
//
//  For a given gkFragment, scan the sequence for a linker.  If found,
//  generate two new mated reads and delete the original read.
//
//

typedef struct {
  char     h_alignA[AS_READ_MAX_NORMAL_LEN + AS_READ_MAX_NORMAL_LEN + 2];
  char     h_alignB[AS_READ_MAX_NORMAL_LEN + AS_READ_MAX_NORMAL_LEN + 2];
  dpCell   h_matrix[AS_READ_MAX_NORMAL_LEN + 1][AS_READ_MAX_NORMAL_LEN + 1];
} dpMatrix;

dpMatrix  *globalMatrix = NULL;


static
int
processMate(gkFragment *fr,
            gkFragment *m1,
            gkFragment *m2,
            char        *linker[AS_LINKER_MAX_SEQS],
            int          search[AS_LINKER_MAX_SEQS],
	    const int   stringent) {

  alignLinker_s  al = {0};

  //  Did we find enough of the linker to do something?  We just need
  //  to throw out the obviously bad stuff.  When we get shorter and
  //  shorter, it's hard to define reasonable cutoffs.
  //
  //  Things that are called good here, but are actually bad, will be
  //  examined in OBT's chimera.  If there are no overlaps spanning,
  //  they'll be trimmed out, usually by being called chimeric.
  //
  //  int  goodAlignment = 0;
  int  bestAlignment = 0;
  int  functionalAlignment = 0;
  int  fractionalAlignment = 0;
  int  minimalAlignment = 0;
  int  foundAlignment = 0;
  int  linkerID      = 0;
  int  mismatches = 0;
  const int LINKER_POSITIVE = 1;
  const int LINKER_NEGATIVE = 0;
  int  allowedToSplit = 0;
  uint32 lSize = 0;
  uint32 rSize = 0;
  int another1, another2;

  int linkerLength = 0;
  const int MAX_MISMATCH_CONSIDERED_FUNCTIONAL = 2;
  const int DIFFERENCE_CONSIDERED_MINIMAL = 15;
  const int DIFFERENCE_CONSIDERED_FRACTIONAL = 25;

  assert(fr->clrBgn < fr->clrEnd);
  if ((m1 == NULL) && (m2 == NULL)) {
    allowedToSplit = 0;
  } else {
    allowedToSplit = 1;
  }
  assert (stringent==0 || stringent==1);

  //  Linker array contains multiple linkers or the forward and reverse of one linker.
  //  Loop tests each linker in the array.
  //  For each linker, find the best alignment to the given read.
  //  On first alignment found, exit the loop.
  //  The later code may recurse, so it might find other linkers in the same read.

  //  Known problem with our loop structure.
  //  We act on the first linker with an alignment, not the best.
  //  There is probably some loss of sensitivity.
  //  However, a read with multiple linker hits is suspect anyway.

  while (linkerID < AS_LINKER_MAX_SEQS && foundAlignment == 0) {
    if (search[linkerID] == TRUE) {
      assert(linker[linkerID] != NULL);
      linkerLength = strlen ( linker[linkerID] );

      char *seq      = fr->gkFragment_getSequence();
      char  stopBase = seq[fr->clrEnd];

      seq[fr->clrEnd] = 0;

      alignLinker(globalMatrix->h_alignA,
                  globalMatrix->h_alignB,
                  linker[linkerID],
                  seq + fr->clrBgn,
                  globalMatrix->h_matrix,
                  &al,
                  FALSE, FALSE, 0, 0);

      seq[fr->clrEnd] = stopBase;

      lSize = al.begJ;
      rSize = al.lenB - al.endJ;
      al.begJ += fr->clrBgn;
      al.endJ += fr->clrBgn;
   
      assert(lSize >= 0);
      assert(lSize <= fr->gkFragment_getSequenceLength());
      assert(rSize >= 0);
      assert(rSize <= fr->gkFragment_getSequenceLength());

      // Minimal =>    Mark linker sequence as possible contaminant.
      // Fractional => Remove the linker sequence.
      // Functional => Try to split the read into 2 mates.
      mismatches = al.alignLen - al.matches;
      if (mismatches <= 5) {
	if (1==stringent) {
	  // Look for a full-length alignment.
	  if ((al.alignLen >=  linkerLength - MAX_MISMATCH_CONSIDERED_FUNCTIONAL ) 
	      && (mismatches <= MAX_MISMATCH_CONSIDERED_FUNCTIONAL )) {
	    minimalAlignment = 1;
	    fractionalAlignment = 1; 
	    functionalAlignment = 1;
	  } 
	} else {
	  // Look for a partial alignment.
	  functionalAlignment = 0;
	  if (al.matches - mismatches > DIFFERENCE_CONSIDERED_MINIMAL) {
	    minimalAlignment = 1;
	    if (al.matches - mismatches > DIFFERENCE_CONSIDERED_FRACTIONAL) {
	      fractionalAlignment = 1; 
	    }
	  }
	}
      }
      foundAlignment = // Exit the loop after we find an alingment.
	minimalAlignment+fractionalAlignment+functionalAlignment;
      if (foundAlignment>0) {
	fprintf(logFile, 
		"ProcessMate(%s,%d) found %d,%d,%d.\n",
		AS_UID_toString(fr->gkFragment_getReadUID()), stringent,
		minimalAlignment,fractionalAlignment,functionalAlignment); 
	fprintf(logFile, "D fragment %d clr=%d,%d root=%s\n", 
		gkpStore->gkStore_getNumFragments() + 1, 
		fr->clrBgn, fr->clrEnd, 
		AS_UID_toString(fr->gkFragment_getReadUID()));
	fprintf(logFile, "D alignLen=%d matches=%d lSize=%d rSize=%d\n", 
		al.alignLen, al.matches, lSize, rSize);
	fprintf(logFile, "D I %3d-%3d %s\n", al.begI, al.endI, globalMatrix->h_alignA);
	fprintf(logFile, "D J %3d-%3d %s\n", al.begJ, al.endJ, globalMatrix->h_alignB);
      } // log
      
    } // search linker
    linkerID++; 
  } // while
  
  if ((0==allowedToSplit) && (1==minimalAlignment)) {
    // This is a recursive search for secondary linker.
    // We found secondary linker.
    // Repair the read? Too complicated!
    // Just delete the read.
    fprintf(logFile, 
            "Secondary linker found. Calling proc will delete the read.\n");
    // Cannot print read ID. This variable is undefined at this point. 
    // AS_UID_toString(fr->gkFragment_getReadUID())); 
    return (LINKER_POSITIVE);
  }  

  if (minimalAlignment == 0) {
    //  No match after trying all possible linkers. 
    //  Signal no change to read.
    return (LINKER_NEGATIVE);
  }

  if (fractionalAlignment == 1) { // or functionalAlignment
    if ((lSize < MATE_MIN_LEN) && (rSize < MATE_MIN_LEN)) {
      //  Linker found but read too short. 
      //  Not enough sequence on either side of the linker to make a mate. 
      fprintf(logFile, 
              "Read %s (len=%d) is nearly all linker. Linker at %d-%d. Read deleted.\n",
              AS_UID_toString(fr->gkFragment_getReadUID()), 
              fr->gkFragment_getSequenceLength(),
              al.begJ, al.endJ);
      // Delete the read.
      fr->gkFragment_setReadUID(AS_UID_undefined());
      fr->gkFragment_setIsDeleted(1);
      return (LINKER_POSITIVE);
    }

    if ((lSize < MATE_MIN_LEN)) {
      //  Linker on the left.
      //  Sufficient sequence left on the right.
      fprintf(logFile, "Trim linker from left side of %s\n",
              AS_UID_toString(fr->gkFragment_getReadUID()));
      // Trim the linker.
      uint32 oldLen = fr->gkFragment_getSequenceLength();
      fr->gkFragment_setLength(rSize + (oldLen - fr->clrEnd));
      char *seq = fr->gkFragment_getSequence();
      char *qlt = fr->gkFragment_getQuality();
      memmove(seq, seq + al.endJ, rSize + (oldLen - fr->clrEnd));
      memmove(qlt, qlt + al.endJ, rSize + (oldLen - fr->clrEnd));
      seq[rSize + (oldLen - fr->clrEnd)] = 0;
      qlt[rSize + (oldLen - fr->clrEnd)] = 0;
      fr->clrBgn = 0;
      fr->clrEnd = (fr->clrEnd < al.endJ) ? 0 : (fr->clrEnd - al.endJ);
      if (fr->maxEnd >= fr->maxBgn) {
         fr->maxBgn = (fr->maxBgn < al.endJ) ? 0 : (fr->maxBgn - al.endJ);
         fr->maxEnd = (fr->maxEnd < al.endJ) ? 0 : (fr->maxEnd - al.endJ);
      }
      fr->vecBgn = 1;
      fr->vecEnd = 0;
      assert(fr->clrEnd <= rSize);
      // Launch recursive search for linker in the newly trimmed read.
      if (0 != processMate(fr, NULL, NULL, linker, search, 0)) {  // low stringency, do not split
        fprintf(logFile, "Found more linker after left trim. Delete %s.\n",
                AS_UID_toString(fr->gkFragment_getReadUID()));
      	// Delete the read.
      	fr->gkFragment_setReadUID(AS_UID_undefined());
      	fr->gkFragment_setIsDeleted(1);
      }
      return (LINKER_POSITIVE);
    }

    if ((rSize < MATE_MIN_LEN)) {
      //  Linker on the right.
      //  Sufficiently long sequence on the left.
      fprintf(logFile, "Trim linker from right side of %s\n",
              AS_UID_toString(fr->gkFragment_getReadUID()));
      // Trim the linker.
      fr->gkFragment_setLength(al.begJ);
      char *seq = fr->gkFragment_getSequence();
      char *qlt = fr->gkFragment_getQuality();
      seq[al.begJ] = 0;
      qlt[al.begJ] = 0;
      fr->clrEnd = MIN(fr->clrEnd, al.begJ);
      fr->maxEnd = MIN(fr->maxEnd, al.begJ);
      fr->vecEnd = MIN(fr->vecEnd, al.begJ);
      // Launch recursive search for linker in the newly trimmed read.
      if (0 != processMate(fr, NULL, NULL, linker, search, 0)) {  // low stringency, do not split
        fprintf(logFile, "Found more linker after right trim. Delete %s.\n",
                AS_UID_toString(fr->gkFragment_getReadUID()));
      	// Delete the read.
      	fr->gkFragment_setReadUID(AS_UID_undefined());
      	fr->gkFragment_setIsDeleted(1);
      }
      return (LINKER_POSITIVE);
    }
  }

  if (1==fractionalAlignment && 0==functionalAlignment) {
    //  We have a fractional alignment.
    //  We do not have a functional alignment, so it must be partial.
    //  The alignment is not on the left or right, so it must be in the middle.
    assert (!(rSize < MATE_MIN_LEN));
    assert (!(lSize < MATE_MIN_LEN));
    //
    //  What shall be done?
    //  We choose not to split this read into mates.
    //  We choose not to carve out the larger of the left & right non-linker.
    //  We choose to delete the read. This is just for simplicity.
    fprintf(logFile, 
            "Partial linker in the middle (position %d-%d).  Delete read %s (len %d).\n",
            al.begJ, al.endJ, 
            AS_UID_toString(fr->gkFragment_getReadUID()), 
            fr->gkFragment_getSequenceLength());
    // Delete the read.
    fr->gkFragment_setReadUID(AS_UID_undefined());
    fr->gkFragment_setIsDeleted(1);
    return (LINKER_POSITIVE);
  }
    
  if (1==functionalAlignment) {
    //  Linker found in the middle with enough sequence to make two mated reads.
    assert (1==minimalAlignment);
    assert (1==fractionalAlignment);
    assert ((lSize >= MATE_MIN_LEN) && (rSize >= MATE_MIN_LEN));
      
    // SPLIT THE READ INTO TWO MATES!
    fprintf(logFile, "Split read %s into Mates.\n",
            AS_UID_toString(fr->gkFragment_getReadUID()));
    
    //  0.  Copy the fragments to new mated fragments
    //      CANNOT just copy fr over m1 -- that nukes seq/qlt pointers!
    //memcpy(m1, fr, sizeof(gkFragment));
    //memcpy(m2, fr, sizeof(gkFragment));

    m1->gkFragment_setType(GKFRAGMENT_NORMAL);
    m2->gkFragment_setType(GKFRAGMENT_NORMAL);

    m1->gkFragment_setLibraryIID(1);
    m2->gkFragment_setLibraryIID(1);
    
    //  1.  Make new UIDs for the two mated reads.  Nuke the old
    //  read.  Make the mates.
    //
    //  WARNING!  See those getLastElemStore() below?  It forces us to
    //  load the gkm1 read before the gkm2 read.
    {
      char  uid[64];
      int   len;
      strcpy(uid, AS_UID_toString(fr->gkFragment_getReadUID()));
      len = strlen(uid);
      uid[len+1] = 0;      
      uid[len] = 'a';
      m1->gkFragment_setReadUID(AS_UID_load(uid));
      m1->gkFragment_setIsDeleted(0);
      uid[len] = 'b';
      m2->gkFragment_setReadUID(AS_UID_load(uid));
      m2->gkFragment_setIsDeleted(0); 
      fr->gkFragment_setReadUID(AS_UID_undefined());
      fr->gkFragment_setIsDeleted(1);
      m1->gkFragment_setMateIID(gkpStore->gkStore_getNumFragments() + 2);
      m2->gkFragment_setMateIID(gkpStore->gkStore_getNumFragments() + 1);	
      m1->gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
      m2->gkFragment_setOrientation(AS_READ_ORIENT_INNIE);
    }
      
    //  2.  Propagate clear ranges.  Math.
    //
    //  m1 is reverse complemented, so the start of m1 is next to the
    //  linker, and the end can extend into low quality sequence at
    //  the start of the read.
    //
    //  lSize - size of the left half of the read, excluding X's
    //  rSize - size of the right half of the read, excluding X's
    //
    //       v clearBeg                   clearEnd v
    //  XXXXXX-------------------[linker]----------XXXXXXXXXXX
    //                   al.begJ ^      ^ al.endJ
    //
    //  
    m1->clrBgn = 0;
    m1->clrEnd = lSize;
    if (fr->maxEnd < fr->maxBgn) {
      m1->maxBgn = fr->maxBgn;
      m1->maxEnd = fr->maxEnd;
    } else {
       m1->maxBgn = 0;
       m1->maxEnd = lSize;
    }
    m1->vecBgn = m1->tntBgn = 1;
    m1->vecEnd = m1->tntEnd = 0;       
    m2->clrBgn = 0;
    m2->clrEnd = (fr->clrEnd < al.endJ) ? 0 : (fr->clrEnd - al.endJ);
    if (fr->maxEnd < fr->maxBgn) {
      m2->maxBgn = fr->maxBgn;
      m2->maxEnd = fr->maxEnd;      
    } else {    
      m2->maxBgn = (fr->maxBgn < al.endJ) ? 0 : (fr->maxBgn - al.endJ);
      m2->maxEnd = (fr->maxEnd < al.endJ) ? 0 : (fr->maxEnd - al.endJ);
    }
    m2->vecBgn = m2->tntBgn = 1;
    m2->vecEnd = m2->tntEnd = 0;
      
    //  3.  Construct new rm1, rm2.  Nuke the linker.  Reverse
    //  complement -- inplace -- the left mate.
    {
      char  *seq = m1->gkFragment_getSequence();
      char  *qlt = m1->gkFragment_getQuality();
      assert(fr->clrBgn >= 0);
      assert(lSize > 0);
      memmove(seq, fr->gkFragment_getSequence(), al.begJ);
      memmove(qlt, fr->gkFragment_getQuality() , al.begJ);
      reverseComplement(seq, qlt, al.begJ);
      m1->gkFragment_setLength(al.begJ);
      seq[al.begJ] = 0;
      qlt[al.begJ] = 0;
      assert(strlen(seq) == al.begJ);
    }
      
    {
      char  *seq = m2->gkFragment_getSequence();
      char  *qlt = m2->gkFragment_getQuality();
      assert(al.endJ >= 0);
      assert(rSize > 0);
      memmove(seq, fr->gkFragment_getSequence() + al.endJ, rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd));
      memmove(qlt, fr->gkFragment_getQuality()  + al.endJ, rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd));
      m2->gkFragment_setLength(rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd));
      seq[rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd)] = 0;
      qlt[rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd)] = 0;
      assert(strlen(seq) == rSize + (fr->gkFragment_getSequenceLength() - fr->clrEnd));
    }
    
    fprintf(logFile, "Mates '%s' (%d-%d) and '%s' (%d-%d) created.\n",
            AS_UID_toString(m1->gkFragment_getReadUID()), 0, al.begJ,
            AS_UID_toString(m2->gkFragment_getReadUID()), al.endJ, al.lenB);

    //  4.  Recurseive search for linker in the left and right mates.
    //      Issue low-stringency search and don't allow any more splitting.

    if (processMate(m1, NULL, NULL, linker, search, 0) ||
	processMate(m2, NULL, NULL, linker, search, 0)) {
      fprintf(logFile, "Found more linker after split. Delete mates %s,%s.\n",
              AS_UID_toString(m1->gkFragment_getReadUID()),
              AS_UID_toString(m2->gkFragment_getReadUID()));
      // The read is already deleted.
      // Must delete the mates.
      m1->gkFragment_setReadUID(AS_UID_undefined());
      m1->gkFragment_setIsDeleted(1);
      m2->gkFragment_setReadUID(AS_UID_undefined());
      m2->gkFragment_setIsDeleted(1);
    } else {
      fprintf(logFile, "Mates (%s,%s) survived recursive search for linker.\n",
              AS_UID_toString(m1->gkFragment_getReadUID()),
              AS_UID_toString(m2->gkFragment_getReadUID()));
    }
  
  return (LINKER_POSITIVE);
  } // if functionalAlignment
  
  //  There was significant alignment.
  //  To reach here means we didn't do anything with it.  
  //  Is this a problem? Not sure.
  
  assert (minimalAlignment == 1);
    
  fprintf(logFile, 
          "Linker detected in '%s' (%d bp, %d match).  Mark %d-%d as possible contaminant.\n",
          AS_UID_toString(fr->gkFragment_getReadUID()), 
          al.alignLen, al.matches,
          al.begJ, al.endJ);
  
  //  Action: mark the aligned region as possible contaminant.
  //  Put the marks in the gatekeeper store (read database).
  //  Hope that the OBT module can check the contaminant region for excessive overlaps.
  
  //  Overload some of the later clear ranges to convey coordinates to downstream processes.
  //  We have 64 bits split into 4 16 bit words.  We want to store the
  //  position of the match in both the linker and the read, as well
  //  as length and number of matches.  Six things.  If we assume the
  //  linker is 255bp or smaller, we can pack.
  //  (This could benefit from having 64-bits of generic unioned data in the fragment record.)
  
  fr->tntBgn = al.begJ;
  fr->tntEnd = al.endJ;
  
  return (LINKER_NEGATIVE);
}


int
detectMates(char *linker[AS_LINKER_MAX_SEQS], int search[AS_LINKER_MAX_SEQS]) {
  gkFragment    fr;
  gkFragment    m1;
  gkFragment    m2;
  int readChanged = 0;

  fr.gkFragment_enableGatekeeperMode(gkpStore);
  m1.gkFragment_enableGatekeeperMode(gkpStore);
  m2.gkFragment_enableGatekeeperMode(gkpStore);

  globalMatrix = (dpMatrix *)safe_malloc(sizeof(dpMatrix));

  int32  lastElem = gkpStore->gkStore_getNumFragments();

  fprintf(stderr, "detectMates()-- from %d to %d\n", 1, lastElem);

  for (int32 thisElem=1; thisElem<=lastElem; thisElem++) {
    if ((thisElem % 1000000) == 0)
      fprintf(stderr, "detectMates()--  at %d\n", thisElem);

    gkpStore->gkStore_getFragment(thisElem, &fr, GKFRAGMENT_QLT);

    if (fr.gkFragment_getIsDeleted()) {
      st.notExaminedForLinker++;  //  because it was deleted already
      continue;
    }

    if (fr.clrBgn >= fr.clrEnd || fr.clrEnd <= 0) {
      // assert(fr.clrBgn < fr.clrEnd);  // REMOVED BY JASON
      // This happens when read begins NN... and user option trim = pair-of-n.
      st.notExaminedForLinker++;  //  because it was deleted already
      continue;
    }
      
    m1.gkFragment_setType(GKFRAGMENT_NORMAL);
    m1.gkFragment_setReadUID(AS_UID_undefined());
    m1.gkFragment_setIsDeleted(1);

    m2.gkFragment_setType(GKFRAGMENT_NORMAL);
    m2.gkFragment_setReadUID(AS_UID_undefined());
    m2.gkFragment_setIsDeleted(1);

    //  If processMate returns true, something changed.  Delete the
    //  original read, and add either the new single read, or the two
    //  mates.
    //
    //  WARNING!  The mates MUST be added in this order, otherwise,
    //  the UID<->IID mapping will be invalid.

    readChanged = processMate(&fr, &m1, &m2, linker, search, 1); // high stringency; split if found
    if (0==readChanged) 
      readChanged = processMate(&fr, &m1, &m2, linker, search, 0); // low stringency; trim if found

    //  Now figure out what happened.
      
    if (readChanged == 0) {
      st.noLinker++;

    } else if (fr.gkFragment_getIsDeleted() == 0) {
      st.partialLinker++;

      assert(AS_UID_isDefined(m1.gkFragment_getReadUID()) == 0);
      assert(AS_UID_isDefined(m2.gkFragment_getReadUID()) == 0);

      gkpStore->gkStore_delFragment(thisElem);
      gkpStore->gkStore_addFragment(&fr);

    } else if ((m1.gkFragment_getIsDeleted() == 0) &&
               (m2.gkFragment_getIsDeleted() == 0)) {
      st.fullLinker++;

      assert(AS_UID_isDefined(fr.gkFragment_getReadUID()) == 0);

      gkpStore->gkStore_delFragment(thisElem);
      gkpStore->gkStore_addFragment(&m1);
      gkpStore->gkStore_addFragment(&m2);

    } else if ((fr.gkFragment_getIsDeleted() == 1) &&
               (m1.gkFragment_getIsDeleted() == 1) &&
               (m2.gkFragment_getIsDeleted() == 1)) {
      st.badLinker++;

      assert(AS_UID_isDefined(fr.gkFragment_getReadUID()) == 0);
      assert(AS_UID_isDefined(m1.gkFragment_getReadUID()) == 0);
      assert(AS_UID_isDefined(m2.gkFragment_getReadUID()) == 0);

      gkpStore->gkStore_delFragment(thisElem);

    } else {
      fprintf(stderr, "ERROR:  linker found, but we failed to handle it.\n");
      exit(1);
    }
  }

  safe_free(globalMatrix);
  globalMatrix = NULL;

  return(0);
}






void
addLibrary(char *libraryName,
           int   insertSize,
           int   insertStdDev,
           int   haveLinker)  {
  gkLibrary   gkl;

  gkl.libraryUID = AS_UID_load(libraryName);

  gkl.forceBOGunitigger          = 1;

  gkl.doNotTrustHomopolymerRuns  = 1;

  gkl.doRemoveDuplicateReads     = 1;
  gkl.doNotQVTrim                = 1;
  gkl.goodBadQVThreshold         = 1;
  gkl.doNotOverlapTrim           = 0;

  if (haveLinker == FALSE) {
    gkl.mean        = 0;
    gkl.stddev      = 0;
    gkl.orientation = AS_READ_ORIENT_UNKNOWN;
  } else {
    gkl.mean        = insertSize;
    gkl.stddev      = insertStdDev;
    gkl.orientation = AS_READ_ORIENT_INNIE;
  }

  gkpStore->gkStore_addLibrary(gkl.libraryUID, &gkl);

  assert(gkpStore->gkStore_getNumLibraries() == 1);
}


//  This is an efficient version of dumpGateKeeperAsFRG() in AS_GKP_dump.c
void
dumpFragFile(char *outName, FILE *outFile) {
  gkFragment        fr;

  GenericMesg       pmesg;
  LibraryMesg       libMesg;
  FragMesg          frgMesg;
  LinkMesg          lnkMesg;

  int               i;

  AS_UID           *frgUID = new AS_UID [gkpStore->gkStore_getNumFragments() + 1];

  //  Dump the format message
  //
  {
    VersionMesg  vmesg;

    AS_MSG_setFormatVersion(2);

    vmesg.version = 2;

    pmesg.m = &vmesg;
    pmesg.t = MESG_VER;

    WriteProtoMesg_AS(outFile, &pmesg);
  }

  //  Exactly one library here.

  {
    gkLibrary  gkl;

    gkpStore->gkStore_getLibrary(1, &gkl);

    frgUID[0] = gkl.libraryUID;

    pmesg.m = &libMesg;
    pmesg.t = MESG_LIB;

    libMesg.action       = AS_ADD;
    libMesg.eaccession   = gkl.libraryUID;
    libMesg.mean         = gkl.mean;
    libMesg.stddev       = gkl.stddev;
    libMesg.source       = NULL;

    libMesg.link_orient.setIsUnknown();

    switch(gkl.orientation) {
      case AS_READ_ORIENT_INNIE:
        libMesg.link_orient.setIsInnie();
        break;
      case AS_READ_ORIENT_OUTTIE:
        libMesg.link_orient.setIsOuttie();
        break;
      case AS_READ_ORIENT_NORMAL:
        libMesg.link_orient.setIsNormal();
        break;
      case AS_READ_ORIENT_ANTINORMAL:
        libMesg.link_orient.setIsAnti();
        break;
      case AS_READ_ORIENT_UNKNOWN:
        libMesg.link_orient.setIsUnknown();
        break;
      default:
        //  Cannot happen, unless someone adds a new orientation to gkFragment.
        assert(0);
        break;
    }

    gkl.gkLibrary_encodeFeatures(&libMesg);

    WriteProtoMesg_AS(outFile, &pmesg);

    gkl.gkLibrary_encodeFeaturesCleanup(&libMesg);
  }

  //  Dump fragments -- as soon as both reads in a mate are defined,
  //  we dump the mate relationship.

  gkStream *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_QLT);

  while (fs->next(&fr)) {
    if (fr.gkFragment_getIsDeleted())
      continue;

    frgUID[fr.gkFragment_getReadIID()] = fr.gkFragment_getReadUID();

    st.fragmentsOutput++;

    pmesg.m = &frgMesg;
    pmesg.t = MESG_FRG;

    //  This code used in AS_GKP_dump.c (dumpFRG).
    frgMesg.action            = fr.gkFragment_getIsDeleted() ? AS_DELETE : AS_ADD;
    frgMesg.eaccession        = fr.gkFragment_getReadUID();
    frgMesg.library_uid       = frgUID[0];
    frgMesg.library_iid       = fr.gkFragment_getLibraryIID();
    frgMesg.plate_uid         = AS_UID_undefined();
    frgMesg.plate_location    = 0;
    frgMesg.type              = AS_READ;
    frgMesg.is_random         = (fr.gkFragment_getIsNonRandom()) ? 0 : 1;
    frgMesg.status_code       = 'G';
    frgMesg.clear_rng.bgn     = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_CLR);
    frgMesg.clear_rng.end     = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_CLR);
    frgMesg.clear_vec.bgn     = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC);
    frgMesg.clear_vec.end     = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC);
    frgMesg.clear_max.bgn     = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_MAX);
    frgMesg.clear_max.end     = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_MAX);
    frgMesg.contamination.bgn = fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_TNT);
    frgMesg.contamination.end = fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_TNT);
    frgMesg.source            = NULL;
    frgMesg.sequence          = fr.gkFragment_getSequence();
    frgMesg.quality           = fr.gkFragment_getQuality();
    frgMesg.hps               = NULL;
    frgMesg.iaccession        = fr.gkFragment_getReadIID();

    WriteProtoMesg_AS(outFile, &pmesg);

    if ((fr.gkFragment_getMateIID() > 0) &&
        (fr.gkFragment_getMateIID() < fr.gkFragment_getReadIID())) {
      st.matesOutput++;
      st.fragmentsOutput--;
      st.fragmentsOutput--;

      pmesg.m = &lnkMesg;
      pmesg.t = MESG_LKG;

      //  The link_orient is not used here.  These should be dumped as
      //  version 2 fragments.

      lnkMesg.action      = AS_ADD;
      lnkMesg.type.setIsMatePair();
      lnkMesg.link_orient.setIsUnknown();
      lnkMesg.frag1       = frgUID[fr.gkFragment_getMateIID()];
      lnkMesg.frag2       = fr.gkFragment_getReadUID();
      lnkMesg.distance    = frgUID[0];

      WriteProtoMesg_AS(outFile, &pmesg);
    }
  }

  delete [] frgUID;

  delete fs;
}


int
main(int argc, char **argv) {
  int       insertSize       = 0;
  int       insertStdDev     = 0;
  char     *libraryName      = 0L;
  int       firstFileArg     = 0;

  char      oPrefix[FILENAME_MAX] = {0};
  char      frgName[FILENAME_MAX] = {0};
  char      gkpName[FILENAME_MAX] = {0};
  char      logName[FILENAME_MAX] = {0};
  char      stsName[FILENAME_MAX] = {0};

  bool      doDeDup          = 1;

  // initialize linker search structure
  // One array stores the character sequences of the linker
  // A boolean array stores which linkers are to be used in the search

  int       haveLinker        = FALSE;
  int       invalidLinkerSeq  = FALSE;
  char     *linker[AS_LINKER_MAX_SEQS] = { NULL };
  int       search[AS_LINKER_MAX_SEQS] = { 0    };
  
  // the first slot of the linker array is the FLX mate pair linker (which is a palindrome)
  char   *linkerFLX     = linker[0] = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";  // palindrome
  // the next two slots are the Titanium linker. It requires two linkers because they are not palindromes
  char   *linkerFIX     = linker[1] = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";    // linker for Titanium reads
  char   *linkerXIF     = linker[2] = "CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA";    // rc of linker for Titanium reads
  // subsequent linkers will be used for future barcoding
  // final linkers are custom, provided by the user, filled in when parsing parameters

  int       bogusOptions[256] = {0};
  int       bogusOptionsLen   = 0;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-insertsize") == 0) {
      insertSize   = atoi(argv[++arg]);
      insertStdDev = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-libraryname") == 0) {
      libraryName = argv[++arg];

    } else if (strcmp(argv[arg], "-clear") == 0) {
      arg++;

      //  If this is the first time we get a -clear switch, set
      //  clearAction to exactly that value.  Later times through,
      //  we'll add in more options.
      //
      if      (strcasecmp(argv[arg], "all") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_ALL) : (clearAction | CLEAR_ALL);
      else if (strcasecmp(argv[arg], "454") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_454) : (clearAction | CLEAR_454);
      else if (strcasecmp(argv[arg], "n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_N) : (clearAction | CLEAR_N);
      else if (strcasecmp(argv[arg], "pair-of-n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_PAIR_N) : (clearAction | CLEAR_PAIR_N);
      else if (strcasecmp(argv[arg], "discard-n") == 0)
        clearAction = (clearSet == 0) ? (CLEAR_DISCARD_N) : (clearAction | CLEAR_DISCARD_N);
      else {
        clearAction = (clearSet == 0) ? (CLEAR_ERRR) : (clearAction | CLEAR_ERRR);
        err++;
      }
      clearSet++;

    } else if (strcmp(argv[arg], "-trim") == 0) {
      arg++;

      if      (strcasecmp(argv[arg], "none") == 0)
        trimAction = TRIM_NONE;
      else if (strcasecmp(argv[arg], "soft") == 0)
        trimAction = TRIM_SOFT;
      else if (strcasecmp(argv[arg], "hard") == 0)
        trimAction = TRIM_HARD;
      else if (strcasecmp(argv[arg], "chop") == 0)
        trimAction = TRIM_CHOP;
      else {
        trimAction = TRIM_ERRR;
        err++;
      }

    } else if (strcmp(argv[arg], "-linker") == 0) {
      arg++;

      if      (strcasecmp(argv[arg], "flx") == 0) {
        search[0]    = TRUE;
        haveLinker   = TRUE;
      }
      else if (strcasecmp(argv[arg], "titanium") == 0) {
        search[1]    = TRUE;
        search[2]    = TRUE;
        haveLinker   = TRUE;
      }
      else {
        int start = AS_LINKER_CUSTOM_OFFSET;
        if (AS_UTL_isValidSequence(argv[arg], strlen(argv[arg]))) {
          linker[start]      = argv[arg];
          search[start++]    = TRUE;
          haveLinker         = TRUE;
        } else {
          invalidLinkerSeq     = TRUE;
          err++;
        }
      }

    } else if (strcmp(argv[arg], "-nodedup") == 0) {
      doDeDup = 0;

    } else if (strcmp(argv[arg], "-output") == 0) {
      strcpy(oPrefix, argv[++arg]);

    } else {
      if (argv[arg][0] == '-') {
        bogusOptions[bogusOptionsLen++] = arg;
        err++;
      } else {
        firstFileArg = arg;
        arg          = argc;
      }
    }

    arg++;
  }

  //  Have a linker but no insert size?  Error.
  if ((haveLinker) && ((insertSize == 0) || (insertStdDev == 0)))
    err++;

  //  Have an insert size but no linker?  Error.
  if ((!haveLinker) && ((insertSize != 0) || (insertStdDev != 0)))
    err++;

  if ((err) || (libraryName == 0L) || (oPrefix[0] == 0) || (firstFileArg == 0)) {
    fprintf(stderr, "usage: %s [opts] -libraryname LIB -output NAME IN.SFF ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -insertsize i d        Mates are on average i +- d bp apart.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -libraryname n         The UID of the library these reads are added to.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -clear all             Use the whole read.\n");
    fprintf(stderr, "  -clear 454             Use the 454 clear ranges as is (default).\n");
    fprintf(stderr, "  -clear n               Use the whole read up to the first N.\n");
    fprintf(stderr, "  -clear pair-of-n       Use the whole read up to the frist pair of Ns.\n");
    fprintf(stderr, "  -clear discard-n       Delete the read if there is an N in the clear range.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  If multiple -clear options are supplied, the intersection is used.  For\n");
    fprintf(stderr, "  'discard-n', the clear range is first computed, then if there is still an\n");
    fprintf(stderr, "  N in the clear range, the read is deleted.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Caution!  Even though the default is '454', when any -clear option is used,\n");
    fprintf(stderr, "  the list of clear ranges to intersect is reset.  To get both '454' and 'n',\n");
    fprintf(stderr, "  BOTH '-clear 454' and '-clear n' must be supplied on the command line.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -trim none             Use the whole read regardless of -clear settings.\n");
    fprintf(stderr, "  -trim soft             OBT and ECR can increase the clear range.\n");
    fprintf(stderr, "  -trim hard             OBT can only shrink the clear range, but ECR can extend (default).\n");
    fprintf(stderr, "  -trim chop             Erase sequence outside the clear range.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  'none' will emit the whole read, and reset clear ranges to cover the whole read.\n");
    fprintf(stderr, "  'soft' will emit the whole read, and leave clear ranges as set.\n");
    fprintf(stderr, "  'hard' is like soft, with the addition of a 'clm' message to stop OBT.\n");
    fprintf(stderr, "  'chop' is like none, but after the read is chopped down to just the clear bases.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -linker [name | seq]   Search for linker, create mated reads.\n");
    fprintf(stderr, "                         Name is one of:\n");
    fprintf(stderr, "                           'flx'      == %s\n",     linkerFLX);
    fprintf(stderr, "                           'titanium' == %s and\n", linkerFIX);
    fprintf(stderr, "                                         %s\n",     linkerXIF);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -nodedup               Do not remove reads that are a perfect prefix of another read.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -output name           Write output to files prefixed with 'name'.  Three files are created:\n");
    fprintf(stderr, "                           name.frg   -- CA format fragments.\n");
    fprintf(stderr, "                           name.log   -- Actions taken; deleted fragments, mate splits, etc.\n");
    fprintf(stderr, "                           name.stats -- Human-readable statistics.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "See http://apps.sourceforge.net/mediawiki/wgs-assembler/index.php?title=Formatting_Inputs\n");
    fprintf(stderr, "\n");

    for (err=0; err<bogusOptionsLen; err++)
      fprintf(stderr, "ERROR:  Unknown option '%s'\n", argv[bogusOptions[err]]);

    if (libraryName == 0L)
      fprintf(stderr, "ERROR:  Need to supply -libraryname.\n");

    if (oPrefix[0] == 0)
      fprintf(stderr, "ERROR:  Need to supply -output.\n");

    if (firstFileArg == 0)
      fprintf(stderr, "ERROR:  Need to supply some SFF files.\n");

    if ((haveLinker) && ((insertSize == 0) ||
                         (insertStdDev == 0)))
      fprintf(stderr, "ERROR:  Have a linker sequence, but no insert size set with -insertsize.\n");

    if ((!haveLinker) && ((insertSize != 0) ||
                          (insertStdDev != 0)))
      fprintf(stderr, "ERROR:  Have an insert size, bu no linker sequence set with -linker.\n");

    if (clearAction == CLEAR_ERRR)
      fprintf(stderr, "ERROR:  Unknown -clear value.\n");

    if (trimAction == TRIM_ERRR)
      fprintf(stderr, "ERROR:  Unknown -trim value.\n");

    if (invalidLinkerSeq == TRUE)
      fprintf(stderr, "ERROR:  Invalid -linker value. It must be one of titanium, flx, or a valid ACGT string.\n");
    
    exit(1);
  }

  {
    int32  oLen = strlen(oPrefix);

    if ((oPrefix[oLen-4] == '.') &&
        (oPrefix[oLen-3] == 'f') &&
        (oPrefix[oLen-2] == 'r') &&
        (oPrefix[oLen-1] == 'g'))
      oPrefix[oLen-4] = 0;
  }

  strcpy(frgName, oPrefix);  strcat(frgName, ".frg");
  strcpy(gkpName, oPrefix);  strcat(gkpName, ".tmpStore");
  strcpy(logName, oPrefix);  strcat(logName, ".log");
  strcpy(stsName, oPrefix);  strcat(stsName, ".stats");

  if (AS_UTL_fileExists(frgName, FALSE, FALSE))
    fprintf(stderr, "ERROR: Output file '%s' exists; I will not clobber it.\n", frgName), exit(1);

  errno = 0;
  logFile = fopen(logName, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open the log file '%s': %s\n", logName, strerror(errno)), exit(1);

  errno = 0;
  FILE *frgFile = fopen(frgName, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open the output file '%s': %s\n", frgName, strerror(errno)), exit(1);

  if (AS_UTL_fileExists(gkpName, TRUE, FALSE)) {
    fprintf(stderr, "ERROR: Temporary Gatekeeper Store '%s' exists; I will not clobber it.\n", gkpName);
    fprintf(stderr, "       If this is NOT from another currently running sffToCA, simply remove this directory.\n");
    exit(1);
  }

  gkpStore = new gkStore(gkpName, TRUE, TRUE);

  addLibrary(libraryName, insertSize, insertStdDev, haveLinker);

  for (int32 file=firstFileArg; file < argc; file++)
    loadSFF(argv[file]);

  if (doDeDup)
    removeDuplicateReads();

  if (haveLinker)
    detectMates(linker, search);

  dumpFragFile(frgName, frgFile);

  gkpStore->gkStore_delete();
  delete gkpStore;

  errno = 0;
  fclose(frgFile);
  if (errno)
    fprintf(stderr, "Failed to close '%s': %s\n", frgName, strerror(errno)), exit(1);

  errno = 0;
  if ((logFile) && (logFile != stderr))
    fclose(logFile);
  if (errno)
    fprintf(stderr, "Failed to close '%s': %s\n", logName, strerror(errno)), exit(1);

  writeStatistics(argv, argc, firstFileArg, frgName, haveLinker, linker, search, stsName);

  return(0);
}
