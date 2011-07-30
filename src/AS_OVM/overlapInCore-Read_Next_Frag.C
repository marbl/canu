
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

static const char *rcsid = "$Id: overlapInCore-Read_Next_Frag.C,v 1.1 2011-07-30 01:16:05 brianwalenz Exp $";

#include "overlapInCore.H"



/* Read the next fragment from fragment stream  stream  and store it in  frag ,
 *  with its quality values in  quality .  Put the read itself in  myRead
 *  and the screen information in  screen .  Set  (* last_frag_read)
 *  to the iid of the fragment read.
 *  If successful, return  TRUE ; if there is no fragment, return  FALSE . */

int
Read_Next_Frag(char frag [AS_READ_MAX_NORMAL_LEN + 1],
               char quality [AS_READ_MAX_NORMAL_LEN + 1],
               gkStream *stream,
               gkFragment *myRead,
               uint32 * last_frag_read,
               uint32 minLibToRead, uint32 maxLibToRead) {
  int  return_val, success, match_ct;
  size_t  i, frag_len;
  char   *seqptr;
  char   *qltptr;

  //  BPW says we don't need to mutex this
  success = stream->next (myRead);

  if  (! success)
    return(0);

  *last_frag_read = myRead->gkFragment_getReadIID ();

  if  (myRead->gkFragment_getIsDeleted ())
    return (DELETED_FRAG);

  if (minLibToRead != 0 && (myRead->gkFragment_getLibraryIID () < minLibToRead || myRead->gkFragment_getLibraryIID() > maxLibToRead)) {
    return (DELETED_FRAG);
  }

  //  We got a read!  Lowercase it, adjust the quality, and extract
  //  the clear region.

  seqptr   = myRead->gkFragment_getSequence();
  qltptr   = myRead->gkFragment_getQuality();
  frag_len = myRead->gkFragment_getSequenceLength();

  for  (i = 0;  i < frag_len;  i ++) {
    frag [i]    = tolower (seqptr [i]);
    quality [i] = qltptr[i] - QUALITY_BASE_CHAR;
  }

  if  (! Ignore_Clear_Range) {
    uint32 clear_start, clear_end;
    myRead->gkFragment_getClearRegion(clear_start, clear_end);

    frag_len = clear_end - clear_start;
    if  (clear_start > 0) {
      memmove (frag, frag + clear_start, frag_len);
      memmove (quality, quality + clear_start, frag_len);
    }
  }

  if  (OFFSET_MASK < frag_len) {
    fprintf (stderr, "ERROR:  Read "F_IID" is too long (%lu) for hash table\n",
             myRead->gkFragment_getReadIID(), frag_len);
    exit (-1);
  }

  frag [frag_len] = '\0';
  quality [frag_len] = '\0';

  return (VALID_FRAG);
}
