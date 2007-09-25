
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignStore_CNS.h"


MultiAlignStoreT *
CreateMultiAlignStoreT(void) {
  MultiAlignStoreT *mas = (MultiAlignStoreT *)safe_calloc(1, sizeof(MultiAlignStoreT));
  mas->maMax = 1024;
  mas->maPtr = (MultiAlignT **)safe_calloc(mas->maMax, sizeof(MultiAlignT *));
  return(mas);
}

void
ClearMultiAlignStoreT(MultiAlignStoreT *mas) {
  int    i;
#ifdef __FreeBSD__
  //  BPW has had problems with ClearMultiAlignStoreT on FreeBSD --
  //  it's just slow.
  int    u=0;
  for (i=0; i<mas->maMax; i++)
    if (mas->maPtr[i])
      u++;
  fprintf(stderr, "ClearMultiAlignT()-- %d used out of max %d.\n",
          u, mas->maMax);
  fprintf(stderr, "ClearMultiAlignT()-- Bypassing!\n");
  return;
#endif
  for (i=0; i<mas->maMax; i++)
    DeleteMultiAlignT(mas->maPtr[i]);
  memset(mas->maPtr, 0, mas->maMax * sizeof(MultiAlignT *));
}

void 
DeleteMultiAlignStoreT(MultiAlignStoreT *mas) {
  ClearMultiAlignStoreT(mas);
  safe_free(mas->maPtr);
  safe_free(mas);
}

void 
SetMultiAlignInStore(MultiAlignStoreT *mas, int index, MultiAlignT *ma) {
  assert(0 <= index);
  if (mas->maMax <= index) {
    uint32 lastMax = mas->maMax;
    while (mas->maMax <= index)
      mas->maMax *= 2;
    mas->maPtr = (MultiAlignT **)safe_realloc(mas->maPtr, mas->maMax * sizeof(MultiAlignT *));
    memset(mas->maPtr + lastMax, 0, (mas->maMax - lastMax) * sizeof(MultiAlignT *));
  }
  //fprintf(stderr, "SetMultiAlignInStore()-- add ma 0x%016p to index=%d\n", ma, index);

  //  Already here?  Don't do anything.
  if (mas->maPtr[index] == ma)
    return;

  //  Look for duplicate pointers -- this is a nice sanity check --
  //  just kind of slow.
#if 0
  if (ma != NULL) {
    int i;
    for (i=0; i<mas->maMax; i++) {
      if (mas->maPtr[i] == ma) {
        fprintf(stderr, "SetMultiAlignInStore()-- Found ma 0x%016p at location %d, wanted to add it to location %d\n",
                ma, i, index);
        assert(0);
      }
    }
  }
#endif

  DeleteMultiAlignT(mas->maPtr[index]);
  mas->maPtr[index] = ma;
}

MultiAlignT *
GetMultiAlignInStore(MultiAlignStoreT *mas, int index) {
  assert(0 <= index);
  if ((mas->maMax <= index) || (mas->maPtr[index] == NULL))
    //  It's perfectly valid to request something not in the store.
    return(NULL);
  return(mas->maPtr[index]);
}

void
RemoveMultiAlignFromStore(MultiAlignStoreT *mas, int index) {
  assert(0 <= index);
  if (index < mas->maMax) {
    //fprintf(stderr, "RemoveMultiAlignFromStore()--  ma 0x%016lx index=%d\n", mas->maPtr[index], index);
    DeleteMultiAlignT(mas->maPtr[index]);
    mas->maPtr[index] = NULL;
  }
}
