
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

MultiAlignT *
CreateMultiAlignTFromIUM(IntUnitigMesg *ium,
                         int localID,
                         int sequenceOnly) {
  MultiAlignT *ma = CreateMultiAlignT();

  //  if localID >= 0, use that for the IMP 'sourceInt' id, otherwise,
  //  use the 'sourceInt' in the input ium->imp.

  if (ium->length != strlen(ium->consensus))
    fprintf(stderr, "Reported Length of IUM %d (%d) doesnt matches strlen (%d)\n",
            ium->iaccession, ium->length, strlen(ium->consensus));
  assert(ium->length == strlen(ium->consensus));
  assert(ium->length == strlen(ium->quality));

  ma->maID          = ium->iaccession;

  ma->consensus = CreateVA_char(ium->length + 1);
  EnableRangeVA_char(ma->consensus, ium->length + 1);
  strcpy(Getchar(ma->consensus,0), ium->consensus);

  ma->quality = CreateVA_char(ium->length + 1);
  EnableRangeVA_char(ma->quality, ium->length + 1);
  strcpy(Getchar(ma->quality,0), ium->quality);

  if (sequenceOnly == 0) {
    IntUnitigPos unitigPos;
    int          deltai, i;

    deltai = 0;
    for (i=0; i<ium->num_frags; i++)
      deltai += ium->f_list[i].delta_length;
    ma->fdelta = CreateVA_int32(deltai);
    ma->f_list = CreateVA_IntMultiPos(ium->num_frags);

    deltai = 0;
    for (i=0; i<GetMultiAlignLength(ma); i++)
      if (ium->consensus[i] == '-')
        deltai++;
    ma->udelta = CreateVA_int32(deltai);
    ma->u_list = CreateVA_IntUnitigPos(1);

    ma->v_list = CreateVA_IntMultiVar(0);

    for (i=0; i<ium->num_frags; i++) {
      IntMultiPos *cfr =  ium->f_list + i;
      IntMultiPos  tmp = *cfr;

      if (localID < 0)
        tmp.sourceInt = cfr->sourceInt;
      else
        tmp.sourceInt = localID++;

      //  Because Getint32() returns NULL for elements out of range,
      //  we need to set the delta pointer after we append, adjusting
      //  back to the first element we appended. we guarantee that the
      //  pointer is NULL if there are no delta values.
      //
      if (cfr->delta_length > 0) {
        for (deltai=0; deltai<cfr->delta_length; deltai++)
          AppendVA_int32(ma->fdelta, cfr->delta + deltai);
        tmp.delta = Getint32(ma->fdelta, GetNumint32s(ma->fdelta) - cfr->delta_length);
      } else {
        tmp.delta = NULL;
      }

      SetIntMultiPos(ma->f_list, i, &tmp);
    }

    assert(ium->num_frags == GetNumIntMultiPoss(ma->f_list));

    //  Make a unitig

    deltai = 0;
    for (i=0; i<GetMultiAlignLength(ma); i++)
      if (ium->consensus[i] == '-')
        Appendint32(ma->udelta, &deltai);
      else
        deltai++;

    unitigPos.type         = AS_OTHER_UNITIG;
    unitigPos.ident        = ium->iaccession;
    unitigPos.position.bgn = 0;
    unitigPos.position.end = GetMultiAlignLength(ma);
    unitigPos.delta_length = GetNumint32s(ma->udelta);
    unitigPos.delta        = Getint32(ma->udelta,0);

    AppendIntUnitigPos(ma->u_list, &unitigPos);
  }  //  end of (sequenceOnly == 0)

  CheckMAValidity(ma);

  return ma;
}
