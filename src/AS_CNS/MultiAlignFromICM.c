
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
CreateMultiAlignTFromICM(IntConConMesg *icm,
                         int localID,
                         int sequenceOnly) {
  MultiAlignT *ma = CreateMultiAlignT();

  //  if localID >= 0, use that for the IMP 'sourceInt' id, otherwise,
  //  use the 'sourceInt' in the input ium->imp.

  if (icm->length != strlen(icm->consensus))
    fprintf(stderr, "Reported Length of ICM %d (%d) doesnt matches strlen (%d)\n",
            icm->iaccession, icm->length, strlen(icm->consensus));
  assert(icm->length == strlen(icm->consensus));
  assert(icm->length == strlen(icm->quality));

  ma->maID          = icm->iaccession;

  ma->consensus = CreateVA_char(icm->length + 1);
  EnableRangeVA_char(ma->consensus, icm->length + 1);
  strcpy(Getchar(ma->consensus,0), icm->consensus);

  ma->quality = CreateVA_char(icm->length + 1);
  EnableRangeVA_char(ma->quality, icm->length + 1);
  strcpy(Getchar(ma->quality,0), icm->quality);

  if (sequenceOnly == 0) {
    IntUnitigPos unitigPos;
    int          deltai, i;

    deltai = 0;
    for (i=0; i<icm->num_pieces; i++)
      deltai += icm->pieces[i].delta_length;
    ma->fdelta = CreateVA_int32(deltai);
    ma->f_list = CreateVA_IntMultiPos(icm->num_pieces);

    ma->udelta = CreateVA_int32(0);
    ma->u_list = CreateVA_IntUnitigPos(icm->num_unitigs);

    ma->v_list = CreateVA_IntMultiVar(icm->num_vars);

    for (i=0; i<icm->num_pieces; i++) {
      IntMultiPos *cfr =  icm->pieces + i;
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

    assert(icm->num_pieces == GetNumIntMultiPoss(ma->f_list));

    for (i=0; i<icm->num_vars; i++) {
      IntMultiVar *cvr =  icm->v_list + i;
      IntMultiVar  tmp = *cvr;

#warning memory leak converting VAR from ICM to MultiAlignT
      tmp.nr_conf_alleles  = strdup(cvr->nr_conf_alleles);
      tmp.weights          = strdup(cvr->weights);
      tmp.var_seq          = strdup(cvr->var_seq);
      tmp.conf_read_iids   = strdup(cvr->conf_read_iids);

      SetIntMultiVar(ma->v_list, i, &tmp);
    }

    //  Make unitigs -- not authentic, since ICM doesn't retain delta values.

    for (i=0; i<icm->num_unitigs; i++) {
      IntUnitigPos  unitigPos = {0};

      unitigPos.type         = icm->unitigs[i].type;
      unitigPos.ident        = icm->unitigs[i].ident;
      unitigPos.position     = icm->unitigs[i].position;
      unitigPos.delta_length = 0;
      unitigPos.delta        = NULL;

      AppendIntUnitigPos(ma->u_list, &unitigPos);
    }
  }  //  end of (sequenceOnly == 0)

  CheckMAValidity(ma);

  return ma;
}
