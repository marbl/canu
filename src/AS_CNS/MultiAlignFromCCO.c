
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

VA_DEF(SnapMultiPos)
VA_DEF(UnitigPos)

MultiAlignT *
CreateMultiAlignTFromCCO(SnapConConMesg *cco,
                         int localID,
                         int sequenceOnly) {
  MultiAlignT *ma = CreateMultiAlignT();

  //  if localID >= 0, use that for the IMP 'sourceInt' id, otherwise,
  //  use the 'sourceInt' in the input ium->imp.

  if (cco->length != strlen(cco->consensus))
    fprintf(stderr, "Reported Length of CCO %d (%d) doesnt matches strlen (%d)\n",
            cco->iaccession, cco->length, strlen(cco->consensus));
  assert(cco->length == strlen(cco->consensus));
  assert(cco->length == strlen(cco->quality));

  ma->maID          = cco->iaccession;

  ma->consensus = CreateVA_char(cco->length + 1);
  EnableRangeVA_char(ma->consensus, cco->length + 1);
  strcpy(Getchar(ma->consensus,0), cco->consensus);

  ma->quality = CreateVA_char(cco->length + 1);
  EnableRangeVA_char(ma->quality, cco->length + 1);
  strcpy(Getchar(ma->quality,0), cco->quality);

  if (sequenceOnly == 0) {
    IntUnitigPos unitigPos;
    int          deltai, i;

    deltai = 0;
    for (i=0; i<cco->num_pieces; i++)
      deltai += cco->pieces[i].delta_length;
    ma->fdelta = CreateVA_int32(deltai);
    ma->f_list = CreateVA_IntMultiPos(cco->num_pieces);

    ma->udelta = CreateVA_int32(0);
    ma->u_list = CreateVA_IntUnitigPos(cco->num_unitigs);

    ma->v_list = CreateVA_IntMultiVar(cco->num_vars);

    for (i=0; i<cco->num_pieces; i++) {
      SnapMultiPos *cfr =  cco->pieces + i;
      SnapMultiPos  tmp = *cfr;

#ifdef AS_ENABLE_SOURCE
      tmp.source = strdup(cfr->source);
#endif

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

      SetSnapMultiPos(ma->f_list, i, &tmp);
    }

    for (i=0; i<cco->num_vars; i++) {
      IntMultiVar *cvr =  cco->vars + i;
      IntMultiVar  tmp = *cvr;

#warning memory leak converting VAR from CCO to MultiAlignT
      tmp.nr_conf_alleles  = strdup(cvr->nr_conf_alleles);
      tmp.weights          = strdup(cvr->weights);
      tmp.var_seq          = strdup(cvr->var_seq);
      tmp.conf_read_iids   = strdup(cvr->conf_read_iids);

      SetIntMultiVar(ma->v_list, i, &tmp);
    }

    //  Make unitigs -- not authentic, since CCO doesn't retain delta values.

    for (i=0; i<cco->num_unitigs; i++) {
      UnitigPos  unitigPos = {0};

      unitigPos.type         = cco->unitigs[i].type;
      unitigPos.eident       = cco->unitigs[i].eident;
      unitigPos.position     = cco->unitigs[i].position;
      unitigPos.delta_length = 0;
      unitigPos.delta        = NULL;

      AppendUnitigPos(ma->u_list, &unitigPos);
    }
  }  //  end of (sequenceOnly == 0)

  CheckMAValidity(ma);

  return ma;
}
