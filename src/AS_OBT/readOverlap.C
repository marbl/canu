
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

#include "readOverlap.H"

OVSoverlap*
readOverlap(OverlapStore *ovsprimary,
            OverlapStore *ovssecondary) {
  static  OVSoverlap  ovlp;
  static  OVSoverlap  ovls;

  static  uint32      okp = 0;
  static  uint32      oks = 0;

  if (ovsprimary && !okp)
    okp = AS_OVS_readOverlapFromStore(ovsprimary, &ovlp, AS_OVS_TYPE_OBT);

  if (ovssecondary && !oks)
    oks = AS_OVS_readOverlapFromStore(ovssecondary, &ovls, AS_OVS_TYPE_OBT);

  //  Read neither, return null
  if (!okp && !oks)
    return(NULL);

  //  Have both, pick one.
  if (okp && oks) {
    if (ovlp.a_iid <= ovls.a_iid) {
      okp = 0;
      return(&ovlp);
    } else {
      oks = 0;
      return(&ovls);
    }
  }

  //  Only one.

  if (okp) {
    okp = 0;
    return(&ovlp);
  }

  if (oks) {
    oks = 0;
    return(&ovls);
  }
}
