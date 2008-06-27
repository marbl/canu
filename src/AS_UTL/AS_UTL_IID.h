
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

// $Id: AS_UTL_IID.h,v 1.2 2008-06-27 06:29:21 brianwalenz Exp $

#ifndef AS_UTL_IID_H
#define AS_UTL_IID_H


typedef int32  CDS_CID_t;
#define CDS_CID_MAX     INT32_MAX
#define F_CID           F_S32
#define F_CIDP          F_S32P

typedef int32  CDS_COORD_t;
#define CDS_COORD_MIN   INT32_MIN
#define CDS_COORD_MAX   INT32_MAX
#define F_COORD         F_S32
#define F_COORDP        F_S32P

typedef uint32 AS_IID;
#define CDS_IID_MAX     UINT32_MAX
#define F_IID           F_U32
#define F_IIDP          F_U32P

static
inline
AS_IID
AS_IID_fromString(char *str, char **nxt) {
  return(strtoull(str, nxt, 10));
};

static
inline
int
AS_IID_isDefined(AS_IID iid) {
  return(iid > 0);
};


#endif  //  AS_UTL_IID_H
