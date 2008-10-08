
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

#ifndef UID_ERROR_H
#define UID_ERROR_H

static const char *rcsid_UID_ERROR_H = "$Id: SYS_UIDerror.h,v 1.8 2008-10-08 22:03:00 brianwalenz Exp $";

void        SYS_UIDhandleAcceptError(int32 err_code);
void        SYS_UIDhandleRegisterError(int32 err_code);
void        SYS_UIDhandleCreateError(int32 err_code);
void        SYS_UIDhandleActivateError(int32 err_code);
void        SYS_UIDerrorMsg(const char* err_str);

#endif




