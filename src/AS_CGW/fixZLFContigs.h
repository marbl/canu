
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
#ifndef FIXZLFCONTIGS_H
#define FIXZLFCONTIGS_H

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"

typedef struct
{
  CDS_CID_t ident;
  LengthT aEnd;
  LengthT bEnd;
} IEPish;
VA_DEF(IEPish)

typedef struct
{
  CDS_CID_t id;
  cds_int32 zlfUsFound;
  VA_TYPE(MultiAlignT) * zlfUMAs;
} ZLFContig;
VA_DEF(ZLFContig)
  
typedef struct
{
  CDS_CID_t id;
  cds_int32 zlfContigsFound;
  VA_TYPE(ZLFContig) * zlfContigs;
} ZLFScaffold;
VA_DEF(ZLFScaffold)

  
void PopulateContigUnitigCoords(ContigT * contig,
                                VA_TYPE(IEPish) * uCoords);
void PopulateContigScaffoldCoords(CIScaffoldT * scaffold,
                                  HashTable_AS * zlfCUOs,
                                  VA_TYPE(IEPish) * cCoords);
void CreateNewContigsFromUnitigs(ContigT * oldContig,
                                 ZLFContig * zlfContig,
                                 VA_TYPE(IEPish) * uCoords);
int FixZLFContigs(VA_TYPE(ZLFScaffold) * zlfScaffolds,
                  int checkScaffolds,
                  int splitAsNeeded);

void WriteIMPToFile(IntMultiPos * imp, FILE * fp);
void ReadIMPFromFile(IntMultiPos * imp, FILE * fp);

void WriteIMPsToFile(VA_TYPE(IntMultiPos) * imps, FILE * fp);
VA_TYPE(IntMultiPos) * ReadIMPsFromFile(FILE * fp);

void WriteIUPToFile(IntUnitigPos * iup, FILE * fp);
void ReadIUPFromFile(IntUnitigPos * iup, FILE * fp);

void WriteIUPsToFile(VA_TYPE(IntUnitigPos) * iups, FILE * fp);
VA_TYPE(IntUnitigPos) * ReadIUPsFromFile(FILE * fp);

void WriteMAToFile(MultiAlignT * ma, FILE * fp);
int ReadMAFromFile(MultiAlignT * ma, FILE * fp);

#endif
