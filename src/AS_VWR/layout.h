
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
/* $Id: layout.h,v 1.5 2008-06-27 06:29:22 brianwalenz Exp $ */

void SetCurrentFile(char *name);
void ReadAssembly(FILE *);
void BuildAssembly(void);

void LayoutAssembly(void);
void GetAssemblyDims(int *min, int *max, int *rows, int *smallest);

void SetAssemblyMap(int xl, int xh ,int rl, int rh,
                    int w, int h, long b, long t, long m);

void  DrawAssembly(MT_OBJECT *canvas);

#define PICK_NEAREST 0
#define PICK_SEGMENT 1
#define PICK_PIECES  2

char *PickAssembly(MT_OBJECT *canvas, int x, int y ,int v, int mode);
void  PickRelease(MT_OBJECT *canvas);

int  NumberOfStyles(void);
int  StartStyleList(int menu);
int  NextUsedStyle(void);
void SetStyleVis(int id, int menu, int vis);
void GetStyle(int id, long *color, int *thick, int *dash, char **name);

typedef void HistoPacket;
typedef void CoverPacket;
typedef void ObjectPacket;

void QueryTableSize(int);
void QueryTableEntry(char *, CoverPacket *, ObjectPacket *);
void QueryTableFree(void);

char *ProcessQuery(char *query, ObjectPacket **obj, CoverPacket **cover);

void SetQueryColor(ObjectPacket *, long);
void ClearQueryColor(ObjectPacket *);

void  DrawHistogram(MT_OBJECT *canvas, HistoPacket *hist,
                    int xmin, int xmax, int *bsize);
int   DrawHistoInfo(MT_OBJECT *bar, HistoPacket *hist,char *label);
int   DrawHistoInfoLabel(MT_OBJECT *bar, HistoPacket *hist, char *label);
void  GetHistoDims(HistoPacket *hist, int *min, int *max);
void  SetHistoBorderMax(HistoPacket *hist, int minw, int minh);

HistoPacket *Object2Hist(ObjectPacket *obj);
HistoPacket *Cover2Hist(CoverPacket *obj);
CoverPacket *Object2Cover(ObjectPacket *obj);

void  DrawCoverage(MT_OBJECT *canvas, CoverPacket *cover, int xmin, int xmax);

void  FreeCoverage(CoverPacket *cover);
void  FreeHistogram(HistoPacket *hist);
void  FreeObjects(ObjectPacket *obj);
