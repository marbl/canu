
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
/* $Id: Polygon.cc,v 1.2 2004-09-23 20:25:23 mcschatz Exp $ */
#include <iostream>

#include "Point.h"
#include "Polygon.h"
#include "ConvexHull.h"
#include "PolygonIntersector.h"

bool Polygon<UnitType>::intersects(const Polygon<UnitType> & other) const
{
  PolygonIntersector<UnitType> intersector;
  return intersector.polygonsIntersect(*this, other);
}
