
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
/* $Id: Segment.h,v 1.1.1.1 2004-04-14 13:52:08 catmandew Exp $ */
#ifndef SEGMENT_H
#define SEGMENT_H


template <class UnitType>
class Segment
{
public:
  Segment(){}

  Segment(const Point<UnitType> & p1, const Point<UnitType> & p2) :
    pp1(p1), pp2(p2)
    {
    }

  void setPoint1(const Point<UnitType> & p1) {pp1 = p1;}
  void setPoint2(const Point<UnitType> & p2) {pp2 = p2;}
  
  const Point<UnitType> & getPoint1() {return pp1;}
  const Point<UnitType> & getPoint2() {return pp2;}

  
  
private:
  Point<UnitType> pp1;
  Point<UnitType> pp2;
};


#endif
