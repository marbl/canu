
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

// $Id: AS_UTL_UID.h,v 1.3 2008-06-16 18:07:43 brianwalenz Exp $

#ifndef AS_UTL_UID_H
#define AS_UTL_UID_H

#include <string.h>

//  The UID interface.
//
//  Historically, UIDs were 64-bit integers.  This was inconvenient
//  for people using "external" data, and for incorporating 454 reads
//  directly.
//
//  To use a UID read from, say, the command line, it must first be
//  converted from a character string to a form recognizable to the
//  assembler:
//
//    AS_UID uid = AS_UID_lookup(character_string);
//
//  This will parse the character_string.  If it is an integer (at
//  most 63-bits), the resulting UID will be "numeric".  Otherwise,
//  the UID is a "string", and the GateKeeper store will be used to
//  convert that string to something the assembler can use.
//
//  If the UID is not known to the store, AS_UID_undefined() is
//  returned, and one should AS_UID_load() the uid into the store.
//
//  To print a UID, use AS_UID_toString(uid) to convert it to a static
//  character string.

typedef struct {
  uint64   isString:1;
  uint64   UID:63;
} AS_UID;


//      isString  UID
//  0)  0         = 0      -- empty UID; invalid
//
//  1)  0         > 0      -- UID is the integer UID
//
//  2)  1        >= 0      -- UID is a pointer to the UID in the store



//  String UIDs absolutely must know who gatekeeper is.  Without it,
//  string UIDs cannot be loaded or converted to "internal" UIDs.
//  Printing with this not set results in string UIDs being printed as
//  very large numbers -- for example, 9223372036854776016.
//
//  This is called whenever a gatekeeper store is opened, but you may
//  want to call it manually, if you have more than one store opened.
//
void     AS_UID_setGatekeeper(void *gkp);


//  Convert to/from a 64-bit integer that is guaranteed unique to this
//  assembly, regardless of if the UID came from an integer or from a
//  character string.
//
static
inline
uint64
AS_UID_toInteger(AS_UID uid) {
  uint64  uii = uid.isString;  //  Need to force it to a 64-bit value
  return((uii << 63) | uid.UID);
};

static
inline
AS_UID
AS_UID_fromInteger(uint64 uidn) {
  AS_UID   uid;
  uint64   mask = 1;
  mask <<= 63;

  uid.isString  = ((uidn & mask) != 0);
  uid.UID       = (uidn & ~mask);

  assert(uid.UID == (uidn & ~mask));
  assert(AS_UID_toInteger(uid) == uidn);
  return(uid);
};



//  Return a (static) character string containing the UID.  It
//  selects, round-robin, a static buffer from a set of 16.  See the
//  implementation.
//
char    *AS_UID_toString(AS_UID uid);


//  Return the undefined / uninitialized UID.
static
inline
AS_UID
AS_UID_undefined(void) {
  AS_UID   uid;
  uid.isString  = 0;
  uid.UID       = 0;
  return(uid);
}

//  True if the uid is undefined / uninitialized.
static
inline
int
AS_UID_isDefined(AS_UID uid) {
  return((uid.UID > 0) || (uid.isString == 1));
};

//  A special case, used only by terminator.  Please do not use.
static
inline
int
AS_UID_isString(AS_UID uid) {
  return(uid.isString);
}

//  Compare two UIDs.
static
inline
int
AS_UID_compare(AS_UID a, AS_UID b) {
  uint64  ai = AS_UID_toInteger(a);
  uint64  bi = AS_UID_toInteger(b);

  if (ai < bi)
    return(-1);
  if (ai > bi)
    return(1);
  return(0);
};



AS_UID  AS_UID_lookup(char *uidstr, char **nxtstr);
AS_UID  AS_UID_load(char *uidstr);


#endif  //  AS_UTL_UID_H
