
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
static char CM_ID[] = "$Id: AS_UTL_HashCommon.c,v 1.2 2004-09-23 20:25:29 mcschatz Exp $";

#include "AS_UTL_HashCommon.h"


/* mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() takes 36 machine instructions, but only 18 cycles on a superscalar
  machine (like a Pentium or a Sparc).  No faster mixer seems to work,
  that's the result of my brute-force search.  There were about 2^68
  hashes to choose from.  I only tested about a billion of those.
*/

/* All one line for RedHat Linux 5.2 */
#define mix(a,b,c) { a -= b; a -= c; a ^= (c>>13);    b -= c; b -= a; b ^= (a<<8);   c -= a; c -= b; c ^= (b>>13);   a -= b; a -= c; a ^= (c>>12);    b -= c; b -= a; b ^= (a<<16);   c -= a; c -= b; c ^= (b>>5);   a -= b; a -= c; a ^= (c>>3);    b -= c; b -= a; b ^= (a<<10);   c -= a; c -= b; c ^= (b>>15); }

/* Hash_AS -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.
The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.
If you are hashing n strings (uint8 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = Hash_AS k[i], len[i], h);
By Bob Jenkins, 1996.  bob_jenkins@compuserve.com.  You may use this
code any way you wish, private, educational, or commercial.  It's free.
See http://ourworld.compuserve.com/homepages/bob_jenkins/evahash.htm
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
*/

uint32 Hash_AS( register uint8 *k,        /* the key */
	  register uint32  length,   /* the length of the key */
	  register uint32  initval)   /* the previous hash, or an arbitrary value */
{
   register uint32 a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9U;  /* the golden ratio; an arbitrary value */
   c = initval;         /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((uint32)k[1]<<8) +((uint32)k[2]<<16) +((uint32)k[3]<<24));
      b += (k[4] +((uint32)k[5]<<8) +((uint32)k[6]<<16) +((uint32)k[7]<<24));
      c += (k[8] +((uint32)k[9]<<8) +((uint32)k[10]<<16)+((uint32)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }
   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((uint32)k[10]<<24);
   case 10: c+=((uint32)k[9]<<16);
   case 9 : c+=((uint32)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((uint32)k[7]<<24);
   case 7 : b+=((uint32)k[6]<<16);
   case 6 : b+=((uint32)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((uint32)k[3]<<24);
   case 3 : a+=((uint32)k[2]<<16);
   case 2 : a+=((uint32)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}





