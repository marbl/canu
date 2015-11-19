
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/utgcns/libNDalign/prefixEditDistance-allocateMoreSpace.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUL-20 to 2015-AUG-05
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "NDalgorithm.H"


//  Allocate another block of 64mb for edits

//  Needs to be at least:
//       52,432 to handle 40% error at  64k overlap
//      104,860 to handle 80% error at  64k overlap
//      209,718 to handle 40% error at 256k overlap
//      419,434 to handle 80% error at 256k overlap
//    3,355,446 to handle 40% error at   4m overlap
//    6,710,890 to handle 80% error at   4m overlap
//  Bigger means we can assign more than one Edit_Array[] in one allocation.

uint32  EDIT_SPACE_SIZE  = 1 * 1024 * 1024;

bool
NDalgorithm::allocateMoreEditSpace(void) {

  //  Determine the last allocated block, and the last assigned block

  int32  b = 0;  //  Last edit array assigned
  int32  e = 0;  //  Last edit array assigned more space
  int32  a = 0;  //  Last allocated block

  while (Edit_Array_Lazy[b] != NULL)
    b++;

  while (Edit_Space_Lazy[a] != NULL)
    a++;

  //  Fill in the edit space array.  Well, not quite yet.  First, decide the minimum size.
  //
  //  Element [0] can access from [-2] to [2] = 5 elements.
  //  Element [1] can access from [-3] to [3] = 7 elements.
  //
  //  Element [e] can access from [-2-e] to [2+e] = 5 + e * 2 elements
  //
  //  So, our offset for this new block needs to put [e][0] at offset...

  int32 Offset = 2 + b;
  int32 Del    = 6 + b * 2;
  int32 Size   = EDIT_SPACE_SIZE;

  while (Size < Offset + Del)
    Size *= 2;

  //  Allocate another block

  Edit_Space_Lazy[a] = new pedEdit [Size];

  //  And, now, fill in the edit space array.

  e = b;

  while ((Offset + Del < Size) &&
         (e < Edit_Space_Max)) {
    Edit_Array_Lazy[e++] = Edit_Space_Lazy[a] + Offset;

    Offset += Del;
    Del    += 2;
  }

  if (e == b) {
    fprintf(stderr, "Allocate_More_Edit_Space()-- ERROR: couldn't allocate enough space for even one more entry!  e=%d\n", e);
    return(false);
  }
  assert(e != b);

  return(true);
  //fprintf(stderr, "WorkArea %d allocates space %d of size %d for array %d through %d\n", thread_id, a, Size, b, e-1);
}

