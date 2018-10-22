
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "files.H"

//  Basic routines to fetch and stash files from an object store.
//  Most of this is done in the executive, but low level fetching
//  of sqStore and ovStore data is done here.
//
//  NOTE that this function is limited in its ability to fetch files.
//  It will ONLY work with seqStore and ovlStore data files:
//     seqStore/blobs.*
//     ovlStore/0000<000>
//
//  Returns false if the file was not fetched (either no object store
//  in use, or the file existed already), true if it was fetched.
//
bool   fetchFromObjectStore(char *filename);
