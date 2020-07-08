
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "types.H"
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
