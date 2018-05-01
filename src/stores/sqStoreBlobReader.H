
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
 *    src/stores/gkStoreBlobReader.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-MAR-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef GKSTOREBLOBREADER_H
#define GKSTOREBLOBREADER_H

#include "objectStore.H"

//  Manages access to blob data.  You need one of these per thread.
//
class sqStoreBlobReader {
public:
  sqStoreBlobReader() {
    _filesMax = 0;
    _files    = NULL;
  };

  ~sqStoreBlobReader() {
    for (uint32 ii=0; ii<_filesMax; ii++)
      AS_UTL_closeFile(_files[ii]);

    delete [] _files;
  };

  FILE      *getFile(const char *storePath, sqRead *read) {
    uint32  file = read->sqRead_mSegm();
    uint32  posn = read->sqRead_mByte();

    if (_filesMax == 0) {
      _filesMax = 8192;                   //  Limited in sqRead->H
      allocateArray(_files, _filesMax);
    }

    while (_filesMax <= file)
      resizeArray(_files, _filesMax, _filesMax, _filesMax * 2, resizeArray_copyData | resizeArray_clearNew);

    if (_files[file] == NULL) {
      char  N[FILENAME_MAX + 1];

      snprintf(N, FILENAME_MAX, "%s/blobs.%04u", storePath, file);

#pragma omp critical
      fetchFromObjectStore(N);   //  Fetch from object store, if needed and possible.

      _files[file] = AS_UTL_openInputFile(N);
    }

    AS_UTL_fseek(_files[file], posn, SEEK_SET);

    return(_files[file]);
  };


  uint32    _filesMax;
  FILE    **_files;      //  One file per blob file.
};


#endif  //  GKSTOREBLOBREADER_H
