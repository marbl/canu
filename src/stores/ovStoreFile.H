
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
 *    Brian P. Walenz beginning on 2016-OCT-24
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_OVSTOREFILE_H
#define AS_OVSTOREFILE_H

#include "AS_global.H"
#include "gkStore.H"

#include "ovOverlap.H"


class ovStoreHistogram;


//  The default, no flags, is to open for normal overlaps, read only.  Normal overlaps mean they
//  have only the B id, i.e., they are in a fully built store.
//
//  Output of overlapper (input to store building) should be ovFileFullWrite.  The specialized
//  ovFileFullWriteNoCounts is used internally by store creation.
//
enum ovFileType {
  ovFileNormal              = 0,  //  Reading of b_id overlaps (aka store files)
  ovFileNormalWrite         = 1,  //  Writing of b_id overlaps
  ovFileFull                = 2,  //  Reading of a_id+b_id overlaps (aka dump files)
  ovFileFullWrite           = 3,  //  Writing of a_id+b_id overlaps
  ovFileFullWriteNoCounts   = 4   //  Writing of a_id+b_id overlaps, omitting the counts of olaps per read
};


class ovFile {
public:
  ovFile(gkStore     *gkpName,
         const char  *name,
         ovFileType   type = ovFileNormal,
         uint32       bufferSize = 1 * 1024 * 1024);
  ~ovFile();

  void    writeBuffer(bool force=false);
  void    writeOverlap(ovOverlap *overlap);
  void    writeOverlaps(ovOverlap *overlaps, uint64 overlapLen);

  void    readBuffer(void);
  bool    readOverlap(ovOverlap *overlap);
  uint64  readOverlaps(ovOverlap *overlaps, uint64 overlapMax);

  void    seekOverlap(off_t overlap);

  //  The size of an overlap record is 1 or 2 IDs + the size of a word times the number of words.
  uint64  recordSize(void) {
    return(sizeof(uint32) * ((_isNormal) ? 1 : 2) + sizeof(ovOverlapWORD) * ovOverlapNWORDS);
  };

  //  For use in conversion, force snappy compression.  By default, it is ENABLED, and we cannot
  //  read older ovb files.
#ifdef SNAPPY
  void    enableSnappy(bool enabled) {
    _useSnappy = enabled;
  };
#endif

  //  Move the stats in our histogram to the one supplied, and remove our data
  void    transferHistogram(ovStoreHistogram *copy);

private:
  gkStore                *_gkp;
  ovStoreHistogram       *_histogram;

  uint32                  _bufferLen;    //  length of valid data in the buffer
  uint32                  _bufferPos;    //  position the read is at in the buffer
  uint32                  _bufferMax;    //  allocated size of the buffer
  uint32                 *_buffer;

#ifdef SNAPPY
  size_t                  _snappyLen;
  char                   *_snappyBuffer;
#endif

  bool                    _isOutput;     //  if true, we can writeOverlap()
  bool                    _isSeekable;   //  if true, we can seekOverlap()
  bool                    _isNormal;     //  if true, 3 words per overlap, else 4
#ifdef SNAPPY
  bool                    _useSnappy;    //  if true, compress with snappy before writing
#endif

  compressedFileReader   *_reader;
  compressedFileWriter   *_writer;

  char                    _prefix[FILENAME_MAX];
  FILE                   *_file;
};


#endif  //  AS_OVSTOREFILE_H
