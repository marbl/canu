
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
 *    kmer/libutil/bitPackedFile.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-01
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAR-29 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-16 to 2014-APR-11
 *      are Copyright 2005-2008,2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "bitPackedFile.H"
#include "AS_UTL_fileIO.H"

//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <errno.h>
//#include <string.h>
#include <fcntl.h>


//  N.B. any read() / write() pair (either order) must have a seek (or
//  a fflush) in between.

bitPackedFile::bitPackedFile(char const *name, uint64 offset, bool forceTruncate) {

  _file = 0;
  _name = new char [strlen(name) + 1];
  strcpy(_name, name);

#ifdef WITH_BZIP2
  _bzFILE = 0L;
  _bzerr  = 0;
  _bzfile = 0L;
#endif

  _bfrmax = 1048576 / 8;
  _bfr    = new uint64 [_bfrmax];
  _pos    = uint64ZERO;
  _bit    = uint64ZERO;

  memset(_bfr, 0, sizeof(uint64) * _bfrmax);

  _inCore         = false;
  _bfrDirty       = false;
  _forceFirstLoad = false;
  _isReadOnly     = false;
  _isBzip2        = false;

  stat_seekInside   = uint64ZERO;
  stat_seekOutside  = uint64ZERO;
  stat_dirtyFlushes = uint64ZERO;

  file_offset        = 0;
  endianess_offset   = 0;
  endianess_flipped  = false;


  //  Try to open the original name -- we don't support compressed
  //  files for rewrite.  We just fail with a can't open message.
  //
  //  To get read/write and create we have to use open(2), as mode
  //  "r+" of fopen(3) will not create.  (Yes, but w+ does, sigh.)
  //
  if (forceTruncate) {
    errno = 0;
    _file = open(_name,
                 O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- failed to open and truncate '%s': %s\n",
              _name, strerror(errno)), exit(1);
  } else if (AS_UTL_fileExists(_name)) {
    errno = 0;
    _file = open(_name,
                 O_RDONLY | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- failed to open '%s': %s\n",
              _name, strerror(errno)), exit(1);
    _isReadOnly = true;
  } else {
    errno = 0;
    _file = open(_name,
                 O_RDWR | O_CREAT | O_LARGEFILE,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
    if (errno)
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- failed to open '%s': %s\n",
              _name, strerror(errno)), exit(1);
  }

  //  Move to the correct position in the file.
  //
  file_offset = offset;
  if (file_offset > 0)
    lseek(_file, file_offset, SEEK_SET);

  //  Deal with endianess.  We write out some bytes (or read back some bytes) to the start of
  //  the file, and then hide them from the user.
  //
  endianess_offset  = 32 + file_offset;
  endianess_flipped = false;

  char    t[16] = { 'b', 'i', 't', 'P', 'a', 'c', 'k', 'e', 'd', 'F', 'i', 'l', 'e', 0, 0, 1 };
  char    c[16] = { 0 };
  uint64  at = uint64NUMBER(0xdeadbeeffeeddada );
  uint64  bt = uint64NUMBER(0x0abeadedbabed8f8);
  uint64  ac = uint64NUMBER(0);
  uint64  bc = uint64NUMBER(0);
  size_t  nr = 0;

  errno = 0;
  nr += read(_file, c, sizeof(char) * 16);
  nr += read(_file, &ac, sizeof(uint64));
  nr += read(_file, &bc, sizeof(uint64));

  if (nr == 0) {
    //  Empty file!  Write the magic number and our endianess check.

    errno = 0;
    write(_file,  t,  sizeof(char) * 16);
    write(_file, &at, sizeof(uint64));
    write(_file, &bt, sizeof(uint64));
    if (errno)
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' failed to write the header: %s\n", _name, strerror(errno)), exit(1);

    return;
  }


  if ((c[0] == 'B') && (c[1] == 'Z') && (c[2] == 'h')) {
#ifdef WITH_BZIP2
    //  Looks like a bzip2 file!

    errno = 0;
    _bzFILE = fopen(_name, "r");
    if (errno) {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- failed to open bzip2 file '%s'\n", _name);
      exit(1);
    }

    _bzerr = 0;
    _bzfile = BZ2_bzReadOpen(&_bzerr, _bzFILE, 0, 0, 0L, 0);
    if ((_bzfile == 0L) || (_bzerr != BZ_OK)) {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- failed to init bzip2 file '%s'\n", _name);
      exit(1);
    }

    BZ2_bzRead(&_bzerr, _bzfile, c,   sizeof(char) * 16);
    BZ2_bzRead(&_bzerr, _bzfile, &ac, sizeof(uint64));
    BZ2_bzRead(&_bzerr, _bzfile, &bc, sizeof(uint64));

    //  XXX  should check bzerr!

    _isReadOnly  = true;
    _isBzip2 = true;
#else
    fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' looks like a bzip2 file, but bzip2 support not available!\n", _name);
    exit(1);
#endif
  }


  //  Check the magic number, decide on an endianess to use.
  //
  if (strncmp(t, c, 16) == 0) {
    if ((at == ac) && (bt == bc)) {
      endianess_flipped = false;
    } else if ((at == uint64Swap(ac)) && (bt == uint64Swap(bc))) {
      endianess_flipped = true;
    } else {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' looked like a bitPackedFile, but failed the endianess check, not opened.\n", _name);
      exit(1);
    }
  } else {
    fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' doesn't appear to be a bitPackedFile, not opened.\n", _name);
    fprintf(stderr, "bitPackedFile::bitPackedFile()-- found ");
    for (uint32 i=0; i<16; i++)
      fprintf(stderr, "%c", isascii(c[i]) ? c[i] : '.');
    fprintf(stderr, " at position "F_X64"\n", file_offset);
    exit(1);
  }

  _forceFirstLoad = true;
  seek(0);
}


bitPackedFile::~bitPackedFile() {
  flushDirty();
  delete [] _bfr;
  delete [] _name;
  close(_file);

#ifdef WITH_BZIP2
  if (_bzFILE)
    fclose(_bzFILE);

  if (_bzfile)
    BZ2_bzReadClose(&_bzerr, _bzfile);
#endif
}



//  If the page is dirty, flush it to disk
//
void
bitPackedFile::flushDirty(void) {

  if (_bfrDirty == false)
    return;

  if (_isReadOnly) {
    fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' is readonly, but is dirty!\n", _name);
    exit(1);
  }

  stat_dirtyFlushes++;

  errno = 0;
  lseek(_file, _pos * sizeof(uint64) + endianess_offset, SEEK_SET);
  if (errno) {
    fprintf(stderr, "bitPackedFile::seek()-- '%s' failed: %s\n",
            _name, strerror(errno));
    exit(1);
  }

  //  If we need to, flip all the words we are going to write
  //
  if (endianess_flipped)
    for (uint32 i=0; i<_bfrmax; i++)
      _bfr[i] = uint64Swap(_bfr[i]);

  //  We should only write bits up to _bit, the position we are
  //  currently at.  However, we don't know if the block is being
  //  flushed because we're totally finished with it, or because we
  //  are moving on to the next block.  If we're done with it, we
  //  want to flush the word that contains _bit, and if we're moving
  //  on to the next one, we'll flush that word again.  So, in
  //  either case, we flush the word that contains _bit.
  //
  errno = 0;
  write(_file, _bfr, sizeof(uint64) * _bfrmax);
  if (errno) {
    fprintf(stderr, "bitPackedFile::write()-- '%s' failed: %s\n",
            _name, strerror(errno));
    exit(1);
  }

  //  And then flip them back
  //
  if (endianess_flipped)
    for (uint32 i=0; i<_bfrmax; i++)
      _bfr[i] = uint64Swap(_bfr[i]);

  _bfrDirty = false;
}



void
bitPackedFile::seekBzip2(uint64 UNUSED(bitpos)) {

#ifdef WITH_BZIP2
  //  All we can do here is check that bitpos is
  //  a) in our current buffer
  //  b) would be in the next buffer once we read it

  uint64  newpos = bitpos >> 6;

  if (_pos + _bfrmax < newpos) {
    //  nope, not in the buffer -- we could probably handle this by just reading and
    //  discarding from the file until we get to the correct bitpos.
    fprintf(stderr, "bitPackedFile::seekBzip2()-- '%s' seek was not contiguous!\n", _name);
    exit(1);
  }

  //  Copy the remaining bits of the current buffer to the start.  Or
  //  not, if this is the first load.

  uint64  lastpos = _bit >> 6;            //  The word we are currently in
  uint64  lastlen = (_bfrmax - lastpos);  //  The number of words left in the buffer

  if (_forceFirstLoad == true) {
    lastpos = 0;
    lastlen = 0;
  } else {
    memcpy(_bfr, _bfr + lastpos, sizeof(uint64) * lastlen);
  }

  //  Update _bit and _pos -- lastlen is now the first invalid word
  //
  _bit  = bitpos & 0x3f;  //  64 * lastlen;
  _pos  = bitpos >> 6;

  //  Fill the buffer

  size_t  wordsread = 0;

  if (_bzfile) {
    _bzerr = 0;
    wordsread = BZ2_bzRead(&_bzerr, _bzfile, _bfr + lastlen, sizeof(uint64) * (_bfrmax - lastlen));
    if (_bzerr == BZ_STREAM_END) {
      //fprintf(stderr, "bitPackedFile::seekBzip2() file ended.\n");
      BZ2_bzReadClose(&_bzerr, _bzfile);
      fclose(_bzFILE);
      _bzfile = 0L;
      _bzFILE = 0L;
    } else if (_bzerr != BZ_OK) {
      fprintf(stderr, "bitPackedFile::seekBzip2() '%s' read failed.\n", _name);
      exit(1);
    }
  }

  //fprintf(stderr, "Filled buffer with %d words!\n", wordsread);

  //  Adjust to make wordsread be the index of the last word we actually read.
  //
  wordsread += lastlen;

  //  Flip all the words we just read, if needed
  //
  if (endianess_flipped)
    for (uint32 i=lastlen; i<wordsread; i++)
      _bfr[i] = uint64Swap(_bfr[i]);

  //  Clear any words that we didn't read (supposedly, because we hit
  //  EOF).
  //
  while (wordsread < _bfrmax)
    _bfr[wordsread++] = uint64ZERO;
#else
  fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s'\n", _name);
  fprintf(stderr, "bitPackedFile::bitPackedFile()-- bzip2 support not present, but still tried to read it??\n");
  exit(1);
#endif
}



void
bitPackedFile::seekNormal(uint64 bitpos) {

  if (_inCore) {
    fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' is in core, but still needed to seek??\n",
            _name);
    exit(1);
  }

  //  Somewhat of a gross hack to allow sequential access backwards.
  //
  //  If the new position (bitpos >> 6) is just before the old
  //  position (_pos), assume that we are being accessed iteratively
  //  backwards and load a full buffer so that the position we want to
  //  access is at the end.
  //
  //  Easy to think of bone-headed ways to break this (e.g., seek to
  //  the second element in a structure, access the first, then access
  //  the third).  Not so easy to think of a logical reason someone
  //  would want to do that.
  //
  if (((bitpos >> 6) < _pos) && (_pos <= (bitpos >> 6) + 32)) {
    _pos = bitpos >> 6;
    if (_pos > _bfrmax)
      _pos = _pos - _bfrmax + 32;
    else
      _pos = 0;
  } else {
    _pos = bitpos >> 6;
  }

  _bit = bitpos - (_pos << 6);


  errno = 0;
  lseek(_file, _pos * 8 + endianess_offset, SEEK_SET);
  if (errno) {
    fprintf(stderr, "bitPackedFile::seekNormal() '%s' seek to pos="F_U64" failed: %s\n",
            _name,
            _pos * 8 + endianess_offset, strerror(errno));
    exit(1);
  }

  errno = 0;
  size_t wordsread = read(_file, _bfr, sizeof(uint64) * _bfrmax);
  if (errno) {
    fprintf(stderr, "bitPackedFile::seekNormal() '%s' read of "F_U64" bytes failed': %s\n",
            _name,
            sizeof(uint64) * _bfrmax,
            strerror(errno));
    exit(1);
  }

  //  Flip all the words we just read, if needed
  //
  if (endianess_flipped)
    for (uint32 i=0; i<wordsread; i++)
      _bfr[i] = uint64Swap(_bfr[i]);

  //  Clear any words that we didn't read (supposedly, because we hit
  //  EOF).
  //
  while (wordsread < _bfrmax)
    _bfr[wordsread++] = uint64ZERO;
}





//  Seeks to bitposition pos in the file, reads in a new block.
//
void
bitPackedFile::seek(uint64 bitpos) {

  //  If we are seeking to somewhere in the current block, don't do a
  //  real seek, just move our position within the block.
  //
  if (_forceFirstLoad == false) {
    uint64 np = bitpos >> 6;

    if ((_pos <= np) && (np <= _pos + _bfrmax - 32)) {
      _bit = bitpos - (_pos << 6);
      stat_seekInside++;
      //fprintf(stderr, "SEEK INSIDE to _bit="F_U64"\n", _bit);
      return;
    }
  }

  if (_inCore) {
    fprintf(stderr, "bitPackedFile::seek()-- file '%s' is in core, but still needed to seek??\n",
            _name);
    exit(1);
  }

  stat_seekOutside++;

  flushDirty();

  if (_isBzip2)
    seekBzip2(bitpos);
  else
    seekNormal(bitpos);

  _forceFirstLoad = false;

  //fprintf(stderr, "SEEK OUTSIDE to _pos="F_U64" _bit="F_U64"\n", _pos, _bit);
}




uint64
bitPackedFile::loadInCore(void) {
  struct stat  sb;

  //  Convert this disk-based, read/write bitPackedFile to memory-based read-only.

  flushDirty();

  fstat(_file, &sb);

  //  The extra 1024 words is to keep seek() from attempting to grab
  //  the next block (there isn't a next block, we've got it all!)
  //  when we're near the end of this block.  We just make the block
  //  a little bigger than it really is.

  delete [] _bfr;

  _bfrmax = sb.st_size / 8 + 1024;
  _bfr    = new uint64 [_bfrmax];
  _pos    = 0;
  _bit    = 0;

  //  Tada!  All we need to do now is load the block!

  _forceFirstLoad = true;

  seek(0);

  _inCore = true;

  return(_bfrmax * 8);
}
