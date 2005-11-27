#include "util++.H"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>


//
//  N.B. any read/write pair (either way) must have a seek (or a fflush) in between.
//

bitPackedFile::bitPackedFile(char const *name, u64bit offset) {

  //  Try to open the original name -- we don't support compressed
  //  files for rewrite.  We just fail with a can't open message.
  //
  //  To get read/write and create we have to use open(2), as mode
  //  "r+" of fopen(3) will not create.
  //
  errno = 0;
  _file = open(name,
               O_RDWR | O_CREAT | O_LARGEFILE,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
  int errno_fopen = errno;

  //  If that fails, try to open the compressed version
  //
  if (errno_fopen) {
    size_t  l = strlen(name);

    char *command = new char [l + 64];
    sprintf(command, "%s.bz2", name);

    if (fileExists(command)) {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- Can't use compressed bitPackedFiles here.\n");
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' not opened.\n", command);
    } else {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- Couldn't open bitPackedFile.\n");
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- %s: %s\n", name, strerror(errno_fopen));
    }

    exit(1);
  }

  _bfrmax = 4096;
  _bfrmax = 1048576 / 8;
  _bfr    = new u64bit [_bfrmax];
  _pos    = u64bitZERO;
  _bit    = u64bitZERO;

  _bfrDirty = false;

  stat_seekInside   = u64bitZERO;
  stat_seekOutside  = u64bitZERO;
  stat_dirtyFlushes = u64bitZERO;

  //  Move to the correct position in the file.
  //
  file_offset = offset;
  lseek(_file, file_offset, SEEK_SET);

  //  Deal with endianess.  We write out some bytes (or read back some bytes) to the start of
  //  the file, and then hide them from the user.
  //
  endianess_offset  = 32 + file_offset;
  endianess_flipped = false;

  char    t[16] = { 'b', 'i', 't', 'P', 'a', 'c', 'k', 'e', 'd', 'F', 'i', 'l', 'e', 0, 0, 1 };
  char    c[16] = { 0 };
  u64bit  a = u64bitNUMBER(0xdeadbeeffeeddada );
  u64bit  b = u64bitNUMBER(0x0abeadedbabed8f8);

  size_t num = read(_file, c, sizeof(char) * 16);
  if (num == 0) {
    //  Empty file!  Write the magic number and our endianess check.
    //
    errno = 0;
    write(_file,  t, sizeof(char) * 16);
    write(_file, &a, sizeof(u64bit));
    write(_file, &b, sizeof(u64bit));
    if (errno)
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' failed to write the header: %s\n", name, strerror(errno)), exit(1);
  } else {
    //  We read something.  Make sure it's correct, then check endianess.
    //
    if (strncmp(t, c, 16) == 0) {
      u64bit ac, bc;
      read(_file, &ac, sizeof(u64bit));
      read(_file, &bc, sizeof(u64bit));

      if ((a == ac) && (b == bc)) {
        endianess_flipped = false;
      } else if ((a == u64bitSwap(ac)) && (b == u64bitSwap(bc))) {
        endianess_flipped = true;
      } else {
        fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' looked like a bitPackedFile, but failed the endianess check, not opened.\n", name);
        exit(1);
      }
    } else {
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- '%s' doesn't appear to be a bitPackedFile, not opened.\n", name);
      fprintf(stderr, "bitPackedFile::bitPackedFile()-- found ");
      for (u32bit i=0; i<16; i++)
        fprintf(stderr, "%c", isascii(c[i]) ? c[i] : '.');
      fprintf(stderr, " at position "u64bitHEX"\n", file_offset);
      exit(1);
    }
  }

  seek(0, true);
}


bitPackedFile::~bitPackedFile() {
  flushDirty();
  delete [] _bfr;
  close(_file);
}




//  If the page is dirty, flush it to disk
//
void
bitPackedFile::flushDirty(void) {

  if (_bfrDirty) {
    stat_dirtyFlushes++;

    errno = 0;
    lseek(_file, _pos * sizeof(u64bit) + endianess_offset, SEEK_SET);
    if (errno) {
      fprintf(stderr, "bitPackedFile::seek() failed: %s\n", strerror(errno));
      exit(1);
    }

    //  We should only write bits up to _bit, the position we are
    //  currently at.  However, we don't know if the block is being
    //  flushed because we're totally finished with it, or because we
    //  are moving on to the next block.  If we're done with it, we
    //  want to flush the word that contains _bit, and if we're moving
    //  on to the next one, we'll flush that word again.  So, in
    //  either case, we flush the word that contains _bit.
    //

    //  If we need to , flip all the words we are going to write
    //
    if (endianess_flipped)
      for (u32bit i=0; i<_bfrmax; i++)
        _bfr[i] = u64bitSwap(_bfr[i]);

    errno = 0;
    write(_file, _bfr, sizeof(u64bit) * _bfrmax);
    if (errno) {
      fprintf(stderr, "bitPackedFile::write() failed: %s\n", strerror(errno));
      exit(1);
    }

    //  And then flip them back
    //
    if (endianess_flipped)
      for (u32bit i=0; i<_bfrmax; i++)
        _bfr[i] = u64bitSwap(_bfr[i]);

    _bfrDirty = false;
  }
}





//  Seeks to bitposition pos in the file, reads in a new block.
//
void
bitPackedFile::seek(u64bit bitpos, bool forceLoad) {

  //  If we are seeking to somewhere in the current block, don't do a
  //  real seek, just move our position within the block.
  //
  if (forceLoad == false) {
    u64bit np = bitpos >> 6;
    
    if ((_pos <= np) && (np <= _pos + _bfrmax - 32)) {
      _bit = bitpos - (_pos << 6);
      stat_seekInside++;
      //fprintf(stderr, "SEEK INSIDE to _bit="u64bitFMT"\n", _bit);
      return;
    }
  }

  stat_seekOutside++;

  flushDirty();

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
  if (((bitpos >> 6) < _pos) && ((bitpos >> 6) + 32 >= _pos)) {
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
    fprintf(stderr, "bitPackedFile::seek() failed: %s\n", strerror(errno));
    exit(1);
  }

  errno = 0;
  size_t wordsread = read(_file, _bfr, sizeof(u64bit) * _bfrmax);
  if (errno) {
    fprintf(stderr, "bitPackedFile::bitPackedFile got %s\n", strerror(errno));
    exit(1);
  }

  //  Flip all the words we just read, if needed
  //
  if (endianess_flipped)
    for (u32bit i=0; i<wordsread; i++)
      _bfr[i] = u64bitSwap(_bfr[i]);

  //  Clear any words that we didn't read (supposedly, because we hit
  //  EOF).
  //
  while (wordsread < _bfrmax)
    _bfr[wordsread++] = u64bitZERO;

  //fprintf(stderr, "SEEK OUTSIDE to _pos="u64bitFMT" _bit="u64bitFMT"\n", _pos, _bit);
}




u64bit
bitPackedFile::getBits(u32bit siz) {
  sync();
  u64bit ret = getDecodedValue(_bfr, _bit, siz);
  _bit += siz;
  return(ret);
}

void
bitPackedFile::putBits(u64bit bits, u32bit siz) {
  sync();
  setDecodedValue(_bfr, _bit, siz, bits);
  _bit += siz;
  _bfrDirty = true;
}

u64bit
bitPackedFile::getNumber(void) {
  sync();
  u64bit siz = 0;
  u64bit ret = getFibonacciEncodedNumber(_bfr, _bit, &siz);
  _bit += siz;
  return(ret);
}

void
bitPackedFile::putNumber(u64bit val) {
  sync();
  u64bit siz = 0;
  setFibonacciEncodedNumber(_bfr, _bit, &siz, val);
  _bit += siz;
  _bfrDirty = true;
}
