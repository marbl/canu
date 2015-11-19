
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sffFile.H"
#include "bitOperations.H"

//  Lots of ths came from AS_GKP_sff.c



sffFile::sffFile() {
  clear();
}

sffFile::sffFile(const char *name) {

  clear();

  strcpy(_filename, name);

  _rb = new readBuffer(name);

  _rb->read(&_header, 31);

  if (_header.magic_number != 0x2e736666) {
    _header.swap_endianess           = 1;
    _header.magic_number             = uint32Swap(_header.magic_number);
    _header.index_offset             = uint64Swap(_header.index_offset);
    _header.index_length             = uint32Swap(_header.index_length);
    _header.number_of_reads          = uint32Swap(_header.number_of_reads);
    _header.header_length            = uint16Swap(_header.header_length);
    _header.key_length               = uint16Swap(_header.key_length);
    _header.number_of_flows_per_read = uint16Swap(_header.number_of_flows_per_read);
  }

  assert(_header.magic_number == 0x2e736666);
  assert(_header.number_of_flows_per_read < SFF_NUMBER_OF_FLOWS_MAX);
  assert(_header.key_length < SFF_KEY_SEQUENCE_MAX);

  _rb->read(_header.flow_chars,   sizeof(char) * _header.number_of_flows_per_read);
  _rb->read(_header.key_sequence, sizeof(char) * _header.key_length);

  _firstReadLocation = _header.header_length;

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.
  //
  if ((_header.index_length > 0) && (_header.index_offset == _header.header_length))
    _firstReadLocation += _header.index_length;

  //  Index
  //
  _index = new sffIndex [_header.number_of_reads];


  for (uint64 i=0; i<_header.number_of_reads; i++) {
    uint64  pos = _rb->tell();

    _rb->read(&_read, 16);

    if (_header.swap_endianess) {
      _read.read_header_length = uint16Swap(_read.read_header_length);
      _read.name_length        = uint16Swap(_read.name_length);
      _read.number_of_bases    = uint32Swap(_read.number_of_bases);
    }

    _index[i]._seqPos = pos;
    _index[i]._seqLen = _read.number_of_bases;
    _index[i]._namLen = _read.name_length;

    pos += _read.read_header_length;
    pos += sizeof(uint16) * _header.number_of_flows_per_read;
    pos += sizeof(uint8)  * _read.number_of_bases;
    pos += sizeof(char)   * _read.number_of_bases;
    pos += sizeof(uint8)  * _read.number_of_bases;

    pos += (_header.number_of_flows_per_read * sizeof(uint16) +
            _read.number_of_bases * sizeof(uint8) +
            _read.number_of_bases * sizeof(char) +
            _read.number_of_bases * sizeof(uint8)) % 8;

    _rb->seek(pos);
  }
  //
  //  Index

  _rb->seek(_firstReadLocation);

  _numberOfSequences = _header.number_of_reads;
}

sffFile::~sffFile() {
  delete    _rb;
  delete [] _index;
}

////////////////////////////////////////

seqFile *
sffFile::openFile(const char *name) {
  struct stat  st;

  //  Open the file, return if it matches the SFF magic_number.

  errno = 0;
  stat(name, &st);
  if (errno)
    return(0L);
  if ((st.st_mode & S_IFREG) == 0)
    return(0L);

  FILE *F = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "sffFile::openFile()- failed to open '%s': %s\n", name, strerror(errno));
    return(0L);
  }

  uint32  magic_number = 0;
  AS_UTL_safeRead(F, &magic_number, "sff magic_number", sizeof(uint32), 1);

  fclose(F);

  if ((magic_number == 0x2e736666) || (uint32Swap(magic_number) == 0x2e736666))
    return(new sffFile(name));

  return(0L);
}



bool
sffFile::getSequence(uint32 iid,
                     char *&h, uint32 &hLen, uint32 &hMax,
                     char *&s, uint32 &sLen, uint32 &sMax) {

  if (iid > _header.number_of_reads)
    return(false);

  memset(&_read, 0, sizeof(sffRead));

  _rb->seek(_index[iid]._seqPos);

  _rb->read(&_read, 16);

  if (_header.swap_endianess) {
    _read.read_header_length = uint16Swap(_read.read_header_length);
    _read.name_length        = uint16Swap(_read.name_length);
    _read.number_of_bases    = uint32Swap(_read.number_of_bases);
    _read.clip_quality_left  = uint16Swap(_read.clip_quality_left);
    _read.clip_quality_right = uint16Swap(_read.clip_quality_right);
    _read.clip_adapter_left  = uint16Swap(_read.clip_adapter_left);
    _read.clip_adapter_right = uint16Swap(_read.clip_adapter_right);
  }

  assert(_read.read_header_length < SFF_NAME_LENGTH_MAX);
  assert(_read.number_of_bases < SFF_NUMBER_OF_BASES_MAX);

  _rb->read(_read.name, sizeof(char) * _read.name_length);
  _read.name[_read.name_length] = 0;

  uint64 pos = _rb->tell();

  pos += _read.read_header_length;
  pos += sizeof(uint16) * _header.number_of_flows_per_read;
  pos += sizeof(uint8)  * _read.number_of_bases;

  _rb->seek(pos);

  _rb->read(_read.bases, sizeof(char) * _read.number_of_bases);
  _read.bases[_read.number_of_bases] = 0;

  return(true);
}





bool
sffFile::getSequence(uint32 iid,
                     uint32 bgn, uint32 end, char *s) {

  if (iid > _header.number_of_reads)
    return(false);

  //  Same as above, mostly.

  return(false);
}



void
sffFile::clear(void) {

  memset(_filename, 0, FILENAME_MAX);
  memset(_typename, 0, FILENAME_MAX);

  strcpy(_typename, "SFF");

  _numberOfSequences = 0;

  _rb = 0L;

  memset(&_header, 0, sizeof(sffHeader));
  memset(&_read,   0, sizeof(sffRead));

  _index = 0L;

  _firstReadLocation = 0;
  _readIID = 0;
}
