#include "sffFile.H"

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
    _header.magic_number             = u32bitSwap(_header.magic_number);
    _header.index_offset             = u64bitSwap(_header.index_offset);
    _header.index_length             = u32bitSwap(_header.index_length);
    _header.number_of_reads          = u32bitSwap(_header.number_of_reads);
    _header.header_length            = u16bitSwap(_header.header_length);
    _header.key_length               = u16bitSwap(_header.key_length);
    _header.number_of_flows_per_read = u16bitSwap(_header.number_of_flows_per_read);
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


  for (u64bit i=0; i<_header.number_of_reads; i++) {
    u64bit  pos = _rb->tell();

    _rb->read(&_read, 16);

    if (_header.swap_endianess) {
      _read.read_header_length = u16bitSwap(_read.read_header_length);
      _read.name_length        = u16bitSwap(_read.name_length);
      _read.number_of_bases    = u32bitSwap(_read.number_of_bases);
    }

    _index[i]._seqPos = pos;
    _index[i]._seqLen = _read.number_of_bases;
    _index[i]._namLen = _read.name_length;

    pos += _read.read_header_length;
    pos += sizeof(u16bit) * _header.number_of_flows_per_read;
    pos += sizeof(u8bit)  * _read.number_of_bases;
    pos += sizeof(char)   * _read.number_of_bases;
    pos += sizeof(u8bit)  * _read.number_of_bases;

    pos += (_header.number_of_flows_per_read * sizeof(u16bit) +
            _read.number_of_bases * sizeof(u8bit) +
            _read.number_of_bases * sizeof(char) +
            _read.number_of_bases * sizeof(u8bit)) % 8;

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

  u32bit  magic_number = 0;
  safeRead(fileno(F), &magic_number, "sff magic_number", sizeof(u32bit));

  fclose(F);

  if ((magic_number == 0x2e736666) || (u32bitSwap(magic_number) == 0x2e736666))
    return(new sffFile(name));

  return(0L);
}



bool
sffFile::getSequence(u32bit iid,
                     char *&h, u32bit &hLen, u32bit &hMax,
                     char *&s, u32bit &sLen, u32bit &sMax) {

  if (iid > _header.number_of_reads)
    return(false);

  memset(&_read, 0, sizeof(sffRead));

  _rb->seek(_index[iid]._seqPos);

  _rb->read(&_read, 16);

  if (_header.swap_endianess) {
    _read.read_header_length = u16bitSwap(_read.read_header_length);
    _read.name_length        = u16bitSwap(_read.name_length);
    _read.number_of_bases    = u32bitSwap(_read.number_of_bases);
    _read.clip_quality_left  = u16bitSwap(_read.clip_quality_left);
    _read.clip_quality_right = u16bitSwap(_read.clip_quality_right);
    _read.clip_adapter_left  = u16bitSwap(_read.clip_adapter_left);
    _read.clip_adapter_right = u16bitSwap(_read.clip_adapter_right);
  }

  assert(_read.read_header_length < SFF_NAME_LENGTH_MAX);
  assert(_read.number_of_bases < SFF_NUMBER_OF_BASES_MAX);

  _rb->read(_read.name, sizeof(char) * _read.name_length);
  _read.name[_read.name_length] = 0;

  u64bit pos = _rb->tell();

  pos += _read.read_header_length;
  pos += sizeof(u16bit) * _header.number_of_flows_per_read;
  pos += sizeof(u8bit)  * _read.number_of_bases;

  _rb->seek(pos);

  _rb->read(_read.bases, sizeof(char) * _read.number_of_bases);
  _read.bases[_read.number_of_bases] = 0;

  return(true);
}





bool
sffFile::getSequence(u32bit iid,
                     u32bit bgn, u32bit end, char *s) {

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
