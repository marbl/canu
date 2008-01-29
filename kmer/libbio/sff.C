#include "sff.H"

//  Lots of ths came from AS_GKP_sff.c


#if 0

inline
u64bit
u64bitSwap(u64bit x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffLLU) | ((x <<  8) & 0xff00ff00ff00ff00LLU);
  x = ((x >> 16) & 0x0000ffff0000ffffLLU) | ((x << 16) & 0xffff0000ffff0000LLU);
  x = ((x >> 32) & 0x00000000ffffffffLLU) | ((x << 32) & 0xffffffff00000000LLU);
  return(x);
}

inline
u32bit
u32bitSwap(u32bit x) {
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);
  return(x);
}

inline
u16bit
u16bitSwap(u16bit x) {
  x = ((x >>  8) & 0x000000ff) | ((x <<  8) & 0x0000ff00);
  return(x);
}

#endif



#if 0

int
Load_SFF(FILE *sff) {

  sffHeader *h  = (sffHeader *)safe_calloc(sizeof(sffHeader), 1);
  sffRead   *r  = (sffRead   *)safe_calloc(sizeof(sffRead),   1);
  int        rn = 0;

  char       encodedsequence[AS_FRAG_MAX_LEN+1] = {0};

  readsff_header(sff, h);

  //  Construct a gkpLibraryRecord for this sff file.  Well, this is
  //  where we'd LIKE to do it, but since the sff doesn't give us any
  //  reasonable way to make a UID from the header, we defer until we
  //  get the first read.  Then, we use the read timestamp, hash and
  //  region to make a library.

  for (rn=0; rn < h->number_of_reads; rn++) {
    readsff_read(sff, h, r);

    //  Set clear ranges
    //
    //  These are base-based.  If either value is 0, that means the
    //  value was not computed.
    //
    //  We have a policy decision here.  If only one of the ranges is
    //  set, we can either ignore both, or set the unset one to the
    //  maximum.  We set it to the maximum.

    int  clq = _read.clip_quality_left;
    int  crq = _read.clip_quality_right;
    int  which;

    assert((_read.clip_quality_left == 0) || (h->key_length <= _read.clip_quality_left));
    assert((_read.clip_adapter_left == 0) || (h->key_length <= _read.clip_adapter_left));

    if (clq == 0)  clq = h->key_length + 1;
    if (crq == 0)  crq = _read.number_of_bases;

    for (which=0; which <= AS_READ_CLEAR_LATEST; which++) {
      gkf.clearBeg[which] = clq - h->key_length - 1;
      gkf.clearEnd[which] = crq - h->key_length;
    }

    if ((_read.clip_quality_left > 0) && (_read.clip_quality_right > 0)) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = _read.clip_quality_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = _read.clip_quality_right - h->key_length;
    } else if (_read.clip_quality_left > 0) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = _read.clip_quality_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = _read.number_of_bases - h->key_length;
    } else if (_read.clip_quality_right > 0) {
      gkf.hasQualityClear = 1;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = 0;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = _read.clip_quality_right - h->key_length;
    } else {
      gkf.hasQualityClear = 0;
      gkf.clearBeg[AS_READ_CLEAR_QLT] = 0;
      gkf.clearEnd[AS_READ_CLEAR_QLT] = 0;

    }

    if ((_read.clip_adapter_left > 0) && (_read.clip_adapter_right > 0)) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = _read.clip_adapter_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = _read.clip_adapter_right - h->key_length;
    } else if (_read.clip_adapter_left > 0) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = _read.clip_adapter_left  - h->key_length - 1;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = _read.number_of_bases - h->key_length;
    } else if (_read.clip_adapter_right > 0) {
      gkf.hasVectorClear  = 1;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = 0;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = _read.clip_adapter_right - h->key_length;
    } else {
      gkf.hasVectorClear = 0;
      gkf.clearBeg[AS_READ_CLEAR_VEC] = 0;
      gkf.clearEnd[AS_READ_CLEAR_VEC] = 0;
    }


    //  Now add the fragment to the store
    //
    gkf.readIID = getLastElemStore(gkpStore->frg) + 1;

    gkf.seqLen = strlen(_read.bases + h->key_length);
    gkf.hpsLen = 0;
    gkf.srcLen = strlen(_read.name);

    gkf.seqOffset = getLastElemStore(gkpStore->seq);
    gkf.qltOffset = getLastElemStore(gkpStore->qlt);
    gkf.hpsOffset = getLastElemStore(gkpStore->hps);
    gkf.srcOffset = getLastElemStore(gkpStore->src);

    setGatekeeperUIDtoIID(gkpStore, gkf.readUID, gkf.readIID, AS_IID_FRG);
    appendIndexStore(gkpStore->frg, &gkf);

    appendStringStore(gkpStore->seq, _read.bases + h->key_length, gkf.seqLen);

    encodeSequenceQuality(encodedsequence,
                          _read.bases + h->key_length,
                          _read.quality + h->key_length);
    appendStringStore(gkpStore->qlt, encodedsequence, gkf.seqLen);

    appendStringStore(gkpStore->hps, NULL,    0);
    appendStringStore(gkpStore->src, _read.name, gkf.srcLen);

    gkpStore->gkp.sffLoaded++;
  }

  fprintf(stderr, "Added %d 454 reads.\n", rn);

  return(0);
}

#endif


void
sffFile::clear(void) {

  //  I used to just memset(this, 0, sizeof(sffFile)) but that
  //  destroys the vtbl crud on i386 BSD.

  memset(_fileName, 0, sizeof(char) * FILENAME_MAX);

  _fileTimeStamp = 0;
  _file = 0L;

  memset(&_header, 0, sizeof(sffHeader));
  memset(&_read, 0, sizeof(sffRead));

  _index = 0L;

  _firstReadLocation = 0;
  _readIID = 0;
}

sffFile::sffFile() {
  clear();
}

sffFile::sffFile(const char *name) {

  clear();

  strcpy(_fileName, name);

  _fileTimeStamp = time(NULL);

  _file    = new readBuffer(name);

  _file->read(&_header, 31);

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

  _file->read(_header.flow_chars,   sizeof(char) * _header.number_of_flows_per_read);
  _file->read(_header.key_sequence, sizeof(char) * _header.key_length);

  _firstReadLocation = _header.header_length;

  //  The spec says the index might be here, however, all files I've
  //  seen have the index at the end of the file.
  //
  if ((_header.index_length > 0) && (_header.index_offset == _header.header_length))
    _firstReadLocation += _header.index_length;

  _file->seek(_firstReadLocation);
}

sffFile::~sffFile() {
  delete _file;
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

void
sffFile::openIndex(u32bit indextypetoload) {

  if (_index)
    return;

  _index = new sffIndex [_header.number_of_reads];

  rewind();

  for (u64bit i=0; i<_header.number_of_reads; i++) {
    u64bit  pos = _file->tell();

    _file->read(&_read, 16);

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

    _file->seek(pos);
  }
}

bool
sffFile::rewind(void) {
  return(_file->seek(_firstReadLocation));
}

bool
sffFile::eof(void) {
  return(_readIID == _header.number_of_reads);
}

bool
sffFile::getSequence(u32bit &hLen, char *&h,
                     u32bit &sLen, char *&s) {

  if (_curIID > _header.number_of_reads)
    return(false);

  u64bit pos = _file->tell();

  memset(&_read, 0, sizeof(sffRead));

  _file->read(&_read, 16);

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

  _file->read(_read.name, sizeof(char) * _read.name_length);
  _read.name[_read.name_length] = 0;

  pos += _read.read_header_length;
  pos += sizeof(u16bit) * _header.number_of_flows_per_read;
  pos += sizeof(u8bit)  * _read.number_of_bases;

  _file->seek(pos);

  _file->read(_read.bases, sizeof(char) * _read.number_of_bases);
  _read.bases[_read.number_of_bases] = 0;

  pos += sizeof(char)  * _read.number_of_bases;
  pos += sizeof(u8bit) * _read.number_of_bases;

  pos += (_header.number_of_flows_per_read * sizeof(u16bit) +
          _read.number_of_bases * sizeof(u8bit) +
          _read.number_of_bases * sizeof(char) +
          _read.number_of_bases * sizeof(u8bit)) % 8;

  _file->seek(pos);

  return(true);
}

bool
sffFile::find(seqIID  iid) {
  if (_index == 0L) {
    fprintf(stderr, "sffFile::find(IID)-- ERROR: '%s' not opened for random access.\n", _fileName);
    return(false);
  }

  if (iid >= _header.number_of_reads) {
    fprintf(stderr, "sffFile::find(IID)-- index of "u32bitFMT" too large for '%s' (only "u32bitFMT" sequences).n",
            iid, _fileName, _header.number_of_reads);
    return(false);
  }

  if (iid != _curIID) {
    _file->seek(_index[iid]._seqPos);
    _curIID = iid;
  }

  return(true);
}

bool
sffFile::find(char   *uid) {
  assert(0);
  return(false);
}

void
sffFile::printDescription(FILE *out, char *name) {
  assert(0);
}
