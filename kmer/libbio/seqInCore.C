
seqInCore::seqInCore(seqIID iid) {
  _idx       = iid;

  _headerLen = 0;
  _headerMax = 0;
  _header    = 0L;

  _seqLen = 0;
  _seqMax = 0;
  _seq    = 0L;
}

seqInCore::seqInCore(seqIID iid, char *hdr, u32bit hdrlen, char *seq, u32bit seqlen) {
  _idx       = iid;

  _headerLen = hdrlen;
  _headerMax = hdrlen;
  _header    = hdr;

  _seqLen = seqlen;
  _seqMax = seqlen;
  _seq    = seq;
}

seqInCore*
seqInCore::copy(void) {
  char *h = new char [_headerLen + 1];
  char *s = new char [_seqLen    + 1];

  memcpy(h, _header, _headerLen + 1);
  memcpy(s, _seq,    _seqLen    + 1);

  return(new seqInCore(_idx, h, _headerLen, s, _seqLen));
}

seqInCore::~seqInCore() {
  delete [] _header;
  delete [] _seq;
}
