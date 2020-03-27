
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

inline
bool
ReadEnd::operator==(ReadEnd const that) const {
  return((readId() == that.readId()) && (read3p() == that.read3p()));
}

inline
bool
ReadEnd::operator!=(ReadEnd const that) const {
  return((readId() != that.readId()) || (read3p() != that.read3p()));
}

inline
bool
ReadEnd::operator<(ReadEnd const that) const {
  if (readId() != that.readId())
    return readId() < that.readId();
  else
    return read3p() < that.read3p();
}



////////////////////////////////////////


inline
void
BestEdgeOverlap::clear(void) {
  _id     = 0;
  _e3p    = false;
  _ahang  = 0;
  _bhang  = 0;
  _evalue = 0;
}

inline
void
BestEdgeOverlap::set(BAToverlap const &olap) {
  _id     = olap.b_iid;
  if (((olap.a_hang <= 0) && (olap.b_hang >= 0)) ||    //  If contained, _e3p just means
      ((olap.a_hang >= 0) && (olap.b_hang <= 0)))      //  the other read is flipped
    _e3p = olap.flipped;
  else                                                 //  Otherwise, means the
    _e3p = olap.BEndIs3prime();                        //    olap is to the 3' end

  _ahang  = olap.a_hang;
  _bhang  = olap.b_hang;
  _evalue = olap.evalue;
}

inline
void
BestEdgeOverlap::set(uint32 id, bool e3p, int32 ahang, int32 bhang, uint32 evalue) {
  _id     = id;
  _e3p    = e3p;
  _ahang  = ahang;
  _bhang  = bhang;
  _evalue = evalue;
}

inline
void
BestEdgeOverlap::setReverse(BestEdgeOverlap *that, uint32 readId, bool read3p) {
  _id     =  readId;
  _e3p    =  read3p;

  if (that->read3p() != read3p) {   //  Set hangs based on agreement of ends.
    _ahang  = -that->_ahang;        //    that->read3p() is the end this edge is coming out of
    _bhang  = -that->_bhang;        //    read3p         is the end this edge is going into
  } else {                          //  if they disagree, this is a 'normal' overlap (both reads
    _ahang  =  that->_bhang;        //  are the same orientation) and hangs negate.
    _bhang  =  that->_ahang;
  }

  _evalue =  that->_evalue;
}
