#include "mcBucket.H"

bool
mcBucket::readBucket(void) {

  _bucketID++;

  if (_bucketID >= _lastBucketID)
    return(false);

  //  Read the number of items in this bucket
  //
  _items     = _IDX->getBits(32);
  _bitsRead += 32;

  if (_items > 0) {
    if (_items > _itemsMax) {
      delete _checks;
      delete _counts;

      _itemsMax = _items;

      _checks   = new u64bit [_itemsMax];
      _counts   = new u64bit [_itemsMax];
    }

    //  Read the checks and counts.
    //
    for (u32bit i=0; i<_items; i++) {
      u64bit v = _DAT->getBits(_chckBits + 1);

      _bitsRead += _chckBits + 1;

      if (v & _firstBit) {
        _checks[i] = v & _chckMask;
        _counts[i] = 1;
      } else {
        _checks[i] = v;

        u32bit shiftAmount = 0;
        bool   moreBits    = true;

        _counts[i] = 0;

        while (moreBits) {
          v = _DAT->getBits(_chckBits + 1);
          _bitsRead   += _chckBits + 1;
          _counts[i]  |=  (v & _chckMask) << shiftAmount;
          shiftAmount += _chckBits;

          moreBits = ((v & _firstBit) == u64bitZERO);
        }
      }
    }
  }

  return(true);
}


void
mcBucket::read(mcMer *mer) {
  u64bit   bucketIDdesired = mer->mer >> _chckBits;

  //fprintf(stderr, "Looking for bucketID 0x%016lx (mer = 0x%016lx)\n", bucketIDdesired, mer->mer);

  //  Do we have the bucket already?
  //
  if (_bucketID == bucketIDdesired)
    return;

  //  Nope, read it in.  We, sadly, have to read everything in between.
  //
  while (_bucketID < bucketIDdesired)
    readBucket();
}


void
mcBucket::scan(mcMer *mer) {
  u64bit   checkDesired = mer->mer & _chckMask;

  mer->count = 0;

  for (u64bit i=0; i<_items; i++) {
    if (checkDesired == _checks[i]) {
      mer->count = (u32bit)_counts[i];
      i = _items;
    }
  }
}
