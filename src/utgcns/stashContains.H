
#ifndef STASH_CONTAINS_H
#define STASH_CONTAINS_H

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"


#include <map>
#include <algorithm>


class readLength {
public:
  uint32    idx;
  int32     len;

  bool operator<(const readLength &that) const {
    return(len < that.len);
  };
};


class savedChildren {
public:
  savedChildren(tgTig *tig) {
    childrenLen = tig->_childrenLen;
    childrenMax = tig->_childrenMax;
    children    = tig->_children;
  };

  uint32      childrenLen;
  uint32      childrenMax;
  tgPosition *children;
};


savedChildren *
stashContains(tgTig  *tig,
              double  maxCov,
              bool    beVerbose = true);

void
unstashContains(tgTig         *tig,
                savedChildren *saved);


#endif // STASH_CONTAINS_H
