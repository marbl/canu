

#ifndef ADJUST_OVERLAPS_H
#define ADJUST_OVERLAPS_H

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "clearRangeFile.H"


bool
adjustNormal(clearRangeFile  *iniClr,
             gkStore         *gkp,
             ovsOverlap      *ovl,
             uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
             uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend);

bool
adjustFlipped(clearRangeFile  *iniClr,
              gkStore         *gkp,
              ovsOverlap      *ovl,
              uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
              uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend);


#endif
