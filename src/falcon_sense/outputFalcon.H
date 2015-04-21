

#ifndef OUTPUT_FALCON_H
#define OUTPUT_FALCON_H

#include "AS_global.H"
#include "gkStore.H"
#include "tgStore.H"

void
outputFalcon(gkStore      *gkpStore,
             tgTig        *tig,
             bool          trimToAlign,
             FILE         *F,
             gkReadData   *readData);


#endif  //  OUTPUT_FALCON_H
