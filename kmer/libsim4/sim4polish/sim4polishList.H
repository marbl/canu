#ifndef SIM4_POLISH_LIST_H
#define SIM4_POLISH_LIST_H

#include "sim4polish.h"

//
//  A list of sim4polishes
//

class sim4polishList {
public:
  sim4polishList();
  sim4polishList(char const *filename);
  ~sim4polishList();

  void          push(sim4polish *p);
  void          remove(u32bit i);

  u32bit        length(void)         { return(len); };

  sim4polish   *operator[](u32bit i) { if (i >= len) return(0L); return(list[i]); };
  sim4polish   *get(u32bit i)        { if (i >= len) return(0L); return(list[i]); };

  void          sortBycDNAIID(void);
  void          sortByGenomicIID(void);

  //  Removes polishes below the specified quality
  void          filterByQuality(u32bit minI, u32bit minC);

private:
  u32bit        len;
  u32bit        max;
  sim4polish  **list;
};


#endif  //  SIM4_POLISH_LIST_H
