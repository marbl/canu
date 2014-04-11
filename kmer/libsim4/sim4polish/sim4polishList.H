#ifndef SIM4_POLISH_LIST_H
#define SIM4_POLISH_LIST_H

#include "sim4polish.H"

//
//  A list of sim4polishes
//

class sim4polishList {
public:
  sim4polishList();
  sim4polishList(char const *filename);
  ~sim4polishList();

  void          push(sim4polish *p);
  void          remove(uint32 i);

  uint32        length(void)         { return(len); };

  sim4polish   *operator[](uint32 i) { if (i >= len) return(0L); return(list[i]); };
  sim4polish   *get(uint32 i)        { if (i >= len) return(0L); return(list[i]); };

  void          sortBycDNAIID(void);
  void          sortByGenomicIID(void);

  //  Removes polishes below the specified quality
  void          filterByQuality(uint32 minI, uint32 minC);

private:
  uint32        len;
  uint32        max;
  sim4polish  **list;
};


#endif  //  SIM4_POLISH_LIST_H
