#ifndef SIM4_POLISH_BUILDER_H
#define SIM4_POLISH_BUILDER_H

#include "sim4polish.h"


//  A list of sim4polishes
//
class sim4polishList {
public:
  sim4polishList();
  ~sim4polishList();
  void          push(sim4polish *p);
  u32bit        length(void)         { return(len); };

  //  pretty useless as I return a pointer to a sim4polishList...
  sim4polish   *operator[](u32bit i) { if (i >= len) return(0L); return(list[i]); };
  sim4polish   *get(u32bit i)        { if (i >= len) return(0L); return(list[i]); };
private:
  u32bit        len;
  u32bit        max;
  sim4polish  **list;
};



class sim4polishBuilder {
public:
  sim4polishBuilder();
  ~sim4polishBuilder();

  void          create(u32bit estid, u32bit estlen,
                       u32bit genid, u32bit genlo, u32bit genhi);

  void          setPolyTails(u32bit pa,
                             u32bit pt);

  void          setESTdefline(char *defline);
  void          setGENdefline(char *defline);

  void          setNumberOfMatches(u32bit nummatches, u32bit nummatchesN);
  void          setPercentIdentity(u32bit id);
  void          setMatchOrientation(char o);
  void          setStrandOrientation(char o);

  void          addExon(u32bit estlo, u32bit esthi,
                        u32bit genlo, u32bit genhi,
                        u32bit nummatches, u32bit nummatchesN, u32bit percentid,
                        char intronorientation);

  void          addExonAlignment(char *estalign,
                                 char *genalign);

  sim4polish   *release(void);
private:
  sim4polish      *it;

  u32bit           exMax;  //  maximum number of exons available
  u32bit           exPos;  //  next exon
  u32bit           exAli;  //  next exon without alignment
  sim4polishExon **ex;
};

#endif  //  SIM4_POLISH_BUILDER_H
