#ifndef SIM4_POLISH_BUILDER_H
#define SIM4_POLISH_BUILDER_H

#include "sim4polish.h"

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
