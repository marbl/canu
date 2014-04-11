#ifndef SIM4_POLISH_BUILDER_H
#define SIM4_POLISH_BUILDER_H

#include "sim4polish.H"

class sim4polishBuilder {
public:
  sim4polishBuilder();
  ~sim4polishBuilder();

  void          create(uint32 estid, uint32 estlen,
                       uint32 genid, uint32 genlo, uint32 genhi);

  void          setPolyTails(uint32 pa,
                             uint32 pt);

  void          setESTdefline(char *defline);
  void          setGENdefline(char *defline);

  void          setNumberOfMatches(uint32 nummatches, uint32 nummatchesN);
  void          setPercentIdentity(uint32 id);
  void          setMatchOrientation(char o);
  void          setStrandOrientation(char o);

  void          addExon(uint32 estlo, uint32 esthi,
                        uint32 genlo, uint32 genhi,
                        uint32 nummatches, uint32 nummatchesN, uint32 percentid,
                        char intronorientation);

  void          addExonAlignment(char *estalign,
                                 char *genalign);

  sim4polish   *release(void);
private:
  sim4polish      *it;

  uint32           exMax;  //  maximum number of exons available
  uint32           exPos;  //  next exon
  uint32           exAli;  //  next exon without alignment
  sim4polishExon **ex;
};

#endif  //  SIM4_POLISH_BUILDER_H
