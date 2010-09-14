#include "sim4polish.H" 

bool
sim4polish::s4p_makeForward(void) {
  if (_matchOrientation == SIM4_MATCH_FORWARD)
    return(false);

  for (u32bit e=0; e < _numExons; e++) {
    u32bit t = _estLen - _exons[e]._estFrom + 1;
    _exons[e]._estFrom = _estLen - _exons[e]._estTo + 1;
    _exons[e]._estTo = t;
  }

  _matchOrientation = SIM4_MATCH_FORWARD;

  return(true);
}


bool
sim4polish::s4p_makeReverse(void) {
  if (_matchOrientation == SIM4_MATCH_COMPLEMENT)
    return(false);

  for (u32bit e=0; e < _numExons; e++) {
    u32bit t = _estLen - _exons[e]._estFrom + 1;
    _exons[e]._estFrom = _estLen - _exons[e]._estTo + 1;
    _exons[e]._estTo = t;
  }

  _matchOrientation = SIM4_MATCH_COMPLEMENT;

  return(true);
}
