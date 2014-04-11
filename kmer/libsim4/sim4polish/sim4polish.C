#include "sim4polish.H" 

bool
sim4polish::s4p_makeForward(void) {
  if (_matchOrientation == SIM4_MATCH_FORWARD)
    return(false);

  for (uint32 e=0; e < _numExons; e++) {
    uint32 t = _estLen - _exons[e]._estFrom + 1;
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

  for (uint32 e=0; e < _numExons; e++) {
    uint32 t = _estLen - _exons[e]._estFrom + 1;
    _exons[e]._estFrom = _estLen - _exons[e]._estTo + 1;
    _exons[e]._estTo = t;
  }

  _matchOrientation = SIM4_MATCH_COMPLEMENT;

  return(true);
}
