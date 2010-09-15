#include "sim4polish.H"

void
sim4polish::s4p_printPolish(FILE *O, sim4polishStyle style) {
  char *str = s4p_polishToString(style);

  fprintf(O, "%s", str);

  delete [] str;
}
