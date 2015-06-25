#include "Display_Alignment.H"
#include "gkStore.H"

#define  DISPLAY_WIDTH   250

//  Show (to  stdout ) the alignment encoded in  delta [0 .. (delta_ct - 1)]
//  between strings  a [0 .. (a_len - 1)]  and  b [0 .. (b_len - 1)] .
//  Capitialize  a  characters for positions at and after  capitalize_start .

void
Display_Alignment(char  *a,    int32   a_len,
                  char  *b,    int32   b_len,
                  int32 *delta,
                  int32  delta_ct,
                  int32  capitalize_start) {

  char  *top = new char [AS_MAX_READLEN + 1];
  char  *bot = new char [AS_MAX_READLEN + 1];

  int32 i          = 0;
  int32 j          = 0;
  int32 top_len    = 0;
  int32 bot_len    = 0;

  for (int32 k = 0;  k < delta_ct;  k++) {
    for (int32 m = 1;  m < abs (delta[k]);  m++) {
      if (i >= capitalize_start)
        top[top_len++] = toupper (a[i++]);
      else
        top[top_len++] = a[i++];
      j++;
    }

    if (delta[k] < 0) {
      top[top_len++] = '-';
      j++;

    } else {
      if (i >= capitalize_start)
        top[top_len++] = toupper (a[i++]);
      else
        top[top_len++] = a[i++];
    }
  }

  while (i < a_len && j < b_len) {
    if (i >= capitalize_start)
      top[top_len++] = toupper (a[i++]);
    else
      top[top_len++] = a[i++];
    j++;
  }
  top[top_len] = '\0';


  i = 0;
  j = 0;

  for (int32 k = 0;  k < delta_ct;  k++) {
    for (int32 m = 1;  m < abs (delta[k]);  m++) {
      bot[bot_len++] = b[j++];
      i++;
    }

    if (delta[k] > 0) {
      bot[bot_len++] = '-';
      i++;
    } else {
      bot[bot_len++] = b[j++];
    }
  }

  while (j < b_len && i < a_len) {
    bot[bot_len++] = b[j++];
    i++;
  }

  bot[bot_len] = '\0';


  for (int32 i=0; (i < top_len) || (i < bot_len); i += DISPLAY_WIDTH) {
    printf("%d\n", i);
    printf("A: ");

    for (int32 j=0;  (j < DISPLAY_WIDTH) && (i+j < top_len);  j++)
      putchar(top[i+j]);

    printf("\n");
    printf("B: ");

    for (int32 j=0; (j < DISPLAY_WIDTH) && (i+j < bot_len); j++)
      putchar (bot[i+j]);

    printf("\n");
    printf("   ");

    for (int32 j=0; (j<DISPLAY_WIDTH) && (i+j < bot_len) && (i+j < top_len); j++)
      if ((top[i+j] != ' ') && (bot[i+j] != ' ') && (tolower(top[i+j]) != tolower(bot[i+j])))
        putchar('^');
      else
        putchar(' ');

    putchar('\n');
  }

  delete [] top;
  delete [] bot;
}
