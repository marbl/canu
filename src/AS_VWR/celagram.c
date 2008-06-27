
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
/* $Id: celagram.c,v 1.7 2008-06-27 06:29:21 brianwalenz Exp $ */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>

#include "button.h"

#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (1)
#endif


#ifndef TRUE
#define TRUE  (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

static int semilog = FALSE;
static double semilog_multiplier = 100;
static double semilog_base = 10.;
static double semilogFac = 43.42945;
static int ordinate_remap( int x) {
  // Apply semi-log transformation
  double xy = ((double)x > 1. ? semilogFac*log((double)x) : 0. );
  return (int) xy;
}


/* Color parameters */

#undef SWITCH_BACKGROUND_COLOR

#ifndef SWITCH_BACKGROUND_COLOR
#define TEXTCOLOR  255,255,255   /* Text */
#define BACKCOLOR  0,0,0         /* Background */
#define SELCOLOR   200,200,100   /* General selection hilite color */
#define BARCOLOR   255,255,0     /* Color for the histogram bars */
#else // SWITCH_BACKGROUND_COLOR
#define TEXTCOLOR  0,0,0         /* Text */
#define BACKCOLOR  255,255,255   /* Background */
#define SELCOLOR   100,100,100   /* General selection hilite color */
#define BARCOLOR   TEXTCOLOR     /* Color for the histogram bars */
#endif // SWITCH_BACKGROUND_COLOR

#define KSELCOLOR  125,125,70    /* Selection color for scrollbar buttons */
#define KNOBCOLOR  50,50,225     /* Scrollbar button & slider color */

/* Color values for above parameters */

long BackGround, TextColor;

long SelectColor, KnobColor, KselColor;

char *HistoTitle;

int RoundDown(int val, int modulus)
{ int xvl;

  xvl = val - (val % modulus);
  if (val < 0 && val % modulus != 0)
    xvl -= modulus;
  return (xvl);
}

int RoundUp(int val, int modulus)
{ if (val % modulus != 0)
    val = RoundDown(val,modulus) + modulus;
  return (val);
}

/* Histogram specific code */

typedef struct {
    int   numpts;
    long   sumofpts;
    int  *datapts;
    float stddev;
    int   bordermax;
} HistoPacket;

#define MAX_TITLE_LENGTH 60

static int HISTO_COMPARE(const void *l, const void *r)
{ int x, y;
  x = *((int *) l);
  y = *((int *) r);
  return (x - y);
}


int GetColumnValue(FILE * file, int * histo_cell, int column)
{
  char line[2048];
  char * lineptr = line;

  // get line from file
  if(fgets(line, 2047, file) == NULL)
    return -1;

  // get column-th number from line
  // get past any initial whitespace
  while(isspace(*lineptr)) lineptr++;
  for(; *lineptr != '\0' && column != 1; lineptr++)
  {
    if(!isspace(*lineptr))
    {
      while(!isspace(*lineptr)) lineptr++;
      column--;
    }
  }

  if(column == 1 && *lineptr != '\0')
  {
    *histo_cell = atoi(lineptr);
    return 1;
  }

  return 0;
}


HistoPacket *ReadHistogram(FILE *file, char **title, int column)
{ HistoPacket *hist;
  int i, *dps;
  long sum;
  double mean, dev;
  int *Histogram, HistoLen, HistoMax;

  if(*title == NULL) {
    char buffer[2048];
    if (fgets(buffer,2047,file) == NULL)
      { fprintf(stderr,"File is empty !\n");
      exit (1);
      }
    buffer[MAX_TITLE_LENGTH] = 0;
    buffer[strlen(buffer)-1] = 0;
    if (strlen(buffer) > MAX_TITLE_LENGTH+1)
      { fprintf(stderr,
                "Title %s is too long (%d > %d chars)\n",
                buffer, (int)strlen(buffer),MAX_TITLE_LENGTH);
      exit (1);
      }
    if ((*title = strdup(buffer)) == NULL)
      { fprintf(stderr,"Out of memory\n");
      exit (1);
      }
  }

  HistoLen  = 0;
  HistoMax  = 100000;
  Histogram = (int *) malloc(sizeof(int)*HistoMax);
  while (GetColumnValue(file,Histogram+HistoLen,column) == 1)
    { HistoLen += 1;
      if (HistoLen >= HistoMax)
        { HistoMax  = (int)(1.4*HistoMax) + 10000;
          Histogram = (int *) realloc(Histogram,sizeof(int)*HistoMax);
          if (Histogram == NULL)
            { fprintf(stderr,"Out of memory\n");
              exit (1);
            }
        }
    }
  if (!feof(file))
    { fprintf(stderr,"Did not recognize number\n");
      exit (1);
    }

  if (HistoLen > 0)
    qsort(Histogram,HistoLen,sizeof(int),HISTO_COMPARE);

  hist = (HistoPacket *) malloc(sizeof(HistoPacket));
  dps  = (int *) malloc(sizeof(int)*(HistoLen+1));
  sum  = 0;
  for (i = 0; i < HistoLen; i++)
    sum += (dps[i] = Histogram[i]);
  if (HistoLen > 0)
    mean = (1.*sum) / HistoLen;
  else
    mean = 0.;
  dev  = 0.;
  for (i = 0; i < HistoLen; i++)
    dev += (dps[i] - mean) * (dps[i] - mean);
  hist->numpts   = HistoLen;
  hist->sumofpts = sum;
  hist->datapts  = dps;
  if (HistoLen > 0)
    hist->stddev = sqrt(dev/HistoLen);
  else
    hist->stddev = 0.;

  return (hist);
}

void FreeHistogram(HistoPacket *hist)
{ free(hist->datapts);
  free(hist);
}

void GetHistoDims(HistoPacket *hist, int *min, int *max)
{ *min = hist->datapts[0];
  *max = hist->datapts[hist->numpts-1];
}


static int TickRounder(double value, int *order, int *zeros)
{ int z, l, o, i;
  double v;

  if (value < 1.) value = 1.;
  z = (int)(log10(value));
  v = value / pow(10.,1.*z);
  if (v > 5)
    { l = 1; z += 1; }
  else if (v > 2)
    l = 5;
  else if (v > 1)
    l = 2;
  else
    l = 1;
  o = 1;
  for (i = 1; i <= z; i++)
    o *= 10;
  *order = o;
  *zeros = z;
  return (l*o);
}

static void AdjustTicks(int *offset, int *lead, int *order, int *zeros)
{ int i, digit;

  digit = (*lead)/(*order);
  if (*offset > 0)
    for (i = *offset; i > 0; i--)
      if (digit == 5)
        { digit = 1; *order *= 10; *zeros += 1; }
      else if (digit == 2)
        digit = 5;
      else
        digit = 2;
  else
    for (i = *offset; i < 0; i++)
      if (digit == 1)
        { if (*zeros == 0)
            { *offset -= i; break; }
          else
            { digit = 5; *order /= 10; *zeros -= 1; }
        }
      else if (digit == 2)
        digit = 1;
      else
        digit = 2;
  *lead = digit * (*order);
}

#define HWIDE    6
#define HWMAX  100
#define HHIGH   20
#define HBORDER 10
#define HTEXTB  18
#define HTITLEB 24

void SetHistoBorderMax(HistoPacket *hist, int minw, int minh)
{ int xmin, xmax, hmax;
  int bmin, bmax;
  int lead, order, zeros;
  int border, llead;
  static char  label[100];

  xmin = hist->datapts[0];
  xmax = hist->datapts[hist->numpts-1];
  border = 0;

  lead = TickRounder((1.*(xmax-xmin)) / ((minw-(2*HBORDER))/HWMAX),
                       &order,&zeros);

readjust:
  llead = lead;
  bmin  = RoundDown(xmin,lead);
  bmax  = RoundDown(xmax,lead) + lead;

  { int i, k, b;

    hmax = 0;
    i = 0;
    while (i < hist->numpts && hist->datapts[i] < bmin)
      i += 1;
    for (k = bmin+lead; k <= bmax; k += lead)
      { b = 0;
        while (i < hist->numpts && hist->datapts[i] < k)
          { b += 1; i += 1; }
        if (b > hmax) hmax = b;
      }
    if (hmax == 0) hmax = 1;
  }

  lead = TickRounder( (1.*hmax) / ((minh-(HBORDER+HTEXTB))/HHIGH),
                       &order, &zeros);
  hmax = RoundUp(hmax,lead);

  { int h, b;

    sprintf(label,"%d",hmax);
    mt_string_size(label,&border,&h,&b);
  }

  lead = TickRounder((1.*(xmax-xmin)) / ((minw-(2*HBORDER - border))/HWMAX),
                       &order,&zeros);

  if (llead != lead) goto readjust;

  hist->bordermax = border;
}


void DrawHistogram(MT_OBJECT *canvas, HistoPacket *hist,
                   int xmin, int xmax, int *bsize)
{ int xl, xh, yl, yh;
  int lst, bot;

  static char  label[100];
  static char *letter[] = { "", "K", "M", "G", "T" };
  FILE *fout = stdout;

  mt_get_extent(canvas,&xl,&xh,&yl,&yh);
  xh -= (xl-1);
  yh -= (yl-1);

  { int xzeros, xlead, xorder;
    int yzeros, ylead, yorder;
    int bmin, bmax, hmax;
    int lborder, hwide;

    lborder = hist->bordermax;

    xlead = TickRounder((1.*(xmax-xmin)) / ((xh-(2*HBORDER+lborder))/HWIDE),
                        &xorder,&xzeros);
    AdjustTicks(bsize,&xlead,&xorder,&xzeros);

    while (1)
      { int adjust;

        bmin  = RoundDown(xmin,xlead);
        bmax  = RoundDown(xmax,xlead) + xlead;
        hwide = (xh - (2*HBORDER+lborder)) / ((bmax-bmin) / xlead);

        if (hwide < 2)
          { adjust = 1;
            AdjustTicks(&adjust,&xlead,&xorder,&xzeros);
            *bsize += 1;
          }
        else if (hwide > HWMAX)
          { adjust = -1;
            AdjustTicks(&adjust,&xlead,&xorder,&xzeros);
            if (adjust < 0)
              *bsize += -1;
            else
              break;
          }
        else
          break;
      }

    { int i, k, b;

      hmax = 0;
      i = 0;
      while (i < hist->numpts && hist->datapts[i] < bmin)
        i += 1;
      for (k = bmin+xlead; k <= bmax; k += xlead)
        { b = 0;
          while (i < hist->numpts && hist->datapts[i] < k)
            { b += 1; i += 1; }
          if (b > hmax) hmax = b;
        }
      if (hmax == 0) hmax = 1;
    }
    if(semilog) {
      hmax = ordinate_remap(1.*hmax);
    }

    ylead = TickRounder( (1.*hmax) / ((yh-(HTITLEB+HTEXTB))/HHIGH),
                         &yorder, &yzeros);
    hmax = RoundUp(hmax,ylead);

    { int w, h ,b;

      sprintf(label,"%d",hmax);
      mt_string_size(label,&w,&h,&b);
      if (bmin <= hist->datapts[0] && bmax > hist->datapts[hist->numpts-1])
        lst = (w+4) + (xh - (hwide*(bmax-bmin)/xlead + w + 4))/2;
      else
        lst = lborder + HBORDER;
      bot = yh - HTEXTB;
    }

#if 0
    if(dump_file != NULL) {
      FILE * fdump = fopen(dump_file,"w");
      for(i=0; i < hist->numpts; i++) {
        fprintf(fdump,"%d %d\n", i,(hist->datapts[i]));
      }
      fclose(fdump);
    }
#endif

    { int i, k, b, y, lft;
      double yfact;
      mt_set_color(mt_get_color(BARCOLOR));
      i = 0;
      yfact = (1.*yh - (HTITLEB+HTEXTB))/hmax;
      while (i < hist->numpts && hist->datapts[i] < bmin)
        i += 1;
      lft = lst;
      for (k = bmin+xlead; k <= bmax; k += xlead)
        { b = 0;
          while (i < hist->numpts && hist->datapts[i] < k)
            { b += 1; i += 1; }
          if(semilog) {
#ifdef SEMILOG_DEBUG
            if( b > 0) { fprintf(fout,"%d %d %d %d\n",
                                 k, (hist->datapts[i]), b, ordinate_remap(b));}
#endif
            b = ordinate_remap(b);
          }
          y = yfact * b;
          if (y == 0 && b != 0) y = 1;
          mt_draw_rect(canvas,lft,bot-y,hwide-1,y);
          lft += hwide;
        }
    }
    fflush(fout);

    { int ticks, lft, plead, qlead;
      int adjust, qmax;
      double yfact;

      mt_set_color(TextColor);

      lft   = lst;
      mt_draw_line(canvas,lft-1,bot,lft+(hwide*(bmax-bmin)/xlead)-1,bot,0,0);

      adjust = 1;
      plead  = xlead;

      while (1)
        { int w, h, b;

          qmax = RoundDown(bmax,xlead);
          if (bmin < 0)
            { if (-10*RoundUp(bmin,xlead) > qmax) qmax = RoundUp(bmin,xlead); }
          if (xzeros % 3 == 2)
            sprintf(label,"%d.%d%s",qmax/(xorder*10),(qmax/xorder)%10,
                                    letter[xzeros/3+1]);
          else if (xzeros % 3 == 1)
            sprintf(label,"%d%s",qmax/(xorder/10),letter[xzeros/3]);
          else
            sprintf(label,"%d%s",qmax/xorder,letter[xzeros/3]);

          mt_string_size(label,&w,&h,&b);
          if (w+4 >= hwide*(xlead/plead))
            AdjustTicks(&adjust,&xlead,&xorder,&xzeros);
          else
            break;
        }
      if (xlead/plead < 10)
        qlead = bmax+1;
      else if (xlead/xorder == 5)
        qlead = xlead/5;
      else
        qlead = xlead/2;

      ticks = bmin;
      while (ticks <= bmax)
        { int w, h, b;

          if (ticks % xlead == 0)
            { if (xzeros % 3 == 2)
                sprintf(label,"%d.%d%s",ticks/(xorder*10),(ticks/xorder)%10,
                                        letter[xzeros/3+1]);
              else if (xzeros % 3 == 1)
                sprintf(label,"%d%s",ticks/(xorder/10),letter[xzeros/3]);
              else
                sprintf(label,"%d%s",ticks/xorder,letter[xzeros/3]);
              mt_string_size(label,&w,&h,&b);
              if (lft >= w/2 && lft <= xh-(w/2))
                { mt_draw_text(canvas,lft-(w/2),bot+16,label);
                  mt_draw_line(canvas,lft,bot+1,lft,bot+4,2,0);
                }
            }
          if (ticks % qlead == 0)
            mt_draw_line(canvas,lft-1,bot+1,lft-1,bot+4,1,0);
          ticks += plead;
          lft += hwide;
        }

      mt_set_color(mt_get_color(BARCOLOR));
      mt_draw_title(canvas,lst+5,HTITLEB-8,HistoTitle);

      mt_set_color(TextColor);
      mt_draw_line(canvas,lst-1,bot,lst-1,HTITLEB,1,0);
      yfact = (1.*yh - (HTITLEB+HTEXTB))/hmax;
      ticks = 0;
      while (ticks <= hmax)
        { int y, w, h, b;

          y = bot - yfact*ticks;
          sprintf(label,"%d",ticks);
          mt_string_size(label,&w,&h,&b);
          mt_draw_text(canvas,lst-(w+4),y+h/2-1,label);
          mt_draw_line(canvas,lst,y,lst-4,y,1,0);
          ticks += ylead;
        }
    }
  }
}

static int printwidth(long value)
{ int width;

  if (value == 0)
    width = 1;
  else if (value < 0)
    width = log10(-1.*value) + 2;
  else
    width = log10(1.*value) + 1;
  return (width);
}

void PrintHistogram(HistoPacket *hist, char *title)
{ int max, alt;
  int width1, width2;

  if(title != NULL) { printf("\nHistogram of: %s\n\n",title);}

  if (hist->numpts == 0)
    { printf("  No data points!\n");
      return;
    }

  width1 = printwidth(hist->sumofpts);

  max = (1.*hist->sumofpts)/hist->numpts;
  if (max < 0) max *= -10.;
  alt = hist->stddev;
  if (alt < 0) alt *= -10.;
  if (alt > max) max = alt;
  alt = hist->datapts[hist->numpts/2];
  if (alt < 0) alt *= -10.;
  if (alt > max) max = alt;
  if (max < 1.)
    width2 = 1;
  else
    width2 = log10(1.*max) + 1;

  printf("  Min = %*d  ",width1,hist->datapts[0]);
  printf("  Mean     = %*.2f  ",width2+3,(1.*hist->sumofpts)/hist->numpts);
  printf("  # of pts.  = %d\n",hist->numpts);

  printf("  Max = %*d  ",width1,hist->datapts[hist->numpts-1]);
  printf("  Median   = %*d\n",width2,hist->datapts[hist->numpts/2]);

  printf("  Sum = %*ld  ",width1,hist->sumofpts);

  printf("  Std Dev. = %*.2f\n",width2+3,hist->stddev);
}

int DrawHistoInfo(MT_OBJECT *bar, HistoPacket *hist)
{ char label[100];
  int border, max, width, alt;
  int w, h, b;

  if (hist->numpts == 0) return (0);

  mt_set_color(mt_black());
  border = 5;

  width = printwidth(hist->sumofpts);

  sprintf(label,"Min = %*d",width,hist->datapts[0]);
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  sprintf(label,"Max = %*d",width,hist->datapts[hist->numpts-1]);
  if (bar != NULL) mt_draw_text(bar,border,23,label);
  sprintf(label,"Sum = %*ld",width,hist->sumofpts);
  if (bar != NULL) mt_draw_text(bar,border,34,label);

  mt_string_size(label,&w,&h,&b);
  border += w+10;

  max = (1.*hist->sumofpts)/hist->numpts;
  if (max < 0) max *= -10.;
  alt = hist->stddev;
  if (alt < 0) alt *= -10.;
  if (alt > max) max = alt;
  alt = hist->datapts[hist->numpts/2];
  if (alt < 0) alt *= -10.;
  if (alt > max) max = alt;
  if (max < 1.)
    width = 1;
  else
    width = log10(1.*max) + 1;

  sprintf(label,"Mean     = %*.2f",width+3,(1.*hist->sumofpts)/hist->numpts);
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  sprintf(label,"Median   = %*d",width,hist->datapts[hist->numpts/2]);
  if (bar != NULL) mt_draw_text(bar,border,34,label);
  sprintf(label,"Std Dev. = %*.2f",width+3,hist->stddev);
  if (bar != NULL) mt_draw_text(bar,border,23,label);

  mt_string_size(label,&w,&h,&b);
  border += w+10;

  sprintf(label,"# of pts.  = %d",hist->numpts);
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  mt_string_size(label,&w,&h,&b);
  max = w;

  border += max+5;

  return (border);
}

/* Geometry parameters */

#define BORDER      2   /* Border relief for all objects */

#define HISTWIDE  300   /* Histogram window dimensions */
#define HISTHIGH  250
#define HISTBUTTON 50
#define SCROLLWIDE 16

/* Scroll bar object */

typedef struct stag {
    MT_OBJECT *bar;                       /* Central slider */
    MT_OBJECT *lft;                       /* Left/Up button */
    MT_OBJECT *rgt;                       /* Right/Down button */
    void (*updater)(long, struct stag *); /* Event handler routine */
    long      data;                       /* Value to pass to action routine */
} MY_SCROLL;

/* Arrow bit mask creation routines */

static MT_OBJECT *lftarrow = NULL;
static MT_OBJECT *rgtarrow;
static MT_OBJECT *uparrow  = NULL;
static MT_OBJECT *dwnarrow;

void CreateHorizontalArrows(void)
{ int size, size2, i;

  if (lftarrow != NULL) return;

  size = SCROLLWIDE-2*(BORDER+2);
  size2 = (size-1)/2;
  lftarrow = mt_new_bitmap(size2+1,size);
  mt_set_color(KnobColor);
  mt_draw_rect(lftarrow,0,0,size2+1,size);
  mt_set_color(mt_white());
  for (i = 0; i <= size2; i++)
    mt_draw_line(lftarrow,size2-i,i,size2-i,(size-1)-i,1,0);

  rgtarrow = mt_new_bitmap(size2+1,size);
  mt_set_color(KnobColor);
  mt_draw_rect(rgtarrow,0,0,size2+1,size);
  mt_set_color(mt_white());
  for (i = 0; i <= size2; i++)
    mt_draw_line(rgtarrow,i,i,i,(size-1)-i,1,0);
}

void CreateVerticalArrows(void)
{ int size, size2, i;

  if (uparrow != NULL) return;

  size = SCROLLWIDE-2*(BORDER+2);
  size2 = (size-1)/2;
  uparrow = mt_new_bitmap(size,size2+1);
  mt_set_color(KnobColor);
  mt_draw_rect(uparrow,0,0,size,size2+1);
  mt_set_color(mt_white());
  for (i = 0; i <= size2; i++)
    mt_draw_line(uparrow,i,size2-i,(size-1)-i,size2-i,1,0);

  dwnarrow = mt_new_bitmap(size,size2+1);
  mt_set_color(KnobColor);
  mt_draw_rect(dwnarrow,0,0,size,size2+1);
  mt_set_color(mt_white());
  for (i = 0; i <= size2; i++)
    mt_draw_line(dwnarrow,i,i,(size-1)-i,i,1,0);
}

/* Scroll Event handlers (for all scrollbars) */

void scrollupdate0(long d, MY_SCROLL *s)  /* Non-slave update */
{ mt_draw((MT_OBJECT *) d); }

void scrollupdate1(long d, MY_SCROLL *s)  /* Slave update: linear */
{ int val0, wid0, min0, max0;
  int val1, wid1, min1, max1;

  val0 = mt_get_scrollbar(s[-1].bar,&wid0,&min0,&max0);
  max0 += wid0;
  val1 = mt_get_scrollbar(s->bar,&wid1,&min1,&max1);
  wid0 = val1;
  if (wid0 > max0-min0) wid0 = max0-min0;
  if (val0+wid0 > max0) val0 = max0-wid0;
  if (val0 < min0) val0 = min0;
  mt_set_scrollbar(s[-1].bar,val0,wid0,min0,max0-wid0);
  mt_draw(s[-1].bar);
  mt_draw((MT_OBJECT *) d);
}

int scrollleft(MT_OBJECT *o, long d)
{ int val, wid, min, max;
  MY_SCROLL *s;

  s = (MY_SCROLL *) d;
  val = mt_get_scrollbar(s->bar,&wid,&min,&max);
  val -= wid/2;
  if (val < min) val = min;
  mt_set_scrollbar(s->bar,val,wid,min,max);
  mt_draw(s->bar);
  s->updater(s->data,s);
  return (0);
}

int scrollright(MT_OBJECT *o, long d)
{ int val, wid, min, max;
  MY_SCROLL *s;

  s = (MY_SCROLL *) d;
  val = mt_get_scrollbar(s->bar,&wid,&min,&max);
  val += wid/2;
  if (val > max) val = max;
  mt_set_scrollbar(s->bar,val,wid,min,max);
  mt_draw(s->bar);
  s->updater(s->data,s);
  return (0);
}

int scrollmove(MT_OBJECT *o, long d)
{ MY_SCROLL *s;
  s = (MY_SCROLL *) d;
  s->updater(s->data,s);
  return (0);
}

/* Scroll bar creation routines */

void MakeHorizontalScroll(MY_SCROLL *s, int x, int y, int w,
                          void (*u)(long, MY_SCROLL *), long d)
{ MT_OBJECT *bit;

  CreateHorizontalArrows();

  s->bar = bit = mt_new_scrollbar(x+SCROLLWIDE,y,
                                  w-2*SCROLLWIDE,SCROLLWIDE,0,1,50,10,0,100);
  mt_set_xscale(bit,0.,1.);
  mt_set_yscale(bit,1.,0.);
  mt_scrollbar_colors(bit,KnobColor,-1,-1,-1,SelectColor);
  mt_set_callback(bit,scrollmove,(long) s);

  s->lft = bit = mt_new_button(x,y,SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(bit,0.,0.);
  mt_set_yscale(bit,1.,0.);
  mt_button_colors(bit,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(bit,mt_bitmap_label(lftarrow));
  mt_set_callback(bit,scrollleft,(long) s);

  s->rgt = bit = mt_new_button(x+w-SCROLLWIDE,y,
                               SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(bit,1.,0.);
  mt_set_yscale(bit,1.,0.);
  mt_button_colors(bit,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(bit,mt_bitmap_label(rgtarrow));
  mt_set_callback(bit,scrollright,(long) s);

  s->updater = u;
  s->data    = d;
}

void MakeVerticalScroll(MY_SCROLL *s, int x, int y, int h,
                        void (*u)(long, MY_SCROLL *), long d)
{ MT_OBJECT *bit;

  CreateVerticalArrows();

  s->bar = bit = mt_new_scrollbar(x,y+SCROLLWIDE,
                                  SCROLLWIDE,h-2*SCROLLWIDE,1,1,50,10,0,100);
  mt_set_xscale(bit,1.,0.);
  mt_set_yscale(bit,0.,1.);
  mt_scrollbar_colors(bit,KnobColor,-1,-1,-1,SelectColor);
  mt_set_callback(bit,scrollmove,(long) s);

  s->lft = bit = mt_new_button(x,y,SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(bit,1.,0.);
  mt_set_yscale(bit,0.,0.);
  mt_button_colors(bit,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(bit,mt_bitmap_label(uparrow));
  mt_set_callback(bit,scrollleft,(long) s);

  s->rgt = bit = mt_new_button(x,y+h-SCROLLWIDE,
                               SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(bit,1.,0.);
  mt_set_yscale(bit,1.,0.);
  mt_button_colors(bit,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(bit,mt_bitmap_label(dwnarrow));
  mt_set_callback(bit,scrollright,(long) s);

  s->updater = u;
  s->data    = d;
}

typedef struct {
    HistoPacket *hist;
    MT_OBJECT   *darea;
    MT_OBJECT   *window;
    MY_SCROLL   sbar[2];
    int         bindelta;
} HistoLoad;

int histo_handler(MT_OBJECT *o, long d)
{ HistoLoad *p;
  int val0, wid0, min, max;

  p = (HistoLoad *) d;
  val0 = mt_get_scrollbar(p->sbar[0].bar,&wid0,&min,&max);
  DrawHistogram(p->darea,p->hist,val0,val0+wid0,&p->bindelta);
  return (0);
}

int histo_info(MT_OBJECT *o, long d)
{ HistoLoad *p;

  p = (HistoLoad *) d;
  DrawHistoInfo(o,p->hist);
  return (0);
}

int cbutton(MT_OBJECT *o, long d)
{ HistoLoad *p;

  p = (HistoLoad *) d;
  FreeHistogram(p->hist);
  mt_free(p->window);
  mt_shutdown();
  free(p);
  return (0);
}

int BinButtonUp(MT_OBJECT *o, long d)
{ HistoLoad *p;

  p = (HistoLoad *) d;
  p->bindelta += 1;
  mt_draw(p->darea);
  return (0);
}

int BinButtonDown(MT_OBJECT *o, long d)
{ HistoLoad *p;

  p = (HistoLoad *) d;
  p->bindelta -= 1;
  mt_draw(p->darea);
  return (0);
}

int DrawBinButton(MT_OBJECT *o, long d)
{ int wlab, wplus, wminus;
  int xl, xh, yl, yh;
  int h, b;

  mt_string_size("Bin Size",&wlab,&h,&b);
  mt_string_size("+",&wplus,&h,&b);
  mt_string_size("-",&wminus,&h,&b);
  mt_get_extent(o,&xl,&xh,&yl,&yh);
  mt_set_color(mt_black());
  mt_draw_text(o,((xh-xl+1) - wlab)/2,h-b,"Bin Size");
  mt_draw_text(o,3,(yh-yl)-b,"-");
  mt_draw_text(o,(xh-xl+1)-(wplus+2),(yh-yl)-b,"+");
  return (0);
}

void GetBinButtonDims(int *x, int *y)
{ int wlab, wplus, wminus;
  int h, b;

  mt_string_size("Bin Size",&wlab,&h,&b);
  mt_string_size("+",&wplus,&h,&b);
  mt_string_size("-",&wminus,&h,&b);
  if (wlab+4 > wplus + wminus + 2*SCROLLWIDE+8)
    *x = wlab + 4;
  else
    *x = wplus + wminus + 2*SCROLLWIDE + 8;
  *y = h + SCROLLWIDE;
}

void CreateBinButton(HistoLoad *p, int binwide, int binhigh)
{ MT_OBJECT *obj;

  CreateHorizontalArrows();

  obj = mt_new_button(binwide/2-SCROLLWIDE,binhigh-SCROLLWIDE,
                      SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(obj,.5,0.);
  mt_set_yscale(obj,0.,0.);
  mt_button_colors(obj,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(obj,mt_bitmap_label(lftarrow));
  mt_set_callback(obj,BinButtonDown,(long) p);

  obj = mt_new_button(binwide/2,binhigh-SCROLLWIDE,
                      SCROLLWIDE,SCROLLWIDE,1,CLICK,NULL);
  mt_set_xscale(obj,.5,0.);
  mt_set_yscale(obj,0.,0.);
  mt_button_colors(obj,KnobColor,KselColor,-1,-1,KselColor,-1);
  mt_set_label(obj,mt_bitmap_label(rgtarrow));
  mt_set_callback(obj,BinButtonUp,(long) p);

  p->bindelta = 0;
}


void Usage(void)
{
  fprintf(stderr,
          "Usage: <-c column> <-t title> <-p> <-s> <-S float,float> Celagram file ...\n"
          "-c column uses the specified column for count data\n"
          "-t string uses the string as the histogram title\n"
          "-s has each accumulated count replaced by 100*log(count)/log(10)\n"
          "-S m,b has each accumulated count replaced by m*log(count)/log(base)\n"
          );
  exit(1);
}


int main(int argc, char *argv[])
{ MT_OBJECT *obj = NULL;
  int i;
  FILE *file = NULL;
  HistoLoad *packet = NULL;
  int textwide, histwide, binwide, binhigh;
  char *title = NULL;
  HistoPacket *hist = NULL;
  int column = 1;
  int first_file = 1;
  int printTheHistograms = FALSE;

  if(argc < 2) {
    Usage();
  }

  { /* Parse the argument list using "man 3 getopt". */
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg &&
	   ((ch = getopt(argc, argv,
			 "c:pst:S:"
                         )) != EOF)) {
      switch(ch) {
      case 'c':
        column = atoi(optarg);
        break;
      case 'p':
        printTheHistograms = TRUE;
        break;
      case 'S':
        sscanf(optarg,"%lf,%lf",&semilog_multiplier,&semilog_base);
      case 's':
        semilog = TRUE;
        semilogFac = semilog_multiplier/log(semilog_base);
        break;
      case 't':
        title = strdup(optarg);
        break;
      default:
        fprintf(stderr,"Unknown option %c\n", ch);
        exit(1);
        break;
      }
    }
  }


  first_file = optind;
  for (i = first_file; i < argc; i++)
    { file = fopen(argv[i],"r");
      if (file == NULL)
        { fprintf(stderr,"Can't open file %s\n",argv[i]);
          exit (1);
        }
      hist = ReadHistogram(file,&title,column);
      fclose(file);
      if (argc - first_file > 1 || printTheHistograms)
        { PrintHistogram(hist,title);
          FreeHistogram(hist);
        }
      // free(title); title = NULL;
    }

  if (argc - first_file > 1 ||  printTheHistograms) exit (0);

  mt_startup();

  mt_current_border(BORDER);
  SelectColor = mt_get_color(SELCOLOR);
  KnobColor   = mt_get_color(KNOBCOLOR);
  KselColor   = mt_get_color(KSELCOLOR);
  TextColor   = mt_get_color(TEXTCOLOR);
  BackGround  = mt_get_color(BACKCOLOR);

  packet = (HistoLoad *) malloc(sizeof(HistoLoad));
  packet->hist  = hist;

  textwide = DrawHistoInfo(NULL,hist);
  GetBinButtonDims(&binwide,&binhigh);
  histwide = textwide + 72 + binwide;
  if (histwide < HISTWIDE) histwide = HISTWIDE;

  if(semilog) {
    char buffer[126];
    sprintf(buffer,"%lf*log_%lf : ", semilog_multiplier, semilog_base);
    HistoTitle = strcat(buffer,title);
  } else {
    HistoTitle = title;
  }

  packet->window = obj =
            mt_new_window(10,10,histwide,
                          HISTHIGH+HISTBUTTON+2*SCROLLWIDE,0,title);
  mt_set_window_bounds(obj,histwide,150,-1,-1);
  SetHistoBorderMax(hist,histwide,150-(HISTBUTTON+2*SCROLLWIDE));

  packet->darea = obj =
       mt_new_frame(0,HISTBUTTON,histwide,HISTHIGH,1,REGULAR,REDRAW_EVENT);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,1.);
  mt_set_outline(obj,0);
  mt_set_border(obj,0);
  mt_frame_colors(obj,BackGround,-1,-1,-1,-1);
  mt_set_callback(obj,histo_handler,(long) packet);
  mt_pop_frame();

  MakeHorizontalScroll(packet->sbar,0,HISTHIGH+HISTBUTTON,histwide,
                                    scrollupdate0,(long) (packet->darea));
  MakeHorizontalScroll(packet->sbar+1,0,HISTHIGH+HISTBUTTON+SCROLLWIDE,histwide,
                                     scrollupdate1,(long) (packet->darea));

  { int xmin, xmax;
    int high, tab;

    GetHistoDims(hist,&xmin,&xmax);
    high = xmax - xmin;
    if (high*.05 < 10)
      { tab = 10;
        if (high < 10) high = 10;
      }
    else
      tab = high*.05;
    mt_set_scrollbar(packet->sbar[0].bar,xmin,high,xmin,xmin);
    mt_set_scrollbar(packet->sbar[1].bar,high,tab,10,high);
  }

  obj = mt_new_frame(0,0,histwide,HISTBUTTON,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_hilite(obj,0);

    obj = mt_new_frame(BORDER+3,BORDER+3,textwide,HISTBUTTON-2*(BORDER+3),
                       1,REGULAR,REDRAW_EVENT);
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_frame_colors(obj,mt_white(),-1,-1,-1,-1);
    mt_set_hilite(obj,0);
    mt_set_border(obj,0);
    mt_set_callback(obj,histo_info,(long) packet);
    mt_pop_frame();

  { int cornerx, cornery;

    cornerx = BORDER + textwide + 3 + (histwide - (textwide + binwide + 53))/2;
    cornery = (HISTBUTTON - binhigh)/2;
    obj = mt_new_frame(cornerx,cornery,binwide,binhigh,1,REGULAR,REDRAW_EVENT);
    mt_set_xscale(obj,.5,0.);
    mt_set_yscale(obj,0.,0.);
    mt_set_hilite(obj,0);
    mt_set_border(obj,0);
    mt_set_outline(obj,0);
    mt_set_callback(obj,DrawBinButton,0);

      CreateBinButton(packet,binwide,binhigh);

    mt_pop_frame();
  }

  { int cornery;

    cornery = (HISTBUTTON - 24)/2;

    obj = mt_new_button(histwide-50,cornery,40,24,1,CLICK,"Quit");
    mt_set_xscale(obj,1.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,cbutton,(long) packet);
  }

  mt_pop_frame();

  mt_pop_frame();

  mt_set_vis(packet->window,1);

  mt_event_loop();

  exit(0);
}
