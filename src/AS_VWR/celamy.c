
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
/* $Id: celamy.c,v 1.4 2005-03-22 19:49:31 jason_miller Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "button.h"
#include "agrep.h"
#include "layout.h"

/* Geometry parameters */

#define FRAMEWIDE 600   /* Main window dimensions */
#define FRAMEHIGH 300
#define TEXTHIGH   48
#define SCROLLWIDE 16
#define BUTTONHIGH 36

#define BORDER      2   /* Border relief for all objects */

#define HISTWIDE  300   /* Histogram window dimensions */
#define HISTHIGH  250
#define HISTBUTTON 50

#define COVERHIGH 150   /* Coverage window dimensions */

/* Color parameters */

#undef SWITCH_BACKGROUND_COLOR

#ifndef SWITCH_BACKGROUND_COLOR
#define TEXTCOLOR  255,255,255   /* Text */
#define BACKCOLOR  0,0,0         /* Background */
#define MATCHCOLOR 255,255,255   /* Line Selection Color */
#else
#define TEXTCOLOR  0,0,0   /* Text */
#define BACKCOLOR  255,255,255         /* Background */
#define MATCHCOLOR 0,0,0   /* Line Selection Color */
#endif

#define SELCOLOR   200,200,100   /* General selection hilite color */
#define KSELCOLOR  125,125,70    /* Selection color for scrollbar buttons */
#define MSELCOLOR  150,150,0     /* Selection color for menu buttons */
#define KNOBCOLOR  50,50,225     /* Scrollbar button & slider color */
#define PICKCOLOR  80,200,80     /* Menu button on color */
#define ERSELCOLOR 255,150,150   /* Error text selection color */ 
#define TXSELCOLOR 0,255,255     /* Normal text selection color */
#define ERRORCOLOR 255,0,0       /* Error message text color */

/* Color values for above parameters */

long BackGround, TextColor, MatchColor, ErrorColor, ErSelColor, TxSelColor;

long SelectColor, KnobColor, KselColor, MselColor, PickColor;

/* Scroll bar object */

typedef struct stag {
    MT_OBJECT *bar;                       /* Central slider */
    MT_OBJECT *lft;                       /* Left/Up button */
    MT_OBJECT *rgt;                       /* Right/Down button */
    void (*updater)(long, struct stag *); /* Event handler routine */
    long      data;                       /* Value to pass to action routine */
} MY_SCROLL;

#define HSCROLL 0             /* Main window scrollbar indices */
#define HZOOM   1
#define VSCROLL 2
#define VZOOM   3

static MY_SCROLL sbars[4];    /* Scroll bars on main window */

static MT_OBJECT *Window;     /* Main window */
static MT_OBJECT *canvas;     /* Layout canvas within main window */
static MT_OBJECT *textpane;   /* Textpane within main window */
static MT_OBJECT *errorpane;  /* Errorpage within main window */

#define LINEMENU 0        /* Main window menu indices */
#define MARKMENU 1
#define LINKMENU 2
#define OVLPMENU 3

static MT_OBJECT  *Menu[4];          /* Main window menus */
static MT_OBJECT **MenuButtons[4];   /* Main window menu buttons */
static int        *MenuIds[4];       /* Main window menu item arrays */

static MT_OBJECT  *ToolMenu;          /* Main window tool menu */

static int MainWidth;                /* Main canvas width in pixels */
static int MainHeight;               /* Main canvas height in pixels */

static char *CelamyHome;	     /* Celamy home file */

/* Main Window canvas area handlers */

static int canvas_remap = 1;  /* 0 ==> no need to reset world<->device map */

void setup_map(void)
{ int val0, wid0, min, max;
  int val1, wid1;
  int xl, xh, yl, yh;

  val0 = mt_get_scrollbar(sbars[HSCROLL].bar,&wid0,&min,&max);
  val1 = mt_get_scrollbar(sbars[VSCROLL].bar,&wid1,&min,&max);
  mt_get_extent(canvas,&xl,&xh,&yl,&yh);
  MainWidth  = (xh-xl)+1;
  MainHeight = (yh-yl)+1;
  SetAssemblyMap(val0+1,val0+wid0,val1,(val1+wid1)-1,
                 (xh-xl)+1,(yh-yl)+1,BackGround,TextColor,MatchColor);
}

int canvas_handler(MT_OBJECT *o, long d)
{ int e, x, y, v, m;
  char *pick;
  void ResizeAllCovers(int);
  void viewreset(void);
  void canvasreset(void);
  void viewsetrange(int,int);
  static int zoom_x, zoom_mode;

  e = mt_get_event(o,&x,&y,&v,&m);
  if (e == REDRAW_EVENT)
    { if (canvas_remap)
        { v = MainWidth;
          setup_map();
          ResizeAllCovers(v != MainWidth);
          viewreset();
          canvasreset();
        }
      canvas_remap = 1;
      DrawAssembly(canvas);
    }
  else if (e == PRESS_EVENT)
    { if (zoom_mode = (y >= MainHeight - 18))
        zoom_x = x;
      else
        { if (m & SHIFT_MODIFIER)
            pick = PickAssembly(canvas,x,y,v,PICK_SEGMENT);
          else if (m & CNTRL_MODIFIER)
            pick = PickAssembly(canvas,x,y,v,PICK_PIECES);
          else
            pick = PickAssembly(canvas,x,y,v,PICK_NEAREST);
          if (pick != NULL)
            { mt_set_text(textpane,pick,0,0);
              mt_draw(textpane);
            }
        }
    }
  else if (e == RELEASE_EVENT)
    if (zoom_mode)
      { double xb, xe;
        int val, wid, min, max;
        int beg, end;

        if (x < zoom_x)
          { xb = x; xe = zoom_x; }
        else
          { xe = x; xb = zoom_x; }
        val = mt_get_scrollbar(sbars[HSCROLL].bar,&wid,&min,&max) + 1;
	if (xe - xb >= 3)
	  { beg = (int)( ( xb / MainWidth ) * wid ) + val;
            end = (int)( ( xe / MainWidth ) * wid ) + val;
            if (beg < 0) beg = 0;
            if (end < 0) end = 0;
          }
        else
          { beg = (int)( ( xb / MainWidth ) * wid ) + val;
            if (v <= 1)
              wid /= 4;
            end = beg + wid;
            beg = beg - wid;
            if (beg == end) end = end+1;
          }
        viewsetrange(beg,end);
      }
    else
      { PickRelease(canvas);
        mt_set_text(textpane,"",0,0);
        mt_draw(textpane);
      }
  return (0);
}

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

int    htabexp;
double alpha;

void scrollupdate2(long d, MY_SCROLL *s)  /* Slave update: geometric */
{ int val0, wid0, min0, max0;
  int val1, wid1, min1, max1;

  val0 = mt_get_scrollbar(s[-1].bar,&wid0,&min0,&max0);
  max0 += wid0;
  val1 = mt_get_scrollbar(s->bar,&wid1,&min1,&max1);
  min0   += htabexp;
  max0   -= htabexp;
  val1 = min1 + (int)((max1-min1) *  pow( (1.*val1-min1)/(max1-min1) , alpha ));
  htabexp = (int)(val1*.025);
  min0   -= htabexp;
  max0   += htabexp;
  wid0    = val1 + 2*htabexp;
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

/* Quit Button Handler */

int qbutton(MT_OBJECT *o, long d)
{ mt_shutdown();
  exit (0);
  return (0);
}

/* Query Button Handler */

int dbutton(MT_OBJECT *o, long d)
{
  ObjectPacket *objp;
  CoverPacket  *cover;
  char *query, *error;
  int   bsel, esel;
  void CreateHistoWin(HistoPacket *, char *, int);
  void CreateCoverWin(CoverPacket *, char *, int);
  void QuReferenceSetup(void);

  QuReferenceSetup();

  query = mt_get_text(textpane,&bsel,&esel);
  error = ProcessQuery(query,&objp,&cover);

  QueryTableFree();

  if (error != NULL)
    { mt_set_text(errorpane,error,0,0);
      mt_draw(errorpane);
      free(error);
    }
  else
    { if (cover != NULL)
        { CreateCoverWin(cover,query,1);
          CreateHistoWin(Cover2Hist(cover),query,1);
        }
      else if (objp != NULL)
        { CreateHistoWin(Object2Hist(objp),query,1);
          CreateCoverWin(Object2Cover(objp),query,1);
          FreeObjects(objp);
        }
    }
  return (0);
}

/* Object Menu Routines */

static int MenuRadix;

void process_state(int menu)
{ int i, count;

  count = StartStyleList(menu);
  for (i = 0; i < count; i++)
    SetStyleVis(MenuIds[menu][i],menu,1-mt_get_button(MenuButtons[menu][i]));
  canvas_remap = 0;
  mt_draw(canvas);
}

int menuitem(MT_OBJECT *frame, long d)
{ char label[20], *name;
  int  id, thick, dash, menu;
  long color;
  MT_OBJECT *button;

  menu = d/MenuRadix;
  if (mt_get_vis(Menu[menu]))
    { mt_set_color(mt_black());
      id = MenuIds[menu][d%MenuRadix];
      GetStyle(id,&color,&thick,&dash,&name);
      sprintf(label,"%s",name);
      mt_draw_text(frame,64,15,label);
      mt_set_color(BackGround);
      mt_draw_rect(frame,20,3,40,14);
      mt_set_color(color);
      switch (menu)
      { case LINEMENU:
         mt_draw_line(frame,23,10,57,10,thick,dash);
         break;
        case MARKMENU:
         mt_draw_line(frame,30,7,30,13,thick,dash);
         mt_draw_line(frame,40,7,40,13,thick,dash);
         mt_draw_line(frame,50,7,50,13,thick,dash);
         break;
        case LINKMENU:
         mt_draw_line(frame,23,10,57,10,thick,dash);
         mt_draw_line(frame,23,10-thick/2,23,17,thick,dash);
         mt_draw_line(frame,57,10-thick/2,57,17,thick,dash);
         break;
        case OVLPMENU:
         mt_draw_rect(frame,30,7,20,6);
         mt_set_color(TextColor);
         mt_draw_line(frame,20,7,50,7,2,0);
         mt_draw_line(frame,30,13,60,13,2,0);
         break;
      }
    }
  else
    { button = MenuButtons[menu][d%MenuRadix];
      mt_set_button(button,1-mt_get_button(button));
      process_state(menu);
    }
  return (0);
}

int menubutton(MT_OBJECT *button, long d)
{ if (!mt_get_vis(Menu[d]))
    { mt_set_button(button,1-mt_get_button(button));
      process_state(d);
    }
  return (0);
}

int menuframe(MT_OBJECT *frame, long d)
{ if (!mt_get_vis(Menu[d]))
    process_state(d);
  return (0);
}

int menunone(MT_OBJECT *frame, long d)
{ if (mt_get_vis(Menu[d]))
    { mt_set_color(mt_black());
      mt_draw_text(frame,30,BORDER+16,"NONE");
    }
  return (0);
}

void makemenu(int d, int x, int y)
{ MT_OBJECT *obj;
  int count, i, id;
  int menu_width;

  MenuRadix = NumberOfStyles();

  count = StartStyleList(d);
  menu_width = 0;
  if (count > 0)
    MenuIds[d] = (int *) malloc(sizeof(int)*count);
  for (i = 0; i < count; i++)
    { int   thick, dash, wide, high, base;
      long  color;
      char *name;

      id = NextUsedStyle();
      MenuIds[d][i] = id;
      GetStyle(id,&color,&thick,&dash,&name);
      mt_string_size(name,&wide,&high,&base);
      if (wide > menu_width) menu_width = wide;
    }

  if (count == 0)
    { menu_width = 80;
      Menu[d] = mt_new_frame(x,y,menu_width,2*(BORDER+1)+20,0,
                             REGULAR,REDRAW_EVENT);
      mt_set_callback(Menu[d],menunone,d);
    }
  else
    { menu_width += 72; 
      Menu[d] = mt_new_frame(x,y,menu_width,2*(BORDER+1)+20*count,0,
                             REGULAR,REDRAW_EVENT);
      mt_set_callback(Menu[d],menuframe,d);
      MenuButtons[d] = (MT_OBJECT **) malloc(sizeof(MT_OBJECT *)*count);
    }
  mt_set_hilite(Menu[d],0);

  mt_current_outline(0);
  for (i = 0; i < count; i++)
    { obj = mt_new_frame(BORDER+1,BORDER+1+20*i,menu_width-2*(BORDER+1),20,1,
                         REGULAR,REDRAW_EVENT);
      mt_frame_colors(obj,-1,mt_medium_grey(),mt_medium_grey(),
                           MselColor,MselColor);
      mt_set_callback(obj,menuitem,MenuRadix*d+i);
        obj = mt_new_button(4,4,12,12,1,TOGGLE,NULL);
        MenuButtons[d][i] = obj;
        mt_button_colors(obj,PickColor,mt_medium_grey(),
                             -1,-1,MselColor,-1);
        mt_set_outline(obj,1);
        mt_set_callback(obj,menubutton,d);
      mt_pop_frame();
    }
  mt_current_outline(1);

  mt_pop_frame();
}

/* Tool menu routines */

void CreateQueryWin(long);
void CreateViewWin(long);
void CreateCanvasWin(long);

typedef void procpt(long);

static char   *ToolLabels[] = { "Query ...", "View ...",
                                "Canvas ...", "" };
static int     ToolIsOn[]   = { 0, 0, 0, 0 };
static procpt *ToolCalls[]  = { CreateQueryWin, CreateViewWin,
                                CreateCanvasWin, };

int toolitem(MT_OBJECT *frame, long d)
{ char label[20];
  int e, x, y, v, m;

  e = mt_get_event(ToolMenu,&x,&y,&v,&m);
  if (mt_get_vis(ToolMenu) && e == REDRAW_EVENT)
    { mt_set_color(mt_black());
      sprintf(label,"%s",ToolLabels[d]);
      mt_draw_text(frame,2,15,label);
    }
  else
    { if (e == PRESS_EVENT) mt_set_vis(ToolMenu,0);
      if (!ToolIsOn[d])
        { ToolIsOn[d] = 1;
          ToolCalls[d](d);
        }
    }
  return (0);
}

void maketool(int x, int y)
{ MT_OBJECT *obj;
  int count, i;
  int menu_width;

  menu_width = 0;
  for (count = 0; ToolLabels[count][0] != '\0'; count++)
    { int   wide, high, base;

      mt_string_size(ToolLabels[count],&wide,&high,&base);
      if (wide > menu_width) menu_width = wide;
    }

  menu_width += 2*(BORDER+3); 
  ToolMenu = mt_new_frame(x,y,menu_width,2*(BORDER+1)+20*count,0,
                         REGULAR,REDRAW_EVENT);
  mt_set_hilite(ToolMenu,0);

  mt_current_outline(0);
  for (i = 0; i < count; i++)
    { obj = mt_new_frame(BORDER+1,BORDER+1+20*i,menu_width-2*(BORDER+1),20,1,
                         REGULAR,REDRAW_EVENT|PRESS_EVENT);
      mt_frame_colors(obj,-1,mt_medium_grey(),mt_medium_grey(),
                           MselColor,MselColor);
      mt_set_callback(obj,toolitem,i);
      mt_pop_frame();
    }
  mt_current_outline(1);

  mt_pop_frame();
}

/* Create main window and all object there in and start up event loop */

int main(int argc, char *argv[])
{ MT_OBJECT *obj;
  FILE *file;
  static char Title[100];
  void QueryInit(void);
  int ch;
  int batch=0;
  int noptargc=argc;
  optarg = NULL;
  while ( ((ch = getopt(argc, argv, "h?b")) != EOF)) {
        switch(ch) {
        case 'b':
          batch = 1;
          break;
        case 'h':
        case '?':
          fprintf(stderr,"  Usage:\n\n");
          fprintf(stderr," %s [-b] file.cam\n",argv[0]);
          fprintf(stderr,"\t\t-b\t(opt)\tbatch process queries from stdin\n");
          exit(1);
        default: 
          fprintf(stderr,"  Incorrect usage. Please try %s -h\n",argv[0]);
          exit(1);
      }
  }

  mt_startup();

  mt_current_border(BORDER);
  SelectColor = mt_get_color(SELCOLOR);
  KnobColor   = mt_get_color(KNOBCOLOR);
  KselColor   = mt_get_color(KSELCOLOR);
  MselColor   = mt_get_color(MSELCOLOR);
  PickColor   = mt_get_color(PICKCOLOR);
  TextColor   = mt_get_color(TEXTCOLOR);
  ErrorColor  = mt_get_color(ERRORCOLOR);
  BackGround  = mt_get_color(BACKCOLOR);
  MatchColor  = mt_get_color(MATCHCOLOR);
  ErSelColor  = mt_get_color(ERSELCOLOR);
  TxSelColor  = mt_get_color(TXSELCOLOR);

  if (optind > 1) { noptargc = argc - optind + 1;}
  if (noptargc < 2)
    { fprintf(stderr,"Usage: Celamy file ...\n");
      exit (1);
    }

  { char *home;

    CelamyHome = getenv("HOME");
    if (CelamyHome != NULL)
      { home = (char *) malloc(strlen(CelamyHome)+50);
        strcpy(home,CelamyHome);
        strcpy(home+strlen(CelamyHome),"/.celquery");
        CelamyHome = home;
      }
  }

  { int i;

    for (i = 1; i < noptargc; i++)
      { file = fopen(argv[i+optind-1],"r");
        if (file == NULL)
          { fprintf(stderr,"Can't open file %s\n",argv[i+optind-1]);
            exit (1);
          }
        if (noptargc > 2)
          SetCurrentFile(argv[i+optind-1]);
        ReadAssembly(file);
      }
  }

  if (argc > 2)
    sprintf(Title,"Celamy: %.30s ...",argv[1+optind]);
  else
    sprintf(Title,"Celamy: %.30s",argv[1]);

  BuildAssembly();
  LayoutAssembly();

  Window = mt_new_window(10,10,FRAMEWIDE+2*SCROLLWIDE,
                         FRAMEHIGH+2*SCROLLWIDE+BUTTONHIGH+TEXTHIGH,0,Title);
  mt_set_window_bounds(Window,360,200,-1,-1);

  canvas = mt_new_frame(0,BUTTONHIGH,FRAMEWIDE,FRAMEHIGH,1,
                        REGULAR,REDRAW_EVENT|PRESS_EVENT|RELEASE_EVENT);
  mt_set_xscale(canvas,0.,1.);
  mt_set_yscale(canvas,0.,1.);
  mt_set_outline(canvas,0);
  mt_set_border(canvas,0);
  mt_frame_colors(canvas,BackGround,-1,-1,-1,-1);
  mt_set_callback(canvas,canvas_handler,-1);
  mt_pop_frame();

  obj = mt_new_frame(0,0,FRAMEWIDE+2*SCROLLWIDE,BUTTONHIGH,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_hilite(obj,0);

  { int cornery;
    int cornerl;
    int midpnt;

    cornery = (BUTTONHIGH - 24)/2;
    cornerl = cornery + 24;
    midpnt  = (FRAMEWIDE + 2*SCROLLWIDE - 210) / 2;

    mt_current_xscale(0.,0.);
    mt_current_yscale(0.,0.);

    maketool(10,cornerl);
    obj = mt_new_button(10,cornery,40,24,1,PRESS,"Tools");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_menu_pair(obj,ToolMenu);

    mt_current_xscale(.5,0.);

    makemenu(LINEMENU,midpnt,cornerl);
    obj = mt_new_button(midpnt,cornery,40,24,1,PRESS,"Lines");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_menu_pair(obj,Menu[LINEMENU]);

    makemenu(MARKMENU,midpnt+50,cornerl);
    obj = mt_new_button(midpnt+50,cornery,40,24,1,PRESS,"Marks");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_menu_pair(obj,Menu[MARKMENU]);

    makemenu(LINKMENU,midpnt+100,cornerl);
    obj = mt_new_button(midpnt+100,cornery,40,24,1,PRESS,"Links");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_menu_pair(obj,Menu[LINKMENU]);

    makemenu(OVLPMENU,midpnt+150,cornerl);
    obj = mt_new_button(midpnt+150,cornery,60,24,1,PRESS,"Overlaps");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_menu_pair(obj,Menu[OVLPMENU]);

    mt_current_xscale(0.,0.);

    obj = mt_new_button(FRAMEWIDE+2*SCROLLWIDE-50,
                         cornery,40,24,1,CLICK,"Quit");
    mt_set_xscale(obj,1.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,qbutton,0);
  }
 
  mt_pop_frame();

  mt_current_hilite(1);

  MakeHorizontalScroll(sbars+HSCROLL,0,FRAMEHIGH+BUTTONHIGH,
                                     FRAMEWIDE,scrollupdate0,(long) canvas);
  MakeHorizontalScroll(sbars+HZOOM,0,FRAMEHIGH+BUTTONHIGH+SCROLLWIDE,
                                   FRAMEWIDE,scrollupdate2,(long) canvas);

  MakeVerticalScroll(sbars+VSCROLL,FRAMEWIDE,BUTTONHIGH,
                                   FRAMEHIGH,scrollupdate0,(long) canvas);
  MakeVerticalScroll(sbars+VZOOM,FRAMEWIDE+SCROLLWIDE,BUTTONHIGH,
                                 FRAMEHIGH,scrollupdate1,(long) canvas);

  obj = mt_new_frame(FRAMEWIDE,FRAMEHIGH+BUTTONHIGH,
                     2*SCROLLWIDE,2*SCROLLWIDE,1,REGULAR,0);
  mt_set_xscale(obj,1.,0.);
  mt_set_yscale(obj,1.,0.);
  mt_set_hilite(obj,0);
  mt_pop_frame();

  { int min, max, rows, smallest;
    int wide, high, tab;

    GetAssemblyDims(&min,&max,&rows,&smallest);
    min -= 1;
    wide = (max-min);
    tab  = (int)(wide*.05);
    mt_set_scrollbar(sbars[HZOOM].bar,wide,tab,2*smallest,wide);
    alpha = log((1.*smallest)/wide) / log(.025);
    if (alpha < 1.) alpha = 1.;
    htabexp = (int)(.025*wide);
    min    -= htabexp;
    max    += htabexp;
    wide   += 2*htabexp;
    mt_set_scrollbar(sbars[HSCROLL].bar,min,wide,min,min);
    if (rows < FRAMEHIGH/10)
      high = FRAMEHIGH / 10;
    else
      high = rows;
    if (rows*.05 < 2)
      tab = 2;
    else
      tab = (int)(rows*.05);
    mt_set_scrollbar(sbars[VZOOM].bar,high,tab,2,high);
    mt_set_scrollbar(sbars[VSCROLL].bar,0,high,0,0);
  }

  obj = mt_new_frame(0,BUTTONHIGH+FRAMEHIGH+2*SCROLLWIDE,
                     FRAMEWIDE+2*SCROLLWIDE,TEXTHIGH,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,1.,0.);
  mt_set_hilite(obj,0);

  { int cornery;

    cornery = (TEXTHIGH - 24)/2;

    obj = mt_new_button(10,cornery,40,24,1,CLICK,"Query");
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_xscale(obj,0.,0.);
    mt_set_callback(obj,dbutton,0);

    obj = mt_new_frame(60,6,FRAMEWIDE+2*SCROLLWIDE-66,TEXTHIGH-12,1,REGULAR,0);
    mt_set_xscale(obj,0.,1.);
    mt_set_hilite(obj,0);
    mt_frame_colors(obj,-1,mt_light_grey(),mt_dark_grey(),-1,-1);

      errorpane = mt_new_textbox(3,3+(TEXTHIGH-18)/2,
                              FRAMEWIDE+2*SCROLLWIDE-72,(TEXTHIGH-18)/2,1,NULL);
      mt_text_colors(errorpane,-1,-1,-1,ErrorColor,TxSelColor);
      mt_set_xscale(errorpane,0.,1.);
      mt_set_border(errorpane,0);
      mt_set_outline(errorpane,0);

      textpane = mt_new_textbox(3,3,FRAMEWIDE+2*SCROLLWIDE-72,
                                    TEXTHIGH-18,1,"\n\r");
      mt_text_colors(textpane,-1,-1,-1,-1,TxSelColor);
      mt_set_xscale(textpane,0.,1.);
      mt_set_border(textpane,0);
      mt_set_outline(textpane,0);
      mt_set_callback(textpane,dbutton,0);

    mt_pop_frame();
  }
  mt_pop_frame();

  mt_pop_frame();

  setup_map();

  mt_set_focus(textpane);

  mt_set_vis(Window,1);

  QueryInit();

if (batch) {
  char query[1001];
  ObjectPacket *objp;
  CoverPacket  *cover;
  char *error;
  void CreateHistoWin(HistoPacket *, char *, int);
  void CreateCoverWin(CoverPacket *, char *, int);
  void QuReferenceSetup(void);

  QuReferenceSetup();

  // read queries from standard in and process, then exit
  while ( fgets(query,1000,stdin) != NULL ) {
    error = ProcessQuery(query,&objp,&cover);

    QueryTableFree();
 
    query[strlen(query)-1] = ':';
    fprintf(stdout,"%s\t",query);
    query[strlen(query)-1] = '\n';

    if (error != NULL)
    { mt_set_text(errorpane,error,0,0);
      mt_draw(errorpane);
      free(error);
    }
    else
    { if (cover != NULL)
        { CreateCoverWin(cover,query,1);
          CreateHistoWin(Cover2Hist(cover),query,1);
          FreeCoverage(cover);
        }
      else if (objp != NULL)
        { CreateHistoWin(Object2Hist(objp),query,1);
          CreateCoverWin(Object2Cover(objp),query,1);
          FreeObjects(objp);
        }
    }
    fprintf(stdout,"\n");
  }
  exit(0);
} 

  mt_event_loop();

  exit(0);
}

/* View & Window Setting routines */

typedef struct {
    MT_OBJECT *widgets[3];
    char      *labels[3];
    int        last_width;
    int        toolno;
    int       (*actionroutine)(MT_OBJECT *, long);
    void      (*fetchroutine)(int *, int *);
} pairbox;

int closepair(MT_OBJECT *o, long d)
{ pairbox *box;
  int xl, xh, yl, yh;

  box = (pairbox *) d;
  mt_get_extent(box->widgets[0],&xl,&xh,&yl,&yh);
  box->last_width = (xh-xl)+1;
  mt_free(box->widgets[0]); 
  ToolIsOn[box->toolno] = 0;
  return (0);
}

int getinteger(MT_OBJECT *o)
{ char *text, *end;
  int bsel, esel, value;

  text  = mt_get_text(o,&bsel,&esel);
  value = strtol(text,&end,10);
  if (*end != '\0' || end == text || value < 0)
    { mt_text_colors(o,-1,-1,-1,-1,ErSelColor);
      if (end == text)
        mt_set_text(o,"?",0,1);
      else if (value < 0)
        mt_set_text(o,NULL,0,strlen(text));
      else
        mt_set_text(o,NULL,end-text,strlen(text));
      mt_draw(o);
      mt_set_focus(o);
      mt_text_colors(o,-1,-1,-1,-1,TxSelColor);
      return (-1);
    }
  else
    return (value);
}

void pairreset(pairbox *box)
{ int no1, no2;
  static char label[50];

  if (box->widgets[0] == NULL || !mt_get_vis(box->widgets[0])) return;

  box->fetchroutine(&no1,&no2);

  sprintf(label,"%d",no1);
  mt_set_text(box->widgets[1],label,0,strlen(label));
  mt_set_focus(box->widgets[1]);
  mt_draw(box->widgets[1]);

  sprintf(label,"%d",no2);
  mt_set_text(box->widgets[2],label,0,0);
  mt_draw(box->widgets[2]);
}

int pairfocus(MT_OBJECT *o, long d)
{ pairbox *box;
  char *text;
  int x, y, v, m, e;

  box = (pairbox *) d;
  e = mt_get_event(o,&x,&y,&v,&m);
  if (v != '\t')
    { box->actionroutine(o,d);
      return (0);
    }
  if (box->widgets[1] == o)
    o = box->widgets[2];
  else
    o = box->widgets[1];
  mt_set_focus(o);
  text = mt_get_text(o,&x,&y);
  mt_set_text(o,NULL,0,strlen(text));
  mt_draw(o);
  return (0);
}

int pairlabels(MT_OBJECT *o, long d)
{ pairbox *box;
  int xl, xh, yl, yh;

  box = (pairbox *) d;
  mt_set_color(mt_black());
  mt_draw_text(o,10,18,box->labels[1]);
  mt_get_extent(o,&xl,&xh,&yl,&yh);
  mt_draw_text(o,(xh-xl)/2+5,18,box->labels[2]);
  return (0);
}

void CreatePairWin(pairbox *box)
{ MT_OBJECT *obj;
  int no1, no2, lwide;
  static char label[50];

  lwide = box->last_width;
  box->widgets[0] = obj = 
       mt_new_window(10,10,lwide,60+BUTTONHIGH,0,box->labels[0]);
  mt_set_window_bounds(obj,180,60+BUTTONHIGH,-1,60+BUTTONHIGH);

  obj = mt_new_frame(0,BUTTONHIGH,lwide,60,1,REGULAR,REDRAW_EVENT);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,1.);
  mt_set_outline(obj,0);
  mt_set_border(obj,0);
  mt_set_callback(obj,pairlabels,(long) box);

    box->widgets[1] = obj = mt_new_textbox(10,22,(lwide-30)/2,22,1,"\t\n\r");
    mt_set_xscale(obj,0.,.5);
    mt_set_callback(obj,pairfocus,(long) box);
    mt_text_colors(obj,-1,-1,-1,-1,TxSelColor);

    box->widgets[2] = obj = mt_new_textbox(20+(lwide-30)/2,22,
                                           (lwide-30)/2,22,1,"\t\n\r");
    mt_set_xscale(obj,.5,.5);
    mt_set_callback(obj,pairfocus,(long) box);
    mt_text_colors(obj,-1,-1,-1,-1,TxSelColor);

    box->fetchroutine(&no1,&no2);

    sprintf(label,"%d",no1);
    mt_set_text(box->widgets[1],label,0,strlen(label));
    mt_set_focus(box->widgets[1]);
    sprintf(label,"%d",no2);
    mt_set_text(box->widgets[2],label,0,0);

  mt_pop_frame();

  obj = mt_new_frame(0,0,lwide,BUTTONHIGH,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_hilite(obj,0);

  { int cornery;

    cornery = (BUTTONHIGH - 24)/2;

    obj = mt_new_button(10,cornery,40,24,1,CLICK,"Set");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,box->actionroutine,(long) box);

    obj = mt_new_button(lwide-50,cornery,40,24,1,CLICK,"Close");
    mt_set_xscale(obj,1.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,closepair,(long) box);
  }
 
  mt_pop_frame();

  mt_pop_frame();

  mt_set_vis(box->widgets[0],1);
}

void viewfetch(int *lft, int *rgt)
{ int val, wid, min, max;

  val  = mt_get_scrollbar(sbars[HSCROLL].bar,&wid,&min,&max);
  *lft = val + 1 + htabexp;
  *rgt = val + wid - htabexp;
}

void viewsetrange(int beg, int end)
{ int val, wid, min, max;
  int zval, zwid, zmin, zmax;

  mt_get_scrollbar(sbars[HZOOM].bar,&zwid,&zmin,&zmax);

  val = mt_get_scrollbar(sbars[HSCROLL].bar,&wid,&min,&max);
  max = max + wid - htabexp;
  min = min + 1 + htabexp;

  if (end - beg < zmin)
    end = beg + zmin;
  if (beg < min)
    beg = min;
  if (end > max)
    end = max;
  if (end - beg < zmin)
    beg = end - zmin;

  val = beg-1;
  wid = end-val;
  min = min-1;
  zval = wid;
  htabexp = (int)(.025*zval);
  min -= htabexp;
  max += htabexp;
  wid = zval + 2*htabexp;
  val -= htabexp;
  if (wid > max-min) wid = max-min;
  if (val+wid > max) val = max-wid;
  if (val < min) val = min;
  mt_set_scrollbar(sbars[HSCROLL].bar,val,wid,min,max-wid);
  zval = zmin + (int)((zmax-zmin) *  pow( (1.*zval-zmin)/(zmax-zmin) , 1./alpha ));
  mt_set_scrollbar(sbars[HZOOM].bar,zval,zwid,zmin,zmax);
  mt_draw(sbars[HSCROLL].bar);
  mt_draw(sbars[HZOOM].bar);
  mt_draw(canvas);
}

int viewset(MT_OBJECT *o, long d)
{ pairbox *box;
  int beg, end;
  int wid, min, max;

  box = (pairbox *) d;

  beg = getinteger(box->widgets[1]);
  if (beg < 0) return (0);
  end = getinteger(box->widgets[2]);
  if (end < 0) return (0);

  if (beg >= end)
    { MT_OBJECT *o;

      mt_get_scrollbar(sbars[HSCROLL].bar,&wid,&min,&max);
      if (end <= max && beg > max)
        o = box->widgets[1];
      else
        o = box->widgets[2];
      mt_text_colors(o,-1,-1,-1,-1,ErSelColor);
      mt_set_text(o,NULL,0,strlen(mt_get_text(o,&beg,&end)));
      mt_draw(o);
      mt_set_focus(o);
      mt_text_colors(o,-1,-1,-1,-1,TxSelColor);
      return (0);
    }
  
  viewsetrange(beg,end);
  return (0);
}

static pairbox ViewBox;

void viewreset(void)
{ pairreset(&ViewBox); }

void CreateViewWin(long d)
{ static char *title = "Set View";
  static char *lab1  = "From:";
  static char *lab2  = "To:";
  static int Firstime = 1;

  if (Firstime)
    { Firstime = 0;
      ViewBox.labels[0] = title;
      ViewBox.labels[1] = lab1;
      ViewBox.labels[2] = lab2;
      ViewBox.last_width = 180;
      ViewBox.toolno    = d;
      ViewBox.fetchroutine  = viewfetch;
      ViewBox.actionroutine = viewset; 
    }
  CreatePairWin(&ViewBox);
}

void canvasfetch(int *lft, int *rgt)
{ int xl, xh, yl, yh;

  mt_get_extent(canvas,&xl,&xh,&yl,&yh);
  *lft = xh - xl + 1;
  *rgt = yh - yl + 1;
}

int canvasset(MT_OBJECT *o, long d)
{ pairbox *box;
  int beg, end;
  static char label[50];

  box = (pairbox *) d;

  beg = getinteger(box->widgets[1]);
  if (beg < 0) return (0);
  end = getinteger(box->widgets[2]);
  if (end < 0) return (0);

  beg += 2*SCROLLWIDE;
  end += 2*SCROLLWIDE + BUTTONHIGH + TEXTHIGH;
  if (beg < 360)
    { beg = 360;
      sprintf(label,"%d",beg);
      mt_set_text(box->widgets[1],label,0,0);
      mt_draw(box->widgets[1]);
    }
  if (end < 200)
    { end = 200;
      sprintf(label,"%d",end);
      mt_set_text(box->widgets[2],label,0,0);
      mt_draw(box->widgets[2]);
    }

  mt_set_window_size(Window,beg,end);
  return (0);
}

static pairbox CanvasBox;

void canvasreset(void)
{ pairreset(&CanvasBox); }

void CreateCanvasWin(long d)
{ static char *title = "Set Canvas";
  static char *lab1  = "Width:";
  static char *lab2  = "Height:";
  static int Firstime = 1;

  if (Firstime)
    { Firstime = 0;
      CanvasBox.labels[0] = title;
      CanvasBox.labels[1] = lab1;
      CanvasBox.labels[2] = lab2;
      CanvasBox.last_width = 180;
      CanvasBox.toolno    = d;
      CanvasBox.fetchroutine  = canvasfetch;
      CanvasBox.actionroutine = canvasset; 
    }
  CreatePairWin(&CanvasBox);
}

/* Color Selection Dialog */

#define COLOR_LEVELS 59
#define COLOR_WIDTH   4

static MT_OBJECT *colorwin;
static MT_OBJECT *colhue[60];
static MT_OBJECT *colsat[60];
static MT_OBJECT *colval[60];
static MT_OBJECT *curcolbut;

static int curhue, cursat, curval;

static int CBarLabelX;
static int CButLabelX;

typedef void (*colorcall)(int,int,int);

static colorcall ColorCallBack;

void HSV_2_RGB(int hue, int val, int sat,
               int *red, int *green, int *blue)
{ int p, q, t, f;

  val = (255*val) / 59;
  if (sat == 0)
    *red = *green = *blue = val;
  else
    { f = hue - (hue/10)*10;
      p = (val*(59 - sat))/59;
      q = (val*(590 - sat*f))/590;
      t = (val*(590 - sat*(10-f)))/590;
      switch (hue/10)
      { case 0:
          *red = val; *green = t; *blue = p; break;
        case 1:
          *red = q; *green = val; *blue = p; break;
        case 2:
          *red = p; *green = val; *blue = t; break;
        case 3:
          *red = p; *green = q; *blue = val; break;
        case 4:
          *red = t; *green = p; *blue = val; break;
        case 5:
          *red = val; *green = p; *blue = q; break;
      }
    }   
}

void SetWheelColors(void)
{ int i, red, green, blue;
  long base;

  for (i = 0; i <= COLOR_LEVELS; i++)
    { HSV_2_RGB(i,curval,cursat,&red,&green,&blue);
      base = mt_get_color(red,green,blue);
      mt_button_colors(colhue[i],base,base,-1,-1,-1,-1);
      mt_draw(colhue[i]);
    }
  for (i = 0; i <= COLOR_LEVELS; i++)
    { HSV_2_RGB(curhue,i,cursat,&red,&green,&blue);
      base = mt_get_color(red,green,blue);
      mt_button_colors(colval[i],base,base,-1,-1,-1,-1);
      mt_draw(colval[i]);
    }
  for (i = 0; i <= COLOR_LEVELS; i++)
    { HSV_2_RGB(curhue,curval,i,&red,&green,&blue);
      base = mt_get_color(red,green,blue);
      mt_button_colors(colsat[i],base,base,-1,-1,-1,-1);
      mt_draw(colsat[i]);
    }
  HSV_2_RGB(curhue,curval,cursat,&red,&green,&blue);
  base = mt_get_color(red,green,blue);
  mt_button_colors(curcolbut,base,base,mt_black(),mt_black(),-1,-1);
  mt_draw(curcolbut);
}

int HuePress(MT_OBJECT *button, long d)
{ mt_set_outline(colhue[curhue],0);
  curhue = d;
  mt_set_outline(colhue[curhue],1);
  SetWheelColors();
  return (0);
}

int SatPress(MT_OBJECT *button, long d)
{ mt_set_outline(colsat[cursat],0);
  cursat = d;
  mt_set_outline(colsat[cursat],1);
  SetWheelColors();
  return (0);
}

int ValPress(MT_OBJECT *button, long d)
{ mt_set_outline(colval[curval],0);
  curval = d;
  mt_set_outline(colval[curval],1);
  SetWheelColors();
  return (0);
}

int ColorLabels(MT_OBJECT *obj, long d)
{ mt_set_color(mt_black());
  mt_draw_text(colorwin,CBarLabelX,26,"Hue");
  mt_draw_text(colorwin,CBarLabelX,58,"Saturation");
  mt_draw_text(colorwin,CBarLabelX,90,"Value");
  mt_draw_text(colorwin,CButLabelX,122,"Selected");
  mt_draw_text(colorwin,CButLabelX,138,"Color");
  return (0);
}

int CloseColorWin(MT_OBJECT *obj, long d)
{ mt_free(colorwin);
  if (d)
    ColorCallBack(curhue,cursat,curval);
  return (0);
}

void CreateColorWin(int hue, int sat, int val, colorcall response)
{ MT_OBJECT *obj;
  int        i, w, h;
  int        txtwide, winwidth;

  ColorCallBack = response;
  curhue = hue;
  cursat = sat;
  curval = val;

  winwidth = (COLOR_LEVELS+1)*COLOR_WIDTH + 25;
  mt_string_size("Saturation",&txtwide,&w,&h);
  winwidth += txtwide;

  colorwin = obj = mt_new_dialog(10,10,winwidth,198,0,"Color Selector");
  mt_set_window_bounds(obj,winwidth,198,winwidth,198);

  mt_current_outline(0);

  obj = mt_new_frame(0,0,winwidth,300,1,REGULAR,REDRAW_EVENT);
  CBarLabelX = winwidth-(txtwide+9);
  mt_set_callback(obj,ColorLabels,winwidth-(txtwide+9));
  mt_set_border(obj,0);

  obj = mt_new_frame(7,8,(COLOR_LEVELS+1)*COLOR_WIDTH + 6,24,1,REGULAR,0);
  mt_frame_colors(obj,-1,mt_medium_grey(),mt_medium_grey(),
                         MselColor,MselColor);
  mt_set_border(obj,2);
  obj = mt_new_frame(2,2,(COLOR_LEVELS+1)*COLOR_WIDTH + 2,20,1,REGULAR,0);
  mt_set_border(obj,0);
  mt_set_outline(obj,1);
    for (i = 0; i <= COLOR_LEVELS; i++)
      { colhue[i] = obj =
           mt_new_button(1+i*COLOR_WIDTH,1,COLOR_WIDTH,18,1,CLICK,NULL);
        mt_set_border(obj,0);
        mt_set_hilite(obj,0);
        mt_set_callback(obj,HuePress,i);
      }
  mt_pop_frame();
  mt_pop_frame();

  obj = mt_new_frame(7,40,(COLOR_LEVELS+1)*COLOR_WIDTH + 6,24,1,REGULAR,0);
  mt_frame_colors(obj,-1,mt_medium_grey(),mt_medium_grey(),
                         MselColor,MselColor);
  mt_set_border(obj,2);
  obj = mt_new_frame(2,2,(COLOR_LEVELS+1)*COLOR_WIDTH + 2,20,1,REGULAR,0);
  mt_set_border(obj,0);
  mt_set_outline(obj,1);
    for (i = 0; i <= COLOR_LEVELS; i++)
      { colsat[i] = obj =
           mt_new_button(1+i*COLOR_WIDTH,1,COLOR_WIDTH,18,1,CLICK,NULL);
        mt_set_border(obj,0);
        mt_set_hilite(obj,0);
        mt_set_callback(obj,SatPress,i);
      }
  mt_pop_frame();
  mt_pop_frame();

  obj = mt_new_frame(7,72,(COLOR_LEVELS+1)*COLOR_WIDTH + 6,24,1,REGULAR,0);
  mt_frame_colors(obj,-1,mt_medium_grey(),mt_medium_grey(),
                         MselColor,MselColor);
  mt_set_border(obj,2);
  obj = mt_new_frame(2,2,(COLOR_LEVELS+1)*COLOR_WIDTH + 2,20,1,REGULAR,0);
  mt_set_border(obj,0);
  mt_set_outline(obj,1);
    for (i = 0; i <= COLOR_LEVELS; i++)
      { colval[i] = obj = 
           mt_new_button(1+i*COLOR_WIDTH,1,COLOR_WIDTH,18,1,CLICK,NULL);
        mt_set_border(obj,0);
        mt_set_hilite(obj,0);
        mt_set_callback(obj,ValPress,i);
      }
  mt_pop_frame();
  mt_pop_frame();

  w = 20 + (COLOR_LEVELS+1)*COLOR_WIDTH;
  curcolbut = mt_new_button((w-150)/2,106,150,40,1,PRESS,NULL);
  CButLabelX = (w-150)/2 + 153;
  mt_set_border(curcolbut,2);
  mt_set_hilite(curcolbut,0);

  mt_current_outline(1);

  obj = mt_new_button((w-140)/2,160,60,24,1,CLICK,"Select");
  mt_set_callback(obj,CloseColorWin,1);
  obj = mt_new_button((w-140)/2+80,160,60,24,1,CLICK,"Cancel");
  mt_set_callback(obj,CloseColorWin,0);

  mt_set_outline(colhue[curhue],1);
  mt_set_outline(colsat[cursat],1);
  mt_set_outline(colval[curval],1);
  SetWheelColors();

  mt_set_vis(colorwin,1);
}

/* Save Dialog */

static MT_OBJECT *savewin;

int CloseSaveWin(MT_OBJECT *obj, long d)
{ void QueryWrite(void);
  if (d) QueryWrite();
  mt_free(savewin);
  return (0);
}

int DrawSaveText(MT_OBJECT *obj, long d)
{ mt_set_color(mt_black());
  mt_draw_text(obj,16,20,"Overwrite existing .celquery file?");
  return (0);
}

void SaveDialog(void)
{ MT_OBJECT *obj;

  savewin = obj = mt_new_dialog(10,10,240,64,0,"Save Dialog");
  mt_set_window_bounds(obj,240,64,240,64);

  obj = mt_new_frame(0,0,240,20,1,VIRTUAL,REDRAW_EVENT);
  mt_set_callback(obj,DrawSaveText,0);
  mt_pop_frame();
  
  obj = mt_new_button(40,30,60,24,1,CLICK,"Write");
  mt_set_callback(obj,CloseSaveWin,1);
  obj = mt_new_button(140,30,60,24,1,CLICK,"Cancel");
  mt_set_callback(obj,CloseSaveWin,0);

  mt_pop_frame();

  mt_set_vis(savewin,1);
}

/* Query Management Window */

static MT_OBJECT *querywin = NULL;
static int query_width  = 400;
static int query_height = BUTTONHIGH + 132;

typedef struct qtag {
  struct qtag *qpred;
  struct qtag *qsucc;
  int          selected;
  int          hue;
  int          sat;
  int          val;
  char        *name;
  char        *query;
  MT_OBJECT   *active;
  MT_OBJECT   *histon;
  MT_OBJECT   *covron;
  MT_OBJECT   *showon;
  MT_OBJECT   *colbox;
  MT_OBJECT   *textbox;
  MT_OBJECT   *frame;
  CoverPacket *cover;
  ObjectPacket *object;
} query_record;

static int querynum = 0;
static query_record *qlist = NULL;

void QuFetchAllNames(void)
{ query_record *q;
  char *name;
  int i, bsel, esel;

  q = qlist;
  for (i = 0; i < querynum; i++)
    { name = mt_get_text(q->textbox,&bsel,&esel);
      if (strcmp(name,q->name) != 0)
        { free(q->name);
          q->name = strdup(name);
        }
      q = q->qsucc;
    }
}

int QuNameEvent(MT_OBJECT *obj, long d)
{ return (0); }

int QueryFileExists(void)
{ struct stat buffer;
  return (stat(CelamyHome,&buffer));
}

void QueryInit(void)
{ FILE *unit;
  query_record *q;
  int   i, c;
  int   sel, hue, sat, val, hon, con, son;
  int   namelen, querylen, size, len;
  char *name = NULL, *query = NULL;
  CoverPacket *cover;
  ObjectPacket *objp;

  unit = fopen(CelamyHome,"r");
  if (unit == NULL) return;
  if (fscanf(unit,"%d %d %d\n",&namelen,&querylen,&size) != 3) goto error;
  name = (char *) malloc(sizeof(char)*namelen+1);
  query = (char *) malloc(sizeof(char)*querylen+1);
  QueryTableSize(size);
  while ((len = fscanf(unit,"%d %d %d %d %d %d %d",
                       &sel,&hue,&sat,&val,&hon,&con,&son)) != EOF)
    { if (len != 7) goto error;
      if (fscanf(unit," %d,",&len) != 1) goto error;
      for (i = 0; i < len; i++)
        if ((c = fgetc(unit)) == EOF)
          goto error;
        else
          name[i] = c;
      name[len] = '\0';
      if (fscanf(unit," %d,",&len) != 1) goto error;
      for (i = 0; i < len; i++)
        if ((c = fgetc(unit)) == EOF)
          goto error;
        else
          query[i] = c;
      query[len] = '\0';
      if (ProcessQuery(query,&objp,&cover) != NULL) goto error;

      q = (query_record *) malloc(sizeof(query_record));
      q->selected = sel;
      if (qlist == NULL)
        qlist = q->qsucc = q->qpred = q;
      else
        { q->qpred = qlist->qpred;
          q->qsucc = qlist;
          q->qpred->qsucc = q->qsucc->qpred = q;
        }
      q->hue    = hue;
      q->sat    = sat;
      q->val    = val;
      q->name   = strdup(name);
      q->query  = strdup(query);
      q->histon = (hon ? Window : NULL);
      q->covron = (con ? Window : NULL);
      q->showon = (son ? Window : NULL);
      q->cover  = cover;
      q->object = objp;
      querynum += 1;

      QueryTableEntry(q->name,cover,objp);
    }
  free(name);
  free(query);
  fclose(unit);
  QueryTableFree();
  return;

error:
  fprintf(stderr,"Invalid query save file: ignored\n");
  while (querynum > 0)
    { q = qlist->qsucc;
      free(q->name);
      free(q->query);
      if (q->cover != NULL)
        FreeCoverage(q->cover);
      if (q->object != NULL)
        FreeObjects(q->object);
      q->qpred->qsucc = q->qsucc;
      q->qsucc->qpred = q->qpred;
      free(q);
      querynum -= 1;
    }
  qlist = NULL;
  free(name);
  free(query);
  fclose(unit);
  QueryTableFree();
  return;
}

void QueryWrite()
{ FILE *unit;
  query_record *q;
  int namelen, querylen, len;

  QuFetchAllNames();
  unit = fopen(CelamyHome,"w");
  if (unit == NULL) return;
  if (querynum == 0)
    fprintf(unit,"%d %d %d\n",0,0,0);
  else
    { q = qlist;
      namelen = querylen = 0;
      do
        { len = strlen(q->name);
          if (len > namelen) namelen = len;
          len = strlen(q->query);
          if (len > querylen) querylen = len;
          q = q->qsucc;
        }
      while (q != qlist);
      fprintf(unit,"%d %d %d\n",namelen,querylen,querynum);
      q = qlist;
      do
        { fprintf(unit,"%d %d %d %d",q->selected,q->hue,q->sat,q->val);
          fprintf(unit," %d",mt_get_button(q->histon));
          fprintf(unit," %d",mt_get_button(q->covron));
          fprintf(unit," %d",mt_get_button(q->showon));
          fprintf(unit," %d,%s %d,%s\n",(int) (strlen(q->name)),q->name,
                                        (int) (strlen(q->query)),q->query);
          q = q->qsucc;
        }
      while (q != qlist);
    }
  fclose(unit);
}

int QuSave(MT_OBJECT *o, long d)
{ if (QueryFileExists() == 0)
    SaveDialog();
  else
    QueryWrite();
  return (0);
}

void TrimStringDraw(MT_OBJECT *obj, int x, int y, char *s, int pix)
{ int v, w, h, b, c, len;

  mt_string_size(s,&w,&h,&b);
  len = strlen(s);
  if (w <= pix)
    { mt_set_color(mt_black());
      mt_draw_text(obj,x,y,s);
      return;
    }

  mt_string_size(" ...",&v,&h,&b);
  pix -= v;
  if (pix < 0) return;
  len = (int)(((1.*pix)/w)*len);

  len -= 1;
  do
    { len += 1;
      c = s[len];
      s[len] = 0;
      mt_string_size(s,&w,&h,&b);
      s[len] = c;
    }
  while (w <= pix);
  
  while (w > pix)
    { len -= 1;
      c = s[len];
      s[len] = 0;
      mt_string_size(s,&w,&h,&b);
      if (w <= pix)
        { mt_set_color(mt_black());
          mt_draw_text(obj,x,y,s);
          mt_draw_text(obj,x+w,y," ...");
        } 
      s[len] = c;
    }
}

static query_record *QColBeingSet;
static MT_OBJECT    *QButBeingSet;

void NewColor(int h, int s, int v)
{ int red, green, blue;
  long base;

  QColBeingSet->hue = h;
  QColBeingSet->sat = s;
  QColBeingSet->val = v;
  HSV_2_RGB(h,v,s,&red,&green,&blue);
  base = mt_get_color(red,green,blue);
  mt_button_colors(QButBeingSet,base,base,
                                -1,-1,-1,-1);
  mt_draw(QButBeingSet);
}

int QColSet(MT_OBJECT *obj, long d)
{ query_record *q;

  q = (query_record *) d;
  QColBeingSet = q;
  QButBeingSet = obj;
  CreateColorWin(q->hue,q->sat,q->val,NewColor);
  return (0);
}

int QuSelect(MT_OBJECT *obj, long d)
{ query_record *q;
  int val;

  q = (query_record *) d;
  val = mt_get_button(obj);
  if (val)
    { mt_set_button(obj,0);
      mt_draw(obj);
    }
  else
    { q->selected = 1;
      q = q->qsucc;
      while ( ! q->selected)
        q = q->qsucc;
      mt_set_button(q->active,1);
      mt_draw(q->active);
      q->selected = 0;
    }
  return (0);
}

int qframelab(MT_OBJECT *obj, long d)
{ query_record *q;
  int xl, xh, yl, yh;

  q = (query_record *) d;
  mt_get_extent(querywin,&xl,&xh,&yl,&yh);
  mt_set_color(mt_black());
  mt_draw_text(obj,220,20,"=");
  TrimStringDraw(obj,232,20,q->query,xh-175);
  return (0);
}

void DrawUnderlinedString(MT_OBJECT *obj, int x, int y, char *string)
{ int w, h, b;

  mt_draw_text(obj,x,y,string);
  mt_string_size(string,&w,&h,&b);
  mt_draw_line(obj,x,y+1,x+w,y+1,1,0);
}

int TitleFrame(MT_OBJECT *obj, long d)
{ long blue;

  blue = mt_get_color(0,0,255);
  mt_set_color(blue);
  DrawUnderlinedString(obj,5,16,"On?");
  DrawUnderlinedString(obj,30,16,"Color");
  DrawUnderlinedString(obj,70,16,"H?");
  DrawUnderlinedString(obj,90,16,"P?");
  DrawUnderlinedString(obj,110,16,"S?");
  DrawUnderlinedString(obj,142,16,"Name");
  DrawUnderlinedString(obj,232,16,"Query");
  return (0);
}
  
void MakeQueryFrame(int pos, query_record *q)
{ MT_OBJECT *obj;
  int xl, xh, yl, yh;
  int red, green, blue;
  long base;

  mt_get_extent(querywin,&xl,&xh,&yl,&yh);
  mt_current_frame(querywin);

  obj = mt_new_frame(0,BUTTONHIGH+28*pos+20,(xh-xl)+1,28,1,
                       REGULAR,REDRAW_EVENT);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_border(obj,0);
  mt_set_outline(obj,0);
  mt_set_callback(obj,qframelab,(long) q);
  q->frame = obj;

    obj = mt_new_button(10,7,14,14,1,RADIO,NULL);
    mt_button_colors(obj,PickColor,mt_medium_grey(),
                         -1,-1,MselColor,-1);
    mt_set_outline(obj,1);
    mt_set_button(obj, ! q->selected);
    mt_set_callback(obj,QuSelect,(long) q);
    q->active = obj;

    obj = mt_new_button(30,3,30,22,1,CLICK,NULL);
    HSV_2_RGB(q->hue,q->val,q->sat,&red,&green,&blue);
    base = mt_get_color(red,green,blue);
    mt_button_colors(obj,base,base,
                         -1,-1,-1,-1);
    mt_set_hilite(obj,0);
    mt_set_callback(obj,QColSet,(long) q);
    q->colbox = obj;

    obj = mt_new_button(70,9,11,11,1,TOGGLE,NULL);
    mt_button_colors(obj,mt_medium_grey(),PickColor,
                         -1,-1,MselColor,-1);
    mt_set_outline(obj,1);
    mt_set_button(obj,q->histon != NULL);
    q->histon = obj;

    obj = mt_new_button(90,9,11,11,1,TOGGLE,NULL);
    mt_button_colors(obj,mt_medium_grey(),PickColor,
                         -1,-1,MselColor,-1);
    mt_set_outline(obj,1);
    mt_set_button(obj,q->covron != NULL);
    q->covron = obj;

    obj = mt_new_button(110,9,11,11,1,TOGGLE,NULL);
    mt_button_colors(obj,mt_medium_grey(),PickColor,
                         -1,-1,MselColor,-1);
    mt_set_outline(obj,1);
    mt_set_button(obj,q->showon != NULL);
    q->showon = obj;

    obj = mt_new_textbox(130,3,80,22,1,"\n\r");
    mt_text_colors(obj,-1,-1,-1,-1,TxSelColor);
    mt_set_text(obj,q->name,0,0);
    q->textbox = obj;
    mt_set_callback(obj,QuNameEvent,(long) 0);

  mt_pop_frame();
}

void QuReferenceSetup(void)
{ query_record *q;

  if (querywin != NULL)
    QuFetchAllNames();
  QueryTableSize(querynum);

  if (querynum > 0)
    { q = qlist;
      do
        { QueryTableEntry(q->name,q->cover,q->object);
          q = q->qsucc;
        }
      while (q != qlist);
    }
}

int QuAdder(MT_OBJECT *obj, long d)
{ char *query, *error;
  int bsel, esel;
  query_record *q;
  CoverPacket  *cover;
  ObjectPacket *objp;

  QuReferenceSetup();

  query = mt_get_text(textpane,&bsel,&esel);
  error = ProcessQuery(query,&objp,&cover);

  QueryTableFree();

  if (error != NULL)
    { mt_set_text(errorpane,error,0,0);
      mt_draw(errorpane);
      free(error);
      return (0);
    }
  
  q = (query_record *) malloc(sizeof(query_record));
  q->selected = 1;
  q->object   = objp;
  q->cover    = cover;
  if (qlist == NULL)
    qlist = q->qsucc = q->qpred = q;
  else
    { q->qpred = qlist->qpred;
      q->qsucc = qlist;
      q->qpred->qsucc = q->qsucc->qpred = q;
    }
  q->hue   = 10;
  q->sat   = 59;
  q->val   = 59;
  q->name  = strdup("");
  q->query = strdup(query);
  q->histon = Window;
  q->covron = Window;
  q->showon = Window;

  MakeQueryFrame(querynum,q);
  querynum += 1;

  mt_draw(q->frame);
  if (querynum > 1)
    QuSelect(q->active,(long) q);

  return (0);
}

int QuDeleter(MT_OBJECT *obj, long d)
{
  int   red, green, blue;
  long  base;
  query_record *q, *p;

  if (querynum == 0) return(0);

  QuFetchAllNames();

  for (q = qlist; ! q->selected; q = q->qsucc)
    ;

  free(q->name);
  free(q->query);
  if (q->object != NULL)
    { ClearQueryColor(q->object);
      DrawAssembly(canvas);
      FreeObjects(q->object);
    }
  if (q->cover != NULL)
    FreeCoverage(q->cover);

  p = q->qsucc;
  if (p == qlist && q != qlist)
    { qlist->selected = 1;
      mt_set_button(qlist->active,0);
      mt_draw(qlist->active);
    }

  while (p != qlist)
    { q->name  = p->name;
      q->query = p->query;
      q->hue   = p->hue;
      q->sat   = p->sat;
      q->val   = p->val;
      q->cover = p->cover;
      q->object = p->object;

      HSV_2_RGB(q->hue,q->val,q->sat,&red,&green,&blue);
      base = mt_get_color(red,green,blue);
      mt_button_colors(q->colbox,base,base,
                                 -1,-1,-1,-1);
      mt_set_button(q->histon,mt_get_button(p->histon));
      mt_set_button(q->covron,mt_get_button(p->covron));
      mt_set_button(q->showon,mt_get_button(p->showon));
      mt_set_text(q->textbox,q->name,0,0);
      mt_draw(q->frame);

      q = p;
      p = p->qsucc;
    }

  q->qsucc->qpred = q->qpred;
  q->qpred->qsucc = q->qsucc;
  mt_free(q->frame);
  free(q);
  querynum -= 1;
  if (querynum == 0)
    qlist = NULL;

  return (0);
}

int QuQuery(MT_OBJECT *obj, long d)
{ query_record *q;
  int  red, green, blue;
  long base;
  void CreateHistoWin(HistoPacket *, char *, int);
  void CreateCoverWin(CoverPacket *, char *, int);

  if (querynum == 0) return (0);

  for (q = qlist; ! q->selected; q = q->qsucc)
    ;

  QuFetchAllNames();

  if (q->cover != NULL)
    { if (mt_get_button(q->covron))
        CreateCoverWin(q->cover,q->query,0);
      if (mt_get_button(q->histon))
        CreateHistoWin(Cover2Hist(q->cover),q->query,1);
    }
  else if (q->object != NULL)
    { if (mt_get_button(q->covron))
        CreateCoverWin(Object2Cover(q->object),q->query,1);
      if (mt_get_button(q->histon))
        CreateHistoWin(Object2Hist(q->object),q->query,0);
      if (mt_get_button(q->showon))
        { HSV_2_RGB(q->hue,q->val,q->sat,&red,&green,&blue);
          base = mt_get_color(red,green,blue);
          SetQueryColor(q->object,base);
          DrawAssembly(canvas);
        }
    }

  return (0);
}

int QuClear(MT_OBJECT *obj, long d)
{ query_record *q;

  if (querynum == 0) return (0);

  for (q = qlist; ! q->selected; q = q->qsucc)
    ;

  if (q->object != NULL)
    { ClearQueryColor(q->object);
      DrawAssembly(canvas);
    }

  return (0);
}

int closequery(MT_OBJECT *o, long d)
{ int xl, xh, yl, yh;
  query_record *q;

  QuFetchAllNames();
  q = qlist;
  if (querynum > 0)
    do
      { if ( ! mt_get_button(q->histon))
          q->histon = NULL;
        if ( ! mt_get_button(q->covron))
          q->covron = NULL;
        if ( ! mt_get_button(q->showon))
          q->showon = NULL;
        q = q->qsucc;
      }
    while (q != qlist);

  mt_get_extent(querywin,&xl,&xh,&yl,&yh);
  query_width  = (xh-xl)+1;
  query_height = (yh-yl)+1;
  ToolIsOn[d] = 0;
  mt_free(querywin); 
  querywin = NULL;
  return (0);
}

void CreateQueryWin(long toolno)
{ MT_OBJECT *obj;
  int i;
  query_record *qr;

  querywin = mt_new_window(10,10,query_width,query_height,0,"Query Manager");
  mt_set_window_bounds(querywin,350,132+BUTTONHIGH,-1,-1);

  obj = mt_new_frame(0,0,query_width,BUTTONHIGH,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_hilite(obj,0);

  { int cornery;

    cornery = (BUTTONHIGH - 24)/2;

    obj = mt_new_button(10,cornery,46,24,1,CLICK,"Add");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,QuAdder,0);

    obj = mt_new_button(66,cornery,46,24,1,CLICK,"Delete");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,QuDeleter,0);

    obj = mt_new_button(122,cornery,46,24,1,CLICK,"Query");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,QuQuery,0);

    obj = mt_new_button(178,cornery,46,24,1,CLICK,"Clear");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,QuClear,0);

    obj = mt_new_button(234,cornery,46,24,1,CLICK,"Save");
    mt_set_xscale(obj,0.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,QuSave,0);

    obj = mt_new_button(query_width-50,cornery,40,24,1,CLICK,"Close");
    mt_set_xscale(obj,1.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,closequery,(long) toolno);
  }
 
  mt_pop_frame();

  obj = mt_new_frame(0,BUTTONHIGH,query_width,20,1,REGULAR,REDRAW_EVENT);
  mt_set_xscale(obj,0.,1.);
  mt_set_border(obj,0);
  mt_set_outline(obj,0);
  mt_set_callback(obj,TitleFrame,(long) 0);
  mt_pop_frame();

  mt_pop_frame();

  qr = qlist;
  for (i = 0; i < querynum; i++)
    { MakeQueryFrame(i,qr);
      qr = qr->qsucc;
    }

  mt_set_vis(querywin,1);
}

/* Histogram Window Management */

typedef struct {
    HistoPacket *hist;
    int          owner;
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
  char label[800];

  p = (HistoLoad *) d;
  DrawHistoInfo(o,p->hist,label);
  return (0);
}

int cbutton(MT_OBJECT *o, long d)
{ HistoLoad *p;

  p = (HistoLoad *) d;
  if (p->owner)
    FreeHistogram(p->hist);
  mt_free(p->window);
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
  
void CreateHistoWin(HistoPacket *hist, char *qmesg, int owner)
{ MT_OBJECT *obj;
  HistoLoad *packet;
  int textwide, histwide, binwide, binhigh;
  char label[800];

  if (hist == NULL) return;

  packet = (HistoLoad *) malloc(sizeof(HistoLoad));
  packet->hist  = hist;
  packet->owner = owner;

  textwide = DrawHistoInfo(NULL,hist,label);
  fprintf(stdout,"%s",label);
  GetBinButtonDims(&binwide,&binhigh);
  histwide = textwide + 72 + binwide;
  if (histwide < HISTWIDE) histwide = HISTWIDE;

  packet->window = obj =
            mt_new_window(10,10,histwide,
                          HISTHIGH+HISTBUTTON+2*SCROLLWIDE,0,qmesg);
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
      tab = (int)(high*.05);
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

    obj = mt_new_button(histwide-50,cornery,40,24,1,CLICK,"Close");
    mt_set_xscale(obj,1.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,cbutton,(long) packet);
  }
 
  mt_pop_frame();

  mt_pop_frame();

  mt_set_vis(packet->window,1);
}

/* Coverage Window Routines */

typedef struct CoverTag {
    struct CoverTag *next, *prev; 
    CoverPacket *cover;
    MT_OBJECT   *darea;
    MT_OBJECT   *window;
    int          owner;
} CoverLoad;

static CoverLoad *CoverList = NULL;

int cover_handler(MT_OBJECT *o, long d)
{ CoverLoad *p;
  int val0, wid0, min, max;

  p = (CoverLoad *) d;
  val0 = mt_get_scrollbar(sbars[HSCROLL].bar,&wid0,&min,&max);
  DrawCoverage(p->darea,p->cover,val0+1,val0+wid0);
  return (0);
}

int cobutton(MT_OBJECT *o, long d)
{ CoverLoad *p;

  p = (CoverLoad *) d;
  if (p->owner)
    FreeCoverage(p->cover);
  mt_free(p->window);
  if (p->prev != NULL)
    p->prev->next = p->next;
  else
    CoverList = p->next;
  if (p->next != NULL)
    p->next->prev = p->prev;
  free(p);
  return (0);
}

void CreateCoverWin(CoverPacket *cover, char *qmesg, int owner)
{ MT_OBJECT *obj;
  CoverLoad *packet;
  int xl = 0, yh = 0;

  packet = (CoverLoad *) malloc(sizeof(CoverLoad));
  packet->cover   = cover;
  packet->owner   = owner;
  packet->next    = CoverList;
  if (CoverList != NULL)
    CoverList->prev = packet;
  packet->prev    = NULL;
  CoverList       = packet;

  packet->window = obj =
           mt_new_window(xl,yh,MainWidth,COVERHIGH+BUTTONHIGH,0,qmesg);
  mt_set_window_bounds(obj,MainWidth,100,MainWidth,-1);

  packet->darea = obj =
        mt_new_frame(0,BUTTONHIGH,MainWidth,COVERHIGH,1,REGULAR,REDRAW_EVENT);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,1.);
  mt_set_outline(obj,0);
  mt_set_border(obj,0);
  mt_frame_colors(obj,BackGround,-1,-1,-1,-1);
  mt_set_callback(obj,cover_handler,(long) packet);
  mt_pop_frame();

  obj = mt_new_frame(0,0,MainWidth,BUTTONHIGH,1,REGULAR,0);
  mt_set_xscale(obj,0.,1.);
  mt_set_yscale(obj,0.,0.);
  mt_set_hilite(obj,0);

  { int cornery;

    cornery = (BUTTONHIGH - 24)/2;

    obj = mt_new_button(MainWidth-50,cornery,40,24,1,CLICK,"Close");
    mt_set_xscale(obj,1.,0.);
    mt_set_yscale(obj,0.,0.);
    mt_button_colors(obj,-1,SelectColor,-1,-1,SelectColor,-1);
    mt_set_callback(obj,cobutton,(long) packet);
  }
 
  mt_pop_frame();

  mt_pop_frame();

  mt_set_vis(packet->window,1);
}

void ResizeAllCovers(int resize)
{ CoverLoad *p;

  for (p = CoverList; p != NULL; p = p->next)
    if (resize)
      mt_set_window_bounds(p->window,MainWidth,-1,MainWidth,-1);
    else
      mt_draw(p->darea);
}
