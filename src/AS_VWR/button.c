
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
/* $Id: button.c,v 1.2 2004-09-23 20:25:30 mcschatz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "wpack.h"

#define FRAME     0  /* Widget Types */
#define BUTTON    1
#define SCROLLBAR 2
#define BITMAP    3
#define TEXTBOX   4

#define WINDOW  0 /* Frame Classes */
#define REGULAR 1
#define VIRTUAL 2

#define CLICK  0  /* Button Classes */
#define PRESS  1
#define TOGGLE 2
#define RADIO  3

#define IDLE      0  /* Window Event States */
#define PRESSED   1
#define SLIDER    2
#define MENUPRESS 3
#define MENUDOWN  4
#define MENUTWICE 5
#define TEXTSEL   6

static long greyf, greyh, greys, cyan;  /* Default widget colors */

typedef struct obj { struct obj   *forw, *back;
                     int           type;
                     int           vis;
                     window_desc   wpwin;
                     struct obj   *frame;
                     int           level;
                     int           xl, xh;
                     int           yl, yh;
                     int           ox, ow;
                     int           oy, oh;
                     int           minx, miny;
                     float         xpf, xof;
                     float         ypf, yof;
                     int           border, outline, hilite;
                     int           pixlab;
                     char         *label;
                     void        (*free_routine)();
                     void        (*draw_routine)();
                     int         (*event_routine)();
                     long          event_data;
                   } MT_OBJECT;

#define REDRAW_EVENT  0x01 /* Event Chord */
#define PRESS_EVENT   0x02
#define RELEASE_EVENT 0x04
#define KEY_EVENT     0x08
#define ENTRY_EVENT   0x10
#define EXIT_EVENT    0x20

#define SHIFT_MODIFIER 0x1
#define CNTRL_MODIFIER 0x2
#define META_MODIFIER  0x4

typedef struct frame { MT_OBJECT     basis;
                       int           class;
                       MT_OBJECT    *ctlbut;
                       int           event_chord;
                       MT_OBJECT    *afor, *abck;
                       long          col_dark;
                       long          col_lite;
                       long          col_up;
                       long          col_dsel;
                       long          col_lsel;
                     } frame;

typedef struct window { frame       basis;
                        MT_OBJECT  *nxtdial;
                        int         state;
                        int         but_dwn;
                        int         slide_delta;
                        long        dwn_time;
                        MT_OBJECT  *menuf;
                        MT_OBJECT  *menub;
                        MT_OBJECT  *object;
                        bitmap_desc backing;
                        int         xdown;
                        int         xlast;
                        int         started;
                        int         recentime;
                        MT_OBJECT  *focus;
                        int         minwide;
                        int         minhigh;
                        int         maxwide;
                        int         maxhigh;
                      } window;

typedef struct textinput { frame    basis;
                           char    *text;
                           int      bufsize;
                           int      bsel, esel;
                           int      drawn;
                           long     lastev;
                           char     keyact[256];
                         } textinput;

typedef struct { MT_OBJECT   basis;
                 int         class;
                 int         status;
		 long        col_dark;
		 long        col_lite;
                 long        col_up;
                 long        col_down;
                 long        col_sel;
                 long        col_text;
                 bitmap_desc map_down;
                 bitmap_desc map_sel;
                 MT_OBJECT  *menu;
               } button;

typedef struct { button      basis;
                 int         xm, ym;
                 int         radius;
               } radio;

typedef struct { MT_OBJECT   basis;
                 int         value, width;
                 int         pxpos,  pxwid;
                 int         minim, maxim;
                 int         vertical;
                 long        col_dark;
                 long        col_lite;
                 long        col_tab;
                 long        col_well;
                 long        col_sel;
               } scrollbar;

/* ATTRIBUTE VARS AND ROUTINES */

static void erase_frame(MT_OBJECT *);
static void erase_button(MT_OBJECT *);
static void erase_scrollbar(MT_OBJECT *);
static void scroll_v2p(scrollbar *);
static void bp_handler(input_event *, window_desc, void *);
static void scale_adjust(MT_OBJECT *, int, int);

static int current_border;  /* Current Attribute Vars */
static int current_outline;
static int current_hilite;
static float current_xpf, current_xof;
static float current_ypf, current_yof;
static MT_OBJECT *current_frame;

static int ev_x, ev_y;  /* Event packet for frame events */
static int ev_modifiers;
static int ev_value;
static int ev_type;
static MT_OBJECT *dialog_stack = NULL;

static font_desc current_font;

static MT_OBJECT nonobject, *no_object;

int mt_get_event(MT_OBJECT *o, int *x, int *y, int *value, int *mods)
{ *x = ev_x - o->xl;
  *y = ev_y - o->yl;
  *value = ev_value;
  *mods  = ev_modifiers;
  return (ev_type);
}

static void init_attributes(void)
{ nonobject.level = -1;
  no_object = &nonobject;
  current_xpf = 0.;
  current_xof = 0.;
  current_ypf = 0.;
  current_yof = 0.;
  current_frame = no_object;
  current_border = 2;
  current_outline = 1;
  current_hilite = 1;
  current_font = font_default();
}

void mt_current_xscale(double xpf, double xof)
{ current_xpf = xpf;
  current_xof = xof;
}

void mt_current_yscale(double ypf, double yof)
{ current_ypf = ypf;
  current_yof = yof;
}

void mt_current_border(int width)
{ current_border = width; }

void mt_current_outline(int is_on)
{ current_outline = is_on; }

void mt_current_frame(MT_OBJECT *fr)
{ current_frame = fr; }

void mt_pop_frame(void)
{ if (current_frame != no_object)
    current_frame = current_frame->frame;
}

void mt_current_hilite(int is_on)
{ current_hilite = is_on; }


/* GENERAL OBJECT ATTRIBUTE SETTING ROUTINES */

static void set_minimums(MT_OBJECT *o)  /* Button specific? */
{ int wid, hgt, base;

  if (o->label != NULL)
    if (o->pixlab)
      { MT_OBJECT *pix;
        pix = (MT_OBJECT *) (o->label);
        wid = pix->ow;
        hgt = pix->oh;
        base = 0;
      }
    else
      font_string_size(current_font,o->label,&wid,&hgt,&base);
  else
    wid = hgt = base = 0;
  o->minx = 2*(o->border + o->outline) + wid;
  o->miny = 2*(o->border + o->outline) + hgt + base;
  if (o->type == SCROLLBAR)
    { o->minx += 2*o->border;
      o->miny += 2*o->border;
    }
}
  
void mt_set_callback(MT_OBJECT *o,
                     int (*call)(MT_OBJECT *p, long d), long data)
{ o->event_routine = call;
  o->event_data    = data;
}

void mt_set_xscale(MT_OBJECT *o, double xpf, double xof)
{ o->xpf = xpf;
  o->xof = xof;
}

void mt_set_yscale(MT_OBJECT *o, double ypf, double yof)
{ o->ypf = ypf;
  o->yof = yof;
}

void mt_set_border(MT_OBJECT *o, int width)
{ o->border = width;
  set_minimums(o);
}

void mt_set_outline(MT_OBJECT *o, int is_on)
{ if (is_on != o->outline)
    { o->outline = is_on;
      if (o->type == SCROLLBAR)
        scroll_v2p((scrollbar *) o);
    }
}

void mt_set_hilite(MT_OBJECT *o, int is_on)
{ o->hilite = is_on; }

char *mt_bitmap_label(MT_OBJECT *o)
{ static char StringEncode[50];
  int i, cnt;
  unsigned long x;

  if (o->type != BITMAP)
    { fprintf(stderr,"Object is not a bitmap (mt_bitmap_label)\n");
      exit (1);
    }
  StringEncode[0] = '\01';
  x = (unsigned long) o;
  for (cnt = 0; x != 0; cnt++)
    x /= 10;
  StringEncode[cnt] = '\0';
  x = (unsigned long) o;
  for (i = cnt; i > 0; i--)
    { StringEncode[i] = (x%10) + '0';
      x /= 10;
    }
  return (StringEncode);
}

static MT_OBJECT *label_test(char *label)
{ unsigned long x;
  int i;

  if (label == NULL) return (NULL);
  if (*label != '\01') return (NULL);
  x = 0;
  for (i = 1; label[i] != '\0'; i++)
    x = 10*x + (label[i] - '0');
  return ((MT_OBJECT *) x);
}

void mt_set_label(MT_OBJECT *o, char *label)
{ MT_OBJECT *bitlabel;

  bitlabel = label_test(label); 
  if (bitlabel != NULL)

    { if (!o->pixlab && o->label != NULL)
        free(o->label);

      if (o->type == BUTTON)
        { button *b;
          bitmap_desc bm;

          b = (button *) o;
          if (o->pixlab)
            { free(b->map_sel);
              free(b->map_down);
            }
          bm = (bitmap_desc) (bitlabel->wpwin);
          b->map_sel  = bm_new(bitlabel->ow,bitlabel->oh);
          b->map_down = bm_new(bitlabel->ow,bitlabel->oh);
          bm_draw_bitmap(b->map_sel,0,0,bm,MODE_SRC);
          bm_draw_bitmap(b->map_down,0,0,bm,MODE_SRC);
          bm_map_pixels(b->map_sel,b->col_up,b->col_sel);
          bm_map_pixels(b->map_down,b->col_up,b->col_down);
        }

      o->pixlab = 1;
      o->label  = (char *) bitlabel;
      set_minimums(o);
      return;
    }

  if (o->pixlab)
    { if (label != NULL)
        o->label = strcpy((char *) malloc(strlen(label)+1),label);
      else
        o->label = NULL;

      if (o->type == BUTTON)
        { bm_free(((button *) o)->map_sel);
          bm_free(((button *) o)->map_down);
        }

      o->pixlab = 0;
      set_minimums(o);
      return;
    }

  if (o->label == NULL && label == NULL)
    return;
  if (o->label != NULL && label != NULL && strlen(label) < strlen(o->label))
    { strcpy(o->label,label);
      set_minimums(o);
      return;
    }
  if (o->label != NULL)
    free(o->label);
  if (label != NULL)
    o->label = strcpy((char *) malloc(strlen(label)+1),label);
  else
    o->label = NULL;
  set_minimums(o);
}

void mt_set_vis(MT_OBJECT *o, int vis)
{ if (vis != o->vis)
    { o->vis = vis;
      if (o->type == FRAME && ((frame *) o)->ctlbut != NULL)
        { window     *root;
          window_desc win;
          int         s;

          root = (window *) (o->frame);
          win  = o->wpwin;
          s = root->state;
          if (vis)
            { if (s == MENUDOWN || s == MENUPRESS || s == MENUTWICE)
                { fprintf(stderr,"Activating menu while one is already down");
                  fprintf(stderr," (mt_set_vis)\n");
                  exit (1);
                }
              root->menuf = o;
              root->menub = ((frame *) o)->ctlbut;
              root->backing = bm_new(o->ow,o->oh);
              wp_win_saveblk(root->backing,0,0,o->ow,o->oh,
                             win,o->xl,o->yl,MODE_SRC);
              o->draw_routine(o,0);
              root->state = MENUDOWN;
            }
          else
            { win_draw_bitmap(win,root->menuf->xl,root->menuf->yl,
                              root->backing,MODE_SRC);
              bm_free(root->backing);
              root->object  = no_object;
              root->state   = IDLE;
            }
        }
      else
        if (vis)
          o->draw_routine(o,0);
        else
          switch (o->type)
          { case TEXTBOX:
            case FRAME:
              erase_frame(o);
              break;
            case BUTTON:
              erase_button(o);
              break;
            case SCROLLBAR:
              erase_scrollbar(o);
              break;
          }
    }
  if (o->type == FRAME && ((frame *) o)->class == WINDOW)
    win_set_vis(o->wpwin,vis);
}

void mt_get_extent(MT_OBJECT *o, int *xl, int *xh, int *yl, int *yh)
{ *xl = o->xl;
  *xh = o->xh;
  *yl = o->yl;
  *yh = o->yh;
}

int mt_get_vis(MT_OBJECT *o)
{ return (o->vis); }

void mt_draw(MT_OBJECT *o)
{ o->draw_routine(o,0); }

void mt_free(MT_OBJECT *o)
{
  mt_set_vis(o,0);
  o->forw->back = o->back;
  o->back->forw = o->forw;
  if (o->free_routine != NULL) o->free_routine(o);
  free(o);
}


/* FRAME OBJECT ROUTINES */

static void free_frame(MT_OBJECT *a)
{ MT_OBJECT *o, *p, *q;
  frame     *f;
  window    *w;

  f = (frame *) a;
  if (f->class == WINDOW)
    { if (dialog_stack != NULL)
        { if (dialog_stack == a)
            dialog_stack = ((window *) a)->nxtdial;
          else
            { fprintf(stderr,"Freeing window that is not at top of");
              fprintf(stderr," dialog stack\n");
              exit (1);
            }
        }
    }
  p = (MT_OBJECT *) ( & (f->afor));
  for (o = p->forw; o != p; o = q)
    { q = o->forw;
      if (o->free_routine != NULL) o->free_routine(o);
      free(o);
    }
  if (f->class == WINDOW)
    win_free(a->wpwin);
  else if (a->type == TEXTBOX)
    { p = a;
      while (p->frame != no_object)
        p = p->frame;
      w = (window *) p;
      if (a == w->focus)
        w->focus = no_object;
    }
}

static void draw_frame(MT_OBJECT *o, int high)
{ MT_OBJECT *p;
  frame     *fr;
  window_desc win;
  int out, cut;
  int uni, all, wid;
  int lite, dark;
  register int xl, xh, yl, yh;
  register int  i,  w;

  fr = (frame *) o;
  if (fr->class != VIRTUAL)
    { xl  = o->xl;
      xh  = o->xh;
      yl  = o->yl;
      yh  = o->yh;
      w   = o->border;
      win = o->wpwin;
      out = o->outline;
      cut = 1-out;

      if (high > 0)
        { lite = fr->col_lsel; dark = fr->col_dsel; }
      else
        { lite = fr->col_lite; dark = fr->col_dark; }
    
      uni = (lite == dark);
      all = (fr->col_lite == fr->col_up);
    
      if (uni && all && (high == 0))
        { mt_set_color(fr->col_lite);
          wp_win_fill(win,xl+out,yl+out,xh+cut,yh+cut,MODE_SRC);
        }
      else
        { if (uni)
            { mt_set_color(lite);
              wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
              wp_win_fill(win,xl+out,yl+out,xl+out+w,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xl+w+out,yl+out,xh+cut-w,yl+out+w,MODE_SRC);
            }
          else
            { mt_set_color(dark);
              wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
              mt_set_color(lite);
              wid = w - cut;
              for (i = out; i <= wid; i++)
                { win_draw_line(win,xl+out,yl+i,xh-1-i,yl+i,MODE_SRC);
                  win_draw_line(win,xl+i,yl+out,xl+i,yh-1-i,MODE_SRC);
                }
            }
    
          if (high == 0)
            { mt_set_color(fr->col_up);
              wp_win_fill(win,xl+out+w,yl+out+w,xh+cut-w,yh+cut-w,MODE_SRC);
            }
        }
      
      if (out && (high == 0))
        { mt_set_color(mt_black());
          win_draw_line(win,xl,yl,xl,yh,MODE_SRC);
          win_draw_line(win,xh,yl,xh,yh,MODE_SRC);
          win_draw_line(win,xl,yl,xh,yl,MODE_SRC);
          win_draw_line(win,xl,yh,xh,yh,MODE_SRC);
        }

      if (high == 0 && (fr->event_chord & REDRAW_EVENT) != 0
                    && o->event_routine != NULL)
        { fr->event_chord &= ~ REDRAW_EVENT;
          ev_type = REDRAW_EVENT;
          o->event_routine(o,o->event_data);
          fr->event_chord |= REDRAW_EVENT;
        }
    }

  if (high == 0)
    { if ((fr->event_chord & REDRAW_EVENT) != 0 && o->event_routine != NULL)
        { fr->event_chord &= ~ REDRAW_EVENT;
          ev_type = REDRAW_EVENT;
          o->event_routine(o,o->event_data);
          fr->event_chord |= REDRAW_EVENT;
        }
      p = (MT_OBJECT *) ( & (fr->afor));
      for (o = p->forw; o != p; o = o->forw)
        if (o->vis)
          o->draw_routine(o,0);
    }
}

static void erase_frame(MT_OBJECT *o)
{ frame *fr; 

  fr = (frame *) o;
  if (fr->class == WINDOW)
    mt_set_color(fr->col_up);
  else if (fr->class == REGULAR)
    mt_set_color(((frame *) o->frame)->col_up);
  if (fr->class != VIRTUAL)
    wp_win_fill(o->wpwin,o->xl,o->yl,o->xh,o->yh,MODE_SRC);
}

MT_OBJECT *mt_new_window(int x, int y, int w, int h, int vis, char *title)
{ frame     *wfr;
  window    *wwn;
  MT_OBJECT *win;

  wwn = (window *) malloc(sizeof(window));
  wfr = (frame *) wwn;
  win = (MT_OBJECT *) wwn;

  wwn->state   = IDLE;
  wwn->object  = no_object;
  wwn->focus   = no_object;
  wwn->but_dwn = -1;
  wwn->started = 0;
  wwn->minwide = 0;
  wwn->minhigh = 0;
  wwn->maxwide = 0x7FFFFFFF;
  wwn->maxhigh = 0x7FFFFFFF;

  win->forw  = win->back = win;

  win->type  = FRAME;
  win->vis   = vis;
  win->wpwin = win_new(x,y,w,h,bp_handler,win,title); 
  win->frame = no_object;
  win->level = 0;

  win->ox = x;
  win->oy = y;
  win->xl = 0;
  win->yl = 0;
  win->xh = (win->ow = w) - 1;
  win->yh = (win->oh = h) - 1;

  win->xpf = 0.;
  win->xof = 1.;
  win->ypf = 0.;
  win->yof = 1.;

  win->border  = 0;
  win->outline = 0;
  win->hilite  = 0;
  win->pixlab  = 0;
  win->label   = NULL;

  win->free_routine  = free_frame;
  win->draw_routine  = draw_frame;
  win->event_routine = NULL;

  wfr->class  = WINDOW;
  wfr->ctlbut = NULL;
  wfr->afor   = wfr->abck = (MT_OBJECT *) ( & (wfr->afor));

  wfr->col_dark = greys;
  wfr->col_lite = greyh;
  wfr->col_up   = greyf;
  wfr->col_dsel = greyh;
  wfr->col_lsel = greys;

  if (vis)
    { win->draw_routine(win,0);
      win_set_vis(win->wpwin,1);
    }

  current_frame = win;
  return (current_frame);
}

MT_OBJECT *mt_new_dialog(int x, int y, int w, int h, int vis, char *title)
{ window    *wwn;
  MT_OBJECT *win;

  win = mt_new_window(x,y,w,h,vis,title);
  wwn = (window *) win;
  wwn->nxtdial = dialog_stack;
  dialog_stack = win;
  return (current_frame);
}

void mt_set_window_bounds(MT_OBJECT *o, int minw, int minh, int maxw, int maxh)
{ window *win;
  int w, h;

  if (o->type != FRAME || ((frame *) o)->class != WINDOW)
    { fprintf(stderr,"Object is not a window (mt_set_window_bounds)\n");
      exit (1);
    }
  win = (window *) o;
  if (minw >= 0) win->minwide = minw;
  if (minh >= 0) win->minhigh = minh;
  if (maxw >= 0) win->maxwide = maxw;
  if (maxh >= 0) win->maxhigh = maxh;

  if (win->minwide > win->maxwide)
    { fprintf(stderr,"Fatal: Window min width (%d) > max width (%d)\n",
                     win->minwide,win->maxwide);
      exit (1);
    }
  if (win->minhigh > win->maxhigh)
    { fprintf(stderr,"Fatal: Window min height (%d) > max height (%d)\n",
                     win->minhigh,win->maxhigh);
      exit (1);
    }

  w = o->xh - o->xl + 1;
  h = o->yh - o->yl + 1;
  if (win->maxwide < w || w < win->minwide ||
      win->maxhigh < h || h < win->minhigh)
    { if (win->maxwide < w)
        w = win->maxwide;
      else if (w < win->minwide)
        w = win->minwide;
      if (win->maxhigh < h)
        h = win->maxhigh;
      else if (h < win->minhigh)
        h = win->minhigh;
      win_set_size(o->wpwin,w,h);
    }
}

void mt_set_window_size(MT_OBJECT *o, int width, int height)
{ window *win;

  if (o->type != FRAME || ((frame *) o)->class != WINDOW)
    { fprintf(stderr,"Object is not a window (mt_set_window_bounds)\n");
      exit (1);
    }
  win = (window *) o;
  win_set_size(o->wpwin,width,height);
  scale_adjust(o,width,height);
  win_set_clip(o->wpwin,0,0,width,height);
  o->xh = width-1;
  o->yh = height-1;
  if (win->state == MENUDOWN)
    win->menuf->vis = 0;
  o->draw_routine(o,0);
  if (win->state == MENUDOWN)
    { wp_win_saveblk(win->backing,0,0,
                     win->menuf->ow,win->menuf->oh,
                     o->wpwin,win->menuf->xl,win->menuf->yl,MODE_SRC);
      win->menuf->vis = 1;
      win->menuf->draw_routine(win->menuf,0);
    }
  win->object = no_object;
}

static MT_OBJECT *new_frame_fill(frame *f, int x, int y, int w, int h,
                                 int vis, int class, int ev_chord)
{ MT_OBJECT *o, *fr;
  MT_OBJECT *p;

  fr = current_frame;
  if (fr == no_object)
    { fprintf(stderr,"Fatal: No current frame (mt_new_frame)\n");
      exit (1);
    }

  o = (MT_OBJECT *) f;

  p = (MT_OBJECT *) ( & (((frame *) fr)->afor));
  o->back = p->back;
  o->forw = p;
  o->back->forw = o->forw->back = o;

  o->type  = FRAME;
  o->vis   = vis;
  o->wpwin = fr->wpwin;
  o->frame = fr;
  o->level = fr->level + 1;

  o->xl = fr->xl + x;
  o->yl = fr->yl + y;
  o->xh = o->xl + (w-1);
  o->yh = o->yl + (h-1);
  o->ow = w;
  o->oh = h;
  o->ox = x;
  o->oy = y;

  o->xpf = current_xpf;
  o->xof = current_xof;
  o->ypf = current_ypf;
  o->yof = current_yof;

  o->border  = current_border;
  o->outline = current_outline;
  o->hilite  = current_hilite;
  o->pixlab  = 0;
  o->label   = NULL;

  o->free_routine  = free_frame;
  o->draw_routine  = draw_frame;
  o->event_routine = NULL;
  f->event_chord   = ev_chord;

  f->class  = class;
  f->ctlbut = NULL;
  f->afor   = f->abck = (MT_OBJECT *) ( & (f->afor));

  f->col_dark = greys;
  f->col_lite = greyh;
  if (class == VIRTUAL)
    f->col_up = ((frame *) fr)->col_up;
  else
    f->col_up = greyf;
  f->col_dsel = greyh;
  f->col_lsel = greys;

  if (o->vis) o->draw_routine(o,0);

  current_frame = o;
  return (current_frame);
}

MT_OBJECT *mt_new_frame(int x, int y, int w, int h,
                        int vis, int class, int ev_chord)
{ frame *f;

  f = (frame *) malloc(sizeof(frame));
  return (new_frame_fill(f,x,y,w,h,vis,class,ev_chord));
}

void mt_frame_colors(MT_OBJECT *o, long fore, long dark, long lite,
                                   long dsel, long lsel)
{ frame *f;

  f = (frame *) o;
  if (o->type != FRAME)
    { fprintf(stderr,"Object is not a frame (mt_frame_colors)\n");
      exit (1);
    }

  if (dark >= 0)
    f->col_dark = dark;
  if (lite >= 0)
    f->col_lite = lite;
  if (fore >= 0)
    f->col_up = fore;
  if (dsel >= 0)
    f->col_dsel = dsel;
  if (lsel >= 0)
    f->col_lsel  = lsel;
}


/* BUTTON ROUTINES */

static void set_radio_geom(radio *r)
{ register int xl, xh, yl, yh;
  MT_OBJECT *o;

  o = (MT_OBJECT *) r;
  xl = o->xl;
  xh = o->xh;
  yl = o->yl;
  yh = o->yh;
  r->xm = (xl+xh)/2;
  r->ym = (yl+yh)/2;
  if (xh-xl > yh-yl)
    r->radius = r->ym-yl;
  else
    r->radius = r->xm-xl;
}

static void draw_radio(MT_OBJECT *o, int high)
{ window_desc win;
  button *b;
  int status, out, cut;
  int wid, hgt, base;
  register int xl, xh, yl, yh;
  register int  i,  w;
  register int xm, ym, rad;

  b = (button *) o;

  w   = o->border;
  win = o->wpwin;
  out = o->outline;
  cut = 1-out;

  status = b->status;

  xm  = ((radio *) b)->xm;
  ym  = ((radio *) b)->ym;
  rad = ((radio *) b)->radius;
  xl  = xm - rad;
  yl  = ym - rad;
  xh  = xm + rad;
  yh  = ym + rad;

  if (high == 0)

    { mt_set_color(b->col_dark);
      for (i = w-cut; i >= out; i--)
        { win_draw_line(win,xm+1,yh-i-1,xh-i,ym,MODE_SRC);
          win_draw_line(win,xh-i,ym,xm+1,yl+i+1,MODE_SRC);
        }
  
      mt_set_color(b->col_lite);
      for (i = w-cut; i >= out; i--)
        { win_draw_line(win,xm,yl+i,xl+i,ym,MODE_SRC);
          win_draw_line(win,xl+i,ym,xm,yh-i,MODE_SRC);
        }

      if (status)
        mt_set_color(b->col_down);
      else
        mt_set_color(b->col_up);
      for (i = w+out; i <= rad; i++)
        { win_draw_line(win,xm,yh-i,xh-i,ym,MODE_SRC);
          win_draw_line(win,xh-i,ym,xm,yl+i,MODE_SRC);
          win_draw_line(win,xm,yl+i,xl+i,ym,MODE_SRC);
          win_draw_line(win,xl+i,ym,xm,yh-i,MODE_SRC);
        }

      if (out)
        { mt_set_color(mt_black());
          win_draw_line(win,xl,ym,xm,yh,MODE_SRC);
          win_draw_line(win,xm,yh,xh,ym,MODE_SRC);
          win_draw_line(win,xh,ym,xm,yl,MODE_SRC);
          win_draw_line(win,xm,yl,xl,ym,MODE_SRC);
        }
    }

  else
    { if (high > 0)
        mt_set_color(b->col_sel);
      else if (status)
        mt_set_color(b->col_down);
      else
        mt_set_color(b->col_up);
      for (i = w+out; i <= rad; i++)
        { win_draw_line(win,xm,yh-i,xh-i,ym,MODE_SRC);
          win_draw_line(win,xh-i,ym,xm,yl+i,MODE_SRC);
          win_draw_line(win,xm,yl+i,xl+i,ym,MODE_SRC);
          win_draw_line(win,xl+i,ym,xm,yh-i,MODE_SRC);
        }
    }
    
  if (o->label != NULL)
    if (o->pixlab)
      { bitmap_desc bm;
        MT_OBJECT *map;
        
        map = (MT_OBJECT *) (o->label);
        if (high > 0)
          bm = b->map_sel;
        else if (status)
          bm = b->map_down;
        else
          bm = (bitmap_desc) (map->wpwin);
        win_draw_bitmap(win,xm-((map->ow-1)/2),ym-((map->oh-1)/2),bm,MODE_SRC);
      }
    else
      { mt_set_color(b->col_text);
        font_string_size(current_font,o->label,&wid,&hgt,&base);
        win_draw_string(win,xm-wid/2+1,ym+hgt/2-1,o->label,
                            current_font,MODE_SRC);
      }
}

static void draw_button(MT_OBJECT *o, int high)
{ window_desc win;
  button *b;
  int status, out, cut;
  int uni, all;
  int wid, hgt, base;
  register int xl, xh, yl, yh;
  register int  i,  w;
  register int xm, ym;

  b = (button *) o;

  xl  = o->xl;
  xh  = o->xh;
  yl  = o->yl;
  yh  = o->yh;
  xm  = (xl+xh)/2;
  ym  = (yl+yh)/2;
  w   = o->border;
  win = o->wpwin;
  out = o->outline;
  cut = 1-out;

  status = b->status;

  if (high == 0)
  
    { uni = (b->col_lite == b->col_dark);
      if (status)
        all = (b->col_lite == b->col_down);
      else
        all = (b->col_lite == b->col_up);
    
      if (uni && all)
        { mt_set_color(b->col_lite);
          wp_win_fill(win,xl+out,yl+out,xh+cut,yh+cut,MODE_SRC);
        }
      else
        { if (uni)
            { mt_set_color(b->col_lite);
              wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
              wp_win_fill(win,xl+out,yl+out,xl+out+w,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xl+w+out,yl+out,xh+cut-w,yl+out+w,MODE_SRC);
            }
          else
            { if (status)
                mt_set_color(b->col_lite);
              else
                mt_set_color(b->col_dark);
              wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
              wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
              if (status)
                mt_set_color(b->col_dark);
              else
                mt_set_color(b->col_lite);
              wid = w - cut;
              for (i = out; i <= wid; i++)
                { win_draw_line(win,xl+out,yl+i,xh-1-i,yl+i,MODE_SRC);
                  win_draw_line(win,xl+i,yl+out,xl+i,yh-1-i,MODE_SRC);
                }
            }
    
          if (status)
            mt_set_color(b->col_down);
          else
            mt_set_color(b->col_up);
          wp_win_fill(win,xl+out+w,yl+out+w,xh+cut-w,yh+cut-w,MODE_SRC);
        }
    
      if (out)
        { mt_set_color(mt_black());
          win_draw_line(win,xl,yl,xl,yh,MODE_SRC);
          win_draw_line(win,xh,yl,xh,yh,MODE_SRC);
          win_draw_line(win,xl,yl,xh,yl,MODE_SRC);
          win_draw_line(win,xl,yh,xh,yh,MODE_SRC);
        }
    }

  else
    { if (high > 0)
        mt_set_color(b->col_sel);
      else if (status)
        mt_set_color(b->col_down);
      else
        mt_set_color(b->col_up);
      wp_win_fill(win,xl+out+w,yl+out+w,xh+cut-w,yh+cut-w,MODE_SRC);
    }

  if (o->label != NULL)
    if (o->pixlab)
      { bitmap_desc bm;
        MT_OBJECT *map;
        
        map = (MT_OBJECT *) (o->label);
        if (high > 0)
          bm = b->map_sel;
        else if (status)
          bm = b->map_down;
        else
          bm = (bitmap_desc) (map->wpwin);
        win_draw_bitmap(win,xm-((map->ow-1)/2),ym-((map->oh-1)/2),bm,MODE_SRC);
      }
    else
      { mt_set_color(b->col_text);
        font_string_size(current_font,o->label,&wid,&hgt,&base);
        win_draw_string(win,xm-wid/2+1,ym+hgt/2-1,o->label,
                            current_font,MODE_SRC);
      }
}

static void erase_button(MT_OBJECT *o)
{ register window_desc win;
  register button *b;
  register radio  *r;
  register int d;

  b = (button *) o;
  win = o->wpwin;
  mt_set_color(((frame *) o->frame)->col_up);
  if (b->class == RADIO)
    { r = (radio *) o;
      d = r->radius;
      wp_win_fill(win,r->xm-d,r->ym-d,r->xm+d+1,r->ym+d+1,MODE_SRC);
    }
  else
    wp_win_fill(win,o->xl,o->yl,o->xh+1,o->yh+1,MODE_SRC);
}

MT_OBJECT *mt_new_button(int x, int y, int w, int h, int vis,
                         int class, char *label)
{ button    *b;
  MT_OBJECT *fr, *o;
  MT_OBJECT *p;

  fr = current_frame;
  if (fr == no_object)
    { fprintf(stderr,"Fatal: No current frame (mt_new_button)\n");
      exit (1);
    }

  if (class == RADIO)
    b = (button *) malloc(sizeof(radio));
  else
    b = (button *) malloc(sizeof(button));
  o = (MT_OBJECT *) b;

  p = (MT_OBJECT *) ( & (((frame *) fr)->afor));
  o->back = p->back;
  o->forw = p;
  o->back->forw = o->forw->back = o;

  o->type  = BUTTON;
  o->vis   = vis;
  o->wpwin = fr->wpwin;
  o->frame = fr;
  o->level = fr->level + 1;

  o->xl = fr->xl + x;
  o->yl = fr->yl + y;
  o->xh = o->xl + (w-1);
  o->yh = o->yl + (h-1);
  o->ow = w;
  o->oh = h;
  o->ox = x;
  o->oy = y;

  o->xpf = current_xpf;
  o->xof = current_xof;
  o->ypf = current_ypf;
  o->yof = current_yof;

  o->border  = current_border;
  o->outline = current_outline;
  o->hilite  = current_hilite;

  b->col_dark = greys;
  b->col_lite = greyh;
  b->col_up   = greyf;
  b->col_down = greyf;
  b->col_sel  = greyh;
  b->col_text = mt_black();

  { MT_OBJECT *bitlabel;
    bitmap_desc bm;

    bitlabel   = label_test(label);
    if (bitlabel == NULL)
      { o->pixlab  = 0;
        if (label != NULL)
          o->label  = strcpy((char *) malloc(strlen(label)+1),label);
        else
          o->label = NULL;
      }
    else
      { o->pixlab = 1;
        o->label  = (char *) bitlabel;
        bm = (bitmap_desc) (bitlabel->wpwin);
        b->map_sel  = bm_new(bitlabel->xh+1,bitlabel->yh+1);
        b->map_down = bm_new(bitlabel->xh+1,bitlabel->yh+1);
        bm_draw_bitmap(b->map_sel,0,0,bm,MODE_SRC);
        bm_draw_bitmap(b->map_down,0,0,bm,MODE_SRC);
        bm_map_pixels(b->map_sel,b->col_up,b->col_sel);
        bm_map_pixels(b->map_down,b->col_up,b->col_down);
      }
  }

  o->free_routine  = NULL;
  if (class == RADIO)
    o->draw_routine  = draw_radio;
  else
    o->draw_routine  = draw_button;
  o->event_routine = NULL;

  b->class  = class;
  b->status = 0;
  b->menu   = NULL;

  if (class == RADIO) set_radio_geom((radio *) b);

  set_minimums(o);
  if ((o->xh - o->xl) + 1 < o->minx)
    { fprintf(stderr,"Warning: ");
      fprintf(stderr,"Button x-dim is < min needed for label & border\n");
      o->xh = (o->xl + o->minx) - 1;
    }
  if ((o->yh - o->yl) + 1 < o->miny)
    { fprintf(stderr,"Warning: ");
      fprintf(stderr,"Button y-dim is < min needed for label & border\n");
      o->yh = (o->yl + o->miny) - 1;
    }

  if (o->vis) o->draw_routine(o,0);

  return (o);
}

int mt_get_button(MT_OBJECT *o)
{ button *b;

  b = (button *) o;
  if (o->type != BUTTON)
    { fprintf(stderr,"Object is not a button (mt_get_button)\n");
      exit (1);
    }
  return (b->status);
}

void mt_set_button(MT_OBJECT *o, int is_down)
{ button *b;

  b = (button *) o;
  if (o->type != BUTTON)
    { fprintf(stderr,"Object is not a button (mt_set_button)\n");
      exit (1);
    }
  b->status = is_down;
}

void mt_button_colors(MT_OBJECT *o, long   up, long down, long dark,
                                    long lite, long  sel, long text)
{ button *b;
  long upc, slc, dnc;

  b = (button *) o;
  if (o->type != BUTTON)
    { fprintf(stderr,"Object is not a button (mt_button_colors)\n");
      exit (1);
    }

  upc = b->col_up;
  slc = b->col_sel;
  dnc = b->col_down;

  if (dark >= 0)
    b->col_dark = dark;
  if (lite >= 0)
    b->col_lite = lite;
  if (up >= 0)
    b->col_up   = up;
  if (down >= 0)
    b->col_down = down;
  if (sel >= 0)
    b->col_sel  = sel;
  if (text >= 0)
    b->col_text = text;

  if (o->pixlab)
    { bm_map_pixels((bitmap_desc) (((MT_OBJECT *) (o->label))->wpwin),
                    upc,b->col_up);
      bm_map_pixels(b->map_sel,slc,b->col_sel);
      bm_map_pixels(b->map_down,dnc,b->col_down);
    }
}

/* Text widget */

static void trim_text(textinput *t)
{ MT_OBJECT *o;
  int wid, hgt, base;
  int max, len;

  o = (MT_OBJECT *) t;
  max = o->xh - (o->xl + 2*(o->border + o->outline) + 2);
  if (max < 0) max = 0;
  len = strlen(t->text);
  while (1)
    { font_string_size(current_font,t->text,&wid,&hgt,&base);
      if (wid <= max) break;
      t->text[--len] = '\0';
    }
  if (t->esel > len) t->esel = len;
  if (t->bsel > len) t->bsel = len;
}

static int pick_text(textinput *t, int x, int *p)
{ int lft, rgt, mid, a;
  int wid, hgt, base;
  MT_OBJECT *o;

  o = (MT_OBJECT *) t;
  x -= (o->xl + o->border + o->outline + 2); 
  lft = strlen(t->text);
  rgt = 0;
  while (rgt < lft)
    { mid = (lft+rgt)/2;
      a = t->text[mid];
      t->text[mid] = '\0';
      font_string_size(current_font,t->text,&wid,&hgt,&base);
      t->text[mid] = a;
      if (wid >= x)
        lft = mid;
      else
        rgt = mid+1;
    }
  if (lft > 0)
    { a = t->text[lft];
      t->text[lft] = '\0';
      font_string_size(current_font,t->text,&wid,&hgt,&base);
      t->text[lft] = a;
      mid = wid-x;
      *p = wid;
      a = t->text[lft-1];
      t->text[lft-1] = '\0';
      font_string_size(current_font,t->text,&wid,&hgt,&base);
      t->text[lft-1] = a;
      if (x-wid < mid)
        { lft -= 1;
          *p = wid;
        }
    }
  else
    *p = 0;
  return (lft);
}

static void draw_text(MT_OBJECT *o, int high)
{ textinput *t;
  window_desc win;
  register int xl, xh, yl, yh;
  int  w, wid, wda, hgt, base;

  t = (textinput *) o;
  draw_frame(o,0);
  trim_text(t);
  w   = o->border + o->outline;
  xl  = o->xl + w;
  xh  = o->xh - w;
  yl  = o->yl + w;
  yh  = o->yh - w;
  win = o->wpwin;
  w = t->text[t->bsel];
  t->text[t->bsel] = '\0';
  font_string_size(current_font,t->text,&wid,&hgt,&base);
  t->text[t->bsel] = w;
  if (t->bsel != t->esel)
    { w = t->text[t->esel];
      t->text[t->esel] = '\0';
      font_string_size(current_font,t->text,&wda,&hgt,&base);
      t->text[t->esel] = w;
      mt_set_color(((frame *) o)->col_lsel);
      wp_win_fill(win,xl+wid+1,yl,xl+wda+2,yl+hgt+base,MODE_SRC);
    }
  mt_set_color(((frame *) o)->col_dsel);
  win_draw_string(win,xl+2,yl+hgt,t->text,current_font,MODE_SRC);
  if (o->hilite && high < 0)
    { if (t->drawn)
        { mt_set_color(((frame *) o)->col_up);
          win_draw_line(win,xl+wid+1,yl+1,xl+wid+1,yl+hgt,MODE_SRC);
          t->drawn = 0;
        }
    }
}

static void toggle_cursor(MT_OBJECT *o)
{ textinput *t;
  window_desc win;
  register int xl, xh, yl, yh;
  int  w, wid, hgt, base;

  t = (textinput *) o;
  w   = o->border + o->outline;
  xl  = o->xl + w;
  xh  = o->xh - w;
  yl  = o->yl + w;
  yh  = o->yh - w;
  win = o->wpwin;
  w = t->text[t->bsel];
  t->text[t->bsel] = '\0';
  font_string_size(current_font,t->text,&wid,&hgt,&base);
  t->text[t->bsel] = w;
  if (t->drawn)
    { mt_set_color(((frame *) o)->col_up);
      win_draw_line(win,xl+wid+1,yl+1,xl+wid+1,yl+hgt,MODE_SRC);
      t->drawn = 0;
    }
  else
    { mt_set_color(((frame *) o)->col_dsel);
      win_draw_line(win,xl+wid+1,yl+1,xl+wid+1,yl+hgt,MODE_SRC);
      t->drawn = 1;
    }
}

MT_OBJECT *mt_new_textbox(int x, int y, int w, int h, int vis, char *screen)
{ textinput *t;
  frame     *f;
  MT_OBJECT *o;

  t = (textinput *) malloc(sizeof(textinput)); 
  f = (frame *) t;
  o = (MT_OBJECT *) t;

  new_frame_fill(f,x,y,w,h,vis,REGULAR,0);

  t->text    = (char *) malloc(sizeof(char)*1024);
  t->bufsize = 1024;
  *(t->text) = '\0';
  t->bsel = t->esel = 0;
  t->lastev   = clock();

  { int i;

    for (i = 0; i < 256; i++)
      t->keyact[i] = 0;
    if (screen != NULL)
      for (i = 0; screen[i] != '\0'; i++)
        t->keyact[(int) screen[i]] = 1;
  }

  f->col_up   = mt_white();
  f->col_lite = greys;
  f->col_dark = greyh;
  f->col_dsel = mt_black();
  f->col_lsel = cyan;

  o->draw_routine = draw_text;
  o->type = TEXTBOX;
  o->hilite = 1;

  current_frame = current_frame->frame;
  return (o);
}

void mt_text_colors(MT_OBJECT *o, long fore, long dark, long lite,
                                  long text, long sel)
{ frame *f;

  f = (frame *) o;
  if (o->type != TEXTBOX)
    { fprintf(stderr,"Object is not a textbox (mt_text_colors)\n");
      exit (1);
    }

  if (dark >= 0)
    f->col_dark = dark;
  if (lite >= 0)
    f->col_lite = lite;
  if (fore >= 0)
    f->col_up = fore;
  if (text >= 0)
    f->col_dsel = text;
  if (sel >= 0)
    f->col_lsel = sel;
}

MT_OBJECT *mt_get_focus(MT_OBJECT *o)
{ MT_OBJECT *p;

  p = o;
  while (p->frame != no_object)
    p = p->frame;
  p = ((window *) p)->focus;
  if (p == no_object)
    return (NULL);
  else
    return (p);
}

void mt_set_focus(MT_OBJECT *o)
{ textinput *t;
  window    *w;
  MT_OBJECT *p;

  t = (textinput *) o;
  if (o->type != TEXTBOX)
    { fprintf(stderr,"Object is not a textbox (mt_get_text)\n");
      exit (1);
    }
  p = o;
  while (p->frame != no_object)
    p = p->frame;
  w = (window *) p;
  if (o != w->focus)
    { if (w->focus != no_object)
        { t = (textinput *) (w->focus);
          if (t->drawn)
            toggle_cursor(w->focus);
          else if (t->bsel != t->esel)
            { t->esel = t->bsel;
              draw_text(w->focus,0);
            }
        }
      w->focus = o;
    }
}

char *mt_get_text(MT_OBJECT *o, int *bsel, int *esel)
{ textinput *t;

  t = (textinput *) o;
  if (o->type != TEXTBOX)
    { fprintf(stderr,"Object is not a textbox (mt_get_text)\n");
      exit (1);
    }
  *bsel = t->bsel;
  *esel = t->esel;
  return (t->text);
}

void mt_set_text(MT_OBJECT *o, char *text, int bsel, int esel)
{ int        len;
  textinput *t;

  t = (textinput *) o;
  if (o->type != TEXTBOX)
    { fprintf(stderr,"Object is not a textbox (mt_set_text)\n");
      exit (1);
    }

  if (text != NULL)
    { len = strlen(text);
      if (len >= t->bufsize)
        { t->bufsize = 1.2*len + 100;
          t->text    = realloc(t->text,sizeof(char)*t->bufsize);
        }
      strcpy(t->text,text);
    }
  else
    len = strlen(t->text);
  if (bsel >= 0)
    t->bsel = bsel;
  if (esel >= 0)
    t->esel = esel;

  if (t->bsel > len || t->esel > len || t->bsel > t->esel)
    { fprintf(stderr,"Text selection [%d,%d] is not valid (text len is %d)\n",
                     t->bsel,t->esel,len);
      exit (1);
    }

  trim_text(t);
}

/* SCROLL object routines */

static void scroll_p2v(scrollbar *s)
{ int pw, dn;
  int lw, hg;
  MT_OBJECT *o;
 
  o = (MT_OBJECT *) s;
  if (s->vertical)
    { lw = o->yl; hg = o->yh; }
  else
    { lw = o->xl; hg = o->xh; }
  pw = (hg + 1) - (lw + 4*o->border + 2*o->outline + s->pxwid);
  dn = s->maxim - s->minim;
  if (pw == 0)
    s->value = s->minim;
  else
    s->value = s->minim
             + (((double) dn) / pw) * (s->pxpos - lw - o->border - o->outline);
}

static void scroll_v2p(scrollbar *s)
{ int pw, dn;
  int lw, hg;
  MT_OBJECT *o;
 
  o = (MT_OBJECT *) s;
  if (s->vertical)
    { lw = o->yl; hg = o->yh; }
  else
    { lw = o->xl; hg = o->xh; }
  pw = (hg + 1) - (lw + 4*o->border + 2*o->outline);
  dn = (s->maxim - s->minim) + s->width;
  if (dn == 0)
    { s->pxwid = pw;
      s->pxpos = lw + o->border + o->outline;
    }
  else
    { s->pxwid = (((double) pw) / dn) * s->width;
      pw -= s->pxwid;
      dn -= s->width;
      if (dn == 0)
        s->pxpos = lw + o->border + o->outline;
      else
        s->pxpos = lw + o->border + o->outline
                 + (((double) pw) / dn) * (s->value - s->minim);
    }
}

static void draw_scrollbar(MT_OBJECT *o, int high)
{ window_desc win;
  scrollbar *s;
  int out, cut, uni;
  int wid;
  register int xl, xh, yl, yh;
  register int  i,  w;
  int    xs, xp;
  int    ys, yp;

  s = (scrollbar *) o;

  xl  = o->xl;
  xh  = o->xh;
  yl  = o->yl;
  yh  = o->yh;
  w   = o->border;
  win = o->wpwin;
  out = o->outline;
  cut = 1-out;

  if (s->vertical)
    { yp = s->pxpos;
      ys = yp + s->pxwid + 2*w - 1;
      xp = xl+(w+out);
      xs = xh-(w+out);
    }
  else
    { xp = s->pxpos;
      xs = xp + s->pxwid + 2*w - 1;
      yp  = yl+(w+out);
      ys  = yh-(w+out);
    }

  if (high == 0)
  
    { uni = (s->col_lite == s->col_dark);

      if (uni)
        { mt_set_color(s->col_lite);
          wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
          wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
          wp_win_fill(win,xl+out,yl+out,xl+out+w,yh+cut-w,MODE_SRC);
          wp_win_fill(win,xl+w+out,yl+out,xh+cut-w,yl+out+w,MODE_SRC);
        }
      else
        { mt_set_color(s->col_lite);
          wp_win_fill(win,xl+out,yh+cut,xh+cut,yh+cut-w,MODE_SRC);
          wp_win_fill(win,xh+cut,yl+out,xh+cut-w,yh+cut,MODE_SRC);
          mt_set_color(s->col_dark);
          wid = w - cut;
          for (i = out; i <= wid; i++)
            { win_draw_line(win,xl+out,yl+i,xh-1-i,yl+i,MODE_SRC);
              win_draw_line(win,xl+i,yl+out,xl+i,yh-1-i,MODE_SRC);
            }
        }
   
      mt_set_color(s->col_well);
      wp_win_fill(win,xl+out+w,yl+out+w,xh+cut-w,yh+cut-w,MODE_SRC);

      if (uni)
        { mt_set_color(s->col_lite);
          wp_win_fill(win,xp,ys+1,xs+1,ys+1-w,MODE_SRC);
          wp_win_fill(win,xs+1,yp,xs+1-w,ys+1,MODE_SRC);
          wp_win_fill(win,xp,yp,xp+w,ys+1-w,MODE_SRC);
          wp_win_fill(win,xp+w,yp,xs+1-w,yp+w,MODE_SRC);
        }
      else
        { mt_set_color(s->col_dark);
          wp_win_fill(win,xp,ys+1,xs+1,ys+1-w,MODE_SRC);
          wp_win_fill(win,xs+1,yp,xs+1-w,ys+1,MODE_SRC);
          mt_set_color(s->col_lite);
          for (i = 0; i < w; i++)
            { win_draw_line(win,xp,yp+i,xs-1-i,yp+i,MODE_SRC);
              win_draw_line(win,xp+i,yp,xp+i,ys-1-i,MODE_SRC);
            }
        }

      mt_set_color(s->col_tab);
      wp_win_fill(win,xp+w,yp+w,xs+1-w,ys+1-w,MODE_SRC);

      if (out)
        { mt_set_color(mt_black());
          win_draw_line(win,xl,yl,xl,yh,MODE_SRC);
          win_draw_line(win,xh,yl,xh,yh,MODE_SRC);
          win_draw_line(win,xl,yl,xh,yl,MODE_SRC);
          win_draw_line(win,xl,yh,xh,yh,MODE_SRC);
        }
    }

  else
    { if (high > 0)
        mt_set_color(s->col_sel);
      else
        mt_set_color(s->col_well);
      if (s->vertical)
        { wp_win_fill(win,xp,yl+out+w,xs+1,yp,MODE_SRC);
          wp_win_fill(win,xp,ys+1,xs+1,yh+cut-w,MODE_SRC);
        }
      else
        { wp_win_fill(win,xl+out+w,yp,xp,ys+1,MODE_SRC);
          wp_win_fill(win,xs+1,yp,xh+cut-w,ys+1,MODE_SRC);
        }
    }
}

static void move_scrollbar(MT_OBJECT *o, int delta)
{ window_desc win;
  scrollbar *s;
  int out, uni;
  register int xl, xh, yl, yh;
  register int  i,  w;
  int    xs, xp;
  int    ys, yp;

  s = (scrollbar *) o;

  xl  = o->xl;
  xh  = o->xh;
  yl  = o->yl;
  yh  = o->yh;
  w   = o->border;
  win = o->wpwin;
  out = o->outline;

  if (s->vertical)
    { yp = s->pxpos;
      ys = yp + s->pxwid + 2*w - 1;
      xp = xl+(w+out);
      xs = xh-(w+out);
    }
  else
    { xp = s->pxpos;
      xs = xp + s->pxwid + 2*w - 1;
      yp  = yl+(w+out);
      ys  = yh-(w+out);
    }

  uni = (s->col_lite == s->col_dark);

  mt_set_color(s->col_sel);
  wp_win_fill(win,xp,yp,xs+1,ys+1,MODE_SRC);

  if (s->vertical)
    { yp += delta; ys += delta; }
  else
    { xp += delta; xs += delta; }

  if (uni)
    { mt_set_color(s->col_lite);
      wp_win_fill(win,xp,ys+1,xs+1,ys+1-w,MODE_SRC);
      wp_win_fill(win,xs+1,yp,xs+1-w,ys+1,MODE_SRC);
      wp_win_fill(win,xp,yp,xp+w,ys+1-w,MODE_SRC);
      wp_win_fill(win,xp+w,yp,xs+1-w,yp+w,MODE_SRC);
    }
  else
    { mt_set_color(s->col_dark);
      wp_win_fill(win,xp,ys+1,xs+1,ys+1-w,MODE_SRC);
      wp_win_fill(win,xs+1,yp,xs+1-w,ys+1,MODE_SRC);
      mt_set_color(s->col_lite);
      for (i = 0; i < w; i++)
        { win_draw_line(win,xp,yp+i,xs-1-i,yp+i,MODE_SRC);
          win_draw_line(win,xp+i,yp,xp+i,ys-1-i,MODE_SRC);
        }
    }
 
  mt_set_color(s->col_tab);
  wp_win_fill(win,xp+w,yp+w,xs+1-w,ys+1-w,MODE_SRC);
}

static void erase_scrollbar(MT_OBJECT *o)
{ register window_desc win;

  win = o->wpwin;
  mt_set_color(((frame *) o->frame)->col_up);
  wp_win_fill(win,o->xl,o->yl,o->xh+1,o->yh+1,MODE_SRC);
}

MT_OBJECT *mt_new_scrollbar(int x, int y, int w, int h, int vert, int vis,
                            int val, int wid, int min, int max)
{ scrollbar *s;
  MT_OBJECT *fr, *o;
  MT_OBJECT *p;

  fr = current_frame;
  if (fr == no_object)
    { fprintf(stderr,"Fatal: No current frame (mt_new_button)\n");
      exit (1);
    }

  s = (scrollbar *) malloc(sizeof(scrollbar));
  o = (MT_OBJECT *) s;

  p = (MT_OBJECT *) ( & (((frame *) fr)->afor));
  o->back = p->back;
  o->forw = p;
  o->back->forw = o->forw->back = o;

  o->type  = SCROLLBAR;
  o->vis   = vis;
  o->wpwin = fr->wpwin;
  o->frame = fr;
  o->level = fr->level + 1;

  o->xl = fr->xl + x;
  o->yl = fr->yl + y;
  o->xh = o->xl + (w-1);
  o->yh = o->yl + (h-1);
  o->ow = w;
  o->oh = h;
  o->ox = x;
  o->oy = y;

  o->xpf = current_xpf;
  o->xof = current_xof;
  o->ypf = current_ypf;
  o->yof = current_yof;

  o->border  = current_border;
  o->outline = current_outline;
  o->hilite  = current_hilite;
  o->pixlab  = 0;
  o->label   = NULL;

  o->free_routine  = NULL;
  o->draw_routine  = draw_scrollbar;
  o->event_routine = NULL;

  s->vertical = vert;
  s->value = val;
  s->width = wid;
  s->minim = min;
  s->maxim = max;
  if (min > val || wid < 0 || val > max)
    { fprintf(stderr,"Scroll bar being set to impossible values:\n");
      fprintf(stderr,"  min(%d) <= val(%d) <= max(%d), wid(%d) >= 0 ???\n",
                     min,val,max,wid);
      exit (1);
    }
  scroll_v2p(s);

  s->col_dark = greys;
  s->col_lite = greyh;
  s->col_tab  = greyf;
  s->col_well = greyf;
  s->col_sel  = greyh;

  set_minimums(o);
  if ((o->xh - o->xl) + 1 < o->minx)
    { fprintf(stderr,"Warning: ");
      fprintf(stderr,"Scrollbar x-dim is < min needed for label & border\n");
      o->xh = (o->xl + o->minx) - 1;
    }
  if ((o->yh - o->yl) + 1 < o->miny)
    { fprintf(stderr,"Warning: ");
      fprintf(stderr,"Scrollbar y-dim is < min needed for label & border\n");
      o->yh = (o->yl + o->miny) - 1;
    }

  if (o->vis) o->draw_routine(o,0);

  return (o);
}

int mt_get_scrollbar(MT_OBJECT *o, int *wid, int *min, int *max)
{ scrollbar *s;

  s = (scrollbar *) o;
  if (o->type != SCROLLBAR)
    { fprintf(stderr,"Object is not a scrollbar (mt_get_scrollbar)\n");
      exit (1);
    }
  *wid = s->width;
  *min = s->minim;
  *max = s->maxim;
  return (s->value);
}

void mt_set_scrollbar(MT_OBJECT *o, int val, int wid, int min, int max)
{ scrollbar *s;

  s = (scrollbar *) o;
  if (o->type != SCROLLBAR)
    { fprintf(stderr,"Object is not a scrollbar (mt_set_scrollbar)\n");
      exit (1);
    }
  if (min > val || wid < 0 || val > max)
    { fprintf(stderr,"Scroll bar begin set to impossible values:\n");
      fprintf(stderr,"  min(%d) <= val(%d) <= max(%d), wid(%d) >= 0 ???\n",
                     min,val,max,wid);
    }
  s->value = val;
  s->width = wid;
  s->minim = min;
  s->maxim = max;
  scroll_v2p(s);
}

void mt_scrollbar_colors(MT_OBJECT *o, long  tab, long well, long dark,
                                       long lite, long  sel)
{ scrollbar *s;

  s = (scrollbar *) o;
  if (o->type != SCROLLBAR)
    { fprintf(stderr,"Object is not a scrollbar (mt_scrollbar_colors)\n");
      exit (1);
    }

  if (dark >= 0)
    s->col_dark = dark;
  if (lite >= 0)
    s->col_lite = lite;
  if (tab >= 0)
    s->col_tab  = tab;
  if (well >= 0)
    s->col_well = well;
  if (sel >= 0)
    s->col_sel  = sel;
}

/* Bitmap object routines */

void draw_bitmap(MT_OBJECT *o, int high)
{ return; }

void free_bitmap(MT_OBJECT *o)
{ bm_free((bitmap_desc) (o->wpwin)); }

MT_OBJECT *mt_new_bitmap(int w, int h)
{ MT_OBJECT *o;

  o = (MT_OBJECT *) malloc(sizeof(MT_OBJECT));
  o->forw = o;
  o->back = o;

  o->type  = BITMAP;
  o->vis   = 0;
  o->wpwin = (window_desc) bm_new(w,h);
  o->frame = NULL;
  o->level = 0;

  o->xl = 0;
  o->yl = 0;
  o->xh = o->xl + (w-1);
  o->yh = o->yl + (h-1);
  o->ow = w;
  o->oh = h;
  o->ox = 0;
  o->oy = 0;

  o->xpf = current_xpf;
  o->xof = current_xof;
  o->ypf = current_ypf;
  o->yof = current_yof;

  o->border  = current_border;
  o->outline = current_outline;
  o->hilite  = current_hilite;
  o->pixlab  = 0;
  o->label   = NULL;

  o->free_routine  = free_bitmap;
  o->draw_routine  = draw_bitmap;
  o->event_routine = NULL;

  return (o);
}

void mt_menu_pair(MT_OBJECT *menubutton, MT_OBJECT *menuframe)
{ button *b;
  frame  *f;
  MT_OBJECT *o, *p, *fr;

  o = menuframe;
  f = (frame *) o;

  if (menubutton->type != BUTTON)
    { fprintf(stderr,"1st arg is not a button (mt_menu_pair)\n");
      exit (1);
    }
  if (menuframe->type != FRAME)
    { fprintf(stderr,"2nd arg is not a frame (mt_menu_pair)\n");
      exit (1);
    }
  if (f->ctlbut != NULL)
    { fprintf(stderr,"frame is already assigned to another button");
      fprintf(stderr," (mt_menu_pair)\n");
      exit (1);
    }

  b = (button *) menubutton;
  b->menu = menuframe;

  f->ctlbut = menuframe;
  o->forw->back = o->back;
  o->back->forw = o->forw;

  fr = o->frame;
  while (fr->frame != no_object)
    fr = fr->frame;
  p = (MT_OBJECT *) ( & (((frame *) fr)->afor));
  o->back = p->back;
  o->forw = p;
  o->back->forw = o->forw->back = o;
  o->frame = fr;
}

void mt_startup()
{ winpack_init();
  greyf = mt_get_color(0xB2,0xB2,0xB2);
  greys = mt_get_color(0x65,0x65,0x65);
  greyh = mt_get_color(0xE1,0xE1,0xE1);
  cyan  = mt_get_color(0x00,0xFF,0xFF);
  init_attributes();
}

void mt_shutdown()
{ winpack_shutdown(); }

static MT_OBJECT *find(MT_OBJECT *bw, int x, int y)
{ register MT_OBJECT *o, *p;
  register radio  *r;
  register int off;

  p = (MT_OBJECT *) ( & (((frame *) bw)->afor));
  for (o = p->back; o != p; o = o->back)
    if (o->vis && o->xl <= x && x <= o->xh && o->yl <= y && y <= o->yh)
      if (o->type == BUTTON)
        if (((button *) o)->class == RADIO)
          { r = (radio *) o; 
            if (y <= r->ym)
              off = r->ym - y;
            else
              off = y - r->ym;
            off = r->radius - off;
            if (r->xm - off <= x && x <= r->xm + off)
              return (o);
          }
        else
          return (o);
      else if (o->type == SCROLLBAR || o->type == TEXTBOX)
        return (o);
      else
        return (find(o,x,y));
  return (bw);
}

static void scale_adjust(MT_OBJECT *bin, int w, int h)
{ MT_OBJECT *o, *p;
  int delta;

  p = (MT_OBJECT *) ( & (((frame *) bin)->afor));

  for (o = p->forw; o != p; o = o->forw)
    { delta = w - bin->ow;
      o->xl = bin->xl + o->ox + o->xpf*delta;
      o->xh = bin->xl + o->ox + o->ow + (o->xpf + o->xof)*delta - 1;

      delta = h - bin->oh;
      o->yl = bin->yl + o->oy + o->ypf*delta;
      o->yh = bin->yl + o->oy + o->oh + (o->ypf + o->yof)*delta - 1;

      if (o->type == FRAME)
        scale_adjust(o,o->xh - o->xl + 1,o->yh - o->yl + 1);

      else if (o->type == BUTTON && ((button *) o)->class == RADIO)
        set_radio_geom((radio *) o);

      else if (o->type == SCROLLBAR)
        scroll_v2p((scrollbar *) o);
    }
}

static void hilite_update(MT_OBJECT *o, MT_OBJECT *p)
{ if (o == p) return;

  if (p->level >= o->level && p->hilite)
    p->draw_routine(p,-1);

  if (o->level > p->level)
    hilite_update (o->frame,p);
  else if (o->level < p->level)
    hilite_update(o,p->frame);
  else
    hilite_update (o->frame,p->frame);

  if (o->level >= p->level && o->hilite)
    o->draw_routine(o,1);
}

static void millisleep(int millisecs)
{ clock_t beg;

  beg = clock();
  while ( (1000.*(clock() - beg))/CLOCKS_PER_SEC < millisecs)
    continue;
}

static void bp_handler(input_event *ev, window_desc win, void *data)
{ int code, x, y, value;
  long when;
  int w, h, click;
  window    *root;
  MT_OBJECT *focus;
  MT_OBJECT *o;
  button    *bu;
  scrollbar *sb;
  static window *curwin;

  if (dialog_stack != NULL)
    { if (dialog_stack->wpwin != win) return; }

  code    = ev->code;
  x       = ev->x;
  y       = ev->y;
  value   = ev->value;
  when    = (long) ev->when;
  root    = (window *) data;

  ev_x         = x;
  ev_y         = y;
  ev_value     = value;
  ev_modifiers = ev->updown_mask;

  if (code != win_resize)
    { if (root->started && when < root->recentime) return;
      root->recentime = when;
      root->started = 1;
    }
  if (code == win_exit)
    curwin = NULL;
  else if (code == win_enter)
    curwin = root;

  if (root->state == MENUDOWN)
    focus = (MT_OBJECT *) (root->menuf);
  else
    focus = (MT_OBJECT *) data;
  switch (root->state)
  { case IDLE:
    case MENUDOWN:
      if (code == win_resize)
        { MT_OBJECT *winobj;

          winobj = (MT_OBJECT *) root;
          w = win_width(win);
          h = win_height(win);
          if (root->maxwide < w || w < root->minwide ||
              root->maxhigh < h || h < root->minhigh)
            { if (w < root->minwide)
                w = root->minwide;
              else if (root->maxwide < w)
                w = root->maxwide;
              if (h < root->minhigh)
                h = root->minhigh;
              else if (root->maxhigh < h)
                h = root->maxhigh;
              win_set_size(win,w,h);
            }
          scale_adjust(winobj,w,h);
          win_set_clip(win,0,0,w,h);
          winobj->xh = w-1;
          winobj->yh = h-1;
          if (root->state == MENUDOWN)
            root->menuf->vis = 0;
          winobj->draw_routine(winobj,0);
          if (root->state == MENUDOWN)
            { wp_win_saveblk(root->backing,0,0,
                             root->menuf->ow,root->menuf->oh,
                             win,root->menuf->xl,root->menuf->yl,MODE_SRC);
              root->menuf->vis = 1;
              root->menuf->draw_routine(root->menuf,0);
            }
          root->object = no_object;
        }
      if (code == win_exit)
        o = no_object;
      else if (root->state == MENUDOWN)
        { o = root->menub;
          if (o->xl > x || x > o->xh || o->yl > y || y > o->yh)
            o = find(focus,x,y);
        }
      else
        o = find(focus,x,y);
      hilite_update(o,root->object);
      root->object = o;
      if (code == loc_but_up)
        { root->but_dwn = -1;
          if (o->type == FRAME)
            { if (o->event_routine != NULL)
                { if ((((frame *) o)->event_chord & RELEASE_EVENT) != 0)
                    { ev_type = RELEASE_EVENT;
                      o->event_routine(o,o->event_data);
                    }
                }
            }
          wp_user_event(root->recentime,win);
        } 
      else if (code == loc_but_dn)
        { root->but_dwn  = value;
          root->dwn_time = when;
          if (o->type == FRAME)
            { if (o->event_routine != NULL)
                { if ((((frame *) o)->event_chord & PRESS_EVENT) != 0)
                    { ev_type = PRESS_EVENT;
                      o->event_routine(o,o->event_data);
                    }
                }
            }
          else if (o->type == BUTTON)
            { bu = (button *) o;
              bu->status = 1 - bu->status;
              if (bu->menu != NULL)
                { if (root->state != MENUDOWN)
                    { root->menuf = bu->menu;
                      root->menub = o;
                      root->backing = bm_new(root->menuf->ow,root->menuf->oh);
                      wp_win_saveblk(root->backing,0,0,
                                     root->menuf->ow,root->menuf->oh,
                                     win,root->menuf->xl,root->menuf->yl,
                                     MODE_SRC);
                      bu->menu->vis = 1;
                      bu->menu->draw_routine(bu->menu,0);
                      root->state = MENUPRESS;
                    }
                  else
                    root->state = MENUTWICE;
                }
              else
                switch (bu->class)
                { case CLICK:
                    root->state = PRESSED;
                    break;
                  case PRESS:
                    root->state = PRESSED;
                    if (o->event_routine != NULL)
                      { ev_type = PRESS_EVENT;
                        o->event_routine(o,o->event_data);
                      }
                    break;
                  case TOGGLE:
                  case RADIO:
                    if (o->event_routine != NULL)
                      { ev_type = PRESS_EVENT;
                        o->event_routine(o,o->event_data);
                      }
                    break;
                }
              o->draw_routine(o,0);
            }
          else if (o->type == SCROLLBAR)
            { int xp, xs;
              int lw, hg, ps;
              int move;

              sb = (scrollbar *) o;
              if (sb->vertical)
                { lw = o->yl + o->border + o->outline;
                  hg = o->yh - (o->border + o->outline);
                  ps = y;
                }
              else
                { lw = o->xl + o->border + o->outline;
                  hg = o->xh - (o->border + o->outline);
                  ps = x;
                }
              xp = sb->pxpos;
              xs = xp + sb->pxwid + 2*o->border - 1;
              if (ps < xp)
                { if (move = (ps >= lw))
                    sb->pxpos = ps;
                }
              else if (ps > xs)
                { if (move = (ps <= hg))
                    sb->pxpos += ps-xs;
                }
              else
                { root->state = SLIDER;
                  root->slide_delta = sb->pxpos - ps;
                }
              if (root->state != SLIDER && move)
                { scroll_p2v(sb);
                  o->draw_routine(o,0);
                  o->draw_routine(o,1);
                  if (o->event_routine != NULL)
                    { ev_type = PRESS_EVENT;
                      o->event_routine(o,o->event_data);
                    }
                }
            }
          else if (o->type == TEXTBOX && o->hilite)
            { textinput *t;
              int xs, del;
              int hgt, wid, base;

              if (o != root->focus)
                { if (root->focus != no_object)
                    { t = (textinput *) (root->focus);
                      if (t->drawn)
                        toggle_cursor(root->focus);
                      else if (t->bsel != t->esel)
                        { t->esel = t->bsel;
                          draw_text(root->focus,0);
                        }
                    }
                  root->focus = o;
                }

              t = (textinput *) o;
              if (t->drawn)
                toggle_cursor(o);
              else if (t->bsel != t->esel)
                { t->esel = t->bsel;
                  draw_text(o,0);
                }
              xs = o->xl + o->border + o->outline + 1;
              t->bsel = t->esel = pick_text(t,x,&del);
              font_string_size(current_font,t->text,&wid,&hgt,&base);
              root->xdown = root->xlast = xs + del;
              toggle_cursor(o);
              root->state = TEXTSEL;
            }
        }
      else if (code == key_dn)
        { o = root->focus;
          if (o != no_object)
            { textinput *t;
              char *u, *v;
              int   a,  b;
              int wid1, wid2, hgt, base;
              int xl, yl, del, rez;

              t = (textinput *) o;

              ev_type = KEY_EVENT;
              rez = 1;
              if (t->keyact[value] && o->event_routine != NULL)
                rez = o->event_routine(o,o->event_data);

              if (rez)
                { font_string_size(current_font,t->text,&wid2,&hgt,&base);
              
                  del = 0;
                  if (t->bsel != t->esel || (value == 0x08 && t->bsel > 0))
                    { if (value == 0x08 && t->bsel == t->esel)
                        t->bsel -= 1;
                      u = t->text + t->bsel;
                      v = t->text + t->esel;
                      while (*v != '\0')
                        *u++ = *v++;
                      *u = '\0';
                      del -= (t->esel - t->bsel);
                      t->esel = t->bsel;
                    }

                  u = t->text + t->bsel;
                  a = *u;
                  *u = '\0';
                  font_string_size(current_font,t->text,&wid1,&hgt,&base);

                  if (value != 0x08)
                    { *u++ = value;
                      while (a != '\0')
                        { b = *u;
                          *u++ = a;
                          a = b;
                        }
                      *u++ = a;
                      *u = '\0';
                      del += 1;
                    }
                  else
                    *u = a;

                  if (t->drawn) toggle_cursor(o);
                  trim_text(t);
                  if (del > 0)
                    font_string_size(current_font,t->text,&wid2,&hgt,&base);
                  xl = o->xl + o->border + o->outline + 2;
                  yl = o->yl + o->border + o->outline;
                  mt_set_color(((frame *) o)->col_up);
                  wp_win_fill(win,xl+wid1-1,yl,xl+wid2,yl+hgt+base,MODE_SRC);
                  mt_set_color(((frame *) o)->col_dsel);
                  win_draw_string(win,xl+wid1,yl+hgt,
                                      t->text + t->bsel,current_font,MODE_SRC);
                  if (value != 0x08) t->bsel += 1;
                  t->esel  = t->bsel;
                }
            }
        }
      if (root->focus != no_object && root->but_dwn < 0)
        { textinput *t;
          t = (textinput *) (o = root->focus);
          if (curwin != root)
            { if (t->drawn) toggle_cursor(o); }
          else if (t->bsel == t->esel)
            { if ( (1000.*(clock() - t->lastev))/CLOCKS_PER_SEC > 500)
                { toggle_cursor(o); 
                  t->lastev = clock();
                }
              millisleep(25);
              wp_user_event(root->recentime,win);
            }
        }
      break;

    case TEXTSEL:
      { textinput *t;
        int yl, xl;
        int wid, hgt, base;

        o  = (MT_OBJECT *) (root->object);
        t  = (textinput *) o;

        yl = o->yl + o->border + o->outline;
        xl = o->xl + o->border + o->outline;

        if (x < xl)
          x = o->xl + o->border + o->outline;
        if (x > o->xh - (o->border + o->outline))
          x = o->xh - (o->border + o->outline);
        if (t->drawn) toggle_cursor(o);

        font_string_size(current_font,t->text,&wid,&hgt,&base);
        base += hgt;
        if (code == loc_but_up)
          { t->esel = pick_text(t,x,&wid);
            x = xl + 1 + wid;
          }
        if (root->xdown < root->xlast)
          if (x < root->xdown)
            { mt_set_color(((frame *) o)->col_up);
              wp_win_fill(win,root->xdown+1,yl,
                              root->xlast+1,yl+base,MODE_SRC);
              mt_set_color(((frame *) o)->col_lsel);
              wp_win_fill(win,x,yl,root->xdown+1,yl+base,MODE_SRC);
            }
          else if (x < root->xlast)
            { mt_set_color(((frame *) o)->col_up);
              wp_win_fill(win,x+1,yl,root->xlast+1,yl+base,MODE_SRC);
            }
          else
            { mt_set_color(((frame *) o)->col_lsel);
              wp_win_fill(win,root->xlast+1,yl,x+1,yl+base,MODE_SRC);
            }
        else
          if (x > root->xdown)
            { mt_set_color(((frame *) o)->col_up);
              wp_win_fill(win,root->xlast,yl,
                              root->xdown,yl+base,MODE_SRC);
              mt_set_color(((frame *) o)->col_lsel);
              wp_win_fill(win,root->xdown,yl,x+1,yl+base,MODE_SRC);
            }
          else if (x > root->xlast)
            { mt_set_color(((frame *) o)->col_up);
              wp_win_fill(win,root->xlast,yl,x,yl+base,MODE_SRC);
            }
          else
            { mt_set_color(((frame *) o)->col_lsel);
              wp_win_fill(win,x,yl,root->xlast+1,yl+base,MODE_SRC);
            }
        mt_set_color(((frame *) o)->col_dsel);
        win_draw_string(win,xl+2,yl+hgt,t->text,current_font,MODE_SRC);

        root->xlast = x;
        if (code == loc_but_up)
          { if (t->bsel == t->esel)
              { mt_set_color(((frame *) o)->col_up);
                win_draw_line(win,xl+wid+1,yl,xl+wid+1,yl+(base-1),MODE_SRC);
              }
            else if (t->bsel > t->esel)
              { x = t->bsel;
                t->bsel = t->esel;
                t->esel = x;
              }
            root->but_dwn = -1;
            root->state   = IDLE;
            wp_user_event(root->recentime,win);
          }
        break;
      }

    case SLIDER:
      { int lw, hg, nw, ps;

        sb = (scrollbar *) root->object;
        o  = root->object;
        if (sb->vertical)
          { lw = o->yl + o->border + o->outline;
            hg = (o->yh + 1) - (3*o->border + o->outline + sb->pxwid);
            ps = y;
          }
        else
          { lw = o->xl + o->border + o->outline;
            hg = (o->xh + 1) - (3*o->border + o->outline + sb->pxwid);
            ps = x;
          }
        nw = ps+root->slide_delta;
        if (nw < lw)
          nw = lw;
        else if (nw > hg)
          nw = hg;
        if (nw != sb->pxpos)
          move_scrollbar(root->object,nw-sb->pxpos);
        sb->pxpos = nw;
        if (code == loc_but_up)
          { root->state = IDLE;
            scroll_p2v(sb);
            if (o->event_routine != NULL)
              { ev_type = RELEASE_EVENT;
                o->event_routine(o,o->event_data);
              }
            root->but_dwn = -1;
            wp_user_event(root->recentime,win);
          }
        break;
      }

    case PRESSED:
      o = find(focus,x,y);
      if (o != root->object)
        { MT_OBJECT *press;

          press = root->object;
          root->object = o;
          if (press->type == BUTTON)
            { bu = (button *) press;
              bu->status = 1 - bu->status;
              if (bu->class == PRESS)
                if (press->event_routine != NULL)
                  { ev_type = RELEASE_EVENT;
                    press->event_routine(press,press->event_data);
                  }
              press->draw_routine(press,0);
            }
          hilite_update(o,press->frame);
          root->state = IDLE;
        }
      else if (code == loc_but_up & root->but_dwn == value)
        { if (o->type == BUTTON)
            { bu = (button *) o;
              switch (bu->class)
              { case CLICK:
                  click = (when - root->dwn_time < 1000);
                  bu->status = 1 - bu->status;
                  o->draw_routine(o,0);
                  if (click)
                    if (o->event_routine != NULL)
                      { ev_type = RELEASE_EVENT;
                        o->event_routine(o,o->event_data);
                      }
                  break;
                case PRESS:
                  bu->status = 1 - bu->status;
                  o->draw_routine(o,0);
                  if (o->event_routine != NULL)
                    { ev_type = RELEASE_EVENT;
                      o->event_routine(o,o->event_data);
                    }
                  break;
              }
            }
          root->but_dwn = -1;
          root->state   = IDLE;
          wp_user_event(root->recentime,win);
        }
      break;

    case MENUPRESS:
    case MENUTWICE:
      o = find(root->menuf,x,y);
      hilite_update(o,root->object);
      root->object = o;
      if (code == loc_but_up & root->but_dwn == value)
        { bu = (button *) root->menub;
          bu->status = 1-bu->status;
          root->menub->draw_routine(root->menub,0);
          if (find(focus,x,y) == root->menub && when - root->dwn_time < 1000)
            if (root->state == MENUPRESS)
              { root->state   = MENUDOWN;
                root->but_dwn = -1; 
                wp_user_event(root->recentime,win);
                break;
              }
          root->menuf->vis = 0;
          win_draw_bitmap(win,root->menuf->xl,root->menuf->yl,
                          root->backing,MODE_SRC);
          bm_free(root->backing);
          if (o != root->menuf || root->state != MENUPRESS)
            if (o->event_routine != NULL)
              { ev_type = RELEASE_EVENT;
                o->event_routine(o,o->event_data);
              }
          root->object  = no_object;
          root->but_dwn = -1;
          root->state   = IDLE;
          wp_user_event(root->recentime,win);
        }
      break;
  }
}

void mt_event_loop()
{ winpack_main_loop(); }

long mt_medium_grey()
{ return (greyf); }

long mt_dark_grey()
{ return (greys); }

long mt_light_grey()
{ return (greyh); }

MT_OBJECT *mt_parent(MT_OBJECT *o)
{ if (o->frame == no_object)
    return (NULL);
  else
    return (o->frame);
}

MT_OBJECT *mt_children(MT_OBJECT *o)
{ frame *f;

  f = (frame *) o;
  if (o->type != FRAME)
    { fprintf(stderr,"Object is not a frame (mt_first_child)\n");
      exit (1);
    }

  if (f->afor == (MT_OBJECT *) ( & (f->afor)))
    return (NULL);
  else
    return (f->afor);
}

MT_OBJECT *mt_sibling(MT_OBJECT *o)
{ frame *f;

  if (o->frame == no_object)
    return (NULL);
  f = (frame *) (o->frame);
  if (f->afor == (MT_OBJECT *) ( & (f->afor)))
    return (NULL);
  else
    return (o->forw);
}
