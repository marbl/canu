
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

#ifndef BUTTON_H
#define BUTTON_H

static const char *rcsid_BUTTON_H = "$Id: button.h,v 1.6 2008-10-08 22:03:02 brianwalenz Exp $";

typedef void MT_OBJECT;

#define CLICK  0
#define PRESS  1
#define TOGGLE 2
#define RADIO  3

#define REGULAR 1
#define VIRTUAL 2

#define REDRAW_EVENT  0x01
#define PRESS_EVENT   0x02
#define RELEASE_EVENT 0x04
#define KEY_EVENT     0x08
#define ENTRY_EVENT   0x10
#define EXIT_EVENT    0x20

#define SHIFT_MODIFIER 0x1
#define CNTRL_MODIFIER 0x2
#define META_MODIFIER  0x4

void mt_current_xscale(double xpf, double xof);
void mt_current_yscale(double ypf, double yof);
void mt_current_border(int width);
void mt_current_outline(int is_on);
void mt_current_hilite(int is_on);

void mt_current_frame(MT_OBJECT *f);
void mt_pop_frame(void);

void mt_set_xscale(MT_OBJECT *o, double xpf, double xof);
void mt_set_yscale(MT_OBJECT *o, double ypf, double yof);
void mt_set_border(MT_OBJECT *o, int width);
void mt_set_outline(MT_OBJECT *o, int is_on);
void mt_set_hilite(MT_OBJECT *o, int is_on);

void mt_set_vis(MT_OBJECT *o, int is_on);
void mt_set_label(MT_OBJECT *o, char *label);
void mt_set_callback(MT_OBJECT *o,
                     int (*callback)(MT_OBJECT *p, long d), long data);

void mt_get_extent(MT_OBJECT *o, int *xl, int *xh, int *yl, int *yh);
int  mt_get_vis(MT_OBJECT *o);
int  mt_get_event(MT_OBJECT *o, int *x, int *y, int *value, int *mods);

void mt_draw(MT_OBJECT *o);
void mt_free(MT_OBJECT *o);

MT_OBJECT *mt_new_dialog(int x, int y, int w, int h, int vis, char *title);
MT_OBJECT *mt_new_window(int x, int y, int w, int h, int vis, char *title);
void  mt_set_window_bounds(MT_OBJECT *o, int minw, int minh,
                                             int maxw, int maxh);
void mt_set_window_size(MT_OBJECT *o, int width, int height);
MT_OBJECT  *mt_new_frame(int x, int y, int w, int h,
                         int vis, int classa, int ev_chord);
void mt_frame_colors(MT_OBJECT *f, long fore, long dark, long lite,
                                   long dsel, long lsel);

MT_OBJECT  *mt_new_textbox(int x, int y, int w, int h, int vis, char *screen);
void mt_text_colors(MT_OBJECT *t, long fore, long dark, long lite,
                                  long text, long sel);
char *mt_get_text(MT_OBJECT *o, int *bsel, int *esel);
void  mt_set_text(MT_OBJECT *o, char *text, int bsel, int esel);
void  mt_set_focus(MT_OBJECT *o);

MT_OBJECT *mt_new_button(int x, int y, int w, int h, int vis,
                         int classa, char *label);
int  mt_get_button(MT_OBJECT *b);
void mt_set_button(MT_OBJECT *b, int is_down);
void mt_button_colors(MT_OBJECT *b, long   up, long down, long dark,
                                    long lite, long  sel, long text);

MT_OBJECT *mt_new_scrollbar(int x, int y, int w, int h, int vert, int vis,
                            int val, int wid, int min, int max);
int  mt_get_scrollbar(MT_OBJECT *s, int *wid, int *min, int *max);
void mt_set_scrollbar(MT_OBJECT *s, int val, int wid, int min, int max);
void mt_scrollbar_colors(MT_OBJECT *s, long  tab, long well, long dark,
                                       long lite, long  sel);

MT_OBJECT *mt_new_bitmap(int w, int h);
char *mt_bitmap_label(MT_OBJECT *b);

void mt_menu_pair(MT_OBJECT *button, MT_OBJECT *menuframe);

void mt_startup(void);
void mt_shutdown(void);
void mt_event_loop(void);

long mt_get_color(int red, int green, int blue);
void mt_set_color(long color);
long mt_black(void);
long mt_white(void);
long mt_medium_grey(void);
long mt_dark_grey(void);
long mt_light_grey(void);

void mt_draw_line (MT_OBJECT *o, int xl, int yl, int xh, int yh,
                                 int thick, int dash);
void mt_draw_text (MT_OBJECT *o, int x, int y, char *text);
void mt_draw_pixel(MT_OBJECT *o, int x, int y);
void mt_draw_rect (MT_OBJECT *o, int xl, int yl, int w, int h);
void mt_string_size(char *text, int *wide, int *height, int *base);
void mt_begin_clip_to(MT_OBJECT *o);
void mt_end_clip_to(MT_OBJECT *o);

/* This is a temporary kludge */
void mt_draw_title(MT_OBJECT *o, int x, int y, char *text);

MT_OBJECT *mt_parent(MT_OBJECT *o);
MT_OBJECT *mt_children(MT_OBJECT *o);
MT_OBJECT *mt_sibling(MT_OBJECT *o);

#endif
