
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
/* $Id: wpack.c,v 1.6 2007-03-09 03:05:58 brianwalenz Exp $ */

#undef DEBUG

/* Minimal Retargetable Window Package
  X-Windows library routines
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#include <X11/Core.h>
#include <X11/cursorfont.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#include "wpack.h"

/* UGLY PATCH BY EWM TO GET A COUPLE OF BUTTON PACKAGE OBJECTS VISIBLE
   TO THE UNDERLYING WPACK.  EVENTUALLY REDO THIS.                     */

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
                     char         *label;
                     void        (*free_routine)(struct obj *);
                     void        (*draw_routine)(struct obj *);
                     void        (*event_routine)(struct obj *);
                     long          event_data;
                   } MT_OBJECT;

#define BITMAP  3

/* END UGLY PATCH */

#define STD_FONT "fixed"

/* the following are for peeking inside core widget structures */
/* because I can't seem to get XtGetValues to work */
#include <X11/IntrinsicP.h>
#include <X11/CoreP.h>

/*---------------------------------------------*/

/* structures underlying the opaque pointers passed outside */

/* eventually need to unify structures for windows and bitmaps */

typedef struct bitmap_struct {          /* bitmap */
  Pixmap pm;
  GC    gc;
  int   width,height;
  int   clip_x1,clip_y1,clip_x2,clip_y2;
  int   origin_x,origin_y;
  int   is_clipped;                     /* If a new clipping region is used */
} bitmap_struct, *bitmap_ptr;

typedef struct {
  Cursor cur;
} cursor_struct, *cursor_ptr;

typedef XFontStruct *font_ptr;          /* font */

typedef struct win_struct {             /* window */
  Widget shell;
  Widget widget;
  Window xwin;
  Cursor cur_cursor;
  GC    gc;
  Pixmap pm;
  event_handler_proc    user_event_proc;
  opaque_pointer        user_data;
  int   clip_x1,clip_y1,clip_x2,clip_y2;
  int   origin_x,origin_y;
  int   is_clipped;                     /* If a new clipping region is used */
  int   is_visible;  
  int   batch_mode;
  int   draw_mode;
  Font  fid;
  struct win_struct *next;
} win_struct, *win_ptr;

/*---------------------------------------------*/

/* macros to cast opaque pointers into internal pointers */

#define WP(win)  ((win_ptr)(win))
#define BP(bm)   ((bitmap_ptr)(bm))
#define FP(font) ((font_ptr)(font))
#define CP(curs) ((cursor_ptr)(curs))

/*---------------------------------------------*/

static XtAppContext app;        /* application context */
static Widget topwidget;        /* top-level shell widget */
static Display *dpy;            /* display pointer */
static Window root;             /* root window */
static int scrn;                /* screen number */
static int depth;               /* its depth */
static int white, black;        /* "white" and "black" pixel values */

static GC blkgc;                /* constant "black" (normal) gc */
static GC whtgc;                /* constant "white" (inverted) gc */
static GC tmpgc;                /* scratch gc for use by pattern raster ops */
static XtTranslations trans_table;

static win_ptr winlist;

static window_desc win_get(Widget);
static void expose(Widget, XEvent *, String *, Cardinal);
static void report(Widget, XEvent *, String *, Cardinal);
static long cvt_updown(unsigned int);
static int cvt_key(XKeyEvent *);
static int ReadBitmap(Display *, Drawable, char *, int *,
                      int *, Pixmap *, int *, int *);
static void resize(window_desc, int, int);
static void win_set_clip_global(window_desc,int,int,int,int);
static void bm_set_clip_global(bitmap_desc,int,int,int,int);

static Colormap my_cmap, their_cmap;
static int current_color;
static int *color_list, *color_own, cofree, cunused;
static unsigned long *unpixels;

/*---------------------------------------------*/

/* table to translate between package drawing modes and X drawing modes */

static int modetable[] = {
  /* clear      */   GXclear,
  /* ~(src|dst) */   GXnor,
  /* ~src&dst   */   GXandInverted,
  /* ~src       */   GXcopyInverted,
  /* src&~dst   */   GXandReverse,
  /* ~dst       */   GXinvert,
  /* src^dst    */   GXxor,
  /* ~(src&dst) */   GXnand,
  /* src&dst    */   GXand,
  /* ~(src^dst) */   GXequiv,
  /* dst        */   GXnoop,
  /* ~src|dst   */   GXorInverted,
  /* src        */   GXcopy,
  /* src|~dst   */   GXorReverse,
  /* src|dst    */   GXor,
  /* set        */   GXset
};

/* Table to translate cursor names to X defaults (2 * index is the X name)   */
char *cursor_names[] = {
   "X_cursor", "arrow", "based_arrow_down", "based_arrow_up", "boat",
   "bogosity", "bottom_left_corner", "bottom_right_corner", "bottom_side",
   "bottom_tee", "box_spiral", "center_ptr", "circle", "clock", "coffee_mug",
   "cross", "cross_reverse", "crosshair", "diamond_cross", "dot", "dotbox",
   "double_arrow", "draft_large", "draft_small", "draped_box", "exchange",
   "fleur", "gobbler", "gumby", "hand1", "hand2", "heart", "icon", "iron_cross",
   "left_ptr", "left_side", "left_tee", "leftbutton", "ll_angle", "lr_angle",
   "man", "middlebutton", "mouse", "pencil", "pirate", "plus", "question_arrow",
   "right_ptr", "right_side", "right_tee", "rightbutton", "rtl_logo",
   "sailboat", "sb_down_arrow", "sb_h_double_arrow", "sb_left_arrow",
   "sb_right_arrow", "sb_up_arrow", "sb_v_double_arrow", "shuttle", "sizing",
   "spider", "spraycan", "star", "target", "tcross", "top_left_arrow",
   "top_left_corner", "top_right_corner", "top_side", "top_tee", "trek",
   "ul_angle", "umbrella", "ur_angle", "watch", "xterm",
   0
};

void print_GC(GC gc)
{ XGCValues gcval;

  XGetGCValues(dpy,gc,
     GCFunction|GCPlaneMask|GCForeground|GCBackground|
                                     GCLineWidth|GCLineStyle,&gcval);
  printf("Fct = %x, PM = %lx, Fore = %lx, Back = %lx, Width = %d, Style = %x\n",
         gcval.function,gcval.plane_mask,gcval.foreground,gcval.background,
         gcval.line_width,gcval.line_style);
}

/*---------------------------------------------------------------------------*/

int color_mapped;

void set_up_colors(int depth)
{ static unsigned long *pixels, planes;
  int i, ncol;
  XColor fill;

  ncol = (1<<depth);

  color_list = (int *) malloc(sizeof(int)*ncol);
  color_own  = (int *) malloc(sizeof(int)*ncol);
  unpixels = (unsigned long *) malloc(sizeof(unsigned long)*ncol);
  pixels = (unsigned long *) malloc(sizeof(unsigned long)*ncol);

  for (cunused = ncol; cunused > 0; cunused--)
    if (XAllocColorCells(dpy,their_cmap,0,&planes,0,unpixels,cunused) != 0)
      break;
  
  for (i = 0; i < ncol; i++)
    color_own[i] = 0;
  for (i = 0; i < cunused; i++)
    color_own[unpixels[i]] = 1;

/*
  printf("There are %d free colors\n",cunused);
  fill.flags = (DoRed | DoGreen | DoBlue);
  printf("The ones in use are:\n");
  for (i = 0; i < ncol; i++)
    if (color_own[i] == 0)
      { fill.pixel = i;
        XQueryColor(dpy,their_cmap,&fill);
        printf("%3d: %d,%d,%d\n",i,fill.red,fill.green,fill.blue);
      }
*/

  if (XAllocColorCells(dpy,my_cmap,0,&planes,0,pixels,ncol) == 0)
    { printf("Serious problem, let Gene know\n"); exit(1); }

  cofree = -1;
  for (i = 0; i < ncol; i++)
    if (color_own[i] == 0)
      { fill.pixel = i;
        XQueryColor(dpy,their_cmap,&fill);
        XStoreColor(dpy,my_cmap,&fill);
        if (i != black && i != white)
          { color_list[i] = cofree;
            cofree = i;
          }
      }
  for (i = ncol-1; i >= 0; i--)
    if (color_own[i] == 1)
      { color_list[i] = cofree;
        cofree = i;
      }

/*
  printf("black = %3d\n",black);
  printf("white = %3d\n",white);
  for (i = cofree; i >= 0; i = color_list[i])
    printf("  %3d\n",i);
*/

  free(pixels);
} 

static void init_callback_table(void);

void winpack_init()
  {
    static int argc = 1;
    static char *argv[] = {"winpack", 0};

    static XtActionsRec actions[] = {
      {"Expose",(XtActionProc)expose},
      {"Report",(XtActionProc)report},
    };
    static String translations = "\
        <Expose>:    Expose()\n\
        <Enter>:     Report()\n\
        <Motion>:    Report()\n\
        <Leave>:     Report()\n\
        <BtnDown>:   Report()\n\
        <BtnUp>:     Report()\n\
        <KeyDown>:   Report()\n\
        <KeyUp>:     Report()\n\
        <Configure>: Report()\n\
        ";

    XGCValues gcval;


    XtToolkitInitialize();
    app = XtCreateApplicationContext();
    XtAppAddActions(app,actions,sizeof(actions)/sizeof(actions[0]));
    trans_table = XtParseTranslationTable(translations);

    dpy = XtOpenDisplay(app,NULL,"winpack","Winpack",NULL,0,&argc,argv);
    if (!dpy)
      { fprintf(stderr,"winpack: can't open X-windows display\n"); exit(1); }
    scrn = DefaultScreen(dpy);
    black = BlackPixel(dpy,scrn);
    white = WhitePixel(dpy,scrn);

    topwidget = XtAppCreateShell(NULL,"winpack",applicationShellWidgetClass,
      dpy,NULL,0);
    depth = topwidget->core.depth;

    color_mapped = 0;
    their_cmap = DefaultColormap(dpy,scrn);

/*
{ Visual *viz;
  XVisualInfo vt, *vta;
  int na, i;

  viz = DefaultVisual(dpy,scrn);
  vt.visualid = XVisualIDFromVisual(viz);
  vta = XGetVisualInfo(dpy,VisualIDMask,&vt,&na);

  printf("Classes:\n");
  printf("  StaticGray = %d\n",StaticGray);
  printf("  StaticColor = %d\n",StaticColor);
  printf("  TrueColor = %d\n",TrueColor);
  printf("  GrayScale = %d\n",GrayScale);
  printf("  PseudoColor = %d\n",PseudoColor);
  printf("  DirectColor = %d\n",DirectColor);

  printf("Got %d visuals\n",na);
  for (i = 0; i < na; i++)
    printf("Class = %d  depth = %d  masks = (%x,%x,%x) mapsize = %d\n",
           vta[i].class,vta[i].depth,vta[i].red_mask,vta[i].green_mask,
           vta[i].blue_mask,vta[i].colormap_size);
  XFree(vta);
}
*/

    if (depth <= 10) {
        my_cmap = XCreateColormap ( dpy, XDefaultRootWindow (dpy),
                                          XDefaultVisual (dpy, scrn),
                                          AllocNone );
        topwidget->core.colormap = my_cmap;
        set_up_colors(depth);
        color_mapped = 1;
    }
    current_color = black;

    gcval.foreground = black;
    gcval.background = white;
    gcval.cap_style = CapButt;
    gcval.join_style = JoinMiter;
    root = RootWindow(dpy,scrn);
    blkgc = XCreateGC(dpy,root,GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
            &gcval);
    tmpgc = XCreateGC(dpy,root,GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
            &gcval);

    gcval.foreground = white;
    gcval.background = black;
    whtgc = XCreateGC(dpy,root,GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
            &gcval);
    init_callback_table();
  }

/*---------------------------------------------*/

void winpack_main_loop()
  {
    XEvent event;

    while (app)  
      {
        XtAppNextEvent(app,&event);
        XtDispatchEvent(&event);
      }
  }

/*---------------------------------------------*/

/* File input callback stuff */
#define CALLBACKTABLESIZE 64

struct {
   int fd;
   file_callback_proc proc;
   XtInputId x_id;
} callback_table[CALLBACKTABLESIZE];

static void init_callback_table(void)
{
   int i;
   for(i = 0; i < CALLBACKTABLESIZE; i++)
      callback_table[i].fd = -1;
}

static int find_callback(int fd)
{
   int i;
   for(i = 0; i < CALLBACKTABLESIZE; i++)
      if(callback_table[i].fd == fd)
	 return i;
   return -1;
}

static int get_callback_slot(void)
{
   int i;
   for(i = 0; i < CALLBACKTABLESIZE; i++)
      if(callback_table[i].fd == -1)
	 return i;
   return -1;
}

/* The procedure that gets the Xt callbacks and then calls the user proc */
static void gen_callback_proc(int fd, int *junque, XtInputId xid)
{
   int i;
   i = find_callback(fd);
   /* <<< COMMENTED OUT >>>
   if(callback_table[i].x_id != xid || fd != *junque) {
      fprintf(stderr, "Discrepancy in gen_callback_proc(): ");
      if(callback_table[i].x_id != xid)
	 fprintf(stderr, "XtInputId mismatch ");
      else
	 fprintf(stderr, "client_data source mismatch ");
      fprintf(stderr, "for fd = %d\n", fd);
   }
    */

   (*callback_table[i].proc)(fd);
}

/* Commented out by E. Myers, Jan. 25, '93

void winpack_monitor_fd(int fd, file_callback_proc user_proc)
{
   int slot;
   if((slot = get_callback_slot()) < 0) {
      fprintf(stderr, "Ran out of space for fd monitoring... fd = %d\n", fd);
      return;
   }
   printf("Using slot %d for fd %d.\n", slot, fd);

   callback_table[slot].fd = fd;
   callback_table[slot].proc = user_proc;
   callback_table[slot].x_id =
      XtAppAddInput(app, fd, XtInputReadMask, gen_callback_proc, (caddr_t) fd);
}
*/

void winpack_unmonitor_fd(int fd)
{
   /* Search for the XtInputId */
   int slot;
   slot = find_callback(fd);
   if(slot < 0) return;
   
   XtRemoveInput(callback_table[slot].x_id);
   callback_table[slot].fd = -1;
}


/*---------------------------------------------*/

void winpack_shutdown(void)
  {
    if (!app) return;

    if (cunused > 0)
      XFreeColors(dpy,their_cmap,unpixels,cunused,0);
    XtDestroyApplicationContext(app);
    topwidget = 0;
    app = 0;
  }

/*---------------------------------------------*/

#include <string.h>
#define streq(a, b) (!strcmp((a), (b)))

cursor_desc cursor_open(filename, maskname)
char *filename, *maskname;
{
   bitmap_ptr bm, mask;
   int x, y, r;
   Cursor cc;
   cursor_ptr c;

   if(*filename == '#') {
      int i;
      /* lookup table */
      i = 0;
      while(cursor_names[i] != 0) {
         if(streq(cursor_names[i], filename+1)) {
            cc = XCreateFontCursor(dpy, i*2);
            break;
         }
         i++;
      }
      if(!cursor_names[i]) return 0;
   }
   else {
      unsigned int rwidth, rheight;
      XColor black, white;

      bm = (bitmap_ptr)XtMalloc(sizeof(bitmap_struct));
      r = ReadBitmap(dpy,root,filename,&bm->width,&bm->height,&bm->pm,&x,&y);
      if (r != BitmapSuccess)
        return 0;
      if(maskname) {
         int xx, yy;
         mask = (bitmap_ptr)XtMalloc(sizeof(bitmap_struct));
         r = ReadBitmap(dpy,root,maskname,&mask->width,&mask->height,&mask->pm,
                 &xx,&yy);
         if (r != BitmapSuccess) {
           XtFree((char *)bm);
           return 0;
         }
      }
      else
         mask = 0;

      XQueryBestCursor(dpy, root, (unsigned)bm->width, (unsigned)bm->height,
                       &rwidth, &rheight);
      if(bm->width != rwidth || bm->height != rheight) {
         XFreePixmap(dpy,bm->pm);
         XtFree((char *)bm);
         return 0;
      }

      black.red = black.green = black.blue = 0;
      black.flags = DoRed | DoGreen | DoBlue;
      black.pixel = 1; /* What the hell */
      white.red = white.green = white.blue = 65535;
      white.flags = DoRed | DoGreen | DoBlue;
      white.pixel = 1;
      if(mask) {
         cc = XCreatePixmapCursor(dpy, bm->pm, mask->pm, &black, &white, x, y);
         XFreePixmap(dpy, mask->pm);
         XtFree((char *)mask);
      }
      else
         cc = XCreatePixmapCursor(dpy, bm->pm, 0, &black, &white, x, y);
      XFreePixmap(dpy, bm->pm);
      XtFree((char *)bm);
   }

   if(cc) {
      c = (cursor_ptr) XtMalloc(sizeof(cursor_struct));
      c->cur = cc;
      return (cursor_desc)c;
   }
   else
      return 0;
}

void cursor_free(c)
cursor_desc c;
{
   if(!c) return;

   XFreeCursor(dpy, CP(c)->cur);
   XtFree((char *)c);
}

void win_set_cursor(w, c)
window_desc w;
cursor_desc c;
{
   if(!w || !c) return;

   WP(w)->cur_cursor = CP(c)->cur;
   XDefineCursor(dpy, WP(w)->xwin, CP(c)->cur);
}

/*---------------------------------------------*/

static void set_win_mode(window_desc win, drawing_mode mode)
  {
    mode &= 0xF;
    if (WP(win)->draw_mode != mode)
      {
        WP(win)->draw_mode = mode;
        XSetFunction(dpy,WP(win)->gc,modetable[mode]);
      }
    XSetForeground(dpy,WP(win)->gc,current_color);
  }

/*---------------------------------------------*/

static void set_bm_mode(bitmap_desc bm, drawing_mode mode)
  {
    mode &= 0xF;
    /****  xx need to cache draw mode for bitmaps like windows ***
    if (WP(win)->draw_mode != mode)
      {
        WP(win)->draw_mode = mode;
        XSetFunction(dpy,WP(win)->gc,modetable[mode]);
      }
    ****/
    XSetFunction(dpy,BP(bm)->gc,modetable[mode]);
    XSetForeground(dpy,BP(bm)->gc,current_color);
  }

/*---------------------------------------------*/

static window_desc win_get(Widget wid)
  {
    win_ptr p;

    for (p = winlist; p; p = p->next)
      if (p->widget == wid)
        return (window_desc) p;
    fprintf(stderr,"winpack internal error: can't map widget to window\n");
    exit(1);
    /*NOTREACHED*/
  }

/*---------------------------------------------*/

static void expose(Widget wid, XEvent *event, String *params, Cardinal n)
  {
    int w, h;
    window_desc win;

    if (event->xexpose.count > 0) return;
    win = win_get(wid);
    win_get_size(win,&w,&h);
    XCopyArea(dpy,WP(win)->pm,XtWindow(wid),blkgc,0,0,w,h,0,0);
  }

/*---------------------------------------------*/

static void report(Widget wid, XEvent *event, String *params, Cardinal n)
  {
    XMotionEvent *me = (XMotionEvent *) event;
    window_desc win = win_get(wid);
    input_event ev;

    if (!WP(win)->user_event_proc)
      return;

    ev.in_window = win;
    ev.when = me->time;
    ev.x = me->x - WP(win)->origin_x;
    ev.y = me->y - WP(win)->origin_y;
    ev.global_x = ev.x;
    ev.global_y = ev.y;
    ev.updown_mask = cvt_updown(me->state);
    ev.dev_num = 0;
    ev.value = 0;

    switch (me->type)
      {
        case EnterNotify:
          ev.code = win_enter;
        break;

        case MotionNotify:
          ev.code = loc_move;
        break;

        case LeaveNotify:
          ev.code = win_exit;
        break;

        case KeyPress:
          ev.code = key_dn;
          ev.value = cvt_key((XKeyEvent*)event);
          if (ev.value < 0)
            return;
        break;

        case KeyRelease:
          ev.code = key_up;
          ev.value = cvt_key((XKeyEvent*)event);
          if (ev.value < 0)
            return;
        break;

        case ButtonPress:
          ev.code = loc_but_dn;
          ev.value = event->xbutton.button;
        break;

        case ButtonRelease:
          ev.code = loc_but_up;
          ev.value = event->xbutton.button;
        break;

        case ConfigureNotify:
          ev.code = win_resize;
          ev.x = event->xconfigure.width;
          ev.y = event->xconfigure.height;
          resize(win, ev.x, ev.y);
          break;
      }

    (*WP(win)->user_event_proc)(&ev,win,WP(win)->user_data);
  }

/*---------------------------------------------*/

static long cvt_updown(unsigned int x_state)
  {
    long result;

    result = 0;
    if (x_state & ShiftMask)   result |= SHIFT_KEY_DOWN;
    if (x_state & ControlMask) result |= CNTRL_KEY_DOWN;
    if (x_state & Mod1Mask)    result |= META_KEY_DOWN;

    return result;
  }

/*---------------------------------------------*/

static int cvt_key(XKeyEvent *ev)
  {
    int n;
    char buf[4];

    n = XLookupString(ev,buf,4,NULL,NULL);
    if (n != 1) return -1;
    return buf[0];
  }

/*---------------------------------------------------------------------------*/

char *event_code_string(input_code cd)
  {
    static char *event_names[] = {
      "EVENT_none",	"loc_but_dn",	"loc_but_up",	"but_dn",
      "but_up",		"loc_move",	"val_move",	"key_dn",
      "key_up",		"win_enter",	"win_exit",	"win_resize",
      "special_event"
    };

    static char *unknown_event_name = "??unknown??";

    switch(cd)
      {
        case EVENT_none:	return event_names[0];
        case loc_but_dn:	return event_names[1];
        case loc_but_up:	return event_names[2];
        case but_dn:		return event_names[3];
        case but_up:		return event_names[4];
        case loc_move:		return event_names[5];
        case val_move:		return event_names[6];
        case key_dn:		return event_names[7];
        case key_up:		return event_names[8];
        case win_enter:		return event_names[9];
        case win_exit:		return event_names[10];
        case win_resize:	return event_names[11];
        case special_event:	return event_names[12];

	default:		return unknown_event_name;
      }
  }

/*---------------------------------------------------------------------------*/

char *event_string(input_event *evt)
  {
    static char buffer[512],str[512];
    static char *null_message = "NULL EVENT";

    if (!evt) return null_message;

    /* clear the buffer */
    buffer[0] ='\0';

    strncat(buffer,event_code_string(evt->code),511);

    sprintf(str,":(%ld,%ld) [0x%lx] ",evt->x, evt->y, (long)evt->in_window);
    strncat(buffer,str,511);

    str[0] = '\0';
    if (evt->updown_mask & SHIFT_KEY_DOWN) strcat(str,"SHIFT ");
    if (evt->updown_mask & CNTRL_KEY_DOWN) strcat(str,"CNTRL ");
    if (evt->updown_mask & META_KEY_DOWN) strcat(str,"META ");
    strncat(buffer,str,511);

    sprintf(str,"dev#%d val=%ld",evt->dev_num,evt->value);
    strncat(buffer,str,511);

    return buffer;
  }


/*---------------------------------------------------------------------------*/

window_desc win_new(x,y,w,h,event_proc,data,title)
  int x,y,w,h;
  event_handler_proc event_proc;
  opaque_pointer data;
  char *title;
  {
    char buf[100];
    static int num;
    win_ptr result;
    XGCValues gcval;
    Arg args[25];
    int nargs;
    XWMHints wmhints; /* [SEH] */


    if (h < 0) h = -h;
    if (w < 0) w = -w;

    result = (win_ptr)XtMalloc(sizeof(win_struct));
    result->user_event_proc = event_proc;
    result->user_data = data;
    result->clip_x1 = 0;
    result->clip_y1 = 0;
    result->clip_x2 = w;
    result->clip_y2 = h;
    result->is_clipped = 0;
    result->origin_x   = 0;
    result->origin_y   = 0;
    result->draw_mode = -1;
    result->fid = 0;

    nargs = 0;
    XtSetArg(args[nargs],XtNx,x);       nargs++;
    XtSetArg(args[nargs],XtNy,y);       nargs++;
    XtSetArg(args[nargs],XtNwidth,w);   nargs++;
    XtSetArg(args[nargs],XtNheight,h);  nargs++;

    num++;
    if (!title) title = "   ";
    result->shell = XtCreatePopupShell(title,
      applicationShellWidgetClass,topwidget,args,nargs);

    sprintf(buf,"win%d",num);
    result->widget = XtCreateManagedWidget(buf,
      widgetClass,result->shell,NULL,0);
    XtOverrideTranslations(result->widget,trans_table);

    XtRealizeWidget(result->shell);

    result->xwin = XtWindow(result->widget);

    /* tell the window manager that we want passive input [SEH] */
    wmhints.flags = InputHint;
    wmhints.input = 1;
    XSetWMHints(dpy,XtWindow(result->shell),&wmhints);

    result->pm = XCreatePixmap(dpy,XtWindow(result->widget),w,h,depth);
    result->cur_cursor = 0;

    gcval.foreground = black; 
    gcval.background = white;
    gcval.cap_style  = CapButt;
    gcval.join_style = JoinMiter;
    result->gc = XCreateGC(dpy,result->xwin,
                           GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
                           &gcval);
    result->is_visible = 0;
    result->batch_mode = 0;
    XFillRectangle(dpy,result->pm,whtgc,0,0,w,h);

    result->next = winlist;
    winlist = result;

    set_win_mode((window_desc)result,MODE_SRC);

    return (window_desc)result;
  }

/*---------------------------------------------*/

void win_free(win)
  window_desc win;
  {
    win_ptr p;

    if (!win) return;

    if (winlist == WP(win))
      winlist = WP(win)->next;
    else
      {
        for (p = winlist; p && p->next != WP(win); p = p->next)
          ;
        if (p)
          p->next = WP(win)->next;
      }
      
    XFreePixmap(dpy,WP(win)->pm);
    XFreeGC(dpy,WP(win)->gc);
    XtDestroyWidget(WP(win)->shell);
  }

/*---------------------------------------------*/

void  win_set_vis(win,v)
  window_desc win;
  int          v;
  {
    if (!win) return;
    if (v && !WP(win)->is_visible)
      {
        XtPopup(WP(win)->shell,XtGrabNone);
        WP(win)->is_visible = 1;
      }
    else if (!v && WP(win)->is_visible)
      {
        XtPopdown(WP(win)->shell);
        WP(win)->is_visible = 0;
      }
    XFlush(dpy);
  }

/*---------------------------------------------*/

int  win_get_vis(win)
  window_desc win;
  {
    if (!win) return 0;
    return WP(win)->is_visible;
  }

/*---------------------------------------------*/

void  win_set_data(win,data)
  window_desc win;
  opaque_pointer data;
  {
    if (!win) return;
    WP(win)->user_data = data; 
  }

/*---------------------------------------------*/

opaque_pointer win_get_data(win)
  window_desc win;
  {
    if (!win) return 0;
    return WP(win)->user_data;
  }
  
/*---------------------------------------------*/

void win_set_handler(win, p)
  window_desc win; 
  event_handler_proc p;
  {
    if (!win) return;
    WP(win)->user_event_proc = p;
  }

/*---------------------------------------------*/

event_handler_proc win_get_handler(win)
  window_desc win;
  {
    if (!win) return 0;
    return WP(win)->user_event_proc;
  }

/*---------------------------------------------*/

void win_set_pos(win,x,y)
  window_desc win;
  int x,y;
  {
    Arg args[2];

    if (!win) return;
    XtSetArg(args[0],XtNx,x);
    XtSetArg(args[1],XtNy,y);
    XtSetValues(WP(win)->shell,args,2);
  }
  
/*---------------------------------------------*/

void win_get_pos(win, x_result, y_result)
  window_desc win;
  int *x_result;
  int *y_result;
  {
    if (!win || !x_result || !y_result) return;
    *x_result = WP(win)->shell->core.x;
    *y_result = WP(win)->shell->core.y;
  }
  
/*---------------------------------------------*/

void win_set_size(win,w,h)
  window_desc win;
  int w,h;
  {
    Pixmap pm;
    Arg args[2];

    if (!win) return;
    win_reset_clip(win);
    XtSetArg(args[0],XtNwidth,w);
    XtSetArg(args[1],XtNheight,h);
    XtSetValues(WP(win)->shell,args,2);

    pm = XCreatePixmap(dpy,XtWindow(WP(win)->widget),w,h,depth);
    XFillRectangle(dpy,pm,whtgc,0,0,w,h);
    XCopyArea(dpy,WP(win)->pm,pm,blkgc,0,0,w,h,0,0);
    XFreePixmap(dpy,WP(win)->pm);
    WP(win)->pm = pm;
    win_reset_clip(win);
  }

/* An internal stub to do the right thing if user resizes the window. */
static void resize(window_desc win, int w, int h)
{
    Pixmap pm;

    pm = XCreatePixmap(dpy,XtWindow(WP(win)->widget),w,h,depth);
    XFillRectangle(dpy,pm,whtgc,0,0,w,h);
    XCopyArea(dpy,WP(win)->pm,pm,blkgc,0,0,w,h,0,0);
    XFreePixmap(dpy,WP(win)->pm);
    WP(win)->pm = pm;
    /* If a clipping region is in effect, don't touch it */
    if(!WP(win)->is_clipped)
       win_reset_clip(win);
}

/*---------------------------------------------*/

void win_get_size(win, w_result, h_result)
  window_desc win;
  int *w_result;
  int *h_result;
  {
    if (!win || !w_result || !h_result) return;
    *w_result = WP(win)->widget->core.width;
    *h_result = WP(win)->widget->core.height;
  }

/*---------------------------------------------*/

int win_width(win)
  window_desc win;
  {
    if (!win) return 0;
    return WP(win)->widget->core.width;
  }

/*---------------------------------------------*/

int win_height(win)
  window_desc win;
  {
    if (!win) return 0;
    return WP(win)->widget->core.height;
  }

/*---------------------------------------------*/

void win_reduce_clip(win,x1,y1,x2,y2)
  window_desc win;
  int x1,y1,x2,y2;
  {
    if (!win) return;

    /* translate clip region into global window coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    /* intersect old clip area with new one */
    if (WP(win)->clip_x1 > x1) x1 = WP(win)->clip_x1;
    if (WP(win)->clip_y1 > y1) y1 = WP(win)->clip_y1;
    if (WP(win)->clip_x2 < x2) x2 = WP(win)->clip_x2;
    if (WP(win)->clip_y2 < y2) y2 = WP(win)->clip_y2;

    /* set the region */
    win_set_clip_global(win,x1,y1,x2,y2);
  }

/*---------------------------------------------*/

void bm_reduce_clip(bm,x1,y1,x2,y2)
  bitmap_desc bm;
  int x1,y1,x2,y2;
  {
    if (!bm) return;

    /* translate clip region into global window coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    /* intersect old clip area with new one */
    if (BP(bm)->clip_x1 > x1) x1 = BP(bm)->clip_x1;
    if (BP(bm)->clip_y1 > y1) y1 = BP(bm)->clip_y1;
    if (BP(bm)->clip_x2 < x2) x2 = BP(bm)->clip_x2;
    if (BP(bm)->clip_y2 < y2) y2 = BP(bm)->clip_y2;

    /* set the region */
    bm_set_clip_global(bm,x1,y1,x2,y2);
  }

/*------------------------------------------------------------*/

static void win_set_clip_global(window_desc win, int x1, int y1, int x2, int y2)
  {
    XRectangle r;

    if (!win) return;

    /* if empty give null list else set clipping rectangle */
    if (x1 > x2 || y1 > y2)
      {
        XSetClipRectangles(dpy,WP(win)->gc,0,0,NULL,0,YXBanded);
      }
    else
      {
        r.x = 0;
        r.y = 0;
        r.width = x2 - x1;
        r.height = y2 - y1;
        XSetClipRectangles(dpy,WP(win)->gc,x1,y1,&r,1,YXBanded);
      }

    WP(win)->clip_x1 = x1; WP(win)->clip_y1 = y1;
    WP(win)->clip_x2 = x2; WP(win)->clip_y2 = y2;
    WP(win)->is_clipped = 1;
  }

void mt_begin_clip_to(MT_OBJECT *o)
{ XRectangle r;
  window_desc win;

  win = o->wpwin;

  /* if empty give null list else set clipping rectangle */
  if (o->xl > o->xh || o->yl > o->yh)
    XSetClipRectangles(dpy,WP(win)->gc,0,0,NULL,0,YXBanded);
  else
    { r.x = o->outline + o->border;
      r.y = r.x;
      r.width = (o->xh - o->xl) + (1 - 2*r.x);
      r.height = (o->yh - o->yl) + (1 - 2*r.x);
      XSetClipRectangles(dpy,WP(win)->gc,o->xl,o->yl,&r,1,YXBanded);
    }

  WP(win)->clip_x1 = o->xl; WP(win)->clip_y1 = o->yl;
  WP(win)->clip_x2 = o->xh; WP(win)->clip_y2 = o->yh;
  WP(win)->is_clipped = 1;
}

/*------------------------------------------------------------*/

void win_set_clip(win,x1,y1,x2,y2)
  window_desc win;
  int x1,y1,x2,y2;
  {

    if (!win) return;

    /* translate clip region into global window coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    /* set the region */
    win_set_clip_global(win,x1,y1,x2,y2);
  }

/*------------------------------------------------------------*/

static void bm_set_clip_global(bitmap_desc bm, int x1, int y1, int x2, int y2)
  {
    XRectangle r;

    if (!bm) return;

    /* if empty give null list else set clipping rectangle */
    if (x1 > x2 || y1 > y2)
      {
        XSetClipRectangles(dpy,BP(bm)->gc,0,0,NULL,0,YXBanded);
      }
    else
      {
        r.x = 0;
        r.y = 0;
        r.width = x2 - x1;
        r.height = y2 - y1;
        XSetClipRectangles(dpy,BP(bm)->gc,x1,y1,&r,1,YXBanded);
      }

    BP(bm)->clip_x1 = x1; BP(bm)->clip_y1 = y1;
    BP(bm)->clip_x2 = x2; BP(bm)->clip_y2 = y2;
    BP(bm)->is_clipped = 1;
  }

/*------------------------------------------------------------*/

void bm_set_clip(bm,x1,y1,x2,y2)
  bitmap_desc bm;
  int x1,y1,x2,y2;
  {

    if (!bm) return;

    /* translate clip region into global window coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    /* set the region */
    bm_set_clip_global(bm,x1,y1,x2,y2);
  }

/*------------------------------------------------------------*/

void win_get_clip(win,x1_result,y1_result,x2_result,y2_result)
  window_desc win;
  int *x1_result,*y1_result,*x2_result,*y2_result;
  {
    if (!win || !x1_result || !y1_result || !x2_result || !y2_result) return;
    
    /* return clip region translated into local coordinates */
    *x1_result = WP(win)->clip_x1 - WP(win)->origin_x; 
    *y1_result = WP(win)->clip_y1 - WP(win)->origin_y;
    *x2_result = WP(win)->clip_x2 - WP(win)->origin_x; 
    *y2_result = WP(win)->clip_y2 - WP(win)->origin_y;
  }

/*------------------------------------------------------------*/

void bm_get_clip(bm,x1_result,y1_result,x2_result,y2_result)
  bitmap_desc bm;
  int *x1_result,*y1_result,*x2_result,*y2_result;
  {
    if (!bm || !x1_result || !y1_result || !x2_result || !y2_result) return;
    
    /* return clip region translated into local coordinates */
    *x1_result = BP(bm)->clip_x1 - BP(bm)->origin_x; 
    *y1_result = BP(bm)->clip_y1 - BP(bm)->origin_y;
    *x2_result = BP(bm)->clip_x2 - BP(bm)->origin_x; 
    *y2_result = BP(bm)->clip_y2 - BP(bm)->origin_y;
  }

/*------------------------------------------------------------*/

void mt_end_clip_to(MT_OBJECT *o)
{ window_desc win;
  int w, h;

  win = o->wpwin;
  if (!win) return;
  win_get_size(win,&w,&h);
  XSetClipMask(dpy,WP(win)->gc,None);
  WP(win)->clip_x1 = 0;      
  WP(win)->clip_y1 = 0;
  WP(win)->clip_x2 = w;
  WP(win)->clip_y2 = h;
  WP(win)->is_clipped = 0;
}

void win_reset_clip(win)
  window_desc win;
  {
    int w, h;

    if (!win) return;
    win_get_size(win,&w,&h);
    XSetClipMask(dpy,WP(win)->gc,None);
    WP(win)->clip_x1 = 0;      
    WP(win)->clip_y1 = 0;
    WP(win)->clip_x2 = w;
    WP(win)->clip_y2 = h;
    WP(win)->is_clipped = 0;
  }

/*------------------------------------------------------------*/

void bm_reset_clip(bm)
  bitmap_desc bm;
  {
    int w, h;

    if (!bm) return;
    w = bm_width(bm);
    h = bm_height(bm);
    XSetClipMask(dpy,BP(bm)->gc,None);
    BP(bm)->clip_x1 = 0;      
    BP(bm)->clip_y1 = 0;
    BP(bm)->clip_x2 = w;
    BP(bm)->clip_y2 = h;
    BP(bm)->is_clipped = 0;
  }

/*---------------------------------------------*/

void win_reset_origin(win)
  window_desc win;
  {
    if (!win) return;
    WP(win)->origin_x = 0;
    WP(win)->origin_y = 0;
  }

/*---------------------------------------------*/

void bm_reset_origin(bm)
  bitmap_desc bm;
  {
    if (!bm) return;
    BP(bm)->origin_x = 0;
    BP(bm)->origin_y = 0;
  }

/*---------------------------------------------*/

void win_set_origin(win, x, y)
  window_desc win;
  int x;
  int y;
  {
    if (!win) return;
    WP(win)->origin_x += x;
    WP(win)->origin_y += y;
  }
  
/*---------------------------------------------*/

void bm_set_origin(bm, x, y)
  bitmap_desc bm;
  int x;
  int y;
  {
    if (!bm) return;
    BP(bm)->origin_x += x;
    BP(bm)->origin_y += y;
  }
  
/*---------------------------------------------*/

void win_restore_origin(win, x, y)
  window_desc win;
  int x;
  int y;
  {
    if (!win) return;
    WP(win)->origin_x = x;
    WP(win)->origin_y = y;
  }
  
/*---------------------------------------------*/

void bm_restore_origin(bm, x, y)
  bitmap_desc bm;
  int x;
  int y;
  {
    if (!bm) return;
    BP(bm)->origin_x = x;
    BP(bm)->origin_y = y;
  }

/*---------------------------------------------*/

void win_get_origin(win, x, y)
  window_desc win;
  int *x;
  int *y;
  {
    if (!win) return;
    if (x) *x = WP(win)->origin_x;
    if (y) *y = WP(win)->origin_y;
  }
  
/*---------------------------------------------*/

void bm_get_origin(bm, x, y)
  bitmap_desc bm;
  int *x;
  int *y;
  {
    if (!bm) return;
    if (x) *x = BP(bm)->origin_x;
    if (y) *y = BP(bm)->origin_y;
  }
  
/*---------------------------------------------*/

int win_trivial_reject(win, x1, y1, x2, y2)
  window_desc win;
  int x1; 
  int y1;
  int x2;
  int y2;
  {
    if (!win) return (1);

    /* translate from local to global coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    return 
      (y1 > WP(win)->clip_y2 && y2 > WP(win)->clip_y2) ||
      (y1 < WP(win)->clip_y1 && y2 < WP(win)->clip_y1) ||
      (x1 > WP(win)->clip_x2 && x2 > WP(win)->clip_x2) ||
      (x1 < WP(win)->clip_x1 && x2 < WP(win)->clip_x1);
  }

/*---------------------------------------------*/

int bm_trivial_reject(bm, x1, y1, x2, y2)
  bitmap_desc bm;
  int x1; 
  int y1;
  int x2;
  int y2;
  {
    if (!bm) return (1);

    /* translate from local to global coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    return 
      (y1 > BP(bm)->clip_y2 && y2 > BP(bm)->clip_y2) ||
      (y1 < BP(bm)->clip_y1 && y2 < BP(bm)->clip_y1) ||
      (x1 > BP(bm)->clip_x2 && x2 > BP(bm)->clip_x2) ||
      (x1 < BP(bm)->clip_x1 && x2 < BP(bm)->clip_x1);
  }

/*---------------------------------------------*/

void win_clear_rect(win,x1,y1,x2,y2) 
  window_desc win;
  int x1,y1,x2,y2;
  {
    int h,w;

    if (!win) return;

    /* translate from local to global coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    if (x2 <= x1) return;
    if (y2 <= y1) return;
    if (x1 > WP(win)->clip_x2) return;
    if (x2 < WP(win)->clip_x1) return;
    if (y1 > WP(win)->clip_y2) return;
    if (y2 < WP(win)->clip_y1) return;

    if (x1 < WP(win)->clip_x1)  x1 = WP(win)->clip_x1;
    if (x2 > WP(win)->clip_x2)  x2 = WP(win)->clip_x2;
    if (y1 < WP(win)->clip_y1)  y1 = WP(win)->clip_y1;
    if (y2 > WP(win)->clip_y2)  y2 = WP(win)->clip_y2;

    w = x2 - x1;  if (w <= 0)  return;
    h = y2 - y1;  if (h <= 0)  return;

    if (!WP(win)->batch_mode)
      XFillRectangle(dpy,WP(win)->xwin,whtgc,x1,y1,w,h);
    XFillRectangle(dpy,WP(win)->pm,whtgc,x1,y1,w,h);
  }

/*---------------------------------------------*/

void bm_clear_rect(bm,x1,y1,x2,y2) 
  bitmap_desc bm;
  int x1,y1,x2,y2;
  {
    int h,w;

    if (!bm) return;

    /* translate from local to global coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    if (x2 <= x1) return;
    if (y2 <= y1) return;
    if (x1 > BP(bm)->clip_x2) return;
    if (x2 < BP(bm)->clip_x1) return;
    if (y1 > BP(bm)->clip_y2) return;
    if (y2 < BP(bm)->clip_y1) return;

    if (x1 < BP(bm)->clip_x1)  x1 = BP(bm)->clip_x1;
    if (x2 > BP(bm)->clip_x2)  x2 = BP(bm)->clip_x2;
    if (y1 < BP(bm)->clip_y1)  y1 = BP(bm)->clip_y1;
    if (y2 > BP(bm)->clip_y2)  y2 = BP(bm)->clip_y2;


    w = x2 - x1;  if (w < 0) return;
    h = y2 - y1;  if (h < 0) return;

    XFillRectangle(dpy,BP(bm)->pm,whtgc,x1,y1,w,h);
  }

/*---------------------------------------------*/

void win_start_batch(win)
  window_desc win;
  {
    if (!win) return;

    WP(win)->batch_mode = 1;
  }

/*---------------------------------------------*/

void win_end_batch(win)
  window_desc win;
  {
    int w, h;

    if (!win) return;
    if (WP(win)->batch_mode) 
      {
        win_get_size(win,&w,&h);
        XCopyArea(dpy,WP(win)->pm,WP(win)->xwin,blkgc,0,0,w,h,0,0);
        WP(win)->batch_mode = 0;
      }
    XFlush(dpy);
  }
/*---------------------------------------------------------------------------*/

void win_draw_line(win,x1,y1,x2,y2,mode)
  window_desc win;
  int x1,y1,x2,y2;
  drawing_mode mode;
  {
    int temp;

    if (!win) return;

    /* translate from local to global coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    if(x2 < x1) {
       temp = x2; x2 = x1; x1 = temp;
       temp = y2; y2 = y1; y1 = temp;
    }
    else if(y2 < y1) {
       temp = x2; x2 = x1; x1 = temp;
       temp = y2; y2 = y1; y1 = temp;
    }

    set_win_mode(win,mode);

#ifdef DEBUG
    printf("WLine: ");
    print_GC(WP(win)->gc);
#endif

    if (!WP(win)->batch_mode)
      XDrawLine(dpy,WP(win)->xwin,WP(win)->gc,x1,y1,x2,y2);
    XDrawLine(dpy,WP(win)->pm,WP(win)->gc,x1,y1,x2,y2);
  }

void mt_draw_line(MT_OBJECT *o, int xl, int yl, int xh, int yh,
                                int thick, int dash)
{ window_desc win;
  bitmap_desc bm;

  xl += o->xl; yl += o->yl;
  xh += o->xl; yh += o->yl;

  if (o->type == BITMAP)
    { bm = (bitmap_desc) (o->wpwin);
      set_bm_mode(bm,MODE_SRC);
      if (dash)
        XSetLineAttributes(dpy,BP(bm)->gc,thick,
                           LineOnOffDash,CapButt,JoinMiter);
      else
        XSetLineAttributes(dpy,BP(bm)->gc,thick,LineSolid,CapButt,JoinMiter);

#ifdef DEBUG
      printf("BLine: ");
      print_GC(BP(bm)->gc);
#endif

      XDrawLine(dpy,BP(bm)->pm,BP(bm)->gc,xl,yl,xh,yh);

      XSetLineAttributes(dpy,BP(bm)->gc,0,LineSolid,CapButt,JoinMiter);
    }
  else
    { win = o->wpwin;
      set_win_mode(win,MODE_SRC);
      if (dash)
        XSetLineAttributes(dpy,WP(win)->gc,thick,
                           LineOnOffDash,CapButt,JoinMiter);
      else
        XSetLineAttributes(dpy,WP(win)->gc,thick,LineSolid,CapButt,JoinMiter);

#ifdef DEBUG
      printf("WLine: ");
      print_GC(WP(win)->gc);
#endif

      if (!WP(win)->batch_mode)
        XDrawLine(dpy,WP(win)->xwin,WP(win)->gc,xl,yl,xh,yh);
      XDrawLine(dpy,WP(win)->pm,WP(win)->gc,xl,yl,xh,yh);

      XSetLineAttributes(dpy,WP(win)->gc,0,LineSolid,CapButt,JoinMiter);
    }
}

/*---------------------------------------------------------------------------*/

void bm_draw_line(bm,x1,y1,x2,y2,mode)
  bitmap_desc bm;
  int x1,y1,x2,y2;
  drawing_mode mode;
  {
    int temp;

    if (!bm) return;

    /* translate from local to global coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    if(x2 < x1) {
       temp = x2; x2 = x1; x1 = temp;
       temp = y2; y2 = y1; y1 = temp;
    }
    else if(y2 < y1) {
       temp = x2; x2 = x1; x1 = temp;
       temp = y2; y2 = y1; y1 = temp;
    }
    set_bm_mode(bm,mode);

#ifdef DEBUG
    printf("BLine: ");
    print_GC(BP(bm)->gc);
#endif

    XDrawLine(dpy,BP(bm)->pm,BP(bm)->gc,x1,y1,x2,y2);
  }

/*---------------------------------------------------------------------------*/

font_desc font_open(from_file)
  char *from_file;
  {
    XFontStruct *result;

    if (!from_file) return 0;
    result = XLoadQueryFont(dpy,from_file);
    return (font_desc)result;
  }

/*---------------------------------------------*/

font_desc font_default(void)
  {
    font_desc result;
    char *fontname;
    
    fontname = XGetDefault(dpy,"winpack","Font");
    if (!fontname) fontname = STD_FONT;
    result = font_open(fontname);
    return result;
  }

/*---------------------------------------------*/

void font_free(font)
  font_desc font;
  {
    if (!font) return;
    XFreeFont(dpy,(XFontStruct*)font);
  }

/*---------------------------------------------*/

void font_string_size(font,str,w_result,h_result,baseline_result)
  font_desc      font;
  char          *str;
  int           *w_result;
  int           *h_result;
  int           *baseline_result;
  {
    if (!font || !str || !w_result || !h_result || !baseline_result)
      return;
    
    *w_result = XTextWidth((XFontStruct*)font,str,strlen(str));
    *h_result = ((XFontStruct*)font)->ascent + ((XFontStruct*)font)->descent;
    *baseline_result = ((XFontStruct*)font)->descent;
  }

void mt_string_size(char *str, int *w_result, int *h_result, int *baseline_result)
  { font_desc font;

    if (!str || !w_result || !h_result || !baseline_result)
      return;
    
    font = font_default();
    *w_result = XTextWidth((XFontStruct*)font,str,strlen(str));
    *h_result = ((XFontStruct*)font)->ascent + ((XFontStruct*)font)->descent;
    *baseline_result = ((XFontStruct*)font)->descent;
  }

/*---------------------------------------------*/

extern int XDrawImageString(
#if NeedFunctionPrototypes
    Display*            /* display */,
    Drawable            /* d */,
    GC                  /* gc */,
    int                 /* x */,
    int                 /* y */,
    const char*         /* string */,
    int                 /* length */
#endif
);    

extern int XDrawString(
#if NeedFunctionPrototypes
    Display*            /* display */,
    Drawable            /* d */,
    GC                  /* gc */,
    int                 /* x */,
    int                 /* y */,
    const char*         /* string */,
    int                 /* length */
#endif
);  

void win_draw_string(window_desc win, int x, int y, char *str, font_desc font, drawing_mode mode)
  {
    int n;

    if (!win || !str || !font) return;

    /* translate from local to global coords */
    x += WP(win)->origin_x; 
    y += WP(win)->origin_y;

    if (WP(win)->fid != FP(font)->fid)
      {
        XSetFont(dpy,WP(win)->gc,FP(font)->fid);
        WP(win)->fid = FP(font)->fid;
      }

    set_win_mode(win,mode);   /* on Dec3100, needed even w/ DrawImageString */

    n = strlen(str);
    if (!WP(win)->batch_mode)
      XDrawString(dpy,WP(win)->xwin,WP(win)->gc,x,y,str,n);
    XDrawString(dpy,WP(win)->pm,WP(win)->gc,x,y,str,n);

#ifdef DEBUG
    printf("Wstring: ");
    print_GC(WP(win)->gc);
#endif
  }

void mt_draw_text(MT_OBJECT *o, int x, int y, char *str)
{ window_desc win;
  bitmap_desc bm;
  font_desc font;
  int n;

  font = font_default();
  if (str == NULL) return;

  x += o->xl; 
  y += o->yl;

  if (o->type == BITMAP)
    { bm = (bitmap_desc) (o->wpwin);
      XSetFont(dpy,BP(bm)->gc,FP(font)->fid);

      set_bm_mode(bm,MODE_SRC);

#ifdef DEBUG
      printf("Bstring: ");
      print_GC(BP(bm)->gc);
#endif

      n = strlen(str);
      XDrawString(dpy,BP(bm)->pm,BP(bm)->gc,x,y,str,n);
    }
  else
    { win = o->wpwin;
      if (WP(win)->fid != FP(font)->fid)
        { XSetFont(dpy,WP(win)->gc,FP(font)->fid);
          WP(win)->fid = FP(font)->fid;
        }

      set_win_mode(win,MODE_SRC);

#ifdef DEBUG
      printf("Wstring: ");
      print_GC(WP(win)->gc);
#endif

      n = strlen(str);
      if (!WP(win)->batch_mode)
        XDrawString(dpy,WP(win)->xwin,WP(win)->gc,x,y,str,n);
      XDrawString(dpy,WP(win)->pm,WP(win)->gc,x,y,str,n);
    }
}

void mt_draw_title(MT_OBJECT *o, int x, int y, char *str)
{ window_desc win;
  bitmap_desc bm;
  static font_desc font = NULL;
  int n;

  if (font == NULL)
    font = font_open("*-helvetica-bold-r-*-*-14-*-*-*-*-*-*-*");
  if (str == NULL) return;

  x += o->xl; 
  y += o->yl;

  if (o->type == BITMAP)
    { bm = (bitmap_desc) (o->wpwin);
      XSetFont(dpy,BP(bm)->gc,FP(font)->fid);

      set_bm_mode(bm,MODE_SRC);

#ifdef DEBUG
      printf("Bstring: ");
      print_GC(BP(bm)->gc);
#endif

      n = strlen(str);
      XDrawString(dpy,BP(bm)->pm,BP(bm)->gc,x,y,str,n);
    }
  else
    { win = o->wpwin;
      if (WP(win)->fid != FP(font)->fid)
        { XSetFont(dpy,WP(win)->gc,FP(font)->fid);
          WP(win)->fid = FP(font)->fid;
        }

      set_win_mode(win,MODE_SRC);

#ifdef DEBUG
      printf("Wstring: ");
      print_GC(WP(win)->gc);
#endif

      n = strlen(str);
      if (!WP(win)->batch_mode)
        XDrawString(dpy,WP(win)->xwin,WP(win)->gc,x,y,str,n);
      XDrawString(dpy,WP(win)->pm,WP(win)->gc,x,y,str,n);
    }
}

/*---------------------------------------------------------------------------*/

void bm_draw_string(bm,x,y,str,font,mode)
  bitmap_desc bm;
  int x,y;
  char *str;
  font_desc font;
  drawing_mode mode;
  {
    int n;

    if (!bm || !str || !font) return;

    /* translate from local to global coords */
    x += BP(bm)->origin_x; 
    y += BP(bm)->origin_y;

    /*** xx need to cache font in bm as in windows (?) ****
    if (WP(win)->fid != FP(font)->fid)
      {
        XSetFont(dpy,WP(win)->gc,FP(font)->fid);
        WP(win)->fid = FP(font)->fid;
      }
    *******/
    XSetFont(dpy,BP(bm)->gc,FP(font)->fid);

    set_bm_mode(bm,mode);     /* on Dec3100, needed even w/ DrawImageString */

    n = strlen(str);
    XDrawString(dpy,BP(bm)->pm,BP(bm)->gc,x,y,str,n);

#ifdef DEBUG
    printf("Bstring: ");
    print_GC(BP(bm)->gc);
#endif
  }

/*---------------------------------------------------------------------------*/

bitmap_desc bm_new(w,h)
  int w,h;
  {
    bitmap_ptr bm;
    XGCValues gcval;

    bm = (bitmap_ptr)XtMalloc(sizeof(bitmap_struct));
    bm->pm = XCreatePixmap(dpy,root,w,h,depth);
    XFillRectangle(dpy,bm->pm,whtgc,0,0,w,h);
    bm->width = w;
    bm->height = h;
    bm->clip_x1 = 0;
    bm->clip_y1 = 0;
    bm->clip_x2 = w;
    bm->clip_y2 = h;
    bm->is_clipped = 0;
    bm->origin_x =  bm->origin_y = 0;
    gcval.foreground = black;
    gcval.background = white;
    gcval.cap_style  = CapButt;
    gcval.join_style = JoinMiter;
    bm->gc = XCreateGC(dpy,bm->pm,
                       GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
                       &gcval);
    return (bitmap_desc)bm;
  }

/*---------------------------------------------*/

bitmap_desc bm_read_file(filename)
  char *filename;
  {
    bitmap_ptr bm;
    int x, y, r;
    XGCValues gcval;

    bm = (bitmap_ptr)XtMalloc(sizeof(bitmap_struct));
    r = ReadBitmap(dpy,root,filename,&bm->width,&bm->height,&bm->pm,&x,&y);
    bm->clip_x1 = 0;
    bm->clip_y1 = 0;
    bm->clip_x2 = bm->width;
    bm->clip_y2 = bm->height;
    bm->is_clipped = 0;
    bm->origin_x =  bm->origin_y = 0;
    gcval.foreground = black;
    gcval.background = white;
    gcval.cap_style  = CapButt;
    gcval.join_style = JoinMiter;
    if (r == BitmapSuccess)
      {
        bm->gc = 
          XCreateGC(dpy,bm->pm,
                    GCForeground|GCBackground|GCCapStyle|GCJoinStyle,
                    &gcval);
        return (bitmap_desc)bm;
      }
    else
      return 0;
  }

/*---------------------------------------------*/

int bm_write_file(bm,filename)
  bitmap_desc bm;
  char *filename;
  {
    int r;

    if (!filename || !bm) return 0;

    r = XWriteBitmapFile(dpy,filename,
	    BP(bm)->pm, BP(bm)->width, BP(bm)->height, 0,0);

    return (r == BitmapSuccess);
  }

/*---------------------------------------------*/

void bm_free(bm)
  bitmap_desc bm;
  {
    XFreePixmap(dpy,BP(bm)->pm);
    XtFree((char*)bm);
  }

/*---------------------------------------------*/

bitmap_desc bm_build(data, w, h)
char *data;
int w, h;
{
   bitmap_ptr bm;
   XGCValues gcval;

   bm = (bitmap_ptr)XtMalloc(sizeof(bitmap_struct));
   bm->pm = XCreatePixmapFromBitmapData(dpy, root, data, w, h, black, white,
               depth);
   if(bm->pm == None) {
      XtFree((char *)bm);
      return 0;
   }
   bm->height = h;
   bm->width = w;
   gcval.foreground = black;
   gcval.background = white;
   bm->gc = XCreateGC(dpy,bm->pm,GCForeground|GCBackground,&gcval);
   bm->clip_x1    = 0;
   bm->clip_y1    = 0;
   bm->clip_x2    = 0;
   bm->clip_y2    = 0;
   bm->is_clipped = 0;
   bm->origin_x   = 0;
   bm->origin_y   = 0;

   return (bitmap_desc)bm;
}

/*---------------------------------------------*/

int bm_width(bm)
  bitmap_desc bm;
  {
    if (!bm) 
      return 0;
    else
      return (BP(bm)->width);
  }

/*---------------------------------------------*/

int bm_height(bm)
  bitmap_desc bm;
  {
    if (!bm) 
      return 0;
    else
      return (BP(bm)->height);
  }

/*---------------------------------------------*/

void bm_get_size(bm,w_result,h_result)
  bitmap_desc  bm;
  int         *w_result;
  int         *h_result;
  {
    if (!bm) return;
    if (w_result) *w_result = BP(bm)->width;
    if (h_result) *h_result = BP(bm)->height;
  }

/*---------------------------------------------*/

void bm_draw_pixel(bm,x,y,mode)
  bitmap_desc bm;
  int x,y,mode;
  {

    if (!bm) return;

    /* translate from local to global coords */
    x += BP(bm)->origin_x; 
    y += BP(bm)->origin_y;

    if (bm && x >= BP(bm)->clip_x1 && x <= BP(bm)->clip_x2 && 
	      y >= BP(bm)->clip_y1 && y <= BP(bm)->clip_y2) 
      { set_bm_mode(bm,mode);

        XDrawPoint(dpy,BP(bm)->pm,BP(bm)->gc,x,y);
      }

  }

/*---------------------------------------------*/

void win_draw_pixel(win,x,y,mode)
  window_desc win;
  int x,y,mode;
  {
    if (!win) return;

    /* translate from local to global coords */
    x += WP(win)->origin_x; 
    y += WP(win)->origin_y;

    if (win && x >= WP(win)->clip_x1 && x <= WP(win)->clip_x2 && 
	       y >= WP(win)->clip_y1 && y <= WP(win)->clip_y2) 
      { set_win_mode(win,mode);

        if (!WP(win)->batch_mode)
	  XDrawPoint(dpy,WP(win)->xwin,WP(win)->gc,x,y);
        XDrawPoint(dpy,WP(win)->pm,WP(win)->gc,x,y);
      }
  }

void mt_draw_pixel(MT_OBJECT *o, int x, int y)
{ window_desc win;
  bitmap_desc bm;

  x += o->xl; 
  y += o->yl;

  if (x >= o->xl && x <= o->xh && 
      y >= o->yl && y <= o->yh) 
    if (o->type == BITMAP)
      { bm = (bitmap_desc) (o->wpwin);
        set_bm_mode(bm,MODE_SRC);
        XDrawPoint(dpy,BP(bm)->pm,BP(bm)->gc,x,y);
      }
    else
      { win = o->wpwin;
        set_win_mode(win,MODE_SRC);
        if (!WP(win)->batch_mode)
          XDrawPoint(dpy,WP(win)->xwin,WP(win)->gc,x,y);
        XDrawPoint(dpy,WP(win)->pm,WP(win)->gc,x,y);
      }
}

/*---------------------------------------------*/

int bm_get_pixel(bm,x,y)
  bitmap_desc bm;
  int x,y;
  {
    XImage *im;
    int v;

    if (!bm) return -1;

    /* translate from local to global coords */
    x += BP(bm)->origin_x; 
    y += BP(bm)->origin_y;

    if (x < 0 || y < 0) return -1;
    if (x >= bm_width(bm) || y >= bm_height(bm)) return 0;

    im = XGetImage(dpy,BP(bm)->pm,x,y,1,1,AllPlanes,XYPixmap);
    v = XGetPixel(im,0,0);
    XDestroyImage(im);
    return (v == black);
  }

void bm_map_pixels(bitmap_desc bm, long old, long newa)
  {
    XImage *im;
    long v;
    int w, h;
    int x, y;

    w = bm_width(bm);
    h = bm_height(bm);

    im = XGetImage(dpy,BP(bm)->pm,0,0,w,h,AllPlanes,XYPixmap);
    for (x = 0; x < w; x++)
      for (y = 0; y < h; y++)
        { v = XGetPixel(im,x,y);
          if (v == old)
            XPutPixel(im,x,y,newa);
        }
    XPutImage(dpy,BP(bm)->pm,BP(bm)->gc,im,0,0,0,0,w,h);
    XDestroyImage(im);
    return;
  }

/*---------------------------------------------*/

int win_get_pixel(win,x,y)
  window_desc win;
  int x,y;
  {
    XImage *im;
    int v;

    if (!win) return -1;

    /* translate from local to global coords */
    x += WP(win)->origin_x; 
    y += WP(win)->origin_y;

    if (x < 0 || y < 0) return -1;
    if (x >= win_width(win) || y >= win_height(win)) return 0;

    im = XGetImage(dpy,WP(win)->pm,x,y,1,1,AllPlanes,XYPixmap);
    v = XGetPixel(im,0,0);
    XDestroyImage(im);
    return (v == black);
  }

/*---------------------------------------------*/

void bm_raster_op(dst,dx,dy,w,h,src,sx,sy,mode)
  bitmap_desc dst,src;
  int dx,dy,w,h,sx,sy;
  drawing_mode mode;
  {
    if (!src || !dst) return;

    /* translate from local to global coords */
    dx += BP(dst)->origin_x; 
    dy += BP(dst)->origin_y;
    sx += BP(src)->origin_x; 
    sy += BP(src)->origin_y;

    set_bm_mode(dst,mode);
    XCopyArea(dpy,BP(src)->pm,BP(dst)->pm,BP(dst)->gc,sx,sy,w,h,dx,dy);

#ifdef DEBUG
    printf("ToBBlt: ");
    print_GC(BP(dst)->gc);
#endif
  }

/*---------------------------------------------*/

void win_raster_op(dst,dx,dy,w,h,src,sx,sy,mode)
  window_desc dst;
  bitmap_desc src;
  int dx,dy,w,h,sx,sy;
  drawing_mode mode;
  {
    if (!src || !dst) return;

    /* translate from local to global coords */
    dx += WP(dst)->origin_x; 
    dy += WP(dst)->origin_y;
    sx += BP(src)->origin_x; 
    sy += BP(src)->origin_y;

    set_win_mode(dst,mode);
    XCopyArea(dpy,BP(src)->pm,WP(dst)->pm,WP(dst)->gc,sx,sy,w,h,dx,dy);

#ifdef DEBUG
    printf("ToWBlt: ");
    print_GC(WP(dst)->gc);
#endif

    if (!WP(dst)->batch_mode)
      XCopyArea(dpy,BP(src)->pm,WP(dst)->xwin,WP(dst)->gc,sx,sy,w,h,dx,dy);
  }

/*---------------------------------------------*/

void wp_win_saveblk(dst,dx,dy,w,h,src,sx,sy,mode)
  bitmap_desc dst;
  window_desc src;
  int dx,dy,w,h,sx,sy,mode;
  {
    if (!src || !dst) return;

    /* translate from local to global coords */
    dx += BP(dst)->origin_x; 
    dy += BP(dst)->origin_y;
    sx += WP(src)->origin_x; 
    sy += WP(src)->origin_y;

    set_bm_mode(dst,mode);
    XCopyArea(dpy,WP(src)->pm,BP(dst)->pm,BP(dst)->gc,sx,sy,w,h,dx,dy);
  }

/*---------------------------------------------*/

void bm_pattern_raster_op(dst,dx,dy,w,h,pat,sx,sy,mode)
  bitmap_desc dst,pat;
  int dx,dy,w,h,sx,sy;
  drawing_mode mode;
  {
    Pixmap pm;

    if (!dst || !pat) return;

    /* translate from local to global coords */
    dx += BP(dst)->origin_x; 
    dy += BP(dst)->origin_y;
    sx += BP(pat)->origin_x; 
    sy += BP(pat)->origin_y;

    pm = XCreatePixmap(dpy,root,w,h,depth);
    XSetFunction(dpy,tmpgc,GXcopy);
    XSetFillStyle(dpy,tmpgc,FillTiled);
    XSetTile(dpy,tmpgc,BP(pat)->pm);
    XSetTSOrigin(dpy,tmpgc,sx,sy);
    XFillRectangle(dpy,pm,tmpgc,0,0,w,h);
    XSetFillStyle(dpy,tmpgc,FillSolid);

    set_bm_mode(dst,mode);
    XCopyArea(dpy,pm,BP(dst)->pm,BP(dst)->gc,sx,sy,w,h,dx,dy);

    XFreePixmap(dpy,pm);

  }

/*---------------------------------------------*/

void win_pattern_raster_op(dst,dx,dy,w,h,pat,sx,sy,mode)
  window_desc dst;
  bitmap_desc pat;
  int dx,dy,w,h,sx,sy;
  drawing_mode mode;
  {
    Pixmap pm;

    if (!dst || !pat) return;

    /* translate from local to global coords */
    dx += WP(dst)->origin_x; 
    dy += WP(dst)->origin_y;
    sx += BP(pat)->origin_x; 
    sy += BP(pat)->origin_y;

    pm = XCreatePixmap(dpy,root,w,h,depth);
    XSetFunction(dpy,tmpgc,GXcopy);
    XSetFillStyle(dpy,tmpgc,FillTiled);
    XSetTile(dpy,tmpgc,BP(pat)->pm);
    XSetTSOrigin(dpy,tmpgc,sx,sy);
    XFillRectangle(dpy,pm,tmpgc,0,0,w,h);
    XSetFillStyle(dpy,tmpgc,FillSolid);

    set_win_mode(dst,mode);
    XCopyArea(dpy,pm,WP(dst)->pm,WP(dst)->gc,sx,sy,w,h,dx,dy);
    if (!WP(dst)->batch_mode)
      XCopyArea(dpy,pm,WP(dst)->xwin,WP(dst)->gc,sx,sy,w,h,dx,dy);

    XFreePixmap(dpy,pm);
  }

/*---------------------------------------------*/

void win_draw_bitmap(win,x,y,bm,mode)
  window_desc win;
  int x,y;
  bitmap_desc bm;
  drawing_mode mode;
  {
    if (!win || !bm) return;

    /* translate from local to global coords */
    x += WP(win)->origin_x; 
    y += WP(win)->origin_y;

    set_win_mode(win,mode);
    if (!WP(win)->batch_mode)
      XCopyArea(dpy,BP(bm)->pm,WP(win)->xwin,WP(win)->gc,
        0,0,BP(bm)->width,BP(bm)->height,x,y);
    XCopyArea(dpy,BP(bm)->pm,WP(win)->pm,WP(win)->gc,
      0,0,BP(bm)->width,BP(bm)->height,x,y);
  }

/*---------------------------------------------*/

void bm_draw_bitmap(on_bm,x,y,bm,mode)
  bitmap_desc on_bm;
  int x,y;
  bitmap_desc bm;
  drawing_mode mode;
  {
    if (!on_bm || !bm) return;

    /* translate from local to global coords */
    x += BP(bm)->origin_x; 
    y += BP(bm)->origin_y;

    set_bm_mode(on_bm,mode);
    XCopyArea(dpy,BP(bm)->pm,BP(on_bm)->pm,BP(on_bm)->gc,
      0,0,BP(bm)->width,BP(bm)->height,x,y);
  }

/*---------------------------------------------*/

void wp_win_fill(win,x1,y1,x2,y2,mode)
  window_desc win;
  int x1,y1, x2,y2, mode;
  {
    int t;

    if (!win) return;

    /* translate from local to global coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    if (x1 > x2) {t = x1; x1 = x2; x2 = t;}
    if (y1 > y2) {t = y1; y1 = y2; y2 = t;}

    set_win_mode(win,mode);

#ifdef DEBUG
    printf("Wfill: ");
    print_GC(WP(win)->gc);
#endif

    if (!WP(win)->batch_mode)
      XFillRectangle(dpy,WP(win)->xwin,WP(win)->gc, x1,y1,x2-x1,y2-y1);
    XFillRectangle(dpy,WP(win)->pm,WP(win)->gc, x1,y1,x2-x1,y2-y1);
  }

void mt_draw_rect(MT_OBJECT *o, int xl, int yl, int w, int h)
{ window_desc win;
  bitmap_desc bm;

  xl += o->xl;
  yl += o->yl;

  if (o->type == BITMAP)
    { bm = (bitmap_desc) (o->wpwin);

#ifdef DEBUG
      printf("Bfill: ");
      print_GC(BP(bm)->gc);
#endif

      set_bm_mode(bm,MODE_SRC);
      XFillRectangle(dpy,BP(bm)->pm,BP(bm)->gc, xl,yl,w,h);
    }
  else
    { win = o->wpwin;
      set_win_mode(win,MODE_SRC);

#ifdef DEBUG
      printf("Wfill: ");
      print_GC(WP(win)->gc);
#endif

      if (!WP(win)->batch_mode)
        XFillRectangle(dpy,WP(win)->xwin,WP(win)->gc, xl,yl,w,h);
      XFillRectangle(dpy,WP(win)->pm,WP(win)->gc, xl,yl,w,h);
    }
}

/*---------------------------------------------*/

void wp_btm_fill(bitmap_desc bm, int x1, int y1, int x2, int y2, int mode)
  {
    int t;

    if (!bm) return;

    /* translate from local to global coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    if (x1 > x2) {t = x1; x1 = x2; x2 = t;}
    if (y1 > y2) {t = y1; y1 = y2; y2 = t;}

    set_bm_mode(bm,mode);

#ifdef DEBUG
    printf("Bfill: ");
    print_GC(BP(bm)->gc);
#endif

    XFillRectangle(dpy,BP(bm)->pm,BP(bm)->gc, x1,y1,x2-x1,y2-y1);
  }

/*---------------------------------------------*/

void win_invert_rect(window_desc win, int x1, int y1, int x2, int y2)
  {
    int t;

    if (!win) return;

    /* translate from local to global coords */
    x1 += WP(win)->origin_x; y1 += WP(win)->origin_y;
    x2 += WP(win)->origin_x; y2 += WP(win)->origin_y;

    if (x1 > x2) {t = x1; x1 = x2; x2 = t;}
    if (y1 > y2) {t = y1; y1 = y2; y2 = t;}

    set_win_mode(win,MODE_XOR);

    if (!WP(win)->batch_mode)
      XFillRectangle(dpy,WP(win)->xwin,WP(win)->gc, x1,y1,x2-x1,y2-y1);
    XFillRectangle(dpy,WP(win)->pm,WP(win)->gc, x1,y1,x2-x1,y2-y1);
  }

/*---------------------------------------------*/

void bm_invert_rect(bm,x1,y1,x2,y2)
  bitmap_desc bm;
  int x1,y1, x2,y2;
  {
    int t;

    if (!bm) return;

    /* translate from local to global coords */
    x1 += BP(bm)->origin_x; y1 += BP(bm)->origin_y;
    x2 += BP(bm)->origin_x; y2 += BP(bm)->origin_y;

    if (x1 > x2) {t = x1; x1 = x2; x2 = t;}
    if (y1 > y2) {t = y1; y1 = y2; y2 = t;}

    set_bm_mode(bm,MODE_XOR);

    XFillRectangle(dpy,BP(bm)->pm,BP(bm)->gc, x1,y1,x2-x1,y2-y1);
  }

/*------------------------------------------------------------
 * fill_bitmap(bm,x,y,nv)
 *
 * one page seed fill program, 1 channel frame buffer version.
 *
 * (adapted from code by:)
 *       Paul Heckbert  13 Sept 1982, 28 Jan 1987
 * "doesn't read each pixel twice like the BASICFILL algorithm in
 *      Alvy Ray Smith, "Tint Fill", SIGGRAPH '79"
 */

/*
 * segment of scan line y for xl<=x<=xr was filled,
 * now explore adjacent pixels in scan line y+dy
 */
struct seg {short y, xl, xr, dy;};
#define SMAX 5000               /* max depth of stack (was 10000)*/

#define PUSH(Y, XL, XR, DY) \
    if (sp<stack+SMAX && Y+(DY)>=0 && Y+(DY)<h) \
    {sp->y = Y; sp->xl = XL; sp->xr = XR; sp->dy = DY; sp++;}

#define POP(Y, XL, XR, DY) \
    {sp--; Y = sp->y+(DY = sp->dy); XL = sp->xl; XR = sp->xr;}


void bm_fill(bm,x, y, nv)
  bitmap_desc bm;       /* pixrect to be filled */
  int x, y;             /* seed point */
  int nv;               /* new pixel value */
  {
    int l, x1, x2, dy;
    int ov;                     /* old pixel value */
    int w, h;   /* bounds */
    XImage *im;
    struct seg stack[SMAX], *sp = stack;
    int rx, ry, rw, rh; /* region to fill within */

    if (!bm) return;

    /* translate from local to global coords */
    x += BP(bm)->origin_x; 
    y += BP(bm)->origin_y;

    w = BP(bm)->width;
    h = BP(bm)->height;

    if (x<0 || x>=w || y<0 || y>=h) return;

    /* work out clipping region */
    if (BP(bm)->is_clipped)
      {
	if (BP(bm)->clip_x1 == BP(bm)->clip_x2 ||
	    BP(bm)->clip_y1 == BP(bm)->clip_y2) return;
        rx = BP(bm)->clip_x1;
        ry = BP(bm)->clip_y1;
        if (rx < 0) rx = 0;
        if (ry < 0) ry = 0;
	rw = BP(bm)->clip_x2 - rx;
	rh = BP(bm)->clip_y2 - ry;
	if (x < rx || y < ry) return;
	if (x >= rx+rw || y >= ry+rh) return;

	x -= rx;
	y -= ry;
      }
    else
      {
	rx = 0; ry = 0;
	rw = w; rh = h;
      }

    nv = nv ? black : white;
    im = XGetImage(dpy,BP(bm)->pm,rx,ry,rw,rh,AllPlanes,ZPixmap);

    ov = XGetPixel(im,x,y);             /* read old value at seed point */
    if (ov==nv) return;
    PUSH(y, x, x, 1);                   /* needed in some cases */
    PUSH(y+1, x, x, -1);                /* seed segment (popped 1st) */

    while (sp>stack) {
        /* pop segment off stack and fill a neighboring scan line */
        POP(y, x1, x2, dy);
        for (x=x1; x>=0 && XGetPixel(im,x,y)==ov; x--)
            XPutPixel(im,x,y,nv);
        if (x>=x1) goto skip;
        l = x+1;
        if (l<x1) PUSH(y, l, x1-1, -dy);                /* leak on left? */
        x = x1+1;
        do {
            for (; x<w && XGetPixel(im,x,y)==ov; x++)
                XPutPixel(im,x,y,nv);
            PUSH(y, l, x-1, dy);
            if (x>x2+1) PUSH(y, x2+1, x-1, -dy);        /* leak on right? */
skip:       for (x++; x<=x2 && XGetPixel(im,x,y)!=ov; x++);
            l = x;
        } while (x<=x2);
    }
    XPutImage(dpy,BP(bm)->pm,blkgc,im,0,0,rx,ry,w,h);
    XDestroyImage(im);
  }

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . -*/

void win_fill(win,x, y, nv)
  window_desc win;      /* window to be filled */
  int x, y;             /* seed point */
  int nv;               /* new pixel value */
  {
    int l, x1, x2, dy;
    int ov;                     /* old pixel value */
    int w, h;   /* bounds */
    XImage *im;
    struct seg stack[SMAX], *sp = stack;
    int cx1, cy1, cx2, cy2;
    int rx, ry, rw, rh; /* region to fill within */

    if (!win) return;

    /* translate from local to global coords */
    x += WP(win)->origin_x; 
    y += WP(win)->origin_y;

    w = win_width(win);
    h = win_height(win);

    if (x<0 || x>=w || y<0 || y>=h) return;

    /* work out clipping region */
    if (WP(win)->is_clipped)
      {
	/* pull out clip region */
	cx1 = WP(win)->clip_x1;
	cy1 = WP(win)->clip_y1;
	cx2 = WP(win)->clip_x2;
	cy2 = WP(win)->clip_y2;

	/* force it to be within bounds */
	if (cx1 < 0) cx1 = 0; if (cy1 < 0) cy1 = 0;
	if (cx2 >= w) cx2 = w-1; if (cy1 >= h) cy1 = h-1;

	/* bail out is everything clipped or seed point is outside */
	if (cx1 >= cx2 || cy1 >= cy2) return;
	if (x < cx1 || y < cy1 || x > cx2 || y < cy2) return;

	/* establish region to bring over from server */
        rx = cx1; 
        ry = cy1;
	rw = cx2 - rx;
	rh = cx2 - ry;

	/* adjust point to coodinates of region */
	x -= rx;
	y -= ry;
      }
    else
      {
	rx = 0; ry = 0;
	rw = w; rh = h;
      }

    nv = nv ? black : white;
    im = XGetImage(dpy,WP(win)->pm,rx,ry,rw,rh,AllPlanes,ZPixmap);

    ov = XGetPixel(im,x,y);             /* read old value at seed point */
    if (ov==nv) return;
    PUSH(y, x, x, 1);                   /* needed in some cases */
    PUSH(y+1, x, x, -1);                /* seed segment (popped 1st) */

    while (sp>stack) {
        /* pop segment off stack and fill a neighboring scan line */
        POP(y, x1, x2, dy);
        for (x=x1; x>=0 && XGetPixel(im,x,y)==ov; x--)
            XPutPixel(im,x,y,nv);
        if (x>=x1) goto skip;
        l = x+1;
        if (l<x1) PUSH(y, l, x1-1, -dy);                /* leak on left? */
        x = x1+1;
        do {
            for (; x<w && XGetPixel(im,x,y)==ov; x++)
                XPutPixel(im,x,y,nv);
            PUSH(y, l, x-1, dy);
            if (x>x2+1) PUSH(y, x2+1, x-1, -dy);        /* leak on right? */
skip:       for (x++; x<=x2 && XGetPixel(im,x,y)!=ov; x++);
            l = x;
        } while (x<=x2);
    }

    XPutImage(dpy,WP(win)->pm,blkgc,im,0,0,rx,ry,rw,rh);
    if (!WP(win)->batch_mode) 
      XPutImage(dpy,WP(win)->xwin,blkgc,im,0,0,rx,ry,rw,rh);
    XDestroyImage(im);
  }

#undef SMAX
#undef PUSH
#undef POP

/*---------------------------------------------------------------------------*/

/*
 *  The code following this point is mostly due to
 *  
 *      Mic Kaczmarczik
 *      User Services Unix Support Group
 *      UT Austin Computation Center
 *      Austin, Texas
 *  
 *      Internet:       mic@emx.utexas.edu
 *      UUCP:           ...!uunet!cs.utexas.edu!ut-emx!mic
 *      BITNET:         mic@utaivc
 *      THENET:         utaivc::mic
 *  
 *  This code reads X Bitmaps, Sun Icons, and Sun Rasterfiles interchangeably.
 *
 *  Arizona mods were to allow use with depth-not-one displays.
 */

/* Copyright, 1987, Massachusetts Institute of Technology */

/*
 *      Code to read bitmaps from disk files. Interprets 
 *      data from X10 and X11 bitmap files and creates
 *      Pixmap representations of files. Returns Pixmap
 *      ID and specifics about image.
 *
 *      Modified for speedup by Jim Becker, changed image
 *      data parsing logic (removed some fscanf()s). 
 *      Aug 5, 1988
 *
 *      Sun rasterfile and Suntools icon support by
 *      Mic Kaczmarczik, November 1988. The rasterfile and
 *      icon support can be used on machines other than Suns.
 *      In particular, it works fine on VAXen and Encores.
 *
 *      Uses code to decode run-length encoded Sun rasterfiles
 *      by Martin Boyer.
 *
 * Note that this file and ../Xmu/RdBitF.c look very similar....  Keep them
 * that way (but don't use common source code so that people can have one 
 * without the other).
 */

#define MAX_SIZE 255
#define Xfree free
#define Xmalloc malloc

static int _ReadBitmapFile(Display *, Drawable, FILE *,
                       int *, int *, Pixmap *, int *, int *);
static int _ReadRasterFile(Display *, Drawable, FILE *,
                       int *, int *, Pixmap *, int *, int *);
static int _ReadIconFile(Display *, Drawable, FILE *,
                     int *, int *, Pixmap *, int *, int *);
static int msb_int(unsigned char *);
static void pixtoX11(unsigned char *, int, int);
static unsigned int 
    decode_image(unsigned char *inpix, unsigned char *utpix, unsigned int lin);

/*
 * ReadBitmap is a front end for different routines to read
 * in bitmap files.  The order in which formats are tried is somewhat
 * arbitrary, but I figure X bitmaps should come first, then Suntools
 * icons (which will probably be used more often than big rasterfiles),
 * and finally rasterfiles.
 */

static
int ReadBitmap(display, d, filename, width, height, btm, x_hot, y_hot)
     Display *display;
     Drawable d;
     char *filename;
     int *width, *height;      /* RETURNED */
     Pixmap *btm;              /* RETURNED */
     int *x_hot, *y_hot;       /* RETURNED */
{
  int           status;
  FILE          *fstream = NULL /*, *fopen() */;

#define RETURN(code) { if (fstream) fclose(fstream); return (code); }

  if ((fstream = fopen(filename, "r")) == NULL)
          RETURN (BitmapOpenFailed);

  /* X11/10 bitmaps */
  status = _ReadBitmapFile(display, d, fstream,
                           width, height, btm, x_hot, y_hot);
  if (status != BitmapFileInvalid)
          RETURN (status);

  /* Suntools icons */
  rewind(fstream);
  status = _ReadIconFile(display, d, fstream, width, height,
                        btm, x_hot, y_hot);
  if (status != BitmapFileInvalid)
          RETURN (status);

  /* Sun rasterfiles */
  rewind(fstream);
  status = _ReadRasterFile(display, d, fstream,
                           width, height, btm, x_hot, y_hot);

  fclose(fstream);
  return (status);
}

/* shared data for the image read/parse logic */
static short hexTable[256];             /* conversion value */
static Bool initialized = False;        /* easier to fill in at run time */


/*
 *      Table index for the hex values. Initialized once, first time.
 *      Used for translation value or delimiter significance lookup.
 */
static void initHexTable(void)
{
    /*
     * We build the table at run time for several reasons:
     *
     *     1.  portable to non-ASCII machines.
     *     2.  still reentrant since we set the init flag after setting table.
     *     3.  easier to extend.
     *     4.  less prone to bugs.
     */
    hexTable['0'] = 0;  hexTable['1'] = 1;
    hexTable['2'] = 2;  hexTable['3'] = 3;
    hexTable['4'] = 4;  hexTable['5'] = 5;
    hexTable['6'] = 6;  hexTable['7'] = 7;
    hexTable['8'] = 8;  hexTable['9'] = 9;
    hexTable['A'] = 10; hexTable['B'] = 11;
    hexTable['C'] = 12; hexTable['D'] = 13;
    hexTable['E'] = 14; hexTable['F'] = 15;
    hexTable['a'] = 10; hexTable['b'] = 11;
    hexTable['c'] = 12; hexTable['d'] = 13;
    hexTable['e'] = 14; hexTable['f'] = 15;

    /* delimiters of significance are flagged w/ negative value */
    hexTable[' '] = -1; hexTable[','] = -1;
    hexTable['}'] = -1; hexTable['\n'] = -1;
    hexTable['\t'] = -1;
        
    initialized = True;
}

/*
 *      read next hex value in the input stream, return -1 if EOF
 */
static int NextInt(FILE *fstream)
{
    int ch;
    int value = 0;
    int gotone = 0;
    int done = 0;
    
    /* loop, accumulate hex value until find delimiter  */
    /* skip any initial delimiters found in read stream */

    while (!done) {
        ch = getc(fstream);
        if (ch == EOF) {
            value       = -1;
            done++;
        } else {
            /* trim high bits, check type and accumulate */
            ch &= 0xff;
            if (isascii(ch) && isxdigit(ch)) {
                value = (value << 4) + hexTable[ch];
                gotone++;
            } else if ((hexTable[ch]) < 0 && gotone)
              done++;
        }
    }
    return value;
}


/*
 * Read X bitmap files (this code is based on X11R3's XReadBitmapFile()).
 */

static int
_ReadBitmapFile(Display *display, Drawable d, FILE *fstream,
                int *width, int *height,
                Pixmap *pixmap, int *x_hot, int *y_hot)
{
    Pixmap pix;                         /* value to return */
    char *data = NULL;                  /* working variable */
    char line[MAX_SIZE];                /* input line from file */
    int size;                           /* number of bytes of data */
    char name_and_type[MAX_SIZE];       /* an input line */
    char *type;                         /* for parsing */
    int value;                          /* from an input line */
    int version10p;                     /* boolean, old format */
    int padding;                        /* to handle alignment */
    int bytes_per_line;                 /* per scanline of data */
    unsigned int ww = 0;                /* width */
    unsigned int hh = 0;                /* height */
    int hx = -1;                        /* x hotspot */
    int hy = -1;                        /* y hotspot */

    /* first time initialization */
    if (initialized == False) initHexTable();

    /* error cleanup and return macro   */
#ifdef RETURN
#undef RETURN
#endif
#define RETURN(code) { if (data) free (data); return code; }

    while (fgets(line, MAX_SIZE, fstream)) {
        if (strlen(line) == MAX_SIZE-1) {
            RETURN (BitmapFileInvalid);
        }
        if (sscanf(line,"#define %s %d",name_and_type,&value) == 2) {
            if (!(type = rindex(name_and_type, '_')))
              type = name_and_type;
            else
              type++;

            if (!strcmp("width", type))
              ww = (unsigned int) value;
            if (!strcmp("height", type))
              hh = (unsigned int) value;
            if (!strcmp("hot", type)) {
                if (type-- == name_and_type || type-- == name_and_type)
                  continue;
                if (!strcmp("x_hot", type))
                  hx = value;
                if (!strcmp("y_hot", type))
                  hy = value;
            }
            continue;
        }
    
        if (sscanf(line, "static short %s = {", name_and_type) == 1)
          version10p = 1;
        else if (sscanf(line,"static unsigned char %s = {",name_and_type) == 1)
          version10p = 0;
        else if (sscanf(line, "static char %s = {", name_and_type) == 1)
          version10p = 0;
        else
          continue;

        if (!(type = rindex(name_and_type, '_')))
          type = name_and_type;
        else
          type++;

        if (strcmp("bits[]", type))
          continue;
    
        if (!ww || !hh)
          RETURN (BitmapFileInvalid);

        if ((ww % 16) && ((ww % 16) < 9) && version10p)
          padding = 1;
        else
          padding = 0;

        bytes_per_line = (ww+7)/8 + padding;

        size = bytes_per_line * hh;
        data = (char *) Xmalloc ((unsigned int) size);
        if (!data) 
          RETURN (BitmapNoMemory);

        if (version10p) {
            char *ptr;
            int bytes;

            for (bytes=0, ptr=data; bytes<size; (bytes += 2)) {
                if ((value = NextInt(fstream)) < 0)
                  RETURN (BitmapFileInvalid);
                *(ptr++) = value;
                if (!padding || ((bytes+2) % bytes_per_line))
                  *(ptr++) = value >> 8;
            }
        } else {
            char *ptr;
            int bytes;

            for (bytes=0, ptr=data; bytes<size; bytes++, ptr++) {
                if ((value = NextInt(fstream)) < 0) 
                  RETURN (BitmapFileInvalid);
                *ptr=value;
            }
        }
    }                                   /* end while */

    if (data == NULL) {
        RETURN (BitmapFileInvalid);
    }

    pix = XCreatePixmapFromBitmapData (display,d,data,ww,hh,black,white,depth);
    if (pix == None) {
        RETURN (BitmapNoMemory);
    }
    *pixmap = pix;
    *width = ww;
    *height = hh;
    if (x_hot) *x_hot = hx;
    if (y_hot) *y_hot = hy;

    RETURN (BitmapSuccess);
}

/*
 * ReadRasterFile -- read a Sun rasterfile into a pixmap
 *
 * Read a Sun rasterfile into a pixmap. Calling sequence is the same as
 * XReadBitmapFile(). Run-length decoding routine courtesy of Martin Boyer,
 * posted to the Sun-Spots mailing list May 1988.  This version has been
 * tested on VAXen and Encores as well as Suns.  No chance to test it
 * on a Cray (yet).
 */

#define RASTERFILE_MAGIC        0x59a66a95
#define OLD                     0
#define STANDARD                1
#define BYTE_ENCODED            2

struct rasterfile_header {
        int     magic;
        int     width;
        int     height;
        int     depth;
        int     length;
        int     type;
        int     maptype;
        int     maplength;
};

static int
_ReadRasterFile(Display *display, Drawable d, FILE *fstream,
                int *width, int *height,
                Pixmap *btm, int *x_hot, int *y_hot)
{
  int hx = -1;
  int hy = -1;

  int   imagebytes, databytes, linewidth;
  char  *image = NULL, *data = NULL;
  struct rasterfile_header ras;

  Pixmap pix;

#ifdef  RETURN
#undef  RETURN
#endif
#define RETURN(code) { if (image) Xfree(image); if (data) Xfree(data); \
                               return (code); }

  /*
   * read in header and convert MSBFirst 4-byte longwords into integers
   */
  if (fread(&ras, sizeof(ras), 1, fstream) < 1) 
          RETURN (BitmapFileInvalid);

  ras.magic = msb_int((unsigned char *)&ras.magic);
  ras.width = msb_int((unsigned char *)&ras.width);
  ras.height = msb_int((unsigned char *)&ras.height);
  ras.depth = msb_int((unsigned char *)&ras.depth);
  ras.length = msb_int((unsigned char *)&ras.length);
  ras.type = msb_int((unsigned char *)&ras.type);
  ras.maptype = msb_int((unsigned char *)&ras.maptype);
  ras.maplength = msb_int((unsigned char *)&ras.maplength);

  /*
   * Check for unsupported rasterfile formats.   Since we don't handle color
   * images, ignore color map if provided.
   */

   if ((ras.magic != RASTERFILE_MAGIC) || (ras.depth != 1))
          RETURN (BitmapFileInvalid);

  if ((ras.type != OLD) && (ras.type != STANDARD) && (ras.type != BYTE_ENCODED))
          RETURN (BitmapFileInvalid);

  if (ras.maplength &&
      (0 == fseeko(fstream, (off_t) (sizeof(ras) + ras.maplength), 0)))
          RETURN (BitmapFileInvalid);

  /*
   * Determine the size of the image buffer and allocate memory for it.
   * If we're provided with a length for the image data, use it, else
   * the data is unencoded and uses the same amount of space as the
   * image buffer.
   */

  linewidth = ((ras.width + 15) / 16) * 2;
  imagebytes = linewidth * ras.height * ras.depth;
  databytes = ras.length ? ras.length : imagebytes;

  if ((image = (char *)Xmalloc(imagebytes)) == NULL)
          RETURN (BitmapNoMemory);

  /*
   * Read in the image.  If the format is BYTE_ENCODED, read into
   * a separate buffer and unpack into the image.
   */

  if (ras.type != BYTE_ENCODED) {
          if (fread(image, 1, imagebytes, fstream) != imagebytes)
                  RETURN (BitmapNoMemory);
  } else {
          if ((data = (char *)Xmalloc(databytes)) == NULL)
                  RETURN (BitmapNoMemory);
          if ((fread(data, 1, databytes, fstream) != databytes) ||
           (decode_image((unsigned char *)data,(unsigned char *)image,databytes)
                    != imagebytes))
                  RETURN (BitmapFileInvalid);
  }
          

  /*
   * Convert the bitmap data to X 11 format (basically unpad lines with
   * an odd line width), then create a pixmap from the data.  It might
   * be more efficient to query the server for the preferred bit order,
   * but I don't pretend to know the right way to do it.
   */

  pixtoX11((unsigned char *) image, ras.width, ras.height);
  pix = XCreatePixmapFromBitmapData (display, d, image,
                   ras.width, ras.height, black, white, depth);
  if (pix == None)
          RETURN (BitmapNoMemory);

  /*
   * Return the pixmap and its size, clean up and return.
   */

  *btm = pix;
  *width = ras.width;
  *height = ras.height;

  if (x_hot)
    *x_hot = hx;
  if (y_hot)
    *y_hot = hy;

  RETURN (BitmapSuccess);
}

/*
 * ReadIconFile -- read a Suntools icon into an X11 pixmap
 */

static int
_ReadIconFile(Display *display, Drawable d, FILE *fstream,
              int *width, int *height,
              Pixmap *btm, int *x_hot, int *y_hot)
{
  int hx = -1;
  int hy = -1;

  int   iconwidth = 64, iconheight = 64, valid_bits = 16, icondepth = 1;
  int   iconbytes, linewidth, nbytes, nitems, c1, c2;
  char  buf[BUFSIZ];
  register char *data = NULL, *cp;

  Pixmap pix;

#ifdef  RETURN
#undef  RETURN
#endif
#define RETURN(code) { if (data) Xfree(data); return (code); } 

  /*
   * The string "Format_version=1" signifies a Sun icon file. If found,
   * the rest of the line contains dimensions, etc.
   */
  for (;;) {
          if (fgets(buf, sizeof(buf), fstream) == NULL)
                  RETURN (BitmapFileInvalid);

          if (0 == strncmp("/* Format_version=1",buf,19))
                  break;
  }

  nitems = sscanf(buf,
"/* Format_version=1, Width=%d, Height=%d, Depth = %d, Valid_bits_per_item=%d",
                  &iconwidth, &iconheight, &icondepth, &valid_bits);

  if (nitems != 4)
          RETURN (BitmapFileInvalid);

  if ((icondepth != 1) || (valid_bits != 16))
          RETURN (BitmapFileInvalid);

  /*
   * seek to the end of the header section (end of the C comment)
   */
  while ((c1 = getc(fstream)) != EOF)
          if ((c1 == '*') && ((c2 = getc(fstream)) == '/'))
                  break;
  if ((c1 == EOF) || (c2 == EOF))
          RETURN (BitmapFileInvalid);

  /*
   * Determine the size of the pixrect-format image and allocate memory for it.
   * A Sun pixrect is composed of raster scan lines, each padded to 16-bit
   * word boundaries.  Pixels within each scan line are composed of 1, 4, or 8
   * bits, laid out in big-endian (e.g. most signifigant bit/byte first)
   * order.  In the case of single-plane images, we can ignore the depth
   * of the image in calculating the width of each scan line.
   */

  linewidth = ((iconwidth + 15) / 16) * 2;
  iconbytes = linewidth * iconheight;
  if ((data = (char *)malloc(iconbytes)) == NULL)
          RETURN (BitmapNoMemory);

  /*
   * Read in the image, reading each integer as an MSBFirst 2-byte short
   * Could of course be done faster.
   */
  nbytes = 0;
  cp = data; 
  while (nbytes < iconbytes) {
          if (fscanf(fstream," 0x%2x%2x", &c1, &c2) != 2)
                  RETURN (BitmapFileInvalid);
          *(cp++) = c1;
          *(cp++) = c2;
          if (((nbytes += 2) < iconbytes) && (fscanf(fstream,",") == EOF))
                  RETURN (BitmapFileInvalid);
  }

          
  /*
   * Convert the data into X 11 format and create a pixmap from it.
   * It might be more efficient to query the server for the preferred
   * bitmap bit order.
   */
  pixtoX11((unsigned char *) data, iconwidth, iconheight);
  pix = XCreatePixmapFromBitmapData (display,d,data,
                                iconwidth,iconheight,black,white,depth);
  if (pix == None)
          RETURN (BitmapNoMemory);


  /*
   * Return the pixmap and its size, clean up and return.
   */

  *btm = pix;
  *width = iconwidth;
  *height = iconheight;

  if (x_hot)
    *x_hot = hx;
  if (y_hot)
    *y_hot = hy;

  RETURN (BitmapSuccess);
}

/*
 * Date:    Mon, 9 May 88 17:17:07 EDT
 * From:    gamin%amadeus.UUCP@larry.mcrcim.mcgill.edu (Martin Boyer)
 * Subject: Re: Format for byte encoded rasterfiles (2)
 *
 * sow@cad.luth.se (Sven-Ove Westberg) and JDEBE@MTUS5.BITNET (John de
 * Beaubien) want to know the format of Suns runlength encoding for
 * rasterfiles.
 *
 * At one point, I wrote a quick filter to decode such files (I didn't want
 * the overhead of the pixrect library).  It took a while to find out (at the
 * time, there was a bug in pr_dump with RT_BYTE_ENCODED) but it paid off,
 * and I got something that is easily twice as fast as
 * /usr/lib/rasfilters/convert.2.
 *
 * The rules are the following:
 *
 * 1.  Encoding is byte per byte (regardless of the depth of the image)
 * 2.  Lengths of less than 3 are not encoded.
 * 3.  Lengths of 3 and more duplicate bytes are encoded as follows:
 *     for n consecutive "YZ" bytes, write out 0x80 (n-1) "YZ"
 *     so:  YZ YZ YZ YZ => 0x80 0x03 0xYZ
 *
 *     exception: 0x80 (alone) => 0x80 0x00    (note missing byte)
 *      and NOT 0x80 0x00 0x80
 *     but 0x80 0x80 => 0x80 0x01 0x80  (follows general rule)
 *
 * Even better, here is a function that does the decoding:
 * (I don't have any for encoding)
 *
 * It worked for at least a year on all sorts of images without a burp.
 * Good luck!
 *
 * Martin Boyer                               amadeus!gamin@mcgill-vision.uucp
 * Institut de recherche d'Hydro-Quebec       sun!sunlegende!amadeus!gamin
 * Varennes, QC, Canada   J0L 2P0             +1 514 652-8136
 */

/*
 * First argument is a pointer to the encoded array of pixels.
 *    Second is a pointer to enough space for the decoded array of pixels,
 *    which will be ras_width (rounded up to a short) * ras_height
 *      * ras_depth long
 *    Third is length of input (in bytes) (ras_length).
 *
 *    Returns length of decoded output, can be used to verify correct
 *    encoding/decoding.
 */

static
unsigned int 
decode_image(unsigned char *inpix, unsigned char *outpix, unsigned int lin)
{
    register unsigned char   value;
    register unsigned int    n;
    unsigned char  *outpix_0;

    outpix_0 = outpix;
    while (lin) {
        if ((value = *inpix++) == 0x80) {
            if ((n = *inpix++) == 0) {
                *outpix++ = 0x80;       /* special: 0x80 0x00 ==> 0x80 */
                lin -= 2;
            } else {
                for (value = *inpix++, n++; n; n--)
                    *outpix++ = value;
                lin -= 3;
            }
        } else {
            *outpix++ = value;
            lin--;
        }
    }
    return (outpix - outpix_0);
}

/*
 * Treat cp[0] through cp[3] as a big-endian 32-bit longword and convert
 * it into an integer.
 */

static int msb_int(unsigned char *cp)
{
        return ((cp[0] << 24) | (cp[1] << 16) | (cp[2] << 8) | cp[3]);
}

/*
 * 256-byte table for quickly reversing the bits in an unsigned 8-bit char,
 * used to convert between MSBFirst and LSBFirst image formats.
 */

static char revtable[256] = { 
            0, -128,   64,  -64,   32,  -96,   96,  -32,
           16, -112,   80,  -48,   48,  -80,  112,  -16,
            8, -120,   72,  -56,   40,  -88,  104,  -24,
           24, -104,   88,  -40,   56,  -72,  120,   -8,
            4, -124,   68,  -60,   36,  -92,  100,  -28,
           20, -108,   84,  -44,   52,  -76,  116,  -12,
           12, -116,   76,  -52,   44,  -84,  108,  -20,
           28, -100,   92,  -36,   60,  -68,  124,   -4,
            2, -126,   66,  -62,   34,  -94,   98,  -30,
           18, -110,   82,  -46,   50,  -78,  114,  -14,
           10, -118,   74,  -54,   42,  -86,  106,  -22,
           26, -102,   90,  -38,   58,  -70,  122,   -6,
            6, -122,   70,  -58,   38,  -90,  102,  -26,
           22, -106,   86,  -42,   54,  -74,  118,  -10,
           14, -114,   78,  -50,   46,  -82,  110,  -18,
           30,  -98,   94,  -34,   62,  -66,  126,   -2,
            1, -127,   65,  -63,   33,  -95,   97,  -31,
           17, -111,   81,  -47,   49,  -79,  113,  -15,
            9, -119,   73,  -55,   41,  -87,  105,  -23,
           25, -103,   89,  -39,   57,  -71,  121,   -7,
            5, -123,   69,  -59,   37,  -91,  101,  -27,
           21, -107,   85,  -43,   53,  -75,  117,  -11,
           13, -115,   77,  -51,   45,  -83,  109,  -19,
           29,  -99,   93,  -35,   61,  -67,  125,   -3,
            3, -125,   67,  -61,   35,  -93,   99,  -29,
           19, -109,   83,  -45,   51,  -77,  115,  -13,
           11, -117,   75,  -53,   43,  -85,  107,  -21,
           27, -101,   91,  -37,   59,  -69,  123,   -5,
            7, -121,   71,  -57,   39,  -89,  103,  -25,
           23, -105,   87,  -41,   55,  -73,  119,   -9,
           15, -113,   79,  -49,   47,  -81,  111,  -17,
           31,  -97,   95,  -33,   63,  -65,  127,   -1,
};

/*
 * Convert pixrect-format data into X 11 format data, which means reversing
 * the bits, and removing the extra byte of padding.  Removing the padding
 * is done by backing up at the end of a pixrect line, if necessary.
 */

static void pixtoX11(unsigned char *data, int width, int height)
{
        register unsigned char *cp, *end, *datap;
        int row, linewidth, unpad = False;

        linewidth = ((width + 7) / 8);
        if (linewidth & 1) {
                linewidth++;
                unpad = True;
        }

        datap = data;
        for (row = 0 ; row < height ; row++) {
                cp = data + (row * linewidth);
                end = cp + linewidth;
                while (cp < end)
                        *(datap++) = revtable[*(cp++)];
                if (unpad)
                        datap--;        /* back up to unpad line        */
        }
}

/*======================================================================*/

void wp_user_event(long time, window_desc win)
{ XMotionEvent event;
  XEvent query;
  long status;
  int  wx, wy;

  while (1)
    { if (XPending(dpy) == 0) break;
      XPeekEvent(dpy,&query);
      if (query.type != MotionNotify) break;
      XNextEvent(dpy,&query);
    }

  { int xa, xb;
    unsigned int mk;
    Window wa, wb;

    XQueryPointer(dpy,WP(win)->xwin,&wa,&wb,&xa,&xb,&wx,&wy,&mk);
  }

  event.type = MotionNotify;
  event.send_event = 1;
  event.display = dpy;
  event.window = WP(win)->xwin;
  event.root = root;
  event.time = time;
  event.x = wx;
  event.y = wy;
  event.x_root = wx;
  event.y_root = wy;
  event.state = 0;
  event.is_hint = 0;
  event.same_screen = 1;
  status = XSendEvent(dpy,WP(win)->xwin,0,0,(XEvent *) (&event));
}

font_desc wp_get_font(char *face, char *wgt, char *slant, char *size)
{ char name[1000];
  sprintf(name,"-*-%s-%s-%s-*-*-%s-*-*-*-*-*-*-*",face,wgt,slant,size);
  return (font_open(name));
}

long mt_get_color(r,g,b) int r, g, b;
{ static XColor fill;
  long c;

  fill.pixel = 0;
  fill.red   = (r << 8);
  fill.green = (g << 8);
  fill.blue  = (b << 8);
  fill.flags = (DoRed | DoGreen | DoBlue);
  XAllocColor(dpy,their_cmap,&fill);
  c = fill.pixel;
  return (c);
}

void mt_set_color(long c)
{ current_color = c; }

long mt_black()
{ return (black); }

long mt_white()
{ return (white); }
