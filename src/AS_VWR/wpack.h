
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
/* $Id: wpack.h,v 1.1.1.1 2004-04-14 13:54:20 catmandew Exp $ */

/*
  Minimal Retargetable Window Package

  Introduction
    This header file defines data structures and routines for "winpack"
    a minimalist retargetable window package developed at the University
    of Arizona.  Winpack is designed to provide a minimal set of window 
    facilities in an easy to learn and use package.  Winpack is supported 
    "on top of" an existing "host" window system (e.g., the X window system) 
    and is designed to be easy to retarget to new host window systems.

    more later...

  Overview of Major Constructs
    As described below, the system supports 5 basic types of objects: 
    windows, bitmaps, fonts, cursors, and events.  Windows and bitmaps 
    are the primary objects of the system.  Drawing can take place on 
    either windows or bitmaps uniformly.  Cursors and fonts are system
    resources that define the image that tracks the pointing device and
    the appearance of textual output respectively.  Finally, events are a 
    record of significant input actions.  Events are queued in the order 
    of their occurrence and handled later by means of an application callback 
    routine attached to each window.
    
    The name of winpack each routine contains a prefix which indicates which 
    type of object it deals with.  Window routines begin with "win_", 
    bitmap routines begin with "bm_", cursor routines begin with "cursor_", 
    and finally, font routines begin with "font_".  In addition, the general 
    routines needed to initialize, run, and terminate the package as a whole 
    begin with "winpack_" and routines to manipulate events begin with "event_".
    
    Windows
      Winpack windows provide a set of (potentially) overlapping drawing 
      surfaces on the workstation screen.  Winpack maintains backing 
      store for all windows, consequently, all overlap relationships
      between windows are hidden from the application.  As a result, each 
      window can be treated as an abstract drawing surface without regard to
      overlap or visibility considerations.  For example, output which is 
      obscured by overlap or iconfication of a window appears automatically 
      when window relationships change or the window is deiconified.  Winpack 
      windows may be explicitly positioned and resized both under program and 
      user control.

    Bitmaps
      Bitmaps provide a drawing surface that does not directly appear on
      the screen (to make the image of a bitmap appear on the screen, it
      must be drawn in some window).   Both bitmaps and windows provide a
      rectangular drawing surface of black and white (1 bit) pixels.  
      Drawing operations on bitmaps and windows include: lines, "raster-op"
      or "bit-blit" operations, text in various fonts,  setting and clearing 
      of individual pixels, and flood fills.  In addition, both bitmaps and 
      windows support clipping to an arbitrary rectangle.
      
    Fonts
      Fonts provide a means to describe text output within windows or bitmaps. 
      The details of font naming and exact appearance depend on the underlying
      host window system.  The X-utilities "xfontsel" and "xlsfont" allow
      you to explore the set of possible fonts available and their naming
      conventions (for use in "open fonts").  xfontsel interactively allows
      you to browse the available fonts and preview their appearance.
      xlsfont gives a list of all fonts available on the system.  Font
      names area a sequence of characteristics separated by dashes, e.g.,
         -adobe-courier-bold-o-normal-*-10-100-75-75-m-0-iso8859-1
      Stars are don't care indicators.

    Cursors
      Cursors can be used to define the image that tracks the pointing device
      within a particular window on the screen.  Cursors can come from a 
      predefined set supported by the underlying host window system and/or may 
      be defined by the application.

    Events
      All inputs within the winpack system are modeled as events.  An
      event record encodes the "what", "when", and "where" of each significant 
      input action.  Winpack works as an event driven system. As events occur, 
      they are queued by the system and passed to an application "callback" 
      routine to be handled.  Note that winpack controls the "main loop" of a 
      program.  The application program plays a passive role responding
      to events as they occur rather than directly managing control flow.  
  
  Output Model
    Graphical output may be created uniformly on either windows or bitmaps.
    Both drawing surfaces provide a rectangular array of black and white
    (1 bit) pixels.  Pixels with the integer value 1 appear as black whereas
    pixels with the integer value 0 appear as white.  All output operations 
    are defined in terms of an integer coordinate system with its origin at 
    the top-left.  For example, in a 64x64 pixel window or bitmap the top-left 
    pixel is found at <0,0> while the bottom right pixel is located at <63,63>.

    The following primitive output operations are provided (actual routine 
    names are prefixed with "win_" or "bm_"):

      draw_pixel    - Set the value of one pixel.
      draw_line     - Draw a 1 pixel wide line between 2 given points.
      draw_string   - Draw a text string at a given point in a given font.
      raster_op     - Perform a bitwise combination of 2 rectangular areas.
      fill          - Perform a flood fill operation.

    In addition to these primitive operations, a number of commonly used 
    shortcut operations are supported.  (Each of these operations can 
    expressed in terms of a more complicated raster_op operation).   
    These include:

      pattern_raster_op   - Perform a raster_op using a repeated pattern.
      draw_bitmap         - Perform a raster_op using a full bitmap.
      clear_rect          - Clear a rectangular region of pixels to white.
      invert_rect         - Invert a rectangular region of pixels.

    Each output operation (except fill)  combines a "source" object (a line, 
    text string image, or rectangular region of pixels)  with a "destination" 
    object (a window or bitmap).  This combination is done in a bitwise 
    fashion.  Each bit (pixel) of the source object is logically combined 
    with each bit (pixel) of the destination using a "drawing mode" function.  
    For example, if drawing is done in "or" mode the following results are
    obtained:

      Source     Destination     Result
      Value      Value           Value
      -------------------------------------------
	0           0        |     0
	0           1        |     1
	1           0        |     1
	1           1        |     1

    while drawing in "and" mode would produce the following results:

      Source     Destination     Result
      Value      Value           Value
      -------------------------------------------
	0           0        |     0
	0           1        |     0
	1           0        |     0
	1           1        |     1

    In general, any of the 16 possible logical operators may be used as 
    a drawing function.  To facilitate use of these functions, special
    symbolic constants have been defined: MODE_SRC and MODE_DST.  These
    values can be combined using the C bitwise operators (& | ~ ^) to create
    the constant needed to specify each of the drawing modes.  For example 
    the "or" mode drawing function is denoted by "MODE_SRC | MODE_DST".  
    The following table gives a synopsis of some of the more useful drawing
    modes:

      MODE_SRC             - copy source pixels over destination.
      MODE_SRC | MODE_DST  - add new black pixels from source ("or" mode).
      MODE_SRC & MODE_DST  - leave black pixels only where there is black 
			     in source ("mask" operation, "and" mode).
      ~MODE_DST            - invert destination pixels.
      MODE_SRC ^ MODE_DST  - invert destination pixels only where there is black
			     in the source.  note: two xor's of the same image
			     leave the destination in its original state.
      ~MODE_SRC & MODE_DST - clear areas in destination that are black in 
			     source.

    In addition, several shortcut constants are defined:

      MODE_PAINT           - same as "or" mode.
      MODE_XOR             - same as "xor" mode.
      MODE_INVERT          - same as ~MODE_DST
      MODE_SET             - set destination pixels to black unconditionally
      MODE_CLEAR           - set destination pixels to white unconditionally

    Note that when lines and text are being draw, only the pixels of the lines
    or characters themselves are part of the source.  These source pixels are 
    always black.  

    xx need to talk about clipping and batching

  Input Model
    All input accepted by the system is represented by a uniform event record
    data structure.  Each significant input action causes an event record to be 
    placed on an internal queue of such events.  During normal processing 
    (e.g., inside winpack_main_loop), the system loops: waiting for a new 
    event, passing that event on to the application by calling a "callback" 
    routine, and then returning to wait for the next event.

    Each event records facts about one significant input action.  This record 
    includes information about:
      * what kind of event occurred.
      * where the locator was pointing when the event occurred
        (which window and the pixel coordinate within the window).
      * when the event occurred.
      * what the status of the modifier keys (shift, control, etc.) was 
	when the event occurred.
      * specific values associated with the particular event type (e.g, button
	number, ASCII code of keypress, etc.)

    The system currently supports 12 types of events:
      loc_but_dn    - a button on the locator (mouse) went down
      loc_but_up    - a button on the locator went up
      but_dn        - some other button went down
      but_up        - some other button went up
      loc_move      - the locator moved
      val_move      - a valuator (knob or slider) moved
      key_dn        - a keyboard key went down
      key_up        - a keyboard key went up
      win_enter     - the cursor entered the bounds of a window
      win_exit      - the cursor exited the bounds of a window
      win_resize    - a window changed size (under user control)
      special_event - "hook" for device dependent extensions

    In addition to the type of input action, each event record also encodes 
    the location of the locator at the time of the action.  This includes 
    the window over which the locator was positioned, as well as the position 
    within that window.  Positions are reported in the local coordinate
    system of the window they appear in (i.e., relative to the top-left 
    corner of the window).  Note: as an exception to the general rule, 
    win_resize events record the new size of the window rather than the 
    locator position.

    Input event records also encode the time at which events occur.  Event
    times are recorded as a number of seconds and microseconds past some
    fixed previous point in time.  These times can be used to determine the 
    interval between two events but may not be suitable for determining wall 
    clock time.  Note: accuracy and resolution of timing varies widely and is 
    never accurate to the nearest microsecond.

    Another component of each event record is an encoding of the status of 
    several modifier keys at the time of the event.  This encoding consists
    of an integer value formed by the bitwise or-ing of any of the values:
      SHIFT_KEY_DOWN  - the shift key was being held down
      CNTRL_KEY_DOWN  - the control key was being held down
      META_KEY_DOWN   - another alternate modifier key was being held down

    In order to support multiple devices of the same type (i.e., 2 or more 
    valuators), events also provide a device number to indicate which
    physical device the event came from.  Note: currently, multiple devices
    of the same type are not supported.

    The final component of an event record is a type specific value.  This 
    value is always a long integer value, but its interpretation varies
    depending on event type as indicated below:
      key up/dn      - ASCII code of the key
      button up/dn   - device dependent button number 
      loc_but up/dn  - locator button number
      val_move       - new valuator value
      loc_move       - ignored
      win enter/exit - ignored
      win_resize     - ignored

    xx need to talk about event dispatch and callback

  Synopsis of Data Structures
    Except for events, objects of each of the types supported by the system 
    are represented by descriptors.  These descriptors are implemented as 
    "opaque" pointers -- pointers to hidden machine dependent data structures 
    that cannot be accessed directly.  The C type names for these descriptor
    include: window_desc, bitmap_desc, cursor_desc, and font_desc.  Input
    events are implemented by the input_event type.

  Synopsis of Routines
    General Routines
      winpack_init       - Initialize the package.  Must be called before any 
			   other routines.
      winpack_shutdown   - Close all windows, release all resources, and exit
			   from winpack_main_loop.  Once winpack has been shut 
			   down, it cannot be restarted.
      winpack_main_loop  - Enter the main event dispatch loop. Does not exit 
			   until winpack_shutdown is called.
  
    Event Manipulation Routines
      event_in_mask      - Determine if the given event code is contained in the
			   given set of events (event "mask").
      event_code_string  - Produce a human readable string corresponding to an
			   event type code.
      event_string       - Produce a human readable string describing an event.
  
    Object Creation & Reclamation Routines
      win_new            - Create a new window with a given size, global screen 
			   position, and event handling callback procedure.  
			   Windows are initially created "invisible" (i.e., 
			   not mapped to the screen or appearing on the screen).
      win_free           - Destroy a window and reclaim any resources associated
			   with it.
      bm_new             - Create a new (all white) bitmap of a given size.
      bm_build           - Create and initialize "static" bitmap from a data 
			   byte array.
      bm_read_file       - Create a new bitmap and read its image from the named
			   file.  Bitmap file formats supported include at least
			   the one produced by bm_write_file.  (Currently this
			   is X-windows bitmap format.  SunView bitmap format
			   is also currently supported.)
      bm_free            - Reclaim any resources associated with a bitmap.  
      cursor_open        - Create or open a cursor.  Cursors can come from a 
			   named standard set, or may be created from a pair of 
			   bitmap files (one for the cursor image and one for 
			   the mask).
      cursor_free        - Reclaim any resources associated with a cursor.
      font_open          - Open a named font.  Font names are host window 
			   system dependent.
      font_default       - Return the default font.
      font_free          - Reclaim any resources associated with a font.
      
    Common Window and Bitmap Query Functions
      *_width            - Return the width of a window or bitmap.
      *_height           - Return the height of a window or bitmap.
      *_get_size         - Return both width and height.
      *_get_pixel        - Return the pixel value at a given coordinate.
  
    Common Window and Bitmap Clipping Manipulation Routines
      *_set_clip         - Set the clipping rectangle
      *_get_clip         - Retrieve the current clipping Rectangle
      *_reset_clip       - Reset the clipping rectangle to the whole window or 
			   bitmap
      *_reduce_clip      - Set the clipping rectangle, but make sure it stays
			   within the current one.
      *_trivial_reject   - Given a bounding box, quickly test if all drawing 
			   within the area of that box would be entirely 
			   clipped away and hence need not be drawn.
  
    Common Window and Bitmap Drawing Routines
      *_draw_line        - Draw a 1 pixel wide line between 2 points.
      *_draw_string      - Draw a text string at a given point in a given font.
      *_draw_pixel       - Set the value of a single pixel.
      *_raster_op        - Perform a bitwise logical combination between 2
			   rectangular regions of pixels.
      *_pattern_raster_op - Perform a raster_op operation where the source 
			   bitmap is taken to be an infinite bitmap constructed 
			   from tiling a pattern bitmap.
      *_draw_bitmap      - Perform a simple copy raster_op from the full extent
			   of a given bitmap.
      *_fill             - Perform an 8-connected region fill starting from a 
			   given seed point.  Fills can be black (filling white
			   regions) or white (filling black regions).
      *_clear_rect       - Clear a rectangular area to white.
      *_invert_rect      - Invert a rectangular area.

    Window Redraw Batching Control Routines
      win_start_batch    - Start batching updates (output does not appear on the
			   screen until win_end_batch is called).
      win_end_batch      - End a batch and force any deferred output to appear 
                           on the screen.  Note: batches do no nest.

    Specialized Window Property Manipulation Routines
      win_set_pos        - Move a window.
      win_get_pos        - Retrieve the global screen position of a window
      win_set_size       - Change the size of a window.
      win_set_vis        - Change whether a window appears on the screen or not.
      win_get_vis        - Inquire whether a window is eligible to appear on the
			   screen.  Note: some or all of the image of a visible 
			   window may not actually appear on the screen if it 
			   is obscured by other windows or iconified by the 
			   user.
      win_set_handler    - Change the event callback routine (event handler) for
			   a given window.
      win_get_handler    - Retrieve the current callback routine for a window.
      win_set_data       - Each window maintains 32 bits of private user data,
			   set the value associated with a given window.
      win_get_data       - Retrieve the private user data associated with a 
			   window.
      win_set_cursor     - Set the cursor to be used when the locator is within
			   a given window.
  
    Specialized Bitmap Operations
      bm_write_file      - Write a bitmap file.
  
    Font Operations
      font_string_size   - Determine the drawing size of a given string to be 
			   drawn in a given font.
  
  Use
    more later ...

  Limitations
    In keeping with its minimalist foundations, the winpack package does 
    not support the following features you might be expecting to see in a 
    window package:
      * color or greyscale
      * lines wider than 1 pixel
      * circles, arcs, and curves
      * support for cut and paste or other "interclient communication"
      * a high level imaging model (e.g., PostScript)

  Contributors
    Winpack was first developed under SunView and the X window systems by
    Scott Hudson.  The second X port (which is the basis of the current 
    version) was done by Gregg Townsend [gmt].  Additional work to support
    cursors, "static" bitmaps, and window resize events under X was done 
    by Shamim Mohamed [sham].

  Revisions
    8/88  - Rework for version 2.0
    11/88 - Port to X-windows
    11/88 - Add additional features
    2/89  - Fix bugs in clipping
    2/89  - Add fonts
    3/89  - Attach windows to events
    4/89  - Numerous updates and new X-windows port [gmt]
    5/89  - Add event masks
    6/89  - Small name changes
    1/90  - Add fix for broken time.h under g++
          - Add invert_rect()
          - Add titles to windows
          - Add set_win_clip() and get_win_clip()
    1/90  - Add trivial_reject()
    5/90  - Added cursor support
            Added window-resize events
            Added static bitmaps
            [sham]
    6/90  - Added write_bitmap_file()
    7/90  - Version 3.0
	  - Reorganized for easier interface to C++ 
	  - Unified drawing behavior of bitmaps and windows
	  - Introduced new naming scheme and #defines for backward compat
	  - Added event to string conversion routines
	  - Various bug fixes
    7/90  - Misc small changes for ver 3.1
    7/90  - Added translation routines ( *_*_origin() )
    12/90 - Added global_x and global_y to events;
    12/90 - Added text_bitmaps and network routines [sham]

  Notes 
    At different times this header file is compiled as "old-fashioned" C, 
    ANSI C, and C++ and must be compatible with all three.

    Bug reports should be sent to Scott Hudson (hudson@cs.arizona.edu).

  Known Problems
    Drawing modes should be defined as an enum
    Fills under X-windows don't work with clipping correctly

  Things left to do in current release (xx):
    Finish top documentation 
    Improve documentation of individual routines
*/
/*----------------------------------------------------------------------*/
/*======================================================================*/
/* Configuration */
/*======================================================================*/

typedef void *pointer;
#define Params(a) a

/*======================================================================*/
/* Pointers to Major Data Types */
/*======================================================================*/

/* "Opaque" pointers to machine dependent structures */

typedef pointer bitmap_desc;
typedef pointer window_desc;
typedef pointer font_desc;
typedef pointer cursor_desc;
typedef pointer opaque_pointer;

/*======================================================================*/
/* Overall Window Package Control */
/*======================================================================*/

/* Start up and shut down the window system */
void            winpack_init      Params((void));
void            winpack_shutdown Params((void));

/* Main event dispatch loop */

void            winpack_main_loop Params((void));

/* Network support - add fd to the select */
typedef void (*file_callback_proc) Params((int fd));
void            winpack_monitor_fd   Params((int fd, file_callback_proc uprc));
void            winpack_unmonitor_fd Params((int fd));

/*======================================================================*/
/* Input Events */
/*======================================================================*/

/* Codes specifying what kind of event something is */

typedef enum INPUT_CODE {
  EVENT_none    = 0x00000000,
  loc_but_dn    = 0x00000001, /* locator (mouse) button went down       */
  loc_but_up    = 0x00000002, /* locator (mouse) button went up         */
  but_dn        = 0x00000004, /* some other button went down            */
  but_up        = 0x00000008, /* some other button went up              */
  loc_move      = 0x00000010, /* locator (mouse) has moved              */
  val_move      = 0x00000020, /* valuator has moved                     */
  key_dn        = 0x00000040, /* keyboard key has went down             */
  key_up        = 0x00000080, /* keyboard key has went up               */
  win_enter     = 0x00000100, /* cursor has entered a window            */
  win_exit      = 0x00000200, /* cursor has exited a window             */
  win_resize    = 0x00000400, /* the underlying window has been resized */
  special_event = 0x00000800, /* code for device dependent extensions   */
  EVENT_all     = 0x0fffffff
} input_code, input_mask;

#define event_in_mask(mask,code)   ((code)&(mask))

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Masks for modifier keys (later this should be an enum) */
typedef enum KEY_STATUS {
  KEY_none       = 0,
  SHIFT_KEY_DOWN = 0x0001,
  CNTRL_KEY_DOWN = 0x0002,
  META_KEY_DOWN  = 0x0004
} key_status;

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Record of an input event */

typedef struct INPUT_EVENT {
  input_code    code;       /* what was the event                            */
  long          x,y;        /* where was the pointing device at the time     */
  long          global_x;   /* a second copy that will always stay on the    */
  long          global_y;   /* window coordinate system                      */
  window_desc   in_window;  /* which window was the pointing device in       */
  long          when;       /* when did the event occur                      */
  key_status    updown_mask;/* what modifier keys were in effect             */
  short         dev_num;    /* which device event is from                    */
  long          value;      /* value associated with event                    
                                values for various kinds of events are:
                                  key_up/dn      ASCII code for key
                                  but_up/dn      button number
                                  loc_but_up/dn  button number
                                  val_move       new value
                                  special_event  depends on dev_num
                            */
} input_event;

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Debugging routines to provide human readable description of event */ 
/*   note: return values point to static storage and are overwritten */
/*         by subsequent calls.                                      */

char *event_code_string Params((input_code cd));
char *event_string      Params((input_event *evt));

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Type for event handler procedures */

typedef void (*event_handler_proc)
  Params((input_event *,window_desc,opaque_pointer));


/*======================================================================*/
/* Drawing Routines (common to bitmaps and windows) */
/*======================================================================*/

/* Clipping control */

void win_set_clip       Params((window_desc win, int x1,int y1, int x2,int y2));
void bm_set_clip        Params((bitmap_desc bm, int x1,int y1, int x2,int y2));

void win_get_clip       Params((window_desc win, 
                                int* x1_result, int *y1_result, 
				int* x2_result, int* y2_result));
void bm_get_clip        Params((bitmap_desc bm, 
                                int* x1_result, int *y1_result, 
				int* x2_result, int* y2_result));

void win_reset_clip     Params((window_desc win));
void bm_reset_clip      Params((bitmap_desc bm));

void win_reduce_clip    Params((window_desc win, int x1,int y1, int x2,int y2));
void bm_reduce_clip     Params((bitmap_desc bm, int x1,int y1, int x2,int y2));

int  win_trivial_reject Params((window_desc win, int x1,int y1, int x2,int y2));
int  bm_trivial_reject  Params((bitmap_desc bm, int x1,int y1, int x2,int y2));

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Translation (origin control) routines */

void win_reset_origin   Params((window_desc win));
void bm_reset_origin    Params((bitmap_desc win));

void win_set_origin     Params((window_desc win, int x, int y));
void bm_set_origin      Params((bitmap_desc win, int x, int y));

void win_get_origin     Params((window_desc win, int *x, int *y));
void bm_get_origin      Params((bitmap_desc win, int *x, int *y));

void win_restore_origin Params((window_desc win, int x, int y));
void bm_restore_origin  Params((bitmap_desc win, int x, int y));

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Drawing modes -- may be combined with bitwise logical operators: & | ^ ~ */

typedef enum DRAWING_MODE {
  MODE_SRC    = 0xc,
  MODE_DST    = 0xa,

  MODE_CLEAR  = (MODE_SRC&~MODE_SRC),
  MODE_SET    = (MODE_SRC|~MODE_SRC),
  MODE_XOR    = (MODE_SRC^MODE_DST)
} drawing_mode;

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/*  Drawing primitives */

void win_draw_line Params((
  window_desc win, 
  int x1, int y1, 
  int x2, int y2, 
  drawing_mode mode));

void bm_draw_line Params((
  bitmap_desc bm, 
  int x1, int y1, 
  int x2, int y2, 
  drawing_mode mode));

void      win_draw_string Params((
  window_desc win, 
  int         x, 
  int         y, 
  char       *str, 
  font_desc   fnt, 
  drawing_mode         mode));

void      bm_draw_string Params((
  window_desc win, 
  int         x, 
  int         y, 
  char       *str, 
  font_desc   fnt, 
  drawing_mode         mode));

void bm_raster_op Params((
  bitmap_desc dst, 
  int dx, 
  int dy, 
  int w, 
  int h, 
  bitmap_desc src, 
  int sx, 
  int sy,
  drawing_mode mode));

void win_raster_op Params((
  window_desc dst, 
  int dx, 
  int dy, 
  int w, 
  int h, 
  bitmap_desc src, 
  int sx, 
  int sy,
  drawing_mode mode));

void bm_pattern_raster_op Params((
  bitmap_desc dst,
  int dx, 
  int dy, 
  int w, 
  int h, 
  bitmap_desc pat, 
  int sx, 
  int sy,
  drawing_mode mode));

void win_pattern_raster_op Params((
  window_desc dst,
  int dx, 
  int dy, 
  int w, 
  int h, 
  bitmap_desc pat, 
  int sx, 
  int sy,
  drawing_mode mode));

void win_draw_bitmap Params((
  window_desc win, 
  int         x, 
  int         y, 
  bitmap_desc bm, 
  drawing_mode         mode));

void bm_draw_bitmap Params((
  bitmap_desc on_bm, 
  int         x, 
  int         y, 
  bitmap_desc from_bm, 
  drawing_mode         mode));

void bm_fill  Params((bitmap_desc bm, int x, int y, int val));
void win_fill Params((window_desc win, int x, int y, int val));

void win_clear_rect  Params((window_desc win, int x1,int y1, int x2,int y2));
void bm_clear_rect   Params((bitmap_desc win, int x1,int y1, int x2,int y2));

void win_invert_rect Params((window_desc win, int x1, int y1, int x2, int y2));
void bm_invert_rect  Params((bitmap_desc win, int x1, int y1, int x2, int y2));

void win_draw_pixel Params((window_desc win, int x, int y, int mode));
void bm_draw_pixel  Params((bitmap_desc bm,  int x, int y, int mode));

/*======================================================================*/
/* Inquiry and Manipulation Functions (common to windows and bitmaps) */
/*======================================================================*/

void win_get_size  Params((window_desc win, int *w_result, int *h_result));
void bm_get_size   Params((bitmap_desc bm, int *w_result, int *y_result));

int  win_width     Params((window_desc win));
int  bm_width      Params((bitmap_desc bm));

int  win_height    Params((window_desc win));
int  bm_height     Params((bitmap_desc bm));

int  bm_get_pixel  Params((bitmap_desc bm, int x, int y));
int  win_get_pixel Params((window_desc win, int x, int y));

void bm_map_pixels Params((bitmap_desc bm, long old, long new));

/*======================================================================*/
/* Window Routines */
/*======================================================================*/

/* Create and destroy windows */

window_desc win_new Params((
  int                 x,
  int                 y,
  int                 width,
  int                 height,
  event_handler_proc  event_proc, 
  opaque_pointer      data,
  char               *title));

void win_free Params((window_desc win));

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Query and manipulate window properties */

void win_set_pos        Params((window_desc win, int x, int y));
void win_get_pos        Params((window_desc win, int *x_result, int *y_result));

void win_set_size       Params((window_desc win, int w, int h));

void           win_set_vis     Params((window_desc win, int v));
int            win_get_vis     Params((window_desc win));

void           win_set_handler Params((window_desc win, event_handler_proc p));
event_handler_proc 
	       win_get_handler Params((window_desc win));

void           win_set_data    Params((window_desc win, opaque_pointer data));
opaque_pointer win_get_data    Params((window_desc win));

void           win_set_cursor  Params((window_desc w, cursor_desc c));

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Redraw batching control */

void win_start_batch    Params((window_desc win));
void win_end_batch      Params((window_desc win));

/*======================================================================*/
/* Bitmap Routines */
/*======================================================================*/

/* Create and destroy bitmaps */

bitmap_desc bm_new         Params((int w, int h));
bitmap_desc bm_read_file   Params((char *file_name));
bitmap_desc bm_build       Params((char *data, int w, int h));
void        bm_free        Params((bitmap_desc bm));

int         bm_write_file  Params((bitmap_desc bm, char *filename));

/*======================================================================*/
/* Cursor Routines   */
/*======================================================================*/

/* Create and destroy cursors */

cursor_desc cursor_open Params((char *filename, char *maskname));
void        cursor_free Params((cursor_desc c));

/*======================================================================*/
/* Font Routines */
/*======================================================================*/

/* Open and release */

font_desc font_open    Params((char *file_name));
font_desc font_default Params((void));
void      font_free    Params((font_desc fnt));

/* Deterimine size of string drawn in font */

void font_string_size Params((
  font_desc  in_font,
  char      *str,
  int       *w_result,
  int       *h_result,
  int       *baseline_result));


/*======================================================================*/
/* Backward Compatability */
/*======================================================================*/

/* 
 * Aliases for various names to maintain backward compatability 
 *   Define the symbol BACKWARD_COMPATIBLE to get access to old (pre v3.0) 
 *   routine names.  
 */

#ifdef BACKWARD_COMPATIBLE

        /*----------*/      /*----------*/
        /* OLD NAME */      /* NEW NAME */
        /*----------*/      /*----------*/
#define clear_rect          win_clear_rect
#define invert_rect         win_invert_rect
#define draw_line           win_draw_line
#define draw_string         win_draw_string
#define set_win_clip        win_set_clip
#define reset_win_clip      win_reset_clip
#define reduce_win_clip     win_reduce_clip
#define get_win_clip        win_get_clip
#define trivial_reject      win_trivial_reject  
#define draw_bitmap         win_draw_bitmap
#define new_window          win_new
#define free_window         win_free
#define set_win_vis         win_set_vis
#define get_win_vis         win_get_vis
#define set_win_data        win_set_data
#define get_win_data        win_get_data
#define set_win_pos         win_set_pos
#define get_win_pos         win_get_pos
#define start_draw          win_start_batch
#define end_draw            win_end_batch

#define put_bitmap_bit      bm_draw_pixel
#define get_bitmap_bit      bm_get_pixel
#define fill_bitmap         bm_fill

#define new_bitmap          bm_new
#define get_bitmap          bm_read_file
#define static_bitmap       bm_build
#define free_bitmap         bm_free
#define bitmap_w            bm_width
#define bitmap_h            bm_height

#define write_bitmap_file(f,b)   bm_write_file((b),(f))
#define raster_op(d,dx,dy,w,h,m,s,sx,sy) \
		 bm_raster_op((d),(dx),(dy),(w),(h),(s),(sx),(sy),(m))
#define pattern_raster_op(d,dx,dy,w,h,m,s,sx,sy) \
		 bm_pattern_raster_op((d),(dx),(dy),(w),(h),(s),(sx),(sy),(m))

#define start_windows       winpack_init
#define main_loop           winpack_main_loop
#define end_windows         winpack_shutdown

#define get_font            font_open
#define default_font        font_default
#define free_font           font_free
#define string_size         font_string_size

#define get_cursor          cursor_open
#define free_cursor         cursor_free
#define set_cursor          win_set_cursor

typedef input_mask          event_mask;  
#define in_evt_mask         event_in_mask

#define EMASK_loc_but_dn    loc_but_dn
#define EMASK_loc_but_up    loc_but_up
#define EMASK_but_dn        but_dn
#define EMASK_but_up        but_up
#define EMASK_loc_move      loc_move
#define EMASK_val_move      val_move
#define EMASK_key_dn        key_dn
#define EMASK_key_up        key_up
#define EMASK_win_enter     win_enter
#define EMASK_win_exit      win_exit
#define EMASK_special_event special_event
#define EMASK_window_resize win_resize
#define EMASK_ALL           EVENT_all
#define EMASK_NONE          EVENT_none

#endif /* BACKWARD_COMPATIBLE */

void wp_user_event(long,window_desc);
void mt_set_color(long);
long mt_get_color(int,int,int);
long mt_white(void);
long mt_black(void);
void wp_win_fill(window_desc,int,int,int,int,int);
void wp_win_saveblk(bitmap_desc,int,int,int,int,window_desc,int,int,int);
