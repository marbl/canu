
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
/* $Id: layout_uint.c,v 1.7 2007-02-08 02:04:56 brianwalenz Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>

#include "button.h"
#include "hash.h"
#include "agrep.h"

#include "AS_global.h"  //  F_S64 and int64, used to be cds.h

#undef  DEBUG_READ
#undef  DEBUG_SORT
#undef  DEBUG_ATTACH
#undef  DEBUG_ROWCON
#undef  DEBUG_PARSE
#undef  DEBUG_QUERY

#ifdef DEBUG_READ
#define DEBUG_LINKIDS
#endif
#ifdef DEBUG_SORT
#define DEBUG_LINKIDS
#endif

#define TRUNCATE_NEGATIVE_COORDINATES
#define COORDSCALE 10

#define MAX_LINE_LEN (128 * 4096)  /* Input line limit */
#define LIST_BUFFER  (16 * 4096)  /* Limit on list lengths */

#define PADDING      25    /* Min # of bps between pieces on same row */
#define RULER_HEIGHT 20    /* Pixel height of ruler at bottom of canvas */
#define PICK_THRESH   3    /* Max # of pixels from target for a pick */


/**** GLOBAL TYPES AND OBJECTS ****/

#define SOLID 0         /* Line style scalar values */
#define DASH  1

#define LINESTYLE 0x1   /* Bitvector set of object types */
#define MARKSTYLE 0x2
#define LINKSTYLE 0x4
#define OVLPSTYLE 0x8
#define ALLSTYLES 0xF

#define SHOWSTYLE 0x10
#define ANYSTYLES 0x1F

/* 44(60) #Style + [52(84),60(100)] #Lines + 8 #Points +
   40(56) #Links + 4(8) #Card_Links  + 20(36) #Overlaps + 1 #Comments
   4(8) #Clusters + 4(8) #Rows.
*/

typedef struct LSTag
  { struct LSTag *next;  /* Link for list of all such records */
    int    idnum;        /* Id number of style */
    char  *label;        /* Label of style */
    long   color;        /* Color */
    int    style;        /* In {SOLID,DASH} */
    int    thick;        /* Thickness (>= 1) */
    int    usage;        /* Type of objects in this style, <= ALLSTYLES */
    int    chord;        /* Current vis of each object type, <= ALLSTYLES */
  } Style;

typedef struct LNTag
  { struct LNTag *next;  /* Link for list of all such records */
    int    idnum;        /* Id number of line */
    int    locale;       /* Locale of first segment in line (for # ops) */
    char  *idnam;        /* Identifier name of line */
    char  *label;        /* Label of line */
    int    rowcon;       /* Row constraint (if >= 0), -1 if none, and -(r+3)
                            if also search selected (where r is row code). */
    int    cluster;      /* Cluster membership (-1 if none) */
    int64    *beg;          /* [*beg,*end) is an list of coords with */
    int64   *end;          /*   style refs inbetween.               */
    int64    attach;       /* Attachment point for links (-1 if none) */
#ifdef DEBUG_ATTACH
    int64  ottach;    /* Unadjusted value */
#endif
    int    row;          /* Current row in which line is placed */
    int64    extent;       /* Max_L ( end[-1] ) over all L s.t. *(L.beg) < *beg */
    struct LKTag *high;  /* Highest link attached to this line */
  } Line;

typedef struct LKTag
  { struct LKTag *next;  /* Link for list of all such records  */
#ifdef DEBUG_LINKIDS
    int    idnum;
#endif
    char  *label;        /* Label of link */
    Line **beg;          /* [*beg,*end) is a list of all lines */
    Line **end;          /*    in this link group.             */
    int    att;          /* Attribute for group.               */
    int64    min;          /* Min attachment coord  */
    int64    max;          /* Max attachment coord  */
    int64   emin;          /* Min fragment position  */
    int64   emax;          /* Max fragment position  */
    int    row;          /* Current row in which group bus bar is placed */ 
    int64    extent;       /* Max_G ( G->max ) over all G s.t. (G.min) < min */
  } Link;

typedef struct OVTag
  { struct OVTag *next;  /* Link for list of all such records */
    Line *id1, *id2;     /* Lines joined by overlap */
    char  *label;        /* Label of overlap */
    int att;             /* Attribute for overlap. */
  } Overlap;

static Style *StyleList;    /* "next"-linked list of all style records */
static Style *SlotsList;    /* "next"-linked list of all available styles */
static Line  *LineList;     /* "next"-linked list of all line records */

static int   NumStyles;     /* Number of styles */
static int   AllStyles;     /* Number of styles including query temps */
static int   NumLines;      /* Number of lines */
static int   NumEvents;     /* Total number of all line segment boundaries */

static Line  **LineIndex;   /* [0..NumLines-1] point to line records. Sorted so
                               that *LineIndex[i].beg <= *LineIndex[i+1].beg */
static Style **StyleIndex;  /* [0..NumStyles-1] point to style records */

static int   NumRowCons;    /* Constraints between rows [0,NumRowCons-1] */
static int   NumClusters;   /* Clusters #'d between [0,NumClusters-1] */

static Line ***RowConIndex;  /* Lines constrained to be in row i =         */
static Line  **RowCons;      /*   { RowCons[j] : j in [RCI[i],RCI[i+1]) }  */
static Line ***ClusterIndex; /* Lines in cluster i =                       */
static Line  **Clusters;     /*   { Clusters[j] : j in [CI[i],CI[i+1]) }   */

static int64   MinPnt, MaxPnt; /* Min and Max x-coords in layout */
static int64   SmallestLen;    /* Length of shortest piece in layout */
static int   Rows, Links;    /* Number of rows needed to layout all pieces */

static Link    *LinkList; /* "next"-linked list of all link records */
static Overlap *OvlpList; /* "next"-linked list of all overlap records */

static int      NumLinks;  /* Number of link sets */
static int      NumOvlps;  /* Number of overlaps */
static int      NumAttach; /* Number of link attachment points */

static Link   **LinkIndex; /* [0..NumLinks-1] point to link records.  Sorted so
                              that *LinkIndex[i].min <= *LinkIndex[i+1].min */
static Link  **ExtentIndex;/* [0..NumLinks-1] point to link records.  Sorted so
                              that *LinkIndex[i].emin <= *LinkIndex[i+1].emin */
static Line  **AttachIndex;/* [0..NumLines-1] point to line records.
                              Sorted so that *LineIndex[i].attach
                                                <= *LineIndex[i+1].attach */


/**** INPUT MODULE ****/
/*       Read in and build initial data structures for a layout */

static char IOBuffer[MAX_LINE_LEN];  /* Input line buffer */

static char *CurrentFile = NULL;

static char *SpaceAdvance(char *s)   /* Skip blanks */
{ while (isspace(*s))
    s += 1;
  return (s);
}

void SetCurrentFile(char *name)
{ CurrentFile = name; }

static void Error(char *msg, int data)  /* Report an error and quit */
{ if (CurrentFile != NULL)
    fprintf(stderr,"In file: %s\n",CurrentFile);
  fprintf(stderr,"    ");
  fprintf(stderr,msg,data);
  fprintf(stderr,"\n");
  if (*IOBuffer != '\0')
    fprintf(stderr,"\t%s",IOBuffer);
  exit (1);
}

/* Check, parse, and build data structures in response to spec file, "file" */

static char *NoInfoString = "<No Info>";

static char *CaptureString(char **scan, int attname)
{ static int   firstime = 1;
  static char *combuffer;
  static int   cptr;

  char *beg, *end, *rez;
  int   a, n;

  if (firstime)
    { firstime = 0;
      cptr = 0;
      combuffer = (char *) malloc(sizeof(char)*LIST_BUFFER);
    }

  beg = *scan + 1;
  while (isspace(*beg))
    beg += 1;
  end = beg;
  while (*end != '\0')
    end += 1;
  *scan = end;
  while (end > beg && isspace(end[-1]))
    end -= 1; 
  if (beg >= end)
    return (NoInfoString);
  else if (attname)
    { end = beg;
      while (!isspace(*end))
        end += 1;
    }
  a = *end;
  *end = '\0';
  n = (end-beg) + 1;
  if (cptr + n > LIST_BUFFER)
    { combuffer = (char *) malloc(sizeof(char)*LIST_BUFFER);
      cptr = 0;
      if (n > LIST_BUFFER)
        Error("Line has too long a comment (%d chars)",n-1);
    }
  rez = combuffer+cptr;
  strncpy(rez,beg,n);
  cptr += n;
  *end = a;
  return (rez);
}

static HashTable *LineTable;
static HashTable *StyleTable;
static HashTable *ClusterTable;
static HashTable *NameTable;

static int *NameToStyle;

static long GetIdentifier(HashTable *table, char *s, char **end, int isdef)
{ char *t;
  long   a, id;
  char  message[MAX_LINE_LEN];

  t = s;
  if (!isdigit(*t++))
    Error("Expecting object name",0);
  while (isalnum(*t))
    t++;
  a = *t;
  *t = '\0';
  id = HashLookup(s,table,&isdef);
  if (id < 0)
    id = HashAdd(s,table,isdef);
  else if (isdef < 0)
#ifdef VWR_TOLERATE_TWICE_DEFINED_NAMES
	{
		char str[1000];
		fprintf(stderr,"Warning: Name '%s' is defined more than once\n",s);
		sprintf(str,"x%s",s);
		id=HashLookup(str,table,&isdef);
		if(id<0)
			id=HashAdd(str,table,isdef);
	}
#else
    { sprintf(message,"Name '%s' is defined twice",s);
      *t = a;
      Error(message,0);
    }
#endif
  *t = a;
  *end = t;
  return (id);
}

static void ScanIdentifier(char *s, char **end)
{ char *t;

  t = s;
  if (!isdigit(*t++))
    Error("Expecting object designator",0);
  while (isalnum(*t))
    t++;
  *end = t;
}

void ReadAssembly(FILE *file)
{ static long   color;
  static int    thick, style;
  static int    hex[256];
  static int64   *int64buffer;
  static int64    xptr;
  static Line **linebuffer;
  static int    lptr;
  static int   *glosty, glocnt, glomax;

  static int FirstCall = 1;

  if (FirstCall)
    { FirstCall = 0;

      { int i;

        for (i = 0; i < 256; i++)    /* Hex digit to int table setup */
          hex[i] = -1;
        for (i = '0'; i <= '9'; i++)
          hex[i] = i-'0';
        for (i = 'A'; i <= 'F'; i++)
          hex[i] = 10 + (i-'A');
        for (i = 'a'; i <= 'f'; i++)
          hex[i] = 10 + (i-'f');
      }

      LineTable = NewHashTable(0x2000-1);
      StyleTable = NewHashTable(0x100-1);
      ClusterTable = NewHashTable(0x400-1);
    
      color = mt_white();         /* Setup default style attributes */
      thick = 1;
      style = SOLID;

      /* Initially build linked lists of each defined object and count their
         number.  Build indices to said objects later.                     */

      StyleList = NULL;
      NumStyles = 0;

      LineList  = NULL;
      NumLines  = 0;
      NumEvents = 0;

      LinkList  = NULL;
      NumLinks  = 0;

      OvlpList  = NULL;
      NumOvlps  = 0;

      xptr = 0;
      int64buffer = (int64 *) malloc(sizeof(int64)*LIST_BUFFER);

      lptr = 0;
      linebuffer = (Line **) malloc(sizeof(Line *)*LIST_BUFFER);

      /* Setup initial category application loop */

      glomax = 10;
      glosty = (int *) malloc(sizeof(int)*glomax);
      glocnt = 1;
      glosty[0] = -1;
    }

  /* For each line, parse, check, and build record for object */

  IOBuffer[MAX_LINE_LEN-2] = '\n';
  while (fgets(IOBuffer,MAX_LINE_LEN,file) != NULL)
    { char *s, *idnam;
      int type, idnum;

#define LINETYPE 0
#define ATTRTYPE 1
#define LINKTYPE 2
#define OVLPTYPE 3

      if (IOBuffer[MAX_LINE_LEN-2] != '\n')
        Error("Line too long",0);

      /* Process line prefix, determining line type and id # (if relevant) */

      s = SpaceAdvance(IOBuffer);
      if (strncmp(s,"ATT:",4) == 0)
        type = ATTRTYPE;
      else if (strncmp(s,"LNK:",4) == 0)
        type = LINKTYPE;
      else if (strncmp(s,"OVL:",4) == 0)
        type = OVLPTYPE;
      else
        { char *end, *t;

          ScanIdentifier(s,&end);
          if (*end != ':')
            Error(": should immediately follow name",0);

          t = SpaceAdvance(end+1);
          if (isdigit(*t)
#ifdef TRUNCATE_NEGATIVE_COORDINATES
              || (*t) == '-'
#endif
              )
            { type = LINETYPE;
              idnum = GetIdentifier(LineTable,s,&end,1);
              idnam = GetLastHash(LineTable);
              s = t;
              goto linelist;
	    }
          else
            { idnum = GetIdentifier(StyleTable,s,&end,1);
              s = t;
              goto attlist;
            }
        }
      s = SpaceAdvance(s+4);
      goto idlist;

    /* Process category line */

    attlist:
      { char *com;

        while (*s != '\0' && *s != '#')
          { if (*s == 'C')
              { int i, red, green, blue;
    
                for (i = 1; i <= 6; i++)
                  if (hex[(int) s[i]] < 0)
                    Error("Illegal color",0);
                red   = hex[(int) s[1]]*16 + hex[(int) s[2]];
                green = hex[(int) s[3]]*16 + hex[(int) s[4]];
                blue  = hex[(int) s[5]]*16 + hex[(int) s[6]];
                color = mt_get_color(red,green,blue);
                s = SpaceAdvance(s+7);
              }
            else if (*s == 'D')
              { style = DASH; s = SpaceAdvance(s+1); }
            else if (*s == 'S')
              { style = SOLID; s = SpaceAdvance(s+1); }
            else if (*s == 'T')
              { char *end;
                thick = strtol(s+1,&end,10);  
                if (end == s+1 || thick <= 0)
                  Error("Illegal thickness",0);
                s = SpaceAdvance(end);
              }
            else
              Error("Unrecognizable drawing attribute",0);
          }
        com = NoInfoString;
        if (*s == '#')
          com = CaptureString(&s,1);
        if (*s != '\0')
          Error("Invalid attribute spec.",0);

        { Style *newa;
          newa = (Style *) malloc(sizeof(Style));
          newa->idnum = idnum;
          newa->label = com;
          newa->color = color;
          newa->style = style;
          newa->thick = thick;
          newa->usage = 0;
          newa->chord = ALLSTYLES;
          newa->next  = StyleList;
          StyleList  = newa;
          NumStyles += 1;
        }
      }
      continue;

    /* Process link, overlap, and category loop lines (similar syntax) */

    idlist:
      { char *end, *t, *com;
        int   att, cnt;

        /* Syntax scan and count */
    
        t = s;
        for (cnt = 0; *t != '\0' && *t != 'A' && *t != '#'; cnt += 1)
          { ScanIdentifier(t,&end);
            if (end == t)
              Error("Expecting name in list",0);
            t = SpaceAdvance(end);
          }
        att = -1;
        if (*t == 'A')
          { att = GetIdentifier(StyleTable,t+1,&end,0);
            if (t+1 == end)
              Error("%c-qualifier must be a style name",*t);
            t = SpaceAdvance(end);
          }
        com = NoInfoString;
        if (*t == '#' && type != ATTRTYPE)
          com = CaptureString(&t,0);
        if (*t != '\0')
          if (type == ATTRTYPE)
            Error("Invalid category loop spec.",0);
          else if (type == LINKTYPE)
            Error("Invalid link set spec.",0);
          else
            Error("Invalid overlap spec.",0);

        /* Build and record for each line type */

        switch (type)
        { case ATTRTYPE:
            if (att >= 0)
              Error("Category loop does not take a style qualifier",0);
            if (glomax < cnt)
              { glomax = cnt + 10;
                glosty = (int *) realloc(glosty,sizeof(int)*glomax);
              }
            glocnt = cnt;
            for (cnt = 0; cnt < glocnt; cnt++)
              { glosty[cnt] = GetIdentifier(StyleTable,s,&end,0);  
                s = SpaceAdvance(end);
              }
            break;

          case LINKTYPE:
            if (lptr + cnt > LIST_BUFFER)
              { linebuffer = (Line **) malloc(sizeof(Line *)*LIST_BUFFER);
                lptr = 0;
                if (cnt > LIST_BUFFER)
                  Error("Link set has too many lines (has %d)",cnt);
              }

            { Link *newa; 

              newa = (Link *) malloc(sizeof(Link));
              newa->beg = linebuffer + lptr;
              while (cnt-- > 0)
                { linebuffer[lptr++] =
                      (Line *) GetIdentifier(LineTable,s,&end,0);
                  s = SpaceAdvance(end);
                }
              newa->end  = linebuffer + lptr;
              if (att < 0)
                newa->att = glosty[0];
              else
                newa->att = att;
              newa->label = com;
#ifdef DEBUG_LINKIDS
              newa->idnum = NumLinks;
#endif
              newa->next = LinkList;
              LinkList  = newa;
              NumLinks += 1;
            }
            break;

          case OVLPTYPE:
            if (cnt != 2)
              Error("Overlap between more than two pieces?",0);

            { Overlap *newa; 
              newa = (Overlap *) malloc(sizeof(Overlap));
              newa->id1 = (Line *) GetIdentifier(LineTable,s,&end,0);
              s = SpaceAdvance(end);
              newa->id2 = (Line *) GetIdentifier(LineTable,s,&end,0);
              if (att < 0)
                newa->att = glosty[0];
              else
                newa->att = att;
              newa->label = com;
              newa->next = OvlpList;
              OvlpList  = newa;
              NumOvlps += 1;
            }
            break;
        }
      }
      continue;

    /* Process piece definitions */

    linelist:
      { char *end, *t, *com;
        int   row, chk, cnt; 
        int64 num;
        int   wasA, wasM;

        /* Syntax scan and count */
    
        t = s;
        cnt = wasA = wasM = 0;
        while (*t != '\0' && *t != 'R' && *t != 'C' && *t != '#')
          { if (cnt == 0 && (*t == 'A' || *t == 'M'))
              Error("First coord of line must be a non-mark coord",0);
            if (*t == 'A')
              { if (wasA)
                  Error("Two consecutive category refs",0);
                wasA = 1;
                wasM = 0;
                t   += 1;
              }
            else if (*t == 'M')
              { wasM = 1;
                wasA = 0;
                t   += 1;
              }
            else
              wasM = wasA = 0;
            if (wasA)
              num = GetIdentifier(StyleTable,t,&end,0);  
            else
              { num = STR_TO_INT64(t,&end,10);  
              if (num < 0) {
#ifndef TRUNCATE_NEGATIVE_COORDINATES
                Error("Coordinates must be non-negative ints",0);
#else
                fprintf(stderr,
                        "Warning: Coordinates must be non-negative ints\n"
                        "Truncating " F_S64 " to zero\n", num);
                num = 0;
#endif                
                }
                cnt += 1;
              }
            if (end == t)
              Error("Invalid line spec.",0);
            t = SpaceAdvance(end);
          }
        if (wasM || wasA)
          Error("Last coord of line must be a non-mark coord",0);

        chk = row = -1;
        while (*t != '\0' && *t != '#')
          { switch (*t)
            { case 'C':
                chk = GetIdentifier(ClusterTable,t+1,&end,0);
                break;
              case 'R':
                row = strtol(t+1,&end,10);
                if (row < 0)
                  Error("Row constraint < 0",0);
                break;
              default:
                Error("Invalid line spec.",0);
            }
            if (t+1 == end)
              Error("%c-qualifier must be an int >= 0",*t);
            t = SpaceAdvance(end);
          }

        com = NoInfoString;
        if (*t == '#')
          com = CaptureString(&t,0);
        if (*t != '\0')
          Error("Invalid line spec.",0);

        /* Enough space for coord list? */

        if (xptr + (2*cnt-1) > LIST_BUFFER)
          { int64buffer = (int64 *) malloc(sizeof(int64)*LIST_BUFFER);
            xptr = 0;
            if (2*cnt-1 > LIST_BUFFER)
              Error("Line has too many points (has %d)",cnt);
          }

        /* Build and record */
        
        { Line *newa;
          int64   ap, num, lst;
      
          newa = (Line *) malloc(sizeof(Line));
          newa->idnum   = idnum;
          newa->idnam   = idnam;
          newa->rowcon  = row;
          newa->attach  = 0;
          newa->cluster = chk;
          newa->beg     = int64buffer + xptr;
          newa->locale  = NumEvents - NumLines;
          ap = num = 0;
          while (cnt-- > 0)
            { lst = num;
              if (*s == 'M')
                num = STR_TO_INT64(s+1,&end,10);
              else
                num = STR_TO_INT64(s,&end,10);  
              if (num < 0 ) {
#ifndef TRUNCATE_NEGATIVE_COORDINATES
                Error("Coordinates must be non-negative",0);
#else                
                fprintf(stderr,
                        "Warning: Coordinates must be non-negative ints\n"
                        "Truncating " F_S64 " to zero\n", num);
                num = 0;
#endif
              }
              if (lst > num)
                Error("Division pts out of order",0);
              if (*s == 'M')
                int64buffer[xptr++] = -(num+1)/COORDSCALE;  // divide by COORDSCALE here
              else
                int64buffer[xptr++] = num/COORDSCALE;
              NumEvents += 1;
              s = SpaceAdvance(end);
              if (cnt > 0)
                if (*s == 'A')
                  { int64buffer[xptr++] = GetIdentifier(StyleTable,s+1,&end,0);
                    s = SpaceAdvance(end);
                  }
                else
                  int64buffer[xptr++] = glosty[(ap++) % glocnt];
            }
          newa->end   = int64buffer + xptr;
          newa->label = com;
          newa->next  = LineList;
          LineList   = newa;
          NumLines  += 1;
        }
      }
    }
	if(NumLines<=0)
	{
		fprintf(stderr,"Nothing to draw? bye\n");
		exit(1);
	}
}

void BuildAssembly(void)
{

  /* Debug output of syntax and build scan */

#ifdef DEBUG_READ
  { Line    *ln, **lp;
    Link    *lk;
    Style   *cs;
    Overlap *ov;
    int64     *p;


    printf("\nFIRST SCAN CONSTRUCTS:\n");

    PrintHashTable(LineTable);
    PrintHashTable(StyleTable);
    PrintHashTable(ClusterTable);
    printf("\n");

    for (ln = LineList; ln != NULL; ln = ln->next)
      { printf("  Line %d: ",ln->idnum);
        for (p = ln->beg; p < ln->end; p++)
          if ((p-ln->beg)%2 == 0)
            if (*p < 0)
              printf("M%u ",-(*p+1));
            else
              printf("%u ",*p);
          else
            printf("(%u) ",*p);
        printf(" (R=%d,C=%d) '%s' '%s' (loc=%d)\n",
               ln->rowcon,ln->cluster,ln->idnam,ln->label,ln->locale);
      }
    printf("\n");
    for (cs = StyleList; cs != NULL; cs = cs->next)
      { printf("  Style %d:",cs->idnum);
        printf(" (C=%lx S=%d T=%d) '%s'\n",
               cs->color,cs->style,cs->thick,cs->label);
      }
    printf("\n");
    for (lk = LinkList; lk != NULL; lk = lk->next)
      { printf("  %d:",lk->idnum);
        for (lp = lk->beg; lp < lk->end; lp++)
          printf(" %d",(int) (*lp));
        printf(" A%d '%s'\n",lk->att,lk->label);
      }
    printf("\n");
    for (ov = OvlpList; ov != NULL; ov = ov->next)
      printf("  Ovlp: %d %d A%d '%s'\n",
             (int) (ov->id1),(int) (ov->id2),ov->att,ov->label);
  }
#endif

  *IOBuffer = '\0';  /* Don't echo line in subsequent error messages */

  { char *name;
  
    name = HashRefButNoDef(StyleTable);
    if (name != NULL)
      { sprintf(IOBuffer+1,"Style '%s' referenced but not defined",name);
        Error(IOBuffer+1,0);
      }
    name = HashRefButNoDef(LineTable);
    if (name != NULL)
      { sprintf(IOBuffer+1,"Line '%s' referenced but not defined",name);
        Error(IOBuffer+1,0);
      }
  }

  /* Fill category index and check that id #s are densely allocated */

  { Style *cs; int i;

    AllStyles  = NumStyles + 3;
    StyleIndex = (Style **) malloc(sizeof(Style *)*(AllStyles + 1));
    StyleIndex += 1;

    for (cs = StyleList; cs != NULL; cs = cs->next)
      StyleIndex[cs->idnum] = cs;
    StyleIndex[-1] = cs = (Style *) malloc(sizeof(Style));
    cs->next  = StyleList;
    StyleList = cs;
    cs->idnum = -1;
    cs->label = "Default";
    cs->color = mt_white();
    cs->style = SOLID;
    cs->thick = 1;
    cs->usage = 0;
    cs->chord = ALLSTYLES;

    cs = (Style *) malloc(sizeof(Style)*(AllStyles - NumStyles));
    for (i = NumStyles; i < AllStyles; i += 1)
      { StyleIndex[i] = cs + (i-NumStyles);
        StyleIndex[i]->idnum = i;
        StyleIndex[i]->usage = 0;
      }
    for (i = NumStyles; i < AllStyles-1; i += 1)
      StyleIndex[i]->next = StyleIndex[i+1];
    StyleIndex[i]->next = NULL;
    SlotsList = StyleIndex[NumStyles];
  }

  { int i, id, isdef;

    NameTable = NewHashTable(0x100-1);

    NameToStyle = (int *) malloc(sizeof(int)*(NumStyles+1));
    isdef = 0;
    for (i = -1; i < NumStyles; i++)
      { id = HashLookup(StyleIndex[i]->label,NameTable,&isdef);
        if (StyleIndex[i]->label == NoInfoString)
          { fprintf(stderr,"Warning: Style %d does not have a name.\n",i);
            fprintf(stderr,"         Queries will not work.\n");
          }
        else if (id >= 0)
          { fprintf(stderr,"Warning: Style name '%s' is used more than once.\n",
                           StyleIndex[i]->label);
            fprintf(stderr,"         Queries will not work.\n");
          }
        else
          id = HashAdd(StyleIndex[i]->label,NameTable,isdef);
        NameToStyle[id] = i;
      }
  }

  /* Fill line index and check that id #s are densely allocated.
     Also note category usages (in relevant style records) and count
     constraint rows and utilized clusters in prepartion of building
     indices for said objects.                                       */

  { int  maxrow, maxchk;
    int  numrow, numchk;

    maxrow = maxchk = -1;
    numrow = numchk = 0;

    LineIndex = (Line **) malloc(sizeof(Line *)*NumLines);

    { Line *ln;
      int64  *p;

      for (ln = LineList; ln != NULL; ln = ln->next)
        { LineIndex[ln->idnum] = ln;
          if (ln->cluster >= 0)
            { if (maxchk < ln->cluster) maxchk = ln->cluster;
              numchk += 1;
            }
          if (ln->rowcon >= 0)
            { if (maxrow < ln->rowcon) maxrow = ln->rowcon;
              numrow += 1;
            }
          for (p = ln->beg+1; p < ln->end; p += 2)
            { if (p[1] >= 0)
                StyleIndex[*p]->usage |= LINESTYLE;
              else
                StyleIndex[*p]->usage |= MARKSTYLE;
            }
        }
    }

    /* Allocate and build row constraint and cluster indices.  For each
       we want to now the set of lines in each row and/or cluster.   */

    NumRowCons   = maxrow+1;
    NumClusters  = maxchk+1;
    RowConIndex  = (Line ***) malloc(sizeof(Line **)*(maxrow+2));
    ClusterIndex = (Line ***) malloc(sizeof(Line **)*(maxchk+2));
    RowCons =  (Line **) malloc(sizeof(Line *)*numrow);
    Clusters = (Line **) malloc(sizeof(Line *)*numchk);

    { int i;
      for (i = 0; i <= maxrow; i++)
        RowConIndex[i] = RowCons;
      for (i = 0; i <= maxchk; i++)
        ClusterIndex[i] = Clusters;
    }

    { Line *ln;

      for (ln = LineList; ln != NULL; ln = ln->next)
        { if (ln->cluster >= 0)
            ClusterIndex[ln->cluster] += 1;
          if (ln->rowcon >= 0)
            RowConIndex[ln->rowcon] += 1;
        }
    }

    { int i;

      for (i = 1; i <= maxrow; i++)
        RowConIndex[i] += (RowConIndex[i-1] - RowCons);
      RowConIndex[maxrow+1] = RowCons + numrow;
      for (i = 1; i <= maxchk; i++)
        ClusterIndex[i] += (ClusterIndex[i-1] - Clusters);
      ClusterIndex[maxchk+1] = Clusters + numchk;
    }

    { Line *ln;
      Line **c;

      for (ln = LineList; ln != NULL; ln = ln->next)
        { if (ln->cluster >= 0)
            { c = ClusterIndex[ln->cluster]-1; 
              *c = ln;
              ClusterIndex[ln->cluster] = c;
            }
          if (ln->rowcon >= 0)
            { c = RowConIndex[ln->rowcon]-1; 
              *c = ln;
              RowConIndex[ln->rowcon] = c;
            }
        }
    }

    /* Build index of link constraints.  Convert int line references
       to pointer references and note style usages.  Build reverse index
       of all links attaching to a given line.                          */

    LinkIndex = (Link **) malloc(sizeof(Overlap *)*NumLinks);

    { int i, j;
      int64 x1, x2;
      Link *lk;
      Line *ln, **lp;

      lk = LinkList;
      for (i = 0; i < NumLinks; i++)
        { LinkIndex[i] = lk;
          for (lp = lk->beg; lp < lk->end; lp++)
            { *lp = LineIndex[(long) (*lp)];
              (*lp)->attach += 1;
            }
          StyleIndex[lk->att]->usage |= LINKSTYLE;
          lk = lk->next;
        }

      NumAttach = 0;
      for (i = 0; i < NumLines; i++)
        { ln = LineIndex[i];
          if (ln->attach > 0)
            NumAttach += 1;
        }

      AttachIndex = (Line **) malloc(sizeof(Line *)*NumAttach);

      j = 0;
      for (i = 0; i < NumLines; i++)
        if ((ln = LineIndex[i])->attach > 0)
          { AttachIndex[j++] = ln;
            x1 = ln->beg[0];
            x2 = ln->end[-1];
	    /* the following monstrosity is needed [sic] 
	       in place of the simpler (x1+x2)/2
	       so that we do not have overflow problems on
	       the sum -ALH */
            ln->attach = x1/2+x2/2+(x1%2+x2%2)/2;
          }
        else
          ln->attach = -1;
    }

    /* Also convert int line references to pointer references and
       note style usages for all overlap constraints.              */

    { Overlap *ov;

      for (ov = OvlpList; ov != NULL; ov = ov->next)
        { ov->id1 = LineIndex[(long) (ov->id1)]; 
          ov->id2 = LineIndex[(long) (ov->id2)]; 
          StyleIndex[ov->att]->usage |= OVLPSTYLE;
        }
    }
  }

  /* Second pass debug output showing indices and quantities inferred
     on this pass.                                                      */

#ifdef DEBUG_READ
  { int      i, j; 
    int64 *p;
    Line    *ln, **lp;
    Link    *lk, **kp;
    Style   *cs;
    Overlap *ov;

    printf("\nSECOND SCAN CONSTRUCTS:\n");
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        printf("  Line %d: ",ln->idnum);
        for (p = ln->beg; p < ln->end; p++)
          if ((p-ln->beg)%2)
            printf("%u ",*p);
          else
            printf("(%u) ",*p);
        printf(" (R=%d,C=%d)",ln->rowcon,ln->cluster);
        printf(" attach=%u\n",ln->attach);
      }
    printf("\n");
    for (i = -1; i < NumStyles; i++)
      { char *name;
        int   id;

        cs = StyleIndex[i];
        printf("  Style %d:",cs->idnum);
        printf(" (C=%lx S=%d T=%d) Usage = %x\n",
               cs->color,cs->style,cs->thick,cs->usage);
        name = cs->label;
        id   = 0;
        id   = HashLookup(name,NameTable,&id);
        printf("\t\t%s -> %d -> %d\n",name,id,NameToStyle[id]);
      }
    printf("\n");
    for (i = 0; i < NumRowCons; i++)
      if (RowConIndex[i+1] > RowConIndex[i])
        { printf("  RowSet %d:\n",i);
          for (lp = RowConIndex[i]; lp < RowConIndex[i+1]; lp++)
            printf("\t%d\n",(*lp)->idnum);
        }
    printf("\n");
    for (i = 0; i < NumClusters; i++)
      if (ClusterIndex[i+1] > ClusterIndex[i])
        { printf("  Cluster %d:\n\t",i);
          for (lp = ClusterIndex[i]; lp < ClusterIndex[i+1]; lp++)
            printf("%d ",(*lp)->idnum);
          printf("\n");
        }
    printf("\n");
    for (i = 0; i < NumLinks; i++)
      { lk = LinkIndex[i];
        printf("  Link %d:",i);
        for (lp = lk->beg; lp < lk->end; lp++)
          printf(" %d",(*lp)->idnum);
        printf(" A%d\n",lk->att);
      }
    printf("\n");
    for (ov = OvlpList; ov != NULL; ov = ov->next)
      printf("  Ovlp: %d %d A%d\n",ov->id1->idnum,ov->id2->idnum,ov->att);
  }
#endif
}


/**** MENU COMMUNICATION ROUTINES ****/
/*       Routines for communicating menu info to user interface and for
         receiving visibility changes from said.                        */

static int ListWalker;
static int WhichList;
static int TypeToMask[] = { LINESTYLE, MARKSTYLE, LINKSTYLE, OVLPSTYLE };

int StartStyleList(int type)
{ int i, count;

  WhichList = TypeToMask[type];
  count = 0;
  for (i = -1; i < NumStyles; i++)
    { if (StyleIndex[i]->usage & WhichList)
        count += 1; 
    }
  ListWalker = -1;
  return (count);
}

int NextUsedStyle(void)
{ while (ListWalker < NumStyles)
    { if (StyleIndex[ListWalker]->usage & WhichList)
        break;
      ListWalker++;
    }
  if (ListWalker >= NumStyles)
    return (ListWalker);
  return (ListWalker++);
}

int NumberOfStyles(void)
{ return (NumStyles); }

void GetStyle(int idnum, long *color, int *thick, int *dash, char **name)
{ *color = StyleIndex[idnum]->color;
  *thick = StyleIndex[idnum]->thick;
  *dash  = StyleIndex[idnum]->style;
  *name  = StyleIndex[idnum]->label;
}

void SetStyleVis(int idnum, int type, int vis)
{ Style *sty, *aty;
  int i, old;

  sty = StyleIndex[idnum];
  old = sty->chord;
  if (vis)
    sty->chord |= TypeToMask[type];
  else
    sty->chord &= ~TypeToMask[type];
  if (old != sty->chord)
    for (i = NumStyles; i < AllStyles; i++)
      { aty = StyleIndex[i];
        if ((aty->usage & SHOWSTYLE) != 0)
          if (aty->next == sty)
            aty->chord = sty->chord;
      }
}


/**** ROW ALLOCATION OF ALL OBJECTS ****/
/*       Routine establishes the drawing locations of all objects.   */

static int ATTACH_COMPARE(const void *l, const void *r)
{ int64 x, y;
  x = (*((Line **) l))->attach;
  y = (*((Line **) r))->attach;
  if (x > y ) return 1;
  if (x == y ) return 0;
  return -1;
}

static int64 AttachMin, AttachSep; 

static void MoveLeft(int i, int64 del)
{ int64   play, space;
  Line *ln;

  ln   = AttachIndex[i];
  play = ln->attach - (ln->beg[0] + AttachMin);
  if (play < del)
    del = play;
  if (i == 0)
    space = AttachSep;
  else
    space = (ln->attach - del) - AttachIndex[i-1]->attach;
  if (space < AttachSep)
    { MoveLeft(i-1,AttachSep-space);
      del = (ln->attach - AttachIndex[i-1]->attach) - AttachSep;
      if (del < 0) del = 0;
    }
  ln->attach -= del;
}

static void MoveRight(int i, int64 del)
{ int64   play, space;
  Line *ln;

  ln   = AttachIndex[i];
  play = (ln->end[-1] - AttachMin) - ln->attach;
  if (play < del)
    del = play;
  if (i+1 == NumAttach)
    space = AttachSep;
  else
    space = AttachIndex[i+1]->attach - (ln->attach + del);
  if (space < AttachSep)
    { MoveRight(i+1,AttachSep-space);
      del = (AttachIndex[i+1]->attach - ln->attach) - AttachSep;
      if (del < 0) del = 0;
    }
  ln->attach += del;
}

static void AttachmentSetup(void)
{ AttachSep = SmallestLen/10;
  AttachMin = SmallestLen/4;

  qsort(AttachIndex,NumAttach,sizeof(Line *),ATTACH_COMPARE);

#ifdef DEBUG_ATTACH
  { int i;
    for (i = 0; i < NumAttach; i++)
      AttachIndex[i]->ottach = AttachIndex[i]->attach;
  }
#endif

  { int i; 
    int64 del;

    for (i = 1; i < NumAttach; i++)
      { del = AttachIndex[i]->attach - AttachIndex[i-1]->attach;
        if (del < AttachSep)
          MoveLeft(i-1,AttachSep-del); 
        del = AttachIndex[i]->attach - AttachIndex[i-1]->attach;
        if (del < AttachSep)
          MoveRight(i,AttachSep-del);
      }
  }

#ifdef DEBUG_ATTACH
  { int   i; 
    int64 xl, xh;
    Line *ln;

    printf("\nATTACHMENT SORT:  min = %d  sep = %d\n",AttachMin,AttachSep);
    for (i = 0; i < NumAttach; i++)
      { ln = AttachIndex[i];
        printf("  %3d: %6d %6d",ln->idnum,ln->ottach,ln->attach);
        xl = ln->beg[0] + AttachMin;
        xh = ln->end[-1] - AttachMin;
        printf(" [%d,%d]",xl,xh);
        if (i > 0 && ln->attach - AttachIndex[i-1]->attach < AttachSep)
          printf(" **Violated**");
        if (ln->attach < xl || ln ->attach > xh) printf(" **Bug**");
        printf("\n");
      }
  }
#endif
}

static int LINE_COMPARE(const void *l, const void *r)
{ int64 x, y;
  x = *((*((Line **) l))->beg);
  y = *((*((Line **) r))->beg);
  if (x > y ) return 1;
  if ( x == y ) return 0;
  return -1;
}

static int LINK_COMPARE(const void *l, const void *r)
{ int64 x, y;
  int rc=-1;
  x = (*((Link **) l))->min;
  y = (*((Link **) r))->min;
  if (x > y ) rc=1;
  if (x == y ) rc=0;
  return rc;
}

static int EXTENT_COMPARE(const void *l, const void *r)
{ int64 x, y;
  x = (*((Link **) l))->emin;
  y = (*((Link **) r))->emin;
  if (x > y ) return 1;
  if (x == y ) return 0;
  return -1;
}

static int ROWCON_COMPARE(const void *l, const void *r)
{ Line *ln1, *ln2;
  ln1 = *((Line **) l);
  ln2 = *((Line **) r);
  if (ln1->rowcon != ln2->rowcon)
    return (ln1->rowcon - ln2->rowcon);
  if ( *(ln1->beg) > *(ln2->beg) ) return 1;
  if ( *(ln1->beg) == *(ln2->beg) ) return 0;
  return -1;
}

void LayoutAssembly(void)
{ Rows = 0;
  Links = 0;
  MinPnt = MaxPnt = 0;

  if (NumLines <= 0) return;

  { int   i; 
    int64 bg, ed;
    Line *ln;

    ln = LineIndex[0];
    MinPnt = *(ln->beg);
    MaxPnt = ln->end[-1];
    SmallestLen = MaxPnt - MinPnt;
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        if (ln->rowcon < 0)
          { if (ln->attach < 0)
              ln->rowcon = NumRowCons+1;
            else
              ln->rowcon = NumRowCons;
          }
        bg = *(ln->beg);
        ed = ln->end[-1];
        if (bg < MinPnt)
          MinPnt = bg;
        if (ed > MaxPnt)
          MaxPnt = ed;
        if (ed-bg < SmallestLen)
          SmallestLen = ed-bg;
      }
    SmallestLen += 1;
  }

  if (NumLinks > 0)
    { int i; 
      int64 x;
      Link  *lk;
      Line **lp;

      AttachmentSetup();

      for (i = 0; i < NumLinks; i++)
        { lk = LinkIndex[i];
          lp = lk->beg;
          lk->min = lk->max = (*lp)->attach;
          for (lp += 1; lp < lk->end; lp++)
            { x = (*lp)->attach;
              if (x < lk->min)
                lk->min = x;
              if (x > lk->max)
                lk->max = x;
            }
        }
    }

#ifdef DEBUG_SORT
      { Link *lk;
        int i;

        printf("\nLINK SORT:\n");
        for (i = 0; i < NumLinks; i++)
          { lk = LinkIndex[i];
            printf("   %d: %6d\n",lk->idnum,lk->min);
          }
      }
#endif
  if (NumLinks > 0)
    { qsort(LinkIndex,NumLinks,sizeof(Link *),LINK_COMPARE);

#ifdef DEBUG_SORT
      { Link *lk;
        int i;

        printf("\nLINK SORT:\n");
        for (i = 0; i < NumLinks; i++)
          { lk = LinkIndex[i];
            printf("   %d: %6d\n",lk->idnum,lk->min);
          }
      }
#endif
    }

  ExtentIndex = (Link **) malloc(sizeof(Overlap *)*NumLinks);
  if (NumLinks > 0)
    { int i; 
      int64 x;
      Link  *lk;
      Line **lp;

      for (i = 0; i < NumLinks; i++)
        { ExtentIndex[i] = lk = LinkIndex[i];
          lp = lk->beg;
          lk->emin = (*lp)->beg[0];
          lk->emax = (*lp)->end[-1];
          for (lp += 1; lp < lk->end; lp++)
            { x = (*lp)->beg[0];
              if (x < lk->emin)
                lk->emin = x;
              x = (*lp)->end[-1];
              if (x > lk->emax)
                lk->emax = x;
            }
        }

      qsort(ExtentIndex,NumLinks,sizeof(Link *),EXTENT_COMPARE);

#ifdef DEBUG_SORT
      { Link *lk;
        int i;

        printf("\nEXTENT SORT:\n");
        for (i = 0; i < NumLinks; i++)
          { lk = ExtentIndex[i];
            printf("   %d: %6d %6d\n",lk->idnum,lk->emin,lk->emax);
          }
      }
#endif
    }

  qsort(LineIndex,NumLines,sizeof(Line *),ROWCON_COMPARE);

#ifdef DEBUG_SORT
  { Line *ln;
    int i;

    printf("\nROWCON SORT:\n");
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        printf("   %3d: %3d %6d\n",ln->idnum,ln->rowcon,*(ln->beg));
      }
  }
#endif

  { Line **base, **finger;
    Line *last, *ln, *fn;
    Line  fakeln;
    int64   fakend; 
    int rowtop;
    int i, crow, mrow, try;

    fakeln.end = (&fakend) + 1; 
    fakeln.beg =  &fakend;
    fakend = MaxPnt + 100;

    if (NumRowCons < 100)
      rowtop = 100;
    else
      rowtop = NumRowCons;
    base   = (Line **) malloc(sizeof(Line *)*rowtop);
    finger = (Line **) malloc(sizeof(Line *)*rowtop);

    if (NumRowCons <= 0)
      Rows = 0;
    else
      Rows = NumRowCons-1;
    for (i = 0; i <= Rows; i++)
      base[i] = &fakeln;

    crow = -1;
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        if (ln->rowcon != crow)
          { crow = ln->rowcon;
            if (crow < NumRowCons)
              mrow = crow;
            else
              mrow = 0;
            finger[mrow] = base[mrow];

#ifdef DEBUG_ROWCON
            printf("\nNEXT POP %d,%d:\n",crow,Rows);
            for (j = 0; j <= Rows; j++)
              { printf("Row %d:",j);
                for (fn = base[j]; fn != &fakeln; fn = fn->next)
                  { printf(" %d[%d,%d]",fn->idnum,fn->beg[0],fn->end[-1]);
                    if (fn->end[-1] > 10000) break;
                  }
                printf("\n");
              }
#endif
          }

        if (crow >= NumRowCons)
          try = 0;
        else
          try = crow;
        while (1)
          { last = NULL;
            fn = finger[try];
            while (fn->end[-1]+PADDING < ln->beg[0])
              { last = finger[try];
                fn = fn->next;
              }
            if (ln->end[-1]+PADDING < fn->beg[0])
              { ln->row = try;
                if (last == NULL)
                  base[try] = ln;
                else
                  last->next = ln;
                ln->next = fn;
                finger[try] = ln;
                break;
              }
            else
              { finger[try] = fn;
                try += 1;
                if (try > Rows)
                  { if (try >= rowtop)
                     { rowtop = 1.2*rowtop + 50;
                       base   = (Line **) realloc(base,sizeof(Line *)*rowtop);
                       finger = (Line **) realloc(finger,sizeof(Line *)*rowtop);
                     }
                    Rows = try;
                    base[Rows] = &fakeln;
                  }
                if (try > mrow)
                  { mrow = try;
                    finger[mrow] = base[mrow];
                  }
              }
          }
      }

#ifdef DEBUG_ROWCON
    printf("\nFINAL POP **,%d:\n",Rows);
    for (j = 0; j <= Rows; j++)
      { printf("Row %d:",j);
        for (fn = base[j]; fn != &fakeln; fn = fn->next)
          { printf(" %d[%d,%d]",fn->idnum,fn->beg[0],fn->end[-1]);
            if (fn->end[-1] > 10000) break;
          }
        printf("\n");
      }
#endif

    free(base);
    free(finger);
    Rows += 1;
  }

  qsort(LineIndex,NumLines,sizeof(Line *),LINE_COMPARE);

#ifdef DEBUG_SORT
  { Line *ln;
    int i;

    printf("\nLINE SORT:\n");
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        printf("   %3d: %6d\n",ln->idnum,*(ln->beg));
      }
  }
#endif

  { int i; 
    int64 max;

    max      = LineIndex[0]->end[-1];
    LineIndex[0]->extent = max;
    LineIndex[0]->high   = NULL;
    for (i = 1; i < NumLines; i++)
      { if (LineIndex[i]->end[-1] > max)
          max = LineIndex[i]->end[-1];
        LineIndex[i]->extent = max;
        LineIndex[i]->high   = NULL;
      }
  }

  if (NumLinks > 0)
    { int i, j, last, linkrows, *rowstuff; 
      int64 bp;
      int64 max;

      last     = 0;
      linkrows = 0;
      max      = LinkIndex[0]->max;
      LinkIndex[0]->extent = max;
      for (i = 1; i < NumLinks; i++)
        { while (LinkIndex[last]->max+PADDING < LinkIndex[i]->min)
            last += 1;
          if (i-last > linkrows)
            linkrows = i-last;
          if (LinkIndex[i]->max > max)
            max = LinkIndex[i]->max;
          LinkIndex[i]->extent = max;
        }
      linkrows += 1;

      rowstuff  = (int *) malloc(sizeof(int)*linkrows);

      for (i = 0; i < linkrows; i++)
        rowstuff[i] = MinPnt-1;
      for (i = 0; i < NumLinks; i++)
        { bp = LinkIndex[i]->min;
          for (j = 0; j < linkrows; j++)
            if (rowstuff[j] < bp)
              { rowstuff[j] = LinkIndex[i]->max + PADDING;
                LinkIndex[i]->row = -(j+1);
                break;
              }
          if (j >= linkrows)
            { fprintf(stderr,"Fatal: overflowed row stack (bug)\n");
              exit (1);
            }
          else if (j >= Links)
            Links = j+1;
        }

      free(rowstuff);
    }

  { int i;
    Link *lk;
    Line **kp;

    for (i = 0; i < NumLinks; i++)
      { lk = LinkIndex[i];
        for (kp = lk->beg; kp < lk->end; kp++)
          if ((*kp)->high == NULL)
            (*kp)->high = lk;
          else if (lk->row < (*kp)->high->row)
            (*kp)->high = lk;
      }
  }

#ifdef DEBUG_SORT
  { Line *ln;
    Link *lk;
    int i;

    printf("\nLINE ALLOCATION:\n");
    for (i = 0; i < NumLines; i++)
      { ln = LineIndex[i];
        printf("  %d: [%d,%d] row = %d",ln->idnum,
               *(ln->beg),ln->end[-1],ln->row);
        if (ln->high != NULL) 
          printf(" high = %d(%d)",ln->high->idnum,ln->high->row);
        printf("\n");
      }

    if (NumLinks > 0)
      { printf("\nLINK ALLOCATION:\n");
        for (i = 0; i < NumLinks; i++)
          { lk = LinkIndex[i];
            printf("  %d: [%d,%d] row = %d\n",i,lk->min,lk->max,lk->row);
          }
      }
  }
#endif
}


/**** DRAWING ROUTINES ****/

/* Return dimensions to interface so it can setup scrollbars */

void GetAssemblyDims(int *min, int *max, int *rows, int *smallest)
{ *min = MinPnt;
  *max = MaxPnt;
  *rows = Rows + Links;
  *smallest = SmallestLen;
}

/* Routines to establish current master to pixel coordinate mapping */

#define XMAP(x)  (((x) - Xoff)*Xfact)
#define YMAP(y)  (((y) - Yoff)*Yfact)

#define XIMAP(x)  ((x)/Xfact + Xoff)
#define YIMAP(y)  ((y)/Yfact + Yoff)

static long   BackColor, TextColor, SeldColor;

static double Xfact, Yfact, Xoff, Yoff;

static int    Can_xmin, Can_xmax, Can_ymin, Can_ymax;
static int    Can_wide, Can_high;

void SetAssemblyMap(int bplow, int bphgh, int rowlow, int rowhgh,
                    int wide, int high, long back, long text, long match)
{ BackColor = back;
  TextColor = text;
  SeldColor = match;
  Can_xmin  = bplow;
  Can_xmax  = bphgh;
  Can_ymin  = rowlow - Links;
  Can_ymax  = rowhgh - Links;
  Can_wide  = wide;
  Can_high  = high;

  high -= RULER_HEIGHT;
  Xfact = wide/((bphgh-bplow)+1.);
  Yfact = high/((rowhgh-rowlow)+1.);
  if ((rowhgh-rowlow)+1 > (Rows+Links))
    Yoff = rowlow - ((((rowhgh-rowlow)-(Rows+Links))/2.) + 1);
  else
    Yoff = rowlow - .5;
  Yoff -= Links;
  Xoff  = bplow;
}

#define BINARY_SEARCH(Size,Predicate,Answer)	\
{ int lft, rgt, mid;				\
  lft = Size;					\
  rgt = 0;					\
  while (rgt < lft)				\
    { mid = (lft+rgt)/2;			\
      if (Predicate)				\
        lft = mid;				\
      else					\
        rgt = mid+1;				\
    }						\
  Answer = lft;					\
}

static int pick_obj = -1;
static int pick_button;
static int pick_line;

void PickRelease(MT_OBJECT *frame)
{ int s;

  s = pick_obj;
  if (s < 0) return;
  pick_obj = -1;
  
  if (!pick_line)
    { Link *lk;
      Line *ln, **kp;

      mt_begin_clip_to(frame);
      lk = LinkIndex[s];
      for (kp = lk->beg; kp < lk->end; kp++)
        { Style *at;
          int64   x1, x2, y;
          int   marks; 
          int64 *xp;

          ln = *kp;
          y  = YMAP(ln->row);
          xp = ln->beg;
          x1 = XMAP(*xp);
          marks = 0;
          if (ln->rowcon >= -1)
            for (xp += 2; xp < ln->end; xp += 2)
              { if (*xp < 0)
                  { marks = 1;
                    continue;
                  }
                at = StyleIndex[xp[-1]]; 
                if ((at->chord & LINESTYLE) == 0) continue;
                mt_set_color(at->color);
                x2 = XMAP(*xp);
                if (x1 < 0) x1 = 0;
                if (x2 > Can_wide) x2 = Can_wide;
                mt_draw_line(frame,x1,y,x2,y,at->thick,at->style);
                x1 = x2;
              }
          if (marks)
            { xp = ln->beg;
              for (xp += 2; xp < ln->end; xp += 2)
                { if (*xp >= 0) continue;
                  at = StyleIndex[xp[-1]]; 
                  if ((at->chord & MARKSTYLE) == 0) continue;
                  mt_set_color(at->color);
                  x2 = XMAP(-(*xp+1));
                  if (0 <= x2 && x2 <= Can_wide)
                    mt_draw_line(frame,x2,y-3,x2,y+3,at->thick,at->style);
                }
            }
        }
      mt_end_clip_to(frame);
    }
}

static int RoundInt(double c)
{ if (c < 0.)
    return ( - ((int) (.5 - c)) );
  else
    return ( (int) (.5 + c) );
}

#define PICK_NEAREST 0
#define PICK_SEGMENT 1
#define PICK_PIECES  2

char *PickAssembly(MT_OBJECT *frame, int64 x, int64 y, int v, int mode)
{ int64 bpzone;
  static char result[100];

  if (pick_obj >= 0) return (NULL);

  x = RoundInt(XIMAP(x));
  y = RoundInt(YIMAP(y));

  bpzone = PICK_THRESH/Xfact;

  { int s, d, i, ans;

    x += bpzone;
    if (y < 0)
      { if (NumLinks > 0)
          BINARY_SEARCH(NumLinks,LinkIndex[mid]->min > x,ans)
        else
          return (NULL);
        x -= bpzone;
        s = -1;
        d = bpzone + 1;
        for (i = ans-1; i >= 0; i--)
          { if (LinkIndex[i]->extent < x-bpzone) break;
            if (LinkIndex[i]->row == y)
              { if (LinkIndex[i]->max < x)
                  { if (d > x - LinkIndex[i]->max)
                      { s = i;
                        d = x - LinkIndex[i]->max;
                      }
                  }
                else if (LinkIndex[i]->min > x)
                  { if (d > LinkIndex[i]->min - x)
                      { s = i;
                        d = LinkIndex[i]->min - x;
                      }
                  }
                else
                  { s = i;
                    d = 0;
                    break;
                  }
              }
          }
        if (s < 0) return (NULL);

        pick_obj    = s;
        pick_line   = 0;
        pick_button = v;

        if (v <= 1)
          { Link *lk;
            Line *ln, **kp;
            int  slow, shgh;

            mt_begin_clip_to(frame);
            lk = LinkIndex[s];
            slow = (*(lk->beg))->beg[0];
            shgh = (*(lk->beg))->end[-1];
            for (kp = lk->beg; kp < lk->end; kp++)
              { Style *at;
                int64   x1, x2, y;
                int   marks;
                int64 *xp;

                ln = *kp;
                y  = YMAP(ln->row);
                xp = ln->beg;
                if (*xp < slow) slow = *xp;
                x1 = XMAP(*xp);
                marks = 0;
                for (xp += 2; xp < ln->end; xp += 2)
                  { if (*xp < 0)
                      { marks = 1;
                        continue;
                      }
                    at = StyleIndex[xp[-1]]; 
                    if ((at->chord & LINESTYLE) == 0) continue;
                    mt_set_color(SeldColor);
                    x2 = XMAP(*xp);
                    if (x1 < 0) x1 = 0;
                    if (x2 > Can_wide) x2 = Can_wide;
                    mt_draw_line(frame,x1,y,x2,y,at->thick,at->style);
                    x1 = x2;
                  }
                if (xp[-2] > shgh) shgh = xp[-2];
                if (marks)
                  { xp = ln->beg;
                    for (xp += 2; xp < ln->end; xp += 2)
                      { if (*xp >= 0) continue;
                        at = StyleIndex[xp[-1]]; 
                        if ((at->chord & MARKSTYLE) == 0) continue;
                        mt_set_color(at->color);
                        x2 = XMAP(-(*xp+1));
                        if (0 <= x2 && x2 <= Can_wide)
                          mt_draw_line(frame,x2,y-3,x2,y+3,at->thick,at->style);
                      }
                  }
              }
            mt_end_clip_to(frame);
            sprintf(result,"%d pieces spanning [%d,%d] len = %d",
                    (int) (lk->end-lk->beg),slow,shgh,shgh-slow);
            return (result);
          }
        else
          return (LinkIndex[s]->label);
      }
    else
      { Line *ln;
        int64   min, max;

        BINARY_SEARCH(NumLines,LineIndex[mid]->beg[0] > x,ans)
        x -= bpzone;
        s = -1;
        d = bpzone + 1;
        for (i = ans-1; i >= 0; i--)
          { ln = LineIndex[i];
            if (ln->extent < x-bpzone) break;
            if (ln->row == y)
              { if (ln->end[-1] < x)
                  { if (d > x - ln->end[-1])
                      { s = i;
                        d = x - ln->end[-1];
                      }
                  }
                else if (ln->beg[0] > x)
                  { if (d > ln->beg[0] - x)
                      { s = i;
                        d = ln->beg[0] - x;
                      }
                  }
                else
                  { s = i;
                    d = 0;
                    break;
                  }
              }
          }
        if (s < 0) return (NULL);

        pick_obj    = s;
        pick_line   = 1;
        pick_button = v;

        ln = LineIndex[s]; 
        min = ln->beg[0];
        if (mode != PICK_PIECES && v <= 1)
          { int64 *xp;

            for (xp = ln->beg+2; xp < ln->end-1; xp += 2)
              if (*xp < 0)
                { max = -(*xp+1);
                  if (mode == PICK_NEAREST &&
                      max - bpzone <= x && x <= max + bpzone)
                    { min = max;
                      break;
                    }
                }
              else
                { max = *xp;
                  if (x <= max) 
                    break;
                  else
                    min = max;
                }
            max = *xp;
          }
        else
          max = ln->end[-1];

        if (v <= 1)
          { if (max < 0)
              sprintf(result,"Mark of %s at " F_S64,ln->idnam,min);
            else if (mode != PICK_PIECES)
              sprintf(result,"Segment %s: [" F_S64 "," F_S64 "] len = " F_S64,
                             ln->idnam,min,max,max-min);
            else
              sprintf(result,"Piece of %s: [" F_S64 "," F_S64 "] len = " F_S64,
                             ln->idnam,min,max,max-min);
            return (result);
          }
        else
          return (LineIndex[s]->label);
      }
  }
}

static int TickRounder(double value, int *order, int *zeros)
{ int z, l, o, i;
  double v;

  if (value < 1.) value = 1.;
  z = log10(value);
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

/* Draw the window [bplow,bphgh] x [rowlow,rowhgh] of the layout in
   the given frame which is wide x high pixels in size with origin (0,0)
   at upper left.                                                       */

void DrawAssembly(MT_OBJECT *frame)
{ int ltop, ktop, atop;

  mt_begin_clip_to(frame);

  { int ticks, zeros, lead, place;
    char  label[20];
    char *letter[] = { "", "K", "M", "G", "T" };

    lead = TickRounder( (1.*(Can_xmax-Can_xmin)) / (Can_wide/50),
                        &place,&zeros);

    mt_set_color(TextColor);

    { int x1, x2;

      if (Can_xmin < MinPnt)
        x1 = XMAP(MinPnt);
      else
        x1 = 0;
      if (Can_xmax > MaxPnt)
        x2 = XMAP(MaxPnt);
      else
        x2 = Can_wide-1;
      mt_draw_line(frame,x1,Can_high-18,x2,Can_high-18,0,0);
    }

    ticks = (((Can_xmin-1) / lead) + 1)*lead;
    if (ticks < 0) ticks = 0;
    while (ticks <= Can_xmax && ticks <= MaxPnt)
      { int x1, w, h, b;

        x1 = XMAP(ticks);
        if (zeros % 3 == 2)
          sprintf(label,"%d.%d%s",ticks/(place*10),(ticks/place)%10,
                                  letter[zeros/3+1]);
        else if (zeros % 3 == 1)
          sprintf(label,"%d%s",ticks/(place/10),letter[zeros/3]);
        else
          sprintf(label,"%d%s",ticks/place,letter[zeros/3]);
        mt_string_size(label,&w,&h,&b);
        if (x1 >= w/2 && x1 <= Can_wide-(w/2))
          { mt_draw_text(frame,x1-(w/2),Can_high-2,label);
            mt_draw_line(frame,x1,Can_high-17,x1,Can_high-13,2,0);
          }
        ticks += lead;
      }
  }

  BINARY_SEARCH(NumLines,LineIndex[mid]->beg[0] > Can_xmax,ltop)

  if (NumLinks > 0)
    BINARY_SEARCH(NumLinks,LinkIndex[mid]->min > Can_xmax,ktop)
  else
    ktop = 0;

  if (NumAttach > 0)
    BINARY_SEARCH(NumAttach,AttachIndex[mid]->attach > Can_xmax,atop)
  else
    atop = 0;

  { int i, high;
    int drops;

    /* Draw drops from linker to line segment */

    high  = Can_high - RULER_HEIGHT;
    drops = 0;
    for (i = atop-1; i >= 0; i--)
      { Line  *ln;
        Style *at;
        int   y1, y2, x;

        ln = AttachIndex[i];
        if (ln->attach < Can_xmin) break;
        at = StyleIndex[ln->high->att];
        if ((at->chord & LINKSTYLE) == 0) continue;
        x  = XMAP(ln->attach);
        y2 = YMAP(ln->high->row);
        y1 = YMAP(ln->row);
        if (y1 > high) y1 = high;
        if (y2 <= y1)
          { mt_set_color(at->color);
            mt_draw_line(frame,x,y2,x,y1,at->thick,at->style);
            drops = 1;
          }
      }

    /* Draw oversized black silhouette around line segments to disconnect
       drop lines, and then restore connections where necessary           */

    if (drops)
      for (i = ltop-1; i >= 0; i--) 
        { Line  *ln; 
          Style *at;
          int64   x1, x2, y;
          int   blot; 
          int64 *xp;
  
          ln = LineIndex[i];
          if (ln->extent < Can_xmin) break;
          if (ln->row >= Can_ymin && ln->row <= Can_ymax)
            { y  = YMAP(ln->row);
              xp = ln->beg;
              x1 = XMAP(*xp);
              blot = 0;
              mt_set_color(BackColor);
              for (xp += 2; xp < ln->end; xp += 2)
  	        { at = StyleIndex[xp[-1]]; 
                  if (*xp < 0)
                    { if ((at->chord & MARKSTYLE) == 0) continue;
                      x2 = XMAP(-(*xp+1));
                      if (0 <= x2 && x2 <= Can_wide)
                        { mt_draw_line(frame,x2,y-4,x2,y+4,at->thick+2,0);
                          if (blot < 2) blot = 2;
                        }
                    }
                  else
                    { if ((at->chord & LINESTYLE) == 0) continue;
                      x2 = XMAP(*xp);
                      if (x1 < 0) x1 = 0;
                      if (x2 > Can_wide) x2 = Can_wide;
                      if (x1 <= x2)
                        mt_draw_line(frame,x1-1,y,x2+1,y,at->thick+4,0);
                      x1 = x2;
                      if (blot < at->thick/2) blot = at->thick/2;
                    }
                }
              if (ln->attach >= 0)
                { x2 = XMAP(ln->attach);
                  at = StyleIndex[ln->high->att];
                  if ((at->chord & LINKSTYLE) != 0)
                    { mt_set_color(at->color);
                      mt_draw_line(frame,x2,y-(blot+2),x2,y,
                                   at->thick,at->style);
                    }
                }
            }
        }

    /* Draw line segments and marks */

    for (i = ltop-1; i >= 0; i--) 
      { Line  *ln; 
        Style *at;
        int64   x1, x2, y;
        int   marks; 
        int64 *xp;

        ln = LineIndex[i];
        if (ln->extent < Can_xmin) break;
        if (ln->row >= Can_ymin && ln->row <= Can_ymax)
          { y  = YMAP(ln->row);
            xp = ln->beg;
            x1 = XMAP(*xp);
            marks = 0;
            for (xp += 2; xp < ln->end; xp += 2)
              { if (*xp < 0)
                  { marks = 1;
                    continue;
                  }
                at = StyleIndex[xp[-1]]; 
                x2 = XMAP(*xp);
                if ((at->chord & LINESTYLE) != 0)
                  { mt_set_color(at->color);
                    if (x1 < 0) x1 = 0;
                    if (x2 > Can_wide) x2 = Can_wide;
                    if (x1 <= x2)
                      mt_draw_line(frame,x1,y,x2,y,at->thick,at->style);
                  }
                x1 = x2;
              }
            if (marks)
              { xp = ln->beg;
                for (xp += 2; xp < ln->end; xp += 2)
                  { if (*xp >= 0) continue;
                    at = StyleIndex[xp[-1]]; 
                    if ((at->chord & MARKSTYLE) == 0) continue;
                    mt_set_color(at->color);
                    x2 = XMAP(-(*xp+1));
                    if (0 <= x2 && x2 <= Can_wide)
                      mt_draw_line(frame,x2,y-3,x2,y+3,at->thick,at->style);
                  }
              }
          }
      }

    /* Draw oversized black silhouettes around all linker horizontal bars */

    mt_set_color(BackColor);
    for (i = ktop-1; i >= 0; i--)
      { Link  *lk; 
        Style *at;
        int   x1, x2, y;

        lk = LinkIndex[i];
        if (lk->extent < Can_xmin) break;
        at = StyleIndex[lk->att];
        if ((at->chord & LINKSTYLE) == 0) continue;
        if (lk->row >= Can_ymin && lk->row <= Can_ymax)
          { y  = YMAP(lk->row);
            x1 = XMAP(lk->min);
            x2 = XMAP(lk->max);
            if (x1 < 0) x1 = 0;
            if (x2 > Can_wide) x2 = Can_wide;
            if (x1 <= x2)
              mt_draw_line(frame,x1,y,x2+1,y,at->thick+4,0);
          }
      }

    /* Draw linker horizontal bars and reconnect to drops */

    for (i = ktop-1; i >= 0; i--)
      { Link  *lk; 
        Line  **kp;
        Style *at;
        int64   x1, x2, y, t;

        lk = LinkIndex[i];
        if (lk->extent < Can_xmin) break;
        at = StyleIndex[lk->att];
        if ((at->chord & LINKSTYLE) == 0) continue;
        if (lk->row >= Can_ymin && lk->row <= Can_ymax)
          { y  = YMAP(lk->row);
            x1 = XMAP(lk->min);
            x2 = XMAP(lk->max);
            if (x1 < 0) x1 = 0;
            if (x2 > Can_wide) x2 = Can_wide;
            mt_set_color(at->color);
            if (x1 <= x2)
              mt_draw_line(frame,x1,y,x2+1,y,at->thick,at->style);
            t = at->thick/2;
            for (kp = lk->beg; kp < lk->end; kp++)
              { int x;
                int yl, yh;

                x = XMAP((*kp)->attach);
                if (x < 0 || x > Can_wide) continue;
                if ((*kp)->high != lk)
                  yl = y - (t + 2);  
                else
                  yl = y - (t);
                yh = y + (t + 2); 
                mt_draw_line(frame,x,yl,x,yh+1,at->thick,at->style);
              }
          }
      }
  }

  mt_end_clip_to(frame);
}

/**** QUERY PARSING ROUTINES ****/

/* Recursive Decent Parser */

typedef struct PTag
  { int          lab;
    struct PTag *lft;
    struct PTag *rgt;
  } PTree;

typedef struct {
    int   numpts;
    int   maxdep;
    int   initpt;
    int  *datapts;
} CoverPacket;

typedef struct {
    int   numobj;
    int  *objects;
    int  *segment;
} ObjectPacket;

typedef struct {
    char *refname;
    CoverPacket *refcover;
    ObjectPacket *refobject;
} result;

static char *Scan;
static int   ParseError;
static char *ErrorString;

static result *QTable;
static int     QLength;

#ifdef DEBUG_PARSE
static int   Pnodes;
#endif

static void SkipBlanks(void)
{ while (isspace(*Scan))
    Scan += 1;
}

static PTree *Node(int label, PTree *left, PTree *right)
{ PTree *p;

  p = (PTree *) malloc(sizeof(PTree));
  p->lab = label;
  p->lft = left;
  p->rgt = right;
#ifdef DEBUG_PARSE
  Pnodes += 1;
#endif
  return (p);
}

static PTree *GetRegExp(int iflag)
{ char  *balp, *error;
  int    posn;
  regexp mach;

  if (*Scan != '"')
    { ErrorString = "Expecting regular expression";
      ParseError = 1;
      return (NULL);
    }
  Scan += 1;
  for (balp = Scan; *balp != '"'; balp += 1)
    { if (*balp == '\\')
        balp += 1; 
      if (*balp == '\0')
        { ErrorString = "Unterminated regular expression";
          ParseError = 1;
          Scan = balp;
          return (NULL);
        }
    }
  *balp = '\0';
  mach = re_parse(iflag,Scan,&error,&posn);
  *balp = '"';
  if (mach == NULL)
    { ErrorString = error;
      ParseError = 1;
      Scan += posn; 
      return (NULL);
    }
  Scan = balp+1;
  SkipBlanks();
  return ((PTree *) mach);
}

static PTree *GetName(void)
{ int a, id;
  char *name;

  name = Scan;
  while (!isspace(*Scan) && *Scan != '\0')
    Scan += 1;
  a = *Scan;
  *Scan = '\0';
  id = 0;
  id = HashLookup(name,NameTable,&id);
  *Scan = a;
  if (id < 0)
    { ErrorString = "Don't recognize name";
      ParseError = 1;
      Scan -= 1;
      return (NULL);
    }
  SkipBlanks();
  return (Node(NameToStyle[id],NULL,NULL));
}

static PTree *GetInt(void)
{ int   num;
  char *end;

  num = strtol(Scan,&end,10);
  if (Scan == end)
    { ErrorString = "Expecting an integer";
      ParseError = 1;
      return (NULL);
    }
  Scan = end;
  SkipBlanks();
  return (Node(num,NULL,NULL));
}

void QueryTableSize(int size)
{ QLength = 0;
  QTable = (result *) malloc(sizeof(result)*size);
}

void QueryTableEntry(char *name, CoverPacket *cov, ObjectPacket *obj)
{ QTable[QLength].refname = name;
  QTable[QLength].refcover = cov;
  QTable[QLength].refobject = obj;
  QLength += 1;
}

void QueryTableFree(void)
{ free(QTable); }

static PTree *GetResult(int obqry)
{ char *name;
  int   a, i;

  name = ++Scan;
  while (!isspace(*Scan) && *Scan != '\0')
    Scan += 1;
  a = *Scan;
  *Scan = '\0';
  for (i = 0; i < QLength; i++)
    if (strcmp(QTable[i].refname,name) == 0)
      { *Scan = a;
        if (obqry && QTable[i].refcover != NULL)
          { ErrorString = "Reference is of wrong type";
            Scan -= 1;
            ParseError = 1;
            return (NULL);
	  }
        SkipBlanks();
        return (Node('$',Node(i,NULL,NULL),NULL));
      }
  *Scan = a;
  ErrorString = "Reference is undefined";
  Scan -= 1;
  ParseError = 1;
  return (NULL);
}

static PTree *GetOrClause(int);

static PTree *GetBaseClause(int obqry)
{ PTree *p;

  if (*Scan == '(')
    { Scan += 1;
      SkipBlanks();
      p = GetOrClause(obqry);
      if (ParseError) return (p);
      if (*Scan != ')')
        { ErrorString = "Closing paren is missing";
          ParseError = 1;
          return (p);
        }
      Scan += 1;
      SkipBlanks();
    }
  else if (*Scan == '$')
    p = GetResult(obqry);
  else
    p = GetName();
  return (p);
}

static PTree *GetSufClause(int obqry)
{ PTree *p;
  int c, nocase;

  p = GetBaseClause(obqry);
  while (!ParseError && (*Scan == '_' || *Scan == '>' ||
                         *Scan == '<' || *Scan == '+' || *Scan == '~'))
    { c = *Scan;
      if (obqry)
        { if (c == '_')
            { ErrorString = "Cannot have _ in @-expression";
              ParseError = 1;
              return (p);
            }
          else if (c == '+')
            { ErrorString = "Cannot have + in @-expression";
              ParseError = 1;
              return (p);
            }
        }
      else
        { if (c == '~')
            { ErrorString = "Cannot have ~ outside an @-expression";
              ParseError = 1;
              return (p);
            }
        }
      Scan += 1;
      if (c == '~')
        { if (*Scan == 'n')
            { c = 'n';
              Scan += 1;
            }
          else if (*Scan == 'c')
            { c = 'c';
              Scan += 1;
            }
          else
            c = 'b';
          nocase = 0;
          if (*Scan == 'i')
            { nocase = 1;
              Scan += 1;
            }
          SkipBlanks();
          p = Node(c,p,GetRegExp(nocase));
        }
      else
        { SkipBlanks();  // think this controls the length in < > expressions
	  {
	    PTree *n=GetInt(); 
	    n->lab = n->lab/COORDSCALE;
	    p = Node(c,p,n);
	  }
        }
    }
  return (p);
}

static PTree *GetPreClause(int obqry)
{ PTree *p;

  if (*Scan == '!')
    { Scan += 1;
      SkipBlanks();
      p = Node('!',GetPreClause(obqry),NULL);
    }
  else if (*Scan == '@')
    { if (obqry)
        { ErrorString = "@'s cannot nest";
          ParseError = 1;
          return (NULL);
        }
      Scan += 1;
      SkipBlanks();
      p = Node('@',GetPreClause(1),NULL);
    }
  else if (*Scan == '#')
    { if (!obqry)
        { ErrorString = "Cannot have # outside an @-expression";
          ParseError = 1;
          return (NULL);
        }
      Scan += 1;
      SkipBlanks();
      p = Node('#',GetPreClause(obqry),NULL);
    }
  else
    p = GetSufClause(obqry);
  return (p);
}

static PTree *GetIfClause(int obqry)
{ PTree *p;
  int code;

  p = GetPreClause(obqry);
  while (!ParseError && *Scan == '?')
    { Scan += 1;
      if (*Scan == 'i')
        code = '0';
      else if (*Scan == 'x')
        code = '1';
      else if (*Scan == 'c')
        code = '2';
      else
        { ErrorString = "Illegal ? operator";
          ParseError = 1;
          return (p);
        }
      Scan += 1;
      SkipBlanks();
      p = Node(code,p,GetPreClause(0));
    }
  return (p);
}

static PTree *GetAndClause(int obqry)
{ PTree *p;

  p = GetIfClause(obqry);
  while (!ParseError && *Scan == '&')
    { Scan += 1;
      SkipBlanks();
      p = Node('&',p,GetIfClause(obqry));
    }
  return (p);
}

static PTree *GetMinClause(int obqry)
{ PTree *p;
  char  *mark;

  p = GetAndClause(obqry);
  while (!ParseError && *Scan == '-')
    { mark = Scan;
      Scan += 1;
      SkipBlanks();
      if (isdigit(*Scan))
        { if (obqry)
            { ErrorString = "Cannot have + in @-expression";
              ParseError = 1;
              Scan = mark;
              return (p);
            }
          else
            p = Node('~',p,GetInt());
        }
      else
        p = Node('&',p,Node('!',GetAndClause(obqry),NULL));
    }
  return (p);
}

static PTree *GetOrClause(int obqry)
{ PTree *p;

  p = GetMinClause(obqry);
  while (!ParseError && *Scan == '|')
    { Scan += 1;
      SkipBlanks();
      p = Node('|',p,GetMinClause(obqry));
    }
  return (p);
}

#ifdef DEBUG_PARSE

static void RecTree(int level, PTree *p)
{ if (p->lft != NULL) RecTree(level+1,p->lft);
  if (p->lft == NULL)
    printf("%*s%d\n",3*level,"",p->lab);
  else
    printf("%*s%c\n",3*level,"",p->lab);
  if (p->lft != NULL && p->rgt != NULL)
    if (isalpha(p->lab))
      printf("%*sReg Exp\n",3*level+3,"");
    else
      RecTree(level+1,p->rgt);
}

static void ShowTree(PTree *p)
{ printf("Parse Tree:\n");
  RecTree(1,p);
  fflush(stdout);
}

#endif

static void FreeTree(PTree *p)
{ if (p->lft != NULL)
    { FreeTree(p->lft);
      if (p->rgt != NULL)
        { if (isalpha(p->lab))
            re_free((regexp) (p->rgt));
          else
            FreeTree(p->rgt);
        }
    }
  free(p);
#ifdef DEBUG_PARSE
  Pnodes -= 1;
#endif
}

/* Evaluate Query */

#ifdef DEBUG_QUERY

static void DumpCoverage(CoverPacket *cov)
{ int i, size, *data;

  size = cov->numpts;
  data = cov->datapts;
  printf("\nCOVERAGE LIST:\n");
  printf("  Size = %d  Max Depth = %d  InitVal = %d\n",
         size,cov->maxdep,cov->initpt);
  fflush(stdout);
  printf("  Values:\n");
  for (i = 0; i < size; i++)
    printf("    %8d\n",data[i]);
  fflush(stdout);
}

#endif

static int EXPAND_COMPARE(const void *l, const void *r)
{ int x, y, xa, ya, d;

  x = *((int *) l);
  y = *((int *) r);
  if (x < 0)
    xa = -(x+1);
  else
    xa = x;
  if (y < 0)
    ya = -(y+1);
  else
    ya = y;
  d = xa - ya;
  if (d != 0)
    return (d);
  else
    return (x - y);
}

static CoverPacket *Expander(CoverPacket *oprnd, int amount)
{ int size, *data, *stack;
  int val, npt, top, pos;
  int i, j;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  stack = (int *) malloc(sizeof(int)*oprnd->maxdep);
  
  amount = 2*amount;
  val = oprnd->initpt;
  npt = 0;
  top = -1;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (data[i] < 0)
        { pos = -(pos+1);
          val -= 1;
          if (top < 0) top = val;
          stack[val] = pos;
        }
      else
        { if (top >= val)
            { if (pos - stack[val] > amount)
                { for (j = val; j <= top; j++)
                    data[npt++] = -(stack[j]+1);
                  top = -1;
                  data[npt++] = pos;
                }
              else if (top == val)
                top = -1;
            }
          else
            data[npt++] = pos;
          val += 1;
        }
    }
  for (j = val; j <= top; j++)
    data[npt++] = -(stack[j]+1);
  data[npt] = 0;

  amount = amount/2;
  for (i = 0; i < npt; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1) + amount;
          if (pos > MaxPnt) pos = MaxPnt;
          data[i] = -(pos+1);
        }
      else
        { pos = pos - amount;
          if (pos < MinPnt) pos = MinPnt;
          data[i] = pos;
        }
    }
  qsort(data,npt,sizeof(int),EXPAND_COMPARE);

  oprnd->numpts = npt;
  free(stack);

#ifdef DEBUG_QUERY
  printf("Expander\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *Contractor(CoverPacket *oprnd, int amount)
{ int size, *data, *stack;
  int val, npt, cnf, pos;
  int i, mxd;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  stack = (int *) malloc(sizeof(int)*oprnd->maxdep);
  
  val = cnf = oprnd->initpt;
  npt = 0;
  mxd = val;
  for (i = 0; i < val; i++)
    stack[i] = - 2*amount;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (data[i] < 0)
        { pos = -(pos+1) - amount;
          val -= 1;
          if (pos > stack[val])
            { while (cnf <= val)
                data[npt++] = stack[cnf++];
              if (pos < MinPnt) pos = MinPnt;
              data[npt++] = -(pos+1);
              if (mxd < cnf) mxd = cnf;
            }
        }
      else
        { stack[val] = pos + amount;
          if (cnf > val) cnf = val;
          val += 1;
        }
    }
  while (cnf <= val)
    { pos = stack[cnf++];
      if (pos > MaxPnt) pos = MaxPnt;
      data[npt++] = pos;
    }
  if (mxd < cnf) mxd = cnf;
  data[npt] = 0;

  oprnd->numpts = npt;
  free(stack);

#ifdef DEBUG_QUERY
  printf("Contractor\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *CovContains(CoverPacket *oprnd, CoverPacket *thresh)
{ int size, *data, *stack;
  int val, pos, cnf;
  int npt, ndp;
  int i;
  int tize, *tata;
  int tpt, tpos, tsoon;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  stack = (int *) malloc(sizeof(int)*oprnd->maxdep);

  tize = thresh->numpts;
  tata = thresh->datapts;
  tpos = data[size-1];
  if (tpos < 0) tpos = -(tpos+1);
  tata[tize] = tpos+1;
  if (thresh->initpt & tize > 0)
    tpt = 1;
  else
    tpt  = 0;
  tpos = tata[tpt];

  npt = cnf = 0;
  val = ndp = oprnd->initpt;
  for (i = 0; i < val; i++)
    stack[i] = -1;
  tsoon = -2;

  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          while (tpos <= pos)
            { tpos = tata[++tpt];
              if (tpos < 0)
                tpos = -(tpos+1);
              else
                tsoon = tata[tpt-2];
            }
          while (cnf < val && stack[cnf] <= tsoon)
            { if (stack[cnf] >= 0)
                data[npt++] = stack[cnf];
              cnf += 1;
            }
          if (cnf == val)
            { data[npt++] = -(pos+1);
              cnf -= 1;
            }
          val -= 1;
        }
      else
        { stack[val++] = pos;
          if (val > ndp) ndp = val;
        } 
    }
  while (tpt < tize-1)
    { tpos = tata[++tpt];
      if (tpos < 0)
        tpos = -(tpos+1);
      else
        tsoon = tata[tpt-2];
    }
  while (cnf < val && stack[cnf] <= tsoon)
    { if (stack[cnf] >= 0)
        data[npt++] = stack[cnf];
      cnf += 1;
    }
  data[npt] = 0;
  oprnd->numpts  = npt;
  oprnd->maxdep  = ndp;
  free(stack);

#ifdef DEBUG_QUERY
  printf("Cover Contains\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *CovInter(CoverPacket *oprnd, CoverPacket *thresh)
{ int size, *data, *stack;
  int val, pos, cnf;
  int npt, ndp;
  int i;
  int tize, *tata;
  int tpt, tpos, tlpos, tlast;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  stack = (int *) malloc(sizeof(int)*oprnd->maxdep);

  tize = thresh->numpts;
  tata = thresh->datapts;
  tpos = data[size-1];
  if (tpos < 0) tpos = -(tpos+1);
  tata[tize] = tpos+1;
  tlpos = tlast = -2;
  tpos = tata[tpt=0];
  if (tpos < 0) tpos = -(tpos+1);

  npt = cnf = 0;
  val = ndp = oprnd->initpt;
  for (i = 0; i < val; i++)
    stack[i] = -1;

  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          while (tpos < pos)
            { tlpos = tlast = tata[tpt];
              if (tlpos < 0) tlpos = -(tlpos+1);
              tpos = tata[++tpt];
              if (tpos < 0) tpos = -(tpos+1);
            }
          if (tlast > 0)
            while (cnf < val)
              { if (stack[cnf] >= 0)
                  data[npt++] = stack[cnf];
                cnf += 1;
              }
          else
            while (cnf < val && stack[cnf] < tlpos)
              { if (stack[cnf] >= 0)
                  data[npt++] = stack[cnf];
                cnf += 1;
              }
          if (cnf == val)
            { data[npt++] = -(pos+1);
              cnf -= 1;
            }
          val -= 1;
        }
      else
        { stack[val++] = pos;
          if (val > ndp) ndp = val;
        } 
    }
  tpt = tize-1;
  tlpos = tlast = tata[tpt];
  if (tlpos < 0) tlpos = -(tlpos+1);
  if (tlast > 0)
    while (cnf < val)
      { if (stack[cnf] >= 0)
          data[npt++] = stack[cnf];
        cnf += 1;
      }
  else
    while (cnf < val && stack[cnf] < tlpos)
      { if (stack[cnf] >= 0)
          data[npt++] = stack[cnf];
        cnf += 1;
      }
  data[npt] = 0;
  oprnd->numpts  = npt;
  oprnd->maxdep  = ndp;
  free(stack);

#ifdef DEBUG_QUERY
  printf("Cover Contains\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *CovIncludes(CoverPacket *oprnd, CoverPacket *thresh)
{ int size, *data;
  int val, pos;
  int npt, ndp;
  int sync, last;
  int i, j;
  int tize, *tata;
  int tpt, tpos, tlpos;

  size  = oprnd->numpts;
  data  = oprnd->datapts;

  tize = thresh->numpts;
  tata = thresh->datapts;
  tpos = data[size-1];
  if (tpos < 0) tpos = -(tpos+1);
  tata[tize] = tpos+1;
  tpos = tpt = -1;
  last = -2;

  npt = 0;
  ndp = 0;
  val = oprnd->initpt;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          val -= 1;
          if (val == 0)
            { while (tpos < pos)
                { tlpos = tpos;
                  tpos = tata[++tpt];
                  if (tpos <= 0) tpos = -(tpos+1);
                }
              if (tata[tpt] < 0 && tlpos <= last)
                { for (j = sync; j <= i; j++)
                    { data[npt++] = data[j];
                      if (data[j] < 0)
                        val -= 1;
                      else
                        val += 1;
                      if (val > ndp) ndp = val;
                    }
                }
            }
        }
      else
        { if (val == 0)
            { sync = i; last = pos; }
          val += 1;
        } 
    }
  data[npt] = 0;
  oprnd->numpts  = npt;
  oprnd->maxdep  = ndp;

#ifdef DEBUG_QUERY
  printf("Cover Contained Within\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *GreaterThan(CoverPacket *oprnd, int thresh)
{ int size, *data, *stack;
  int val, pos, cnf;
  int npt, ndp;
  int i;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  stack = (int *) malloc(sizeof(int)*oprnd->maxdep);
  npt = 0;
  val = ndp = cnf = oprnd->initpt;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          while (cnf < val && pos - stack[cnf] > thresh)
            data[npt++] = stack[cnf++];
          if (cnf == val)
            { data[npt++] = -(pos+1);
              cnf -= 1;
            }
          val -= 1;
        }
      else
        { stack[val++] = pos;
          if (val > ndp) ndp = val;
        } 
    }
  for (i = cnf; i < val; i++)
    data[npt++] = stack[i];
  data[npt] = 0;
  oprnd->numpts  = npt;
  oprnd->maxdep  = ndp;
  free(stack);

#ifdef DEBUG_QUERY
  printf("Greater Than\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *LessThan(CoverPacket *oprnd, int thresh)
{ int size, *data;
  int val, pos;
  int npt, ndp;
  int sync, last;
  int i, j;

  size  = oprnd->numpts;
  data  = oprnd->datapts;
  npt = 0;
  ndp = 0;
  val = oprnd->initpt;
  last = - thresh;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          val -= 1;
          if (val == 0 && pos - last < thresh)
            { for (j = sync; j <= i; j++)
                { data[npt++] = data[j];
                  if (data[j] < 0)
                    val -= 1;
                  else
                    val += 1;
                  if (val > ndp) ndp = val;
                }
            }
        }
      else
        { if (val == 0)
            { sync = i; last = pos; }
          val += 1;
        } 
    }
  data[npt] = 0;
  oprnd->numpts  = npt;
  oprnd->maxdep  = ndp;

#ifdef DEBUG_QUERY
  printf("Less Than\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *Negate(CoverPacket *oprnd)
{ int size, *data,*negdata;
  int val;
  int i;
  int icase=0;
  

  data = oprnd->datapts;
  size = oprnd->numpts;
  val  = oprnd->initpt;
  /*
    Check cases:
             |----       ----|  Negate takes two intervals to one
             |----    ----   |  Negate takes two intervals to two
             |   ----    ----|  Negate takes two intervals to two
             |  ----    ---- |  Negate takes two intervals to three
	     |               |  Negate takes no intervals to one
	     |---------------|  Negate takes one interval to none
  */

  if(size==0){
    icase = 6;
  } else {
    if ( data[0] == MinPnt ) {
      if ( data[size-1] == -(MaxPnt+1) ) {
	if(size==2){
	  icase = 5;
	} else {
	  icase = 1;
	}
      } else {
        icase = 2;
      }
    } else {
      if ( data[size-1] == -(MaxPnt+1) ) {
        icase = 3;
      } else {
        icase = 4;
      }
    }
  }

  switch ( icase ) {
  case 1:
    negdata = (int *) malloc((size-2)*sizeof(int));
    for (i = 0; i < size-2; i++)
      { negdata[i] = -(data[i+1]+1);
      }
    oprnd->numpts = size -2;
    break;
  case 2:
    negdata = data;
    for (i = 0; i < size-2; i++)
      { negdata[i] = -(data[i+1]+1);
      }
    negdata[size-2]=data[size-1];
    negdata[size - 1 ] = - (MaxPnt+1);
    oprnd->numpts = size;
    break;
  case 3:
    negdata = data;
    oprnd->numpts = size;
    for (i = size-1; i > 0; i--)
      { negdata[i] = -(data[i-1]+1);
      }
    negdata[0] = MinPnt;
    oprnd->datapts= negdata;
    break;
  case 4:
    negdata =  (int *) malloc((size+2)*sizeof(int));
    negdata[size] = - (data[size-1]+1);
    negdata[size+1] = - (MaxPnt+1);
    oprnd->numpts = size+2;
    for (i = size-1; i > 0; i--)
      { negdata[i] = -(data[i-1]+1);
      }
    negdata[0] = MinPnt;
    oprnd->datapts= negdata;
    break;
  case 5:
    negdata = data; // sloppy but harmless -- not sure we want to have it be null
    oprnd->numpts = 0;
    break;
  case 6:
    negdata = (int *) malloc(2*sizeof(int));
    negdata[0]= MinPnt;
    negdata[1]= -(MaxPnt+1);
    oprnd->numpts=2;
    break;
  default:
    break;
  }

  oprnd->datapts=negdata;

  if(icase!=5){
    oprnd->maxdep = 1;
  }else{
    oprnd->maxdep = 0;
  }

  if ( (icase == 1 || icase == 4 || icase == 6) && data!=NULL ) free(data);
  oprnd->maxdep = 1;
 
#ifdef DEBUG_QUERY
  printf("Negate\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *Divide(CoverPacket *oprnd, int divisor)
{ int size, *data;
  int pos, val;
  int cval, nval;
  int i, s, mxdp;

  data = oprnd->datapts;
  size = oprnd->numpts;
  val  = oprnd->initpt;
  cval = val / divisor;
  mxdp = cval;
  oprnd->initpt = cval;
  s = 0;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { pos = -(pos+1);
          val -= 1;
        }
      else
        val += 1;
      nval = val / divisor;
      if (nval < cval)
        data[s++] = -(pos+1);
      else if (nval > cval)
        data[s++] = pos;
      cval = nval;
      if (cval > mxdp) mxdp = cval;
    }
  oprnd->maxdep = mxdp;
  oprnd->numpts = s;

#ifdef DEBUG_QUERY
  printf("Divide\n");
  DumpCoverage(oprnd);
#endif

  return (oprnd);
}

static CoverPacket *Intersection(CoverPacket *left, CoverPacket *right)
{ int lsize, rsize;
  int *ldata, *rdata, *data;
  int lval, rval, nval, cval, mval, omval;
  int lpos, rpos, npos;
  int l, r, s, smax;

  lsize = left->numpts;
  rsize = right->numpts;
  lpos  = left->datapts[lsize-1];
  rpos  = right->datapts[rsize-1];
  if (lpos < 0) lpos = -(lpos+1);
  if (rpos < 0) rpos = -(rpos+1);

  if (rpos > lpos)
    { CoverPacket *temp;
      temp  = right;
      right = left;
      left  = temp;
    }

  lsize = left->numpts;
  rsize = right->numpts;
  ldata = left->datapts;
  rdata = right->datapts;
  lval  = left->initpt;
  rval  = right->initpt;

  if (lval < rval)
    omval = mval = cval = lval;
  else
    omval = mval = cval = rval;
  data = (int *) malloc(sizeof(int)*(lsize+rsize+1));
  left->initpt = cval;

  l = r = s = smax = 0;
  lpos = ldata[l];
  rpos = rdata[r];
  if (lpos < 0) lpos = -(lpos+1);
  if (rpos < 0) rpos = -(rpos+1);
  while (r < rsize)
    { 
      if (rpos > lpos)
        { 
          if (ldata[l] < 0)
            lval -= 1; 
          else
            lval += 1; 
          npos = lpos;
          lpos = ldata[++l];
          if (lpos < 0) lpos = -(lpos+1);
        }
      else
        { if (rdata[r] < 0)
            rval -= 1; 
          else
            rval += 1; 
          npos = rpos;
          rpos = rdata[++r];
          if (rpos < 0) rpos = -(rpos+1);
        }
      if (lval < rval)
        nval = lval;
      else
        nval = rval;
      if (nval < cval)
        data[s++] = -(npos+1);
      else if (nval > cval)
        data[s++] = npos;
      if ( s > 0 && (s+1)%2==1 ) {
        // check for 0 length interval
        int d0=data[s-2];
        int d1=data[s-1];
        if ( d0 < 0 ) d0 = -(d0+1);
          else d1 = -(d1+1);
        if ( d0 == d1 ) {
          //fprintf(stderr,"BARK (zero length coverage interval): %d,%d\n",data[s-2],data[s-1]);
          // back out the interval
//Following test not fully guaranteed -- Aaron and Karin had different
// proposals -- take your pick
// Karin:           
//	  if ( r < rsize -1 )  mval = omval;
// Aaron:
	  if(s==smax+1) mval = omval;
          s = s-2;
        }
      } 
      cval = nval;
      if (cval > mval) {
        //fprintf(stderr,"changing max depth from %d to %d\n",mval,cval);
        omval = mval; // archive this
        mval = cval;
      }
    }

  while (l < lsize)
    { if (ldata[l] < 0)
        lval -= 1;
      else
        lval += 1;
      npos = lpos;
      lpos = ldata[++l];
      if (lpos < 0) lpos = -(lpos+1);
      if (lval < rval)
        nval = lval;
      else
        nval = rval;
      if (nval < cval)
        data[s++] = -(npos+1);
      else if (nval > cval)
        data[s++] = npos;
      if ( s > 0 && (s+1)%2==1 ) {
        // check for 0 length interval
        int d0=data[s-2];
        int d1=data[s-1];
        if ( d0 < 0 ) d0 = -(d0+1);
          else d1 = -(d1+1);
        if ( d0 == d1 ) {
          // fprintf(stderr,"BARK (zero length coverage interval): %d,%d\n",data[s-2],data[s-1]);
          // back out the interval
//Following test not fully guaranteed -- Aaron and Karin had different
// proposals -- take your pick
// Karin:           
//        if ( l < lsize - 1 ) mval = omval;
// Aaron:
	  if(s==smax+1) mval = omval;
          s = s-2;
        }
      } 
      cval = nval;
      if (cval > mval) {
         omval = mval;
         mval = cval;
      }
    }
 
  data[s] = 0;
  left->numpts = s;
  left->maxdep = mval;
  left->datapts = data;
  free(ldata);
  free(rdata);
  free(right);

#ifdef DEBUG_QUERY
  printf("Intersection\n");
  DumpCoverage(left);
#endif

  return (left);
}


static CoverPacket *Union(CoverPacket *left, CoverPacket *right)
{ int lsize, rsize;
  int *ldata, *rdata, *data;
  int lval, rval, nval, cval, mval;
  int lpos, rpos, npos;
  int l, r, s;

  lsize = left->numpts;
  rsize = right->numpts;
  lpos  = left->datapts[lsize-1];
  rpos  = right->datapts[rsize-1];
  if (lpos < 0) lpos = -(lpos+1);
  if (rpos < 0) rpos = -(rpos+1);

  if (rpos > lpos)
    { CoverPacket *temp;
      temp  = right;
      right = left;
      left  = temp;
    }

  lsize = left->numpts;
  rsize = right->numpts;
  ldata = left->datapts;
  rdata = right->datapts;
  lval  = left->initpt;
  rval  = right->initpt;

  if (lval < rval)
    mval = cval = rval;
  else
    mval = cval = lval;
  data = (int *) malloc(sizeof(int)*(lsize+rsize+1));
  left->initpt = cval;

  l = r = s = 0;
  lpos = ldata[l];
  rpos = rdata[r];
  if (lpos < 0) lpos = -(lpos+1);
  if (rpos < 0) rpos = -(rpos+1);
  while (r < rsize)
    { if (rpos > lpos)
        { if (ldata[l] < 0)
            lval -= 1; 
          else
            lval += 1; 
          npos = lpos;
          lpos = ldata[++l];
          if (lpos < 0) lpos = -(lpos+1);
        }
      else
        { if (rdata[r] < 0)
            rval -= 1; 
          else
            rval += 1; 
          npos = rpos;
          rpos = rdata[++r];
          if (rpos < 0) rpos = -(rpos+1);
        }
      if (lval < rval)
        nval = rval;
      else
        nval = lval;
      if (nval < cval)
        data[s++] = -(npos+1);
      else if (nval > cval)
        data[s++] = npos;
      cval = nval;
      if (cval > mval) mval = cval;
    }

  while (l < lsize)
    { if (ldata[l] < 0)
        lval -= 1;
      else
        lval += 1;
      npos = lpos;
      lpos = ldata[++l];
      if (lpos < 0) lpos = -(lpos+1);
      if (lval < rval)
        nval = rval;
      else
        nval = lval;
      if (nval < cval)
        data[s++] = -(npos+1);
      else if (nval > cval)
        data[s++] = npos;
      cval = nval;
      if (cval > mval) mval = cval;
    }

  data[s] = 0;
  left->numpts = s;
  left->maxdep = mval;
  left->datapts = data;
  free(ldata);
  free(rdata);
  free(right);

#ifdef DEBUG_QUERY
  printf("Union\n");
  DumpCoverage(left);
#endif

  return (left);
}

/* POSSIBILITY:  Object operators -- instead of just coverage concept
     apply to segments & links.  Interconversion is still fuzzy, but
     would have &, |, !, >, <, and -.  The ops + ~ and _ are strictly
     coverage concepts.  Pure object concepts would be (a) label or
     comment matches r.e., (b) contained in an interval of a coverage
     set or contains an interval of a coverage set.
*/

static int ISA_idx, ISA_min, ISA_max, ISA_event;
static char *ISA_name, *ISA_label;

static int InsideSearch(int min, int max, CoverPacket *cover)
{ int size, *data;
  int l, r, m, v;

  size = cover->numpts;
  data = cover->datapts;

  l = 0;
  r = size;
  while (l < r)
    { m = (l+r) / 2;
      v = data[m];
      if (v < 0) v = -(v+1);
      if (v < max)
        l = m+1;
      else
        r = m;
    }
  if (r == 0) return (cover->initpt);
  v = data[r-1];
  return (v > 0 && v <= min);
}

static int ContainSearch(int min, int max, CoverPacket *cover)
{ int size, *data;
  int l, r, m, v;

  size = cover->numpts;
  data = cover->datapts;

  l = 0;
  r = size;
  while (l < r)
    { m = (l+r) / 2;
      v = data[m];
      if (v < 0) v = -(v+1);
      if (v < min)
        l = m+1;
      else
        r = m;
    }
  if (r+1 >= size) return (0);
  v = data[r+1];
  if (v < 0) return (-(v+1) <= max);
  if (r+2 >= size) return (0);
  v = data[r+2];
  return (-(v+1) <= max);
}

static int IntersectSearch(int min, int max, CoverPacket *cover)
{ int size, *data;
  int l, r, m, v;

  size = cover->numpts;
  data = cover->datapts;

  if (size == 0) return (cover->initpt);

  l = 0;
  r = size;
  while (l < r)
    { m = (l+r) / 2;
      v = data[m];
      if (v < 0) v = -(v+1);
      if (v < min)
        l = m+1;
      else
        r = m;
    }
  if (r >= size) return (data[r-1] > 0);
  v = data[r];
  if (v >= 0) return (v < max);
  v = -(v+1);
  if (v > min) return (1);
  if (r+1 >= size) return (0);
  v = data[r+1];
  return (v < max);
}

static int IsAnEval(PTree *p)
{ if (p->lft != NULL)
    switch (p->lab)
    { case '&':
        return (IsAnEval(p->lft) && IsAnEval(p->rgt));
      case '|':
        return (IsAnEval(p->lft) || IsAnEval(p->rgt));
      case '!':
        return ( ! IsAnEval(p->lft));
      case '>':
        return ( IsAnEval(p->lft) && (ISA_max - ISA_min > p->rgt->lab));
      case '<':
        return ( IsAnEval(p->lft) && (ISA_max - ISA_min < p->rgt->lab));
      case 'n':
        return ( IsAnEval(p->lft) && ISA_name != NULL &&
                      re_match((regexp) (p->rgt),ISA_name,0));
      case 'c':
        return ( IsAnEval(p->lft) && re_match((regexp) (p->rgt),ISA_label,0));
      case 'b':
        if ( IsAnEval(p->lft) )
          return (re_match((regexp) (p->rgt),ISA_label,0) ||
                  ISA_name != NULL && re_match((regexp) (p->rgt),ISA_name,0) );
        else
          return (0);
      case '$':
        if (((char *) (p->rgt))[ISA_event] != 0)
          return (1);
        else
          return (0);
      case '#':
        if (ISA_event >= 0 & ((char *) (p->lft))[ISA_event] != 0)
          return (1);
        else
          return (0);
      case '0':
        return (IsAnEval(p->lft) &&
                   InsideSearch(ISA_min,ISA_max,(CoverPacket *) (p->rgt)));
      case '1':
        return (IsAnEval(p->lft) &&
                   IntersectSearch(ISA_min,ISA_max,(CoverPacket *) (p->rgt)));
      case '2':
        return (IsAnEval(p->lft) &&
                   ContainSearch(ISA_min,ISA_max,(CoverPacket *) (p->rgt)));
    }
  else
    return (p->lab == ISA_idx);
}

static int IsAnInstance(PTree *p, int idx, int xmn, int xmx,
                           char *name, char *label, int evn)
{ ISA_idx   = idx;
  ISA_min   = xmn;
  ISA_max   = xmx;
  ISA_name  = name;
  ISA_label = label;
  ISA_event = evn;
  return (IsAnEval(p));
}

void FreeObjects(ObjectPacket *opk)
{ free(opk->objects);
  free(opk->segment);
  free(opk);
}

#ifdef DEBUG_QUERY

static void DumpObjects(ObjectPacket *opk)
{ int i, size, *ject, *segm;

  size = opk->numobj;
  ject = opk->objects;
  segm = opk->segment;
  printf("\nOBJECT LIST: Size = %d\n",size);
  fflush(stdout);
  printf("  Objects:\n");
  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      printf("    %8d link\n",ject[i]);
    else
      printf("    %8d line seg %d\n",ject[i],segm[i]);
  fflush(stdout);
}

#endif

static char *ObjectChord(ObjectPacket *oprnd)
{ int size, *ject, *segm;
  char *chord;
  int i;

  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  chord = (char *) malloc(sizeof(char)*((NumEvents-NumLines)+NumLinks));
  chord += NumLinks;

  for (i = -NumLinks; i < NumEvents-NumLines; i++)
    chord[i] = 0;

  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      chord[-ject[i]] = 1;
    else
      chord[LineIndex[ject[i]]->locale+segm[i]/2] = 1;

  return (chord);
}

static char *LinkExpand(ObjectPacket *oprnd)
{ int size, *ject, *segm;
  char *chord;
  int i, n, m, ev;

  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  chord = (char *) malloc(sizeof(char)*(NumEvents-NumLines));

  for (i = 0; i < NumEvents-NumLines; i++)
    chord[i] = 0;

  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      { Link *lk;
        Line *ln, **lp;

        lk = LinkIndex[ject[i]];
        for (lp = lk->beg; lp < lk->end; lp++)
          { ln = *lp;
            ev = ln->locale;
            for (m = 0; ln->beg + m < ln->end-1; m = n)
              { n = m+2;
                while (ln->beg[n] < 0)
                  n += 2;
                chord[ev+m/2] = 1;
              }
          }
      }

#ifdef DEBUG_QUERY
  printf("\nLink # chord\n");
  for (i = 0; i < NumEvents; i++)
    if (chord[i])
      printf("    %4d\n",i);
#endif

  FreeObjects(oprnd);

  return (chord);
}

static ObjectPacket *LeafEvent(PTree *);
static CoverPacket  *InterpTree(PTree *);

static void SetUpTree(PTree *p)
{ ObjectPacket *set;
  CoverPacket  *ref;

  if (p->lft != NULL)
    switch (p->lab)
    { case '&': case '|':
        SetUpTree(p->lft);
        SetUpTree(p->rgt);
        return;
      case '!': case '>': case '<':
      case 'n': case 'c': case 'b':
        SetUpTree(p->lft);
        return;
      case '$':
        p->rgt = (PTree *) ObjectChord(QTable[p->lft->lab].refobject);
        return;
      case '#':
        set = LeafEvent(p->lft);
        FreeTree(p->lft);
        p->lft = (PTree *) LinkExpand(set);
        return;
      case '0': case '1': case '2':
        SetUpTree(p->lft);
        ref = Negate(Negate(InterpTree(p->rgt)));
        FreeTree(p->rgt);
        p->rgt = (PTree *) ref;
        return;
    }
}

static void CleanUpTree(PTree *p)
{ CoverPacket *ref;
  char        *chord;

  if (p->lft != NULL)
    switch (p->lab)
    { case '&': case '|':
        CleanUpTree(p->lft);
        CleanUpTree(p->rgt);
        return;
      case '!': case '>': case '<':
      case 'n': case 'c': case 'b':
        CleanUpTree(p->lft);
        return;
      case '0': case '1': case '2':
        CleanUpTree(p->lft);
        ref = (CoverPacket *) (p->rgt);
        free(ref->datapts);
        free(ref);
        p->rgt = NULL;
        return;
      case '$':
        chord = (char *) (p->rgt);
        free(chord-NumLinks);
        p->rgt = NULL;
        return;
      case '#':
        chord = (char *) (p->lft);
        free(chord);
        p->lft = NULL;
        return;
    }
}

static int *LineCache, *LinkCache;
static int *QueryMap, *MarkupMap;

static void BuildColorMemory(void)
{ int i, m, ev;
  Line *ln;
  static int Executed = 0;

  if (Executed) return;
  Executed = 1;

  m = 0;
  for (i = 0; i < NumLines; i++)
    { ln = LineIndex[i];
      m += (ln->end - ln->beg) / 2;
    }
  LineCache = (int *) malloc(sizeof(int)*m);
  for (i = 0; i < NumLines; i++)
    { ln = LineIndex[i];
      ev = ln->locale;
      for (m = 1; ln->beg + m < ln->end; m += 2)
        LineCache[ev++] = ln->beg[m];
    }
  LinkCache = (int *) malloc(sizeof(int)*NumLinks);
  for (i = 0; i < NumLinks; i++)
    LinkCache[i] = LinkIndex[i]->att;

  QueryMap = (int *) malloc(sizeof(int)*(AllStyles+1));
  QueryMap += 1;
  MarkupMap = (int *) malloc(sizeof(int)*(AllStyles+1));
  MarkupMap += 1;
}

static int AllocateMarkup(int style, long color)
{ Style *sty, *bty;
  int   i, m;
  Line *ln;
  Link *lk;

  if (SlotsList == NULL)
    { for (i = NumStyles; i < AllStyles; i++)
        MarkupMap[i] = 0;
      for (i = 0; i < NumLines; i++)
        { ln = LineIndex[i];
          for (m = 1; ln->beg + m < ln->end; m += 2)
            MarkupMap[ln->beg[m]] = 1;
        }
      for (i = 0; i < NumLinks; i++)
        { lk = ExtentIndex[i];
          MarkupMap[lk->att] = 1;
        }
      for (i = NumStyles; i < AllStyles; i++)
        if (MarkupMap[i] == 0)
          { StyleIndex[i]->next = SlotsList;
            StyleIndex[i]->usage = 0;
            SlotsList = StyleIndex[i];
          }
    }

  if (SlotsList == NULL)
    { Style *cs;

      StyleIndex = (Style **) realloc(StyleIndex-1,
                                      sizeof(Style *)*(2*AllStyles+1));
      StyleIndex += 1;
      cs = (Style *) malloc(sizeof(Style)*AllStyles);
      for (i = AllStyles; i < 2*AllStyles; i++)
        { StyleIndex[i] = cs + (i-AllStyles);
          StyleIndex[i]->idnum = i;
          StyleIndex[i]->usage = 0;
        }
      for (i = AllStyles; i < 2*AllStyles-1; i += 1)
        StyleIndex[i]->next = StyleIndex[i+1];
      StyleIndex[i]->next = NULL;
      SlotsList = StyleIndex[AllStyles];
      AllStyles = 2*AllStyles;

      QueryMap = (int *) realloc(QueryMap-1,sizeof(int)*(AllStyles+1));
      QueryMap += 1;
      MarkupMap = (int *) realloc(MarkupMap-1,sizeof(int)*(AllStyles+1));
      MarkupMap += 1;
    }

  bty = StyleIndex[style];
  sty = SlotsList;
  SlotsList = sty->next;
  sty->next  = bty;
  sty->color = color;
  sty->thick = bty->thick;
  sty->style = bty->style;
  sty->usage = SHOWSTYLE;
  sty->chord = bty->chord;
  return (sty->idnum);
}

void SetQueryColor(ObjectPacket *oprnd, long color)
{
  int size, *ject, *segm;
  int i, n, x;
  Line *ln;
  Link *lk;

  BuildColorMemory();
  
  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  for (i = -1; i < NumStyles; i++)
    QueryMap[i] = -1;

  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      { lk = LinkIndex[ject[i]];
        x = lk->att;
        if (x >= NumStyles)
          x = LinkCache[ject[i]];
        if (QueryMap[x] < 0)
          QueryMap[x] = AllocateMarkup(x,color);
        lk->att = QueryMap[x];
      }
    else
      { ln = LineIndex[ject[i]];
        n = segm[i]+2;
        while (ln->beg[n] < 0)
          n += 2;
        n -= 1;
        x = ln->beg[n];
        if (x >= NumStyles)
          x = LineCache[ln->locale + n/2];
        if (QueryMap[x] < 0)
          QueryMap[x] = AllocateMarkup(x,color);
        ln->beg[n] = QueryMap[x];
      }
}

void ClearQueryColor(ObjectPacket *oprnd)
{ int size, *ject, *segm;
  int i, n;
  Line *ln;
  Link *lk;

  BuildColorMemory();
  
  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      { lk = LinkIndex[ject[i]];
        lk->att = LinkCache[ject[i]];
      }
    else
      { ln = LineIndex[ject[i]];
        n = segm[i]+2;
        while (ln->beg[n] < 0)
          n += 2;
        n -= 1;
        ln->beg[n] = LineCache[ln->locale + n/2];
      }
}

static ObjectPacket *LeafEvent(PTree *p)
{ ObjectPacket *rez;
  int *ject, *segm;
  int ev, i, s, n, m, ty;
  Line *ln;
  Link *lk;

  SetUpTree(p);

  s = 0;
  for (i = 0; i < NumLines; i++)
    { ln = LineIndex[i];
      ev = ln->locale;
      for (m = 0; ln->beg + m < ln->end-1; m = n)
        { n = m+2;
          while (ln->beg[n] < 0)
            n += 2;
          ty = ln->beg[n-1];
          if (ty >= NumStyles)
            ty = StyleIndex[ty]->next->idnum;
          if (IsAnInstance(p,ty,ln->beg[m],ln->beg[n],
                             ln->idnam,ln->label,ev+m/2))
            s += 1;
        }
    }
  for (i = 0; i < NumLinks; i++)
    { lk = LinkIndex[i];
      ty = lk->att;
      if (ty >= NumStyles)
        ty = StyleIndex[ty]->next->idnum;
      if (IsAnInstance(p,ty,lk->emin,lk->emax,NULL,lk->label,-(i+1)))
        s += 1;
    }

  ject = (int *) malloc(sizeof(int)*(s+1));
  segm = (int *) malloc(sizeof(int)*(s+1));

  s = 0;
  for (i = 0; i < NumLines; i++)
    { ln = LineIndex[i];
      ev = ln->locale;
      for (m = 0; ln->beg + m < ln->end-1; m = n)
        { n = m+2;
          while (ln->beg[n] < 0)
            n += 2;
          ty = ln->beg[n-1];
          if (ty >= NumStyles)
            ty = StyleIndex[ty]->next->idnum;
          if (IsAnInstance(p,ty,ln->beg[m],ln->beg[n],
                             ln->idnam,ln->label,ev+m/2))
            { ject[s] = i;
              segm[s++] = m;
            }
        }
    }
  for (i = 0; i < NumLinks; i++)
    { lk = LinkIndex[i];
      ty = lk->att;
      if (ty >= NumStyles)
        ty = StyleIndex[ty]->next->idnum;
      if (IsAnInstance(p,ty,lk->emin,lk->emax,NULL,lk->label,-(i+1)))
        { ject[s] = i;
          segm[s++] = -1;
        }
    }
  ject[s] = -1;

  rez = (ObjectPacket *) malloc(sizeof(ObjectPacket));
  rez->numobj  = s;
  rez->objects = ject;
  rez->segment = segm;

#ifdef DEBUG_QUERY
  printf("LeafEvent\n");
  DumpObjects(rez);
#endif

  CleanUpTree(p);

  return (rez);
}

CoverPacket *Object2Cover(ObjectPacket *oprnd)
{ int size, *ject, *segm, *data;
  int val, mdp;
  int i, s, n;
  CoverPacket *rez;

  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  data = (int *) malloc(sizeof(int)*(2*size+1));

  s = 0;
  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      { Link *lk;
        lk = LinkIndex[ject[i]];
        data[s++] = lk->emin;
        data[s++] = - (lk->emax + 1);
      }
    else
      { Line *ln;
        ln = LineIndex[ject[i]];
        n = segm[i];
        data[s++] = ln->beg[n];
        n += 2;
        while (ln->beg[n] < 0)
          n += 2;
        data[s++] = - (ln->beg[n] + 1);
      }
  data[s] = 0;

  qsort(data,2*size,sizeof(int),EXPAND_COMPARE);

  rez = (CoverPacket *) malloc(sizeof(CoverPacket));

  val = mdp = 0;
  for (i = 0; i < 2*size; i++)
    { if (data[i] < 0)
        val -= 1;
      else
        val += 1;
      if (val > mdp) mdp = val;
    }

  rez->numpts  = 2*size;
  rez->maxdep  = mdp;
  rez->initpt  = 0;
  rez->datapts = data;

#ifdef DEBUG_QUERY
  printf("Object2Cover\n");
  DumpCoverage(rez);
#endif

  return (rez);
}

CoverPacket *CopyOfCover(CoverPacket *oprnd)
{ int size, i, *data;
  CoverPacket *rez;

  size = oprnd->numpts;

  data = (int *) malloc(sizeof(int)*(size+1));

  for (i = 0; i < size; i++)
    data[i] = oprnd->datapts[i];
  data[size] = 0;

  rez = (CoverPacket *) malloc(sizeof(CoverPacket));

  rez->numpts  = size;
  rez->maxdep  = oprnd->maxdep;
  rez->initpt  = oprnd->initpt;
  rez->datapts = data;

#ifdef DEBUG_QUERY
  printf("CopyOfCover\n");
  DumpCoverage(rez);
#endif

  return (rez);
}

  /* Query Evaluation */

typedef struct {
    int64   numpts;
    int64   sumofpts;
    int64  *datapts;
    float stddev;
    int64   bordermax;
} HistoPacket;

static int64 *Histogram;

static int HISTO_COMPARE(const void *l, const void *r)
{ int64 x, y;
  x = *((int64 *) l);
  y = *((int64 *) r);
  return (x - y);
}

static HistoPacket *MakeHistogram(int64 len)
{ HistoPacket *hist;

  if (len == 0)
    hist = NULL;
  else
    { int i; 
      int64 sum, *dps;
      double mean, dev;

      qsort(Histogram,len,sizeof(int64),HISTO_COMPARE);

      hist = (HistoPacket *) malloc(sizeof(HistoPacket));
      dps  = (int64 *) malloc(sizeof(int64)*len);
      sum  = 0;
      for (i = 0; i < len; i++)
        sum += (dps[i] = Histogram[i]);
      mean = (1.*sum) / len;
      dev  = 0.;
      for (i = 0; i < len; i++)
        dev += (dps[i] - mean) * (dps[i] - mean);
      hist->numpts   = len;
      hist->sumofpts = sum;
      hist->datapts  = dps;
      hist->stddev   = sqrt(dev/len);
    }

  return (hist);
}

HistoPacket *Object2Hist(ObjectPacket *oprnd)
{ int size, *ject, *segm;
  int i, n, m, len;

  size = oprnd->numobj;
  ject = oprnd->objects;
  segm = oprnd->segment;

  len = 0;
  for (i = 0; i < size; i++)
    if (segm[i] < 0)
      { Link *lk;
        lk = LinkIndex[ject[i]];
        Histogram[len++] = lk->emax - lk->emin;
      }
    else
      { Line *ln;
        ln = LineIndex[ject[i]];
        n  = segm[i];
        m  = n+2;
        while (ln->beg[m] < 0)
          m += 2;
        Histogram[len++] = ln->beg[m] - ln->beg[n];
      }

  return (MakeHistogram(len));
}

HistoPacket *Cover2Hist(CoverPacket *cover)
{ int size, *data; 
  int i, val, pos, start, len;

  size = cover->numpts;
  data = cover->datapts;
  val  = cover->initpt;
  start = -1;
  len = 0;
  for (i = 0; i < size; i++)
    { pos = data[i];
      if (pos < 0)
        { val -= 1;
          pos = -(pos+1);
          if (val == 0 && start >= 0)
            Histogram[len++] = pos - start;
        }
      else
        { if (val == 0)
            start = pos;
          val += 1;
        }
    }

  return (MakeHistogram(len));
} 

static CoverPacket *InterpTree(PTree *p)
{ CoverPacket *cov;
  ObjectPacket *obj;

  if (p->lft != NULL)
    switch (p->lab)
    { case '&':
        return (Intersection(InterpTree(p->lft),InterpTree(p->rgt)));
      case '|':
        return (Union(InterpTree(p->lft),InterpTree(p->rgt)));
      case '!':
        return (Negate(InterpTree(p->lft)));
      case '>':
        return (GreaterThan(InterpTree(p->lft),p->rgt->lab));
      case '<':
        return (LessThan(InterpTree(p->lft),p->rgt->lab));
      case '_':
        return (Divide(InterpTree(p->lft),p->rgt->lab));
      case '+':
        return (Expander(InterpTree(p->lft),p->rgt->lab));
      case '~':
        return (Contractor(InterpTree(p->lft),p->rgt->lab));
      case '@':
        obj = LeafEvent(p->lft);
        cov = Object2Cover(obj);
        FreeObjects(obj);
        return (cov);
      case '0':
        return (CovIncludes(InterpTree(p->lft),
                            Negate(Negate(InterpTree(p->rgt)))));
      case '1':
        return (CovInter(InterpTree(p->lft),
                         Negate(Negate(InterpTree(p->rgt)))));
      case '2':
        return (CovContains(InterpTree(p->lft),
                            Negate(Negate(InterpTree(p->rgt)))));
      case '$':
        { CoverPacket *cov;
          ObjectPacket *obj;
          cov = QTable[p->lft->lab].refcover;
          obj = QTable[p->lft->lab].refobject;
          if (cov != NULL)
            return (CopyOfCover(cov));
          else
            return (Object2Cover(obj));
        }
    }
  else
    { obj = LeafEvent(p);
      cov = Object2Cover(obj);
      FreeObjects(obj);
      return (cov);
    }
}

/* Process Query: Top level */

char *ProcessQuery(char *query, ObjectPacket **objR, CoverPacket **covrR)
{ PTree *p;
  char  *estring = NULL;
  static int Firstime = 1;

  Scan = query;
  ParseError = 0;
  SkipBlanks();
  p = GetOrClause(0);
  if (!ParseError)
    { if (*Scan == ')')
        { ErrorString = "Missing matching left paren";
          ParseError = 1;
        }
      else if (*Scan != '\0')
        { ErrorString = "Expecting binary operator";
          ParseError = 1;
        }
    }

  *objR = NULL;
  *covrR = NULL;
  if (ParseError)
    { int len, eln;

      len = Scan-query;
      eln = strlen(ErrorString);
      if (len >= eln+1)
        { estring = malloc(len);
          sprintf(estring,"%*s%s ^",len-(eln+1),"",ErrorString);
        }
      else
        { estring = malloc(len+eln+3);
          sprintf(estring,"%*s^ %s",len,"",ErrorString);
        }
    }
  else if (NumLines > 0)
    { 
#ifdef DEBUG_PARSE
      ShowTree(p);
      fflush(stdout);
#endif

      if (Firstime)
        { Firstime = 0;
          Histogram = (int64 *) malloc(sizeof(int64)*(NumEvents+NumLinks));
        }

      if (p->lab == '@')
        *objR = LeafEvent(p->lft);
      else
        *covrR = InterpTree(p);
    }

#ifdef DEBUG_QUERY
  if (*objR != NULL)
    DumpObjects(*objR);
  if (*covrR != NULL)
    DumpCoverage(*covrR);
#endif

  if (p != NULL)
    FreeTree(p);
#ifdef DEBUG_PARSE
  printf("Memleft = %d\n",Pnodes);
  fflush(stdout);
#endif

  if (ParseError)
    return (estring);
  else
    return (NULL);
}

static int RoundDown(int val, int modulus)
{ int xvl;

  xvl = val - (val % modulus);
  if (val < 0 && val % modulus != 0)
    xvl -= modulus;
  return (xvl);
}

static int RoundUp(int val, int modulus)
{ if (val % modulus != 0)
    val = RoundDown(val,modulus) + modulus;
  return (val);
}

void FreeHistogram(HistoPacket *hist)
{ free(hist->datapts);
  free(hist);
}

void GetHistoDims(HistoPacket *hist, int *min, int *max)
{ *min = hist->datapts[0];
  *max = hist->datapts[hist->numpts-1];
}

#define HWIDE    6
#define HWMAX  100
#define HHIGH   20
#define HBORDER 10
#define HTEXTB  18

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
                   int64 xmin, int64 xmax, int *bsize)
{ int xl, xh, yl, yh;
  int64 lst, bot;

  static char  label[100];
  static char *letter[] = { "", "K", "M", "G", "T" };

  mt_get_extent(canvas,&xl,&xh,&yl,&yh);
  xh -= (xl-1);
  yh -= (yl-1);
  
  { int xzeros, xlead, xorder;
    int yzeros, ylead, yorder;
    int64 bmin, bmax, hmax;
    int64 lborder, hwide;

    lborder = hist->bordermax;

    xlead = TickRounder((COORDSCALE*1.*(xmax-xmin)) / ((xh-(2*HBORDER+lborder))/HWIDE),
                        &xorder,&xzeros);
    AdjustTicks(bsize,&xlead,&xorder,&xzeros);

    while (1)
      { int adjust;
      
        bmin  = RoundDown(COORDSCALE*xmin,xlead);
        bmax  = RoundDown(COORDSCALE*xmax,xlead) + xlead;
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

    ylead = TickRounder( (1.*hmax) / ((yh-(HBORDER+HTEXTB))/HHIGH),
                         &yorder, &yzeros);
    hmax = RoundUp(hmax,ylead);

    { int w, h ,b;

      sprintf(label,F_S64,hmax);
      mt_string_size(label,&w,&h,&b);
      if (bmin <= hist->datapts[0] && bmax > hist->datapts[hist->numpts-1])
        lst = (w+4) + (xh - (hwide*(bmax-bmin)/xlead + w + 4))/2;
      else
        lst = lborder + HBORDER;
      bot = yh - HTEXTB;
    }

    { int i, k, b, y, lft;
      double yfact;

      mt_set_color(mt_get_color(255,255,0));
      i = 0;
      yfact = (1.*yh - (HBORDER+HTEXTB))/hmax; 
      while (i < hist->numpts && hist->datapts[i] < bmin)
        i += 1;
      lft = lst;
      for (k = bmin+xlead; k <= bmax; k += xlead)
        { b = 0;
          while (i < hist->numpts && hist->datapts[i] < k)
            { b += 1; i += 1; }
          y = yfact * b;
          if (y == 0 && b != 0) y = 1;
          mt_draw_rect(canvas,lft,bot-y,hwide-1,y);
          lft += hwide;
        }
    }

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

      mt_draw_line(canvas,lst-1,bot,lst-1,HBORDER,1,0);
      yfact = (1.*yh - (HBORDER+HTEXTB))/hmax; 
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

static int printwidth(int value)
{ int width;

  if (value == 0)
    width = 1;
  else if (value < 0)
    width = log10(-1.*value) + 2;
  else
    width = log10(1.*value) + 1;
  return (width);
}

int DrawHistoInfoLabel(MT_OBJECT *bar, HistoPacket *hist, char *outlabel)
{ char label[100]; 
  int border, max, width, alt;
  int w, h, b;
  char *outlabelptr=outlabel;

  mt_set_color(mt_black());
  border = 5;

  width = printwidth(hist->sumofpts);

  sprintf(label,"Min = %*" F_S64P,width,COORDSCALE*(hist->datapts[0]));
  outlabelptr+=sprintf(outlabelptr,"%*" F_S64P "\t",width,COORDSCALE*(hist->datapts[0]));
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  sprintf(label,"Max = %*" F_S64P,width,COORDSCALE*(hist->datapts[hist->numpts-1]));
  outlabelptr+= sprintf(outlabelptr,"%*" F_S64P "\t",width,COORDSCALE*(hist->datapts[hist->numpts-1]));
  if (bar != NULL) mt_draw_text(bar,border,23,label);
  sprintf(label,"Sum = %*" F_S64P,width,COORDSCALE*(hist->sumofpts));
  outlabelptr+=sprintf(outlabelptr,"%*" F_S64P "\t",width,COORDSCALE*(hist->sumofpts));
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

  sprintf(label,"Mean     = %*.2f",width+3,COORDSCALE*(1.*hist->sumofpts)/hist->numpts);
  outlabelptr+=sprintf(outlabelptr,"%*.2f\t",width+3,COORDSCALE*(1.*hist->sumofpts)/hist->numpts);
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  sprintf(label,"Median   = %*" F_S64P,width,COORDSCALE*hist->datapts[hist->numpts/2]);
  outlabelptr+=sprintf(outlabelptr,"%*" F_S64P "\t",width,COORDSCALE*hist->datapts[hist->numpts/2]);
  if (bar != NULL) mt_draw_text(bar,border,34,label);
  sprintf(label,"Std Dev. = %*.2f",width+3,COORDSCALE*hist->stddev);
  outlabelptr+=sprintf(outlabelptr,"%*.2f\t",width+3,COORDSCALE*hist->stddev);
  if (bar != NULL) mt_draw_text(bar,border,23,label);

  mt_string_size(label,&w,&h,&b);
  border += w+10;

  sprintf(label,"%% of total = %.2f",
          (100.*hist->sumofpts) / (MaxPnt - MinPnt));
  outlabelptr+=sprintf(outlabelptr,"%.2f\t",
          (100.*hist->sumofpts) / (MaxPnt - MinPnt));
  if (bar != NULL) mt_draw_text(bar,border,12,label);
  mt_string_size(label,&w,&h,&b);
  max = w;
  sprintf(label,"# of pts.  = " F_S64,hist->numpts);
  outlabelptr+=sprintf(outlabelptr,F_S64,hist->numpts);
  if (bar != NULL) mt_draw_text(bar,border,23,label);
  mt_string_size(label,&w,&h,&b);
  if (w > max) max = w;

  border += max+5;

  return (border);
}

void FreeCoverage(CoverPacket *cover)
{ free(cover->datapts);
  free(cover);
}

void DrawCoverage(MT_OBJECT *canvas, CoverPacket *cover, int xmin, int xmax)
{ int xl, xh, yl, yh;

  mt_get_extent(canvas,&xl,&xh,&yl,&yh);
  xh -= xl;
  yh -= yl;

#define CHIGH   20
#define CBORDER 10
  
  { int yzeros, ylead, yorder;
    int hmax;

    hmax = cover->maxdep;

    ylead = TickRounder( (1.*hmax) / ((yh-2*CBORDER)/HHIGH),
                         &yorder, &yzeros);
    hmax = RoundUp(hmax,ylead) + 1;

    { int i, d, e;
      int s, x, v;
      double xfact, yfact;

#define XCMAP(x) ((int) (((x) - xmin)*xfact))
#define YCMAP(y) ((int) ((hmax - (y))*yfact))

      xfact = (1.*xh)/(xmax-xmin);
      yfact = (1.*yh)/(hmax+1);

      mt_set_color(mt_get_color(0,125,75));
      mt_draw_line(canvas,0,YCMAP(0),xh+1,YCMAP(0),3,0);
      mt_set_color(mt_get_color(255,255,0));

      i = 0;
      e = d = cover->initpt;
      s = xmin;
      while (i < cover->numpts)
        { v = cover->datapts[i++];
          if (v < 0)
            { x = -(v+1); d -= 1; }
          else
            { x = v; d += 1; }
          while (cover->datapts[i] == v)
            { i += 1;
              if (v < 0)
                d -= 1;
              else
                d += 1;
            }
          if (x >= xmin)
            { if (x > xmax)
                { mt_draw_line(canvas,XCMAP(s),YCMAP(e),
                                      XCMAP(xmax)+1,YCMAP(e),1,0);
                  break;
                }
              else
                mt_draw_line(canvas,XCMAP(s),YCMAP(e),XCMAP(x)+1,YCMAP(e),1,0);
              mt_draw_line(canvas,XCMAP(x),YCMAP(e),XCMAP(x),YCMAP(d),1,0);
              s = x;
            }
          e = d;
        }
      if (x < xmax)
        mt_draw_line(canvas,XCMAP(s),YCMAP(e),XCMAP(xmax)+1,YCMAP(e),1,0);
    }
  }

}
