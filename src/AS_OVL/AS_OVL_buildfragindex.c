
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
/* ************************************************
   Module:
   Description:

   Programmer:  A. Delcher
   Written:  27 October 1998
   File:  ~adelcher/Assembly/Overlaps/buildfragindex.c

   Read fragment data (from the argv[1].urc) and build a
   short fragment index from it, argv[1].sfr and a sequence
   index argv[1].seq.
   Assumptions:

*************************************************/

/* RCS info
 * $Id: AS_OVL_buildfragindex.c,v 1.1.1.1 2004-04-14 13:52:27 catmandew Exp $
 */

#include "AS_global.h"
#include  "AS_OVL_delcher.h"
#include  "AS_PER_ReadStruct.h"
#include  "AS_PER_fragStore.h"

#define  ALPHABET_SIZE       4
#define  DONT_KNOW_CHAR      'n'
#define  MAX_COMMENT_LEN     4096
#define  MAX_FRAG_LEN        2048
#define  MAX_LINE_LEN        2048

#define  DEBUG               1
#define  LONGEST_EXTENSION  20    // Number of characters added to index name
                                  //   to designate a separate file.
#define  SHOW_STATS          0


typedef  uint64  Distance_t;

typedef  enum Frag_Type  {READ, GUIDE, REREAD}
  Frag_Type_t;

typedef  struct Fixed_Info
  {
   int  Accession;
   Frag_Type_t  Type;
   time_t  Entry_Time;  /* int32 seconds from sometime in 1970 */
   int  Frag_Len;
   uint64  Sequence_Start;
   uint64  Quality_Start;
   int  Clear_Begin, Clear_End;
   uint64  Screen_Start;
   int  Mate;
   Distance_t  Distance;
   int  ReRead;
   uint64  Guide_Start;
  }  Fixed_Info_t;





int  Get_Field
    (char *, char *, void *);
int  Get_String_Field
    (char *, char *, int, int);
int  Get_Two_Fields
    (char *, char *, void *, void *);
int  Read_Frag_Mesg
    (Fixed_Info_t *, char *, char *, char *, int *);
void  Skip_Rest_Of_Mesg
    (char []);
void  Echo_Rest_Of_Mesg
    (FILE *fp, char Line []);



#ifdef SFR
void  outputFragmentDigest(FILE *, int x, Fixed_Info_t *Frag, char *Source);
#endif

FILE  * Fixed_fp, 
      * Seq_fp, 
      * Quality_fp, 
      *Frg_fp,
      *Sfr_fp;


int  main  (int argc, char * argv [])

  {
   Fixed_Info_t  Frag;
   char  Quality [MAX_FRAG_LEN], 
         Sequence [MAX_FRAG_LEN],
         Source   [MAX_FRAG_LEN];
   char  * File_Name = NULL;
   int64  Curr_Frag_Num;
   int  i, Skip;
   int start, end;
   char *interval, *repeats;
   ReadStructp myReadStruct;
   FragStoreHandle myFragStoreHandle;

   if  (argc < 2)
       {
        fprintf (stderr, "USAGE:  %s <frag-index-name>  < input \n", argv [0]);
        exit (EXIT_FAILURE);
       }

   File_Name = (char *) Safe_malloc (strlen (argv [1]) + LONGEST_EXTENSION);

   myFragStoreHandle = createFragStore(argv[1], "george", 1);
   myReadStruct = new_ReadStruct();  

   strcpy (File_Name, argv [1]);              // one string per line
   strcat (File_Name, ".urc");                 // clear-range only
   Frg_fp = File_Open (File_Name, "r");
   if(Frg_fp == NULL){
     fprintf(stderr,"%s: couldn't open input file %s\n",
	     argv[1], File_Name);
     exit(1);
   }

   strcpy (File_Name, argv [1]);              // one string per line
   strcat (File_Name, ".sfr");                 // clear-range only
   Sfr_fp = File_Open (File_Name, "w");
   if(Sfr_fp == NULL){
     char buffer[2048];
     sprintf(buffer,"%s: couldn't open output file %s",
	     argv[0], File_Name);
     fprintf(stderr,"%s\n",buffer);
     perror(buffer);
     exit(1);
   }

   strcpy (File_Name, argv [1]);              // one string per line
   strcat (File_Name, ".seq");                 // clear-range only
   Seq_fp = File_Open (File_Name, "w");
   if(Seq_fp == NULL){
     char buffer[2048];
     sprintf(buffer,"%s: couldn't open output file %s",
	     argv[0], File_Name);
     fprintf(stderr,"%s\n",buffer);
     perror(buffer);
     exit(1);
   }


   Curr_Frag_Num = 1;
   while  (Read_Frag_Mesg (& Frag, Sequence, Quality, Source, & Skip))
     {
      if  (Skip)
          continue;
      for  (i = Frag . Clear_Begin;  i < Frag . Clear_End;  i ++)
        fputc (Sequence [i], Seq_fp);
      fputc ('\n', Seq_fp);

#ifdef SFR
      outputFragmentDigest(Sfr_fp,Curr_Frag_Num, &Frag, Source);
#endif
      setAccID_ReadStruct(myReadStruct, Frag.Accession);
      setReadIndex_ReadStruct(myReadStruct, Curr_Frag_Num);
      setReadType_ReadStruct(myReadStruct, ReadType_read);
      setMateIndex_ReadStruct(myReadStruct, Frag.Mate);
      setMateDistanceIndex_ReadStruct(myReadStruct, Frag.Distance);
      setSequence_ReadStruct(myReadStruct, Sequence, Sequence);
      setSource_ReadStruct(myReadStruct, Source);
      setEntryTime_ReadStruct(myReadStruct, Frag.Entry_Time);
      setClearRegion_ReadStruct(myReadStruct, Frag.Clear_Begin, Frag.Clear_End);

      appendFragStore(myFragStoreHandle, myReadStruct);
#ifdef SFR
      outputFragmentDigest(Sfr_fp,Curr_Frag_Num++, &Frag, Source);
#endif


      Curr_Frag_Num++;

     }

   closeFragStore(myFragStoreHandle);
   fclose (Seq_fp);
   fclose (Frg_fp);


   return  0;
  }



int  Read_Frag_Mesg
    (Fixed_Info_t * Frag, char * Sequence, char * Quality, char * Source, int * Skip)

//  Read next fragment message from stardard input and store
//  its fixed info in  Frag , its sequence in  Sequence , and
//  its quality value characters in  Quality .  Set  (* Skip)
//  true if this message should be ignored (currently if a distance
//  message).

  {
   char  Ch, Line [MAX_LINE_LEN], String [MAX_COMMENT_LEN];
   int  Len, Num;

   while  (fgets (Line, MAX_LINE_LEN, Frg_fp) != NULL)
     {
      Len = strlen (Line);
      assert (Line [Len - 1] == '\n');
      if  (Len == 1)
          continue;                       // Skip blank lines
      if  (Line [0] != '{')
          {
           Line [Len - 1] = '\0';
           fprintf (stderr, "ERROR:  Unexpected Line:\n  \"%s\"\n",
                   Line);
           exit (EXIT_FAILURE);
          }
      if  (strncmp (Line + 1, "DST", 3) == 0)
          {
           //  Distance record--skip for now
           Echo_Rest_Of_Mesg (Sfr_fp,Line);
           (* Skip) = TRUE;
           return  TRUE;
          }

      if  (strncmp (Line + 1, "SFG", 3) != 0)
          {
           Line [Len - 1] = '\0';
           fprintf (stderr,
                "ERROR:  Unexpected message type, skipping message:\n  \"%s\"\n",
                Line);
           Skip_Rest_Of_Mesg (Line);
           (* Skip) = TRUE;
           return  TRUE;
          }

      if  (! Get_Field ("act", "%c", & Ch))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }
      switch  (Ch)
        {
         case  'A' :
           break;
         case  'D' :
           // Disallow for now
         default :
           fprintf (stderr, "ERROR:  Unexpected  act  value \"%c\"\n", Ch);
           Line [0] = '\0';
           Skip_Rest_Of_Mesg (Line);
           (* Skip) = TRUE;
           return  TRUE;
        }

      if  (! Get_Field ("acc", "%d", & Num))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }
      Frag -> Accession = Num;

      if  (! Get_Field ("typ", "%c", & Ch))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }
      switch  (Ch)
        {
         case  'R' :
         case  'G' :
         case  'E' :
           Frag -> Type = Ch;
           break;
         default :
           fprintf (stderr, "ERROR:  Unexpected  typ  value \"%c\"\n", Ch);
           Line [0] = '\0';
           Skip_Rest_Of_Mesg (Line);
        }

      if  (! Get_String_Field ("src", Source, MAX_COMMENT_LEN, TRUE))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }
      // Ignore source comment string for now.

      if  (! Get_Field ("etm", "%d", & (Frag -> Entry_Time)))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }

      if  (! Get_String_Field ("seq", Sequence, MAX_FRAG_LEN, FALSE))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }

      Frag->Frag_Len = strlen(Sequence);

      if  (! Get_String_Field ("qlt", Quality, MAX_FRAG_LEN, FALSE))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }

      assert(Frag->Frag_Len == strlen(Quality));

      if  (! Get_Two_Fields ("clr", "%d,%d", & (Frag -> Clear_Begin),
                             & (Frag -> Clear_End)))
          {
           (* Skip) = TRUE;
           return  TRUE;
          }
      /* A HACK!! */
      { 
	char goop[4096];
	if(Get_String_Field ("scn", goop, 4096, FALSE) == 0)
          assert(0);
      }
      
      if  (Frag -> Clear_Begin < 0
              || Frag -> Clear_End <= Frag -> Clear_Begin
              || Frag -> Clear_End > strlen (Sequence))
          {
           fprintf (stderr, "ERROR:  Bad clear range \"%d,%d\"\n",
                    Frag -> Clear_Begin, Frag -> Clear_End);
           Skip_Rest_Of_Mesg (Line);
           (* Skip) = TRUE;
           return  TRUE;
          }

      (Frag -> Distance) = 0L;
      (Frag -> Mate)     = 0L;

      if  (! 
	   Get_Field ("mat", F_UID, & (Frag -> Mate))
	   )
          {
           (* Skip) = TRUE;
           return  TRUE;
          }

      if  (!
	   Get_Field ("dst", F_UID, & (Frag -> Distance))
	   )
          {
           (* Skip) = TRUE;
           return  TRUE;
          }

      Skip_Rest_Of_Mesg (Line);

      (* Skip) = FALSE;
      return  TRUE;
     }

   return  FALSE;
  }



int  Get_Field
    (char * Tag, char * Format, void * Value)

//  Read next line of  Frg_fp  and check that the code at beginning is
//  Tag  followed by a colon.  If so, read the next value into  (* Value)
//  using the format in string  Format .  Return whether successful.

  {
   char  Line [MAX_LINE_LEN];
   int  Len, Tag_Len;

   if  (fgets (Line, MAX_LINE_LEN, Frg_fp) == NULL)
       return  FALSE;
   Len = strlen (Line);
   assert (Line [Len - 1] == '\n');
   Tag_Len = strlen (Tag);
   if  (strncmp (Line, Tag, Tag_Len) != 0 || Len <= Tag_Len
          || Line [Tag_Len] != ':')
       {
        Line [Len - 1] = '\0';
        fprintf (stderr,
              "ERROR:  Expected  %s:  field, skipping message:\n  \"%s\"\n",
              Tag, Line);
        Skip_Rest_Of_Mesg (Line);
        return  FALSE;
       }

   if  (sscanf (Line + Tag_Len + 1, Format, Value) != 1)
       {
        Line [Len - 1] = '\0';
        fprintf (stderr,
              "ERROR:  Failed to read  %s:  field, skipping message:\n  \"%s\"\n",
              Tag, Line);
        Skip_Rest_Of_Mesg (Line);
        return  FALSE;
       }

   return  TRUE;
  }



int  Get_String_Field
    (char * Tag, char * Value, int Max_Len, int PreserveNewLines)

//  Read next line of  Frg_fp  and check that the code at beginning is
//  Tag  followed by a colon.  Then read subsequent lines, and concatenate
//  them onto  Value  until a line beginning with a  .  is read.
//  Fail if string is  Max_Len  or longer.
//  Return whether successful.

  {
   char  Line [MAX_LINE_LEN];
   int  Len, Tag_Len, Value_Len;

   if  (fgets (Line, MAX_LINE_LEN, Frg_fp) == NULL)
       return  FALSE;
   Len = strlen (Line);
   assert (Line [Len - 1] == '\n');
   Tag_Len = strlen (Tag);
   if  (strncmp (Line, Tag, Tag_Len) != 0 || Len <= Tag_Len
          || Line [Tag_Len] != ':')
       {
        Line [Len - 1] = '\0';
        fprintf (stderr,
              "ERROR:  Expected  %s:  field, skipping message:\n  \"%s\"\n",
              Tag, Line);
        Skip_Rest_Of_Mesg (Line);
        return  FALSE;
       }

   Value [0] = '\0';
   Value_Len = 0;
   while  (fgets (Line, MAX_LINE_LEN, Frg_fp) != NULL)
     {
      Len = strlen (Line);
      assert (Line [Len - 1] == '\n');
      if  (strncmp (Line, ".\n", 2) == 0)
          return  TRUE;
      if  (strcmp (Line, "}\n") == 0)
          {
           fprintf (stderr, "ERROR:  Unexpected end of message\n");
           return FALSE;
          }
      if(PreserveNewLines == FALSE){ /* Truncate the newline */
	Line [Len - 1] = '\0';
	assert (Value_Len + Len - 1 < Max_Len);
	strcat (Value, Line);
	Value_Len += Len - 1;
      }else{  /* Don't truncate the newline */
	assert (Value_Len + Len  < Max_Len);
	strcat (Value, Line);
	Value_Len += Len;
      }
     }

   fprintf (stderr, "ERROR:  Unexpected end of file\n");
   exit (EXIT_FAILURE);

   return  FALSE;
  }



int  Get_Two_Fields
    (char * Tag, char * Format, void * Val1, void * Val2)

//  Read next line of  Frg_fp  and check that the code at beginning is
//  Tag  followed by a colon.  If so, read the next values into  (* Val1)
//  and  (* Val2)  using the format in string  Format .
//  Return whether successful.

  {
   char  Line [MAX_LINE_LEN];
   int  Len, Tag_Len;

   if  (fgets (Line, MAX_LINE_LEN, Frg_fp) == NULL)
       return  FALSE;
   Len = strlen (Line);
   assert (Line [Len - 1] == '\n');
   Tag_Len = strlen (Tag);
   if  (strncmp (Line, Tag, Tag_Len) != 0 || Len <= Tag_Len
          || Line [Tag_Len] != ':')
       {
        Line [Len - 1] = '\0';
        fprintf (stderr,
              "ERROR:  Expected  %s:  field, skipping message:\n  \"%s\"\n",
              Tag, Line);
        Skip_Rest_Of_Mesg (Line);
        return  FALSE;
       }

   if  (sscanf (Line + Tag_Len + 1, Format, Val1, Val2) != 2)
       {
        Line [Len - 1] = '\0';
        fprintf (stderr,
              "ERROR:  Failed to read  %s:  field, skipping message:\n  \"%s\"\n",
              Tag, Line);
        Skip_Rest_Of_Mesg (Line);
        return  FALSE;
       }

   return  TRUE;
  }



void  Skip_Rest_Of_Mesg
    (char Line [])

//  Read lines until up through next end-of-message line, i.e.,
//  a line beginning with  "}" .  Don't read anything if  Line
//  is already that line.

  {
   while  (Line [0] != '}')
     {
      if  (fgets (Line, MAX_LINE_LEN, Frg_fp) == NULL)
          return;
     }

   return;
  }

void  Echo_Rest_Of_Mesg
    (FILE *fp, char Line [])

//  Read lines until up through next end-of-message line, i.e.,
//  a line beginning with  "}" .  Don't read anything if  Line
//  is already that line.

  {
    fputs(Line,fp);
   while  (Line [0] != '}')
     {
      if  (fgets (Line, MAX_LINE_LEN, Frg_fp) == NULL)
          return;
      fputs(Line,fp);

     }

   return;
  }


#ifdef SFR
void  outputFragmentDigest(FILE *fp, int id, Fixed_Info_t *Frag, char *Source){

    fprintf(fp, "{SFR\n");
    fprintf(fp, "act:A\n");
    fprintf(fp, "fid:%d\n", id );
    fprintf(fp, "len:%d\n", Frag->Frag_Len);
    fprintf(fp, "src:\n%s.\n", Source);
    fprintf(fp, "mat:" F_UID "\n", Frag->Mate );
    fprintf(fp, "dst:" F_UID "\n", Frag->Distance);
    fprintf(fp, "}\n");
}
#endif

