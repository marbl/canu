
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
/*********************************************************************
   Module:  AS_OVL
   Description:  Functions to access binary overlap store
   Assumptions:  
 *********************************************************************/


/* RCS info
 * $Id: OlapStoreOVL.c,v 1.1.1.1 2004-04-14 13:52:43 catmandew Exp $
 * $Revision: 1.1.1.1 $
*/


#include  "OlapStoreOVL.h"
#include  "AS_OVL_delcher.h"

float32  Expand_Quality
    (int q)

//  Convert  q  from condensed form to floating point equivalent
//  and return the result

  {
   return  q / 1000.0;
  }



void  Free_OVL_Store
    (OVL_Store_t * store)

//  Free memory associated with  store

  {
   free (store);

   return;
  }



void  Free_OVL_Stream
    (OVL_Stream_t * stream)

//  Free memory associated with  stream

  {
   if  (stream -> fp != NULL)
       fclose (stream -> fp);
   if  (stream -> offset_buff != NULL)
       free (stream -> offset_buff);
   free (stream);

   return;
  }



int16  Get_Int_Quality
    (int q)

//  Convert  q  from form in overlap data structure to integer
//  and return the result

  {
   return  q;
  }



void  Init_OVL_Stream
    (OVL_Stream_t * stream, uint32 first, uint32 last, OVL_Store_t * store)

//  Initialize  stream  to read overlaps for fragments  first .. last
//  from overlap store  store .

  {
   char  filename [MAX_FILENAME_LEN];
   size_t  file_position;
   uint32  i, header [3];
   int  len, new_file_index;

   if  (store == NULL || store -> name == NULL || store -> offset_fp == NULL)
       {
        fprintf (stderr,
                 "ERROR:  In Init_OVL_Stream overlap store has not been opened\n");
        exit (EXIT_FAILURE);
       }
   if  (first < 1 || first > store -> max_frag)
       {
        fprintf (stderr,
                 "ERROR:  Init_OVL_Stream first = %u not in range 1 .. %u\n",
                 first, store -> max_frag);
        exit (EXIT_FAILURE);
       }
   if  (last < 1 || last > store -> max_frag)
       {
        fprintf (stderr,
                 "ERROR:  Init_OVL_Stream last = %u not in range 1 .. %u\n",
                 last, store -> max_frag);
        exit (EXIT_FAILURE);
       }
   if  (last < first)
       {
        fprintf (stderr,
                 "ERROR:  Init_OVL_Stream last = %u is less than first = %u\n",
                 last, first);
        exit (EXIT_FAILURE);
       }

   rewind (store -> offset_fp);
   Safe_fread (header, sizeof (uint32), 3, store -> offset_fp);
   assert (store -> max_frag = header [0]);
   assert (header [1] == sizeof (Short_Olap_Data_t));
   assert (store -> frags_per_file = header [2]);

   stream -> start_id = first;
   stream -> stop_id = last;
   len = 2 + last - first;

   stream -> offset_buff = (uint32 *) Safe_realloc (stream -> offset_buff,
                                                    len * sizeof (uint32));
   CDS_FSEEK (store -> offset_fp,
              (off_t) (first * sizeof (uint32)), SEEK_CUR);
   Safe_fread (stream -> offset_buff, sizeof (uint32), len, store -> offset_fp);

   for  (i = first;  i <= last;  i ++)
     {
      if  (i % store -> frags_per_file == 0
             || stream -> offset_buff [i - first]
                   != stream -> offset_buff [1 + i - first])
          break;
     }

   stream -> curr_id = i;
   new_file_index = (int) ceil ((double) i / store -> frags_per_file);
   if  (stream -> curr_file_index != new_file_index
            || stream -> store != store)
       {
        stream -> curr_file_index = new_file_index;
        if  (stream -> fp != NULL)
            fclose (stream-> fp);
        sprintf (filename, "%s/data%02d.olap",
                 store -> name, stream -> curr_file_index);
        stream -> fp = Local_File_Open (filename, "rb");
       }
   file_position = stream -> offset_buff [i - first] * sizeof (Short_Olap_Data_t);
   CDS_FSEEK (stream -> fp, (off_t) file_position, SEEK_SET);
   stream -> curr_offset = stream -> offset_buff [i - first];
   stream -> store = store;

   return;
  }



uint32  Last_Frag_In_OVL_Store
    (OVL_Store_t * store)

//  Returns the iid of the highest numbered fragment in overlap
//  store  store .  The store must have been opened previously.

  {
   return  store -> max_frag;
  }



uint32 *  Load_Frag_Offsets
    (OVL_Store_t * store)

//  Allocate memory to hold the offset values in ovl store  store ,
//  load them, and return a pointer to them.  Store must have
//  been opened already.

  {
   uint32  * offset;
   uint32  num_frags;

   num_frags = 2 + store -> max_frag;
     // empty spot for non-existent frag 0 and one spot for fragment
     // off the end
   offset = (uint32 *) Safe_malloc (num_frags * sizeof (uint32));
   CDS_FSEEK (store -> offset_fp,
              (off_t) (3 * sizeof (uint32)), SEEK_SET);
     // skip over header
   Safe_fread (offset, sizeof (uint32), num_frags, store -> offset_fp);

   return  offset;
  }



FILE *  Local_File_Open
    (const char * filename, const char * mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (filename, mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR %d:  Could not open file  %s \n",
                 errno, filename);
        perror (strerror (errno));
        exit (EXIT_FAILURE);
       }

   return  fp;
  }



OVL_Store_t *  New_OVL_Store
    (void)

//  Allocate and return a pointer for a new overlap store.

  {
   OVL_Store_t  * store;

   store = (OVL_Store_t *) Safe_malloc (sizeof (OVL_Store_t));
   store -> name = NULL;
   store -> offset_fp = NULL;

   return  store;
  }



OVL_Stream_t *  New_OVL_Stream
    (void)

//  Allocate and return a pointer for a new overlap stream.

  {
   OVL_Stream_t  * stream;

   stream = (OVL_Stream_t *) Safe_malloc (sizeof (OVL_Stream_t));
   stream -> store = NULL;
   stream -> offset_buff = NULL; // For realloc of PC/Linux/gcc !
   stream -> start_id = stream -> stop_id = stream -> curr_id = 0;
   stream -> fp = NULL;
   stream -> curr_file_index = -1;   // To avoid random matches

   return  stream;
  }



int  Next_From_OVL_Stream
    (Long_Olap_Data_t * olap, OVL_Stream_t * stream)

//  Store in  (* olap) the next overlap from  stream .
//  Return  TRUE  if successful;  FALSE , otherwise.

  {
   OVL_Store_t  * store;
   Short_Olap_Data_t  short_olap;
   uint32  i;

   if  (stream -> curr_id > stream -> stop_id)
        return  FALSE;

   store = stream -> store;

   if  (fread (& short_olap, sizeof (Short_Olap_Data_t), 1, stream -> fp) == 1)
       {
        olap -> a_iid = stream -> curr_id;
        olap -> b_iid = short_olap . b_iid;
        olap -> flipped = short_olap . flipped;
        olap -> a_hang = short_olap . a_hang;
        olap -> b_hang = short_olap . b_hang;
        olap -> orig_erate = short_olap . orig_erate;
        olap -> corr_erate = short_olap . corr_erate;
        stream -> curr_offset ++;
        for  (i = stream -> curr_id;  i <= stream -> stop_id;  i ++)
          {
           if  (i % store -> frags_per_file == 0
                  || stream -> curr_offset
                        < stream -> offset_buff [1 + i - stream -> start_id])
               break;
          }
        stream -> curr_id = i;
        return  TRUE;
       }

   if  (stream -> curr_id % store -> frags_per_file == 0)
       {
        char  filename [MAX_FILENAME_LEN];
        size_t  file_position;

        if  (stream -> curr_id == stream -> stop_id)
            return  FALSE;

        stream -> curr_id ++;
        stream -> curr_file_index ++;
        fclose (stream -> fp);
        
        sprintf (filename, "%s/data%02d.olap",
                 store -> name, stream -> curr_file_index);
        stream -> fp = Local_File_Open (filename, "rb");

        for  (i = stream -> curr_id;  i <= stream -> stop_id;  i ++)
          {
           if  (i % store -> frags_per_file == 0
                  || stream -> offset_buff [i - stream -> start_id]
                        != stream -> offset_buff [1 + i - stream -> start_id])
               break;
          }

        stream -> curr_id = i;
        file_position = stream -> offset_buff [i - stream -> start_id]
                          * sizeof (Short_Olap_Data_t);
        CDS_FSEEK (stream -> fp, (off_t) file_position, SEEK_SET);
        stream -> curr_offset = stream -> offset_buff [i - stream -> start_id];

        return  Next_From_OVL_Stream (olap, stream);
       }

   fprintf (stderr, "ERROR reading OVL store \"%s\"\n", store -> name);
   fprintf (stderr, "  curr_id = %u  file_index = %d  curr_offset = %u\n",
            stream -> curr_id, stream -> curr_file_index, stream -> curr_offset);
   exit (EXIT_FAILURE);
   
   return  FALSE;
  }



int  Open_OVL_Store
    (OVL_Store_t * store, const char * path)

//  Open the overlap store at  path  for read-only access and
//  make  store  refer to it.   store must have been previously
//  allocated by  New_OVL_Store .  Return  0  if successful; -1
//  otherwise.

  {
   char  filename [MAX_FILENAME_LEN];
   OVL_Store_ID_t  id;
   uint32  header [3];

   assert (strlen (path) + strlen (OFFSET_FILENAME) < MAX_FILENAME_LEN);
   strcpy (filename, path);
   strcat (filename, OFFSET_FILENAME);

   store -> name = strdup (path);
   store -> offset_fp = Local_File_Open (filename, "rb");
   Safe_fread (header, sizeof (uint32), 3, store -> offset_fp);
   store -> max_frag = header [0];
   id . header = header [1];
   assert (id . tag . version == OVL_STORE_VERSION);
   assert (id . tag . record_size == sizeof (Short_Olap_Data_t));
   store -> frags_per_file = header [2];

   return  0;
  }



int  Shrink_Quality
    (float32 q)

//  Convert  q  to a discrete, integral form that uses less space

  {
   double  x = (1000.0 * q + 0.5);

   if  (x > MAX_ERATE)
       return  MAX_ERATE;
     else
       return  (int) x;
  }
