
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
 * $Id: AS_CGB_io.h,v 1.4 2005-03-22 19:48:29 jason_miller Exp $
 *
 * Module: AS_CGB_io.h
 * Description: Header file for the code that reads and writes the 
 * check point data.
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_IO_INCLUDE
#define AS_CGB_IO_INCLUDE

void output_the_chunks
(/* Input Only*/
 MesgWriter WriteMesg_AS,
 const Tfragment frags[],
 const Tedge     edges[],
 const VA_TYPE(char) fragsrc[],
 const TChunkFrag    chunkfrags[],
 const TChunkMesg    thechunks[],
 const VA_TYPE(char) chunkseqs[],
 const VA_TYPE(char) chunkquas[],
 const VA_TYPE(char) chunksrc[],
 const int analysis_flag,
 const int output_fom_messages,
 const float global_fragment_arrival_rate,
 const int fragment_count_target,
 const char * const Graph_Store_File_Prefix
);

void convert_the_chunks_to_IUM
(/* Input Only*/
 const Tfragment        * frags,
 const Tedge            * edges,
 const VA_TYPE(char)    * fragsrc,
 const TChunkFrag       * chunkfrags,
 const TChunkMesg       * thechunks,
 const VA_TYPE(char)    * chunkseqs,
 const VA_TYPE(char)    * chunkquas,
 const VA_TYPE(char)    * chunksrc,
 const int                analysis_level,
 const float              global_fragment_arrival_rate,
 /* Output Only */
 VA_TYPE(IntMultiPos)   * the_imps,
 VA_TYPE(IntUnitigMesg) * the_iums,
 VA_TYPE(char)   * the_imp_source,
 VA_TYPE(char)   * the_ium_source
);

void output_the_IUM_to_file
(/* Input Only*/
 MesgWriter                  WriteMesg_AS,
 const VA_TYPE(char)    *    fragsrc,
 const VA_TYPE(char)    *    chunksrc,
 VA_TYPE(IntMultiPos)   *    the_imps,
 VA_TYPE(IntUnitigMesg) *    the_iums,
 const int                   fragment_count_target,
 const char * const          Graph_Store_File_Prefix
);


#endif // AS_CGB_IO_INCLUDE
