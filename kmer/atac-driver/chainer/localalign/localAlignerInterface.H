/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
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
**************************************************************************/


/*********************************************************************
 *		  Confidential -- Do Not Distribute                  *
 *	    Copyright (c) 2003 Celera Genomics Corporation.          *
 *			 All Rights Reserved.                        *
 *********************************************************************/

#ifndef GF_SYN_LOCALALIGNERINTERFACE_H
#define GF_SYN_LOCALALIGNERINTERFACE_H

void syntenicSegments
(
 FILE * outfile,
 char const * const Aseq, int const Astart, int const Astop, // substring of Aseq
 char const * const Bseq, int const Bstart, int const Bstop, // substring of Bseq
 float const erate
 );

int iterate_Local_Overlap
(
 int& seg_abpos,
 int& seg_bbpos,
 int& seg_alen,
 int& seg_blen,
 float& seg_error
 );

int Print_Local_Overlap4(
 FILE *outfile,
 char const * selfid,
 char const * parentid,
 char const * a_uid,
 char const * b_uid,
 int const a_scaf_pos, // Lowest position in string "A" for the alignment.
 int const b_scaf_pos, // Lowest position in string "B" for the alignment.
 int const alen,       // Length of the alignment box for string "A".
 int const blen,       // Length of the alignment box for string "B".
 int const sequence_reversed_match, /* bool */
 int const a_lft_seed_len,
 int const b_lft_seed_len,
 int const a_rht_seed_len,
 int const b_rht_seed_len
 );


#endif // GF_SYN_LOCALALIGNERINTERFACE_H
