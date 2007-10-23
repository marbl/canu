
/**************************************************************************
 * Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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
/*************************************************
* Module:  ReadIUMs.cc
* Description:
*   Read the contents of an IUM/CGB file
* 
*    Programmer:  E. Venter
*       Started:  15 October 2007
* 
* Assumptions:
* 
* Notes:
*
*************************************************/

/* RCS info
 * $Id: ReadIUMs.cc,v 1.3 2007-10-23 14:42:29 eliv Exp $
 * $Revision: 1.3 $
*/

static const char CM_ID[] = "$Id: ReadIUMs.cc,v 1.3 2007-10-23 14:42:29 eliv Exp $";

//  System include files

#include"AS_BOG_UnitigGraph.hh"
#include"AS_BOG_MateChecker.hh"

using std::cout;
using std::endl;
using AS_BOG::UnitigGraph;
using AS_BOG::LongestHighIdent;

//  Local include files
extern "C" {
}

int  main
    (int argc, char * argv [])

{
   // Get path/name of file from command line
   const char* IUM_File       = argv[1];
   const char* GKP_Store_Path = argv[2];

   AS_BOG::MateChecker mateChecker;
   int numFrgsInGKP = mateChecker.readStore(GKP_Store_Path);

   LongestHighIdent bog( 15 );
   UnitigGraph tigGraph( &bog );

   tigGraph.readIUMsFromFile( IUM_File, numFrgsInGKP );
   AS_BOG::UnitigsConstIter iter = tigGraph.unitigs->begin();
   for(; iter != tigGraph.unitigs->end(); iter++)
   {
       cout << **iter << endl;
   }

   return  0;
}
