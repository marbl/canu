
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
 * $Id: ReadIUMs.cc,v 1.5 2007-12-05 23:46:58 brianwalenz Exp $
 * $Revision: 1.5 $
 */

static const char CM_ID[] = "$Id: ReadIUMs.cc,v 1.5 2007-12-05 23:46:58 brianwalenz Exp $";

//  System include files

#include"AS_BOG_UnitigGraph.hh"
#include"AS_BOG_MateChecker.hh"

using std::cout;
using std::endl;
using namespace AS_BOG;

//  Local include files
extern "C" {
}

int
main (int argc, char * argv []) {
    // Get path/name of file from command line
    const char* IUM_File       = argv[1];
    const char* GKP_Store_Path = argv[2];
    const char* OVL_Store_Path = argv[3];

    argc = AS_configure(argc, argv);

    MateChecker mateChecker;
    int numFrgsInGKP = mateChecker.readStore(GKP_Store_Path);

    OverlapStore* ovlStore = AS_OVS_openOverlapStore(OVL_Store_Path);

    BOG_Runner bogRunner(numFrgsInGKP);
    LongestHighIdent* bog = new LongestHighIdent( 1.5 );
    bogRunner.push_back( bog );
    bogRunner.processOverlapStream(ovlStore, GKP_Store_Path);

    UnitigGraph tigGraph( bog );

    tigGraph.readIUMsFromFile( IUM_File, numFrgsInGKP );
    BestEdgeCounts cnts = tigGraph.countInternalBestEdges();
    std::cerr << std::endl << "Overall best edge counts: dovetail " << cnts.dovetail
              << " oneWayBest " << cnts.oneWayBest << " neither " << cnts.neither
              << std::endl << std::endl;

    /*   AS_BOG::UnitigsConstIter iter = tigGraph.unitigs->begin();
         for(; iter != tigGraph.unitigs->end(); iter++)
         {
         cout << **iter << endl;
         }
    */

    delete bog;
    return  0;
}
