
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
#include <iostream>

/*************************************************************************/
/* Local include files */
/*************************************************************************/

extern "C" {
   #include <unistd.h>
   #include "AS_global.h"
   #include "SYS_UIDclient.h"
}

#include "AS_ARD_database.hh"
#include "IDBOutput.hh"
#include "DBTextOutput.hh"

#ifdef SYBASE
   #include "BCPOutput.hh"
   using AS_ARD::BCPOutput;
#endif

using AS_ARD::AS_ARD_database;
using AS_ARD::IDBOutput;
using AS_ARD::DBTextOutput;

int main(int argc, char ** argv)
{
   char * asmFilename       = NULL;
   char * prefix            = NULL;
   char * server            = NULL;
   char * user              = NULL;
   char * pass              = NULL;
   char * bcp               = NULL;
   char * database          = NULL;
  
   AS_ARD_database * asmDB  = NULL;
   IDBOutput * out          = NULL;
   uint64 uidStart          = 0; 
  
   // parse command line
   {
      int ch,errflg=FALSE;
      while (!errflg && ((ch = getopt(argc, argv, "a:E:s:n:o:S:U:B:D:")) != EOF))
      {
         switch(ch) 
         {
            case 'a':
               asmFilename = optarg;
               break;
            case 'E':       
               SYS_UIDset_euid_server(optarg);
               break;
            case 's':
               uidStart = strtoul(optarg, NULL, 10);
               break;
            case 'o':
               prefix = optarg;
               break;
            case 'n':
               SYS_UIDset_euid_namespace(optarg);
               break;
            case 'B':
               bcp = optarg;
               break;
            case 'S':
               server = optarg;
               break;
            case 'U':
               user = optarg;
               break;
            case 'D':
               database = optarg;
               break;
            default:
   	        fprintf(stderr,"Unrecognized option -%c\n",optopt);
	           errflg++;
            break;
         }
      }
   }
   #ifdef SYBASE
      // check command line requirements
      if(asmFilename == NULL || prefix == NULL || server == NULL || database == NULL || user == NULL)
      {
         std::cerr << "Usage: " << argv[0] << " -a asmFilename -o prefix -S server_name -D database_name -U user_name [-B bcp_path -E EUID server -n EUID namespace -s startEUID]\n";            
         exit(-1);
      }

      // inputs a password without echoing
      // this function is deprecated, need another method
      char * password = getpass("Please enter your password:");
      
      out = new BCPOutput(prefix, server, database, user, password, bcp);
      
      // clear password
      memset(password, 0, strlen(password));
   #else
      if(asmFilename == NULL)
      {
         std::cerr << "Usage: " << argv[0] << " -a asmFilename [ -E EUID server -n EUID namespace -s startEUID]\n";            
         exit(-1);
      }

      out = new DBTextOutput();
   #endif
   
   asmDB = new AS_ARD_database(out);
   assert(asmDB != NULL);
  
   // read asm file into data structures
   {
      UIDserver * uids = UIDserverInitialize(100, uidStart);
      FILE * fi = fopen(asmFilename, "r");
      assert(fi != NULL);
    
      if (asmDB->LoadDatabaseFromASMFile(fi, uids) != true) {
         delete asmDB;
         assert(0);
      }
      fclose(fi);
   }
      
   delete asmDB;
  
   return 0;
}
