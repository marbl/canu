
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
#ifdef SYBASE

#ifndef Sybase_H
#define Sybase_H

#include <iostream>
#include "IDBOutput.hh"
#include "IDBConnection.hh"

extern "C" {
   #include <ctpublic.h>
   #include "AS_UTL_Hash.h"
}

namespace AS_ARD {
   class Sybase : public IDBConnection {
      private:
         // disable copy constructor
         Sybase(Sybase &);
      
         CS_COMMAND * sendCommand(CS_CHAR * command);
         static const int DB_SWITCH_MSG_NUMBER     = 69;
         static const int DB_SWITCH_MSG_SEVERITY   = 22;

      protected:                  
         bool sqlCommand(const char * command);
         
         uint64 getCount(const char * tableName);
         uint64 getLast(const char * columnName, const char * tableName);
         bool populateHash(
                     HashTable_AS * hash, 
                     const char * keyColumn, 
                     const char * valColumn, 
                     const char * tableName,
                     uint64 assemblyID);

         void checkError(int32 ret, const char * str);
                  
         static CS_RETCODE CS_PUBLIC csmsg_cb(CS_CONTEXT *context, CS_CLIENTMSG *msg);
         static CS_RETCODE CS_PUBLIC clientmsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_CLIENTMSG *msg);
         static CS_RETCODE CS_PUBLIC servermsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_SERVERMSG *msg);

         CS_CONTEXT     *context;     /* Context structure     */
         CS_CONNECTION  *connection;  /* Connection structure. */

      public:
      
         Sybase( 
            const char * _server, 
            const char * _user,
            const char * _password,
            const char * _database);
         ~Sybase();
   };
}; 

#endif // Sybase_H
#endif //SYBASE
