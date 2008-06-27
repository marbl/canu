
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
#ifndef IDBConnection_HH
#define IDBConnection_HH

extern "C" {
   #include "AS_global.h"
   #include "AS_UTL_Hash.h"
}

namespace AS_ARD {
   class IDBConnection {
      private:
         // disable copy constructor
         IDBConnection(IDBConnection &);

      public:
         IDBConnection() {};
         virtual ~IDBConnection() {};

         static const int MAX_STR_LEN              = 4096;

         virtual bool sqlCommand(const char * command) = 0;
         virtual uint64 getCount(const char * tableName) = 0;
         virtual uint64 getLast(const char * columnName, const char * tableName) = 0;

         virtual void checkError(int32 ret, const char * str) = 0;
         virtual bool populateHash(
                     HashTable_AS * hash,
                     const char * keyColumn,
                     const char * valColumn,
                     const char * tableName,
                     uint64 assemblyID) = 0;
   };
};

#endif // IDBConnection_HH
