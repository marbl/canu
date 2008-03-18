
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

#include <iostream>
#include <string>

#include "Sybase.hh"

using AS_ARD::Sybase;

#define EX_CTLIB_VERSION    CS_VERSION_150
#define EX_BLK_VERSION      BLK_VERSION_150 

CS_RETCODE CS_PUBLIC Sybase::csmsg_cb(CS_CONTEXT *context, CS_CLIENTMSG *msg) {
   std::cerr << "CSMSG (LIB) ERROR CALLBACK " << msg->msgstring << std::endl;
   
   return CS_SUCCEED;
}

CS_RETCODE CS_PUBLIC Sybase::clientmsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_CLIENTMSG *msg) {
   std::cerr << "CLIENT ERROR CALLBACK " << msg->msgstring << std::endl;
   
   return CS_SUCCEED;
}

CS_RETCODE CS_PUBLIC Sybase::servermsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_SERVERMSG *msg) {      
   // ignore db switching messages
   if (CS_NUMBER(msg->msgnumber) == DB_SWITCH_MSG_NUMBER && 
         CS_SEVERITY(msg->msgnumber) == DB_SWITCH_MSG_SEVERITY) {
      return CS_SUCCEED;
   }
   
   std::cerr << "SERVER ERROR CALLBACK " << msg->text << " from server " << msg->svrname << std::endl;
   
   return CS_SUCCEED;
}


Sybase::Sybase( 
            const char * _server,
            const char * _user,
            const char * _password,
            const char * _database
            ) {

   assert (_server != NULL);
   assert (_user != NULL);
   assert (_password != NULL);
   assert (_database != NULL);

   char        cmd[MAX_STR_LEN];
      
   context = (CS_CONTEXT *) NULL;
   CS_RETCODE ret = cs_ctx_alloc(EX_CTLIB_VERSION, &context);
   checkError(ret, "context allocation failed");
   ret = ct_init(context, EX_CTLIB_VERSION);
   checkError(ret, "ct_init failed");

   // setup error/msg handlers
   ret = cs_config(context, CS_SET, CS_MESSAGE_CB, (CS_VOID *)&csmsg_cb, CS_UNUSED, NULL);
   checkError(ret, "cs_config(CS_MESSAGE_CB) failed");
   ret = ct_callback(context, NULL, CS_SET, CS_CLIENTMSG_CB, (CS_VOID *) &clientmsg_cb);
   checkError(ret,"ct_callback for client messages failed");
   ret = ct_callback(context, NULL, CS_SET, CS_SERVERMSG_CB,(CS_VOID *) &servermsg_cb);
   checkError(ret,"ct_callback for server messages failed");

   // connect to server
   ret = ct_con_alloc(context, &connection);
   checkError(ret, "ct_con_alloc() failed");

   // set credentials
   ret = ct_con_props(connection, CS_SET, CS_USERNAME, (CS_VOID *)_user, strlen(_user), NULL);
   checkError(ret, "Could not set user name");

   ret = ct_con_props(connection, CS_SET, CS_PASSWORD, (CS_VOID *)_password, strlen(_password), NULL);
   checkError(ret, "Could not set password");

   ret = ct_connect(connection, (CS_CHAR *)_server, strlen(_server));
   checkError(ret, "Could not connect!");

   // switch to appropriate DB
   sprintf(cmd, "use %s", _database);
   assert(sqlCommand(cmd));

}

CS_COMMAND * Sybase::sendCommand(CS_CHAR * command) {
   CS_COMMAND     *cmd = NULL;
   CS_RETCODE      ret = 0;
      
   ret = ct_cmd_alloc(connection, &cmd);
   checkError(ret, "ct_cmd_alloc() failed");

   ret = ct_command(cmd, CS_LANG_CMD,
            command,
            CS_NULLTERM, CS_UNUSED); 
   checkError(ret, "ct_command() failed");

   // send command
   ret = ct_send(cmd);
   checkError(ret, "ct_send() failed");   

   return cmd;
}

bool Sybase::populateHash(
   HashTable_AS * hash, 
   const char * keyColumn, 
   const char * valColumn, 
   const char * tableName,
   uint64 assemblyID)
{
   CS_BIGINT   valCol = 0;
   CS_CHAR     keyCol[32];
   CS_INT      result_type = 0;
   CS_INT      count;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_COMMAND *cmd = NULL;
   CS_DATAFMT  column[2];
   
   char        command[Sybase::MAX_STR_LEN];
   std::string theTable = tableName;
   std::transform(theTable.begin(), theTable.end(), theTable.begin(), tolower);
   
   sprintf(command, "SELECT %s, %s FROM %s WHERE %s_AssemblyID = %d\n", keyColumn, valColumn, tableName, theTable.c_str(), assemblyID);
   cmd = sendCommand(command);

   memset(keyCol, 0, 32);
   while ((result_ret = ct_results(cmd, &result_type)) == CS_SUCCEED)
   {
      switch (static_cast<int>(result_type))
      {
         case CS_ROW_RESULT:
            column[0].datatype = CS_CHAR_TYPE;
            column[0].format = CS_FMT_NULLTERM;
            column[0].count = 1;
            column[0].maxlength = 32;
            column[0].locale = NULL;
            
            column[1].datatype = CS_BIGINT_TYPE;
            column[1].format = CS_FMT_UNUSED;
            column[1].count = 1;
            column[1].locale = NULL;
            
            ret = ct_bind(cmd, 1, &column[0], &keyCol, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 2, &column[1], &valCol, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");
            
            while(((ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, &count)) == CS_SUCCEED) 
                   || (ret == CS_ROW_FAIL)) {
               if( ret == CS_ROW_FAIL ) {
                  std::cerr << "Error retrieving row" << std::endl;
               }
               
               InsertInHashTable_AS(hash, strtoull(keyCol, NULL, 10), 0, (uint64)valCol, 0);
            }

            break;
         case CS_END_DATA:
         case CS_CMD_DONE:
            break;
         case CS_CMD_FAIL:
         case CS_CMD_SUCCEED:     
         default:
            ret = ct_cmd_drop(cmd);
            checkError(CS_FAIL, "ct_results returned unexpected result type");
            break;
      }
   }
   switch(static_cast<int>(result_ret)) {
      case CS_END_RESULTS:
         break;
      case CS_FAIL:
      default:
         ret = ct_cmd_drop(cmd);
         checkError(CS_FAIL, "ct_results() returned CS_FAIL.");
         break;
   }
   ret = ct_cmd_drop(cmd);
   checkError(ret, "cmd drop failed.");

   return true;
}

uint64 Sybase::getCount(const char * tableName) {
   CS_BIGINT   count = 0;
   CS_INT      result_type = 0;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_COMMAND *cmd = NULL;
   CS_DATAFMT  column;
   
   char        command[MAX_STR_LEN];
   
   sprintf(command, "SELECT COUNT(*) FROM %s\n", tableName);
   cmd = sendCommand((CS_CHAR *)command);

   while ((result_ret = ct_results(cmd, &result_type)) == CS_SUCCEED)
   {
      switch (static_cast<int>(result_type))
      {
         case CS_ROW_RESULT:
            column.datatype = CS_BIGINT_TYPE;
            column.format = CS_FMT_UNUSED;
            column.count = 1;
            column.locale = NULL;

            ret = ct_bind(cmd, 1, &column, &count, NULL, NULL);
            checkError(ret,"ct_bind() for COUNT(*) failed");
            
            // fetch the actual falue
            ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, NULL);
            if (static_cast<int>(ret) != CS_SUCCEED) {
               checkError(CS_FAIL, "ct_fetch failed");
            }
            
            // fetch one more time should mean we are at the end
            ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, NULL);
            if (static_cast<int>(ret) != CS_END_DATA) {
                 checkError(CS_FAIL, "ct_fetch failed");
            }
            break;
         case CS_END_DATA:
         case CS_CMD_DONE:
            break;
         case CS_CMD_FAIL:
         case CS_CMD_SUCCEED:     
         default:
            ret = ct_cmd_drop(cmd);
            checkError(CS_FAIL, "ct_results returned unexpected result type");
            break;
      }
   }
   switch(static_cast<int>(result_ret)) {
      case CS_END_RESULTS:
         break;
      case CS_FAIL:
      default:
         ret = ct_cmd_drop(cmd);
         checkError(CS_FAIL, "ct_results() returned CS_FAIL.");
         break;
   }
   ret = ct_cmd_drop(cmd);
   checkError(ret, "cmd drop failed.");

   return static_cast<uint64>(count);
}

uint64 Sybase::getLast(const char * columnName, const char * tableName) {
   CS_BIGINT   val = 0;
   CS_INT      result_type = 0;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_COMMAND *cmd = NULL;
   CS_DATAFMT  column;

   char        command[MAX_STR_LEN];
   
   sprintf(command, "SELECT MAX(%s) FROM %s\n", columnName, tableName);
   cmd = sendCommand((CS_CHAR *)command);

   while ((result_ret = ct_results(cmd, &result_type)) == CS_SUCCEED)
   {
      switch (static_cast<int>(result_type))
      {
         case CS_ROW_RESULT:
            column.datatype = CS_BIGINT_TYPE;
            column.format = CS_FMT_UNUSED;
            column.count = 1;
            column.locale = NULL;

            ret = ct_bind(cmd, 1, &column, &val, NULL, NULL);
            checkError(ret,"ct_bind() for MAX(*) failed");
            
            // fetch the actual falue
            ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, NULL);
            if (static_cast<int>(ret) != CS_SUCCEED) {
               checkError(CS_FAIL, "ct_fetch failed");
            }
            
            // fetch one more time should mean we are at the end
            ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, NULL);
            if (static_cast<int>(ret) != CS_END_DATA) {
                 checkError(CS_FAIL, "ct_fetch failed");
            }
            break;
         case CS_END_DATA:
         case CS_CMD_DONE:
            break;
         case CS_CMD_FAIL:
            std::cerr << "The error I got is " << result_type << std::endl;            
            break;
         case CS_CMD_SUCCEED:     
         default:
            ret = ct_cmd_drop(cmd);
            checkError(CS_FAIL, "ct_results returned unexpected result type");
            break;
      }
   }
   switch(static_cast<int>(result_ret)) {
      case CS_END_RESULTS:
         break;
      case CS_FAIL:
      default:
         ret = ct_cmd_drop(cmd);
         checkError(CS_FAIL, "ct_results() returned CS_FAIL.");
         break;
   }
   ret = ct_cmd_drop(cmd);
   checkError(ret, "cmd drop failed.");

   return static_cast<uint64>(val);
}

bool Sybase::sqlCommand(const char * command) {
   CS_COMMAND * cmd = sendCommand((CS_CHAR *)command);
   CS_INT result_type;
   CS_RETCODE ret;
   bool success = false;

//std::cerr << "Running command is " << command << std::endl;   
   // check results
   while ((ret = ct_results(cmd, &result_type)) == CS_SUCCEED) {
      switch (static_cast<int>(result_type)) {
         case CS_CMD_SUCCEED:
            success = true;
            break;
         case CS_ROW_RESULT:
            std::cerr << "Cannot run interactive command using function sqlCommand, please use sqlQueryCommand instead" << std::endl;
            break;
         case CS_CMD_FAIL:
            /*
            ** The server encountered an error while
            ** processing our command. These errors will be
            ** displayed by the server-message callback that
            ** we installed earlier.
            */
            break;
         case CS_CMD_DONE:
            /*
            ** The logical command has been completely
            ** processed.
            */
            break;
         default:
            break;
      } 
   }
   
   ret = ct_cmd_drop(cmd);
   checkError(ret, "ct_cmd_drop failed");
   
   return success;
}

void Sybase::checkError(int32 ret, const char * str) {
   if (static_cast<int>(ret) != CS_SUCCEED) {
      std::cerr << "Fatal error: " << str << std::endl;

      if (connection != NULL) {
         ct_close(connection, CS_FORCE_EXIT);
         ct_con_drop(connection);
      }
      
      if (context != NULL) {
         ct_exit(context, CS_FORCE_EXIT);
         cs_ctx_drop(context);
      }
      exit(1);
   }
}

Sybase::~Sybase() {
   ct_close(connection, CS_UNUSED);
   ct_exit(context, CS_UNUSED);
   
   ct_con_drop(connection);
   cs_ctx_drop(context);
}

#endif //SYBASE
