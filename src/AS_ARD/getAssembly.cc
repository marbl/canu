
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
/*************************************************************************/
/* Local include files */
/*************************************************************************/
#ifdef SYBASE
extern "C" {
   #include <ctpublic.h>
   #include <unistd.h>
}
#endif

extern "C" {   
   #include "AS_global.h"
   #include "SYS_UIDclient.h"
   #include "AS_MSG_pmesg.h" 
   #include "AS_UTL_Hash.h"
}

#include <iostream>
#include <string>

#ifdef SYBASE
#define EX_CTLIB_VERSION    CS_VERSION_150
#define EX_BLK_VERSION      BLK_VERSION_150 
#define MAX_STR_LEN 4096
#define DB_SWITCH_MSG_NUMBER 69
#define DB_SWITCH_MSG_SEVERITY 22

CS_CONTEXT     *context;     /* Context structure     */
CS_CONNECTION  *connection;  /* Connection structure. */
FILE * file;

void checkError(CS_RETCODE ret, const char * str) {
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
      exit(-1);
   }
}

CS_RETCODE CS_PUBLIC csmsg_cb(CS_CONTEXT *context, CS_CLIENTMSG *msg) {
   std::cerr << "CSMSG (LIB) ERROR CALLBACK " << msg->msgstring << std::endl;
   
   return CS_SUCCEED;
}

CS_RETCODE CS_PUBLIC clientmsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_CLIENTMSG *msg) {
   std::cerr << "CLIENT ERROR CALLBACK " << msg->msgstring << std::endl;
   
   return CS_SUCCEED;
}

CS_RETCODE CS_PUBLIC servermsg_cb(CS_CONTEXT *context, CS_CONNECTION *connection, CS_SERVERMSG *msg) {      
   // ignore db switching messages
   if (CS_NUMBER(msg->msgnumber) == DB_SWITCH_MSG_NUMBER && 
         CS_SEVERITY(msg->msgnumber) == DB_SWITCH_MSG_SEVERITY) {
      return CS_SUCCEED;
   }
   
   std::cerr << "SERVER ERROR CALLBACK " << msg->text << " from server " << msg->svrname << std::endl;
   
   return CS_SUCCEED;
}

CS_COMMAND * sendCommand(CS_CHAR * command) {
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

bool sqlCommand(CS_CHAR * command) {
   CS_COMMAND * cmd = sendCommand(command);
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

void getMDIs() {
   CS_CHAR EUID[32];
   CS_INT  CIID;
   CS_FLOAT mean;
   CS_FLOAT std;
   CS_INT min;
   CS_INT max;
   
   CS_COMMAND *cmd = NULL;
   CS_INT      result_type = 0;
   CS_INT      count;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_DATAFMT  column[6];
   
   char        command[MAX_STR_LEN];

   sprintf(command, "SELECT mdi_EUID, mdi_CIID, mdi_mea, mdi_std, mdi_min, mdi_max FROM MDI\n");
std::cerr << "Running command " << command << std::endl;
   cmd = sendCommand(command);

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

            column[1].datatype = CS_INT_TYPE;
            column[1].format = CS_FMT_UNUSED;
            column[1].count = 1;
            column[1].locale = NULL;
            
            column[2].datatype = CS_FLOAT_TYPE;
            column[2].format = CS_FMT_UNUSED;
            column[2].count = 1;
            column[2].locale = NULL;

            column[3].datatype = CS_FLOAT_TYPE;
            column[3].format = CS_FMT_UNUSED;
            column[3].count = 1;
            column[3].locale = NULL;

            column[4].datatype = CS_INT_TYPE;
            column[4].format = CS_FMT_UNUSED;
            column[4].count = 1;
            column[4].locale = NULL;

            column[5].datatype = CS_INT_TYPE;
            column[5].format = CS_FMT_UNUSED;
            column[5].count = 1;
            column[5].locale = NULL;
           
            ret = ct_bind(cmd, 1, &column[0], &EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 2, &column[1], &CIID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 3, &column[2], &mean, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 4, &column[3], &std, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 5, &column[4], &min, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 6, &column[5], &max, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            
            
            SnapMateDistMesg mdi;
            GenericMesg gen;
            gen.m = &mdi;
            gen.t = MESG_MDI;
            
            while(((ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, &count)) == CS_SUCCEED) 
                   || (ret == CS_ROW_FAIL)) {
               if( ret == CS_ROW_FAIL ) {
                  std::cerr << "Error retrieving row" << std::endl;
               }
               mdi.erefines = AS_UID_lookup(EUID, NULL);
               mdi.irefines = CIID;                  
               mdi.mean = mean;
               mdi.stddev = std;
               mdi.min = min;
               mdi.max = max;
               mdi.num_buckets = 0;
               mdi.histogram = NULL;
               
               WriteProtoMesg_AS(file,  &gen);
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
}

void getAFGs(std::string utgEUIDs) {
   CS_CHAR afg_EUID[32];
   CS_INT  afg_CIID;
   CS_CHAR afg_mst;
   CS_BIT afg_chi;
   CS_BIT afg_cha;
   CS_SMALLINT afg_clr1;
   CS_SMALLINT afg_clr2;
   
   CS_COMMAND *cmd = NULL;
   CS_INT      result_type = 0;
   CS_INT      count;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_DATAFMT  column[7];
   
   char        command[MAX_STR_LEN];

   sprintf(command, "SELECT afg_EUID, afg_CIID, afg_mst, afg_chi, afg_cha, afg_clr1, afg_clr2 from UTG join MPS on mps_utg_MSG_ID = utg_MSG_ID join AFG on afg_MSG_ID = mps_afg_MSG_ID where utg_MSG_ID IN (%s)\n", utgEUIDs.c_str());
std::cerr << "Running command " << command << std::endl;
   cmd = sendCommand(command);

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

            column[1].datatype = CS_INT_TYPE;
            column[1].format = CS_FMT_UNUSED;
            column[1].count = 1;
            column[1].locale = NULL;
            
            column[2].datatype = CS_CHAR_TYPE;
            column[2].format = CS_FMT_UNUSED;
            column[2].count = 1;
            column[2].maxlength = 1;
            column[2].locale = NULL;

            column[3].datatype = CS_BIT_TYPE;
            column[3].format = CS_FMT_UNUSED;
            column[3].count = 1;
            column[3].locale = NULL;

            column[4].datatype = CS_BIT_TYPE;
            column[4].format = CS_FMT_UNUSED;
            column[4].count = 1;
            column[4].locale = NULL;

            column[5].datatype = CS_SMALLINT_TYPE;
            column[5].format = CS_FMT_UNUSED;
            column[5].count = 1;
            column[5].locale = NULL;

            column[6].datatype = CS_SMALLINT_TYPE;
            column[6].format = CS_FMT_UNUSED;
            column[6].count = 1;
            column[6].locale = NULL;
           
            ret = ct_bind(cmd, 1, &column[0], &afg_EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 2, &column[1], &afg_CIID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 3, &column[2], &afg_mst, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 4, &column[3], &afg_chi, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 5, &column[4], &afg_cha, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 6, &column[5], &afg_clr1, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 7, &column[6], &afg_clr2, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            
            
            AugFragMesg afg;
            GenericMesg gen;
            gen.m = &afg;
            gen.t = MESG_AFG;
            
            while(((ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, &count)) == CS_SUCCEED) 
                   || (ret == CS_ROW_FAIL)) {
               if( ret == CS_ROW_FAIL ) {
                  std::cerr << "Error retrieving row" << std::endl;
               }
               afg.eaccession = AS_UID_lookup(afg_EUID, NULL);
               afg.iaccession = afg_CIID;   
               afg.mate_status = (MateStatType)afg_mst;
               afg.chimeric = afg_chi;
               afg.chaff = afg_cha;
               afg.clear_rng.bgn = afg_clr1;
               afg.clear_rng.end = afg_clr2;
               
               WriteProtoMesg_AS(file,  &gen);
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
}

void getUTGs(std::string utgEUIDs) {
   CS_CHAR EUID[32];
   CS_INT  CIID;
   CS_CHAR src[255];
   CS_FLOAT mhp;
   CS_FLOAT cov;
   CS_CHAR sta;
   CS_INT len;
   CS_BIT utg_for;
   CS_INT nfr;
   CS_CHAR afg_EUID[32];
   CS_CHAR mps_typ;
   CS_CHAR mps_src[255];
   CS_INT mps_pos1;
   CS_INT mps_pos2;
   CS_CHAR mps_del[1000];
   
   CS_COMMAND *cmd = NULL;
   CS_INT      result_type = 0;
   CS_INT      count;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_DATAFMT  column[15];
   
   int seen_frg = 0;
   
   char        command[MAX_STR_LEN];

   sprintf(command, "SELECT utg_EUID, utg_CIID, utg_src, utg_mhp, utg_cov, utg_sta, utg_len, utg_for, utg_nfr, afg_EUID, mps_type, mps_src, mps_pos1, mps_pos2, mps_del from UTG join MPS on mps_utg_MSG_ID = utg_MSG_ID join AFG on afg_MSG_ID = mps_afg_MSG_ID where utg_MSG_ID IN (%s) ORDER BY utg_MSG_ID ASC\n", utgEUIDs.c_str());
std::cerr << "Running command " << command << std::endl;
   cmd = sendCommand(command);

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

            column[1].datatype = CS_INT_TYPE;
            column[1].format = CS_FMT_UNUSED;
            column[1].count = 1;
            column[1].locale = NULL;

            column[2].datatype = CS_CHAR_TYPE;
            column[2].format = CS_FMT_UNUSED;
            column[2].count = 1;
            column[2].maxlength = 255;
            column[2].locale = NULL;

            column[3].datatype = CS_FLOAT_TYPE;
            column[3].format = CS_FMT_UNUSED;
            column[3].count = 1;            
            column[3].locale = NULL;

            column[4].datatype = CS_FLOAT_TYPE;
            column[4].format = CS_FMT_UNUSED;
            column[4].count = 1;            
            column[4].locale = NULL;

            column[5].datatype = CS_CHAR_TYPE;
            column[5].format = CS_FMT_UNUSED;
            column[5].count = 1;
            column[5].maxlength = 1;
            column[5].locale = NULL;

            column[6].datatype = CS_INT_TYPE;
            column[6].format = CS_FMT_UNUSED;
            column[6].count = 1;
            column[6].locale = NULL;

            column[7].datatype = CS_BIT_TYPE;
            column[7].format = CS_FMT_UNUSED;
            column[7].count = 1;
            column[7].locale = NULL;

            column[8].datatype = CS_INT_TYPE;
            column[8].format = CS_FMT_UNUSED;
            column[8].count = 1;
            column[8].locale = NULL;

            column[9].datatype = CS_CHAR_TYPE;
            column[9].format = CS_FMT_UNUSED;
            column[9].count = 1;
            column[9].maxlength = 32;
            column[9].locale = NULL;
            
            column[10].datatype = CS_CHAR_TYPE;
            column[10].format = CS_FMT_UNUSED;
            column[10].count = 1;
            column[10].maxlength = 1;
            column[10].locale = NULL;

            column[11].datatype = CS_CHAR_TYPE;
            column[11].format = CS_FMT_UNUSED;
            column[11].count = 1;
            column[11].maxlength = 255;
            column[11].locale = NULL;

            column[12].datatype = CS_INT_TYPE;
            column[12].format = CS_FMT_UNUSED;
            column[12].count = 1;
            column[12].locale = NULL;

            column[13].datatype = CS_INT_TYPE;
            column[13].format = CS_FMT_UNUSED;
            column[13].count = 1;
            column[13].locale = NULL;

            column[14].datatype = CS_CHAR_TYPE;
            column[14].format = CS_FMT_UNUSED;
            column[14].count = 1;
            column[14].maxlength = 1000;
            column[14].locale = NULL;
           
            ret = ct_bind(cmd, 1, &column[0], &EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 2, &column[1], &CIID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 3, &column[2], &src, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 4, &column[3], &mhp, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 5, &column[4], &cov, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 6, &column[5], &sta, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 7, &column[6], &len, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 8, &column[7], &utg_for, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 9, &column[8], &nfr, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 10, &column[9], &afg_EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 11, &column[10], &mps_typ, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 12, &column[11], &mps_src, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 13, &column[12], &mps_pos1, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 14, &column[13], &mps_pos2, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 15, &column[14], &mps_del, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            
            
            SnapUnitigMesg utg;
            utg.f_list = NULL;
            utg.num_frags = 0;
            
            GenericMesg gen;
            gen.m = &utg;
            gen.t = MESG_UTG;
                        
            while(((ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, &count)) == CS_SUCCEED) 
                   || (ret == CS_ROW_FAIL)) {
               if( ret == CS_ROW_FAIL ) {
                  std::cerr << "Error retrieving row" << std::endl;
               }               

               if (seen_frg == 0) {
                  // init the utg
                  for (int i = 0; i < utg.num_frags; i++) {
                     delete[] utg.f_list[i].delta;
                  }
                  if (utg.f_list != NULL) {
                     delete[] utg.f_list;
                  }
                  utg.num_frags = nfr;
                  utg.num_vars = 0;
                  utg.v_list = NULL;
                  utg.f_list = new SnapMultiPos[nfr];

std::cerr<< "The UID I got is " << EUID << std::endl;                  
                  utg.eaccession = AS_UID_lookup(EUID, NULL);
std::cerr<< "The converted form is I got is " << AS_UID_toString(utg.eaccession) << std::endl;
                  utg.iaccession = CIID;
                  #ifdef AS_ENABLE_SOURCE
                     utg.source = src;
                  #endif
                  utg.coverage_stat = cov;
                  utg.status = (UnitigStatus) sta;
                  utg.length = len;
                  utg.consensus = "";
                  utg.quality = "";
                  utg.forced = utg_for;
                  
                  memset(EUID, 0, 32);
                  memset(src, 0, 255);
               }
               
               utg.f_list[seen_frg].eident = AS_UID_lookup(afg_EUID, NULL);
               utg.f_list[seen_frg].type = (FragType)mps_typ;
               #ifdef AS_ENABLE_SOURCE
                  utg.f_list[seen_frg].source = strdup(mps_src);
               #endif
               utg.f_list[seen_frg].position.bgn = mps_pos1;
               utg.f_list[seen_frg].position.end = mps_pos2;
               utg.f_list[seen_frg].delta_length = 0;
               utg.f_list[seen_frg].delta = NULL;               
               
               // count the spaces, that is our size of delta               
//std::cerr << "From the server delta is " << mps_del << std::endl;
               std::string delta = mps_del;
               if (delta.length() != 0) {                  
                  std::string::size_type start = 0;
                  std::string::size_type pos = 0;
                  utg.f_list[seen_frg].delta_length = 1;
                  
                  while ((pos = delta.find(" ", pos)) != std::string::npos) {
                     pos++;
                     utg.f_list[seen_frg].delta_length++;
                  }
//std::cerr << "The delta length is " << utg.f_list[seen_frg].delta_length << std::endl;
                  utg.f_list[seen_frg].delta = new int32[utg.f_list[seen_frg].delta_length];
                  
                  pos = 0;
                  int i = 0;
                  while ((pos = delta.find(" ", pos)) != std::string::npos) {
                     //std::cerr << "Found delta of " << delta.substr(start, pos-start-1) << std::endl;
                     utg.f_list[seen_frg].delta[i] = strtol(delta.substr(start, pos-start-1).c_str(), NULL, 10);
                     
                     pos++;
                     i++;
                     start = pos;
                  }
                  //std::cerr << "Found delta of " << delta.substr(start, delta.length()-start-1) << std::endl;
                  utg.f_list[seen_frg].delta[i] = strtol(delta.substr(start, delta.length()-start).c_str(), NULL, 10);
               }
               memset(afg_EUID, 0, 32);               
               memset(mps_src, 0, 255);
               memset(mps_del, 0, 1000);
               seen_frg++;
               
std::cerr << "The count is " << utg.num_frags << std::endl;
std::cerr << "The seen count is " << seen_frg << std::endl;
               if (seen_frg == utg.num_frags) {
std::cerr << "Outputting when count is " << seen_frg << std::endl;
                  WriteProtoMesg_AS(file, &gen);
                  seen_frg = 0;
               }
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
}

void getCCOs(std::string ccoEUIDs) {
   CS_CHAR EUID[32];
   CS_INT  CIID;
   CS_CHAR pla;
   CS_INT len;
   CS_BIT cco_for;
   CS_INT npc;
   CS_INT nou;
   CS_INT nvr;
   CS_CHAR afg_EUID[32];
   CS_CHAR mps_typ;
   CS_CHAR mps_src[255];
   CS_INT mps_pos1;
   CS_INT mps_pos2;
   CS_CHAR mps_del[1000];
   CS_CHAR utg_EUID[32];
   CS_CHAR ups_typ;
   CS_INT ups_pos1;
   CS_INT ups_pos2;
   CS_CHAR ups_del[1000];
   
   CS_COMMAND *cmd = NULL;
   CS_INT      result_type = 0;
   CS_INT      count;
   CS_RETCODE  ret = 0;
   CS_RETCODE  result_ret = 0;   
   CS_DATAFMT  column[19];
   
   int seen_frg = 0;
   int seen_var = 0;
   int seen_utg = 0;
   
   char        command[MAX_STR_LEN];

   sprintf(command, " select cco_EUID, cco_CIID, cco_pla, cco_len, cco_for, cco_npc, cco_nou, cco_nvr, afg_EUID, cco_mps_type, cco_mps_src, cco_mps_pos1, cco_mps_pos2, cco_mps_del, utg_EUID, ups_type, ups_pos1, ups_pos2, ups_del from CCO left outer join CCO_MPS on cco_mps_cco_MSG_ID = cco_MSG_ID left outer join UPS on ups_cco_MSG_ID = cco_MSG_ID left outer join VAR on var_cco_MSG_ID = cco_MSG_ID left outer join AFG on afg_MSG_ID = cco_mps_mid left outer join UTG on utg_MSG_ID = ups_utg_MSG_ID where cco_MSG_ID IN (%s) ORDER BY cco_MSG_ID ASC\n", ccoEUIDs.c_str());
std::cerr << "Running command " << command << std::endl;
   cmd = sendCommand(command);

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

            column[1].datatype = CS_INT_TYPE;
            column[1].format = CS_FMT_UNUSED;
            column[1].count = 1;
            column[1].locale = NULL;

            column[2].datatype = CS_CHAR_TYPE;
            column[2].format = CS_FMT_UNUSED;
            column[2].count = 1;
            column[2].maxlength = 1;
            column[2].locale = NULL;

            column[3].datatype = CS_INT_TYPE;
            column[3].format = CS_FMT_UNUSED;
            column[3].count = 1;            
            column[3].locale = NULL;

            column[4].datatype = CS_BIT_TYPE;
            column[4].format = CS_FMT_UNUSED;
            column[4].count = 1;            
            column[4].locale = NULL;

            column[5].datatype = CS_INT_TYPE;
            column[5].format = CS_FMT_UNUSED;
            column[5].count = 1;
            column[5].maxlength = 1;
            column[5].locale = NULL;

            column[6].datatype = CS_INT_TYPE;
            column[6].format = CS_FMT_UNUSED;
            column[6].count = 1;
            column[6].locale = NULL;

            column[7].datatype = CS_INT_TYPE;
            column[7].format = CS_FMT_UNUSED;
            column[7].count = 1;
            column[7].locale = NULL;

            column[8].datatype = CS_CHAR_TYPE;
            column[8].format = CS_FMT_UNUSED;
            column[8].count = 1;
            column[8].maxlength = 32;
            column[8].locale = NULL;
            
            column[9].datatype = CS_CHAR_TYPE;
            column[9].format = CS_FMT_UNUSED;
            column[9].count = 1;
            column[9].maxlength = 1;
            column[9].locale = NULL;

            column[10].datatype = CS_CHAR_TYPE;
            column[10].format = CS_FMT_UNUSED;
            column[10].count = 1;
            column[10].maxlength = 255;
            column[10].locale = NULL;

            column[11].datatype = CS_INT_TYPE;
            column[11].format = CS_FMT_UNUSED;
            column[11].count = 1;
            column[11].locale = NULL;

            column[12].datatype = CS_INT_TYPE;
            column[12].format = CS_FMT_UNUSED;
            column[12].count = 1;
            column[12].locale = NULL;

            column[13].datatype = CS_CHAR_TYPE;
            column[13].format = CS_FMT_UNUSED;
            column[13].count = 1;
            column[13].maxlength = 1000;
            column[13].locale = NULL;

            column[14].datatype = CS_CHAR_TYPE;
            column[14].format = CS_FMT_UNUSED;
            column[14].count = 1;
            column[14].maxlength = 32;
            column[14].locale = NULL;
            
            column[15].datatype = CS_CHAR_TYPE;
            column[15].format = CS_FMT_UNUSED;
            column[15].count = 1;
            column[15].maxlength = 1;
            column[15].locale = NULL;

            column[16].datatype = CS_INT_TYPE;
            column[16].format = CS_FMT_UNUSED;
            column[16].count = 1;
            column[16].locale = NULL;

            column[17].datatype = CS_INT_TYPE;
            column[17].format = CS_FMT_UNUSED;
            column[17].count = 1;
            column[17].locale = NULL;

            column[18].datatype = CS_CHAR_TYPE;
            column[18].format = CS_FMT_UNUSED;
            column[18].count = 1;
            column[18].maxlength = 1000;
            column[18].locale = NULL;

            ret = ct_bind(cmd, 1, &column[0], &EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 2, &column[1], &CIID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 3, &column[2], &pla, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");

            ret = ct_bind(cmd, 4, &column[3], &len, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 5, &column[4], &cco_for, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 6, &column[5], &npc, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 7, &column[6], &nou, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 8, &column[7], &nvr, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 9, &column[8], &afg_EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 10, &column[9], &mps_typ, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 11, &column[10], &mps_src, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 12, &column[11], &mps_pos1, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 13, &column[12], &mps_pos2, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 14, &column[13], &mps_del, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 15, &column[14], &utg_EUID, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 16, &column[15], &ups_typ, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 17, &column[16], &ups_pos1, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 18, &column[17], &ups_pos2, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            ret = ct_bind(cmd, 19, &column[18], &ups_del, NULL, NULL);
            checkError(ret,"ct_bind() for * failed");            

            SnapConConMesg ctg;
            ctg.pieces = NULL;
            ctg.num_pieces = 0;
            ctg.vars = NULL;
            ctg.num_vars = 0;
            ctg.unitigs = NULL;
            ctg.num_unitigs = 0;
            
            GenericMesg gen;
            gen.m = &ctg;
            gen.t = MESG_CCO;
                        
            while(((ret = ct_fetch(cmd, CS_UNUSED, CS_UNUSED, CS_UNUSED, &count)) == CS_SUCCEED) 
                   || (ret == CS_ROW_FAIL)) {
               if( ret == CS_ROW_FAIL ) {
                  std::cerr << "Error retrieving row" << std::endl;
               }               

               if (seen_frg == 0) {
                  // init the contig
                  for (int i = 0; i < ctg.num_pieces; i++) {
                     delete[] ctg.pieces[i].delta;
                  }
                  if (ctg.pieces != NULL) {
                     delete[] ctg.pieces;
                  }
                  for (int i = 0; i < ctg.num_unitigs; i++) {
                     delete[] ctg.unitigs[i].delta;
                  }
                  if (ctg.unitigs != NULL) {
                     delete[] ctg.unitigs;
                  }
                  ctg.num_pieces = npc;  
                  ctg.pieces = new SnapMultiPos[npc];
                  
                  ctg.num_unitigs = nou;
                  ctg.unitigs = new UnitigPos[nou];
                  
                  //ctg.num_vars = nvr;
                  //ctg.vars = new IntMultiVar[nvr];
                  ctg.num_vars = 0;
                  ctg.vars = NULL;
                  
                  ctg.eaccession = AS_UID_lookup(EUID, NULL);
                  ctg.iaccession = CIID;
                  ctg.placed = (ContigPlacementStatusType) pla;
                  ctg.length = len;
                  ctg.consensus = "";
                  ctg.quality = "";
                  ctg.forced = cco_for;
                  
                  memset(EUID, 0, 32);

                  seen_frg = seen_utg = seen_var = 0;
                  seen_var = ctg.num_vars;                  
               }
               
               if (seen_frg < ctg.num_pieces || seen_frg == 0) {
                  ctg.pieces[seen_frg].eident = AS_UID_lookup(afg_EUID, NULL);
                  ctg.pieces[seen_frg].type = (FragType)mps_typ;
                  #ifdef AS_ENABLE_SOURCE
                     ctg.pieces[seen_frg].source = strdup(mps_src);
                  #endif
                  ctg.pieces[seen_frg].position.bgn = mps_pos1;
                  ctg.pieces[seen_frg].position.end = mps_pos2;
                  ctg.pieces[seen_frg].delta_length = 0;
                  ctg.pieces[seen_frg].delta = NULL;               
                  
                  std::string delta = mps_del;
                  if (delta.length() != 0) {                  
                     std::string::size_type start = 0;
                     std::string::size_type pos = 0;
                     ctg.pieces[seen_frg].delta_length = 1;
                     
                     while ((pos = delta.find(" ", pos)) != std::string::npos) {
                        pos++;
                        ctg.pieces[seen_frg].delta_length++;
                     }
                     ctg.pieces[seen_frg].delta = new int32[ctg.pieces[seen_frg].delta_length];
                     
                     pos = 0;
                     int i = 0;
                     while ((pos = delta.find(" ", pos)) != std::string::npos) {
                        ctg.pieces[seen_frg].delta[i] = strtol(delta.substr(start, pos-start-1).c_str(), NULL, 10);
                        
                        pos++;
                        i++;
                        start = pos;
                     }
                     ctg.pieces[seen_frg].delta[i] = strtol(delta.substr(start, delta.length()-start).c_str(), NULL, 10);
                  }
                  memset(afg_EUID, 0, 32);               
                  memset(mps_src, 0, 255);
                  memset(mps_del, 0, 1000);

                  seen_frg++;
               }
               
               if (seen_utg < ctg.num_unitigs || seen_utg == 0) {
                  ctg.unitigs[seen_utg].eident = AS_UID_lookup(utg_EUID, NULL);                 
                  ctg.unitigs[seen_utg].type = (UnitigType)ups_typ;
                  ctg.unitigs[seen_utg].position.bgn = ups_pos1;
                  ctg.unitigs[seen_utg].position.end = ups_pos2;
                  ctg.unitigs[seen_utg].delta_length = 0;
                  ctg.unitigs[seen_utg].delta = NULL;               


                  std::string delta = ups_del;
                  if (delta.length() != 0) {                  
                     std::string::size_type start = 0;
                     std::string::size_type pos = 0;
                     ctg.unitigs[seen_utg].delta_length = 1;
                     
                     while ((pos = delta.find(" ", pos)) != std::string::npos) {
                        pos++;
                        ctg.unitigs[seen_utg].delta_length++;
                     }
                     ctg.unitigs[seen_utg].delta = new int32[ctg.unitigs[seen_utg].delta_length];
                     
                     pos = 0;
                     int i = 0;
                     while ((pos = delta.find(" ", pos)) != std::string::npos) {
                        ctg.unitigs[seen_utg].delta[i] = strtol(delta.substr(start, pos-start-1).c_str(), NULL, 10);
                        
                        pos++;
                        i++;
                        start = pos;
                     }
                     ctg.unitigs[seen_utg].delta[i] = strtol(delta.substr(start, delta.length()-start).c_str(), NULL, 10);
                  }
                  memset(utg_EUID, 0, 32);
                  memset(ups_del, 0, 1000);                  
                  
                  seen_utg++;
               }
               
               if (seen_var < ctg.num_vars || seen_var == 0) {
               }               

std::cerr << "The seen_frg is " << seen_frg << std::endl;
std::cerr << "The seen_utg is " << seen_utg << std::endl;
std::cerr << "The seen_var is " << seen_var << std::endl;
               if (seen_frg == ctg.num_pieces && seen_utg == ctg.num_unitigs && seen_var == ctg.num_vars) {
std::cerr << "Outputting " << std::endl;
                  WriteProtoMesg_AS(file,  &gen);
                  seen_frg = seen_utg = seen_var = 0;
               }               
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
}
#endif

int main(int argc, char ** argv)
{
   char * asmFilename       = NULL;
   char * prefix            = NULL;
   char * server            = NULL;
   char * user              = NULL;
   char * pass              = NULL;
   char * bcp               = NULL;
   char * database          = NULL;
  
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

      // read asm file into data structures
      UIDserver * uids = UIDserverInitialize(100, uidStart);

      file = fopen(asmFilename, "w");
      assert(file != NULL);

      CS_RETCODE  ret = 0;
   
      char        command[MAX_STR_LEN];
      
      context = (CS_CONTEXT *) NULL;
      ret = cs_ctx_alloc(EX_CTLIB_VERSION, &context);
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
      ret = ct_con_props(connection, CS_SET, CS_USERNAME, (CS_VOID *)user, strlen(user), NULL);
      checkError(ret, "Could not set user name");
   
      ret = ct_con_props(connection, CS_SET, CS_PASSWORD, (CS_VOID *)password, strlen(password), NULL);
      checkError(ret, "Could not set password");
   
      ret = ct_connect(connection, (CS_CHAR *)server, strlen(server));
      checkError(ret, "Could not connect!");
   
      // switch to appropriate DB
      sprintf(command, "use %s", database);   
   std::cerr << "Running command " << command << std::endl;
      assert(sqlCommand(command));
      
      getMDIs();
      getAFGs("12526464, 12526465");
      getUTGs("12526464, 12526465");
      getCCOs("5715228");
    
      fclose(file);
      return 0;
   #else
      std::cerr << "Only works on sybase right now" << std::endl;
      assert(0);
   #endif
}

