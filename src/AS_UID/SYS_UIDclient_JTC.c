
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

static const char *rcsid = "$Id: SYS_UIDclient_JTC.c,v 1.9 2008-10-08 22:03:00 brianwalenz Exp $";

#include <string.h> // for memcpy
#include <assert.h> // for assert
#include <curl/curl.h>
#include <curl/types.h>
#include <curl/easy.h>

#include "SYS_UIDcommon.h"


struct JTC_GUIDMemoryStruct {
  char *memory;
  size_t size;
};

const char JTC_GUID_URL[] = "http://guid.jtc.jcvsf.org:8080/guid/GuidClientServer";
char * guidServerNames = NULL;
char * guidNamespace = NULL;

#define JTC_GUID_REQUEST_URL_MAX_SIZE 2048
#define JTC_GUID_HTTP_RESPONSE_MAX_SIZE 4096
#define JTC_GUID_NUM_BUFFER_SIZE 100
#define JTC_GUID_SEPARATOR ","
#define JTC_GUID_MAX_NUM_URL 10

void SYS_UIDset_euid_server(const char * servers)
{
  guidServerNames=strdup(servers);
  assert(guidServerNames != NULL);
}

void SYS_UIDset_euid_namespace(const char * namespaceName)
{
	if (namespaceName == NULL) { return; }

  	guidNamespace=strdup(namespaceName);
  	assert(guidNamespace != NULL);
}


int findGuidStartFromHttpString(char* httpString) {
  /* the first position after the 2nd return character */
  int returnCount = 0;
  int i;
  for (i=0;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
    if (httpString[i] == '\n') {
      returnCount++;
    }
    if (returnCount == 2) {
      return i+1;
    }
  }
  /* error */
  return 0;
}


int findGuidEndFromHttpString(char* httpString) {
  /* the first position before the 3nd return character */
  int returnCount = 0;
  int i;
  for (i=0;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
    if (httpString[i] == '\n') {
      returnCount++;
    }
    if (returnCount == 3) {
      return i-1;
    }
  }
  /* error */
  return 0;
}


size_t
JTC_GUIDWriteMemoryCallback(void *ptr, size_t size, size_t nmemb, void *data)
{
  register int realsize = size * nmemb;
  struct JTC_GUIDMemoryStruct *mem = (struct JTC_GUIDMemoryStruct *)data;

  mem->memory = (char *)safe_realloc(mem->memory, mem->size + realsize + 1);
  memcpy(&(mem->memory[mem->size]), ptr, realsize);
  mem->size += realsize;
  mem->memory[mem->size] = 0;

  return realsize;
}


uint64 getGUIDBlock(int guidRequestSize)
{
  uint64 guidStart = 0;

  CURL *curl_handle;
  int curl_response = -1;
  /* we only handle a limited max of servers that can be specified */
  //char guidRequest[JTC_GUID_REQUEST_URL_MAX_SIZE];
  char *currentURL = NULL;
  char guidRequest[JTC_GUID_REQUEST_URL_MAX_SIZE];
  char guidBuffer[JTC_GUID_REQUEST_URL_MAX_SIZE];
  char httpResponse[JTC_GUID_HTTP_RESPONSE_MAX_SIZE];
  char guidNumResponse[JTC_GUID_NUM_BUFFER_SIZE]; /* needs to comfortably fit a long as a string */
  int guidNumLength = 0;
  struct JTC_GUIDMemoryStruct chunk;
  int i;
  int guidPositionStart = 0;
  int guidPositionEnd = 0;

  chunk.memory=NULL; /* we expect realloc(NULL, size) to work */
  chunk.size = 0;    /* no data at this point */

  if (guidServerNames == NULL)
  {
    // Have to copy because a static string lives in a read-only section on alpha
    strncpy(guidBuffer, JTC_GUID_URL, JTC_GUID_REQUEST_URL_MAX_SIZE);
  }
  else
  {
    strncpy(guidBuffer, guidServerNames, JTC_GUID_REQUEST_URL_MAX_SIZE);
  }

  currentURL = strtok(guidBuffer, ",");

  while (currentURL != NULL && curl_response != CURLE_OK)
  {
		sprintf(guidRequest, "%s?Request=GET&Size=%d", currentURL, guidRequestSize);

		if (guidNamespace != NULL) {
			sprintf(guidRequest, "%s&Namespace=%s", guidRequest, guidNamespace);
		}

		currentURL = strtok(NULL, ",");

		curl_global_init(CURL_GLOBAL_ALL);

		/* init the curl session */
		curl_handle = curl_easy_init();

		/* specify URL to get */
		curl_easy_setopt(curl_handle, CURLOPT_URL, guidRequest);

		/* send all data to this function  */
		curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, JTC_GUIDWriteMemoryCallback);

		/* we pass our 'chunk' struct to the callback function */
		curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, (void *)&chunk);

		/* get it! */
		curl_response = curl_easy_perform(curl_handle);

	    /* cleanup curl stuff */
	    curl_easy_cleanup(curl_handle);

		if (curl_response != CURLE_OK) {
			safe_free(chunk.memory);
			continue;
		}

		if (chunk.size >= JTC_GUID_HTTP_RESPONSE_MAX_SIZE) {
		  safe_free(chunk.memory);
		  continue;//return 0;
		}

		memcpy(httpResponse, chunk.memory, chunk.size);
		httpResponse[chunk.size] = '\0';

		/* HTTP response of this form is assumed:   */
		/* <html><title>SUCCESS</title><body>       */
		/* <h2>Guid Start:</h2>                     */
		/* 1089045040000                            */
		/* </body></html>                           */
		if (strncmp(httpResponse+13,"SUCCESS",7)==0) {
	      guidPositionStart = findGuidStartFromHttpString(httpResponse);
		  guidPositionEnd = findGuidEndFromHttpString(httpResponse);
		  guidNumLength = guidPositionEnd - guidPositionStart;
		  if (guidPositionStart == 0 || guidPositionEnd <= guidPositionStart) {
		    /* error */
		    safe_free(chunk.memory);
		    continue;//return 0;
		  }
		  memcpy(guidNumResponse, httpResponse + guidPositionStart, guidNumLength);
		  for (i=guidNumLength;i<JTC_GUID_NUM_BUFFER_SIZE;i++) {
		    guidNumResponse[i] = '\0';
		  }
		  guidStart = strtoull(guidNumResponse,NULL,10);
		} else {
		  /* error */
		  safe_free(chunk.memory);
		  continue;//return 0;
		}

		safe_free(chunk.memory);
  }

  return guidStart;
}
