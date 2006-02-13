
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
static char CM_ID[] = "$Id: AS_PER_fragDB.c,v 1.6 2006-02-13 22:16:31 eliv Exp $";
/*************************************************************************
 Module:  AS_PER_fragDB
 Description:
     This module defines the interface and implementation of the Assembler
 Fragment Store, as implemented in an Oracle DB.

 Assumptions:
      Oracle DB
 Document:
      FragmentStore.rtf

 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_ReadStruct.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_fragStore.h"
#include "AS_PER_fragStore_private.h"
#include "AS_PER_SafeIO.h"
#include "AS_PER_fragDB.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_Var.h"
#include "PrimitiveVA.h"

#define AS_ORACLEx

#ifdef AS_ORACLE
#include <ocidfn.h>
#include <oratypes.h>
#include <ociap.h>
#include <oci.h>
#endif

const char *login_id = (char *) "drsistg";
const char *passwd = (char *) "welcome";

void printError(void);

#ifdef AS_ORACLE
void initializeOCI(int desiredNumRowsToRead);

void GetFragFromDB(int32 index,               /* The readIndex of the fragment we want to retrieve */
				   ShortFragRecord *fixed,    /* Address of structure to hold fixed length info */
				   void *blob,        /* Address of the  blob */
				   int *blobLength);   /* the size of the blob in bytes */


static OCIEnv * mEnvhp;
static OCISvcCtx * mSvchp;
static OCIServer * mSrvhp;
static OCIError * mErrhp;
static OCISession * mAuthp;
static OCISession * mUsrhp;
#endif

#define SEQUENCE_MASK 0xC0
#define QUALITY_MASK  (~ SEQUENCE_MASK)

#define GET_SEQUENCE(ch) (((ch) & SEQUENCE_MASK))
#define GET_QUALITY(ch) ((ch) & QUALITY_MASK)

#define SEQ_A 00
#define SEQ_C 64
#define SEQ_T 128
#define SEQ_G 192
#define SEQ_N 4
#define QUALITY_MAX 60
#define QUALITY_N 63
#define QUALITY_99 62
#define EOS '\0'        

static const char SeqChars[] = {
  'A',
  'C',
  'T',
  'G',
  'N'
};

#ifdef AS_ORACLE

#define LOBBING 1
#define MAX_READ_ROWS 10000

float deleted[MAX_READ_ROWS], 
  readType[MAX_READ_ROWS], 
  hasQuality[MAX_READ_ROWS],
  numScreenMatches[MAX_READ_ROWS],
  spare1[MAX_READ_ROWS],
  clearRegionStart[MAX_READ_ROWS],
  clearRegionEnd[MAX_READ_ROWS];
double accID[MAX_READ_ROWS];
float readIndex[MAX_READ_ROWS],

sb2 indicator[MAX_READ_ROWS];
ub2 rlenp[MAX_READ_ROWS];
ub2 rcodep[MAX_READ_ROWS];

OCIStmt *mStmtp;
OCILobLocator *mLobSequence[MAX_READ_ROWS];
int startIndex = 0, endIndex = 0;
int numRowsToRead;

#endif

#if 0
int decodeSequenceQuality
( char *encoded, int encodedLength, char *sequence, char *quality, 
  uint hasQuality){
  char *s = sequence;
  char *q = quality;
  char *e = encoded;   /*** NOTE -- encoded value is NOT a null terminated string!!!! */

  int i;

  *q = '\0';

  for(i = 0; i < encodedLength; i++){
      *s = SeqChars[ GET_SEQUENCE(*e)/ SEQ_C ]; // 0-3
    if(hasQuality){
      *q = GET_QUALITY((*e));
      if(*q == QUALITY_N){
	*s = SeqChars[SEQ_N];
	*q = 0;
      }
      *q += '0';
      q++;
    }

    s++;
    e++;
  }
  if(hasQuality)
    *q = EOS;

  *s = EOS;

  assert(strlen(sequence) == encodedLength);
#ifdef DEBUG
  fprintf(stderr,"decode seq = %s\nqual = %s\n",
	  sequence, (hasQuality?quality:""));
#endif
  return 0;
}
#endif
/***********************************************************************************
 * Function unloadFragRecord
 *    Utility routine for unloading a FragRecord.  used by getFragStore and nextFragStream.
 ***********************************************************************************/
void unloadFragBlobs(char *blob, int32 blobLength, 
		     FragRecord *fr, int32 getFlags){
  VLSTRING_SIZE_T sequenceLength, sourceLength;
  uint8 checkSum, recordChecksum;
  char encodeBuffer[MAX_SEQUENCE_LENGTH];
  uint16 localeLength=0, screenLength;
  char *cursor = blob;

  sourceLength = *(VLSTRING_SIZE_T *)cursor;
  cursor += sizeof(VLSTRING_SIZE_T);
  checkSum = *(char *)cursor;
  cursor++;

  switch(fr->frag.readType){
  case AS_READ:
  case AS_B_READ:
  case AS_EXTR:
  case AS_TRNR:
    localeLength = 0;
    break;
  case AS_STS:
  case AS_EBAC:
    localeLength = 1 + sizeof(uint64);
    break;
  case AS_UBAC:
  case AS_FBAC:
    localeLength = 1 + sizeof(uint64) + 2 * sizeof(uint32);
    break;
  default:
    assert(0);
  }

 screenLength = fr->frag.numScreenMatches * sizeof(IntScreenMatch);

#ifdef DEBUG
 fprintf(stderr,"* screenLength = %u localeLength = %u\n",
	 screenLength, localeLength);
#endif
  if((getFlags & FRAG_S_SOURCE) ||
     fr->frag.numScreenMatches > 0 || 
     localeLength != 0){

//    fprintf(stderr,"* Reading source field...screen matches %d sourceLength %u sm length %d \n",
//	    fr->frag.numScreenMatches, sourceLength, fr->frag.numScreenMatches * sizeof(IntScreenMatch));

    /******* Move the sourceBlob into the FragRecord *******/
    memcpy(fr->source, cursor, sourceLength);
    recordChecksum = checkSumBlob(cursor, sourceLength);
	
    assert(recordChecksum == checkSum);

    cursor += sourceLength;

    if(localeLength > 0){
      size_t offset = sourceLength - localeLength + 1;
//      fprintf(stderr,"* Reading Locale at offset " F_SIZE_T "\n",    offset);
      memcpy(&fr->localeID, fr->source + offset, sizeof(uint64));
      offset += sizeof(uint64);
      if(localeLength > sizeof(uint64) + 1){
	memcpy(&fr->localePosStart, fr->source + offset, sizeof(uint32));
	offset += sizeof(int32);
	memcpy(&fr->localePosEnd, fr->source + offset, sizeof(uint32));
      }
    }
      if(screenLength > 0){
	size_t offset = sourceLength - screenLength - localeLength + 1;
        /*
        fprintf(stderr,"* Reading screenMatches at offset " F_SIZE_T "\n",
                offset);
        */
	memcpy(fr->matches, fr->source + offset, fr->frag.numScreenMatches * sizeof(IntScreenMatch));
	{ // What the hell, link the suckers together
	  int i;
	  for(i = 0; i < fr->frag.numScreenMatches - 1; i++){
	    //fr->matches[i].next = fr->matches + i + 1;
	    fr->matches[i].next = &(fr->matches[i + 1]);
	  }
	  fr->matches[i].next = NULL;
	}
      }
  }
  if(getFlags & FRAG_S_SEQUENCE)
  {
	memcpy(&sequenceLength, cursor, sizeof(VLSTRING_SIZE_T));
	
	// sequenceLength = *(VLSTRING_SIZE_T *)cursor;
	cursor += sizeof(VLSTRING_SIZE_T);
	checkSum = *(char *)cursor;
	cursor++;
    /**** Copy the sequence into the encode buffer ****/
    memcpy(encodeBuffer, cursor, sequenceLength);
    recordChecksum = checkSumBlob(cursor, sequenceLength);
    assert(recordChecksum == checkSum);
	
    encodeBuffer[sequenceLength] = '\0';
    decodeSequenceQuality(encodeBuffer, sequenceLength, fr->sequence, fr->quality, fr->frag.hasQuality);
  }
}

/***********************************************************************/
/* getFragDB
   Random Access Read.
   The data is read into the previously allocated ReadStruct rs.
   The types of data read can be restricted using the getFlags.
*/
	
int getFragDB(int64 index, int32 getFlags, ReadStructp rs, int blockSize){
#ifdef AS_ORACLE
  FragRecord *fr = (FragRecord *)rs;
  int32 sequenceLength, sourceLength;
  int32 blobLength;
  char blob[16000];
  static int OCIInitialized = 0;
  
  // fprintf(stderr,"* getFragDB " F_S64 "\n", index);
  if (!OCIInitialized)
  {
	initializeOCI(blockSize);
	OCIInitialized = 1;
  }
  GetFragFromDB(index, &fr->frag, blob, &blobLength);
  unloadFragBlobs( blob, blobLength, fr, getFlags);

#ifdef DEBUG
 fprintf(stderr,"* GetFragDB " F_S64 " with id " F_IID " src: %s  seq: %s \n qu: %s\n",
		 index, fr->frag.readIndex, fr->source, fr->sequence, fr->quality);
#endif
#endif
 return(0);
}


#ifdef AS_ORACLE
void initializeOCI(int desiredNumRowsToRead)
{
  const char *cstring = (const char *)"dev_drsi";
  OCIDefine *defn1p = (OCIDefine *) 0, *defn2p = (OCIDefine *) 0, 
	*defn3p = (OCIDefine *) 0, *defn4p = (OCIDefine *) 0,
	*defn5p = (OCIDefine *) 0, *defn6p = (OCIDefine *) 0,
	*defn7p = (OCIDefine *) 0, *defn8p = (OCIDefine *) 0,
	*defn9p = (OCIDefine *) 0, *defn10p = (OCIDefine *) 0,
	*defn11p = (OCIDefine *) 0, *defn12p = (OCIDefine *) 0;
  char strSel[200];
  OCIBind  *bnd1p = NULL;
  char tableName[200] = "short_fragment_t1";

  if (desiredNumRowsToRead > 0 && desiredNumRowsToRead <= MAX_READ_ROWS)
	numRowsToRead = desiredNumRowsToRead;
  else
	numRowsToRead = MAX_READ_ROWS;

  fprintf(stderr, "numRowsToRead = %d\n", numRowsToRead);
  
  if (OCIInitialize(0, (dvoid *)0,(dvoid * (*)(dvoid *, size_t)) 0,
					(dvoid * (*)(dvoid *, dvoid *, size_t))0,
					(void (*)(dvoid *, dvoid *)) 0 ))
  {
	/* fprintf (stderr, "FAILED TO INITIALIZE OCI"); */
  }

  /* Inititialize the OCI Environment */
  if (OCIEnvInit(&mEnvhp, (ub4) OCI_DEFAULT,
				 (size_t) 0, (dvoid **) 0 ))
  {
	fprintf (stderr, "FAILED TO INITIALIZE ENVIRONMENT");
  }
  /* Allocate a service handle */
  if (OCIHandleAlloc(mEnvhp, (dvoid**) &mSvchp,
					 (ub4) OCI_HTYPE_SVCCTX, (size_t) 0, (dvoid **) 0))
  {
	fprintf (stderr, "FAILED IN OCIHandleAlloc on svchp");
  }
  /* Allocate an error handle */
  if (OCIHandleAlloc(mEnvhp, (dvoid**) &mErrhp,
					 (ub4) OCI_HTYPE_ERROR, (size_t) 0, (dvoid **) 0))
  {
	fprintf (stderr, "FAILED: OCIHandleAlloc() on errhp\n");
  }
  /* Allocate a server handle */
  if (OCIHandleAlloc(mEnvhp, (dvoid **) &mSrvhp,
					 (ub4) OCI_HTYPE_SERVER, (size_t) 0, (dvoid **) 0))
  {
	fprintf (stderr, "FAILED: OCIHandleAlloc() on srvhp\n");
  }
  if (OCIHandleAlloc(mEnvhp, (dvoid **) &mAuthp,
					 (ub4) OCI_HTYPE_SESSION, (size_t) 0, (dvoid **) 0))
  {
	fprintf (stderr, "FAILED: OCIHandleAlloc() on authp\n");
  }
  /* ATTACHES TO THE DEFAULT SERVER */
  if (OCIServerAttach(mSrvhp, mErrhp, (text *) cstring,
					  (sb4) strlen((char *)cstring), (ub4) OCI_DEFAULT))
  {
	fprintf (stderr, "FAILED: OCIServerAttach()\n");
  }
  /* Set the server handle in the service handle */
  if (OCIAttrSet((dvoid *) mSvchp, (ub4) OCI_HTYPE_SVCCTX,
				 (dvoid *) mSrvhp, (ub4) 0, (ub4) OCI_ATTR_SERVER, mErrhp))
  {
	fprintf (stderr, "FAILED: OCIAttrSet() server attribute\n");
  }
  if (OCIAttrSet((dvoid *) mAuthp, (ub4) OCI_HTYPE_SESSION,
				 (dvoid *) login_id, (ub4) strlen((char *) login_id),
				 (ub4) OCI_ATTR_USERNAME, mErrhp))
  {
	fprintf (stderr, "FAILED: OCIAttrSet() userid\n");
  }
  if (OCIAttrSet((dvoid *) mAuthp, (ub4) OCI_HTYPE_SESSION,
				 (dvoid *) passwd, (ub4) strlen((char *) passwd),
				 (ub4) OCI_ATTR_PASSWORD, mErrhp))
  {
	fprintf (stderr, "FAILED: OCIAttrSet() passwd\n");
  }

  if (OCISessionBegin(mSvchp, mErrhp, mAuthp, (ub4) OCI_CRED_RDBMS, 
					  (ub4) OCI_DEFAULT))
  {
	fprintf (stderr, "FAILED: OCIAttrSet() passwd\n");
  }

  fprintf(stderr, "logged on.\n");

  /* Set the authentication handle in the Service handle */
  if (OCIAttrSet((dvoid *) mSvchp, (ub4) OCI_HTYPE_SVCCTX,
				 (dvoid *) mAuthp, (ub4) 0, (ub4) OCI_ATTR_SESSION, mErrhp))
  {
	fprintf (stderr, "FAILED: OCIAttrSet() session\n");
  }



/* statement prep */
  {
	FILE *in;
	
	in = fopen("tableName.txt", "r");
	if (in == NULL)
	  fprintf( stderr, "Could not open tableName.txt\n");
	else
	  fgets( tableName, 200, in);
	fclose(in);
  }
  fprintf( stderr, "*** Using table: %s\n", tableName);
	  

  // note: selecting deleted twice because spare1 did not make it into the database
  if (numRowsToRead > 1)
	sprintf(strSel, "Select deleted, readType, hasQuality, numScreenMatches, deleted, clearRegionStart, clearRegionEnd, accID, readIndex, time_t, seq from %s where readindex >= :1 and readindex <= :2", tableName);
  else
	sprintf(strSel, "Select deleted, readType, hasQuality, numScreenMatches, deleted, clearRegionStart, clearRegionEnd, accID, readIndex, time_t, seq from %s where readindex = :1", tableName);

  fprintf(stderr, "Query: %s\n", strSel);
  
  if (OCIHandleAlloc((dvoid*)mEnvhp, (dvoid**) &mStmtp, 
					 (ub4)OCI_HTYPE_STMT, (CONST size_t) 0, (dvoid **) 0))
  {
	fprintf(stderr, "Failed to alloc statement handle");
  }
  
  if (OCIStmtPrepare(mStmtp, mErrhp, (unsigned char*)strSel, 
					 (ub4)strlen(strSel), OCI_NTV_SYNTAX,
					 OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to prepare statement ");
  }

  if (OCIBindByPos (mStmtp, 
					&bnd1p, 
					mErrhp, 
					1,
					(dvoid *) &startIndex, 
					(sb4) sizeof(startIndex), 
					SQLT_INT,
					(dvoid *) 0, 
					(ub2 *) 0, 
					(ub2) 0, 
					(ub4) 0, 
					(ub4 *) 0, 
					OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to bind index\n");
	printError();
  }

  if (numRowsToRead > 1)
	if (OCIBindByPos (mStmtp, 
					  &bnd1p, 
					  mErrhp, 
					  2,
					  (dvoid *) &endIndex, 
					  (sb4) sizeof(endIndex), 
					  SQLT_INT,
					  (dvoid *) 0, 
					  (ub2 *) 0, 
					  (ub2) 0, 
					  (ub4) 0, 
					  (ub4 *) 0, 
					  OCI_DEFAULT))
	{
	  fprintf(stderr, "Failed to bind index\n");
	  printError();
	}
  
/* bindings */
  if ( OCIDefineByPos(mStmtp, 
					  &defn1p, 
					  mErrhp, 
					  1, 
					  (dvoid*) deleted, 
					  sizeof(deleted[0]), 
					  SQLT_FLT, 
					  (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep, 
					  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 1");
  }

  if ( OCIDefineByPos(mStmtp, &defn2p, mErrhp, 2, 
					  (dvoid*) readType, sizeof(readType[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 2");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn3p, mErrhp, 3, 
					  (dvoid*) hasQuality, sizeof(hasQuality[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 3");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn4p, mErrhp, 4, 
					  (dvoid*) numScreenMatches, sizeof(numScreenMatches[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 4");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn5p, mErrhp, 5, 
					  (dvoid*) spare1, sizeof(spare1[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 5");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn6p, mErrhp, 6, 
					  (dvoid*) clearRegionStart, sizeof(clearRegionStart[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 6");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn7p, mErrhp, 7, 
					  (dvoid*) clearRegionEnd, sizeof(clearRegionEnd[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 7");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn8p, mErrhp, 8, 
					  (dvoid*) accID, sizeof(accID[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 8");
  }
  
  if ( OCIDefineByPos(mStmtp, &defn9p, mErrhp, 9, 
					  (dvoid*) readIndex, sizeof(readIndex[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 9");
  }  

  if ( OCIDefineByPos(mStmtp, &defn10p, mErrhp, 10, 
					  (dvoid*) entryTime, sizeof(entryTime[0]), 
					  SQLT_FLT, (void*) indicator, 
					  (ub2*) &rlenp, 
					  (ub2*) rcodep,  OCI_DEFAULT))
  {
	fprintf(stderr, "Failed to define by pos 10");
  }  

/* now do LOBs */
  if (LOBBING)
  {
	int i;

	// initialize all we could ever want to use
	for (i = 0; i < MAX_READ_ROWS; i++)
	  if ( OCIDescriptorAlloc(mEnvhp, (dvoid**) &mLobSequence[i],
							  (ub4) OCI_DTYPE_LOB, (size_t) 0, 
							  (dvoid **) 0)	!= OCI_SUCCESS )
	  {
		fprintf(stderr, "Failed to allocate LOB sequence descriptor");
	  }
	
	if ( OCIDefineByPos(mStmtp, &defn11p, mErrhp, 11, (dvoid*) mLobSequence,
						(sb4)-1, (ub2)SQLT_BLOB, 0, 0, 0, OCI_DEFAULT))
	{
	  fprintf(stderr, "Define failed for lob 11\n");
	  printError();
	}
  }
  // fprintf(stderr, "\n\n");
}

// all the initialization stuff needs to move outta here
// also, we are not cleaning up after ourselves memory-wise

void GetFragFromDB(int32 index,                 /* The readIndex of the fragment we want to retrieve */
				   ShortFragRecord *fixed,      /* Address of structure to hold fixed length info */
				   void *sequenceBlob,          /* Address of the sequence blob */
				   int *sequenceBlobLength)     /* the size of the sequenceBlob in bytes */
{
  sword empno, sal, deptno;
  sword len, len2, rv, dsize, dsize2;
  sb4   enamelen, joblen, deptlen;
  sb2   sal_ind, job_ind;
  sb2   db_type, db2_type;
  sb1   name_buf[20], name2_buf[20];
  text  *cp, *ename, *job, *dept;
  char col_name[256], data_type[256];
  sb4 status;
  ub4 amtp = 5000000; // approx 5mb
  int indexInArray;

/*  fprintf(stderr, "requested index = %d\n", index);
	fprintf(stderr, "currently startIndex = %d and endIndex = %d\n", startIndex, endIndex); */

// we need to determine how to set endIndex if we don't get back a full set of reads

  if (index < startIndex || index > endIndex)
  {
	startIndex = index;
	endIndex = startIndex + numRowsToRead - 1;

	// fprintf(stderr, "at execute, startIndex = %d and endIndex = %d\n", startIndex, endIndex);
	
/* execution */
	if ( (status = OCIStmtExecute(mSvchp, 
								  mStmtp, 
								  mErrhp, 
								  numRowsToRead, 
								  0, 
								  0, 
								  0, 
								  OCI_EXACT_FETCH))
		 //                       OCI_DEFAULT))
		 && (status != OCI_NO_DATA) )
	{
	  fprintf(stderr, "Failed to Execute Statement");
	  printError();
	}
  }
  
  // now fetch results - for now everything is returned in the call above
  /*
	status = OCIStmtFetch(mStmtp, 
	mErrhp, 
	(ub4) READ_ROWS, 
	(ub4) OCI_FETCH_NEXT, 
	OCI_DEFAULT);
	
	if (status == OCI_SUCCESS || status == OCI_SUCCESS_WITH_INFO)
  */
  
  // indexInArray is where the row corresponding to the current index lives in the arrays we have
  indexInArray = index - startIndex;

  fixed->accID = (long int) accID[indexInArray];
  fixed->deleted = (int) deleted[indexInArray];
  fixed->readType = (char) readType[indexInArray];
  fixed->hasQuality = (int) hasQuality[indexInArray];
  fixed->numScreenMatches = (int) numScreenMatches[indexInArray];
  fixed->spare1 = 0;  // not in db, get rid of?
  fixed->clearRegionStart = (int) clearRegionStart[indexInArray];
  fixed->clearRegionEnd = (int) clearRegionEnd[indexInArray];
  fixed->readIndex = (int) readIndex[indexInArray];
  
  if (LOBBING)
  {
	int lobStatus;
	//for (i = 0; i < numRowsToRead; i++)
	{
	  if ( OCILobGetLength(mSvchp, mErrhp, mLobSequence[indexInArray], (uint *) sequenceBlobLength) )
	  {
		fprintf(stderr, "Failed to read sequence length");
		printError();
	  }
//	  fprintf(stderr, "indexInArray: %d, sequence length: %d\n", indexInArray, *sequenceBlobLength);
	  
	  // fprintf(stderr, "about to do OCILobRead\n");
	  
	  lobStatus = OCILobRead( mSvchp, mErrhp, mLobSequence[indexInArray], &amtp, 1,
							  (dvoid*) sequenceBlob,
							  (ub4) *sequenceBlobLength, (dvoid*)0,
							  (sb4 (*)(dvoid *, CONST dvoid *, ub4, ub1 ))0,
							  (ub2) 0, (ub1) SQLCS_IMPLICIT );
	  if (lobStatus != OCI_SUCCESS)
	  {
		fprintf(stderr, "indexInArray: %d, lob read failed!\n", indexInArray);
		printError();
		assert(0);
	  }
	}
  }
  else
  {
	fprintf(stderr, "status: %d, OCI_SUCCESS: %d\n", status, OCI_SUCCESS);
	printError();
  }
}

void printError(void)
{
  text  msgbuf[512];
  sb4   errcode = 0;

  OCIErrorGet((dvoid *) mErrhp, (ub4) 1, (text *) NULL, &errcode,
			  msgbuf, (ub4) sizeof(msgbuf), (ub4) OCI_HTYPE_ERROR);
  fprintf( stderr, "ERROR CODE = %d\n", errcode);
  fprintf( stderr, "ERROR TEXT = %s\n", msgbuf);
  
  return;
}

#endif
