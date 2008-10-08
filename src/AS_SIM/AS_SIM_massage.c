
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

static const char *rcsid = "$Id: AS_SIM_massage.c,v 1.5 2008-10-08 22:03:00 brianwalenz Exp $";

/*
   Module:      FASTA to 3-code CA Messages Converter.
   Description:
       Massage is a small filter (i.e. reads stdin and produces stdout).

       The input is expected to be the output of a call to "frag" with
       the -F option set.  That is, the input is a FASTA formatted file
       of sequences, where the header file has one of the 3 forms --
       #, #f, or #r -- indicating that the next sequence is either a single
       read, a forward read of a mated pair, or the reverse read of a mated
       pair, respectively.  Moreover, any forward read of pair is always
       immediately followed by the reverse read of the pair.  One can also
       safely assume that every input line is less than 80chars (but the
       code checks and protects just in case).

       The output is a sequence of 3-code Celera Assembler messages for the
       fragments.  See the document on Celera Assembler prototype conventions
       for the exact format of such records.

       The filter expects 2 integer command-line arguments.  The first is
       the accession number of the distance message describing the spatial
       relationships between pairs, and the second is the accession number
       at which to start sequentially allocating such numbers to fragemnts.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h> /* man 3 getopt */
#include <assert.h>
#include <math.h>

#include "AS_global.h"
#include "AS_UTL_rand.h"

double drand48();
void   srand48();

#undef  DEBUG           /* Causes extra diagnostic stuff to get printed */
//#define DEBUG
#define INPUTMAX  2048  /* Maximum poly spec. line length (checked) */
#define NOT_PRESENT 0L

static char *ProgName;  /* Program name 'poly' */

static int Uniform = FALSE;
static int BacEnds = 0;
static int UBacs = 0;
static int LBacs = 0;
static int shreddedFBacs = 0;
static int Bactigs = 0;
static int fullBacs = 0;
static CDS_UID_t   InitialFragId = 0; /* Initial frag accession # */
static CDS_UID_t   FragId = 0;    /* Current frag accession # */
static CDS_UID_t   DistId = 0;    /* Distance message accession # */
static int64 LBacOffset = 0;
static int MinPrefixUnclear = 0;
static int MaxPrefixUnclear = 0;
static int MinSuffixUnclear = 0;
static int MaxSuffixUnclear = 0;

static FragMesg OutFrgMesg;
static GenericMesg OutMesg, OutLnkMesg;
static LinkMesg LnkMesg;
int FalseMate = 0;

#define MAX_SEQ_LENGTH 2048*2048
#define MAX_SOURCE_LENGTH 2048*2048
static char Source[MAX_SOURCE_LENGTH] = {0};
static char Sequence[MAX_SEQ_LENGTH] = {0};
static char Quality[MAX_SEQ_LENGTH] = {0};
static FILE *bacfile = NULL;
static FILE *qltfile = NULL;
int64 BegInterval = 0, EndInterval = 0;
int BTrim = 0, ETrim = 0;
int SequenceLength = 0;
CDS_UID_t Locale = 0;
int Reversed = 0;
int numBactigs = 0;
InternalBactigMesg *bacLengths = NULL;

static time_t tp = CREATION_TIME; /* int32 seconds from sometime in 1970 */

void process(void);

/* ****************************************************** */
void usage(char *ProgName)
{ fprintf(stderr,"\n\t*** Usage is: %s [-q qltfile.fasta] [-UNrelufsbc] "
		  " [PrefixMinUnclear PrefixMaxUnclear SuffixMinUnclear SuffixMaxUnclear]\n",
		  ProgName);
 exit(EXIT_FAILURE);
}

/* ****************************************************** */


/* Output the tail of a fragment record */

#define SEGLEN 70
char dummyQuality[2 * (SEGLEN + 1)];
char BACQuality[2 * (SEGLEN + 1)];

void generateQuality(int FragLen, char *quality)
{ register int i;
 int offset;

 offset = 0;
 for (i = 0; i < FragLen; i += SEGLEN){
   int rand = GetRand_AS(0,SEGLEN - 1, TRUE);
   sprintf(quality + offset,"%.*s",SEGLEN,dummyQuality+ rand);
   offset += SEGLEN;
 }
}

void copyBACQuality(int FragLen, char *quality){
  int remainder = FragLen%SEGLEN;
  int offset = 0;
  int i;
  for(i = 0; i < FragLen/SEGLEN; i++){
    memcpy(quality + offset, BACQuality, SEGLEN);
    offset += SEGLEN;
    //    strcpy(quality + offset,"\n");
    //    offset++;
  }
  if(remainder){
    strncpy(quality + offset, BACQuality, remainder);
    offset += remainder;
    //    strcpy(quality + offset,"\n");
  }
  quality[offset] = '\0';
}

void copyQuality(int FragLen, char *quality){
  int remainder = FragLen%SEGLEN;
  int offset = 0;
  int i;
  for(i = 0; i < FragLen/SEGLEN; i++){
    int rand = GetRand_AS(0,SEGLEN - 1, TRUE);
    memcpy(quality + offset, dummyQuality + rand, SEGLEN);
    offset += SEGLEN;
    //    strcpy(quality + offset,"\n");
    //    offset++;
  }
  if(remainder){
    int rand = GetRand_AS(0,SEGLEN - 1, TRUE);

    strncpy(quality + offset,dummyQuality + rand, remainder);
    offset += remainder;
    //    strcpy(quality + offset,"\n");
  }
  quality[offset] = '\0';
}

void getFastaQuality(FILE *fin, int length, char *quality)
{
  char *line = quality;
  int c;
  int got=0;
  c = fgetc(fin), *line = c;
  if (*line == '>' || *line == '\n') { // ignore the comment line of the fasta file.
	while ( c = fgetc(fin), *line = c, c != '\n') {}
  } else {
	ungetc(c,fin);
  }
  while (got<length) {
    c = toupper(fgetc(fin));
    *line = c;
    if(*line == ' ') {
      fprintf(stderr,"Got %d characters from file (expecting %d) \n",got,length);
	  assert(*line != ' ');
    }
    if (c == EOF) break;
    if (*line != '\n') {
      got++;
      line++;
    }
  }
  *line = '\0';
}

int getBacLengths(FILE *bacin)
{
  int icnt, cntBactigs;
  char buffer[256];

  fgets(buffer, 256, bacin);
  cntBactigs = atoi(buffer);

  // fprintf( stderr, "cntBactigs = %d\n", cntBactigs);

  bacLengths = (InternalBactigMesg *) malloc (cntBactigs * sizeof(InternalBactigMesg));
  if (bacLengths == NULL)
  {
	fprintf( stderr, "Failed to malloc space for bacLengths\n");
	exit(1);
  }

  for (icnt = 0; icnt < cntBactigs; icnt++)
  {
	CDS_UID_t accession;
	int length, numMatched;

	fgets(buffer, 256, bacin);
	numMatched = sscanf(buffer, F_UID " %d", &accession, &length);
	assert(numMatched == 2);
	bacLengths[icnt].eaccession = accession;
	bacLengths[icnt].length = length;
  }
  return cntBactigs;
}

void mate_link(CDS_UID_t frg1, CDS_UID_t frg2){
  LnkMesg.frag1 = frg1;
  LnkMesg.frag2 = frg2;
  LnkMesg.distance = DistId;
  LnkMesg.entry_time = tp;
  /*   if(FalseMate){
	   fprintf(stdout,"# Following Mate link is FALSE\n");
	   }*/
  WriteProtoMesg_AS(stdout,&OutLnkMesg);

}


void outputBAC(int type, CDS_UID_t distacnum, CDS_UID_t locacnum, CDS_UID_t sequence_id, int64 BacGenomeBegin, int64 BacGenomeEnd)
{

  GenericMesg outMesg;
  BacMesg bacMesg;
  char buffer[200];

  bacMesg.action = AS_ADD;
  bacMesg.elength = distacnum;
  bacMesg.ebac_id = locacnum;
  bacMesg.eseq_id = sequence_id;
  bacMesg.iseq_id = 0;
  bacMesg.type = type;
  bacMesg.entry_time = tp; // reproducible timestamp for regression tests
  bacMesg.num_bactigs = numBactigs;
  bacMesg.bactig_list = bacLengths;
  bacMesg.ilength = 0;

  if (type != AS_EBAC)
  {
	sprintf(buffer, "BAC genome coords: [" F_S64 "," F_S64
                "]\n", BacGenomeBegin, BacGenomeEnd);
	bacMesg.source = buffer;
  }
  else
  {
	sprintf(buffer, "BAC ends\n");
	bacMesg.source = buffer;
  }

  outMesg.m = &bacMesg;
  outMesg.t = MESG_BAC;

  // fprintf( stderr, "in outputBAC, about to writeProtoMesg_AS\n");

  WriteProtoMesg_AS(stdout,&outMesg);

  // fprintf( stderr, "in outputBAC, after writeProtoMesg_AS\n");

  fflush(NULL);
}

void output(int direct)
{
  if(direct == -2)
	return;

  fflush(stderr);
  assert(SequenceLength > 0);

  OutFrgMesg.clear_rng.bgn = BTrim;
  OutFrgMesg.clear_rng.end = SequenceLength - ETrim;

  assert( OutFrgMesg.source != NULL );
  assert( OutFrgMesg.clear_rng.bgn >= 0);
  assert( OutFrgMesg.clear_rng.end <= SequenceLength);

#ifdef DEBUG
  fprintf(stderr,"* Output SequenceLength:%d clr:" F_COORD "," F_COORD
          " clrL:%d gen:(" F_S64 "," F_S64 ") genL:" F_S64 "\n\t%s",
		  SequenceLength,
		  OutFrgMesg.clear_rng.bgn, OutFrgMesg.clear_rng.end,
		  abs(OutFrgMesg.clear_rng.bgn - OutFrgMesg.clear_rng.end),
		  BegInterval, EndInterval,
		  abs(EndInterval - BegInterval),
		  OutFrgMesg.source);
#endif

  if (BacEnds)
	OutFrgMesg.type = AS_EBAC;
  else if (LBacs)
	OutFrgMesg.type = AS_LBAC;
  else if (UBacs)
	OutFrgMesg.type = AS_UBAC;
  else if (shreddedFBacs)
	OutFrgMesg.type = AS_FBAC;
  else if (Bactigs)
	OutFrgMesg.type = AS_BACTIG;
  else if (fullBacs)
	OutFrgMesg.type = AS_FULLBAC;
  else
	OutFrgMesg.type = AS_READ;

  if (qltfile == NULL ) {
	if (LBacs || UBacs || shreddedFBacs || Bactigs || fullBacs)
	  copyBACQuality(SequenceLength, OutFrgMesg.quality);
	else
	  copyQuality(SequenceLength, OutFrgMesg.quality);
  } else {
	getFastaQuality(qltfile, SequenceLength, OutFrgMesg.quality);
	assert(strlen(OutFrgMesg.sequence) == strlen(OutFrgMesg.quality));
  }

  if(SequenceLength < 50)
	fprintf(stderr," Fragment length is wierd: %d\n",
			SequenceLength);

  {
	int i;
	char *c;
	for(i = 0, c = Sequence; i < SequenceLength;c++){
	  if(*c != '\n')
		i++;
	  if(!isalpha(*c)){
		fprintf(stderr,"* Non Alpha character 0x%x at pos %d \n",
				*c, i);
	  }
	}
  }
  WriteProtoMesg_AS(stdout,&OutMesg);

  if(direct == -1){
	mate_link(FragId -1 , FragId);
  }
  FragId += 1;
}



/* Read stdin line at a time, transforming to uniform
   assembler message format. */


void process(void)
{
  int  forward = 0, length = 0;
  //char *end = NULL, *start = NULL;
  char *start = NULL;
  char buffer[INPUTMAX] = {0};
  int didOutputBAC = 0;
  int64 BacGenomeBegin = -1, BacGenomeEnd = -1;  // where the BAC actually is in the genome

  //  FILE *test = fopen("test.out","w");
  forward = -2;  /* State variable: no record in progress */
  FalseMate = 0;
  buffer[INPUTMAX-2] = '\n';

  while (fgets(buffer,INPUTMAX,stdin) != NULL){
	//	  fprintf(test,"%s",buffer);
	//	  fflush(test);
	if (buffer[INPUTMAX-2] != '\n')
	{ fprintf(stdout,"\n\t*** %s: Line over %d chars long\n",
			  ProgName,INPUTMAX-2);
	exit(1);
	}
	if (*buffer != '>')  /* sequence line, echo to output */
	{
	  int bufLength = strlen(buffer);


	  // New and 'improved' protoio handles this differently than its predecessor
	  // Strip newline from end of string.
	  if(length < 50 && bufLength < 10)
		fprintf(stderr,"* Wired line length %d read %s\n",
				bufLength, buffer);

	  if(*buffer != '\n'){
		assert(buffer[bufLength -1] == '\n');
		buffer[--bufLength] = '\0';
		strcat(Sequence, buffer);
		length += bufLength;
	  }
	  SequenceLength = length;
	}
	else
	{
	  char *interval = strstr(buffer + 1,"[");
#ifdef DEBUG
	  fprintf(stderr,"* found a header interval = %x\n",
			  interval);
#endif
	  output(forward);      // Output the previous fragment, if any

	  assert(NULL != interval);
	  {
	    int64 b,e;
	    int tmpB, tmpE;
	    int numMatched = sscanf(interval+1, F_S64 "," F_S64 "]",
				    &b, &e);
	    assert(numMatched ==2);

	    if (LBacs)
	      {
		b += LBacOffset;
		e += LBacOffset;
	      }

	    // The clear range is the the entire fragment length - a uniformly distributed
	    // piece at either end in the range of 0-50 base pairs.
	    // How much to trim from fragment to make clear range
	    tmpB = GetRand_AS(MinPrefixUnclear,MaxPrefixUnclear, Uniform);
	    BTrim =  MAX(0,tmpB); // Trim from frag prefix
	    tmpE= GetRand_AS(MinSuffixUnclear,MaxSuffixUnclear, Uniform); // Trim from frag suffix
	    ETrim =  MAX(0,tmpE);

	    assert(BTrim >= 0);
	    assert(ETrim >= 0);

	    assert(b < e);

	    if(strstr(interval,"REVERSED")){
	      BegInterval = e - BTrim;
	      EndInterval = b + ETrim;
	      Reversed = 1;
	    }else{
	      BegInterval = b + BTrim;
	      EndInterval = e - ETrim;
	      Reversed = 0;
	    }
	    if(strstr(interval + 1, "False")){
	      FalseMate = 1;
#ifdef DEBUG
	      fprintf(stderr,"* Detected false mate \n%s\n", buffer);
#endif
	    }
	  }

	  // look for bactig_id
	  interval = strstr(buffer + 1,"bactig_id:");
	  if (interval)
	  {
		CDS_UID_t id;
		int numMatched = sscanf( interval + strlen("bactig_id:"),
                                         F_UID, &id);
		assert(numMatched == 1);
		OutFrgMesg.ebactig_id = id;
	  }
	  else
	  {
		OutFrgMesg.locale_pos.bgn = 0;
		OutFrgMesg.locale_pos.end = 0;
	  }

	  // look for locpos
	  interval = strstr(buffer + 1,"locpos[");
	  if (interval)
	  {
		int64 b,e;
		int numMatched = sscanf( interval + strlen(                                                          "locpos["), F_S64 "," F_S64 "]", &b, &e);
		assert(numMatched == 2);
		OutFrgMesg.locale_pos.bgn = b;
		OutFrgMesg.locale_pos.end = e;
	  }
	  else
	  {
		OutFrgMesg.locale_pos.bgn = 0;
		OutFrgMesg.locale_pos.end = 0;
	  }

	  // look for genomepos
	  interval = strstr(buffer + 1,"genomepos[");
	  if (interval)
	  {
		int64 b,e;
		int numMatched = sscanf( interval + strlen("genomepos["),
                                         F_S64 "," F_S64 "]", &b, &e);
		assert(numMatched == 2);
		BacGenomeBegin = b;
		BacGenomeEnd = e;
	  }
	  else
	  {
		BacGenomeBegin = -1;
		BacGenomeEnd = -1;
	  }

	  // New protoio handles this differently than before, we don't need the newline
	  OutFrgMesg.entry_time = tp;

	  start = buffer + 1;
	  while (isspace(*start))
		start += 1;
	  {
	    char *end = NULL;
            CDS_UID_t fragID = STR_TO_UID(start,&end,10); /* Move past the name */
		char direction = end[0];
		// LOCALE_OFFSET puts us into the UID space above 10^9
		// InitialFragID gives us a unique part of this space
		fragID += InitialFragId;
		OutFrgMesg.eaccession = FragId;

		switch(direction)
		{
		  case (int)' ':
			if(BacEnds)
			{
			  OutFrgMesg.elocale = LOCALE_OFFSET + FragId;
			  outputBAC( AS_ENDS, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id, BacGenomeBegin, BacGenomeEnd);
			}
			else if (UBacs || Bactigs)
			{
			  if (!didOutputBAC)
			  {
                            /*
                            fprintf( stderr, "AS_UBAC: %d, DistId: " F_UID ", OutFrgMesg.elocale: " F_UID ", OutFrgMesg.eseq_id: " F_UID "\n",
                                     AS_UBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id);
                            */
				outputBAC( AS_UBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id, BacGenomeBegin, BacGenomeEnd);
				didOutputBAC = 1;
			  }
			}
			else if (LBacs)
			{
			  if (!didOutputBAC)
			  {
                            /*
                            fprintf( stderr, "AS_UBAC: %d, DistId: " F_UID ", OutFrgMesg.elocale: " F_UID ", OutFrgMesg.eseq_id: " F_UID "\n",
                                     AS_UBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id);
                            */
				outputBAC( AS_LBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id, BacGenomeBegin, BacGenomeEnd);
				didOutputBAC = 1;
			  }
			}
			else if (shreddedFBacs || fullBacs)
			{
			  if (!didOutputBAC)
			  {
				//fprintf( stderr, "AS_FBAC: %d, DistId: " F_UID ", OutFrgMesg.elocale: " F_UID ", OutFrgMesg.eseq_id: " F_UID "\n",
				//	AS_FBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id);
				outputBAC( AS_FBAC, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id, BacGenomeBegin, BacGenomeEnd);
				didOutputBAC = 1;
			  }
			}
			forward = 0;  /* single record in progress */
			break;
		  case (int)'f':
			if (BacEnds)
			{
			  OutFrgMesg.elocale = LOCALE_OFFSET + FragId;  // Locale = LOCALE_OFFSET + FragId;
			  outputBAC( AS_ENDS, DistId, OutFrgMesg.elocale, OutFrgMesg.eseq_id, BacGenomeBegin, BacGenomeEnd);
			}
			forward = 1;  /* forward read of pair in progress */
			break;
		  case (int)'r':
			forward = -1;  /* reverse read of pair in progress */
			break;
		  default:
		  {
			fprintf(stdout,"\n\t*** %s: Unparsable FASTA header\n\t\t%s\n",
					ProgName,buffer);
			exit(1);
		  }
		}

		sprintf(OutFrgMesg.source,
			F_UID "%c %s\n[" F_S64 "," F_S64 "]",
			fragID, direction, (FalseMate == 1?"False Mate":""),
                        BegInterval, EndInterval);
		FalseMate = 0;


		length = 0;
		*Sequence = '\0';
	  }
	}
  }
  output(forward);
  //  fclose(test);
}


/* >>>> MAIN / TOP <<<< */

int main(int argc, char *argv[]){
  char *end1, *end2;
  CDS_UID_t elocale = 0, eseq_id = 0, ebactig_id = 0;
  ProgName = argv[0];

/* Process program flags. */

  { /* Parse the argument list using "man 3 getopt". */
	int ch,errflg=0;
	optarg = NULL;

	Uniform       = 0;   // Gaussian is default
	BacEnds       = 0;   // Mate Links is default
	LBacs         = 0;   // Mate Links is default
	UBacs         = 0;   // Mate Links is default
	shreddedFBacs = 0;   // Mate Links is default
	Bactigs       = 0;   // Mate Links is default
	fullBacs      = 0;   // Mate Links is default
	bacfile = NULL;
	qltfile = NULL;
	while (!errflg && ((ch = getopt(argc, argv, "UNq:relu:f:sb:c:q")) != EOF)){
	  switch(ch) {

		case 'U':
		  Uniform = 1;
		  break;
		case 'N':
		  Uniform = 0;
		  break;

		case 'r':
		  fprintf(stderr,"**** Generating Reads ****\n");
		  DistId = STR_TO_UID(argv[optind++], &end1, 10);
		  break;
		case 'e':
		  fprintf(stderr,"**** Generating BacEnds ****\n");
		  BacEnds = 1;
		  DistId = STR_TO_UID(argv[optind++], &end1, 10);
		  break;
		case 'l':
		  fprintf(stderr,"**** Generating lightly shotgunned BACs ****\n");
		  LBacs = 1;
		  DistId = STR_TO_UID(argv[optind++], &end1, 10);
		  elocale = STR_TO_UID(argv[optind++], &end1, 10);
		  LBacOffset = STR_TO_UID(argv[optind++], &end1, 10);
		  fprintf( stderr, "LBacOffset = " F_S64 "\n",
                           LBacOffset);
		  exit(1);
		  break;
		case 'u':
		  fprintf(stderr,"**** Generating shredded UBACs ****\n");
		  UBacs = 1;
		  DistId =  STR_TO_UID(argv[optind++], &end1, 10);
		  elocale = STR_TO_UID(argv[optind++], &end1, 10);
		  eseq_id = STR_TO_UID(argv[optind++], &end1, 10);
		  fprintf(stderr,"**** Reading bactig lengths from file ****\n");
		  if( (bacfile = fopen(optarg, "r")) == NULL ) {
			fprintf(stderr, "Bac lengths file %s could not be opened\n", optarg);
            usage(ProgName);
		  }
		  numBactigs = getBacLengths(bacfile);
		  fclose(bacfile);
		  break;
		case 'f':
		  fprintf(stderr,"**** Generating shredded FBACs ****\n");
		  shreddedFBacs = 1;
		  DistId =  STR_TO_UID(argv[optind++], &end1, 10);
		  elocale = STR_TO_UID(argv[optind++], &end1, 10);
		  eseq_id = STR_TO_UID(argv[optind++], &end1, 10);
		  fprintf(stderr,"**** Reading bac length from file ****\n");
		  if( (bacfile = fopen(optarg, "r")) == NULL ) {
			fprintf(stderr, "Bac lengths file %s could not be opened\n", optarg);
            usage(ProgName);
		  }
		  numBactigs = getBacLengths(bacfile);
		  fclose(bacfile);
		  break;
		case 'b':
		  fprintf(stderr,"**** Generating BACtigs ****\n");
		  Bactigs = 1;
		  DistId =  STR_TO_UID(argv[optind++], &end1, 10);
		  elocale = STR_TO_UID(argv[optind++], &end1, 10);
		  eseq_id = STR_TO_UID(argv[optind++], &end1, 10);
		  fprintf(stderr,"**** Reading bactig lengths from file ****\n");
		  if( (bacfile = fopen(optarg, "r")) == NULL ) {
			fprintf(stderr, "Bac lengths file %s could not be opened\n", optarg);
            usage(ProgName);
		  }
		  numBactigs = getBacLengths(bacfile);
		  fclose(bacfile);
		  break;
		case 'c':
		  fprintf(stderr,"**** Generating full BACs ****\n");
		  fullBacs = 1;
		  DistId =  STR_TO_UID(argv[optind++], &end1, 10);
		  elocale = STR_TO_UID(argv[optind++], &end1, 10);
		  eseq_id = STR_TO_UID(argv[optind++], &end1, 10);
		  fprintf(stderr,"**** Reading bactig lengths from file ****\n");
		  if( (bacfile = fopen(optarg, "r")) == NULL ) {
			fprintf(stderr, "Bac lengths file %s could not be opened\n", optarg);
            usage(ProgName);
		  }
		  numBactigs = getBacLengths(bacfile);
		  fclose(bacfile);
		  break;

		case 'q':
		  fprintf(stderr,"**** Reading quality values from file ****\n");
		  if( (qltfile = fopen(optarg,"r")) == NULL ) {
			fprintf(stderr,"Input quality file %s could not be opened\n",optarg);
            usage(ProgName);
		  }
		  break;
		case 'P':
		  fprintf(stderr, "-P is depricated; the default now.\n");
		  break;
		case '?':
		  fprintf(stderr,"Unrecognized option -%c",optopt);
		default :
		  errflg++;
	  }
	}
  }
  if((argc - optind < 2 || argc > 20 ))
  {
	fprintf(stderr,"* argc = %d optind = %d\n",
			argc, optind);
	usage(ProgName);
  }

  // Read rest of input params
  InitialFragId = STR_TO_UID(argv[optind++], &end2, 10);
  FragId = InitialFragId;

  if (*end1 != '\0' || *end2 != '\0')
  { fprintf(stdout,"\n\t*** %s: Arguments must be integers\n",ProgName);
  exit(1);
  }

  if(argc <= optind){
    MinPrefixUnclear = 0;
    MaxPrefixUnclear = 50;
    MinSuffixUnclear = 0;
    MaxSuffixUnclear = 50;
  }else{
    MinPrefixUnclear = atoi(argv[optind++]);
    MaxPrefixUnclear = atoi(argv[optind++]);
    MinSuffixUnclear = atoi(argv[optind++]);
    MaxSuffixUnclear = atoi(argv[optind++]);
  }

  assert(MinPrefixUnclear >= 0 && MinPrefixUnclear <= MaxPrefixUnclear);
  assert(MinSuffixUnclear >= 0 && MinSuffixUnclear <= MaxSuffixUnclear);

  fprintf(stderr,"* Massage with DistID = " F_UID
	  " and FragID = " F_UID " clear ranges (%d,%d) (%d,%d)\n",
		  DistId, FragId,
		  MinPrefixUnclear,
		  MaxPrefixUnclear,
		  MinSuffixUnclear,
		  MaxSuffixUnclear);

  // Initialize dummyQuality
  {
    int i = 0;
    srand48(1);
    for(i = 0; i < 2 * (SEGLEN + 1); i++){
      int rand = GetRand_AS(0,19,TRUE);
      dummyQuality[i] = '0' + 40 + rand;
    }
  }

  // Initialize BACQuality
  {
    int i = 0;
    for(i = 0; i < 2 * (SEGLEN + 1); i++){
      BACQuality[i] = '0' + 10;
    }
  }


#ifdef NEVER
  time(&tp); /* int32 seconds from sometime in 1970 */
  tp = 0;  /* A reproducible timestamp for regression tests. */
#endif

  OutFrgMesg.action = AS_ADD;
  OutFrgMesg.source = Source;
  OutFrgMesg.quality = Quality;
  OutFrgMesg.sequence = Sequence;
  OutMesg.m = &OutFrgMesg;
  OutMesg.t = MESG_FRG;


  OutLnkMesg.m = &LnkMesg;
  OutLnkMesg.t = MESG_LKG;
  LnkMesg.link_orient = AS_INNIE;

  if (BacEnds)
	LnkMesg.type = AS_BAC_GUIDE;
  else if (LBacs)
  {
	OutFrgMesg.type = AS_LBAC;
	OutFrgMesg.elocale = elocale;
	OutFrgMesg.eseq_id = eseq_id;
  }
  else if (UBacs)
  {
	OutFrgMesg.type = AS_UBAC;
	OutFrgMesg.elocale = elocale;
	OutFrgMesg.eseq_id = eseq_id;
	OutFrgMesg.ebactig_id = ebactig_id;
  }
  else if (shreddedFBacs)
  {
	OutFrgMesg.type = AS_FBAC;
	OutFrgMesg.elocale = elocale;
	OutFrgMesg.eseq_id = eseq_id;
	OutFrgMesg.ebactig_id = ebactig_id;
  }
  else if (Bactigs)
  {
	OutFrgMesg.type = AS_BACTIG;
	OutFrgMesg.elocale = elocale;
	OutFrgMesg.eseq_id = eseq_id;
	OutFrgMesg.ebactig_id = ebactig_id;
  }
  else if (fullBacs)
  {
	OutFrgMesg.type = AS_FULLBAC;
	OutFrgMesg.elocale = elocale;
	OutFrgMesg.eseq_id = eseq_id;
	OutFrgMesg.ebactig_id = ebactig_id;
  }
  else
	LnkMesg.type = AS_MATE;
  LnkMesg.action = AS_ADD;


/* Pipe stdin, transform, and write to stdout */

  process();

  exit(0);
}


/*
   if(BacEnds)
   OutFrgMesg.type = AS_EBAC;
   else if (LBacs)
   OutFrgMesg.type = AS_LBAC;
   else if (UBacs)
   OutFrgMesg.type = AS_UBAC;
   else if (shreddedFBacs)
   OutFrgMesg.type = AS_FBAC;
   else if (Bactigs)
   OutFrgMesg.type = AS_BACTIG;
   else if (fullBacs)
   OutFrgMesg.type = AS_FULLBAC;
   else
   OutFrgMesg.type = AS_READ;
*/
