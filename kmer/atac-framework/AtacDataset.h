/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Ross Lippert
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
**************************************************************************/


#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

#include "AtacMatch.h"
#include "AtacSequenceSource.h"
#include "AtacGenomicAxis.h"
#include "AtacComment.h"
#include "AtacAttribute.h"

#ifndef _INCLUDE_ATACDATASET
#define _INCLUDE_ATACDATASET

#ifndef _USING_NAMESPACE_STD_
#define _USING_NAMESPACE_STD_
using namespace std;
#endif

/**
 * AtacDataset
 **/
class AtacDataset : public AtacRecord
{
 public:
  AtacDataset();
  AtacDataset(const string file_name, bool read_primary_left);
  ~AtacDataset();

  // version information
  string getFormat() const;
  string getVersion() const;
  bool formatIsSupported() const;

  /* Set the input file */
  void setInputFile(const string fileNameIn);

  /* Set the output file */
  void setOutputFile(const string fileNameOut);

  /* Read all the global data, i.e. non-Matches, scanning the entire file */
  void readGlobalData();

  /* Read all records from the input file */
  void readRecords();

  /* Write all global data to the output file, i.e. non Match records */
  void writeGlobalData();

  /* Write all records to the output file */
  void writeRecords();

  /* Add an allocated sequence source to the dataset by pointer */
  void addSequenceSource(AtacSequenceSource *source);

  /* Get a sequence source by pointer */
  AtacSequenceSource *getSequenceSource(string id);

  /* Add an allocated genomic axis to the dataset by pointer */
  void addGenomicAxis(AtacGenomicAxis *axis);

  /* Does the specified genomic axis exist? */
  bool hasGenomicAxis(string id) const;

  /* Get a genomic axis by pointer */
  AtacGenomicAxis *getGenomicAxis(string id);

  /* Add any implied genomic axes that don't exist yet */
  void addImplicitGenomicAxes(AtacMatch *match);

  /* Add an allocated match to the dataset by pointer */
  void addMatch(AtacMatch *match);

  /* Remove a match from the dataset by pointer (does not delete) */
  void removeMatch(AtacMatch *match);

  /* Return a pointer to a match, specified by its id string */
  AtacMatch *getMatch(string id);

  /* Return a constant reference to a vector of match pointers  */
  const vector<AtacMatch*> &getMatches();

  /* Add a comment record */
  void addComment(AtacComment *comment);

  /* Add a global attribute key-value pair */
  void putAttribute(AtacAttribute attribute);
  void putAttribute(string key, string value);
  void putAttribute(string key, int value);
  void putAttribute(string key, int64 value);
  void putAttribute(string key, double value);

  /* True if attribute exists for key */
  bool hasAttribute(string key);

  /* add a global attribute by key */
  string getAttributeString(string key);

  /* 
   * Read a single record in from the current location of the file, advancing
   * the file pointer to the next record in the file. Returns a newly allocated
   * instance of an AtacRecord, which the client is responsible for deleting. Does
   * not maintain any state information about the record. Client is responsible for
   * adding record information to the dataset.
   *
   * Returns 0 if end of file has been reached.
   *
   * This is the lowest level parsing method, and doesn't maintain any state.
   */
  AtacRecord *readNewNextRecord();

  /*
   * Read the next match from the file and copy into match parameter, skipping over
   * any other types of Records in the file. Does not maintain any state information about the record. Client is responsible for
   * adding record information to the dataset.
   *
   * Returns 0 if end of file has been reached.
   *
   * This is the preferred method for parsing streaming data, and should be called
   * after calling readGlobalData().
   */
  bool readNextMatch(AtacMatch &match);

  /* Write a single record to the end of the output file */
  void writeNextRecord(AtacRecord *record);

 private:
  /* Same as parameterless version, except that true input  value forces it to filter out Matches */
  AtacRecord *readNewNextRecord(bool filterMatches);

  /* Read all the global attributes, scanning the entire file */
  void readAttributes();

  /* Write all the global attributes out to the end of the current output file */
  void writeAttributes();

  /* Read all comments from the input file */
  void readComments();

  /* Write all comments to the output file */
  void writeComments();

  /* Read all matches from the input file */
  void readMatches();

  /* Create implicit genomic axes and sequence sources */
  void computeImplicitRecords();

  /* Determine the matches children from parent links */
  void computeMatchChildren();

  /* Write all sequence sources to the output file */
  void writeSequenceSources();

  /* Write all genomic axes to the output file */
  void writeGenomicAxes();

  /* Write all matches to the output file */
  void writeMatches();

  /* Evaluate any hard-coded information in the records */
  void evaluateRecords();

  /* Evaluate any hard-coded information in the attributes */
  void evaluateAttributes();


  /************************ member variables ***********************/

  string _format;
  string _version;

  FILE *_inStream;
  FILE *_outStream;

  string _fileNameIn;
  string _fileNameOut;

  bool _validating;

  int _nmatches;
  int _nsources;
  int _naxes;

  map<string, string, less<string> > _keysToValues;
  map<string, AtacMatch*, less<string> > _idsToMatches;
  map<string, AtacSequenceSource*, less<string> > _idsToSources;
  map<string, AtacGenomicAxis*, less<string> > _idsToAxes;

  vector<AtacComment*> _comments;
  vector<AtacMatch*> _matches;
  vector<AtacSequenceSource*> _sources;
  vector<AtacGenomicAxis*> _axes;

};

#endif
