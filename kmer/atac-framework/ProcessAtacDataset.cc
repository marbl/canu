// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Ross Lippert
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#include <string>
#include <time.h>
#include "AtacMatch.h"
#include "AtacComment.h"
#include "AtacAttribute.h"
#include "AtacDataset.h"

using namespace std;

/**
 * Parse the ATAC file "fileNameIn" and output it to "fileNameOut", after compiling
 * statistics and filtering out any bad data. This version doesn't preserve the 
 * ordering of comments and global variables
 *
 * STREAMING VERSION
 **/
void processAtacFileStreaming(const string fileNameIn, const string fileNameOut)
{
  cout << "Processing " << fileNameIn << " by streaming, without preserving order..." << endl;
  int nmatches = 0;

  // Create new dataset instance
  AtacDataset *dataset = new AtacDataset();

  // Open the input and output files for reading
  dataset->setInputFile(fileNameIn);
  dataset->setOutputFile(fileNameOut);

  // Read global information, making a complete scan but ignoring Match records
  cout << "Reading global information..." << endl;
  dataset->readGlobalData();

  // Write global information back out at the beginning of the file
  cout << "Writing global information..." << endl;
  dataset->writeGlobalData();

  // Read matches until end of the input file, starting from beginning of file
  cout << "Reading matches and writing out to " << fileNameOut << "..." << endl;
  AtacMatch match;
  while (dataset->readNextMatch(match))
  {
    if (match.isValid())
    {
      // filter out matches with odd starting position
      if (match.getStart(1) % 2 == 0)
      {
        dataset->writeNextRecord(&match);
        nmatches++;
      }
    } 
    else
      cout << "Bad Match: " << match << endl;
  }

  // add a global attribute at the end of the file to record the number of matches
  cout << "Writing additional global attribute..." << endl;
  AtacAttribute attribute("matchCount", nmatches);
  dataset->putAttribute(attribute);
  dataset->writeNextRecord(&attribute);
  cout << "Finished writing " << fileNameOut << endl;
}
/**
 * Parse the ATAC file "fileNameIn" and output it to "fileNameOut", after compiling
 * statistics and filtering out any bad data. This version preserves the ordering
 * of comments and global variables
 * PRESERVE ORDER VERSION
 **/
void processAtacFilePreserveOrder(const string fileNameIn, const string fileNameOut)
{
  cout << "Processing " << fileNameIn << " by streaming, preserving order..." << endl;
  int nmatches = 0;

  // Create new dataset instance
  AtacDataset *dataset = new AtacDataset();

  // Open the input and output files
  dataset->setInputFile(fileNameIn);
  dataset->setOutputFile(fileNameOut);

  // Read global information, excluding Matches. Makes complete scan of file.
  cout << "Reading global information..." << endl;
  dataset->readGlobalData();

  // Read records until end of the input file, starting from beginning of file
  cout << "Reading records and writing out to " << fileNameOut << "..." << endl;
  AtacRecord *record;
  while ((record = dataset->readNewNextRecord()) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == 'M')
      {
        AtacMatch *match = (AtacMatch*)record;
        // filter out matches with odd starting position
        if (match->getStart(1) % 2 == 0)
          dataset->writeNextRecord(match);
        nmatches++;
      }
      else if (record->getClassCode() == 'S')
      {
        AtacSequenceSource *source = (AtacSequenceSource*)record;
        dataset->writeNextRecord(source);
      }
      else if (record->getClassCode() == 'G')
      {
        AtacGenomicAxis *axis= (AtacGenomicAxis*)record;
        dataset->writeNextRecord(axis);
      }
      else if (record->getClassCode() == '#')
      {
        AtacComment *comment = (AtacComment*)record;
        dataset->writeNextRecord(comment);
      }
      else if (record->getClassCode() == '/')
      {
        AtacAttribute *attribute= (AtacAttribute*)record;
        dataset->writeNextRecord(attribute);
      }
    } 
    else
      cout << "Bad Record: " << *record << endl;

    delete record;
  }

  // add a global attribute at the end of the file to record the number of matches
  cout << "Writing additional global attribute..." << endl;
  AtacAttribute attribute("matchCount", nmatches);
  dataset->putAttribute(attribute);
  dataset->writeNextRecord(&attribute);
  cout << "Finished writing " << fileNameOut << endl;
}


/**
 * Parse the ATAC file "fileNameIn" and output it to "fileNameOut", after compiling
 * statistics and filtering out any bad data
 * IN-CORE VERSION
 **/
void processAtacFileInCore(const string fileNameIn, const string fileNameOut)
{
  cout << "Processing " << fileNameIn << " in-core..." << endl;

  // Create new AtacDataset
  AtacDataset *dataset = new AtacDataset();

  // Open the input and output files
  dataset->setInputFile(fileNameIn);
  dataset->setOutputFile(fileNameOut);

  // Read records until end of the input file, starting from beginning of file
  cout << "Reading all records into memory..." << endl;
  dataset->readRecords();

  // Scan through the matches and check for validity
  cout << "Checking for bad data..." << endl;
  const vector<AtacMatch*> &matches = dataset->getMatches();

  // To get a copy of the vector that you can modify, do the following
  // vector<AtacMatch*> mutableMatches = dataset->getMatches();

  int64 nmatches = matches.size();
  int i;
  for (i=0; i<nmatches; i++)
  {
    if (matches[i]!= 0)
    {
      if (!matches[i]->isValid())
        cout << "Bad Match: " << matches[i] << endl;
    }
  }

  // check for matches without parents
  for (i=0; i<nmatches; i++)
  {
    if (matches[i] != 0)
    {
      cout << "Parent: " << matches[i]->getParentId() << endl;
      if (matches[i]->getParentId() == "")
        cout << "Match without parent: " << *matches[i] << endl;
    }
  }

  // print out children relationships
  cout << "Matches with children:" << endl;
  for (i=0; i<nmatches; i++)
  {
    AtacMatch *aMatch = matches[i];
    if (aMatch != 0 && aMatch->hasChildren())
    {
      cout << "Match " <<  aMatch->getId() << ": ";
      vector<string> *children = aMatch->getChildren();
      int nchildren = children->size();
      for (int j=0; j<nchildren; j++)
        cout << " " << (*children)[j];
      cout << endl;
    }
  }

  // Examine a specific match
  AtacMatch *aMatch = dataset->getMatch("m3");
  if (aMatch != 0)
  {
    cout << "Match 3: " << *aMatch << endl;

    // Get the GenomicAxis and SequenceSource information
    // NOTE: this will work in the streaming versions a well
    AtacGenomicAxis *axis1 = dataset->getGenomicAxis(aMatch->getGenomicAxisId(0));
    AtacGenomicAxis *axis2 = dataset->getGenomicAxis(aMatch->getGenomicAxisId(1));
    AtacSequenceSource *source1 = dataset->getSequenceSource(axis1->getSequenceSourceId());
    AtacSequenceSource *source2 = dataset->getSequenceSource(axis2->getSequenceSourceId());
    cout << "Axis 1 File: " << source1->getFileName() << endl;
    cout << "Axis 1 Ordinal: " << axis1->getOrdinal() << endl;
    cout << "Axis 2 File: " << source2->getFileName() << endl;
    cout << "Axis 2 Ordinal: " << axis2->getOrdinal() << endl;

    // Remove the Match
    // NOTE: this nulls out the pointer in the vector returned by getMatches(), only
    cout << "Remove Match 3..." << endl;
    dataset->removeMatch(aMatch);
  }

  // add a global attribute to record the number of matches, and some test attributes
  cout << "Storing additional global attributes..." << endl;
  char nmatchesString[32];
  sprintf(nmatchesString, "%d", nmatches);
  dataset->putAttribute("matchCount", nmatchesString);
  dataset->putAttribute("matchCount2", nmatches);
  // Following line wouldn't compile under linux -- Jason
  //dataset->putAttribute("testLong", 12345678912345L); 
  dataset->putAttribute("testFloat", 1.2345);

  // Write all the data 
  cout << "Outputting to " << fileNameOut << endl;
  cout << "Writing records..." << endl;
  dataset->writeRecords();
  cout << "Finished writing " << fileNameOut << endl;
}

/**
 * Parse the ATAC file "fileNameIn" and output it to "fileNameOut" for
 * Timing test.
 **/
void timeAtacFile(const string fileNameIn, const string fileNameOut)
{
  cout << "Processing " << fileNameIn << " in-core..." << endl;

  // Create new AtacDataset
  AtacDataset *dataset = new AtacDataset();

  // Open the input file for reading
  dataset->setInputFile(fileNameIn);
  if (!dataset->formatIsSupported())
  {
    cout << "File format " << dataset->getFormat() << " version " << dataset->getVersion() << " not supported." << endl;
    exit(1);
  }

  // Open the output file for reading
  dataset->setOutputFile(fileNameOut);

  // Read records until end of the input file, starting from beginning of file
  cout << "Reading all records into memory..." << endl;
  time_t startTime = time(0);
  dataset->readRecords();
  time_t endTime = time(0);
  int loadTime = endTime - startTime;
  cout << "Loaded in " << loadTime << " seconds." << endl;

  // calculate average number of children
  cout << "Average number of children for matches with children:";
  int nmatches = dataset->getMatches().size();
  int nparents = 0;
  int nchildren = 0;
  int i;
  for (i=0; i<nmatches; i++)
  {
    AtacMatch *aMatch = dataset->getMatches()[i];
    if (aMatch != 0 && aMatch->hasChildren())
    {
      nparents++;
      vector<string> *children = aMatch->getChildren();
      nchildren += children->size();
    }
  }
  double average = 0.0;
  if (nparents > 0)
    average = (double)nchildren / nparents;
  cout << average << endl;

  cout << "Outputting to " << fileNameIn << endl;

  // Write all the data
  cout << "Writing comments..." << endl;
  startTime = time(0);
  dataset->writeRecords();
  endTime = time(0);
  int saveTime = endTime - startTime;
  cout << "Saved in " << saveTime << " seconds." << endl;
}

/**
 * Parse the ATAC file "fileNameIn" and output it to "fileNameOut" for
 * Timing test.
 **/
void testGlobals(const string fileNameIn)
{
  cout << "Processing " << fileNameIn << " readGlobalData()..." << endl;

  // Create new AtacDataset
  AtacDataset *dataset = new AtacDataset();

  // Open the input file for reading
  dataset->setInputFile(fileNameIn);
  if (!dataset->formatIsSupported())
  {
    cout << "File format " << dataset->getFormat() << " version " << dataset->getVersion() << " not supported." << endl;
    exit(1);
  }

  // Read records until end of the input file, starting from beginning of file
  cout << "Reading globals..." << endl;
  time_t startTime = time(0);
  dataset->readGlobalData();
  time_t endTime = time(0);
  int loadTime = endTime - startTime;
  cout << "Globals Loaded in " << loadTime << " seconds." << endl;
}


/*
 * Process a test ATAC file and output the results
 */

int main (int argc, char ** argv)
{
  using namespace std;
  if (argc == 2 )
  {
    testGlobals(argv[1]);
/*
    cout << "Usage:" << endl;
    cout << "  ProcessAtacDataset                               (run example test programs)" << endl;
    cout << "  ProcessAtacDataset <input-file> <output-file>    (run timing test)" << endl;
*/
  }
  else if (argc > 2)
  {
    cout << "Starting ProcessAtacDataset timing test program..." << endl << endl;
    timeAtacFile(argv[1], argv[2]);
  }
  else
  {
    cout << "Starting ProcessAtacDataset test program..." << endl << endl;
    processAtacFileStreaming("testIn.atac", "testOutStreaming.atac");
    processAtacFilePreserveOrder("testIn.atac", "testOutPreserveOrder.atac");
    processAtacFileInCore("testIn.atac", "testOutInCore.atac");
  }
  return 0;
}
