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


#include "AtacDataset.h"
#include <assert.h>
#include <string>
#include <ctype.h>

static string supportedFormat = "atac";
static string supportedVersion = "1.0";

AtacDataset::AtacDataset ()
{
  _format = supportedFormat;
  _version = supportedVersion;
  _inStream = 0;
  _outStream = 0;
}

AtacDataset::~AtacDataset ()
{
  if (_inStream != 0)
    fclose(_inStream);
  if (_outStream != 0)
    fclose(_outStream);
}

string AtacDataset::getFormat() const { return _format; }
string AtacDataset::getVersion() const { return _version; }
bool AtacDataset::formatIsSupported() const
{
  return getFormat() == supportedFormat && getVersion() == supportedVersion;
}

void AtacDataset::setInputFile(const string fileNameIn)
{
  _fileNameIn = fileNameIn;
  _inStream = fopen(_fileNameIn.c_str(), "r");
  if (_inStream == 0)
  {
    cout << "Unable to open input file: " << fileNameIn << ", exiting..." << endl;
    exit(1);
  }

  // read the first line to determine file format
  readNewNextRecord();
  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::setOutputFile(const string fileNameOut)
{
   _fileNameOut = fileNameOut;
   _outStream = fopen(_fileNameOut.c_str(), "w");
  if (_outStream == 0)
  {
    cout << "Unable to open output file: " << fileNameOut << ", exiting..." << endl;
    exit(1);
  }
}

void AtacDataset::readAttributes()
{
  char line[1024];
  char *key;
  char *equals;
  char *newline;
  char *value;

  // read line into buffer
  char *result;
  while ((result = fgets(line, 1024, _inStream)) != NULL)
  {
    char classCode = line[0];
    if (classCode == '/')
    {
      // cout << "Attribute: " << line;
      equals = index(line+1, '=');
      newline = index(line+1, '\n');
      *newline = '\0';
      value = equals + 1;
      *equals = '\0';
      key = line+1;
      string keyString = key;
      string valueString = value;
      putAttribute(keyString, valueString);
    }
  }
  evaluateAttributes();

  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::writeAttributes()
{
   map<string, string, less<string> >::iterator i;

   for (i=_keysToValues.begin(); i != _keysToValues.end(); i++)
   {
      string key = (*i).first;
      string value = (*i).second;
      // cout << "Write: " << key << ", " << value << "\n";
      fprintf(_outStream, "/%s=%s\n", key.c_str(), value.c_str());
   }
}

void AtacDataset::putAttribute(AtacAttribute attribute)
{
   string key = attribute.getKey();
   string value = attribute.getValue();
   _keysToValues[key] = value;
}

void AtacDataset::evaluateAttributes()
{
  // create an AtacSequenceSource if specified
  if (hasAttribute("assemblyFilePrefix1") && hasAttribute("assemblyId1") &&
      hasAttribute("assemblyFilePrefix2") && hasAttribute("assemblyId2"))
  {
    string fileNamePrefix1 = getAttributeString("assemblyFilePrefix1");
    string fileName1 = fileNamePrefix1 + ".fasta";
    string id1 = getAttributeString("assemblyId1");
    string fileNamePrefix2 = getAttributeString("assemblyFilePrefix2");
    string fileName2 = fileNamePrefix2 + ".fasta";
    string id2 = getAttributeString("assemblyId2");

    AtacSequenceSource *source1 = new AtacSequenceSource(id1, fileName1, "FASTA", "DNA", 0, 0);
    source1->setImplicit(true);
    addSequenceSource(source1);

    AtacSequenceSource *source2 = new AtacSequenceSource(id2, fileName2, "FASTA", "DNA", 0, 0);
    source2->setImplicit(true);
    addSequenceSource(source2);
  }
}

void AtacDataset::putAttribute(string key, string value)
{
   AtacAttribute attribute(key, value);
   putAttribute(attribute);
}

void AtacDataset::putAttribute(string key, int value)
{
   char valueBuffer[32];
   sprintf(valueBuffer, "%d", value);
   string valueString = valueBuffer;
   AtacAttribute attribute(key, valueString);
   putAttribute(attribute);
}

void AtacDataset::putAttribute(string key, int64 value)
{
   char valueBuffer[64];
   sprintf(valueBuffer, "%lld", value);
   string valueString = valueBuffer;
   AtacAttribute attribute(key, valueString);
   putAttribute(attribute);
}

void AtacDataset::putAttribute(string key, double value)
{
   char valueBuffer[32];
   sprintf(valueBuffer, "%g", value);
   string valueString = valueBuffer;
   AtacAttribute attribute(key, valueString);
   putAttribute(attribute);
}
string AtacDataset::getAttributeString(string key)
{
   return _keysToValues[key];
}

bool AtacDataset::hasAttribute(string key)
{
   return _keysToValues.find(key) != _keysToValues.end();
}


void AtacDataset::readComments()
{
  // Read records until end of the input file, starting from beginning of file
  AtacRecord *record;
  while ((record = readNewNextRecord()) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == '#')
      {
        AtacComment *comment = (AtacComment*)record;
        addComment(comment);
      }
    }
  }

  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::writeComments()
{
   vector<AtacComment*>::iterator i;

   // write out version information
   fprintf(_outStream, "! format %s %s\n", supportedFormat.c_str(), supportedVersion.c_str());
   for (i=_comments.begin(); i != _comments.end(); i++)
   {
      //cout << "Comment Pointer: " << *i << "\n";
      (*i)->write(_outStream);
   }
}

void AtacDataset::readMatches()
{
  // Read records until end of the input file, starting from beginning of file
  AtacRecord *record;
  while ((record = readNewNextRecord()) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == 'M')
      {
        AtacMatch *match = (AtacMatch*)record;
        addMatch(match);
      }
    }
  }

  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::computeMatchChildren()
{
   vector<AtacMatch*>::iterator i;

   for (i=_matches.begin(); i != _matches.end(); i++)
   {
     if ((*i) != 0)
     {
       AtacMatch *child = *i;
       AtacMatch *parent = getMatch(child->getParentId());
       if (parent != 0)
       {
         parent->addChild(child->getId());
       }
     }
   }
}

void AtacDataset::computeImplicitRecords()
{
/*
   vector<AtacMatch*>::iterator i;

   for (i=_matches.begin(); i != _matches.end(); i++)
   {
     if ((*i) != 0)
     {
       AtacMatch *child = *i;
       AtacMatch *parent = getMatch(child->getParentId());
       if (parent != 0)
       {
         parent->addChild(child->getId());
       }
     }
   }
*/
}

void AtacDataset::writeMatches()
{
   vector<AtacMatch*>::iterator i;

   for (i=_matches.begin(); i != _matches.end(); i++)
   {
     if ((*i) != 0 && !(*i)->isImplicit())
       (*i)->write(_outStream);
   }
}

void AtacDataset::writeSequenceSources()
{
   vector<AtacSequenceSource*>::iterator i;

   // cout << "writing sequence sources..." << "\n";
   for (i=_sources.begin(); i != _sources.end(); i++)
   {
     if ((*i) != 0 && !(*i)->isImplicit())
       (*i)->write(_outStream);
   }
}

void AtacDataset::writeGenomicAxes()
{
   vector<AtacGenomicAxis*>::iterator i;

   // cout << "writing genomic axes..." << "\n";
   for (i=_axes.begin(); i != _axes.end(); i++)
   {
     if ((*i) != 0 && !(*i)->isImplicit())
       (*i)->write(_outStream);
   }
}

void AtacDataset::readGlobalData()
{
  // Read records until end of the input file, starting from beginning of file
  AtacRecord *record;
  while ((record = readNewNextRecord(true)) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == 'S')
      {
        AtacSequenceSource *source = (AtacSequenceSource*)record;
        addSequenceSource(source);
      }
      else if (record->getClassCode() == 'G')
      {
        AtacGenomicAxis *axis = (AtacGenomicAxis*)record;
        addGenomicAxis(axis);
      }
      else if (record->getClassCode() == '#')
      {
        AtacComment *comment = (AtacComment*)record;
        addComment(comment);
      }
      else if (record->getClassCode() == '/')
      {
        AtacAttribute *attribute= (AtacAttribute*)record;
        putAttribute(*attribute);
      }
      else
        delete record;
    }
    else
    {
      cout << "Bad Record: " << endl;
      delete record;
    }
  }
  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::readRecords()
{
  // Read records until end of the input file, starting from beginning of file
  AtacRecord *record;
  while ((record = readNewNextRecord()) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == 'M')
      {
        AtacMatch *match = (AtacMatch*)record;
        addMatch(match);
      }
      else if (record->getClassCode() == 'S')
      {
        AtacSequenceSource *source = (AtacSequenceSource*)record;
        addSequenceSource(source);
      }
      else if (record->getClassCode() == 'G')
      {
        AtacGenomicAxis *axis = (AtacGenomicAxis*)record;
        addGenomicAxis(axis);
      }
      else if (record->getClassCode() == '#')
      {
        AtacComment *comment = (AtacComment*)record;
        addComment(comment);
      }
      else if (record->getClassCode() == '/')
      {
        AtacAttribute *attribute= (AtacAttribute*)record;
        putAttribute(*attribute);
      }
      else
        delete record;
    }
    else
    {
      cout << "Bad Record: " << endl;
      delete record;
    }
  }
  evaluateRecords();

  // reset stream
  fclose(_inStream);
  _inStream = fopen(_fileNameIn.c_str(), "r");
}

void AtacDataset::evaluateRecords()
{
  evaluateAttributes();
  computeMatchChildren();
}

void AtacDataset::writeRecords()
{
  writeComments();
  writeAttributes();
  writeSequenceSources();
  writeGenomicAxes();
  writeMatches();
}

void AtacDataset::writeGlobalData()
{
  writeComments();
  writeAttributes();
  writeSequenceSources();
  writeGenomicAxes();
}

void AtacDataset::addComment(AtacComment *comment)
{
  _comments.push_back(comment);
}

void AtacDataset::addSequenceSource(AtacSequenceSource *source)
{
  source->setRowIndex(_nsources++);
  _sources.push_back(source);
  string id = source->getId();
  _idsToSources[id] = source;
}

void AtacDataset::addGenomicAxis(AtacGenomicAxis *axis)
{
  axis->setRowIndex(_naxes++);
  _axes.push_back(axis);
  string id = axis->getId();
  _idsToAxes[id] = axis;
}

AtacSequenceSource *AtacDataset::getSequenceSource(string id)
{
  return _idsToSources[id];
}

bool AtacDataset::hasGenomicAxis(string id) const
{
  return _idsToAxes.find(id) != _idsToAxes.end();
}

AtacGenomicAxis *AtacDataset::getGenomicAxis(string id)
{
  return _idsToAxes[id];
}

void AtacDataset::addMatch(AtacMatch *match)
{
  match->setRowIndex(_nmatches++);
  _matches.push_back(match);
  string id = match->getId();
  _idsToMatches[id] = match;

  // calculate any implied genomic axes
  addImplicitGenomicAxes(match);
}

void AtacDataset::addImplicitGenomicAxes(AtacMatch *match)
{
  // check if genomic axis exists, if not, create one
  for (int which=0; which<2; which++)
  {
    string gaid = match->getGenomicAxisId(which);
    if (!hasGenomicAxis(gaid))
    {
      string ssid = ".";
      int idot;
      if ((idot = gaid.find(":", 0)) != string::npos)
      {
        // cout << "creating genomic axis " << gaid << endl;
        ssid = gaid.substr(0, idot);
        string ordinalString = gaid.substr(idot+1, 32);
        int ordinal = atoi(ordinalString.c_str());
        AtacGenomicAxis *newAxis = new AtacGenomicAxis(".", gaid, ".", 0, ssid, 0L, 0L, 0L, ordinal, "", "");
        newAxis->setImplicit(true);
        addGenomicAxis(newAxis);
      }
    }
  }
}

void AtacDataset::removeMatch(AtacMatch *match)
{
  string id = match->getId();
  _idsToMatches.erase(id);
  int rowIndex = match->getRowIndex();
  // null out the element, rather than erase, to avoid compacting
  _matches[rowIndex] = 0;
}

AtacMatch *AtacDataset::getMatch(string id)
{
  return _idsToMatches[id];
}

const vector<AtacMatch*> &AtacDataset::getMatches()
{
  return _matches;
}

bool AtacDataset::readNextMatch(AtacMatch &match)
{
  AtacRecord *record;
  while ((record = readNewNextRecord()) != 0)
  {
    if (record->isValid())
    {
      if (record->getClassCode() == 'M')
      {
        match = *(AtacMatch*)record;
    	delete record;
        return 1;
      }
    }
    delete record;
  }

  return 0;
}

AtacRecord *AtacDataset::readNewNextRecord()
{
  return readNewNextRecord(false);
}

AtacRecord *AtacDataset::readNewNextRecord(bool filterMatches)
{
  static char line[1024];

  static char typeId[128];
  static char selfId[128];
  static char parentId[128];

  static char axisId1[128];
  static int64 start1;
  static int64 length1;
  static int orient1;

  static char axisId2[128];
  static int64 start2;
  static int64 length2;
  static int orient2;

  static char fileName[512];
  static char fileType[128];
  static char elemType[128];
  static int64 width;
  static int trim;

  static int64 parentStart;
  static char sequenceSourceId[128];
  static int64 start;
  static int64 length;
  static int64 offset;
  static int ordinal;
  static char ordinalString[512];
  static char name[512];
  static char uid[512];

  static char commandString[128];
  static char valueString1[128];
  static char valueString2[128];

  // read line into buffer
  char *result;
  while ((result = fgets(line, 1024, _inStream)) != NULL)
  {
    // cout << "Result: " << line;
    char classCode = line[0];
    if (classCode == 'M')
    {
      if (filterMatches)
        continue;

      // find defline, if it exists, and null out '>'
      char *anglebracket = index(line+1, '>');
      if (anglebracket != 0)
      {
        char *newline = index(line+1, '\n');
        *newline = '\0';
      }

      sscanf(line+1, "%s %s %s %s %lld %lld %d %s %lld %lld %d", 
	     typeId, selfId, parentId, 
	     axisId1, &start1, &length1, &orient1,
	     axisId2, &start2, &length2, &orient2);

      string tid = typeId;
      string sid = selfId;
      string pid = parentId;
      string aid1 = axisId1;
      string aid2 = axisId2;

      AtacMatch *newMatch = new AtacMatch(dotToEmpty(tid), sid, dotToEmpty(pid), dotToEmpty(aid1), start1, length1, orient1, dotToEmpty(aid2), start2, length2, orient2);

      // cout << *newMatch << "\n";

      if (anglebracket != 0)
      {
        string defline = anglebracket+1;
      	newMatch->setDefline(defline);
      }

      if (orient1 > 1 || orient1 < -1 || orient2 > 1 || orient2 < -1 || length1 < 0 || length2 < 0)
      {
        // cout << "invalid\n";
        newMatch->invalidate();
      }


      // match created, return it
      return newMatch;
    }
    else if (classCode == 'S')
    {
      // find defline, if it exists, and null out '>'
      char *anglebracket = index(line+1, '>');
      if (anglebracket != 0)
      {
        char *newline = index(line+1, '\n');
        *newline = '\0';
      }

      sscanf(line+1, "%s %s %s %s %lld %d", selfId, fileName, fileType, elemType, &width, &trim); 

      string sid = selfId;
      string fn = fileName;
      string ft = fileType;
      if (ft.size() == 0)
        ft = "FASTA";
      string et = elemType;
      if (et.size() == 0)
        et = "DNA";

      AtacSequenceSource *newSource = new AtacSequenceSource(sid, fn, ft, et, width, trim);

      // cout << *newSource << "\n";

      if (anglebracket != 0)
      {
        string defline = anglebracket+1;
      	newSource->setDefline(defline);
      }

      // sequence source created, return it
      return newSource;
    }
    else if (classCode == 'G')
    {
      // find defline, if it exists, and null out '>'
      char *anglebracket = index(line+1, '>');
      if (anglebracket != 0)
      {
        char *newline = index(line+1, '\n');
        *newline = '\0';
      }

      sscanf(line+1, "%s %s %s %lld %s %lld %lld %s %s %s", typeId, selfId, parentId, &parentStart, sequenceSourceId, &start, &length, &offset, &ordinalString, name, uid); 

      string tid = typeId;
      string sid = selfId;
      string pid =  parentId;
      string ssid = sequenceSourceId;

      // this is a hack to support the old CGF format, which lacks an "ordinal" column
      bool newFormat = isdigit(ordinalString[0]);
      ordinal = newFormat ? atoi(ordinalString) : 0;
      string nameString = newFormat ? name : ordinalString;
      string uidString = newFormat ? uid : name;
      string n = nameString;
      string u = uidString;

      AtacGenomicAxis *newAxis = new AtacGenomicAxis(dotToEmpty(tid), sid, dotToEmpty(pid), parentStart, dotToEmpty(ssid), start, length, offset, ordinal, dotToEmpty(n), dotToEmpty(u));

      // cout << *newSource << "\n";

      if (anglebracket != 0)
      {
        string defline = anglebracket+1;
      	newAxis->setDefline(defline);
      }

      // sequence source created, return it
      return newAxis;
    }
    else if (classCode == '#')
    {
       char *newline = index(line+1, '\n');
       *newline = '\0';
       string text = line+1;
       AtacComment *newComment = new AtacComment(text);

       return newComment;
    }
    else if (classCode == '/')
    {
      char *key;
      char *equals;
      char *newline;
      char *value;

      // cout << "Attribute: " << line;
      equals = index(line+1, '=');
      newline = index(line+1, '\n');
      *newline = '\0';
      value = equals + 1;
      *equals = '\0';
      key = line+1;
      string keyString = key;
      string valueString = value;
      AtacAttribute *newAttribute = new AtacAttribute(key, value);

      return newAttribute;
    }
    else if (classCode == '!')
    {
      sscanf(line+1, "%s %s %s", commandString, valueString1, valueString2);
      string command = commandString;
      string value1 = valueString1;
      string value2 = valueString2;
      if (command == "format")
      {
        _format = value1;
        _version = value2;
      }
    }
  }

  // end of file reached, return null value
  return 0;
}

void AtacDataset::writeNextRecord(AtacRecord *record)
{
  record->write(_outStream);
}
