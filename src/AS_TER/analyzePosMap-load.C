
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: analyzePosMap-load.C,v 1.1 2012-12-05 01:13:23 brianwalenz Exp $";

#include "analyzePosMap.H"



char                 frgLabel[NUM_FRG_LABELS] = { 'D', 'C', 'd', 'u', 'r', 's', 'R', 'p', 'S', 'o' };
uint32               frgLabelToInt[256]       = { 0 };

map<string,uint32>   libNames;
vector<libEntry>     libDat;

vector<uint32>       frgMate;
vector<uint32>       frgLibrary;

map<string,uint32>   frgNames;
map<string,uint32>   utgNames;
map<string,uint32>   ctgNames;
map<string,uint32>   scfNames;

vector<frgEntry>     frgDat;
vector<utgEntry>     utgDat;
vector<ctgEntry>     ctgDat;
vector<scfEntry>     scfDat;

vector<string>       frgNam;
vector<string>       utgNam;
vector<string>       ctgNam;
vector<string>       scfNam;


FILE *
openFile(char *prefix, char *mapname) {
  char  fileName[FILENAME_MAX];

  sprintf(fileName, "%s.posmap.%s", prefix, mapname);

  errno = 0;
  FILE *file = fopen(fileName, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", fileName, strerror(errno)), exit(1);

  return(file);
}




void
readLibraries(char *prefix) {
  FILE    *file     = openFile(prefix, "libraries");
  uint32   lineMax  = 1048576;
  char    *line     = new char [1048576];

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   libName(W[0]);
    uint32   libId = libNames.size();

    libNames[libName] = libId;

    libDat.push_back(libEntry());

    libDat[libId].initialize(W[0], W(1), W(2));

    fgets(line, lineMax, file);
    chomp(line);
  }

  delete [] line;
}



void
readFrags(char *prefix) {
  FILE    *file     = openFile(prefix, "frags");
  uint32   lineMax  = 1048576;
  char    *line     = new char [1048576];
  uint32   count    = 0;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   frgName(W[0]);
    string   libName(W[5]);

    uint32   frgId = frgNames.size();
    uint32   libId = libNames[libName];

    frgNames[frgName] = frgId;

    frgNam.push_back(frgName);
    frgDat.push_back(frgEntry());

    frgMate.push_back(0);  //  All are initially unmated.
    frgLibrary.push_back(libId);

    frgEntry &frg = frgDat.back();

    frg.len = W(2) - W(1);
    frg.typ = 'f';
    frg.sta = '?';

    if      (W[3][0] == 'p')
      frg.sta = '?';  //  Set later when reading frgscf, frgdeg and frgutg
    else if (W[3][0] == 'd')
      frg.sta = 'D';  //  Deleted
    else if (W[3][0] == 'c')
      frg.sta = 'C';  //  Singleton (chaff)
    else
      assert(0);

    assert(frgDat.size()  == frgId + 1);  //  Already added the object
    assert(frgMate.size() == frgId + 1);
    assert(frgNam.size( ) == frgId + 1);

    if ((++count % 1000000) == 0)
      fprintf(stderr, ".");
    if ((count % 10000000) == 0)
      fprintf(stderr, F_U32, count);

    fgets(line, lineMax, file);
    chomp(line);
  }

  fprintf(stderr, "\n");
  delete [] line;
}




void
readMates(char *prefix) {
  FILE    *file     = openFile(prefix, "mates");
  uint32   lineMax  = 1048576;
  char    *line     = new char [1048576];
  uint32   count    = 0;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   fr1Name(W[0]);
    string   fr2Name(W[1]);

    assert(frgNames.count(fr1Name) > 0);
    assert(frgNames.count(fr2Name) > 0);

    uint32   fr1Id = frgNames[fr1Name];
    uint32   fr2Id = frgNames[fr2Name];

    frgMate[fr1Id] = fr2Id;  //  Now they're mated!
    frgMate[fr2Id] = fr1Id;

    if ((++count % 1000000) == 0)
      fprintf(stderr, ".");
    if ((count % 10000000) == 0)
      fprintf(stderr, F_U32, count);

    fgets(line, lineMax, file);
    chomp(line);
  }

  fprintf(stderr, "\n");
  delete [] line;
}



void
createScf(char   *prefix,
          char   *posmap,
          char    status) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  scfEntry scf;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   scfName(W[0]);

    uint32   scfId = scfNames.size();

    scfNames[scfName] = scfId;

    assert(scfId == scfDat.size());

    scfNam.push_back(scfName);
    scfDat.push_back(scfEntry());

    scfEntry &scf = scfDat.back();

    scf.len = W(1);
    scf.typ = 's';
    scf.sta = status;

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}


void
createCtg(char   *prefix,
          char   *posmap,
          char    status) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  ctgEntry ctg;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   ctgName(W[0]);

    uint32   ctgId = ctgNames.size();

    ctgNames[ctgName] = ctgId;

    assert(ctgId == ctgDat.size());

    ctgNam.push_back(ctgName);
    ctgDat.push_back(ctgEntry());

    ctgEntry &ctg = ctgDat.back();

    ctg.len = W(1);
    ctg.typ = 'c';
    ctg.sta = status;

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}


void
createUtg(char   *prefix,
          char   *posmap,
          char    status) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  utgEntry utg;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    string   utgName(W[0]);

    uint32   utgId = utgNames.size();

    utgNames[utgName] = utgId;

    assert(utgId == utgDat.size());

    utgNam.push_back(utgName);
    utgDat.push_back(utgEntry());

    utgEntry &utg = utgDat.back();

    utg.len = W(1);
    utg.typ = 'u';
    utg.sta = status;

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}





uint32
parsePosMapLine(splitToWords        &W,
                mapEntry            &mpe,
                map<string,uint32>  &objIdMap,
                map<string,uint32>  &parIdMap) {
  string   objName(W[0]);
  string   parName(W[1]);

  uint32   objId = objIdMap[objName];
  uint32   parId = parIdMap[parName];

  mpe.con = parId;
  mpe.bgn = W(2);
  mpe.end = W(3);
  mpe.ori = W[4][0];

  return(objId);
}


char
parseUtgType(splitToWords  &W) {

  assert(W.numWords() == 6);

  if (W.numWords() != 6)
    return('X');

  if (strcmp(W[5], "unique") == 0)    return('u');
  if (strcmp(W[5], "rock") == 0)      return('r');
  if (strcmp(W[5], "stone") == 0)     return('s');
  if (strcmp(W[5], "pebble") == 0)    return('p');
  if (strcmp(W[5], "singleton") == 0) return('S');
  if (strcmp(W[5], "other") == 0)     return('o');

  fprintf(stderr, "WARNING: unknown unitig type '%s'\n", W[5]);
  return('X');
}





void
loadCtgScf(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   ctgId = parsePosMapLine(W, mpe, ctgNames, scfNames);
    uint32   scfId = mpe.con;

    ctgDat[ctgId].scf = mpe;
    scfDat[scfId].ctg.push_back(ctgId);

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
loadUtgScf(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   utgId = parsePosMapLine(W, mpe, utgNames, scfNames);
    uint32   scfId = mpe.con;

    utgDat[utgId].scf = mpe;
    scfDat[scfId].utg.push_back(utgId);

    utgDat[utgId].sta = parseUtgType(W);  //  Status is reset (to the same thing) in utgctg below.

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}

void
loadUtgCtg(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   utgId = parsePosMapLine(W, mpe, utgNames, ctgNames);
    uint32   ctgId = mpe.con;

    utgDat[utgId].ctg = mpe;
    ctgDat[ctgId].utg.push_back(utgId);

    utgDat[utgId].sta = parseUtgType(W);  //  Status was already set by utgscf above.

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
loadUtgDeg(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   utgId = parsePosMapLine(W, mpe, utgNames, ctgNames);
    uint32   ctgId = mpe.con;
    uint32   scfId = ctgId;

    utgDat[utgId].ctg = mpe;
    ctgDat[ctgId].scf = mpe;

    assert(utgDat[utgId].sta == '?');  //  Unitigs are not divided
    assert(ctgDat[ctgId].sta == 'd');
    assert(scfDat[scfId].sta == 'd');

    utgDat[utgId].sta = 'd';

    ctgDat[ctgId].utg.push_back(utgId);
    scfDat[scfId].ctg.push_back(ctgId);

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}


void
loadFrgDeg(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   frgId = parsePosMapLine(W, mpe, frgNames, ctgNames);
    uint32   ctgId = mpe.con;

    frgDat[frgId].ctg = mpe;
    ctgDat[ctgId].frg.push_back(frgId);

    frgDat[frgId].sta = 'd';

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
loadFrgScf(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   frgId = parsePosMapLine(W, mpe, frgNames, scfNames);
    uint32   scfId = mpe.con;

    frgDat[frgId].scf = mpe;
    scfDat[scfId].frg.push_back(frgId);

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}

void
loadFrgCtg(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   frgId = parsePosMapLine(W, mpe, frgNames, ctgNames);
    uint32   ctgId = mpe.con;

    frgDat[frgId].ctg = mpe;
    ctgDat[ctgId].frg.push_back(frgId);

    frgDat[frgId].sta = 'p';

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
loadFrgUtg(char    *prefix,
           char    *posmap) {
  FILE   *file    = openFile(prefix, posmap);
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  mapEntry mpe;

  fgets(line, lineMax, file);
  chomp(line);

  while (!feof(file)) {
    splitToWords W(line);

    uint32   frgId = parsePosMapLine(W, mpe, frgNames, utgNames);
    uint32   utgId = mpe.con;

    frgDat[frgId].utg = mpe;
    utgDat[utgId].frg.push_back(frgId);

    assert(frgDat[frgId].sta != 'D');  //  Can't be deleted!
    assert(frgDat[frgId].sta != 'C');  //  Can't be singleton!

    //  Fragment in a degenerate, nothing further we need to know (or can know).
    if (frgDat[frgId].sta == 'd') {
    }

    //  Status should be either 's' or '?'.  If '?' the fragment was not seen in a contig or
    //  scaffold and so it must be an unplaced surrogate.  Set the mapping fields to generate an
    //  error.
    else if (frgDat[frgId].sta == '?') {
      //frgDat[frgId].con = UINT32_MAX;
      //frgDat[frgId].bgn = 0;
      //frgDat[frgId].end = UINT32_MAX;
      //frgDat[frgId].len = UINT32_MAX;
      //frgDat[frgId].ori = 'z';
      //frgDat[frgId].typ = 'f';
      frgDat[frgId].sta = 'R';
    }

    //  Otherwise, the fragment was placed somewhere, and it inherits the unitig status code.
    else {
      if (frgDat[frgId].sta != 'p')
        fprintf(stderr, "ERROR: fragment "F_U32" has status %c\n", frgId, frgDat[frgId].sta);
      assert(frgDat[frgId].sta == 'p');

      frgDat[frgId].sta = utgDat[utgId].sta;
    }

    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
promoteDegenerates(char *prefix, char *posmap) {
}





//  After all fragments are placed on scaffolds, this will read the frgutg posmap to determine which
//  fragment are unplaced in surrogates.
void
labelFragments(char *prefix) {
  FILE   *file    = openFile(prefix, "frgutg");
  uint32  lineMax = 1048576;
  char   *line    = new char [lineMax];

  fgets(line, lineMax, file);
  chomp(line);
 
  while (!feof(file)) {
    splitToWords W(line);

    string   frgName(W[0]);
    string   utgName(W[1]);

    assert(frgNames.count(frgName) > 0);
    assert(utgNames.count(utgName) > 0);

    uint32   frgId = frgNames[frgName];
    uint32   utgId = utgNames[utgName];



    fgets(line, lineMax, file);
    chomp(line);
  }

  fclose(file);
  delete [] line;
}



void
loadPosMap(char *prefix, char *gkpName) {

  //  Initialize

  for (uint32 i=0; i<NUM_FRG_LABELS; i++)
    frgLabelToInt[frgLabel[i]] = i;

  //  Read library information.

  fprintf(stderr, "Reading library information.\n");
  {
    readLibraries(prefix);
  }

  //  Read fragment information.

  fprintf(stderr, "Reading fragment information.\n");
  {
    readFrags(prefix);
  }

  //  Read mate information.

  fprintf(stderr, "Reading mate information.\n");
  {
    readMates(prefix);
  }

  //  Read scaffold information.

  fprintf(stderr, "Reading assembly information (scf, ctg, utg).\n");
  {
    createScf(prefix, "deglen", 'd');  //  degen scaffolds (promoted from degen contigs)
    createScf(prefix, "scflen", 'p');  //  real scaffolds

    createCtg(prefix, "deglen", 'd');  //  degen contigs (yup, same file as the degen scaffold)
    createCtg(prefix, "ctglen", 'p');  //  real contigs

    createUtg(prefix, "utglen", '?');  //  real and degen unitigs; type set during utgctg and/or utgscf
  }

  fprintf(stderr, "Reading contig/unitig/fragment -> scaffold information.\n");
  {
    loadCtgScf(prefix, "ctgscf");

    loadUtgScf(prefix, "utgscf");
    loadUtgCtg(prefix, "utgctg");

    loadUtgDeg(prefix, "utgdeg");  //  Implicitly creates ctgscf for the degenerate contigs
    loadFrgDeg(prefix, "frgdeg");  //  MUST be before frgctg, frgutg

    loadFrgScf(prefix, "frgscf");
    loadFrgCtg(prefix, "frgctg");
    loadFrgUtg(prefix, "frgutg");  //  MUST be last.
  }

  //  Surrogates are also a special case.  These unitigs have multiple mappings to contigs and scaffolds.
  //  They are not handled.




  //  Fill out the contig layout in each scaffold.  For each contig, annotate the scaffold it is in
  //  with its location.  Then, for each scaffold, invert that list into a list of gaps.

  {
    for (uint32 ci=0; ci<ctgDat.size(); ci++) {
      if (ctgDat[ci].sta != 'p')
        continue;

      uint32 si = ctgDat[ci].scf.con;

      if (scfDat[si].base == NULL) {
        scfDat[si].base = new char [scfDat[si].len];

        memset(scfDat[si].base, 0, sizeof(char) * scfDat[si].len);
      }

      for (uint32 i=ctgDat[ci].scf.bgn; i<ctgDat[ci].scf.end; i++) {
        assert(ctgDat[ci].scf.bgn <  scfDat[si].len);
        assert(ctgDat[ci].scf.end <= scfDat[si].len);

        scfDat[si].base[i] = 1;
      }
    }
  }

  fprintf(stderr, "Loaded "F_SIZE_T" scaffolds.\n",   scfNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" contigs.\n",     ctgNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" unitigs.\n",     utgNames.size());
  fprintf(stderr, "Loaded "F_SIZE_T" fragments.\n",   frgNames.size());

  {
    uint32   fD=0, fC=0, fd=0, fu=0, fr=0, fs=0, fR=0, fp=0, fS=0, fo=0;

    for (uint32 i=0; i<frgDat.size(); i++)
      if      (frgDat[i].sta == 'D')  fD++;
      else if (frgDat[i].sta == 'C')  fC++;
      else if (frgDat[i].sta == 'd')  fd++;
      else if (frgDat[i].sta == 'u')  fu++;
      else if (frgDat[i].sta == 'r')  fr++;
      else if (frgDat[i].sta == 's')  fs++;
      else if (frgDat[i].sta == 'R')  fR++;
      else if (frgDat[i].sta == 'p')  fp++;
      else if (frgDat[i].sta == 'S')  fS++;
      else if (frgDat[i].sta == 'o')  fo++;
      else {
        fprintf(stderr, "frag "F_U32" status %c\n", i, frgDat[i].sta);
        assert(0);
      }

    fprintf(stderr, "  deleted                      "F_U32"\n", fD);
    fprintf(stderr, "  singleton chaff              "F_U32"\n", fC);
    fprintf(stderr, "  in degenerate contig         "F_U32"\n", fd);
    fprintf(stderr, "  in unique unitig             "F_U32"\n", fu);
    fprintf(stderr, "  in rock unitig               "F_U32"\n", fr);
    fprintf(stderr, "  in stone unitig              "F_U32"\n", fs);
    fprintf(stderr, "  in stone unitig, unresolved  "F_U32"\n", fR);
    fprintf(stderr, "  in pebble unitig             "F_U32"\n", fp);
    fprintf(stderr, "  in singleton unitig          "F_U32"\n", fS);
    fprintf(stderr, "  in other unitig              "F_U32"\n", fo);
  }
}

