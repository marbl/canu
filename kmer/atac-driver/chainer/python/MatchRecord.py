#!/usr/bin/env python

import sys, os, copy, string, tempfile

class AtacRow:
    "A general ATAC row object"
    def __init__(self,line):
        self.kind = line[0]
        splitline = line[1:].split(">")
        self.fixed = splitline[0].split()
        if(len(splitline)>1):
            self.defline = splitline[1]
        else:
            self.defline = ""
        # end if
    # end if
# end class

class MatchRecord:
    """Class representing an exact match."""

    #def __init__ (self, *args):
    #print args
    #    if(args):
    # (line,) = args
    # # print " parse line= " + line
    # self.fromString(line)
    # # end if
    # end def

    def sameAs(self,other):
        return (
            (self.x_orientation == other.x_orientation) and
            (self.x_scaf_uid == other.x_scaf_uid) and
            (self.x_start == other.x_start) and
            (self.y_orientation == other.y_orientation) and
            (self.y_scaf_uid == other.y_scaf_uid) and
            (self.y_start == other.y_start) and
            (self.x_length == other.x_length) and
            (self.y_length == other.y_length)
            )
    # end def

    def isInsideBox(self, one, two):
        # We need to modify this because the matches are not points.
        dxone = self.x_start - one.x_start
        dxtwo = self.x_start - two.x_start
        dyone = self.y_start - one.y_start
        dytwo = self.y_start - two.y_start
        flag = (
            # (self.x_orientation == one.x_orientation) and
            # (self.x_orientation == two.x_orientation) and
            # (self.y_orientation == one.y_orientation) and
            # (self.y_orientation == two.y_orientation) and
            (self.x_scaf_uid == one.x_scaf_uid) and
            (self.x_scaf_uid == two.x_scaf_uid) and
            (self.y_scaf_uid == one.y_scaf_uid) and
            (self.y_scaf_uid == two.y_scaf_uid) and
            (((dxone > 0) and ( dxtwo < 0)) or ((dxone < 0) and ( dxtwo > 0))) and
            (((dyone > 0) and ( dytwo < 0)) or ((dyone < 0) and ( dytwo > 0)))
            )
        return flag
    # end def

    def inSameRunAs(self,x):
        return self.runid == x.runid # same parent
    # end def
    
    def copy(self):
        other = copy.copy(self)
        return other

    def convertFromAtacMatchFormat(self,line):
        fields = line.split()
        if(line[0] == 'M'):
            self.rowtype = fields[0]
            self.subtype = fields[1]
            self.matchid = fields[2]
            self.runid = fields[3]
            self.x_scaf_uid = fields[4]
            self.x_start = int(fields[5])
            self.x_length = int(fields[6])
            self.x_orientation = int(fields[7])
            self.y_scaf_uid = fields[8]
            self.y_start = int(fields[9])
            self.y_length = int(fields[10])
            self.y_orientation = int(fields[11])
            #self.mismatches = int(fields[12])
        elif(line[0] == '-'):
            orientation = (fields[0][1:]=='f')
            self.rowtype = 'M'
            self.subtype = 'x'
            self.matchid = '.' # "BMX"+str(lineCount)
            self.runid = '.'
            self.x_scaf_uid = assemblyId1 + ":" + fields[2]
            self.x_start = int(fields[3])
            self.x_length = int(fields[4])
            self.x_orientation = 1
            self.y_scaf_uid = assemblyId2 + ":" + fields[6]
            self.y_start = int(fields[7])
            self.y_length = int(fields[8])
            self.y_orientation = (-1,1)[orientation] # A cvm variant (flag ? x : y).
    # end def

    def __init__ (self, line, *args):
        #print args
        #if(args):
        #   (line,) = args
        if '>' in line:
            (line1, line2) = line.split('>')
        else:
            line1 = line
            line2 = ""
        try:
            self.convertFromAtacMatchFormat(line1)
        except IndexError:
            sys.stderr.write("MatchRecord-- IndexError: line did not split correctly: %s\n"
                             % line1)
            raise
        except ValueError:
            sys.stderr.write("MatchRecord-- ValueError: line did not unpack correctly: %s\n"
                             % line1)
            raise
        self.extend = {}
        extensions = line2.split('/')
        self.identifier = extensions[0].strip()
        for argpair in extensions[1:]:
            if '=' in argpair:
                (key,value) = argpair.split('=')
                self.extend[key] = value.strip()
        return

    def __str__ (self):
        extension = " >" + self.identifier
        for key in self.extend:
            extension += ' /' + key + '=' + str(self.extend[key])
        if(len(extension)<3):
            extension = ""
        return "%s %s %s %s %s %d %d %d %s %d %d %d %s" % (
            self.rowtype, self.subtype, self.matchid, self.runid,
            self.x_scaf_uid, self.x_start, self.x_length, self.x_orientation,
            self.y_scaf_uid, self.y_start, self.y_length, self.y_orientation,
            extension
            )
    # end def

# end class


def convertBrianRecordFormat( inpfile, outfile, assemblyId1, assemblyId2):
    "Convert the match record format from Brian's to atac format."
    lineCount = 0
    for line in inpfile:
        lineCount += 1
        if(lineCount % 100000 == 0):
            print >>sys.stderr, "lineCount=%d" % lineCount
        FB = line.split()
        orientation = (FB[0][1:]=='f')
        FM = MatchRecord("M x . . . 0 0 0 . 0 0 0 0\n")
        FM.x_orientation = 1
        FM.matchid       = "BMX"+str(lineCount)
        FM.x_scaf_uid    = assemblyId1 + ":" + FB[2]
        FM.x_start       = int(FB[3])
        FM.x_length      = int(FB[4])
        FM.y_orientation = (-1,1)[orientation] # A cvm variant (flag ? x : y).
        FM.y_scaf_uid    = assemblyId2 + ":" + FB[6]
        FM.y_start       = int(FB[7])
        FM.y_length      = int(FB[8])
        FM.identifier    = ""
        FM.extend        = {}
        #FM.mismatches   = 0
        print >>outfile, FM
    # end for
    print >>sys.stderr, "convertRecordFormat done: lineCount=%d" % lineCount
    outfile.seek(0)
    return
# end def


def sortInXorderAP( inpfile, outfile):
    # (x_scaf_uid, x_start, x_length, y_scaf_uid, y_start, y_length)
    InXOrderAP = '-k 1,1 -k 2,2 -k 5,5 -k 6n -k 7nr -k 8nr -k 9,9 -k 10n -k 11nr -k 12nr'
    # Use -u to remove the palindromes.
    # Use -k 7nr -k 11nr to remove abutting contained matches.

    if (not os.path.exists('./tmp')):
        os.system('mkdir ./tmp')
    # end if

    inpfile.seek(0)
    outfile.seek(0)
    inpfile.flush()
    outfile.flush()

    ierr = os.system("sync;sync;sync")
    assert(ierr == 0)
    ierr = os.system("sort -T ./tmp %s %s > %s" %
                     (InXOrderAP, inpfile.name, outfile.name));
    assert(ierr == 0)
    ierr = os.system("sync;sync;sync")
    assert(ierr == 0)
    inpfile.seek(0)
    outfile.seek(0)

    return
# end def

def sortInYorderAP( inpfile, outfile):
    # (x_scaf_uid, x_start, x_length, y_scaf_uid, y_start, y_length)
    InYOrderAP = '-k 1,1 -k 2,2 -k 9,9 -k 10n -k 11nr -k 12nr -k 5,5 -k 6n -k 7nr -k 8nr'
    # Use -u to remove the palindromes.
    # Use -k 7nr -k 11nr to remove abutting contained matches.

    if (not os.path.exists('./tmp')):
        os.system('mkdir ./tmp')
    # end if

    inpfile.seek(0)
    outfile.seek(0)
    inpfile.flush()
    outfile.flush()

    ierr = os.system("sync;sync;sync")
    assert(ierr == 0)
    ierr = os.system("sort -T ./tmp %s %s > %s" %
                     (InYOrderAP, inpfile.name, outfile.name));
    assert(ierr == 0)
    ierr = os.system("sync;sync;sync")
    assert(ierr == 0)
    inpfile.seek(0)
    outfile.seek(0)

    return
# end def

def sortInXorderPP( inpname, outfile):
    # (x_win, ywin, x_scaf_uid, y_scaf_uid, x_start, y_start, x_length, y_length)
    assert(1)
# end def

def sortInYorderPP( inpname, outfile):
    # (y_win, x_win, y_scaf_uid, x_scaf_uid, y_start, x_start, y_length, x_length)
    assert(1)
# end def
