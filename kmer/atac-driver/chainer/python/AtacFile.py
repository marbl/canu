#!/usr/bin/env python
# Looking in /usr/local/ir/bin on the Compaqs for the correct Python interpreter.
# export PYTHONPATH=${PYTHONPATH}:$WORK/cds/IR/COMPASS/src/AtacPipeline

"""
Extensive documentation for the Python language is available at
http://www.python.org.
"""

import os, sys, time, getopt, tempfile
import MyFile
import MatchRecord

class AtacFile:

    # The data flow is a pipeline augmented by two read-only indexed
    # FASTA files.
    
    def __init__( self, runName):
        "You must supply a atac file called runName.atac."
        self.runName = runName
        self.comments = []
        self.metacommands = []
        self.globals = {} 
        self.tableformat = {}
        self.tabledata = {}
        self.matches = MyFile.myfile()
        self.runs    = MyFile.myfile()
        
        fp = open(runName,"r")
        for line in fp:
            self.atac_file_parse_line(line)

    def atac_file_parse_line( self, line):
        line = line.strip()
        if(not line):
            return
        # end if
        linetype = line[0]
        if(linetype == '#'):
            # Just a comment: squirrel away or ignore
            self.comments.append(line)
            return
        elif(linetype == '!'):
            self.metacommands.append(line)
            return
        elif(linetype == '/'):
            # Add to the globals dictionary
            (key,value) = line[1:].split('=')
            self.globals[key] = value.strip()
            return
        elif(linetype == '@'):
            list = line[1:].split()
            name = list[0]
            self.tableformat[name] = list[1:]
            self.tabledata[name] = [] # an empty list
            return
        elif(linetype == 'M'):
            fields = line.split()
            if(fields[1] == 'r'):
                print >>self.runs, line
            else:
                print >>self.matches, line
            return
        elif(line == ''):
            pass
        else:
            print >>sys.stderr, "The offending line:"
            print >>sys.stderr, line
            assert(0)
        # end if
    # end def

    def checkpoint(self, filename):
        self.globals["modificationDate"] = time.asctime()
        fp = open(filename,"w")
        for line in self.metacommands:
            print >>fp, line
        for line in self.comments:
            print >>fp, line
        # Output the globals in lexigraphical order.
        list = []
        for key in self.globals:
            list.append("/" + key + "=" + str(self.globals[key]))
        list.sort()
        for line in list:
            print >>fp, line
        self.matches.seek(0)
        for line in self.matches:
            fp.write(line)
        self.matches.seek(0)
        self.runs.seek(0)
        for line in self.runs:
            fp.write(line)
        self.runs.seek(0)
        fp.close()
