#!/usr/bin/env python
# Must look in /usr/local/ir/bin on the Compaqs for the correct Python interpreter.
# export PYTHONPATH=${PYTHONPATH}:$WORK/cds/IR/COMPASS/src/AtacPipeline

"""
Extensive documentation for the Python language is available at
http://www.python.org.
"""

import os, sys, re, tempfile

def main(glist):
    for inpname in glist:

        if(0):
            inpfile = open(inpname,'r')
            tmpname = tempfile.mktemp(".tmp")
            tmpfile = open(tmpname,'w')
            pattern = re.compile(r"^M [gl] ")
            for line in inpfile:
                if(pattern.search(line)):
                    print >>tmpfile, line,
            tmpfile.close()
            os.system("celagram -c 7 -t 'gapped match lengths' %s" % (tmpname,))

        if(0):
            inpfile = open(inpname,'r')
            tmpname = tempfile.mktemp(".tmp")
            tmpfile = open(tmpname,'w')
            pattern = re.compile(r"^M x ")
            for line in inpfile:
                if(pattern.search(line)):
                    print >>tmpfile, line,
            tmpfile.close()
            os.system("celagram -c 7 -t 'exact match lengths' %s" % (tmpname,))

        if(0):
            inpfile = open(inpname,'r')
            tmpname = tempfile.mktemp(".tmp")
            tmpfile = open(tmpname,'w')
            pattern = re.compile(r"^M u ")
            for line in inpfile:
                if(pattern.search(line)):
                    print >>tmpfile, line,
            tmpfile.close()
            os.system("celagram -c 7 -t 'ungapped match lengths' %s" % (tmpname,))

        inpfile = open(inpname,'r')
        tmpname = tempfile.mktemp(".tmp")
        tmpfile = open(tmpname,'w')
        # pattern = re.compile(r"^M\s*[xu]\s")
        pattern = re.compile(r"^M [xu] ")
        for line in inpfile:
            if(pattern.search(line)):
                print >>tmpfile, line,
        tmpfile.close()
        os.system("celagram -c 7 -t '%s ungapped match lengths' %s" % (inpname,tmpname))


        inpfile = open(inpname,'r')
        tmpname = tempfile.mktemp(".tmp")
        tmpfile = open(tmpname,'w')
        pattern = re.compile(r"^M\s*r\s")
        for line in inpfile:
            if(pattern.search(line)):
                print >>tmpfile, line,
        tmpfile.close()
        os.system("celagram -c 7  -t '%s spans in 1st assembly' %s" % (inpname,tmpname))
        os.system("celagram -c 11 -t '%s spans in 2nd assembly' %s" % (inpname,tmpname))


        
if __name__ == '__main__':
    #glist = [ "humR27vsB31-V2.atac", "humB31vsVAN-V1.atac", "humB31vsSC-V3.atac", ]

    main(sys.argv[1:])
