#!/usr/bin/env python

import os, sys, time, tempfile
import MyFile
import MatchRecord
import AtacFile

def xorIntervals( inpname, outname):
    # not tested yet
    leftPicket = 0
    rghtPicket = 0
    inpfile = open(inpname,"r")
    outfile = open(outname,"w")
    for line in inpfile:
        fields = line.split()
        newstart = int(fields[0])
        newend   = int(fields[1])
        assert(leftPicket <= newstart)
        rghtSide = min(newstart,rghtPicket)
        if rghtSide > leftPicket:
            # interval has positive length
            print >>outfile, leftPicket, rghtSide
        leftPicket = max(newstart,  min(rightPicket,newend))
        rghtPicket = max(leftPicket,max(rightPicket,newend))


def findUniformCoverageIntervals(inpfile,outfile):
    # The input records are ("E", id, position, coverage_increment).
    # The output records are ("C", id, start_position, length, coverage_level).
    inpfile.seek(0)
    outfile.seek(0)
    oldaxis = None; oldposition = 0; cov = 0
    for line in inpfile:
        ( recordtype, newaxis, newposition, newchange) = line.split()
        if(recordtype == "E"):
            newposition = int(newposition)
            newchange = int(newchange)
            if(newaxis != oldaxis and cov != 0): print >>sys.stderr, "Woops"
            len = newposition - oldposition
            if(cov>0 and len>0):
                print >>outfile, "C", oldaxis, oldposition, len, cov;
            cov += newchange;
            assert(cov >= 0)
            oldaxis = newaxis; oldposition = newposition;
    assert(cov == 0)
    outfile.flush()


def findCoverageIntervals( inpfile, outfile, processFirstAxis):
    # The input file is an ATAC matches file.
    # The output file is an ATAC coverage intervals file.
    inpfile.seek(0)
    outfile.seek(0)
    t0 = time.time()
    tmpfile3 = MyFile.myfile()
    for line in inpfile:
        if(line[0]=="M"):
            fields = line.split()
            if(fields[1]=="u" or fields[1]=="x"):
                if(processFirstAxis):
                    axis = fields[4]
                    bgn = int(fields[5])
                    end = bgn+int(fields[6])
                else:
                    axis = fields[8]
                    bgn = int(fields[9])
                    end = bgn+int(fields[10])
                print >>tmpfile3, "E", axis,bgn,1
                print >>tmpfile3, "E", axis,end,-1
    tmpfile3.close()
    tmpname = tempfile.mktemp()
    cmd = "sort -T . -k 1,1 -k 2,2 -k 3n -k 4nr  %s > %s" % (tmpfile3.name, tmpname)
    print >>sys.stderr, cmd
    iret = os.system(cmd); assert(iret==0)
    print >>sys.stderr,"time elapsed is ", (time.time() - t0)
    tmpfile4 = open(tmpname)
    t0 = time.time()
    findUniformCoverageIntervals( tmpfile4, outfile)
    print >>sys.stderr,"time elapsed is ", (time.time() - t0)
    tmpfile4.close()
    os.system("rm -f " + tmpname)
    outfile.seek(0)
    

def applyOneKeepMask( inpfile, outfile, keepMaskFile, processFirstAxis):
    # Note that the following merge-like control structure is
    # influenced by the function property of keep intevals to matches.

    debug = 0
    inpfile.seek(0)
    outfile.seek(0)
    keepMaskFile.seek(0)

    # Put the first valid match record into FM.  Each input ATAC match
    # record produces zero, one or more output ATAC matches.
    FM = None;
    ma = None; ms = None; me = None
    qa = None; qs = None; ql = None;

    # the set of masking intervals, using the q variables and iline
    maskiter = iter(keepMaskFile)

    # the set of masked matches using the m variables and mline
    inpiter = iter(inpfile)

    iline = None
    mline = None

    last_matchid = None; subcount = 0
    
    try:   # StopIteration exception from either iterator gets us out
        while 1:
            if(iline == None):
                iline = maskiter.next()
                (subtype, qa, qs, ql, cov, ) = iline.split()
                assert(subtype=='C')
                cov= int(cov)
                if(cov != 1):
                    iline = None
                    continue
                qs = int(qs)
                ql = int(ql)
                qe = qs + ql

            if(mline == None):
                mline = inpiter.next()
                if(mline[0] != 'M'):
                    # not a match record, so just pass it through
                    print >>outfile, mline,
                    mline = None
                    continue
                FM = MatchRecord.MatchRecord(mline)
                assert(FM.subtype == "u" or FM.subtype == "x")
                if(processFirstAxis):
                    ma = FM.x_scaf_uid
                    ms = FM.x_start       # match start
                    me = ms + FM.x_length # match end
                else:
                    ma = FM.y_scaf_uid
                    ms = FM.y_start       # match start
                    me = ms + FM.y_length # match end

            # holding valid iline and mline data now
        
            if not(ma==qa):
                # not on same axis, need to get a new one
                if(ma < qa):
                    mline = None
                else:
                    iline = None

            elif not( (ms < qe) and (qs < me) ):
                # we are not overlapping, need to get a new one of them
                if(ms < qs):
                    mline = None
                else:
                    iline = None

            else:
                # processing for overlaps
                FT = FM.copy()
                mx = max(ms,qs)
                mn = min(me,qe)
                trimFromStart = mx - ms
                trimFromEnd =   me - mn
                trimmedLength = mn - mx
                if( FT.x_orientation == FT.y_orientation):
                    FT.x_start += trimFromStart
                    FT.y_start += trimFromStart
                else:
                    if(processFirstAxis):
                        FT.x_start += trimFromStart
                        FT.y_start += trimFromEnd
                    else:
                        FT.y_start += trimFromStart
                        FT.x_start += trimFromEnd
                FT.x_length = trimmedLength
                FT.y_length = trimmedLength
                if debug:
                    print >>sys.stdout, "# trimmed "
                    print >>sys.stdout, FT

                # We must insure that the match identifier is still unique.
                if last_matchid == FM.matchid :
                    subcount += 1
                else:
                    subcount = 0
                # print >>sys.stderr, last_matchid, FM.matchid, subcount
                last_matchid = FM.matchid

                if(subcount > 0):
                    if processFirstAxis :
                        FT.matchid = FT.matchid + "x" + str(subcount)
                    else:
                        FT.matchid = FT.matchid + "y" + str(subcount)
                
                print >>outfile, FT
                # we need to get a new one
                if(qe < me):
                    iline = None
                else:
                    mline = None

    except StopIteration:
        # If there are any left over non-match lines, then output them!
        for mline in inpiter:
            if(mline[0] != "M"):
                print >>outfile, mline,

def applyBothKeepMasks( inpfile, outfile ):

    # Maybe we can think of a masking implementation where each ATAC match
    # is treated atomicly.  Assume that the keep mask intervals are sorted
    # by start postition.  Assume that the ATAC matches are sorted by start
    # postion.  Assert that all keep mask intervals are non-overlapping and
    # were cut from only one ATAC match.  Thus the mapping from keep mask
    # intervals is a function.  Note that this requires that we do not
    # coalesce abutting keep mask intervals that originate from multiple
    # matches.  Note this still allows an ATAC match to overlap more than
    # one keep mask interval.  Ignore all keep mask intervals with zero
    # length their creation has tie breaking problems.  See notes on 2003
    # Jul 29.

    debug = 0
    debugnum = 0
    inpfile.seek(0)
    outfile.seek(0)


    # Apply the keepMask for the first axis.
    # Make the sorted the keep mask intervals for the first axis.
    processFirstAxis = 1
    keepMaskFile = MyFile.myfile()
    tmpfile2 = inpfile
    tmpfile3 = MyFile.myfile()
    tmpfile4 = MyFile.myfile()

    findCoverageIntervals( inpfile, keepMaskFile, processFirstAxis)
    if debug:
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in keepMaskFile: print >>debugfile, line,
            
    MatchRecord.sortInXorderAP(tmpfile2,tmpfile3)
    if debug:
        #tmpfile2.seek(0)
        #debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        #for line in tmpfile2: print >>debugfile, line,
        tmpfile3.seek(0)
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in tmpfile3: print >>debugfile, line,
        
    applyOneKeepMask( tmpfile3, tmpfile4, keepMaskFile, processFirstAxis)
    if debug:
        tmpfile4.seek(0)
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in tmpfile4: print >>debugfile, line,
        
    # Apply the keepMask for the second axis.
    # Make the sorted the keep mask intervals for the second axis.
    processFirstAxis = 0
    keepMaskFile = MyFile.myfile()
    tmpfile2 = tmpfile4
    tmpfile3 = MyFile.myfile()
    tmpfile4 = outfile

    findCoverageIntervals( inpfile, keepMaskFile, processFirstAxis)
    if debug:
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in keepMaskFile: print >>debugfile, line,


    MatchRecord.sortInYorderAP(tmpfile2,tmpfile3)
    if debug:
        #tmpfile2.seek(0)
        #debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        #for line in tmpfile2: print >>debugfile, line,
        tmpfile3.seek(0)
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in tmpfile3: print >>debugfile, line,

    applyOneKeepMask( tmpfile3, tmpfile4, keepMaskFile, processFirstAxis)
    if debug:
        tmpfile4.seek(0)
        debugnum += 1; debugfile = open("debugfile.%d" % debugnum, "w")
        for line in tmpfile4: print >>debugfile, line,


def main( inpfile, outfile):
    applyBothKeepMasks( inpfile, outfile)

    # Should we check if the first and last characters of the masked
    # matches are matching?

    # Should we compute the percent identity in this module?



# Allow each module to have its own main for testing.
if __name__ == '__main__':
    inpname = sys.argv[1]
    outname = sys.argv[2]
    inpfile = open(inpname)
    outfile = open(outname,"w")
    main(inpfile, outfile)
# end if
