#!/usr/bin/env python
# Looking in /usr/local/ir/bin on the Compaqs for the correct Python interpreter.
# export PYTHONPATH=${PYTHONPATH}:$WORK/cds/IR/COMPASS/src/AtacPipeline

"""
Extensive documentation for the Python language is available at
http://www.python.org.
"""

import sys
import MyFile
import MatchRecord

def cvm(f,x,y):
    # A cvm variant (flag ? x : y) = (y,x)[f]
    if f :
        return x
    else:
        return y
    # end if
# end def

def createSignedEnumeration(inpfile):
    outfile = MyFile.myfile()
    p = 1
    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            FM = MatchRecord.MatchRecord(line)
            forwardX = FM.x_orientation
            forwardY = FM.y_orientation
            srank = cvm(forwardX == forwardY, p, -p)
            p += 1
            FM.extend['srank'] = srank
            print >>outfile, FM
        # end if
    # end while
    return outfile
# end def

def findPerfectRuns ( inpfile, maxJump, runIdPrefix ):
    outfile = MyFile.myfile()
    left = None
    runid = 1
    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            right = MatchRecord.MatchRecord(line)
            pr = int(right.extend['srank'])
            del(right.extend['srank'])
            if(left != None):
                maxGapInXandY = 0
                if(left.x_scaf_uid == right.x_scaf_uid
                   and
                   left.y_scaf_uid == right.y_scaf_uid ):
                    # Find the maximum of the gap in x and y axis.
                    
                    x_rs = right.x_start
                    x_re = x_rs + right.x_length
                    x_ls = left.x_start
                    x_le = x_ls + left.x_length
                    assert(x_rs < x_re)
                    assert(x_ls < x_le)
                    # All matches are positive length.
                    x_gapLeftBeforeRight = x_rs - x_le
                    x_gapRightBeforeLeft = x_ls - x_re
                    assert(not(x_gapLeftBeforeRight>0 and x_gapRightBeforeLeft>0))
                    x_gap = max(x_gapLeftBeforeRight,x_gapRightBeforeLeft)
                    # x_gap == 0 is abutting
                    # x_gap < 0  is overlapping

                    y_rs = right.y_start
                    y_re = y_rs + right.y_length
                    y_ls = left.y_start
                    y_le = y_ls + left.y_length
                    assert(y_rs < y_re)
                    assert(y_ls < y_le)
                    y_gapLeftBeforeRight = y_rs - y_le
                    y_gapRightBeforeLeft = y_ls - y_re
                    assert(not(y_gapLeftBeforeRight>0 and y_gapRightBeforeLeft>0))
                    y_gap = max(y_gapLeftBeforeRight,y_gapRightBeforeLeft)
                    # y_gap == 0 is abutting
                    # y_gap < 0  is overlapping

                    maxGapInXandY = max(x_gap,y_gap)

                    if 1:
                        # Check the sorting of the matches.
                        
                        sorted_by_x = (x_ls <= x_rs)
                        sorted_by_y = (y_ls <= y_rs)

                        if(not(sorted_by_x or sorted_by_y)):
                            print >>sys.stderr, "bad sorting in findPerfectRuns"
                            print >>sys.stderr, left
                            print >>sys.stderr, right
                        assert(sorted_by_x or sorted_by_y)
                        dovetail_in_x = (x_ls <= x_rs) and (x_le <= x_re)
                        dovetail_in_y = (y_ls <= y_rs) and (y_ls <= y_re)
                        if(sorted_by_x and not(dovetail_in_x)):
                            print >>sys.stderr, "contained in x in findPerfectRuns"
                            print >>sys.stderr, left
                            print >>sys.stderr, right
                        if(sorted_by_y and not(dovetail_in_y)):
                            print >>sys.stderr, "contained in y in findPerfectRuns"
                            print >>sys.stderr, left
                            print >>sys.stderr, right
                # endif
                if(
                    (left.x_scaf_uid != right.x_scaf_uid) or  # check first axis id
                    (left.y_scaf_uid != right.y_scaf_uid) or  # check second axis id
                    (maxGapInXandY > maxJump) or
                    (pr != lastpr + 1)  # Using the signed rank NOT the run id !!!!
                    ):
                    runid += 1
                # end if
            # end if
            lastpr = pr
            right.runid = "%s%d" % (runIdPrefix,runid,)   # Assign the run id in the same slot as the signed rank.
            print >>outfile, right
            left = right
        # end if
    # end for
    return outfile
# end def

def formPerfectRuns ( inpfile, firstSort, secondSort, maxJump, runIdPrefix ):
    inpfile.seek(0)
    step = 0
    print >>sys.stderr, 'formPerfectRuns step=' + str(step)
    step += 1

    tmpfile = MyFile.myfile()
    firstSort( inpfile, tmpfile)
    
    print >>sys.stderr, 'formPerfectRuns step=' + str(step)
    step += 1
    outfile = createSignedEnumeration(tmpfile)
    
    print >>sys.stderr, 'formPerfectRuns step=' + str(step)
    step += 1
    tmpfile = MyFile.myfile()
    secondSort( outfile, tmpfile)
    
    print >>sys.stderr, 'formPerfectRuns step=' + str(step)
    step += 1
    outfile = findPerfectRuns( tmpfile, maxJump, runIdPrefix)
    
    return outfile
# end def

def runsAsMatches(inpfile):

    outfile = MyFile.myfile()
    lastF = None
    firstF = None
    runFill = 0
    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            curF = MatchRecord.MatchRecord(line)
            if ((lastF == None) or (curF.runid != lastF.runid)):
                if ((lastF != None) and (firstF.x_scaf_uid != lastF.x_scaf_uid)):
                    print >>sys.stderr, firstF
                    print >>sys.stderr, lastF
                # end if
                assert((lastF==None) or (firstF.x_scaf_uid == lastF.x_scaf_uid))
                assert((lastF==None) or (firstF.y_scaf_uid == lastF.y_scaf_uid))
                if (None != lastF):
                    x1 = firstF.x_start
                    x2 = lastF.x_start
                    startX = cvm(x1 < x2, x1, x2)
                    x1 += firstF.x_length
                    x2 += lastF.x_length
                    endX = cvm(x1 > x2, x1, x2)
                    y1 = firstF.y_start
                    y2 = lastF.y_start
                    startY = cvm( y1 < y2, y1, y2)
                    y1 += firstF.y_length
                    y2 += lastF.y_length
                    endY = cvm(y1 > y2, y1, y2)
                    lastF.subtype = 'r'
                    lastF.matchid = lastF.runid
                    lastF.runid = "."   # the agreed NULL value
                    lastF.x_start = startX
                    lastF.y_start = startY
                    lastF.x_length = endX - startX
                    lastF.y_length = endY - startY
                    lastF.runFill = runFill
                    print >>outfile, lastF
                # end if
                firstF = curF
                runFill = 0
            # end if
            runFill += curF.x_length
            lastF = curF
        # end if
    # end for
    
    if (None != lastF):
        x1 = firstF.x_start
        x2 = lastF.x_start
        startX = cvm( x1 < x2, x1, x2)
        x1 += firstF.x_length
        x2 += lastF.x_length
        endX = cvm( x1 > x2, x1, x2)
        y1 = firstF.y_start
        y2 = lastF.y_start
        startY = cvm( y1 < y2, y1, y2)
        y1 += firstF.y_length
        y2 += lastF.y_length
        endY = cvm( y1 > y2, y1, y2)
        lastF.subtype =  'r'
        lastF.matchid = lastF.runid
        lastF.runid = "."  # the agreed NULL value
        lastF.x_start = startX
        lastF.y_start = startY
        lastF.x_length = endX - startX
        lastF.y_length = endY - startY
        lastF.runFill = runFill
        print >>outfile, lastF
    # end if
    return outfile
# end def

def main(inpname, outname, maxJump, runIdPrefix):
    inpfile = open(inpname)
    tempdata1 = formPerfectRuns(inpfile,
                               MatchRecord.sortInXorderAP,
                               MatchRecord.sortInYorderAP,
                               int(maxJump),
                               runIdPrefix
                               )
    tempdata2 = runsAsMatches( tempdata1)
    tempdata1.link(outname+".matches")
    tempdata2.link(outname+".runs")

def usage():
    print >>sys.stderr, "No usage statement yet"
    print >>sys.stderr, "prog -i inpname -o outname [-j int] [-p string]"
    
def mainTry():
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:vj:p:")
        # , ["help", "inpname", "outname=", "maxJump=", "runIdPrefix="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)
    verbose = False
    inpname = None
    outname = None
    maxJump = 100000
    runIdPrefix = "r"
    for o, a in opts:
        if o == "-v":
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--outname"):
            outname = a
        if o in ("-i", "--inpname"):
            inpname = a
        if o in ("-j","--maxJump",):
            try:
                maxJump = int(a)
            except:
                print >>sys.stderr, "maxJump must be an integer"
                sys.exit(3)
        if o in ("p", "--runIdPrefix",):
            runIdPrefix = a
    # ...
    if(inpname != None and outname != None):
        usage()
        sys.exit(2)

    main(inpname, outname, maxJump, runIdPrefix)

if __name__ == '__main__':
    inpname = sys.argv[1]
    outname = sys.argv[2]
    maxJump = int(sys.argv[3])
    runIdPrefix = sys.argv[4]

    main(inpname, outname, maxJump, runIdPrefix)
