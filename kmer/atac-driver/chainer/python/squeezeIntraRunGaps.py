#!/usr/bin/env python

import sys
import string
import MatchRecord
#import MyFile
import IdxStore
#import localAlignerInterface
#import shelve

# True=1
False=0

        
theIsolatedSNPcount = 0
completefillednotXY = 0
completefilledXnotY = 0
completefilledYnotX = 0
completefilledXandY = 0

def analyzeGap(x,y,left,right,outfile,maxgap,margin):
    global theIsolatedSNPcount
    global completefillednotXY
    global completefilledXnotY
    global completefilledYnotX
    global completefilledXandY
    solidThreshold=20
    inter_run_gap_count = 0
    x_chCount = {}
    y_chCount = {}
    x_notACGT = 0
    y_notACGT = 0

    lp = 0  # We should modify the match instead!
    rp = 0
    x_pos = 0
    x_len = 0
    y_pos = 0
    y_len = 0
    if( (left.x_scaf_uid == right.x_scaf_uid) and
        (left.y_scaf_uid == right.y_scaf_uid) and
        (left.runid == right.runid) ):
        # (left.sindex + 1 == right.sindex) ): # This is obsolete in Russell's file format.
        # sys.stderr.write("Intra-run gap\n")
        left_forward  = (left.x_orientation  == left.y_orientation)
        right_forward = (right.x_orientation == right.y_orientation)
        if( left_forward != right_forward): sys.stderr.write("Bad orientations in run\n")
        assert(left_forward == right_forward)

        sorted_by_x = (left.x_start <= right.x_start)
        dovetail_in_x = sorted_by_x and (left.x_start+left.x_length <= right.x_start+right.x_length)
        sorted_by_y = (left.y_start <= right.y_start)
        dovetail_in_y = sorted_by_y and (left.y_start+left.y_length <= right.y_start+right.y_length)
        if(not(sorted_by_x or sorted_by_y)):
           print >>sys.stderr, "bad sorting in runs"
           print >>sys.stderr, left
           print >>sys.stderr, right
        assert(sorted_by_x or sorted_by_y)
        # This concept of sorted allows neggaps but not containmant in both axes.

        if(sorted_by_x and not dovetail_in_x):
           print >>sys.stderr, "sorted_by_x and not dovetail_in_x"
           print >>sys.stderr, left
           print >>sys.stderr, right
        
        if(sorted_by_y and not dovetail_in_y):
           print >>sys.stderr, "sorted_by_y and not dovetail_in_y"
           print >>sys.stderr, left
           print >>sys.stderr, right
        
        if(not((not left_forward) or (sorted_by_x and sorted_by_y))):
           print >>sys.stderr, "bad sorting in runs"
           print >>sys.stderr, left
           print >>sys.stderr, right
        if(not((left_forward) or (not(sorted_by_x and sorted_by_y)))):
           print >>sys.stderr, "bad sorting in runs"
           print >>sys.stderr, left
           print >>sys.stderr, right
        assert((not left_forward) or (sorted_by_x and sorted_by_y))
        assert((left_forward) or (not(sorted_by_x and sorted_by_y)))

        if(sorted_by_x): # Sorted by x positions.
            #print "Sorted by X"
            x_pos = left.x_start + left.x_length # Start of the intra-run gap.
            x_len = right.x_start - x_pos        # Length of the intra-run gap.
            if(left_forward):
                y_pos = left.y_start + left.y_length
                y_len = right.y_start - y_pos
            else:
                y_pos = right.y_start + right.y_length
                y_len = left.y_start - y_pos
            # end if
        else: # Assume sorted by y positions
            #print "Sorted by Y"
            y_pos = left.y_start + left.y_length
            y_len = right.y_start - y_pos
            if(left_forward):
                x_pos = left.x_start + left.x_length
                x_len = right.x_start - x_pos
            else:
                x_pos = right.x_start + right.x_length
                x_len = left.x_start - x_pos
            # end if
        # end if
        # print "Left %d,%d   Right %d,%d  Width %d,%d" % (x_pos,y_pos,x_pos+x_len,y_pos+y_len,x_len,y_len)

        assert(left.x_start >= 0)
        assert(left.x_length > 0)
        assert(left.y_start >= 0)
        assert(left.y_length > 0)
        assert(right.x_start >= 0)
        assert(right.x_length > 0)
        assert(right.y_start >= 0)
        assert(right.y_length > 0)

        # Trim the intra-run neggaps to become proper gaps.
        if(0 and (x_len < 0 or y_len < 0)):
            sys.stderr.write("neggap x_uid= %s x_pos= %d x_len= %d y_uid= %s y_pos= %d y_len= %d\n" %
                             (left.x_scaf_uid,x_pos,x_len,left.y_scaf_uid,y_pos,y_len))
            trim_len = max( -x_len,  -y_len)
            # Increase the intra-run gap length:
            x_len += trim_len; y_len += trim_len
            # Decrease the right-hand match length:
            right.x_length -= trim_len
            right.y_length -= trim_len
            # Adjust the right-hand gap ending position:
            if(left_forward):
                right.x_start += trim_len
                right.y_start += trim_len
            else:
                if(sorted_by_x):
                    right.x_start += trim_len
                    y_pos -= trim_len
                else: # assume sorted_by_y
                    right.y_start  += trim_len
                    x_pos -= trim_len
                # end if
            #end if
            sys.stderr.write("newgap x_uid= %s x_pos= %d x_len= %d y_uid= %s y_pos= %d y_len= %d\n" %
                             (left.x_scaf_uid,x_pos,x_len,left.y_scaf_uid,y_pos,y_len))
        # end if
        assert(right.x_length > 0)
        assert(right.y_length > 0)
        # We now have a proper intra-run gap segment between two match segments.

        x_substring = ""
        if(x_len > 0):
            x_substring = string.upper(
                x.getStringFromFasta( sorted_by_x, left.x_scaf_uid, x_pos, x_len));
        # end if
        y_substring = ""
        if(y_len > 0):
            y_substring = string.upper(
                y.getStringFromFasta( sorted_by_y, left.y_scaf_uid, y_pos, y_len));
        # end if
        if(x_len > 0 and not(x_len == len(x_substring))):
            sys.stderr.write("x string lengths mismatch asked=%d got=%d\n" % (x_len,len(x_substring)))
            sys.stderr.write("x_uid= %s x_pos= %d x_len= %d y_uid= %s y_pos= %d y_len= %d\n" %
                                       (left.x_scaf_uid,x_pos,x_len,left.y_scaf_uid,y_pos,y_len))
            print >>sys.stderr, "left match"
            print >>sys.stderr, left
            print >>sys.stderr, "right match"
            print >>sys.stderr, right
        # end if
        if(y_len > 0 and not(y_len == len(y_substring))):
            sys.stderr.write("y string lengths mismatch asked=%d got=%d\n" % (y_len,len(y_substring)))
            sys.stderr.write("x_uid= %s x_pos= %d x_len= %d y_uid= %s y_pos= %d y_len= %d\n" %
                                       (left.x_scaf_uid,x_pos,x_len,left.y_scaf_uid,y_pos,y_len))
            print >>sys.stderr, "left match"
            print >>sys.stderr, left
            print >>sys.stderr, "right match"
            print >>sys.stderr, right
        # end if
        
        assert(x_len < 0 or x_len == len(x_substring))
        assert(y_len < 0 or y_len == len(y_substring))
        assert(lp == 0)
        assert(rp == 0)

        # Next we extend the raw matches to squeeze the intra-run gaps
        # with exactly matching sequence.
        if( lp+rp < x_len and lp+rp < y_len ):
            while(lp+rp < x_len and lp+rp < y_len):
                # modify lp
                x_ch = x_substring[lp]; y_ch = y_substring[lp];
                is_a_match = (x_ch==y_ch) and \
                             (x_ch=="A" or x_ch=="C" or x_ch=="G" or x_ch=="T")
                if(is_a_match):
                    lp += 1
                else:
                    break
                # end if
            # end while
            while(lp+rp < x_len and lp+rp < y_len):
                # modify rp
                x_ch = x_substring[-1-rp]; y_ch = y_substring[-1-rp];
                is_a_match = (x_ch==y_ch) and \
                             (x_ch=="A" or x_ch=="C" or x_ch=="G" or x_ch=="T")
                if(is_a_match):
                    rp += 1
                else:
                    break
                # end if
            # end while
        # end if

        # Next we extend the raw matches to squeeze the intra-run gaps.
        # Each mismatch character must be padded on both sides by
        # "solidThreshold" characters from {A,C,G,T}.
        if( x_len > lp+rp and y_len > lp+rp ):
            lq = lp; solid = solidThreshold; tentativeSNPCount = 0
            while(lq+rp < x_len and lq+rp < y_len):
                x_ch = x_substring[lq]; y_ch = y_substring[lq];
                is_a_match = (x_ch==y_ch) and \
                             (x_ch=="A" or x_ch=="C" or x_ch=="G" or x_ch=="T")
                if(solid >= solidThreshold):
                    lp = lq
                    theIsolatedSNPcount += tentativeSNPCount
                    tentativeSNPCount = 0
                    if(is_a_match):
                        solid += 1
                        lp += 1
                    else:
                        solid = 0
                        tentativeSNPCount = 1
                        tentativeSNPposition = lq
                    # end if
                else:
                    if(is_a_match):
                        solid += 1
                    else: # a second mismatch within 20 bp
                        break
                    # end if
                # end if
                lq += 1
            # end while
            if(lq+rp == x_len and lq+rp == y_len):
                lp = lq
                theIsolatedSNPcount += tentativeSNPCount
                tentativeSNPCount = 0
            # end if
            rq = rp; solid = solidThreshold; tentativeSNPCount = 0
            while( lp+rq < x_len and lp+rq < y_len ):
                x_ch = x_substring[-1-rq]; y_ch = y_substring[-1-rq];
                is_a_match = (x_ch==y_ch) and \
                             (x_ch=="A" or x_ch=="C" or x_ch=="G" or x_ch=="T")
                if(solid >= solidThreshold):
                    rp = rq
                    theIsolatedSNPcount += tentativeSNPCount
                    tentativeSNPCount = 0
                    if(is_a_match):
                        solid += 1
                        rp += 1
                    else:
                        solid = 0
                        tentativeSNPCount = 1
                    # end if
                else:
                    if(is_a_match):
                        solid += 1
                    else: # a second mismatch within 20 bp
                        break
                    # end if
                # end if
                rq += 1
            # end while
            if( lp+rq == x_len and lp+rq == y_len ):
                rp = rq
                theIsolatedSNPcount += tentativeSNPCount
                tentativeSNPCount = 0
            # end if
        # end if
            
        # Next we close any remaining intra-run gaps that can form
        # ungapped alignments of a specified high quality.
        # Currently we have hard coded that there must be
        # 5 or less mismatches or
        # better than 95% identity
        # in the intrarun gap remaining after the previous gap closing.
        assert(lp >= 0)
        assert(rp >= 0)
        if( x_len == y_len and x_len > lp+rp):
            lq = lp; mismatchCount = 0
            while(lq+rp < x_len and lq+rp < y_len):
                x_ch = x_substring[lq]; y_ch = y_substring[lq];
                is_a_match = (x_ch==y_ch) and \
                             (x_ch=="A" or x_ch=="C" or x_ch=="G" or x_ch=="T")
                if(not is_a_match): mismatchCount += 1
                lq += 1
            # end while
            if(mismatchCount <= 5 or
               mismatchCount <= 0.05*(x_len-lp-rp) ):
                lp=lq
                # sys.stderr.write("# Closed gap by jumping\n");
            # end if
        # end if

        if(0):
            sys.stderr.write( "rawX: %s\n" % x_substring)
            sys.stderr.write( "rawY: %s\n" % y_substring)
            sys.stderr.write( "squX: %s\n" % x_substring[:lp]+ \
                              string.lower(x_substring[lp:x_len-rp]) + \
                              x_substring[x_len-rp:x_len])
            sys.stderr.write( "squY: %s\n" % y_substring[:lp]+ \
                              string.lower(y_substring[lp:y_len-rp]) + \
                              y_substring[y_len-rp:y_len])
            sys.stderr.write( "sqeX: %s\n" % x_substring[lp:x_len-rp])
            sys.stderr.write( "sqeY: %s\n" % y_substring[lp:y_len-rp])
            sys.stderr.write( "x_seg=(%s %s %d %d)" % \
                              (left.x_orientation, left.x_scaf_uid, x_pos+lp, x_len-rp-lp))
            sys.stderr.write( "y_seg=(%s %s %d %d)\n" % \
                              (left.y_orientation, left.y_scaf_uid, y_pos+lp, y_len-rp-lp))

        if(lp+rp > x_len and x_len >= 0):
            sys.stderr.write("overfilledX ")
            sys.stderr.write( "x_seg=(%s %s %d %d) " % \
                              (left.x_orientation, left.x_scaf_uid, x_pos+lp, x_len-rp-lp))
            sys.stderr.write( "y_seg=(%s %s %d %d)\n" % \
                              (left.y_orientation, left.y_scaf_uid, y_pos+lp, y_len-rp-lp))
        if(lp+rp > y_len and y_len >= 0):
            sys.stderr.write("overfilledY ")
            sys.stderr.write( "x_seg=(%s %s %d %d) " % \
                              (left.x_orientation, left.x_scaf_uid, x_pos+lp, x_len-rp-lp))
            sys.stderr.write( "y_seg=(%s %s %d %d)\n" % \
                              (left.y_orientation, left.y_scaf_uid, y_pos+lp, y_len-rp-lp))
        if(lp+rp < x_len and lp+rp < y_len): completefillednotXY += 1
        if(lp+rp == x_len and x_len < y_len): completefilledXnotY += 1
        if(lp+rp == y_len and y_len < x_len): completefilledYnotX += 1
        if(lp+rp == x_len and x_len == y_len): completefilledXandY += 1

        # Print out abutting intervals to fill gaps.
        if(lp>0):
            left_fill = left.copy()
            if(left_forward):
                left_fill.subtype = "L"
                left_fill.x_start = x_pos
                left_fill.x_length = lp
                left_fill.y_start = y_pos
                left_fill.y_length = lp
                left_fill.matchid = left_fill.matchid + "L"
                
            else:
                if(sorted_by_x):
                    left_fill.subtype = "L"
                    left_fill.x_start = x_pos
                    left_fill.x_length = lp
                    left_fill.y_start = y_pos+y_len-lp
                    left_fill.y_length = lp
                    left_fill.matchid = left_fill.matchid + "L"
                else: # assume sorted_by_y
                    left_fill.subtype = "L"
                    left_fill.x_start = x_pos+x_len-lp
                    left_fill.x_length = lp
                    left_fill.y_start = y_pos
                    left_fill.y_length = lp
                    left_fill.matchid = left_fill.matchid + "L"
                # end if
            # end if
            # outfile.write(str(left_fill))
            print >>outfile, left_fill
        # end if
        if(rp>0):
            right_fill = right.copy()
            if(left_forward):
                right_fill.subtype = "R"
                right_fill.x_start = x_pos+x_len-rp
                right_fill.x_length = rp
                right_fill.y_start = y_pos+y_len-rp
                right_fill.y_length = rp
                right_fill.matchid = right_fill.matchid + "R"
            else:
                if(sorted_by_x):
                    right_fill.subtype = "R"
                    right_fill.x_start = x_pos+x_len-rp
                    right_fill.x_length = rp
                    right_fill.y_start = y_pos
                    right_fill.y_length = rp
                    right_fill.matchid = right_fill.matchid + "R"
                else: # assume sorted_by_y
                    right_fill.subtype = "R"
                    right_fill.x_start = x_pos
                    right_fill.x_length = rp
                    right_fill.y_start = y_pos+y_len-rp
                    right_fill.y_length = rp
                    right_fill.matchid = right_fill.matchid + "R"
                # end if
            # end if
            #outfile.write(str(right_fill))
            print >>outfile, right_fill
        # end if
        
        if(0): # Start gap composition diagnostics.
            if( (x_len > lp+rp) or (y_len > lp+rp) ):
                for ch in x_substring[lp:x_len-rp]:
                    if(not(ch=='A' or ch=='C' or ch=='G' or ch=='T')):
                        x_notACGT += 1
                    try:
                        x_chCount[ch] += 1
                    except KeyError:
                        x_chCount[ch] = 1
                for ch in  y_substring[lp:y_len-rp]:
                    if(not(ch=='A' or ch=='C' or ch=='G' or ch=='T')):
                        y_notACGT += 1
                    try:
                        y_chCount[ch] += 1
                    except KeyError:
                        y_chCount[ch] = 1
                if(1 or x_notACGT > 0 or y_notACGT > 0):
                    sys.stderr.write("Ncounts %d %d\n" % (x_notACGT,y_notACGT))
                    sys.stderr.write("x_gap_len= %d y_gap_len= %d\n" % (x_len-lp-rp,y_len-lp-rp))
                    sys.stderr.write("x_seg=(%s %s %d %d)\n" % \
                                     (left.x_orientation, left.x_scaf_uid, x_pos+lp, x_len-lp-rp))
                    sys.stderr.write("y_seg=(%s %s %d %d)\n" % \
                                     (left.y_orientation, left.y_scaf_uid, y_pos+lp, y_len-lp-rp))
                    #                sys.stderr.write("x_chCount= ")
                    #                sys.stderr.write(x_chCount)
                    #                sys.stderr.write("y_chCount= ")
                    #                sys.stderr.write( y_chCount);
    else:
        # sys.stderr.write("Inter-run gap\n")
        inter_run_gap_count += 1
    # sys.stderr.write("done\n")
    squeezed = lp+rp
    return (inter_run_gap_count,
            squeezed,x_len-squeezed,y_len-squeezed,x_notACGT,y_notACGT)
# end def

def mainLoop( inpfile, outfile, xIdx, yIdx):

    maxgap = 100000 # This should be set by an ATAC global.
    margin = 20     # This should be set by an ATAC global.

    countLines = 0
    inter_run_gap_count_total = 0
    closed_gap_count_total = 0
    squeezed_total = 0
    x_len_total = 0
    y_len_total = 0
    x_nonACGT_total = 0
    y_nonACGT_total = 0

    inpfile.seek(0)
    inpfileIter = iter(inpfile)

    left = None
    for line in inpfileIter:
        if(line[0] == 'M'):
            left = MatchRecord.MatchRecord(line)
            # outfile.write(str(left))
            print >>outfile, left
            countLines += 1
            break;

    for line in inpfileIter:
        if(line[0] == 'M'):
            right = MatchRecord.MatchRecord(line)

            if( countLines % 10000 == 0):
                sys.stderr.write("countLines=%d\n" % countLines)
            (inter_run_gap_count,squeezed,x_len,y_len,x_notACGT,y_notACGT) \
              = analyzeGap(xIdx,yIdx,left,right, outfile, maxgap, margin)
            inter_run_gap_count_total += inter_run_gap_count
            squeezed_total += squeezed
            x_len_total += x_len
            y_len_total += y_len
            x_nonACGT_total += x_notACGT
            y_nonACGT_total += y_notACGT
            if(x_len == 0 and y_len == 0): closed_gap_count_total += 1

            # Output the record which was possibly trimmed.
            #outfile.write(str(right))
            print >>outfile, right
            countLines += 1
            left = right
        # end if
    # end for
    
    sys.stderr.write(
        "countLines %d inter_run_gap_count %d closed_gap_count %d squeezed %d x_len %d y_len %d x_nonACGT %d y_nonACGT %d\n" %
        (countLines,inter_run_gap_count_total,closed_gap_count_total,
         squeezed_total,x_len_total,y_len_total,x_nonACGT_total,y_nonACGT_total))

    sys.stderr.write("theIsolatedSNPcount = %d\n" % theIsolatedSNPcount)
    sys.stderr.write("completefillednotXY = %d\n" % completefillednotXY)
    sys.stderr.write("completefilledXnotY = %d\n" % completefilledXnotY)
    sys.stderr.write("completefilledYnotX = %d\n" % completefilledYnotX)
    sys.stderr.write("completefilledXandY = %d\n" % completefilledXandY)

# end def

# Allow each module to have its own main for testing.
if __name__ == '__main__':
    
    inpname = sys.argv[1]
    outname = sys.argv[2]
    xname = sys.argv[3]
    yname = sys.argv[4]
    assemblyId1 = sys.argv[5]
    assemblyId2 = sys.argv[6]

    #    mismatches = checkExactMatches( x, y, inpfile)
    #    sys.stderr.write("mismatches = %d\n" % mismatches)

    xIdx = IdxStore.IdxStore(xname,assemblyId1)
    yIdx = IdxStore.IdxStore(yname,assemblyId2)

    inpfile = open(inpname)
    outfile = open(outname,"w")
    mainLoop( inpfile, outfile, xIdx, yIdx)
    outfile.close()
# end if
