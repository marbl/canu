#!/usr/bin/env python

import sys
import string
import MatchRecord
import IdxStore
import localAlignerInterface
import halign
#import shelve

# True=1
False=0

def analyzeGap(x,y,left,right,outfile,maxgap,erate,margin):
    inter_run_gap_count = 0

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
        if( left_forward != right_forward): sys.stderr.write("Bad orientations\n")
        assert(left_forward == right_forward)

        sorted_by_x = (left.x_start <= right.x_start) and \
                      (left.x_start+left.x_length <= right.x_start+right.x_length)
        sorted_by_y = (left.y_start <= right.y_start) and \
                      (left.y_start+left.y_length <= right.y_start+right.y_length)
        if(not(sorted_by_x or sorted_by_y)):
           print >>sys.stderr, "bad sorting in runs"
           print >>sys.stderr, left
           print >>sys.stderr, right
        assert(sorted_by_x or sorted_by_y)
        # This concept of sorted allows neggaps but not containmant in both axes.
        
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


        if( 1
            and 0 < x_len and 0 < y_len
            and x_len < maxgap and y_len < maxgap
            ):
            #sys.stderr.write("About to call local aligner with %d margins\n" % margin);
            #sys.stderr.write("# left = %s\n" % str(left))
            #sys.stderr.write("# right= %s\n" % str(right))
            #sys.stderr.write("x_len=%d y_len=%d\n" % (x_len, y_len) );

            # sys.stderr.write("x_substring=%s\n" % x_substring );
            # sys.stderr.write("y_substring=%s\n" % y_substring );

            # Why two orientation flags?  We want the output matches
            # to be in the same sorted order as the left and right
            # matches.
            
            parent_x_start = x_pos - margin
            parent_y_start = y_pos - margin
            parent_x_length = x_len + 2*margin
            parent_y_length = y_len + 2*margin
            #print >>sys.stderr, "parent_x_start=%d" % parent_x_start
            #print >>sys.stderr, "parent_y_start=%d" % parent_y_start
            #print >>sys.stderr, "parent_x_length=%d" % parent_x_length
            #print >>sys.stderr, "parent_y_length=%d" % parent_y_length

            x_seq = ""
            if(x_len > 0):
                x_seq = string.upper(
                    x.getStringFromFasta( sorted_by_x, left.x_scaf_uid,
                                          parent_x_start, parent_x_length));
            # end if
            y_seq = ""
            if(y_len > 0):
                y_seq = string.upper(
                    y.getStringFromFasta( sorted_by_y, left.y_scaf_uid,
                                          parent_y_start, parent_y_length));
            # end if

            if 0:
                print >>outfile, "# STARTED localAlignerInterface.syntenicSegments"
                print >>outfile, "# left = %s" % str(left)
                print >>outfile, "# right= %s" % str(right)
                print >>sys.stderr, "x_seq="+x_seq
                print >>sys.stderr, "len(x_seq)=",len(x_seq)
                print >>sys.stderr, "y_seq="+y_seq
                print >>sys.stderr, "len(y_seq)=",len(y_seq)
                outfile.flush()

            try:
                localAlignerInterface.syntenicSegments(
                    outfile,
                    x_seq, 0, parent_x_length,
                    y_seq, 0, parent_y_length,
                    erate,
                    )

                FM = left
                parent_id = FM.matchid
                #FM.x_orientation = sorted_by_x
                #FM.y_orientation = sorted_by_y

                # Why two orientations and not just a flipped flag?
                # Because we want the resulting matches to come out in
                # the same sorted order as the input matches.

                ii = 0
                for segment in iter(localAlignerInterface.iterateSegments,None):
                    #print >>outfile, segment
                    (bgn1,bgn2,len1,len2,fid) = segment
                    assert(len1 >= 0)
                    assert(len2 >= 0)
                    assert(bgn1 >= 0)
                    assert(bgn2 >= 0)
                    if(not(bgn1 + len1 <= parent_x_length)):
                        print >>sys.stdout,"# warn(not(bgn1 + len1 <= parent_x_length))"
                        print >>sys.stdout,"# bgn1=%d len1=%d parent_x_length=%d" % (bgn1,len1,parent_x_length)
                        print >>sys.stdout,"# left = %s" % str(left)
                        print >>sys.stdout,"# right= %s" % str(right)
                        print >>sys.stdout,"# bgn1,bgn2,len1,len2=", bgn1,bgn2,len1,len2
                        #print >>sys.stdout,"# xseq=%s" % x_seq[bgn1:bgn1+len1]
                        #print >>sys.stdout,"# yseq=%s" % y_seq[bgn2:bgn2+len2]
                        len1 = parent_x_length - bgn1
                        print >>sys.stdout, "# Change len1 = %d" % len1
                    if(not(bgn2 + len2 <= parent_y_length)):
                        print >>sys.stdout,"# warn(not(bgn2 + len2 <= parent_y_length))"
                        print >>sys.stdout,"# bgn2=%d len2=%d parent_y_length=%d" % (bgn2,len2,parent_y_length)
                        print >>sys.stdout,"# left = %s" % str(left)
                        print >>sys.stdout,"# right= %s" % str(right)
                        print >>sys.stdout,"# bgn1,bgn2,len1,len2=", bgn1,bgn2,len1,len2
                        #print >>sys.stdout,"# xseq=%s" % x_seq[bgn1:bgn1+len1]
                        #print >>sys.stdout,"# yseq=%s" % y_seq[bgn2:bgn2+len2]
                        len2 = parent_y_length - bgn2
                        print >>sys.stdout,"# Change len2 = %d" % len2
                    if (len1 == 0):
                        print >>sys.stdout,"# warn(len1 == 0)"
                        print >>sys.stdout,"# bgn1,bgn2,len1,len2=", bgn1,bgn2,len1,len2
                        print >>sys.stdout,"# bgn1=%d len1=%d parent_x_length=%d" % (bgn1,len1,parent_x_length)
                        continue
                    if (len2 == 0):
                        print >>sys.stdout,"# warn(len2 == 0)"
                        print >>sys.stdout,"# bgn1,bgn2,len1,len2=", bgn1,bgn2,len1,len2
                        print >>sys.stdout,"# bgn2=%d len2=%d parent_y_length=%d" % (bgn2,len2,parent_y_length)
                        continue
                    assert(bgn1 >= 0)
                    assert(bgn2 >= 0)
                    assert(len1 > 0)
                    assert(len2 > 0)
                    assert(bgn1 + len1 <= parent_x_length);
                    assert(bgn2 + len2 <= parent_y_length);
                    # Filter by a minimum length? say four bp.
                    ii += 1
                    FM.subtype = 'l'
                    FM.matchid = parent_id + 'l' + str(ii)
                    # FM.runid = parent_id
                    child_x_start = parent_x_start + (parent_x_length-bgn1-len1,bgn1)[sorted_by_x]
                    child_y_start = parent_y_start + (parent_y_length-bgn2-len2,bgn2)[sorted_by_y]
                    child_x_length = len1
                    child_y_length = len2

                    #FM.identifier = " %f" % fid # CMM BEWARE
                    FM.x_start = child_x_start
                    FM.y_start = child_y_start
                    FM.x_length = child_x_length
                    FM.y_length = child_y_length
                    #FM.extend['fid'] = str(fid)
                    #print >>outfile, FM

                    # Here we call the dedasher.
                    #assert(len1 > 0)
                    #assert(len2 > 0)
                    #assert(bgn1 >= 0)
                    #assert(bgn2 >= 0)
                    #assert(bgn1+len1 <= parent_x_length)
                    #assert(bgn2+len2 <= parent_y_length)


                    #print >>sys.stderr, "# x_seq=%s" % x_seq
                    #print >>sys.stderr, "# y_seq=%s" % y_seq
                    #print >>sys.stderr, "# bgn1,bgn2,len1,len2=", bgn1,bgn2,len1,len2
                    #print >>sys.stderr, "# xseq=%s" % x_seq[bgn1:bgn1+len1]
                    #print >>sys.stderr, "# yseq=%s" % y_seq[bgn2:bgn2+len2]

                    halign.halignStart(x_seq[bgn1:bgn1+len1],y_seq[bgn2:bgn2+len2])
                    outfile.flush()
                    for hsegment in iter(halign.halignDedash,None):
                        #print >>outfile, segment
                        (bgn1h,bgn2h,len1h,len2h,nmat) = hsegment
                        # Filter by a minimum length? say four bp.
                        ii += 1
                        FM.subtype = 'u'
                        FM.matchid = parent_id + 'a' + str(ii)
                        # FM.runid = parent_id
                        FM.x_start = child_x_start + (child_x_length-bgn1h-len1h,bgn1h)[sorted_by_x]
                        FM.y_start = child_y_start + (child_y_length-bgn2h-len2h,bgn2h)[sorted_by_y]
                        FM.x_length = len1h
                        FM.y_length = len2h

                        assert(len1h == len2h)
                        mismatches = 0
                        for ic in range(len1h):
                            if(x_seq[bgn1+bgn1h+ic] != y_seq[bgn2+bgn2h+ic]):
                                mismatches += 1
                        FM.extend['mm'] = str(mismatches)
                        #FM.identifier = ""  # BEWARE
                        print >>outfile, FM
                    halign.halignFree()

                # localAlignerInterface.free()
                # print >>outfile,"# FINISHED localAlignerInterface.syntenicSegments"
            except RuntimeError:
                print >>outfile, "# NOTE syntenicSegments failed between these records"
                #print >>outfile, "# STARTED localAlignerInterface.syntenicSegments"
                #print >>outfile, "# left = %s" % str(left)
                #print >>outfile, "# right= %s" % str(right)
                print >>sys.stderr, "NOTE syntenicSegments failed in fillIntraRunGaps for:"
                print >>sys.stderr, "x_seq="+x_seq
                print >>sys.stderr, "len(x_seq)=",len(x_seq)
                print >>sys.stderr, "y_seq="+y_seq
                print >>sys.stderr, "len(y_seq)=",len(y_seq)
        # end if

    else:
        # sys.stderr.write("Inter-run gap\n")
        inter_run_gap_count += 1
    # sys.stderr.write("done\n")
    return (inter_run_gap_count,)
# end def


def mainLoop( inpfile, outfile, xIdx, yIdx, maxgap, erate):

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
            print >>outfile, left
            countLines += 1
            break;

    for line in inpfileIter:
        if(line[0] == 'M'):
            newRight = MatchRecord.MatchRecord(line)
            if( newRight.subtype == 'u' ):
                right = newRight
                if( countLines % 10000 == 0):
                    sys.stderr.write("countLines=%d\n" % countLines)

                (inter_run_gap_count,) =  analyzeGap( xIdx, yIdx, left,right, outfile, maxgap, erate, margin)
                inter_run_gap_count_total += inter_run_gap_count

                # Output the record which was possibly trimmed.
                print >>outfile, right
                countLines += 1
                left = right
        # end if
    # end for
    
    sys.stderr.write(
        "countLines %d inter_run_gap_count %d \n" % (countLines,inter_run_gap_count_total)
        )


# end def

import AtacFile
import MyFile

def main( inpname, outname):
    obj = AtacFile.AtacFile(inpname)
    assemblyId1 = obj.globals['assemblyId1']
    assemblyId2 = obj.globals['assemblyId2']
    assemblyFilePrefix1 = obj.globals['assemblyFilePrefix1']
    assemblyFilePrefix2 = obj.globals['assemblyFilePrefix2']

    if(not obj.globals.has_key('fillIntraRunGapsErate')):
        obj.globals['fillIntraRunGapsErate'] = 0.10
    if(not obj.globals.has_key('fillIntraRunGapsMaxGap')):
        obj.globals['fillIntraRunGapsMaxGap'] = 100000
    fillIntraRunGapsErate = float(obj.globals['fillIntraRunGapsErate'])
    fillIntraRunGapsMaxGap = int(obj.globals['fillIntraRunGapsMaxGap'])
    
    # mismatches = checkExactMatches( x, y, inpfile)
    # sys.stderr.write("mismatches = %d\n" % mismatches)

    xIdx = IdxStore.IdxStore(assemblyFilePrefix1,assemblyId1)
    yIdx = IdxStore.IdxStore(assemblyFilePrefix2,assemblyId2)

    tempfile = MyFile.myfile()
    mainLoop( obj.matches, tempfile, xIdx, yIdx,
              fillIntraRunGapsMaxGap, fillIntraRunGapsErate)
    obj.matches = tempfile
    obj.checkpoint(outname)


# Allow each module to have its own main for testing.
if __name__ == '__main__':
    inpname = sys.argv[1]
    outname = sys.argv[2]
    main( inpname, outname)
# end if




