#!/usr/bin/env python

import sys
import MyFile
import MatchRecord
        
def cvm(f,x,y):
    # A cvm variant (flag ? x : y) = (x,y)[f]
    if f :
        return x
    else:
        return y
    # end if
# end def
    
def coalesceMatches ( inpfile, outfile, needs_to_share_diagonal ):
    "Coalesce overlapping and abutting matches within the same run."
    
    firstF = None
    lastF = None
    
    lastLX = -3
    lastLY = -4
    lastForward = 0

    lowHitPX = None
    lowHitPY = None
    hghHitPX = None
    hghHitPY = None
    inpfile.seek(0)
    outfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            curF = MatchRecord.MatchRecord(line)
            px = curF.x_start
            nx = curF.x_length
            py = curF.y_start
            ny = curF.y_length
            assert(px >= 0)
            assert(nx >= 0)
            assert(py >= 0)
            assert(ny >= 0)
            if (not (not needs_to_share_diagonal or nx == ny)):
                print >>sys.stderr, 'Bombed on:'
                print >>sys.stderr, str(curF)
                print >>sys.stderr, 'needs_to_share_diagonal=' + str(needs_to_share_diagonal)
                print >>sys.stderr, 'nx=' + str(nx) + '  ny=' + str(ny)
            # end if
            assert((hghHitPX == None or (not needs_to_share_diagonal) or nx == ny))
            forward = (curF.x_orientation == curF.y_orientation)
            lx = px
            ly = cvm( forward, py, py + ny)
            rx = px + nx
            ry = cvm( forward, py + ny, py)

            overlapping = ((lastF != None) and
                           (curF.x_scaf_uid == lastF.x_scaf_uid) and
                           (curF.y_scaf_uid == lastF.y_scaf_uid) and
                           (((lx >= lowHitPX and lx <= hghHitPX) and
                             (ly >= lowHitPY and ly <= hghHitPY)) or
                            ((rx >= lowHitPX and rx <= hghHitPX) and
                             (ry >= lowHitPY and ry <= hghHitPY))))
            on_diagonal = ((forward == lastForward) and
                           ((lx - lastLX) == ((ly - lastLY) * cvm(forward, 1, -1))))
            # print >>sys.stdout, lastF, curF
            # print >>sys.stdout, lx,rx,ly,ry
            # print >>sys.stdout, lowHitPX,hghHitPX,lowHitPY,hghHitPY
            # print >>sys.stdout, "overlapping=",overlapping
            # print >>sys.stdout, "on_diagonal=",on_diagonal

            lowMerPX = px
            lowMerPY = py
            hghMerPX = px + nx
            hghMerPY = py + ny
            if (not (overlapping and (not needs_to_share_diagonal or on_diagonal))):
                if (firstF != None):
                    # if (lastF == None or firstF.runid != lastF.runid):
                    # end if
                    firstF.subtype = ('g','u')[needs_to_share_diagonal]
                    firstF.x_start = lowHitPX
                    firstF.y_start = lowHitPY
                    firstF.x_length = hghHitPX - lowHitPX
                    firstF.y_length = hghHitPY - lowHitPY
                    print >>outfile, firstF
                # end if
                firstF = curF
                lowHitPX = lowMerPX
                lowHitPY = lowMerPY
                hghHitPX = hghMerPX
                hghHitPY = hghMerPY
            # end if
            lowHitPX = cvm(lowHitPX < lowMerPX, lowHitPX, lowMerPX)
            lowHitPY = cvm(lowHitPY < lowMerPY, lowHitPY, lowMerPY)
            hghHitPX = cvm(hghHitPX > hghMerPX, hghHitPX, hghMerPX)
            hghHitPY = cvm(hghHitPY > hghMerPY, hghHitPY, hghMerPY)

            lastLX = lx
            lastLY = ly
            lastForward = forward
            lastF = curF
        # end if
    # end for

    
    if (firstF != None):
        firstF.subtype = ('g','u')[needs_to_share_diagonal]
        firstF.x_start = lowHitPX
        firstF.y_start = lowHitPY
        firstF.x_length = hghHitPX - lowHitPX
        firstF.y_length = hghHitPY - lowHitPY
        print >>outfile, firstF

    return
# end def

def trimMatchOverlapsInX(inpfile,outfile, trim_subtype):
    "Trim the match overlaps with respect to the X assembly."
    overlaps=0
    abuts=0
    posgaps=0
    contained = 0
    trimmed = 0
    left = None

    picket = 0
    # For each genomic axis we scan left to right using this picket
    # position to annihilating any part of the current match to the
    # left of this picket.
    
    inpfile.seek(0)
    for line in iter(inpfile):
        if(line[0] == 'M'):
            right = MatchRecord.MatchRecord(line)
            if( right.subtype != trim_subtype):
                print >>outfile, line,
                continue
            if( left == None or
                #left.x_scaf_uid < right.x_scaf_uid):
                left.x_scaf_uid != right.x_scaf_uid):
                picket = 0
            else:
                assert(left != None)
                assert(right != None)
                if(left.x_scaf_uid > right.x_scaf_uid):
                    print >>sys.stderr, "sequence ids out of x sorted order"
                    print >>sys.stderr, left
                    print >>sys.stderr, right
                assert(left.subtype == right.subtype)
                assert(left.x_scaf_uid == right.x_scaf_uid)
                if(not(left.x_start <= right.x_start)):
                    print >>sys.stderr, "trimMatchOverlapsInX: Woops not sorted anymore!"
                    print >>sys.stderr, left
                    print >>sys.stderr, right
                    #assert(0)

                thisbgn = right.x_start
                thisend = right.x_start + right.x_length
                if(picket < thisend):
                    gaplen = thisbgn - picket
                    if(gaplen > 0):
                        posgaps += 1
                    if(gaplen == 0):
                        abuts += 1
                    if(gaplen < 0):
                        overlaps += 1
                        trimmed -= gaplen
                        right.x_start  -= gaplen # modify the match
                        right.x_length += gaplen
                        right.y_length += gaplen
                        if(right.x_orientation == right.y_orientation):
                            right.y_start  -= gaplen # modify the match
                else:
                    # picketed region contains right.
                    #print >>sys.stderr, "trimMatchOverlapsInX: Contained"
                    #print >>sys.stderr, left
                    #print >>sys.stderr, right
                    contained += 1
                    right = None # remove this match
            if(right != None):
                print >>outfile, right
                newpicket = right.x_start + right.x_length
                assert(picket < newpicket)
                picket = newpicket
                left = right
        else:
            print >>outfile, line,
    print >>sys.stderr, "trimMatchOverlapsInX:\n",
    print >>sys.stderr, "#posgaps, #abuts, #overlaps, #contained, bp_trimmed= %d %d %d %d %d\n" \
          % (posgaps, abuts, overlaps, contained, trimmed, )
    return

def trimMatchOverlapsInY(inpfile,outfile, trim_subtype):
    "Trim the match overlaps with respect to the Y assembly."
    overlaps=0
    abuts=0
    posgaps=0
    contained = 0
    trimmed = 0
    left = None

    picket = 0
    # For each genomic axis we scan left to right using this picket
    # position to annihilating any part of the current match to the
    # left of this picket.
    
    inpfile.seek(0)
    for line in iter(inpfile):
        if(line[0] == 'M'):
            right = MatchRecord.MatchRecord(line)
            if( right.subtype != trim_subtype):
                print >>outfile, line,
                continue
            if( left == None or
                #left.y_scaf_uid < right.y_scaf_uid):
                left.y_scaf_uid != right.y_scaf_uid):
                picket = 0
            else:
                assert(left != None)
                assert(right != None)
                if(left.y_scaf_uid > right.y_scaf_uid):
                    print >>sys.stderr, "sequence ids out of y sorted order"
                    print >>sys.stderr, left
                    print >>sys.stderr, right
                assert(left.subtype == right.subtype)
                assert(left.y_scaf_uid == right.y_scaf_uid)
                if(not(left.y_start <= right.y_start)):
                    print >>sys.stderr, "trimMatchOverlapsInY: Woops not sorted anymore!"
                    print >>sys.stderr, left
                    print >>sys.stderr, right
                    #assert(0)

                thisbgn = right.y_start
                thisend = right.y_start + right.y_length
                if(picket < thisend):
                    gaplen = thisbgn - picket
                    if(gaplen > 0):
                        posgaps += 1
                    if(gaplen == 0):
                        abuts += 1
                    if(gaplen < 0):
                        overlaps += 1
                        trimmed -= gaplen
                        right.y_start  -= gaplen # modify the match
                        right.y_length += gaplen
                        right.x_length += gaplen
                        if(right.x_orientation == right.y_orientation):
                            right.x_start  -= gaplen # modify the match
                else:
                    # picketed region contains right.
                    #print >>sys.stderr, "trimMatchOverlapsInY: Contained"
                    #print >>sys.stderr, left
                    #print >>sys.stderr, right
                    contained += 1
                    right = None # remove this match
            if(right != None):
                print >>outfile, right
                newpicket = right.y_start + right.y_length
                assert(picket < newpicket)
                picket = newpicket
                left = right
        else:
            print >>outfile, line,
    print >>sys.stderr, "trimMatchOverlapsInY:\n",
    print >>sys.stderr, "#posgaps, #abuts, #overlaps, #contained, bp_trimmed= %d %d %d %d %d\n" \
          % (posgaps, abuts, overlaps, contained, trimmed, )
    return


def trimMatchOverlapsInBoth(inpfile,outfile,trim_subtype):
    gp = MyFile.myfile()
    MatchRecord.sortInXorderAP(inpfile,gp)
    # The following coalescing assumes perfect runs.
    hp = MyFile.myfile()
    coalesceMatches( gp, hp, ((trim_subtype == 'x') or (trim_subtype == 'u')) )
    gp = MyFile.myfile()
    trimMatchOverlapsInX(hp,gp,trim_subtype)
    hp = MyFile.myfile()
    MatchRecord.sortInYorderAP(gp,hp)
    trimMatchOverlapsInY(hp,outfile,trim_subtype)
    return

def main(inpname, outname, trim_subtype):
    inpfile = open(inpname)
    outfile = open(outname,"w")
    trimMatchOverlapsInBoth(inpfile,outfile,trim_subtype)
    
if __name__ == '__main__':
    inpname = sys.argv[1]
    outname = sys.argv[2]
    trim_subtype = sys.argv[3]
    main(inpname, outname, trim_subtype)
