#!/usr/bin/env python

import sys, os, copy, string, tempfile
import MatchRecord

def countMisMatches( inpfile, outfile, xIdx, yIdx):
    inpfile.seek(0)
    outfile.seek(0)
    mismatches = 0
    for line in inpfile:
        if(line[0] == "M"):
            mr = MatchRecord.MatchRecord(line)
            if(mr.subtype == 'u'):
                x_substring = string.upper(
                    xIdx.getStringFromFasta( (mr.x_orientation==1), mr.x_scaf_uid,
                                             mr.x_start, mr.x_length))
                y_substring = string.upper(
                    yIdx.getStringFromFasta( (mr.y_orientation==1), mr.y_scaf_uid,
                                             mr.y_start, mr.y_length))

                len_x_substring = len(x_substring)
                len_y_substring = len(y_substring)
                assert(len_x_substring == len_y_substring)
                mismatches = 0
                for ii in range(len_x_substring):
                    if(x_substring[ii] != y_substring[ii]):
                        mismatches += 1
                mr.extend['mm'] = mismatches

                if 0:
                    prefixmismatches = 0
                    for ii in range(len_x_substring):
                        if(x_substring[ii] != y_substring[ii]):
                            prefixmismatches += 1
                        else:
                            break
                    mr.extend['mmPrefix'] = prefixmismatches
                    suffixmismatches = 0
                    for ii in range(len_x_substring):
                        if(x_substring[len_x_substring-ii-1] != y_substring[len_x_substring-ii-1]):
                            suffixmismatches += 1
                        else:
                            break
                    mr.extend['mmSuffix'] = suffixmismatches

                print >>outfile, mr

                #sys.stderr.write("Inpfile: <%s>\n" % line)
                #sys.stderr.write("x_orientation=%s x_scaf_uid=%s x_start=%s x_length=%s\n" %
                # (x_orientation, x_scaf_uid, x_start, x_length));
                #sys.stderr.write("y_orientation=%s y_scaf_uid=%s y_start=%s y_length=%s\n" %
                # (y_orientation, y_scaf_uid, y_start, y_length));
                # print >>outfile, "X <%s>"% x_substring
                # print >>outfile, "Y <%s>"% y_substring
            else:
                print >>outfile, line,
        else:
            print >>outfile, line,
    return
# end def

import IdxStore
import AtacFile
import MyFile

def oldmain():
    inpname = sys.argv[1]
    outname = sys.argv[2]
    xname = sys.argv[3]
    yname = sys.argv[4]
    assemblyId1 = sys.argv[5]
    assemblyId2 = sys.argv[6]

    xIdx = IdxStore.IdxStore(xname,assemblyId1)
    yIdx = IdxStore.IdxStore(yname,assemblyId2)

    inpfile = open(inpname)
    outfile = open(outname,"w")
    countMisMatches( inpfile, outfile, xIdx, yIdx)
    outfile.close()


def main(inpname, outname):

    obj = AtacFile.AtacFile(inpname)
    xname = obj.globals["assemblyFilePrefix1"]
    yname = obj.globals["assemblyFilePrefix2"]
    assemblyId1 = obj.globals["assemblyId1"]
    assemblyId2 = obj.globals["assemblyId2"]

    xIdx = IdxStore.IdxStore(xname,assemblyId1)
    yIdx = IdxStore.IdxStore(yname,assemblyId2)

    inpfile = obj.matches
    outfile = MyFile.myfile()
    countMisMatches( inpfile, outfile, xIdx, yIdx)
    obj.matches = outfile
    obj.checkpoint(outname)
    outfile.close()

# Allow each module to have its own main for testing.
if __name__ == '__main__':
    inpname = sys.argv[1]
    outname = sys.argv[2]
    main(inpname,outname)
# end if
