#!/usr/bin/env python

# dedashMatches.py /prod/IR02/synteny/mus-vs-rat/mouse_celera_R13_chr_20030210-vs-rat_celera_R1_chr_20030507-V3.atac.t20.l100.br.squeezed.filled.coalesced mus-vs-rat.out /prod/IR05/GENOMES/mouse_celera_R13_chr_20030210 /prod/IR05/GENOMES/rat_celera_R1_chr_20030507 MR13 RR1

# dedashMatches.py mouse_celera_R13_chr_20030210-vs-rat_celera_R1_chr_20030507-V3.atac.t20.l100.br.squeezed.filled.coalesced mus-vs-rat.out mouse_celera_R13_chr_20030210 rat_celera_R1_chr_20030507 MR13 RR1


import sys
import string
import time
import MatchRecord
import IdxStore
import halign

#import shelve

class dedasher:
    def __init__(self,xstr,ystr):
        pass
    def __iter__(self):
        return iter([1])

x = 3
def suba():
    global x
    x = 7
def subb():
    global x
    x -= 1
    return (x,None)[x == 0]
def subc():
    suba()
    it = iter(subb,None)
    for y in it:
        print y

def main( inpfile, outfile, xIdx, yIdx):
    inpfile.seek(0)
    outfile.seek(0)
    lineCount = 0
    t0 = time.time()
    for line in inpfile:
        lineCount += 1
        if((lineCount % 10000)==0):
            print >>sys.stderr, "lineCount=",lineCount," time=",time.time()-t0
        if(line[0] == 'M'):
            FM = MatchRecord.MatchRecord(line)
            if(FM.subtype == 'g'):
                parentid = FM.matchid
                parent_x_forward = (FM.x_orientation == 1)
                parent_y_forward = (FM.y_orientation == 1)
                parent_x_start = FM.x_start
                parent_y_start = FM.y_start
                parent_x_length = FM.x_length
                parent_y_length = FM.y_length

                # Why two orientations and not just a flipped flag?
                # Because we want the resulting matches to come out in
                # the same sorted order as the input matches.
                
                x_substring = string.upper(
                    xIdx.getStringFromFasta( parent_x_forward,
                                          FM.x_scaf_uid, FM.x_start, FM.x_length));
                y_substring = string.upper(
                    yIdx.getStringFromFasta( parent_y_forward,
                                          FM.y_scaf_uid, FM.y_start, FM.y_length));
                ii = 0
                # Here we call the dedasher.
                halign.halignStart(x_substring,y_substring)
                for segment in iter(halign.halignDedash,None):
                    #print >>outfile, segment
                    (bgn1,bgn2,len1,len2,nmat) = segment
                    # Filter by a minimum length? say four bp.
                    ii += 1
                    FM.subtype = 'u'
                    FM.matchid = parentid + 'u' + str(ii)
                    # FM.runid = parentid
                    FM.x_start = parent_x_start + (parent_x_length-bgn1-len1,bgn1)[parent_x_forward]
                    FM.y_start = parent_y_start + (parent_y_length-bgn2-len2,bgn2)[parent_y_forward]
                    FM.x_length = len1
                    FM.y_length = len2
                    assert(len1 == len2)
                    mismatches = 0
                    for ic in range(len1):
                        if(x_seq[bgn1+ic] != y_seq[bgn2+ic]):
                            mismatches += 1
                    FM.extend['mm'] = str(mismatches)
                    FM.identifier = ""  # BEWARE
                    print >>outfile, FM
                # halign.halignFree()
            else:
                print >>outfile, line,
        else:
            print >>outfile, line,


def oldmain():    
    inpname = sys.argv[1]
    outname = sys.argv[2]
    xIndexName = sys.argv[3]
    yIndexName = sys.argv[4]
    assemblyId1 = sys.argv[5]
    assemblyId2 = sys.argv[6]

    #    mismatches = checkExactMatches( x, y, inpfile)
    #    sys.stderr.write("mismatches = %d\n" % mismatches)

    inpfile = open(inpname)
    outfile = open(outname,"w")
    xIdx = IdxStore.IdxStore(xIndexName,assemblyId1)
    yIdx = IdxStore.IdxStore(yIndexName,assemblyId2)

    main( inpfile, outfile, xIdx, yIdx)
    outfile.close()

import AtacFile
import MyFile

def newmain():
    inpname = sys.argv[1]
    outname = sys.argv[2]

    obj = AtacFile.AtacFile(inpname)
    xname = obj.globals["assemblyFilePrefix1"]
    yname = obj.globals["assemblyFilePrefix1"]
    assemblyId1 = obj.globals["assemblyId1"]
    assemblyId2 = obj.globals["assemblyId2"]


    xIdx = IdxStore.IdxStore(xname,assemblyId1)
    yIdx = IdxStore.IdxStore(yname,assemblyId2)

    inpfile = obj.matches
    outfile = MyFile.myfile()
    main( inpfile, outfile, xIdx, yIdx)
    obj.matches = outfile
    obj.checkpoint(outname)
    outfile.close()
    
    
# Allow each module to have its own main for testing.
if __name__ == '__main__':
    newmain()
# end if
