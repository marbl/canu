#!/usr/bin/env python

"""
The environmental variable "PYTHONPATH" is a colon separated list
of directories of imported Python modules (*.py) and C/C++ shared
libraries (*.so for Unix or *.dll for Windows).

Written by Clark Mobarry, Applied Biosystems, 2002-2004.
"""

"""
Known issues:
(1) I need to remove from parameters from the output: /inpname=, /outname=,.
(4) The checkpointing scheme assumes that a previous existing checkpoint file is GOOD.
"""

import os, sys, time, getopt, tempfile
import MyFile
import MatchRecord
import AtacFile
import IdxStore
import UniqueFilter
import PerfectRuns
import TrimMatchOverlaps
import squeezeIntraRunGaps
import localAlignerInterface
import fillIntraRunGaps

#import dedashMatches

STDERR=sys.stderr
STDOUT=sys.stdout

def die(message):
    print >>STDERR, message
    os.exit(1)

def cvm(f,x,y):
    # A cvm variant (flag ? y : x) = (x,y)[f]
    if f :
        return x
    else:
        return y
    # end if
# end def
    
def carp(x):
    print >>STDERR, x
# end if

def ftsize(path):
    return os.stat(path)[6]
# end def

class GlobalParam:
    def __init__(self,line):
        pass
    def __str__ (self):
        return "/%s=%s" % (self._key,self._value)
    def get(self):
        return (self._key,self._value)
    def put(self,key,value):
        (self._key,self._value) = (key,value)
    

def usage (*_):
    print >>STDERR, "Usage: atacdriver.py matchFilePrefix"
# end def

def filterByMatchLength( inpfile, outfile, minimum_length):
    "Only keep matches that are long enough."
    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            FM = MatchRecord.MatchRecord(line)
            if (FM.x_length >= minimum_length and
                FM.y_length >= minimum_length ):
                print >>outfile, FM
            # end if
        # end if
    # end for
# end def


def onlyKeepLongRuns ( inpfile, outname, lengthThreshold ):
    outfile = MyFile.myfile()
    rejectsfile = MyFile.myfile()
    
    FL = None
    store = []
    lenInMatches = 0
    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            FM = MatchRecord.MatchRecord(line)
            SL = FM.x_length
            if FL != None and FL.runid != FM.runid :
                for x in store:
                    print >>rejectsfile, x
                # end for
                store = []
                lenInMatches = SL
            else:
                lenInMatches += SL
            # end if

            if lenInMatches < lengthThreshold:
                store.append(FM)
            else:
                for x in store:
                    print >>outfile, x
                # end for
                store = []
                print >>outfile, FM
            # end if
            FL = FM
        # end if
    # end for
    rejectsfile.close()
    return outfile
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
                print >>STDERR, 'Bombed on:'
                print >>STDERR, str(curF)
                print >>STDERR, 'needs_to_share_diagonal=' + str(needs_to_share_diagonal)
                print >>STDERR, 'nx=' + str(nx) + '  ny=' + str(ny)
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
            # print >>STDOUT, lastF, curF
            # print >>STDOUT, lx,rx,ly,ry
            # print >>STDOUT, lowHitPX,hghHitPX,lowHitPY,hghHitPY
            # print >>STDOUT, "overlapping=",overlapping
            # print >>STDOUT, "on_diagonal=",on_diagonal

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



# Note that if record has an initial rank for the X and Y sorting,
# then re-sorting and box recovery are simplified.

# Resorting becomes making the inital sparse ranking dense then a
# scattering to the destination.

#  NOTE THAT outname is unused here.

def boxRecovery( inpfile, rawfile, outname):
    inpfile.seek(0)
    rawfile.seek(0)
    outfile = MyFile.myfile()

    rawfileIter = iter(rawfile)
    
    # This is a modified merge operation?
    # The two input files must be sorted the same manner.
    leftMatch = None
    for line in inpfile:
        if(line[0] == 'M'):
            rightMatch = MatchRecord.MatchRecord(line)
            if( leftMatch != None and leftMatch.inSameRunAs(rightMatch) ):
                # print >>STDERR, "In same run leftMatch=", leftMatch, " rightMatch=", rightMatch
                for rawline in rawfileIter:
                    if( rawline[0] == 'M'):
                        rawMatch = MatchRecord.MatchRecord(rawline)
                        if(rawMatch.sameAs(rightMatch)):
                            print >>outfile, rightMatch
                            break
                        else:
                            # print "Inside run rawMatch=", rawMatch
                            if(rawMatch.isInsideBox(leftMatch,rightMatch)):
                                print >>outfile, rawMatch
                            # end if
                        # end if
                    # end if
                # end for
                # We should die here if there is no rawMatch that matched the rightMatch ...
            else:
                # print >>STDERR, "Between runs leftMatch=", leftMatch, " rightMatch=", rightMatch
                for rawline in rawfileIter:
                    if( rawline[0] == 'M'):
                        rawMatch = MatchRecord.MatchRecord(rawline)
                        if(rawMatch.sameAs(rightMatch)):
                            print >>outfile, rightMatch
                            break
                        else: 
                            # print >>STDERR, "Discard rawMatch=", rawMatch
                            pass
                        # end if
                    # end if
                # end for
                # We should die here if there is no rawMatch that matched the rightMatch ...
                # Discard raw Matches until it is ge to the right match.
            # end if
            leftMatch = rightMatch
        # end if
    # end for
    return outfile
# end def


class AtacDriver(AtacFile.AtacFile):
    def runOld(self):
        self.globals['atacAlgorithmVersion'] = str(17)
        print >>STDERR, "runName = %s\n" % self.runName

        # The ATAC globals used by this script:
        opt_t = int(self.globals['globalMatchMinSize'])
        opt_l = int(self.globals['globalPerfectRunMinLen'])
        maxdiff = int(self.globals['globalPerfectRunMaxGapLen'])

        assemblyId1 = self.globals['assemblyId1']
        assemblyId2 = self.globals['assemblyId2']

        assemblyFile1 = self.globals['assemblyFile1']
        assemblyFile2 = self.globals['assemblyFile2']

        boxRecoveryOn = 0  # Deprecated for same species comparisons 2003/09/09.
        if(self.globals.has_key("boxRecoveryOn")):
            boxRecoveryOn = int(self.globals['boxRecoveryOn'])
            
        t0 = time.time()

        assemblyIdx1 = IdxStore.IdxStore(assemblyFile1,assemblyId1)
        assemblyIdx2 = IdxStore.IdxStore(assemblyFile2,assemblyId2)
        rawfile = None
        
        ###################################################################
        # Setup for checkpointing scheme.        
        redo = 0
        keep = 0
        step = 0
        if(self.globals.has_key("ckpKeep")):
            keep = int(self.globals['ckpKeep'])
        ckpName = "AllDone"
        ###################################################################

        print >>STDERR, 'Keep step=' + str(keep)
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)

        outprefix = self.runName

        step += 1
        print >>STDERR, 'At uniqueFilter, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            if(not(self.globals.has_key('uniqueFilterOn') and self.globals['uniqueFilterOn']=="0")):
                print >>STDERR, 'Running UniqueFilter'
                outfile = MyFile.myfile()
                UniqueFilter.main( self.matches, outfile)
                self.matches = outfile
                outprefix += '.uniq'
                self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At filterByMatchLength, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'Running filterByMatchLength'
            outfile = MyFile.myfile()
            filterByMatchLength( self.matches, outfile, opt_t)
            self.matches = outfile
            outprefix += '.t' + str(opt_t)
            self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At trimMatchOverlaps, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, "Start trimming for bp one-to-one-ness"
            tempdata = MyFile.myfile()
            TrimMatchOverlaps.trimMatchOverlapsInBoth(self.matches,tempdata,'u')
            self.matches = tempdata
            print >>STDERR, "Finished trimming for bp one-to-one-ness"
            outprefix += '.trim'
            self.checkpoint(outprefix)

        if( boxRecoveryOn == 1 ):
            # For box recovery later ... but what if we start from a checkpoint?
            rawfile = self.matches

        step += 1
        print >>STDERR, 'At formPerfectRuns, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'from ' + outprefix + ' making ' + outprefix + '.p6'
            tempdata = PerfectRuns.formPerfectRuns(self.matches,
                                                   MatchRecord.sortInXorderAP,
                                                   MatchRecord.sortInYorderAP,
                                                   maxdiff,
                                                   'r')
            self.matches = tempdata
            outprefix += ".p6"
        # end if

        step += 1
        print >>STDERR, 'At onlyKeepLongRuns, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'from ' + outprefix + ' making ' + outprefix + '.l' + str(opt_l)
            tempdata = onlyKeepLongRuns( self.matches, outprefix, opt_l)
            self.matches = tempdata
            outprefix += '.l' + str(opt_l)
            self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At formPerfectRuns, step=' + str(step) 
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'Heal the perfect runs'
            tempdata = PerfectRuns.formPerfectRuns(self.matches,
                                       MatchRecord.sortInYorderAP,
                                       MatchRecord.sortInXorderAP, maxdiff, 'r')
            self.matches = tempdata
            outprefix += '.pr'
            self.checkpoint(outprefix)

        if(boxRecoveryOn == 1): 

            # This is a box recovery step.
            step += 1
            print >>STDERR, 'At boxRecovery, step=' + str(step) 
            print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
            if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
                redo = 1
                print >>STDERR, 'from ' + outprefix + ' making ' + outprefix + '.br'
                print >>STDERR, "Make sorted raw matches"
                outfile = MyFile.myfile()
                MatchRecord.sortInXorderAP( rawfile, outfile)
                rawfile = outfile
                print >>STDERR, "perform box recovery"
                tempdata = boxRecovery( self.matches, rawfile, outprefix)
                self.matches = tempdata
                outprefix += '.br'
                self.checkpoint(outprefix)
            # end if

            step += 1
            print >>STDERR, 'At formPerfectRuns, step=' + str(step)
            print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
            if (redo or ( (keep < step) and not self.globals.has_key(ckpName))):
                print >>STDERR, "form perfect runs"
                redo = 1
                print >>STDERR, 'from ' + outprefix + ' to ' + outprefix + '.p6'
                tempdata = PerfectRuns.formPerfectRuns(self.matches,
                                       MatchRecord.sortInXorderAP,
                                       MatchRecord.sortInYorderAP, maxdiff, 'r')
                self.matches = tempdata
                outprefix += '.pr'
                self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At squeezeIntraRunGaps, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'from ' + outprefix + ' to ' + outprefix + '.sq'
            tempdata = MyFile.myfile()
            squeezeIntraRunGaps.mainLoop(
                self.matches,
                tempdata,
                assemblyIdx1, assemblyIdx2)
            tempy = MyFile.myfile()
            # Beware the current match subtypes are 'x', 'L', and 'R'!
            coalesceMatches( tempdata, tempy, 1)
            self.matches = tempy
            outprefix += '.sq'
            self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At TrimMatchOverlaps, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, "Start trimming for bp one-to-one-ness"
            tempdata = MyFile.myfile()
            TrimMatchOverlaps.trimMatchOverlapsInBoth(self.matches,tempdata,'u')
            self.matches = tempdata
            outprefix += '.trim'
            print >>STDERR, "Finished trimming for bp one-to-one-ness"

        step += 1
        print >>STDERR, 'At RunsAsMatches, step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            self.runs = PerfectRuns.runsAsMatches( self.matches)
            outprefix += '.runs'
            self.checkpoint(outprefix)
        # end if

        if(self.globals.has_key('fillIntraRunGapsOn') and self.globals['fillIntraRunGapsOn']=="1" ):
        
            # Next comes the DNA sequence dependent stuff.
            step += 1
            print >>STDERR, 'At fillIntraRunGaps, step=' + str(step)
            print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
            if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
                redo = 1
                print >>STDERR, "fill the intrarun gaps"
                if(not self.globals.has_key('fillIntraRunGapsErate')):
                    self.globals['fillIntraRunGapsErate'] = 0.10
                if(not self.globals.has_key('fillIntraRunGapsMaxGap')):
                    self.globals['fillIntraRunGapsMaxGap'] = 100000
                fillIntraRunGapsErate = float(self.globals['fillIntraRunGapsErate'])
                fillIntraRunGapsMaxGap = int(self.globals['fillIntraRunGapsMaxGap'])
                tempdata = MyFile.myfile()
                fillIntraRunGaps.mainLoop(self.matches, tempdata,
                                          assemblyIdx1, assemblyIdx2,
                                          fillIntraRunGapsMaxGap, fillIntraRunGapsErate)
                self.matches = tempdata
                outprefix += '.fill'
                self.checkpoint(outprefix)

            step += 1
            print >>STDERR, 'At TrimMatchOverlaps, step=' + str(step)
            print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
            if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
                redo = 1
                print >>STDERR, "trim the overlaps"
                tempdata = MyFile.myfile()
                TrimMatchOverlaps.trimMatchOverlapsInBoth(self.matches,tempdata,'u')
                self.matches = tempdata
                outprefix += '.trim'
                self.checkpoint(outprefix)

        # end if "fillIntraRunGapsOn"

    # end def
        
# end class


class localExecutable :
    def __init__(self, name):
        self.name = name
    def run(self,argline,inpfile,outfile):
        cmd = "%s %s %s %s" % (self.name,argline,inpfile,outfile)
        print >>STDERR,"cmd =", cmd
        iret = os.system(cmd)
        assert(iret == 0)


def main(runName):

    if (0):
        junkfile = sys.stderr
        localAlignerInterface.hello(junkfile)
        return

    print >>STDERR, "I read the output from MatchExtender.\n";
    print >>STDERR, "Start reading checkpoint ...."

    t0 = time.time()
    obj = AtacDriver(runName)

    t1 = time.time()
    print >>STDERR, "Read checkpoint in %d seconds." % (t1-t0)
    t0=t1

    # The following are required:
    assert(obj.globals.has_key('assemblyId1'))
    assert(obj.globals.has_key('assemblyId2'))
    assert(obj.globals.has_key('assemblyFile1'))
    assert(obj.globals.has_key('assemblyFile2'))

    assemblyId1 = obj.globals['assemblyId1']
    assemblyId2 = obj.globals['assemblyId2']

    assemblyFile1 = obj.globals["assemblyFile1"]
    assemblyFile2 = obj.globals["assemblyFile2"]

    assert(os.path.exists(assemblyFile1))
    assert(os.path.exists(assemblyFile2))

    if(not os.path.exists(assemblyFile1+".idxStore")):
        IdxStore.createIndexedFasta( assemblyFile1, assemblyId2)

    if(not os.path.exists(assemblyFile2+".idxStore")):
        IdxStore.createIndexedFasta( assemblyFile2, assemblyId2)

    assert(os.path.exists(assemblyFile1+".idxStore"))
    assert(os.path.exists(assemblyFile2+".idxStore"))

    if not obj.globals.has_key('matchesFile'):
        print >>STDERR, "We need to make the raw matches."
            
    if(not obj.globals.has_key('rawMatchMerSize')):
        obj.globals['rawMatchMerSize'] = 20
    if(not obj.globals.has_key('rawMatchMerMaxDegeneracy')):
        obj.globals['rawMatchMerMaxDegeneracy'] = 1
    if(not obj.globals.has_key('rawMatchMinFillSize')):
        obj.globals['rawMatchMinSize'] = obj.globals['rawMatchMerSize']

    # Many 2*rawMatchMerSize-1 matches are due to isolated single
    # nucleotide mutations in otherwise perfect repeats.
    if(not obj.globals.has_key('globalMatchMinSize')):
        obj.globals['globalMatchMinSize'] = 2*int(obj.globals['rawMatchMerSize'])
    if(not obj.globals.has_key('globalPerfectRunMinLen')):
        obj.globals['globalPerfectRunMinLen'] = 100
    if(not obj.globals.has_key('globalPerfectRunMaxGapLen')):
        obj.globals['globalPerfectRunMaxGapLen'] = 100000

    if(not obj.globals.has_key('intraRunGapIsolatedMismatchLen')):
        obj.globals['intraRunGapIsolatedMismatchLen'] = 20

    obj.runOld()

    t1 = time.time()
    print >>STDERR, "Ran in %d seconds." % (t1-t0)
    t0=t1

    obj.checkpoint(runName+".ckpLast")
    t1 = time.time()
    print >>STDERR, "Wrote checkpoint in %d seconds." % (t1-t0)

if __name__ == '__main__':
    main(sys.argv[1])
