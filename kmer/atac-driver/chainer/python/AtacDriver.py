#!/usr/bin/env python

"""
Extensive documentation for the Python language is available at
http://www.python.org.

The environmental variable "PYTHONPATH" is a colon separated list
of directories of imported Python modules (*.py) and C/C++ shared
libraries (*.so for Unix or *.dll for Windows).  For example, using
the bash shell command in my environment 
"PYTHONPATH=${PYTHONPATH}:$WORK/cds/IR/COMPASS/src/AtacPipeline; export PYTHONPATH"
initializes the PYTHONPATH properly.

This script was tested using Python 2.3.3.  Check using "python
-V" to determine which version you are using.  There are no known
problems with the planned future features of Python.

If you have problems that look platform specific, then search for
"sys.platform" to find computer and OS platform dependencies in the
script.

Written by Clark Mobarry, Applied Biosystems, 2002-2004.
"""

"""
Known issues:
(1) I need to remove from parameters from the output: /inpname=, /outname=,.
(2) Should I remove the requirement of the input parameters: /bindir=, /bridir=,
(3) /assemblyFilePrefix1 is a misnomer.
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

import fillIntraRunGaps
#import dedashMatches
import countMisMatches

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

    def noiseFilter(self):
        pass
    def globalRunFinder(self):
        pass
    def matchRecovery(self):
        pass
    def extractMatchExtension(self):
        pass
    def ungappedMatchExtension(self):
        pass
    def gapFillerPerGap(self):
        pass
    
    def gatekeeper(self):
        # Validates the initial atac checkpoint file.
        # Verifies that the data and executable files are all present.
        # Indexes the fasta files if necessary.
        pass
    # end def

    def rawMatchFinder(self):
        # Calls Brian's meryl and searchGENOMEexactly
        pass

    def chainer(self):
        self.noiseFilter()
        self.globalRunFinder()
        self.matchRecovery()

    def gapFiller(self):
        self.extractMatchExtension()
        self.ungappedMatchExtension()
        self.gapFillerPerGap()

    def terminator(self):
        pass

    def runNew(self):

        self.gatekeeper()
        # The gatekeeper requires a run parameters file and the two FASTA files.
        # The gatekeeper insures that the FASTA files are indexed.
        self.checkpoint()

        self.rawMatchFinder()
        # The rawMatchFinder takes the two indexed FASTA files and produces the
        # initial matches to be chained. 
        self.checkpoint()
        
        self.chainer()
        self.checkpoint()
        
        self.gapFiller()
        self.checkpoint()

    def runOld(self):
        self.globals['atacAlgorithmVersion'] = str(17)
        print >>STDERR, "runName = %s\n" % self.runName

        # The ATAC globals used by this script:
        opt_t = int(self.globals['globalMatchMinSize'])
        opt_l = int(self.globals['globalPerfectRunMinLen'])
        maxdiff = int(self.globals['globalPerfectRunMaxGapLen'])

        assemblyId1 = self.globals['assemblyId1']
        assemblyId2 = self.globals["assemblyId2"]

        assemblyFilePrefix1 = self.globals['assemblyFilePrefix1']
        assemblyFilePrefix2 = self.globals['assemblyFilePrefix2']

        boxRecoveryOn = 0  # Deprecated for same species comparisons 2003/09/09.
        if(self.globals.has_key("boxRecoveryOn")):
            boxRecoveryOn = int(self.globals["boxRecoveryOn"])
            
        t0 = time.time()

        assemblyIdx1 = IdxStore.IdxStore(assemblyFilePrefix1,assemblyId1)
        assemblyIdx2 = IdxStore.IdxStore(assemblyFilePrefix2,assemblyId2)
        rawfile = None
        
        ###################################################################
        # Setup for checkpointing scheme.        
        redo = 0
        keep = 0
        step = 0
        if(self.globals.has_key("ckpKeep")):
            keep = int(self.globals["ckpKeep"])
        ckpName = "AllDone"
        ###################################################################

        print >>STDERR, 'Keep step=' + str(keep)
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)

        outprefix = self.runName

        step += 1
        print >>STDERR, 'At step=' + str(step)
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
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'Running filterByMatchLength'
            outfile = MyFile.myfile()
            filterByMatchLength( self.matches, outfile, opt_t)
            #self.globals[ckpName]=1
            self.matches = outfile
            outprefix += '.t' + str(opt_t)
            self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At step=' + str(step)
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
        print >>STDERR, 'At step=' + str(step)
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
            #self.globals[ckpName]=1
            outprefix += ".p6"
        # end if

        step += 1
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, 'from ' + outprefix + ' making ' + outprefix + '.l' + str(opt_l)
            tempdata = onlyKeepLongRuns( self.matches, outprefix, opt_l)
            #self.globals[ckpName]=1
            self.matches = tempdata
            outprefix += '.l' + str(opt_l)
            self.checkpoint(outprefix)

        step += 1
        print >>STDERR, 'At step=' + str(step) 
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
            #self.globals[ckpName]=1

        if(boxRecoveryOn == 1): 

            # This is a box recovery step.
            step += 1
            print >>STDERR, 'At step=' + str(step) 
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
                #self.globals[ckpName]=1
            # end if

            step += 1
            print >>STDERR, 'At step=' + str(step)
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
                #self.globals[ckpName]=1

        step += 1
        print >>STDERR, 'At step=' + str(step)
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
            #self.globals[ckpName]=1

        step += 1
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, "Start trimming for bp one-to-one-ness"
            tempdata = MyFile.myfile()
            TrimMatchOverlaps.trimMatchOverlapsInBoth(self.matches,tempdata,'u')
            self.matches = tempdata
            outprefix += '.trim'
            #self.checkpoint(outprefix)
            print >>STDERR, "Finished trimming for bp one-to-one-ness"

        step += 1
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            self.runs = PerfectRuns.runsAsMatches( self.matches)
            outprefix += '.runs'
            self.checkpoint(outprefix)
            #self.globals[ckpName]=1
        # end if

        if(self.globals.has_key('fillIntraRunGapsOn') and self.globals['fillIntraRunGapsOn']=="1" ):
        
            # Next comes the DNA sequence dependent stuff.
            step += 1
            print >>STDERR, 'At step=' + str(step)
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
                fillIntraRunGaps.mainLoop( self.matches, tempdata,
                                           assemblyIdx1, assemblyIdx2,
                                           fillIntraRunGapsMaxGap, fillIntraRunGapsErate)
                #self.globals[ckpName]=1
                self.matches = tempdata
                outprefix += '.fill'
                self.checkpoint(outprefix)

            step += 1
            print >>STDERR, 'At step=' + str(step)
            print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
            if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
                redo = 1
                print >>STDERR, "trim the overlaps"
                tempdata = MyFile.myfile()
                TrimMatchOverlaps.trimMatchOverlapsInBoth(self.matches,tempdata,'u')
                #self.globals[ckpName]=1
                self.matches = tempdata
                outprefix += '.trim'
                self.checkpoint(outprefix)

        # end if "fillIntraRunGapsOn"

        step += 1
        print >>STDERR, 'At step=' + str(step)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)
        if (redo or ((keep < step) and not self.globals.has_key(ckpName))):
            redo = 1
            print >>STDERR, "count the number of substitutions"
            tempdata = MyFile.myfile()
            countMisMatches.countMisMatches(self.matches, tempdata, assemblyIdx1, assemblyIdx2)
            #self.globals[ckpName]=1
            self.matches = tempdata
            outprefix += '.mm'
            self.checkpoint(outprefix)
        print >>STDERR, 'Time elapsed=' + str(time.time()-t0)

    # end def
        
# end class

def countMers(id, mersize, merlimit, fastaName, MERYLdir, merylexec ):
    """Check that meryl is finished for each of the inputs"""

    root1 = id + ".k" + str(mersize)
    target1 = os.path.join(MERYLdir,root1)
    if not (os.path.exists( target1 + ".mcdat") and os.path.exists( target1 + ".mcidx")):
        cmd = merylexec
        cmd += " -B -C "
        cmd += " -m " + str(mersize)
        cmd += " -s " + fastaName
        cmd += " -o " + target1
        cmd += " -memory 1500"
        cmd += " -stats " + target1 + ".stats"
        print >>STDOUT, cmd
        iret = os.system(cmd)
        if(iret != 0):
            os.unlink( target1 + ".mcidx")
            os.unlink( target1 + ".mcdat")
            print >>STDERR, "Failed to count mers in " + id
            os.exit(1)

    root2 = id + ".k" + str(mersize) + ".le" + str(merlimit)
    target2 = os.path.join(MERYLdir,root2)
    if not (os.path.exists(target2+".mcdat") and os.path.exists(target2+".mcidx")) :
        cmd = merylexec
        cmd += " -v "
        cmd += " -M lessthanorequal " + str(merlimit)
        cmd += " -s " + target1
        cmd += " -o " + target2
        cmd += " -stats " + target2 + ".stats"
        print >>STDOUT, cmd
        iret = os.system(cmd)
        if(iret != 0):
            os.unlink(target2 + ".mcidx")
            os.unlink(target2 + ".mcdat")
            print >>STDERR, "Failed to count mers lessthanorequal " + merlimit + " in " + id
            os.exit(1)

    return root2

# end def 

def findSegmentsList(numsegments,leaff,seq,prefix1,):
    
    segmentsList = []

    segmentIDs = [0,]

    
    partitionsName = prefix1 + ".partitions"
    cmd = leaff + " -F " + seq + " --partition " + str(numsegments) + " > " + partitionsName
    #CMM The stdout redirection is a Unix-ism.
    iret = os.system(cmd)
    if iret != 0:
        die("leaff partitioning failed")

    FTEST = file(partitionsName,"r")
    ftestiter = iter(FTEST)
    numsegments = int(ftestiter.readline()) # over-write if necessary
    assert(numsegments >= 1)

    segmentID  = 0
    for line in ftestiter:
        pieces = line.split()
        prefix2 = prefix1 + "-segment-" + str(segmentID)

        # the first piece looks like 2](1545909615)
        # after the first piece 1(243615958) 2(199344050) ...
        segment = []
        for piece in pieces[1:]:
            #CMM if ($piece =~ m/(\d+)\(\d+\)/) {
            #CMM     segments += "$1\n";
            #CMM else:
            #CMM     die("Error parsing segment: "+piece)
            segment.append(int(piece.split("(")[0]))
        # end if

        segmentIDs.append(segmentID)
        segmentID += 1

        segmentsList.append(segment)
    # end for
    FTEST.close()

    return segmentsList
# end def


def segmented_search_routine(
    seatac,
    seq1,
    seq2,
    prefix1,
    tmpdir,
    mersize,
    minfill,
    maxgap,
    segmentsList,
    numthreads = 1,
    segmentIDtorun = None,
    buildOnly = None,
    filtername = None,
    filteropts = None,
    ):
    
    for segmentID in xrange(len(segmentsList)):
        #  Now, for each segment that hasn't run, run it.
        prefix2 = prefix1 + "-segment-" + str(segmentID)

        if segmentIDtorun != None and segmentID != segmentIDtorun : continue
        if buildOnly != None and os.path.exists(prefix2 + ".table") : continue

        for segment in segmentsList[segmentID:segmentID+1]:
            prefix2 = prefix1 + "-segment-" + str(segmentID)
            S = file(prefix2,"w")
            print >>STDOUT, "segment", segmentID, segment
            for x in segment:
                print >>S, x
            S.close()
        # end for

        if not (os.path.exists(prefix2 + ".stats") and
                os.path.exists(prefix2 + ".matches")):
            # The matches file is what we want, but Brian only guarantees that
            # the matches file is complete if the stats file was written.
            cmd  = seatac
            cmd += " -verbose";
            cmd += " -mersize " + str(mersize)
            cmd += " -minlength " + str(minfill)
            cmd += " -maxgap " + str(maxgap)
            cmd += " -numthreads " + str(numthreads)
            cmd += " -cdna " + str(seq1)
            cmd += " -genomic " + str(seq2)
            if os.path.exists(prefix1 + ".include.mcdat"):
                cmd += " -only " + prefix1 + ".include"
            if os.path.exists(prefix1 + ".exclude.mcdat"):
                cmd += " -mask " + prefix1 + ".exclude"
            if not os.path.exists(prefix2):
                die("The %d segment file is missing." % segmentID)
            cmd += " -use " + prefix2
            cmd += " -output " + prefix2 + ".matches"
            cmd += " -stats " + prefix2 + ".stats"
            if buildOnly != None :
                cmd += " -buildonly " + prefix2 + ".table"
            if filtername != None :
                cmd += " -filtername " + filtername
            if filtername != None and filteropts != None :
                cmd += " -filteropts " + filteropts

            #  Prevent me from overwriting a run in progress
            if os.path.exists( prefix2 + ".matches") :
                die("WARNING:  Matches already exist!  Exiting!")
            # end if

            F = file(prefix2 + ".cmd", "w");
            F.writelines(cmd)
            F.close()

            iret = os.system(cmd)
            if iret != 0 :
                os.unlink(prefix2 + ".matches.crash")
                os.path.rename(prefix2 + ".matches", prefix2 + ".matches.crash")
                os.unlink(prefix2 + ".stats.crash")
                os.path.rename(prefix2 + ".stats", prefix2 + ".stats.crash")
                die("ERROR: Failed to run " + prefix2)
            # end if

        # end if
    # end for

    #  End early if the segment id to run is defined.
    #
    if segmentIDtorun != None :
        print >>STDERR, "Terminating execution because a specific segmentID was supplied."
        os.exit(0);
    # end if

    #
    #  Join and sort the matches
    #
    if not os.path.exists(prefix1 + ".matches.sorted"):
        mfiles = ""

        #  Check that each search finished, and build a list of all the match files.
        #
        for segmentID in xrange(len(segmentsList)):
            prefix2 = prefix1 + "-segment-" + str(segmentID)
            if os.path.exists( prefix2 + ".stats") and os.path.exists( prefix2 + ".matches"):
                mfiles += prefix2 + ".matches "
            else:
                die(prefix2 + ".matches failed to complete.")
            # end if
         # end for

        #  Ok, all the matches look good, so we can remove the tables.
        for segmentID in xrange(len(segmentsList)):
            prefix2 = prefix1 + "-segment-" + str(segmentID)
            if os.path.exists(prefix2 + ".table"):
                os.unlink( prefix2 + ".table")
        # end for

        cmd = "cat " + mfiles + " | sort -y -T " + tmpdir
        #CMM stdout redirection is a UNIX-ism.
        cmd += " -k 3n -k 7n -k 4n -k 8n "
        cmd += " > " + prefix1 + ".matches.sorted"
        #CMM stdout redirection is a UNIX-ism.
        iret = os.system(cmd)
        if iret != 0 :
            die("Failed to " + cmd)
        # end if

        if 0:
            cmd = "bzip2 " + prefix1 + ".matches.sorted"
            iret = os.system(cmd)
            if iret != 0 :
                die("Failed to " + cmd)
            # end if

        for segmentID in xrange(len(segmentsList)):
            prefix2 = prefix1 + "-segment-" + str(segmentID)
            if os.path.exists(prefix2 + ".matches"):
                os.unlink(prefix2 + ".matches")
        # end for
        
    # end if
    return prefix1 + ".matches.sorted"
# end def

def makeRawMatches(
    assemblyId1,
    assemblyId2,
    seq1,
    seq2,
    MERYLdir,
    Bridir,
    bindir,
    mersize,
    minfill,
    merlimit,
    maxgap,
    numsegments,
    numthreads,
    merylOnly = None,
    segmentIDtorun = None,
    buildOnly = None,
    filtername = None,
    filteropts = None,
    ):

    print >>STDOUT, "makeRawMatches"
    print >>STDOUT, "assemblyId1 = ", assemblyId1
    print >>STDOUT, "assemblyId2 = ", assemblyId2
    print >>STDOUT, "seq1        = ", seq1
    print >>STDOUT, "seq2        = ", seq2

    print >>STDOUT, 'mersize     = ', mersize
    print >>STDOUT, 'merlimit    = ', merlimit
    print >>STDOUT, 'minfill     = ', minfill
    print >>STDOUT, 'maxgap      = ', maxgap
    print >>STDOUT, 'numsegments = ', numsegments
    print >>STDOUT, 'numthreads  = ', numthreads
    print >>STDOUT, 'merylOnly   = ', merylOnly

    print >>STDOUT, ' MERYLdir = ', MERYLdir
    print >>STDOUT, ' Bridir = ', Bridir
    print >>STDOUT, ' bindir = ', bindir

    if segmentIDtorun: print >>STDOUT, 'segmentIDtorun = ', segmentIDtorun
    if buildOnly: print >>STDOUT, 'buildOnly   = ', buildOnly
    if filtername: print >>STDOUT, 'filtername = ' + filtername
    if filteropts: print >>STDOUT, 'filteropts = ' + filteropts
    
    leaff   = os.path.join(bindir,'leaff')
    meryl   = os.path.join(bindir,'meryl')
    existDB = os.path.join(bindir,'existDB')
    seatac  = os.path.join(bindir,'seatac')

    if not os.path.exists(leaff)  : raise ("Can't find " + leaff + "\n")
    if not os.path.exists(meryl)  : raise ("Can't find " + meryl + "\n")
    if not os.path.exists(existDB): raise ("Can't find " + existsDB + "\n")
    if not os.path.exists(seatac) : raise ("Can't find " + seatac + "\n")
    # Brian Walenz would also check if the file is executable.
    
    if not MERYLdir  : raise ("Unset MERYLdir?'\n")
    if not Bridir   : raise ("Unset Bridir?'\n")
    
    if not os.path.isdir(MERYLdir) : raise ("Can't find the MERYLdir '" + MERYLdir + "'\n")
    if not os.path.isdir(Bridir) : os.mkdir(Bridir)

    tmpdir = os.path.join(Bridir,"tmp")
    if os.path.exists(tmpdir) and not os.path.isdir(tmpdir):
        os.unlink(tmpdir)
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    mercount1 = countMers(assemblyId1, mersize, merlimit, seq1, MERYLdir, meryl)
    mercount2 = countMers(assemblyId2, mersize, merlimit, seq2, MERYLdir, meryl)
    matches = assemblyId1 + 'vs' + assemblyId2 + '.k' + str(mersize) + '.u' + str(merlimit) + '.f' + str(minfill) + '.g' + str(maxgap)

    prefix1 = os.path.join(Bridir,matches)

    if not os.path.exists(prefix1 + '.mask.done'):
        minFile = 'min.' + mercount1 + '.' + mercount2
        source1 = os.path.join(MERYLdir, mercount1)
        source2 = os.path.join(MERYLdir, mercount2)
        target3 = os.path.join(Bridir, minFile)
        if (not os.path.exists(target3 + '.mcdat')):
            print >>STDERR, 'Finding the min count between ' + mercount1 + ' and ' + mercount2 + "."
            cmd = meryl
            cmd += ' -M min '
            cmd += ' -s ' + source1
            cmd += ' -s ' + source2
            cmd += ' -o ' + target3
            cmd += ' -stats ' + target3 + '.stats'
            iret = os.system(cmd)
            if iret != 0 :
                os.unlink(target3 + '.mcidx')
                os.unlink(target3 + '.mcdat')
                die('Failed to find the min count between ' + mercount1 + ' and ' + mercount2)


        if not os.path.exists(target3 + '.mcdat') : raise ("Failed to make the mask?\n")
        includeSize = os.path.getsize(target3 + '.mcdat')
        excludeSize = os.path.getsize(source1 + '.mcdat') - os.path.getsize(target3 + '.mcdat')
        print >>STDERR, 'includeSize is about ',includeSize
        print >>STDERR, 'excludeSize is about ',excludeSize
        if includeSize < excludeSize:
            os.rename(target3 + '.mcidx', prefix1 + '.include.mcidx')
            os.rename(target3 + '.mcdat', prefix1 + '.include.mcdat')

        else:
            if (not os.path.exists(prefix1 + '.exclude.mcdat')):
                print >>STDERR, "Finding 'exclude' mers!"
                cmd = meryl
                cmd += ' -M xor'
                cmd += ' -s ' + os.path.join(MERYLdir, assemblyId1 + '.k' + str(mersize))
                cmd += ' -s ' + target3
                cmd += ' -o ' + prefix1 + '.exclude'
                cmd += ' -stats ' + prefix1 + '.exclude.stats'
		print >>STDERR, cmd
                iret = os.system(cmd)
                if iret != 0:
                    os.unlink(prefix1 + '.exclude.mcidx')
                    os.unlink(prefix1 + '.exclude.mcdat')
                    raise ("Failed to make exclude mers!\n")

            if os.path.exists(prefix1 + '.exclude.mcdat'):
                os.unlink(target3 + '.mcdat')
                os.unlink(target3 + '.mcidx')
            else:
                raise ("Failed to find exclude mers?\n")


        F = open(prefix1 + '.mask.done', 'w')
        F.close()

    if merylOnly == 1 : sys.exit(0)

    #  This is the segmented search routine.
    #  assemblyId1 is used as the "query" sequences
    #  assemblyId2 is used for the table
    segmentsList = findSegmentsList(numsegments,leaff,seq2,prefix1)

    print >>STDOUT, segmentsList

    matchesfile = segmented_search_routine(
        seatac,
        seq1, seq2,
        prefix1,
        tmpdir,
        mersize,
        minfill,
        maxgap,
        segmentsList,
        numthreads,
        segmentIDtorun,
        buildOnly,
        filtername,
        filteropts,
        )
    return matchesfile
# end def

class localExecutable :
    def __init__(self, name):
        if 1:
            self.name = name
        else:
            self.name = "%s.%d" % ( name, os.getpid())
            print >>STDERR, self.name
            cmd = "cp `which %s` ./%s" % (name, self.name)
            #CMM the "which" is UNIX-centric.
            #CMM the "./" is UNIX-centric.
            print >>STDERR,"cmd =", cmd
            iret = os.system(cmd)
            assert(iret == 0)
    def __del__(self):
        print >>STDERR, "I almost deleted ", self.name
        #os.system("rm -f %s" % self.name)
        #CMM the "rm -f" is UNIX-centric.
    def run_old(self,argline,inpfile,outfile):
        inpfile.seek(0)
        outfile.seek(0)
        inpfile.flush()
        outfile.flush()
        cmd = "%s %s %s %s" % (self.name,argline,inpfile.name,outfile.name)
        print >>STDERR,"cmd =", cmd
        iret = os.system(cmd)
        assert(iret == 0)
        inpfile.seek(0)
        outfile.seek(0)
    def run(self,argline,inpfile,outfile):
        cmd = "%s %s %s %s" % (self.name,argline,inpfile,outfile)
        print >>STDERR,"cmd =", cmd
        iret = os.system(cmd)
        assert(iret == 0)


def main(runName):

    print >>STDERR, "Start reading checkpoint ...."
    t0 = time.time()
    obj = AtacDriver(runName)
    # Record some debugging information:
    obj.globals["sys.platform"]=str(sys.platform)
    obj.globals["sys.path"]=str(sys.path)
    obj.globals["time.asctime()"]=time.asctime()
    #obj.globals["os.defpath"]=str(os.defpath)
    #obj.globals["sys.version"]=str(sys.version)
    #obj.globals["os.environ"]=str(os.environ)

    if obj.globals.has_key("samespecies"):
        if obj.globals["samespecies"] == "parms1":
            if not obj.globals.has_key("heavyChainsOn"):
                obj.globals["heavyChainsOn"]      = "1"
            if not obj.globals.has_key("matchExtenderOn"):
                obj.globals["matchExtenderOn"]    = "1"
            if not obj.globals.has_key("uniqueFilterOn"):
                obj.globals["uniqueFilterOn"]     = "1"
            if not obj.globals.has_key("fillIntraRunGapsOn"):
                obj.globals["fillIntraRunGapsOn"] = "1"
            if not obj.globals.has_key("numsegments"):
                obj.globals["numsegments"]        = "1"

    heavyProg = None
    matchExtenderProg = None
    if sys.platform == "cygwin" :
        bindir = obj.globals['bindir']  # path to the binaries (hack! hack!)
        heavyProg = localExecutable(os.path.join(bindir,"heavyChainsDriver.exe"))
        matchExtenderProg = localExecutable(os.path.join(bindir,"MatchExtender.exe"))
    else:
        bindir = obj.globals['bindir']  # path to the binaries (hack! hack!)
        print >>STDERR, "obj.globals = ", obj.globals
        print >>STDERR, "bindir = ", bindir
        heavyProg = localExecutable(os.path.join(bindir,"heavyChainsDriver"))
        matchExtenderProg = localExecutable(os.path.join(bindir,"MatchExtender"))
    
    if not obj.globals.has_key('bridir'):
        obj.globals['bridir'] = obj.globals['atacdir']

    t1 = time.time()
    print >>STDERR, "Read checkpoint in %d seconds." % (t1-t0)
    t0=t1

    # The following are required:
    assert(obj.globals.has_key('assemblyId1'))
    assert(obj.globals.has_key('assemblyId2'))
    assemblyId1 = obj.globals["assemblyId1"]
    assemblyId2 = obj.globals["assemblyId2"]

    # Read the assemblies.atai file's mapping of nicknames to datasource.
    genomeAliases = {}
    if obj.globals.has_key("genomeAliases"):
        print >>STDERR, "Reading aliases."
        gafile = file(obj.globals["genomeAliases"],"r")
        for line in gafile:
            if(line[0:1]=="S"):
                (linetype, nickname, source) = line.split()
                genomeAliases[nickname] = source
    else:
        print >>STDERR, "NOT READING ALIASES!"

    # Give the assemblies.atai priority over the atac file.

    print "Id1, Id2 = ", assemblyId1, assemblyId2

    if genomeAliases.has_key(assemblyId1):
        obj.globals["assemblyFilePrefix1"] = genomeAliases.get(assemblyId1)
        print >>STDERR, assemblyId1, genomeAliases.get(assemblyId1)
    if genomeAliases.has_key(assemblyId2):
        obj.globals["assemblyFilePrefix2"] = genomeAliases.get(assemblyId2)
        print >>STDERR, assemblyId2, genomeAliases.get(assemblyId2)
        
    assert(obj.globals.has_key("assemblyFilePrefix1"))
    assert(obj.globals.has_key("assemblyFilePrefix2"))
    assemblyFilePrefix1=obj.globals["assemblyFilePrefix1"]
    assemblyFilePrefix2=obj.globals["assemblyFilePrefix2"]

    print >>STDERR, "Looking for ", assemblyFilePrefix1+".idxStore"
    if(os.path.exists(assemblyFilePrefix1+".idxStore")):
        print >>STDERR, "Looking for ", assemblyFilePrefix1+".seqStore"
        assert(os.path.exists(assemblyFilePrefix1+".seqStore"))
    else:
        assert(obj.globals.has_key('assemblyFilePrefix1'))
        assert(os.path.exists(assemblyFilePrefix1))
        IdxStore.createIndexedFasta( assemblyFilePrefix1, assemblyId1, 1)
    # end if

    print >>STDERR, "Looking for ", assemblyFilePrefix2+".idxStore"
    if(os.path.exists(assemblyFilePrefix2+".idxStore")):
        print >>STDERR, "Looking for ", assemblyFilePrefix2+".seqStore"
        assert(os.path.exists(assemblyFilePrefix2+".seqStore"))
    else:
        assert(obj.globals.has_key('assemblyFilePrefix2'))
        assert(os.path.exists(assemblyFilePrefix2))
        IdxStore.createIndexedFasta( assemblyFilePrefix2, assemblyId2, 1)
    # end if

    # CMM Since we made a private copy of the assembly fasta files, we need to
    # clean up the private copies of the assemblies.

    if not obj.globals.has_key('matchesFile'):
        print >>STDERR, "We need to make the raw matches."
    
        seq1 = assemblyFilePrefix1
        assert(os.path.exists(seq1))
        seq2 = assemblyFilePrefix2
        assert(os.path.exists(seq2))

        # CMM " -merylonly " # only run the meryl components
        # CMM " -samespecies " # use magic values for same species
        # CMM " -crossspecies " # use guesses for different species

        obj.globals['matchesFile'] = makeRawMatches(
            assemblyId1,
            assemblyId2,
            seq1,
            seq2,
            obj.globals['meryldir'],  # path to the MERYL directory
            obj.globals['bridir'],    # path to the BRI directory
            obj.globals['bindir'],    # path to the binaries (hack!)
            int(obj.globals.get('rawMatchMerSize',20)),
            int(obj.globals.get('rawMatchMinFillSize',
                                obj.globals.get('rawMatchMerSize',20))),
            int(obj.globals.get('rawMatchMerMaxDegeneracy',1)),
            int(obj.globals.get('rawMatchAllowedSubstutionBlockSize',0)),
            obj.globals.get('seatac.numsegments',1), # defaults to 1
            obj.globals.get('seatac.numthreads',1),  # defaults to 1
            )

    print >>STDERR, "OK made raw matches ", obj.globals['matchesFile']

    if(obj.globals.has_key('matchesFile')):
        print >>STDERR, "Have matchesFile!"

        inpname = obj.globals['matchesFile']
        outname = runName+".ckpFirst"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            print >>STDERR, "Begin reading: " + inpname
            assert(os.path.exists(inpname))
            # print >>STDERR, "len(obj.matches)=",len(obj.matches)
            print >>STDERR, 'from ' + inpname + ' reading'
            inpfile = open(inpname,'r')
            outfile = MyFile.myfile()
            MatchRecord.convertBrianRecordFormat(inpfile, outfile,
                                                 assemblyId1,
                                                 assemblyId2,
                                                 )
            obj.matches = outfile
            obj.runs = MyFile.myfile()
            del(obj.globals['matchesFile'])
            print >>STDERR, "Finished reading: " + inpname

            print >>STDERR, "Converted matches in %d seconds." % (t1-t0)
            t0=t1
            obj.checkpoint(outname)
            t1 = time.time()
            print >>STDERR, "Wrote checkpoint in %d seconds." % (t1-t0)
            t0=t1
        else:
            print >>STDERR, "MAKE NEW obj!"
            obj = AtacDriver(outname)
    else:
        die("No matches file")

    if( obj.globals.has_key('heavyChainsOn') and obj.globals['heavyChainsOn']=="1"):
        # Set as default parameters:
        if(not obj.globals.has_key('heavyMaxJump')):
            obj.globals['heavyMaxJump'] = 100000
        if(not obj.globals.has_key('heavyMinFill')):
            obj.globals['heavyMinFill'] = 100
        inpname = runName + ".ckpBeforeHeavyChains"
        outname = runName + ".ckpAfterHeavyChains"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            obj.checkpoint(inpname)
            heavyProg.run(
                "-g /assemblyId1=%s -g /assemblyId2=%s -g /heavyMaxJump=%d -g /heavyMinFill=%d" % (
                obj.globals['assemblyId1'],
                obj.globals['assemblyId2'],
                int(obj.globals['heavyMaxJump']),
                float(obj.globals['heavyMinFill'])
                ),
                inpname, outname)
            obj = AtacDriver(outname)
        else:
            obj = AtacDriver(outname)

    elif (obj.globals.has_key('chainGlobalOn') and obj.globals['chainGlobalOn']=="1"):
        # /work/assembly/floreald/ASM/src/Ross/chain-global
        # /work/assembly/floreald/ASM/src/break-chains
        inpname = runName + ".ch"
        outname = inpname + ".M30.dp.runs"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            obj.checkpoint(inpname)
            cmd="chain-global %s -M 30 -p DP > %s.M30.dp 2> %s.M30.dp.errs" % (inpname,inpname,inpname)
	    #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            cmd="break-chains %s.M30.dp -D 0 -M 10 -p DPR | grep -v 'M r ' > %s.M30.dp.runs 2> %s.M30.dp.runs.errs" % (inpname, inpname, inpname)
            #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            obj = AtacDriver(outname)
        else:
            obj = AtacDriver(outname)


    elif (obj.globals.has_key('chainConsvOn') and obj.globals['chainConsvOn']=="1"):
        # /work/assembly/floreald/ASM/src/Ross/chain-consv
        # /work/assembly/floreald/ASM/src/break-chains
        inpname = runName + ".ch"
        outname = inpname + ".cons.runs"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            obj.checkpoint(inpname)
            cmd="chain-consv %s -p CS > %s.cons 2> %s.cons.errs" % (inpname,inpname,inpname)
            #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            cmd="break-chains %s.cons -diffrun -D 0 -M 10 -p TMP > %s.cons.runs.tmp 2> %s.cons.runs.tmp.errs" % (inpname,inpname,inpname)
	    #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            cmd="break-chains %s.cons.runs.tmp -D 0 -M 10 -p CSR | grep -v 'M r ' > %s.cons.runs 2> %s.cons.runs.errs" % (inpname,inpname,inpname)
       	    #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            obj = AtacDriver(outname)
        else:
            obj = AtacDriver(outname)

    elif (obj.globals.has_key('chainGreedyOn') and obj.globals['chainGreedyOn']=="1"):
        # /work/assembly/floreald/ASM/src/Ross/chain-greedy
        # /work/assembly/floreald/ASM/src/break_chains
        inpname = runName + ".ch"
        outname = inpname + ".greedy.runs"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            obj.checkpoint(inpname)
            cmd = "chain-greedy %s -p GR -M 10 -W 500 > %s.greedy 2> %s.greedy.errs" % (inpname, inpname, inpname)
	    #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            cmd = "break-chains %s.greedy -D 0 -M 10 -p GRR | grep -v 'M r ' > %s.greedy.runs 2> %s.greedy.runs.errs" % (inpname, inpname, inpname)
            #CMM The stdout redirection is a Unix-ism.
            print >>STDERR,"cmd= ", cmd
            iret = os.system(cmd)
            assert(iret == 0)
            obj = AtacDriver(outname)
        else:
            obj = AtacDriver(outname)
    else:
        print >>STDERR, "No filtering via chaining used."


    if (obj.globals.has_key('matchExtenderOn') and obj.globals['matchExtenderOn']=="1"):
        inpname = runName + ".ckpBeforeMatchExtender"
        outname = runName + ".ckpAfterMatchExtender"
        if(not( os.access(inpname,os.F_OK) and os.access(outname,os.F_OK))):
            obj.checkpoint(inpname)
            matchExtenderProg.run(
                "",
                inpname, outname)
            obj = AtacDriver(outname)
        else:
            obj = AtacDriver(outname)

            
    if(not obj.globals.has_key('rawMatchMerSize')):
        obj.globals['rawMatchMerSize'] = 20
    if(not obj.globals.has_key('rawMatchMerMaxDegeneracy')):
        obj.globals['rawMatchMerMaxDegeneracy'] = 1
    if(not obj.globals.has_key('rawMatchMinFillSize')):
        obj.globals['rawMatchMinSize'] = obj.globals['rawMatchMerSize']

    if(not obj.globals.has_key('globalMatchMinSize')):
        obj.globals['globalMatchMinSize'] = 2*int(obj.globals['rawMatchMerSize'])
        # Many 2*rawMatchMerSize-1 matches are due to isolated single
        # nucleotide mutations in otherwise perfect repeats.
    if(not obj.globals.has_key('globalPerfectRunMinLen')):
        obj.globals['globalPerfectRunMinLen'] = 100
    if(not obj.globals.has_key('globalPerfectRunMaxGapLen')):
        obj.globals['globalPerfectRunMaxGapLen'] = 100000

    if(not obj.globals.has_key('intraRunGapIsolatedMismatchLen')):
        obj.globals['intraRunGapIsolatedMismatchLen'] = 20

    print >>STDERR, "assemblyId1 = ", obj.globals['assemblyId1']

    obj.runOld()
    t1 = time.time()
    print >>STDERR, "Ran in %d seconds." % (t1-t0)
    t0=t1

    del obj.globals['bridir']
    del obj.globals['inpname']
    del obj.globals['outname']
    
    obj.checkpoint(runName+".ckpLast")
    t1 = time.time()
    print >>STDERR, "Wrote checkpoint in %d seconds." % (t1-t0)

if __name__ == '__main__':
    main(sys.argv[1])
