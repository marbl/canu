#!/usr/bin/env python
import sys, os, copy, tempfile, cStringIO

# from __future__ import generators  # Necessary before Python 2.3

class myfile(file):
    "A temporary anonymous file"
    def __init__(self):
        filename = tempfile.mktemp()
        #print >>sys.stderr, "myfile: creating " + filename
        file.__init__(self,filename,"w+")
    def __del__(self):                 
        #print >>sys.stderr, "myfile: deleting " + self.name
        self.close()
        os.system("rm -f " + self.name)
    def link(self,othername):
        #print >>sys.stderr, "myfile: linking %s to %s" % ( self.name, othername)
        self.flush()
        os.system("ln -f %s %s" % (self.name, othername))

class ListLikeFileIter:
    # See http://www.python.org/peps/pep-0234.html
    # for file iterators.
    def __init__(self,filename):
        self._filename = filename
        self._fileptr = open(self._filename,"r")
        self._fileIter = iter(self._fileptr.readline,"")
    def __del__(self):
        self._fileptr.close()
    def next(self):
        line = self._fileIter.next()
        if line:
            return line
        else:
            raise StopInteration
        # end if
    def __getitem__(self,ii):
        # For files, the list location ii is ignored.
        # line = self._fileptr.readline()
        line = self._fileIter.next()
        if line:
            return line
        else:
            raise IndexError
        # end if
    # end def
    
class ListLikeFile:
    # See Mark Lutz, Programming Python, edition 1, page 18 and page 128.
    def __init__(self):
        #self._filename = tempfile.mktemp()
        #self._fileptr = open(self._filename,"w")
        self._fileptr = cStringIO.StringIO()
        #self._list = []
    def __del__(self):
        self._fileptr.close()
        #pass
    def __iter__(self):
        self._fileptr.flush()
        return iter(cStringIO.StringIO(self._fileptr.getvalue()))
        #return iter(self._fileptr)
        #return ListLikeFileIter(self._filename)
        return iter(self._list)
    def write(self,x):
        self._fileptr.write(x)
        #self._list.append(x)
    # end def
# end class

def tester():
    x = ListLikeFile()
    print >>x, 4
    print >>x, 5

    xi = iter(x)
    print "test 1i"
    for i in xi: print i,

    print >>x, 6
    print >>x, 7
    print "test 2i"
    for i in xi: print i,

    xj = iter(x)
    print "test 3j"
    for i in xj: print i,

    print >>x, 8
    print >>x, 9
    
    print "test 3j"
    for i in xj: print i,

    print "test 3i"
    for i in xi: print i,

    xk = iter(x)
    print "test 3k"
    for i in xk: print i,

    x = None

if __name__ == '__main__':
    tester()
