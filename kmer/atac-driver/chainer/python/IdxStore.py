import os, sys, DNA

#######################################################
# Begin class methods
#######################################################

def createIndexedFasta( prefix, nickname, makeseqstore):
  __doc__ = "Usage prefix, nickname, makeseqstore"
  # This is a class method (as opposed to an object method).
  # This method creates an indexed FASTA file on disk.

  print >>sys.stderr, "Opening %s to write %s.seqStore and %s.idxStore" % (prefix, prefix, prefix, )

  the_uid = None
  defline = None
  seqline = None
  linenumber = 0
  cur_offset = 0

  FASTA    = file( prefix, "r")
  IDXSTORE = file( prefix + ".idxStore", "w")
  if(makeseqstore):
    SEQSTORE = file(prefix + ".seqStore", "w")
  for line in FASTA:
    linenumber += 1
    line = line.strip()
    if(line[0:1] == ">"):
      # Clear current data to make space for new data.
      if( the_uid != None ):
        # The first time thru the_uid is equal to None.
        assert(defline != None)
        assert(seqline != None)

        # If we are using a database, then this might be the place to register the data.
        # uid2defline[the_uid] = defline
        # uid2seqline[the_uid] = seqline

        def_offset = cur_offset
        def_length = len(defline)
        if(makeseqstore): print >>SEQSTORE, defline
        # Store the FASTA defline.
        cur_offset += def_length + 1
        # remember the UNIX newline inserted by print.
        seq_offset = cur_offset
        seq_length = len(seqline)
        if(makeseqstore): print >>SEQSTORE, seqline
        # Store the FASTA sequence.
        cur_offset += seq_length + 1;
        # remember the UNIX newline inserted by print.

        print >>IDXSTORE, the_uid, def_length, seq_length, def_offset, seq_offset
      # end if

      # Now process the new data.
      defline = line
      the_uid = line.split()[0][1:]
      seqline = ""   # Clear any accumulated sequence.
    else:
      seqline += line  # Accumulate more DNA sequence
    # end if
  # end for

  if(the_uid != None):
    # Now make sure that accumulated data makes it to disk.

    # If we are using a database, then this might be the place to register the data.
    # uid2defline[the_uid] = defline
    # uid2seqline[the_uid] = seqline

    def_offset = cur_offset
    def_length = len(defline)
    if(makeseqstore): print >>SEQSTORE, defline
    # Store the FASTA defline.
    cur_offset += def_length + 1;
    # remember the UNIX newline inserted by print.
    seq_offset = cur_offset
    seq_length = len(seqline)
    if(makeseqstore): print >>SEQSTORE, seqline
    # Store the FASTA sequence.
    cur_offset += seq_length + 1;
    # remember the UNIX newline inserted by print.

    print >>IDXSTORE, the_uid, def_length, seq_length, def_offset, seq_offset
  # end if

  FASTA.close()
  IDXSTORE.close()
  if(makeseqstore):
    SEQSTORE.close()
  # end if
# end def

#######################################################
# End class methods
#######################################################

class IdxStore:
  __doc__ = "Class for fast access to multiFASTA files."
  
  #######################################################
  # Begin instance methods
  #######################################################

  def __init__(self,prefix,*optargs):
    __doc__ = "Create an instance of the class"
    if(optargs):
        self.nickname = optargs[0]
    else:
        self.nickname = None
    # end if
    self.uid2iid = {}  # declare an empty mapping
    self.iid2uid = []  # empty
    self.iid2def_length = []
    self.iid2seq_length = []
    self.iid2def_offset = []
    self.iid2seq_offset = []

    filename = prefix + ".idxStore"
    sys.stderr.write( "Openning index file %s for reading\n" % filename);
    idxstore = file(filename, "r");
    # assert(fp != NULL);

    # print dir(self)
    the_iid=0;
    while 1:
      line = idxstore.readline();
      if not line: break
      # sys.stderr.write("idxline %s\n" % line)
      cols = line.split()
      the_uid=cols[0];
      def_length=int(cols[1]); seq_length=int(cols[2]);
      def_offset=eval(cols[3]); seq_offset=eval(cols[4])
      if(self.nickname):
          self.uid2iid[self.nickname + ':' + str(the_iid)] = the_iid; # hashed
      self.uid2iid[the_uid] = the_iid; # hashed
      self.iid2uid.append(the_uid); # vector
      self.iid2def_length.append(def_length) # vector
      self.iid2seq_length.append(seq_length)
      self.iid2def_offset.append(def_offset)
      self.iid2seq_offset.append(seq_offset)
      the_iid += 1;

    #close(idxstore);
    sys.stderr.write("Read index idxStore\n");
    # print self.uid2iid

    filename = prefix + ".seqStore"
    message = "Openning sequence file %s for reading.\n" % filename
    sys.stderr.write(message)
    self.seqstore = file(filename, "r");
    # assert(seqstore != NULL);
    return
  

  def getStringFromFasta(self, forward, scaf_uid, start, length ):
    try:
      scaf_iid = self.uid2iid[scaf_uid]
    except KeyError:
      scaf_iid = eval(scaf_uid)
      if(scaf_iid < 1000000000):
        # sys.stderr.write("Using scaf_uid as the index\n")
        pass
      else:
        sys.stderr.write("scaf_uid=<%s> is invalid.\n" % scaf_uid)
        return ""

    seq_length = self.iid2seq_length[scaf_iid];
    seq_offset = self.iid2seq_offset[scaf_iid];
    # print >>sys.stderr, "seq_length, seq_offset, start =", seq_length, seq_offset, start
    # print >>sys.stderr, "seek to offset =", seq_offset+start
    self.seqstore.seek(seq_offset+start, 0); # from the beginning of file
    substring = self.seqstore.read(length)

    if(not forward):
      # sys.stderr.write("Taking reversecomplement\n")
      try:
          substring = DNA.DNA(substring).reversecomplement()
      except KeyError:
          sys.stderr.write("KeyError in DNA.DNA.reversecomplement()\n")
          sys.stderr.write("The query %d %s %d %d\n" % (forward,scaf_uid,start,length))
          sys.stderr.write("%s\n" % substring)
    #else:
      #sys.stderr.write("Leave as is\n")
    return substring
# end class

def convertIndexToUID (  x_prefix, y_prefix, inpfile, outname, assemblyId1, assemblyId2 ):
    outfile = myfile()
    
    DefLines = file(x_prefix, 'r')
    the_x_uid = {} # Declare an empty dictionary
    ii = 0
    for line in DefLines:
        # A valid idxStore format
        (ga_uid, sln, cln, sst, cst) = line.split()
        the_x_uid[assemblyId1+":"+str(ii)] = ga_uid
        ii += 1
    # end for
    DefLines.close()

    DefLines = file(y_prefix, 'r')
    the_y_uid = {} # Declare an empty dictionary
    ii = 0
    for line in DefLine:
        # A valid idxStore format
        (ga_uid, sln, cln, sst, cst) = line.split()
        the_y_uid[assemblyId2+":"+str(ii)] = ga_uid
        ii += 1
    # end for
    DefLines.close()

    inpfile.seek(0)
    for line in inpfile:
        if(line[0] == 'M'):
            FM = MatchRecord.MatchRecord(line)
            FM.x_scaf_uid = the_x_uid[FM.x_scaf_uid]
            FM.y_scaf_uid = the_y_uid[FM.y_scaf_uid]
            print >>outfile, FM
        # end if
    # end for
    outfile.finished()
    return outfile
# end def

#if __name__ == '__main__':
    # main(sys.argv[1],sys.argv[2])
    #main()
