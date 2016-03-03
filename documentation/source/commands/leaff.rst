leaff
~~~~~~

::

  usage: leaff [-f fasta-file] [options]
  
  SOURCE FILES
     -f file:     use sequence in 'file' (-F is also allowed for historical reasons)
     -A file:     read actions from 'file'
  
  SOURCE FILE EXAMINATION
     -d:          print the number of sequences in the fasta
     -i name:     print an index, labelling the source 'name'
  
  OUTPUT OPTIONS
     -6 <#>:      insert a newline every 60 letters
                    (if the next arg is a number, newlines are inserted every
                    n letters, e.g., -6 80.  Disable line breaks with -6 0,
                    or just don't use -6!)
     -e beg end:  Print only the bases from position 'beg' to position 'end'
                    (space based, relative to the FORWARD sequence!)  If
                    beg == end, then the entire sequence is printed.  It is an
                    error to specify beg > end, or beg > len, or end > len.
     -ends n      Print n bases from each end of the sequence.  One input
                    sequence generates two output sequences, with '_5' or '_3'
                    appended to the ID.  If 2n >= length of the sequence, the
                    sequence itself is printed, no ends are extracted (they
                    overlap).
     -C:          complement the sequences
     -H:          DON'T print the defline
     -h:          Use the next word as the defline ("-H -H" will reset to the
                    original defline
     -R:          reverse the sequences
     -u:          uppercase all bases
  
  SEQUENCE SELECTION
     -G n s l:    print n randomly generated sequences, 0 < s <= length <= l
     -L s l:      print all sequences such that s <= length < l
     -N l h:      print all sequences such that l <= % N composition < h
                    (NOTE 0.0 <= l < h < 100.0)
                    (NOTE that you cannot print sequences with 100% N
                     This is a useful bug).
     -q file:     print sequences from the seqid list in 'file'
     -r num:      print 'num' randomly picked sequences
     -s seqid:    print the single sequence 'seqid'
     -S f l:      print all the sequences from ID 'f' to 'l' (inclusive)
     -W:          print all sequences (do the whole file)
  
  LONGER HELP
     -help analysis
     -help examples
