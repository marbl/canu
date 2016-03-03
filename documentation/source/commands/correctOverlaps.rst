correctOverlaps
~~~~~~

::

  ERROR: no input gatekeeper store (-G) supplied.
  ERROR: no input overlap store (-O) supplied.
  ERROR: no input read corrections file (-c) supplied.
  ERROR: no output erates file (-o) supplied.
  USAGE:  correctOverlaps [-d <dna-file>] [-o <ovl_file>] [-q <quality>]
              [-x <del_file>] [-F OlapFile] [-S OlapStore]
              [-c <cgb_file>] [-e <erate_file>
             <gkpStore> <CorrectFile> <lo> <hi>
  
  Recalculates overlaps for frags  <lo> .. <hi>  in
   <gkpStore>  using corrections in  <CorrectFile> 
  
  Options:
  -e <erate-file>  specifies binary file to dump corrected erates to
                   for later updating of olap store by  update-erates 
  -F             specify file of sorted overlaps to use (in the format
                 produced by  get-olaps
  -o <ovl_file>  specifies name of file to which OVL messages go
  -q <quality>   overlaps less than this error rate are
                 automatically output
  -S             specify the binary overlap store containing overlaps to use
