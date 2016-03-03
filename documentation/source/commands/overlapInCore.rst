overlapInCore
~~~~~~

::

  * No kmer length supplied; -k needed!
  ERROR:  No output file name specified
  USAGE:  overlapInCore [options] <gkpStorePath>
  
  -b <fn>     in contig mode, specify the output file
  -c          contig mode.  Use 2 frag stores.  First is
              for reads; second is for contigs
  -G          do partial overlaps
  -h <range>  to specify fragments to put in hash table
              Implies LSF mode (no changes to frag store)
  -I          designate a file of frag iids to limit olaps to
              (Contig mode only)
  -k          if one or two digits, the length of a kmer, otherwise
              the filename containing a list of kmers to ignore in
              the hash table
  -l          specify the maximum number of overlaps per
              fragment-end per batch of fragments.
  -m          allow multiple overlaps per oriented fragment pair
  -M          specify memory size.  Valid values are '8GB', '4GB',
              '2GB', '1GB', '256MB'.  (Not for Contig mode)
  -o          specify output file name
  -P          write protoIO output (if not -G)
  -r <range>  specify old fragments to overlap
  -t <n>      use <n> parallel threads
  -u          allow only 1 overlap per oriented fragment pair
  -w          filter out overlaps with too many errors in a window
  -z          skip the hopeless check
  
  --maxerate <n>     only output overlaps with fraction <n> or less error (e.g., 0.06 == 6%)
  --minlength <n>    only output overlaps of <n> or more bases
  
  --hashbits n       Use n bits for the hash mask.
  --hashstrings n    Load at most n strings into the hash table at one time.
  --hashdatalen n    Load at most n bytes into the hash table at one time.
  --hashload f       Load to at most 0.0 < f < 1.0 capacity (default 0.7).
  
  --maxreadlen n     For batches with all short reads, pack bits differently to
                     process more reads per batch.
                       all reads must be shorter than n
                       --hashstrings limited to 2^(30-m)
                     Common values:
                       maxreadlen 2048->hashstrings  524288 (default)
                       maxreadlen  512->hashstrings 2097152
                       maxreadlen  128->hashstrings 8388608
  
  --readsperbatch n  Force batch size to n.
  --readsperthread n Force each thread to process n reads.
  
