use DBI;

$JCVI = 1;
$HELPTEXT =
  qq~Request a whole-genome shotgun assembly using Celera Assembler.

 usage: carun [options] <frg>
 
 inputs: 
  <frg>       Sequence file (.frg) generated e.g. by pullfrag.

 general options:
  -d <dir>          Use <dir> as the working directory. 
                    (default: /usr/local/aserver/var/assembly/<...> )
  -p <prefix>       Use <prefix> as the output prefix. (default: 'asm')
  -alias <a>        Identify the request as <a> on Assembly Server Console
  -maxCopy <dir>    Copy max output to <dir> (default: dir=maxCopy)
  -medCopy <dir>    Copy med output to <dir> (default: dir=medCopy) 
  -minCopy <dir>    Copy min output to <dir> (default: dir=minCopy)
  -noCopy           Do not make local copy. (default: medCopy)
  -[no]notify       Send email upon completion/error and create BITS case 
                    for errors. (default: notify)
  -s <specfile>     Read options from the specifications file <specfile>.
  -test             Run in debug mode (same as -D test -nonotify)  
  -version          Outputs the Celera Assembler Version
  -fields           Outputs the specfile fields
  -help	            This text
 
 CA options: 
  -e <percent>      Unitigger assumed percent error [0.00,0.06] (default: 0.015)
  -g <len>          Assign genome length for unitigger  
  -j <lbound>       Set scaffolder A-stat low bound (default: 1)
  -k <hbound>       Set scaffolder A-stat high bound (default: 5)
  -[no]ubs          Enable unitigger bubble smoothing (default: enabled)  
      
 Genomic Assembly Pipeline: 
 https://intranet.jcvi.org/cms/SE/GAP
 Tracking assembly requests:  
 http://assemblyconsole.tigr.org/
~;

1;
