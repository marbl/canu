#!/usr/local/bin/perl 
# $Id: recall.pl,v 1.1 2007-07-03 19:53:59 moweis Exp $
#
# recall: Recall the consensus 
#
# Written by Marwan Oweis
# Designed by Martin Shumway
# The Institure of Genomics Research, 2007

# Program locations
my $BIN     = '/usr/local/bin';
my $COMMON  = '/usr/local/common';

# Needed applications
my %APPS = (
    GET     => "$BIN/getCoverage",
    TR      => "$BIN/trSlice",
    SPLIT   => "$BIN/splitSlice",
    SLCTG   => "$BIN/slice2contig",
    MERG    => "$BIN/mergeSlice",
    CTG2FST => "$COMMON/contig2fasta",
    SPTFST  => "$COMMON/splitfasta",
    NUCMER  => "$COMMON/nucmer",
    SHWALNS => "$COMMON/show-aligns",
    FNDSNPS => "$COMMON/findSNPs"
);
           
my @MY_DEPENDS = 
( 
    'TIGR::Foundation',
    values %APPS     
);

my $MY_VERSION = " 1.3 (Build " . (qw/$Revision: 1.1 $/)[1] . ")";
my $HELP = q~
  recall - Recall the consensus sequence of a contig.
  
  Synoposis: For each contig within the input contig file, it will recall the
  consensus using the quality data per a given ambiguity mode.  Optionally, 
  using the feature information, recall will not recall the consensus for the 
  portion of the contigs. 
  
  Usage:
    recall [options] [feat_filter] [coord_filter] <prefix>

  Input (in directory):
    <prefix>.contig
    <prefix>.qual    
    <prefix>.feat (required with -feat/-Feat option)
    
  Output:
    <prefix>.recall.contig
    <prefix>.recall.qc
    <prefix>.recall.feat (if -feat/-Feat specified)
    
             
  Options:
  --------
    -a <mode>    Ambiguity mode for recall (default: 3)
                 (0=No ambiguity, 1=Minimal Ambiguity, 2=Annotation, 3=Conic)
    -D <dbname>  Use the ambiguity mode for database 'dbname'
    -[no]metrics Generate metrics data file: <prefix>.recall.qc (default: enabled)
    -keep        Keep temp generated files.
                 Output in directory: 'recall/' (default: disabled)
    -o <prefix>  Use the given <prefix> for the output files.
                 
   Feature Filter Options (choose one):
    -----------------------------------
    -feat [tag]  Specify non-recalled features using comma-separated XML
                 feature names of <prefix>.feat file. 'nuu' may be used for the
                 non-unique unitigs tag.  
                 An update .feat file is generated with changes to contig 
                 positions.  This is the only behavior if no tag is provided.                
    -Feat [tag]  Opposite of -feat, the supplied feature names will be recalled 
                 the only features recalled.  See -feat for more info.
               
   Coordinate Filter Options (choose one):
   ---------------------------------------                 
    -cf <file>   Recall all excluding coordinates per file. (See below for file format) 
    -CF <file>   Opposite of -cf, only recall those coordinates per file. 
                 (See below for file format) 
        
   Coordinate File Format
   ----------------------
   A list of coordinates can be specified so as to include or exclude a range 
   from recall. The file format for the coordinates: 
     
        <file>      := <lines>
        <lines>     := <line><lines>|<line>
        <line>      := <contig_id> <partition>
        <partition> := <digit>,<partition>|<digit>
        <digit>     := Integers greater than or equal to 1
        
   The range is inclusive.
   Example coordinates.txt:      
      24353542 100-400,600-900,1050-$
      
      If given with the 'cf' flag, then Contigs 24353542 will be recalled 
      except for the coordinates 100-400 and 600-900 and 1050 until the end of 
      the contig. With the 'CF' flag, the opposite will be done: only 
      coordinates 100-400, 600-900 and 1050 until the end of the contig will be 
      recalled. 
     
                             
See also SliceTools info: 
http://intranet.tigr.org/software_docs/SliceTools.html
~;

use strict;
use bigint;
use TIGR::Foundation;
use IO::File;
use DBI;
use File::Copy;
use XML::Simple;                # To parse .feat file
use Statistics::Descriptive;    # For statistics

my $NO_AMIBIGUITY          = 0;
my $CHURCH_WATER_CODE      = 1; #a.k.a. minimal
my $ANNOTATION_CODE        = 2;
my $CONIC_CODE             = 3;
my $DEFAULT_AMBIGUITY_CODE = $CONIC_CODE; 
my $DEFAULT_OUTPUT_DIR     = 'recall';
my $DEFAULT_DEBUG_LEVEL    = 1;
my $DEFAULT_LOGFILE        = 'recall.log';
my $AMBIGUITIES_UPPER      = 'YRWSKMDVHB'; #List of Amibiguity Codes
my $AMBIGUITIES            = $AMBIGUITIES_UPPER . lc $AMBIGUITIES_UPPER;
my $NUU_TAG                = 'CA_NONUNIQUE_UNITIG';
my $EXCLUDE                = 0;  #for -feat, -cf options
my $INCLUDE                = 1;  #for -Feat, -CF option

my $STATE_NEW                      = 0;
my $STATE_FOUNDLEFT                = 1;
my $STATE_FOUNDRIGHT               = 2;

my %MASK_OPTS = 
    (
        0 => 'Exclude',
        1 => 'Include'
    );

my $dbh = undef;

#TIGR::Foundation connect
my $tf = new TIGR::Foundation() 
    or die("Can't instantiate TIGR::Foundation object.");
$tf->setHelpInfo($HELP);
$tf->addDependInfo(@MY_DEPENDS);
$tf->setVersionInfo($MY_VERSION);
$tf->setLogFile($DEFAULT_LOGFILE);

my $prefix          = undef;
my $outputPrefix    = undef;
my $ambCode         = undef;
my $debug           = undef;
my $metrics         = 1;
my $keep            = undef;
my @contigs         = ();
my %contigFeatPartition = (); # key: contigId, val: comma-seperated positions to partition the contig
my %startHash       = ();     # key: contigId, val: EXCLUDE - exclude the first range from recall
                              #                     INCLUDE - include the first range in recall
my $featureset      = undef;
my %contigToFeatHash = undef;
my %featHash        = ();
my $featMode        = undef;  #depends on -feat/-Feat options
my %coordHash       = ();
my $coorMode        = undef;  #depends on -feat/-Feat options
my $coord_file      = undef;
my $qualfile        = undef;
my $origcon         = undef;
my $featfile        = undef;
my $newcon          = undef;
my @contigFileList  = ();      #contig files to merge
my %partitionHash   = ();
 
sub runCommand ($) {
    my $cmd = shift;
    $tf->logLocal("running command: '$cmd'\n",2);
    my $bad = $tf->runCommand("$cmd");
    if ($bad) {
        $tf->bail( "The command '$cmd' was not successful, error code $bad");
    }    
}

sub validateRangeList($$) {
    my $contigId = shift;
    my $rangeList = shift;
    
    my @partitionPoints = split (/\,/,$rangeList);
    my $length = scalar @partitionPoints;
    if ( $length % 2 == 1 ) {
        $tf->bail ('Invalid coordinates for contig $contigId.  Partition coordinates must be in pairs (e.g. 400-500,1000-$)');
    }
    for ( my $i=1 ; $i < $length ; $i++ ) {
        if ( $partitionPoints[$i] eq '$' and $i < $length-1 ) {
            $tf->bail ('Invalid coordinates for contig $contigId.  The end of range character \'$\' may only be used as the last partition point.');
        } elsif  ( $partitionPoints[$i] eq '$' ) {
            next;
        }
        if ( $partitionPoints[$i] < $partitionPoints[$i-1] ) {
            $tf->bail ("Invalid coordinates for contig $contigId.  Partition points should be in incremental order: '" . $partitionPoints[$i-1] . ',' . $partitionPoints[$i] . "'")
        }    
    }
}

sub coordSetup() {
    
    my $coord_fh = new IO::File("<$coord_file") 
    	or $tf->bail("Cannot open $coord_file: $!");
    	
    my $line = undef;
    while ($line = <$coord_fh>) {
        chomp ($line);
        if ( $line =~ /^(\S+)\s+(\S+)/ ) {            
            my $contigId = $1;
            my $rangeList = $2;
            $rangeList =~ s/-/,/g;
            $rangeList =~ s/\$/inf/g;
            validateRangeList($contigId,$rangeList);
            $coordHash{$contigId} = $rangeList;            
        } elsif ( $line =~ /^\s*$/ ) {
            next;
        } else {
            $tf->bail("Invalid format in coordinates file '$coord_file': '$line'");
        }
    }
    
    close $coord_fh;        
}    


sub getOptions() {

   my @featArr = ();
   my @featArrExclude = ();
   my @featArrInclude = ();
   my $project = undef;
   my $coord_exclude_file = undef;
   my $coord_include_file = undef;
   
   my $result = $tf->TIGR_GetOptions(
        'a=i',      \$ambCode,
        'metrics!', \$metrics,
        'keep',     \$keep,
        'feat:s',   \@featArrExclude,
        'Feat:s',   \@featArrInclude,
        'D=s',      \$project,
        'cf=s',     \$coord_exclude_file,
        'CF=s',     \$coord_include_file,
        'o=s',      \$outputPrefix,
    );

    #checking the validity of the options and creating objects
    if ( ( !defined($result) ) || ( $result eq "" ) ) {
        $tf->bail(
            "The options could not be read. Please refer to recall usage via the -h option."
        );
    }


    if (defined $ambCode and defined $project )  {
        $tf->bail(
            "The -a and -D options are mutually exclusive.
            Please use one or the other (or neither), to specify the recall ambiguity mode."
        );
    } elsif ( defined $project) {
        $ambCode = getAmbiguityCode($project);
    } elsif ( !defined $ambCode) {  #if nothing defined
        $ambCode = $DEFAULT_AMBIGUITY_CODE;
    } else {  #$ambCode is defined
        $tf->bail(
            "Please choose a valid amibiguity mode (see help text).  '$ambCode' is not valid."
        ) if ( $ambCode ne $NO_AMIBIGUITY and
               $ambCode ne $CHURCH_WATER_CODE and 
               $ambCode ne $CONIC_CODE and 
               $ambCode ne $ANNOTATION_CODE );        
    }
        
    $tf->logLocal("Using ambcode = $ambCode\n",1);

    if ( scalar @featArrExclude > 0 and
         scalar @featArrInclude > 0 ) {
        $tf->bail(
            "The -feat and -Feat options are mutually exclusive, use one or the other (or neither)."
        );
    } elsif ( scalar @featArrExclude > 0 ) {
        @featArr = @featArrExclude;
        $featMode = $EXCLUDE;
    } elsif ( scalar @featArrInclude > 0 ) {
        @featArr = @featArrInclude;
        $featMode = $INCLUDE;
    }

    if ( defined $coord_exclude_file and 
         defined $coord_include_file ) {
        $tf->bail(
            "The -cf and -CF options are mutually exclusive, use one or the other (or neither)."
        );
    } elsif ( defined $coord_exclude_file ) {
        $coord_file = $coord_exclude_file;
        $coorMode = $EXCLUDE;
    } elsif ( defined $coord_include_file ) {
        $coord_file = $coord_include_file;
        $coorMode = $INCLUDE;
    }
    
    my $argvSize = scalar(@ARGV);   
    if ( $argvSize == 0 ) {
        $tf->bail("Please provide a genome assembly prefix.");
    } else {
        mkdir $DEFAULT_OUTPUT_DIR;    
        chdir $DEFAULT_OUTPUT_DIR;    
        $prefix = $ARGV[0];
    }
    
    
    #If features are defined
    if ( scalar @featArr > 0 ) {
        #In case features are provided comma separated
        @featArr = split(/,/,join(',',@featArr));
    
        #recall takes the codeword nuu or NUU to ignore NUU features 
        for ( my $index = 0 ; $index <= $#featArr ; $index++) {    
            if (uc $featArr[$index] eq 'NUU') {
                $featArr[$index] = $NUU_TAG;
            }
        }
    
        %featHash = map {$_, 1} @featArr;
    }
    
    # Update debug level for logging if it is not specified by the user
    if ( !defined $tf->getDebugLevel() ) {
        $tf->setDebugLevel($DEFAULT_DEBUG_LEVEL);
    }
    
    $debug = $tf->getDebugLevel();
    
    $tf->logLocal("Prefix: $prefix",2);
}


sub checkInputFiles() {

    $qualfile = "$prefix" . ".qual";
    $origcon = "$prefix" . ".contig";
    $featfile = "$prefix" . ".feat";
    #Required Files
    if ( !-r "../$origcon" ) {
        $tf->bail("Missing required file: $origcon");             
    } elsif ( !-r "../$qualfile" ){
        $tf->bail("Missing required file: $origcon");             
    } else {
        symlink("../$qualfile","$qualfile");
        symlink("../$origcon","$origcon");
    }    
    
    #If -feat option is defined then, a .feat file is needed to determine where
    #the feature regions are.
    if ( defined %featHash and !-r "../$featfile" ) {
        $tf->bail("For '-feat' option,  a $prefix.feat file is required.");
    } elsif ( defined %featHash ) {
        symlink("../$featfile","$featfile");
    }    

    #If -cf or -CF options are defined
    if ( defined $coorMode and !-r "../$coord_file" ) {
        $tf->bail("Unable to read coordinates file: $coord_file.");
    } elsif ( defined $coorMode ) {
        symlink("../$coord_file","$coord_file");
    }    
}

sub getTotalConsensusCalls($) {
    
    my $contigFile = shift;
    
    runCommand("grep bases $contigFile | awk '{print \$3}' > $contigFile.bases");
    
    my $bases_fh = new IO::File "<$contigFile.bases"
    	or $tf->bail("Cannot open $contigFile.bases: $!");    
    my @bases = <$bases_fh>;
    close $bases_fh;
    
    my $totalBases = 0;
    
    foreach (@bases) {
    	$totalBases += $_;
    }
    
    return $totalBases;
}

sub getAmbiguityCode($) {
    my $prefix = shift;
    
    my $ambCode = undef;
    init_db_connection() unless defined $dbh;
    
    my $query = "select ambiguity_method from genomes where db = '$prefix'";

    $tf->logLocal( "Executing the query $query", 3 );

    my $qh = $dbh->prepare($query)
      or bail("Cannot prepare $query: " . $dbh->errstr );

    defined $qh->execute()
      or bail("Database query \'$query\' failed: " . $dbh->errstr );

    my @row             = $qh->fetchrow();
    my $ambCode = $row[0];
    $tf->logLocal( "DB Ambiguity Code=$ambCode", 3 );
    $qh->finish();
    $qh = undef;
    
    $ambCode = $CONIC_CODE if ($ambCode eq 'conic');
    $ambCode = $CHURCH_WATER_CODE if ($ambCode eq 'minimal');
    $ambCode = $ANNOTATION_CODE if ($ambCode eq 'annotation');
    
    $ambCode = $DEFAULT_AMBIGUITY_CODE unless (defined $ambCode and $ambCode ne '');
    
    return $ambCode;    
}

#Name:   init_db_connection
#Input:  none
#Output: none
#Usage:  This function initializes a connection to the asdb database.
sub init_db_connection() {
    $tf->logLocal( "Establishing a connection to the database", 3 );

    # Try to connect to the database server
    $dbh = DBI->connect(
        "dbi:Sybase:server=SYBTIGR;packetSize=8092",
        'access', 'access',
        {
            PrintError => 0,
            RaiseError => 0,
            AutoCommit => 1
        }
    );
    $dbh->do("use common")
      or $tf->bail("Failed to open database common'");
    $dbh->{InactiveDestroy} = 1;
}

my $xs = undef;
sub featSetup($) {
    my $featfile = shift;
    
    my $totalPartitionSize = 0;    
                        
    $tf->logLocal("Features to filter: ".join(',',keys %featHash),2);
    $tf->logLocal("Features filter mode: ".$MASK_OPTS{$featMode},2);
    
    #Read <prefix>.feat file, extract Features
    $xs = new XML::Simple(forcearray => 1, searchpath => ".");
    $featureset = eval { $xs->XMLin($featfile) };
    if ($@) {
      $tf->bail("ERROR: Cannot read XML input file, $featfile: $@");
    }

    #For each contig with Features, 
    #check to see any of the features are to be filtered
    foreach my $contig (@{$featureset->{Contig}}) {
        my %fiveToThrees = (); # key: end5 position, val: end3 position
        $contigToFeatHash{$contig->{Id}} = $contig->{Feature};
        
            print "Contig $contig->{Id} has differences\n";
            
            print "Has following features:\n";
            
            foreach my $feat (@{$contigToFeatHash{$contig->{Id}}}) {
                my $location  = (@{$feat->{Location}})[0];
                #The End5 and End3 elements are the ungapped positions
                #(positions not counting gaps) and are 1-based.
                my $end5 = $location->{End5} <= $location->{End3} ? $location->{End5}: $location->{End3};
                my $end3 = $location->{End5} >= $location->{End3} ? $location->{End5}: $location->{End3};
                                
                print $feat->{Class} . ' ' . $feat->{Type} . ' ' . $end5 . ' ' . $end3 . "\n";
            }
        
        foreach my $feature (@{$contig->{Feature}}) {
            if ( defined $featHash{$feature->{Type}} ) { #Contig has Feature to Filter
                my $location  = (@{$feature->{Location}})[0];
                #The End5 and End3 elements are the ungapped positions
                #(positions not counting gaps) and are 1-based.
                my $end5 = $location->{End5};
                my $end3 = $location->{End3};
                
                #Make end5 the smaller value of the two values.
                if ($end3 < $end5 ) {
                    my $tmp = $end3;
                    $end3 = $end5;
                    $end5 = $tmp;
                }
                if ( defined $fiveToThrees{$end5} && $fiveToThrees{$end5} > $end3) {
                    $tf->logWarning("Feature found but inclosed within another Feature");
                }
                else {                
                    $fiveToThrees{$end5} = $end3;
                }                                             
            }
        }
        
        #If the contig has features to filter        
        if ( scalar keys( %fiveToThrees ) > 0 ) {
            #Create a mask of the partitions needed of the slice file
            #Basically, suppose we have the following set of features to filter 
            #for a contig:
            # (5',3') : (800,1000), (2000,3000), (700,900), (900,1800), (2100,2500)
            #The result of the mask would be: (700,1800), (2000,3000)
            #
            #NOTE: For these pos values, the 3' position is not inclusive, since
            #the pos values correspond to the spaces not the bases. Meaning, that a
            #filtered feature at position (500,550), 
            #this referring to the bases numbered from 500 to 549 (1-based).          
            my @maskedFives = ();
            my @maskedThrees = ();
            my $count = 0;
            foreach my $newFive (sort {$a <=> $b} keys %fiveToThrees) {
                if ( $count == 0 || $newFive > $maskedThrees[$count-1]) {
                    $maskedFives[$count] = $newFive;
                    $maskedThrees[$count] = $fiveToThrees{$newFive};
                    $count++;
                }
                elsif ( $newFive < $maskedThrees[$count-1]) {
                    if ( $fiveToThrees{$newFive} > $maskedThrees[$count-1] ) {
                        $maskedThrees[$count-1] = $fiveToThrees{$newFive};
                    }
                }            
            }
             
            #Creating the partition string to be passed to splitSlice       
            #If for this contig, there were no NRFs then this loop will not run,
            #since 'maskedFives' will be empty.
            my $partition = undef;
            my $contigPartitionSize = 0;
            for( $count = 0 ; $count < scalar(@maskedFives) ; $count++ ) {
                my $fivePrime  = $maskedFives[$count]; 
                my $threePrime = $maskedThrees[$count]; #3' is not inclusive
                $partition .= $fivePrime . ',' . $threePrime;
                $contigPartitionSize += $threePrime - $fivePrime; 
                if ( $count < scalar(@maskedFives) - 1) {
                     $partition .= ',';
                }
            }
                            
            $contigFeatPartition{$contig->{Id}} = $partition;
            $totalPartitionSize += $contigPartitionSize;
            
            #Get the first position of the contigparition. 
            #If it is 1, then
            #the NRF starts at the beginning of the contig
            #otherwise, it starts somewhere later.
            #
            #This info is needed, since that when the files are split
            #using splitSlice we know that the output files will alternate
            #between slices of a NRF and those without.
            #So once we know whether the first file is a NRF or not,
            #we will know the remaining files for this contig.
            #
            #For example,
            #Suppose for contig 12345 of length 10, there is a NRF at [5,8) 
            #We would then use splitSlice, to partition at 5 and 8, so there 
            #would be three files outputted:
            #   12345_0.slice (positions 1-4)
            #   12345_1.slice (positions 5-7)
            #   12345_2.slice (positions 8-9)
            #By knowing that the first position is not 1, we know that the first slice
            #file will be a non NRF, so it should be sent to trSlice, then the second file
            #is NRF, and should not be sent, etc...    
            my @positions = split(/,/, $partition);
            #0 - means that even files will have NRFs, 1-means that odd files will have NRFs
            $tf->logLocal( "Contig: $contig->{Id} has partitions: $contigFeatPartition{$contig->{Id}}",2);
        }
        else {
            $tf->logLocal( "Contig: $contig->{Id} does not have any NRFs",2);
        }
    }    
    return $totalPartitionSize;    
}

sub max($$) {
    my $left = shift;
    my $right = shift;
    
    return $left if ( $left >= $right);
    return $right;
}

sub min($$) {
    my $left = shift;
    my $right = shift;
    
    return $left if ( $left <= $right);
    return $right;
}

sub leftEnd($) {
    my $input = shift;    
    my @inputArr = split(/ /,$input);    
    return $inputArr[0];
}

sub rightEnd($) {
    my $input = shift;    
    my @inputArr = split(/ /,$input);    
    return $inputArr[1];
}

sub excludeExclude ($$) {
    my $input1_ref = shift;
    my $input2_ref = shift;
    my @pairs1 = @$input1_ref;
    my @pairs2 = @$input2_ref;
    my $pairs1length = scalar @pairs1;
    my $pairs2length = scalar @pairs2;
        
    my @partitionArray = ();
    my $first = 0;
    my $second = 0;
    my $state = $STATE_NEW;
    my $currentRight;
    while ( $first < $pairs1length or $second < $pairs2length ) {
        if ( $state == $STATE_FOUNDRIGHT ) {
             if ( $currentRight eq 'first' ) {
                if ( $second == $pairs2length ) {
                    last if ( rightEnd($pairs1[$first]) == inf);
                    push @partitionArray, rightEnd($pairs1[$first]);
                    $first++;
                    $state = $STATE_NEW;
                } elsif ( rightEnd($pairs2[$second]) <= rightEnd($pairs1[$first]) ) {
                    last if ( rightEnd($pairs1[$first]) == inf);
                    $second++;                                    
                } elsif ( leftEnd($pairs2[$second]) <= rightEnd($pairs1[$first]) ) {
                    $first++;
                    $currentRight = 'second';
                } else {
                    push @partitionArray, rightEnd($pairs1[$first]);
                    $first++;
                    $state = $STATE_NEW;                                          
                }
             }
             elsif ( $currentRight eq 'second' ) {
                if ( $first == $pairs1length ) {
                    last if ( rightEnd($pairs2[$second]) == inf);
                    push @partitionArray, rightEnd($pairs2[$second]);
                    $second++;
                    $state = $STATE_NEW;
                } elsif ( rightEnd($pairs1[$first]) <= rightEnd($pairs2[$second])) {
                    last if ( rightEnd($pairs2[$second]) == inf);
                    $first++;                                                    
                } elsif ( leftEnd($pairs1[$first]) <= rightEnd($pairs2[$second]) ) {
                    $second++;
                    $currentRight = 'first';
                } else {
                    push @partitionArray, rightEnd($pairs2[$second]);
                    $second++;
                    $state = $STATE_NEW;                                          
                }
             }
        }                    
        elsif ( $state == $STATE_FOUNDLEFT ) {
            if ( $second == $pairs2length ) {
                last if ( rightEnd($pairs1[$first]) == inf);
                push @partitionArray, rightEnd($pairs1[$first]);
                $first++;
                $state = $STATE_NEW;                     
            } elsif ( $first == $pairs1length) {
                last if ( rightEnd($pairs2[$second]) == inf);
                push @partitionArray, rightEnd($pairs2[$second]);
                $second++;                                    
                $state = $STATE_NEW;                     
            } elsif ( rightEnd($pairs1[$first]) > rightEnd($pairs2[$second]) ) {
                last if ( rightEnd($pairs1[$first]) == inf);
                $currentRight = 'first';
                $second++;
                $state = $STATE_FOUNDRIGHT;             
            } else {
                last if ( rightEnd($pairs2[$second]) == inf);
                $currentRight = 'second';
                $first++;
                $state = $STATE_FOUNDRIGHT;
            }   
        } elsif ( $state == $STATE_NEW )  {
            if ( $second == $pairs2length ) {
                push @partitionArray, leftEnd($pairs1[$first]);
                $state = $STATE_FOUNDLEFT;                     
            } elsif ( $first == $pairs1length) {
                push @partitionArray, leftEnd($pairs2[$second]);
                $state = $STATE_FOUNDLEFT;                     
            } elsif ( leftEnd($pairs1[$first]) < leftEnd($pairs2[$second]) ) {
                push @partitionArray, leftEnd($pairs1[$first]);
                if ( rightEnd($pairs1[$first]) <= leftEnd($pairs2[$second]) ) {
                    push @partitionArray, rightEnd($pairs1[$first]);
                    $first++;
                } else {
                    $state = $STATE_FOUNDLEFT;
                }
            } else {
                push @partitionArray, leftEnd($pairs2[$second]);
                if ( rightEnd($pairs2[$second]) <= leftEnd($pairs1[$first]) ) {
                    push @partitionArray, rightEnd($pairs2[$second]);
                    $second++;
                } else {
                    $state = $STATE_FOUNDLEFT;
                }
            }                           
        }
    }
    return @partitionArray;
}

sub includeInclude ($$) {
    my $input1_ref = shift;
    my $input2_ref = shift;
    my @pairs1 = @$input1_ref;
    my @pairs2 = @$input2_ref;
    my $pairs1length = scalar @pairs1;
    my $pairs2length = scalar @pairs2;
    
    my @partitionArray = ();
    my $first = 0;
    my $second = 0;
    while ( $first < $pairs1length or $second < $pairs2length ) {
        if ( $second == $pairs2length or
             $first == $pairs1length) {
            last;
        } elsif ( rightEnd($pairs1[$first]) <= leftEnd($pairs2[$second]) ) {                
            $first++;
        } elsif ( rightEnd($pairs2[$second]) <= leftEnd($pairs1[$first]) ) {                
            $second++;
        } else {
            push @partitionArray, max(leftEnd($pairs1[$first]),leftEnd($pairs2[$second]));
            if ( min(rightEnd($pairs1[$first]),rightEnd($pairs2[$second])) == inf ) {
                last;
            }
            else {
                push @partitionArray, min(rightEnd($pairs1[$first]),rightEnd($pairs2[$second]));
            }
            $first++;
            $second++;
        }                           
    }
    return @partitionArray;
}

sub excludeInclude($$) {
    my $input1_ref = shift;
    my $input2_ref = shift;
    my @excludePairs = @$input1_ref;
    my @includePairs = @$input2_ref;
    my $excludelength = scalar @excludePairs;
    my $includelength = scalar @includePairs;
    
    my @partitionArray = ();
    my $exclude = 0;
    my $include = 0;
    my $found = 0;
    my $excludeRightEndPoint = undef;
    my $check = undef;
    while ( $exclude < $excludelength or $include < $includelength ) {
        if ( $found == 0 ) {
            ##############
            ##  exclude:  ...
            ##  include:  [...] [...]      
            ##############
            if ( $exclude == $excludelength ) {
                push @partitionArray, leftEnd($includePairs[$include]);
                last if ( rightEnd($includePairs[$include]) == inf);
                push @partitionArray, rightEnd($includePairs[$include]);
                $include++;            
            } 
            ##############
            ##  exclude:  [...] [...]
            ##  include:  -     
            ##############            
            elsif ( $include == $includelength ) {
                last;
            } 
            ##############
            ##  exclude:       [...]
            ##  include:  [...]      
            ##############                        
            elsif ( rightEnd($includePairs[$include]) <  leftEnd($excludePairs[$exclude])) {
                push @partitionArray, leftEnd($includePairs[$include]);
                push @partitionArray, rightEnd($includePairs[$include]);
                $include++;                    
            } 
            ##############
            ##  exclude:  [...] 
            ##  include:        [...]     
            ##############                                    
            elsif ( leftEnd($includePairs[$include]) >  rightEnd($excludePairs[$exclude]) ) {
                $exclude++;
            } 
            ##############
            ##  exclude:    [... 
            ##  include:  [...     
            ##############                                                
            elsif ( leftEnd($includePairs[$include]) < leftEnd($excludePairs[$exclude]) ) {
                push @partitionArray, leftEnd($includePairs[$include]);
                push @partitionArray, leftEnd($excludePairs[$exclude]);
                ##############
                ##  exclude:    [...] 
                ##  include:  [...]     
                ##############                                                
                if ( rightEnd($includePairs[$include]) <= rightEnd($excludePairs[$exclude]) ) {
                    $include++;
                } else {
                    $found = 1;
                    $excludeRightEndPoint = rightEnd($excludePairs[$exclude]);
                    $exclude++;
                }
                
            }
            ##############
            ##  exclude:  [... 
            ##  include:    [...        
            ##############                                                
            elsif ( leftEnd($excludePairs[$exclude]) <  leftEnd($includePairs[$include]) ) {
                if ( rightEnd($includePairs[$include]) <= rightEnd($excludePairs[$exclude]) ) {
                    $include++;
                } else {
                    $found = 1;
                    $excludeRightEndPoint = rightEnd($excludePairs[$exclude]);
                    $exclude++;                
                }        
            }             
        }
        elsif ( $found == 1 ) {  #leftEnd($excludePairs[$exclude]) >  leftEnd($includePairs[$include])
            if ( $exclude == $excludelength ) {
                push @partitionArray, $excludeRightEndPoint;
                last if ( rightEnd($includePairs[$include]) == inf);
                push @partitionArray, rightEnd($includePairs[$include]);
                $include++;
                $found = 0;            
            } 
            elsif ( $include == $includelength ) {
                last;
            } 
            elsif ( leftEnd($excludePairs[$exclude]) <= rightEnd($includePairs[$include]) ) {
                push @partitionArray, $excludeRightEndPoint;
                push @partitionArray, leftEnd($excludePairs[$exclude]);
                if ( rightEnd($includePairs[$include]) <= rightEnd($excludePairs[$exclude]) ) {
                    $found = 0;
                    $include++;
                } else {
                    $excludeRightEndPoint = rightEnd($excludePairs[$exclude]);
                    $exclude++;
                }                            
            } else {
                push @partitionArray, $excludeRightEndPoint;
                last if ( rightEnd($includePairs[$include]) == inf);
                push @partitionArray, rightEnd($includePairs[$include]);
                $include++;
                $found = 0;            
            }
        } 
    }
    return @partitionArray;
}
sub convergeFilters() {
    foreach my $contig (@contigs) {
         if ( exists $contigFeatPartition{$contig} or exists $coordHash{$contig}) {
            $tf->logLocal("Contig: $contig has the following feat range ($MASK_OPTS{$featMode}): $contigFeatPartition{$contig}\n",2) if ( exists $contigFeatPartition{$contig});
            $tf->logLocal("Contig: $contig has the following coord range ($MASK_OPTS{$coorMode}): $coordHash{$contig}\n",2) if ( exists $coordHash{$contig});
            
            my @feats = split (/,/,$contigFeatPartition{$contig});
            my @coords = split (/,/,$coordHash{$contig});
            my @featPairs = ();
            my @coordPairs = ();
            for ( my $i = 0 ; $i < scalar @feats; $i+=2 ) {
                push @featPairs, "$feats[$i] $feats[$i+1]";
            }
            for ( my $i = 0 ; $i < scalar @coords; $i+=2 ) {
                push @coordPairs, "$coords[$i] $coords[$i+1]";
            }

            my $featlength = scalar @featPairs;
            my $coordlength = scalar @coordPairs;
            if ( !defined $featMode and defined $coorMode ) {
                $partitionHash{$contig} = [split(/ /,join(' ',@coordPairs))];
                if ( leftEnd($coordPairs[0]) > 1 ) {
                    $startHash{$contig} = ($coorMode == $EXCLUDE)? $INCLUDE : $EXCLUDE;                    
                }  else {
                    $startHash{$contig} = $coorMode;
                }
            } elsif ( defined $featMode and !defined $coorMode ) {
                $partitionHash{$contig} = [split(/ /,join(' ',@featPairs))];
                if ( leftEnd($featPairs[0]) > 1 ) {
                    $startHash{$contig} = ($featMode == $EXCLUDE)? $INCLUDE : $EXCLUDE;                    
                }  else {
                    $startHash{$contig} = $featMode;
                }
            } elsif ( $featMode == $EXCLUDE and $coorMode == $EXCLUDE ) {  # -feat -cf
                $partitionHash{$contig} = [excludeExclude(\@featPairs,\@coordPairs)]; 
                if ( min(leftEnd($featPairs[0]),leftEnd($coordPairs[0])) > 1 ) {
                    $startHash{$contig} = $INCLUDE;
                }  else {
                    $startHash{$contig} = $EXCLUDE;                    
                }
            } elsif ( $featMode == $INCLUDE and $coorMode == $INCLUDE ) {  # -Feat -CF
                $partitionHash{$contig} = [includeInclude(\@featPairs,\@coordPairs)];
                if ( max(leftEnd($featPairs[0]),leftEnd($coordPairs[0])) > 1 ) {
                    $startHash{$contig} = $EXCLUDE;                    
                }  else {
                    $startHash{$contig} = $INCLUDE;
                }
            } elsif ( $featMode == $EXCLUDE and $coorMode == $INCLUDE ) {  # -feat -CF
                $partitionHash{$contig} = [excludeInclude(\@featPairs,\@coordPairs)];                
                if ( @{$partitionHash{$contig}}[0] > 1 ) {
                    $startHash{$contig} = $EXCLUDE;                    
                }  else {
                    $startHash{$contig} = $INCLUDE;
                }
            } elsif ( $featMode == $INCLUDE and $coorMode == $EXCLUDE ) {  # -Feat -cf
                $partitionHash{$contig} = [excludeInclude(\@coordPairs,\@featPairs)];
                if ( @{$partitionHash{$contig}}[0] > 1 ) {
                    $startHash{$contig} = $EXCLUDE;                    
                }  else {
                    $startHash{$contig} = $INCLUDE;
                }
            }
            pop @{$partitionHash{$contig}} if ( ${$partitionHash{$contig}}[$#{$partitionHash{$contig}}] == inf);
            $tf->logLocal("Final Partition Points = " . join(',',@{$partitionHash{$contig}}) . ", start with $MASK_OPTS{$startHash{$contig}}",2);
         }                      
    }
}


sub getCoverage() {
    #Get coverage is called with the following options:
    #   -1 (split files are 1-based)
    #   -C (Circular contigs)
    #   -l -S (output one file per contig)
    my $getopts = "-1 -C -l -S $prefix > getCoverage.log";
    my $getcmd = "$APPS{GET} $getopts";    
    runCommand($getcmd);
}

sub readContigs($) {
    my $origcon = shift;
    #Open the contig file to read in the contig IDs
	my $contig_fh = new IO::File "<$origcon"
	  or $tf->bail("Cannot open $origcon: $!");    
    @contigs = ();
    while (<$contig_fh>) {
        if ( /^##(\S+)\s*/ ) {
            my $contig = $1;
            chomp($contig);
            push(@contigs,$contig);
        }
    }
    close $contig_fh;
}

sub trSlice($$) {
    my $inputSlice = shift;
    my $outputSlice = shift;
    my $tropts = "-i 1 -g 1 -x 1 -y 1 -c -a $ambCode -o $outputSlice $inputSlice";
    my $trcmd = "$APPS{TR} $tropts";            
    runCommand($trcmd);            
}

sub recallContig ($) {
    my $contig = shift;
    
    #if the contig has associated Partition means
    #that it has non-recalled Features, then split the file
    if ( defined $contigFeatPartition{$contig}) {
        $tf->logLocal("Running $APPS{SPLIT} for contig with non-recalled features",3);
        my $splitopts = "-p $contigFeatPartition{$contig} $contig.slice";
        my $splitcmd = "$APPS{SPLIT} $splitopts";        
        runCommand($splitcmd);             

        my $partitionCounter = 0;
        my $state = $startHash{$contig};
        #splitSlice outputs files in the format $contig_$count.slice
        #We'll loop until none exist
        my $recallSliceFiles = '';            
        while ( -e "${contig}_${partitionCounter}.slice") {
            my $inputSliceFile = "${contig}_${partitionCounter}.slice"; 
            my $outputSliceFile = "${contig}_${partitionCounter}.recall.slice"; 
            if ( $state == $INCLUDE )  {
                $tf->logLocal("${contig}_${partitionCounter}.slice will be recalled",3);
                trSlice($inputSliceFile,$outputSliceFile);
                $state = $EXCLUDE;
            }
            #If Exclude, then file is just copied over to the new format (no recall)
            else {
                $tf->logLocal("${contig}_${partitionCounter}.slice will not be recalled",3);
                copy "$inputSliceFile", "$outputSliceFile"; 
                $state = $INCLUDE;
            }
            $recallSliceFiles .= " $outputSliceFile ";
            $partitionCounter++; 
        }
        
        #Merge the outputted slice files using mergeSlice
        my $mergeopts = "$recallSliceFiles > $contig.recall.slice";
        my $mergecmd = "$APPS{MERG} $mergeopts";
        runCommand($mergecmd);          
    }
    #No non-recalled features exist for this contig, so run it as is (no splitSlice or mergeSlice needed)
    else {
        trSlice("$contig.slice","$contig.recall.slice");
    }
    
}

sub convertSliceToContig($) {
    my $contig = shift;
    #Convert the contig slice file to contig format
    my $tocontigopts = "-i $contig -C ${contig}.recall.slice >> ${contig}.recall.contig";
    my $tocontigcmd = "$APPS{SLCTG} $tocontigopts";        
    runCommand($tocontigcmd);
    push @contigFileList, "${contig}.recall.contig";
}

sub generateOutput() {
    #Cat all the individual contig contig files to one contig file
    $newcon = "$prefix" . ".recall.contig";
    my $contigFileStr = join (' ',@contigFileList);     
    runCommand("cat $contigFileStr > $newcon");        
}

sub recallCleanup() {
    #Remove all intermediary files
    rename("$newcon","$newcon.tmp");
    rename("$origcon","$origcon.tmp");
    runCommand("rm *.contig");        
    runCommand("rm *.slice");
    runCommand("rm $qualfile");
    runCommand("rm $featfile") if ( defined %featHash );    
    #Rename output contig 
    rename("$newcon.tmp","$newcon");
    rename("$origcon.tmp","$origcon");
}

sub publishOutput() {
    $outputPrefix = "$prefix" if ( !defined $outputPrefix);
    copy "$prefix.recall.contig", "../$outputPrefix.recall.contig";
    copy "$prefix.recall.qc", "../$outputPrefix.recall.qc" if ($metrics);
    copy "$prefix.recall.feat", "../$outputPrefix.recall.feat" if (defined %featHash);
}

sub getMetrics($$) { 
    my $totalPartitionSize = shift;
    my $keep = shift;
        
    #Get the total number of bases recalled
    my $totalConsensusCalls = getTotalConsensusCalls($origcon);
    $totalConsensusCalls -= $totalPartitionSize if ( defined %featHash );
    
    my $baseDifferences = Statistics::Descriptive::Full->new();
    #For each contig, create a fasta file, and split the file to get 
    #a seq file for each contig
    foreach my $contigFile ($origcon, $newcon) {
        runCommand("$APPS{CTG2FST} $contigFile > $contigFile.fasta");
        mkdir "${contigFile}Fastas";
        chdir "${contigFile}Fastas";            
        runCommand("$APPS{SPTFST} ../$contigFile.fasta > splitfasta.$contigFile.log");
        chdir '..';
    }
    
    #Get the total number of Ambiguity codes
    runCommand("grep -v '>' $origcon.fasta > origConsensus.txt");
    
    #Open the consensus
    my $origCons_fh = new IO::File "<origConsensus.txt"
      or $tf->bail("Cannot open origConsensus.txt: $!");    
    my @origCons = <$origCons_fh>;
    close $origCons_fh;
    my $consensus = join ('',@origCons);
    uc $consensus;
    my $totalAmbiguities = $consensus =~ tr/YRWSKMDVHB//;
            
    #Create directory to hold the differences of each contig before and after recall
    mkdir 'seqDiffs';            
    chdir 'seqDiffs';
    my %contigDiffs   = (); #key: contig Id that has differences, value: number of bases that differ
    
    my $totalDiffs    = 0;
    my $totalAmbIntro = 0;
    my $totalAmbReslv = 0;
    my $totalNonAmbs  = 0;
    
    #For each contig, see if there are differences due to recall, if so,
    #-Run nucmer to prodcuce the delta file
    #-Use show-aligns to create an alignment file per the deltas
    #-Feed the alignment file to findSNPS, to list the differences (one per line),
    # This is outputted as a .snps file.
    #-Create qc file to indicate the changes 
    foreach my $contig (@contigs) {
        
        #system is used since TIGR::Foundations runCommand interrupts 
        #a 1 exit status as an error (and therefore exists). 
        #'diff' however returns a '1' in the event that the two input 
        #files differ, which is not an error condition. 
        system("diff ../${origcon}Fastas/$contig.seq ../${newcon}Fastas/$contig.seq > $contig.diff");
        #Check if outputted diff file has size > 0 (i.e. there are differences)
        if ( -s "$contig.diff" > 0 ) {
            print "Contig $contig has differences\n";
            
            print "Has following features:\n";
            
            my %featArrHash = ();
            foreach my $feat (@{$contigToFeatHash{$contig}}) {
                $featArrHash{$feat->{Name}} = $feat;
                my $location  = (@{$feat->{Location}})[0];
                #The End5 and End3 elements are the ungapped positions
                #(positions not counting gaps) and are 1-based.
                my $end5 = $location->{End5};
                my $end3 = $location->{End3};                
                print $feat->{Class} . ' ' . $feat->{Type} . ' ' . $end5 . ' ' . $end3 . "\n";
            }
             
            runCommand("$APPS{NUCMER} ../${origcon}Fastas/$contig.seq ../${newcon}Fastas/$contig.seq -p $contig > nucmer.log 2>&1");
            runCommand("$APPS{SHWALNS} $contig.delta $contig $contig > $contig.aligns");
            runCommand("$APPS{FNDSNPS} $contig.aligns > $contig.snps");
            #Open the snps to read the number of lines (i.e. number of differences)
        	my $snps_fh = new IO::File "<$contig.snps"
        	  or $tf->bail("Cannot open $contig.snps: $!");
        	#Count the number of lines (i.e. the number of differences)    
            my @snpDiffs = <$snps_fh>;
            close $snps_fh;
    
            my $numBaseDifferences = 0;     
            my $ambiguitiesResolved = 0;           
            my $ambiguitiesIntroduced = 0;
            my $notAmbiguity = 0;
            my $currentIncrease = 0;
            print "Checking differences:\n";           
            foreach my $snp (@snpDiffs) {
                my @snpInfo  = split (/,/,$snp);
                next if ( $snpInfo[4] != 0);  #Only check the first chain (0) (the best match)
                my $origBase = uc $snpInfo[15];
                my $newBase  = uc $snpInfo[16];
                my $length1 = length $origBase;
                my $length2 = length $newBase;
                
                #TODO if ( $length1 != $length2) then warning!!!! and ?skip it?                 
                if ( $origBase =~ /\.+/ && $newBase =~ /\.+/ ) {
                    #Indiciates a change from (for example) T to t (t: T or Gap)
                    #Not considered a change
                    next;
                }
    
                $numBaseDifferences++;
                
                #If there are amibiguities in the new seq file
                if ( $origBase =~ /[ACTG]/ && $newBase =~ /[$AMBIGUITIES]/) {
                    $ambiguitiesIntroduced++;
                }
                #If there are amibiguities in the old seq file
                elsif ( $origBase =~ /[$AMBIGUITIES]/ && $newBase =~ /[ATCG]/) {
                    $ambiguitiesResolved++;
                }
                else {
                    $notAmbiguity++;
                }
                
                ###################
                # UPDATE Feat file
                ###################
                #check for updates to features
                #col 15-Polymorphism type: I=insertion in query sequence, D=deletion from query sequence                
                my $polyType = uc $snpInfo[14];
                my $diffStartPos  = $snpInfo[10];  #Start index of the difference
                $polyType =~ s/^\s+//;
                $polyType =~ s/\s+$//;
                print "Difference: '$polyType'\n";           
                
                if ( $polyType eq 'I' ) { #Insert
                    print "Insertion for contig $contig\n";        
                    foreach my $feat (values %featArrHash) {
                        my $location  = (@{$feat->{Location}})[0];
                        #The End5 and End3 elements are the ungapped positions
                        #(positions not counting gaps) and are 1-based.
                        my $end5 = $location->{End5};
                        my $end3 = $location->{End3};                
                        if ( $end3 < $diffStartPos ) {
                            next;
                        } elsif ( $end5 > $diffStartPos ) {
                            print "Insertion (1) for contig $contig and feat " . $feat->{Name} . "\n";
                            $location->{End5} += $length1;
                            $location->{End3} += $length1;
                        } else {
                            print "Insertion (2) for contig $contig and feat " . $feat->{Name} . "\n";
                            $location->{End3} += $length1;
                        }
                    }                    
                } elsif ( $polyType eq 'D' ) { #Delete
                    print "Deletion for contig $contig\n";        
                    foreach my $feat (keys %featArrHash) {
                        my $location  = (@{$feat->{Location}})[0];
                        #The End5 and End3 elements are the ungapped positions
                        #(positions not counting gaps) and are 1-based.
                        my $end5 = $location->{End5};
                        my $end3 = $location->{End3};                
                        if ( $end3 < $diffStartPos ) {
                            next;
                        } elsif ( $end5 > $diffStartPos ) {
                            print "Deletion (1) for contig $contig and feat " . $feat->{Name} . "\n";
                            $location->{End5} -= $length1;
                            $location->{End3} -= $length1;
                        } else {
                            print "Deletion (2) for contig $contig and feat " . $feat->{Name} . "\n";
                            $location->{End3} -= $length1;
                        }
                    }                    
                } else {    #No change to features
                    print "Nothing '$polyType' for contig $contig\n";        
                    next;
                }
                             
            }
    
            #Maintain hash of each contig with differences (zero diff contigs not kept in hash)
            $contigDiffs{$contig} = $numBaseDifferences;
    
            $totalDiffs += $numBaseDifferences;
            $totalAmbIntro += $ambiguitiesIntroduced;
            $totalAmbReslv += $ambiguitiesResolved;
            $totalNonAmbs  += $notAmbiguity;
            
            #Total bases differed
            $baseDifferences->add_data($numBaseDifferences);                
        }
    }

    if ( defined %featHash ) {
        #output new feat file
        my $outputFeat = $xs->XMLout($featureset,xmldecl => "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<!DOCTYPE FeatureSet SYSTEM \"http://aserver.tigr.org/aserver/feature.dtd\">");
        my $feat_fn = "../$prefix.recall.feat";
        my $feat_fh = new IO::File ">$feat_fn"
          or $tf->bail("Cannot open $feat_fn: $!");
        $feat_fh->print($outputFeat);
        close $feat_fh;        
    }    
    
        
    chdir '..';
    if ( !defined $keep ) {
        runCommand("rm -rf ${origcon}Fastas ${newcon}Fastas seqDiffs $origcon.fasta $newcon.fasta");
    }
    
    #Create qc file
    my $qc_fn = "$prefix.recall.qc";
    my $qc_fh = new IO::File ">$qc_fn"
      or $tf->bail("Cannot open $qc_fn: $!");        
    my $count = $baseDifferences->count();
    my $sum   = $baseDifferences->sum();
    my $mean  = $baseDifferences->mean();
    my $stdDv = $baseDifferences->standard_deviation();
    my $perConsensusCallsChanged = $totalDiffs / $totalConsensusCalls;
    my $numOfContigs = scalar @contigs;
    my $perContigsRecalled = $count/$numOfContigs;
     
    $qc_fh->print ("[RecallStats]\n");        
    $qc_fh->print ("TotalConsensusCallsConsidered=$totalConsensusCalls\n");
    $qc_fh->print ("TotalNRFCallsNotConsidered=$totalPartitionSize\n");
    $qc_fh->print ("TotalConsensusCallChanges=$totalDiffs\n");
    $qc_fh->printf ("PercentConsensusCallChanges=%0.2f\n",$perConsensusCallsChanged);
    $qc_fh->print ("TotalAmibiguityConsensusCalls=$totalAmbiguities\n");        
    $qc_fh->print ("TotalAmibiguitiesIntroduced=$totalAmbIntro\n");
    $qc_fh->print ("TotalAmibiguitiesResolved=$totalAmbReslv\n");
    $qc_fh->print ("TotalNumberContigs=$numOfContigs\n");
    $qc_fh->print ("TotalNumberContigsRecalled=$count\n");
    $qc_fh->printf ("PercentContigsRecalled=%0.2f\n",$perContigsRecalled);        
    $qc_fh->printf ("MeanBaseDifferencePerDifferingContig=%0.2f\n",$mean);        
    $qc_fh->printf ("StdDevDifferencePerDifferingContig=%0.2f\n",$stdDv);        
    close $qc_fh;
}

sub finalCleanup() {
    chdir '..';
    move 'recall/recall.log', 'recall.log'; 
    runCommand("rm -rf recall");
}

MAIN: 
{

    ### ######################################
    ### Parse program options
    ### ######################################     

    getOptions();

    ### ######################################
    ### Check required files for execution
    ### ######################################     

    checkInputFiles();
    readContigs($origcon);  #list of contigs

    ### ######################################
    ### 0) Coordinates Setup 
    ### (only if $coorMode is defined)
    ### ######################################
         
    print "Coordinates Setup ..." ;
    if ( defined $coorMode ) {
        coordSetup();    
        print " completed\n";
    } else {
        print " skipped\n";
    }         

    ### ######################################
    ### 1) Feature Setup 
    ### (only if %featHash is defined)
    ### ######################################
         
    print "Feature Setup ..." ;
    my $totalPartitionSize = undef;  # Count of the number of bases in the non recalled regions.
                                 # Used to determine the number of bases recalled.
                                 # If all regions are recalled, then this will stay as 0.
                         
    if ( defined %featHash ) {
        $totalPartitionSize = featSetup($featfile);
        print " completed\n";
    } else {
        print " skipped\n";
    }    

    ### ######################################
    ### 1) Converge Filters
    ### ######################################

    if ( scalar keys %contigFeatPartition > 0 or 
         scalar keys %coordHash > 0 ) {
             convergeFilters();
    } 
    
    ### ######################################
    ### 2) Run getCoverage
    ### ######################################     

    print "Generate slice files ..." ;
    getCoverage();
    print " completed\n";
        
    ### ######################################
    ### 3) Recall each contig
    ### ######################################
     
    print "Recall each contig ..." ;
    foreach my $contig (@contigs) {

        ### ######################################
        ### a) Recall the contig, taking into account filter options
        ### ######################################            

        recallContig($contig);

        ### ######################################
        ### b) Convert generated slice file to contig
        ### ######################################                    

        convertSliceToContig($contig);

    }
    print " completed\n";
    
    ### ######################################
    ### 4) Merge generated files to one output
    ### ######################################     

    print "Generate output ..." ;
    generateOutput();
    print " completed\n";
    
    ### ######################################
    ### 5) Cleanup temoprary files (unless keep requested)
    ### ######################################     

    print "Initial Cleanup ..." ;
    if ( !defined $keep ) {
        recallCleanup();
        print " completed\n";
    } else {
        print " skipped\n";
    }         
    
    ### ######################################
    ### 7) Generate metrics (if requested)
    ### ######################################         

    print "Generate metrics ..." ;
    if ( $metrics ) {
        getMetrics($totalPartitionSize,$keep);
        print " completed\n";
    } else {
        print " skipped\n";
    }         
    
    ### ######################################
    ### 6) Publish outputted files (.recall.contig, .recall.qc)
    ### ######################################     

    print "Publish output ..." ;
    publishOutput();
    print " completed\n";
    
    ### ######################################
    ### 8) Final cleanup (unless keep requested)
    ### ######################################         

    print "Final cleanup ..." ;
    if ( !defined $keep ) {
        finalCleanup();
        print " completed\n";
    } else {
        print " skipped\n";
    }         
}
