#!/bin/bash

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

set -e             #  -e          - exit on error
#set -x            #  -x          - echo each line as it is executed
set -o pipefail    #  -o pipefail - fail on errors in pipes


initialize_region() {
    region_type=`cat dnanexus-job.json | jq .region | sed s/,//g | sed s/\"//g | awk -F ":" '{print $1}'`

    if [ $region_type == "aws" ]; then
       export DNANEXUS_REGION="aws"

    elif [ $region_type == "azure" ]; then
       export DNANEXUS_REGION="azure"

    else
       echo "ERROR: region type '$region_type' is unknown, must be 'aws' or 'azure'."
       exit 1
    fi

    echo "."
    echo "INITIALIZE REGION to DNANEXUS_REGION='$DNANEXUS_REGION'"
    echo "."
}


#  Fetch a subscript and run it.
#  See Execution.pm; search for fetch_and_run.
fetch_and_run() {

    echo "."
    echo "FETCH_AND_RUN"
    echo "."

    echo "."
    echo "ENVIRONMENT"
    echo "."

    env

    echo "."
    echo "PARAMETERS"
    echo "."
    echo "output_folder         '${output_folder}'"
    echo "script_path           '${script_path}'"
    echo "script_name           '${script_name}'"
    echo "."
    echo "canu_iteration        '${canu_iteration}'"
    echo "canu_iteration_max    '${canu_iteration_max}'"
    echo "."
    echo "DX_ARRAY_ID           '${DX_ARRAY_ID}'"
    echo "DX_PROJECT_CONTEXT_ID '${DX_PROJECT_CONTEXT_ID}'"
    echo "."

    initialize_region

    #  Give up if the iteration is defined and larger than allowed.  Canu
    #  should be aborting before this job is ever submitted, and so this
    #  _probably_ does nothing.

    if [ "x$canu_iteration" != "x" -a $canu_iteration -gt $canu_iteration_max ]; then
        echo "Error: tried $canu_iteration times, giving up."
        exit 1
    fi

    #  If we're the executive, dump the iteration count into our 'defaults'
    #  file, then grab the script and run it.
    #
    #  Otherwise, we're a compute job.
    #  Fetch $output_folder/$script_path/$script_name into
    #                       $script_path/$script_name and run it.
    #
    #  Note that some $script_name are e.g., 'scripts/2-sort.sh',
    #  and we must make the 'scripts' directory before we can
    #  download the file.

    if [ $script_name = "canu-executive.sh" ] ; then
        canuPath=`which canu`
        canuPath=`dirname $canuPath`

        echo "canuIteration=$canu_iteration"        >  $canuPath/canu.defaults
        echo "canuIterationMax=$canu_iteration_max" >> $canuPath/canu.defaults

        echo "."
        echo "Fetch '$DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh'."
        echo "."

        dx download $DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh -o canu-executive.sh

        echo "."
        echo "Execute 'canu-executive.sh'."
        echo "."

        bash canu-executive.sh

        exit $?

    else
        echo "."
        echo "Fetch '$script_name' into '`dirname $script_path`/'."
        echo "."

        mkdir -p `dirname $script_path/$script_name`
        cd $script_path

        dx download $DX_PROJECT_CONTEXT_ID:$output_folder/$script_path/$script_name -o $script_name

        echo "."
        echo "Execute '$script_name' with DX_ARRAY_ID=$DX_ARRAY_ID."
        echo "."

        export DX_ARRAY_ID=$DX_ARRAY_ID
        bash $script_name

        exit $?
    fi
}




main() {
    echo "."
    echo "MAIN"
    echo "."

    echo "."
    echo "ENVIRONMENT"
    echo "."

    env

    echo "."
    echo "PARAMETERS"
    echo "."
    echo "read_type           '${read_type}'"
    echo "output_folder       '${output_folder}'"
    echo "output_prefix       '${output_prefix}'"
    echo "genome size         '${genome_size}'"
    echo "parameters          '${parameters}'"
    echo "."
    echo "hapa_name           '${hapa_name}'"
    echo "hapb_name           '${hapb_name}'"
    echo "."
    echo "project             '$DX_PROJECT_CONTEXT_ID'"
    echo "."

    echo "."
    echo "DNANEXUS-JOB.JSON"
    echo "."

    cat dnanexus-job.json

    initialize_region

    #  If output contigs exist, do nothing.

    if [ -s ${output_prefix}.contigs.fasta ]; then
      echo "Output '${output_prefix}.contigs.fasta' exists.  Stop."
      exit 0
    fi

    #  Parse the read_files into objects we can give canu.

    for ll in ${!read_files[@]} ; do
      rf=${read_files[$ll]}
      rn=${read_files_name[$ll]}
      rl=`dx-jobutil-parse-link "$rf"`

      echo "File $ll is '$rf' -> '$rl' = '$rn'."
      read_file_list="$read_file_list \"dnanexus:$rl=$rn\""
    done

    #  Parse hapa_files into objects we can give canu.

    for ll in ${!hapa_files[@]} ; do
      rf=${hapa_files[$ll]}
      rn=${hapa_files_name[$ll]}
      rl=`dx-jobutil-parse-link "$rf"`

      echo "File $ll is '$rf' -> '$rl' = '$rn'."
      hapa_file_list="$hapa_file_list \"dnanexus:$rl=$rn\""
    done

    if [ x$hapa_name != x ] ; then
        hapa_option="-haplotype${hapa_name} $hapa_file_list"
    fi

    #  Parse hapb_files into objects we can give canu.

    for ll in ${!hapb_files[@]} ; do
      rf=${hapb_files[$ll]}
      rn=${hapb_files_name[$ll]}
      rl=`dx-jobutil-parse-link "$rf"`

      echo "File $ll is '$rf' -> '$rl' = '$rn'."
      hapb_file_list="$hapb_file_list \"dnanexus:$rl=$rn\""
    done

    if [ x$hapb_name != x ] ; then
        hapb_option="-haplotype${hapb_name} $hapb_file_list"
    fi

    #  Build a nice script to run canu.

    echo  > canu-executive.sh "#!/bin/sh"
    echo >> canu-executive.sh ""
    echo >> canu-executive.sh "canu -haplotype -p ${output_prefix} -d . \\"
    echo >> canu-executive.sh "     executiveMemory=14 \\"             #  Linked to instanceType in
    echo >> canu-executive.sh "     executiveThreads=8 \\"             #  dxapp.json (14 and 8)
    echo >> canu-executive.sh "     objectStore=DNANEXUS \\"
    echo >> canu-executive.sh "     objectStoreClient=`which dx` \\"
    echo >> canu-executive.sh "     objectStoreClientUA=`which ua` \\"
    echo >> canu-executive.sh "     objectStoreClientDA=`which da` \\"
    echo >> canu-executive.sh "     objectStoreNameSpace=$output_folder \\"
    echo >> canu-executive.sh "     objectStoreProject=$DX_PROJECT_CONTEXT_ID \\"
    echo >> canu-executive.sh "     genomeSize=${genome_size} \\"
    echo >> canu-executive.sh "     $parameters \\"
    echo >> canu-executive.sh "     $hapa_option \\"
    echo >> canu-executive.sh "     $hapb_option \\"
    echo >> canu-executive.sh "     ${read_type} ${read_file_list}"
    echo >> canu-executive.sh ""
    echo >> canu-executive.sh "exit \$?"
    echo >> canu-executive.sh ""

    #  Remove any existing file before uploading, so we don't end up with
    #  multiple versions of the script.

    if dx describe --name $DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh > /dev/null 2>&1 ; then
        dx rm --recursive $DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh
    fi

    echo "."
    echo "Upload '$DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh'."
    echo "."

    dx upload --wait --parents --path $DX_PROJECT_CONTEXT_ID:$output_folder/canu-executive.sh canu-executive.sh

    echo "."
    echo "Execute 'canu-executive.sh'."
    echo "."

    bash canu-executive.sh

    exit $?
}
