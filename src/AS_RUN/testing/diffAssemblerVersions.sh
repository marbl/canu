#!/usr/local/bin/bash

VERSION_1=${1}
VERSION_2=${2}

VERSION_1_Listing="/tmp/deleteme_VERSION_1_wgs-assembler_listing.txt"
VERSION_2_Listing="/tmp/deleteme_VERSION_2_wgs-assembler_listing.txt"
DIFFFile="/tmp/deleteme_VERSION_1_VERSION_2_wgs-assembler_diff.txt"

if [ -z ${VERSION_1} ] || [ ! -d ${VERSION_1} ] || [ -z ${VERSION_2} ] || [ ! -d ${VERSION_2} ] ; then
  echo "Please specify wgs-assembler version 1 root dir as parameter 1"
  echo "Please specify wgs-assembler version 2 root dir as parameter 2"
  return
fi

cwd=`pwd`

# get a directory listing of each to identify absent files in one or the other
cd ${VERSION_1}
find . -print | gawk 'BEGIN{FS="/"}{if($2=="doc"||$2=="example"||$2=="inc"||$2=="scripts"||$2=="src"||$2=="testcases")print $0}' | sort > ${VERSION_1_Listing}

cd ${cwd}
cd ${VERSION_2}
find . -print | gawk 'BEGIN{FS="/"}{if($2=="doc"||$2=="example"||$2=="inc"||$2=="scripts"||$2=="src"||$2=="testcases")print $0}' | sort > ${VERSION_2_Listing}

# compare the listings
diff ${VERSION_1_Listing} ${VERSION_2_Listing} > ${DIFFFile}

# "<" prefixed lines are files only in VERSION_1 directory
echo -e "\nFiles only in ${VERSION_1} directory:\n"
egrep "^< " ${DIFFFile} | cut -f 2 -d ' '

# ">" prefixed lines are files only in version 1 directory
echo -e "\nFiles only in ${VERSION_2} directory:\n"
egrep "^> " ${DIFFFile} | cut -f 2 -d ' '
rm -f ${DIFFFile}

echo -e "\n\nFile diffs:"
echo -e "  '>' is in ${VERSION_1}"
echo -e "  '<' is in ${VERSION_2}\n"

cd ${cwd}
for line in `comm -12 ${VERSION_1_Listing} ${VERSION_2_Listing} | cut -f 2- -d '/'` ; do
  cvsCheck=${line%/*}
  cvsCheck=${cvsCheck##*/}
  if [ -f ${VERSION_1}/${line} ] && [ ${cvsCheck} != "CVS" ] ; then
    delta=`diff -w -q ${VERSION_1}/${line} ${VERSION_2}/${line}`
    if [ ! -z "${delta}" ] ; then
      echo -e "DIFFERENCE: ${VERSION_1}/${line} vs ${VERSION_2}/${line}\n"
      diff -w ${VERSION_1}/${line} ${VERSION_2}/${line}
      echo -e "\n\n"
    fi
  fi
done

rm -f ${VERSION_1_Listing}
rm -f ${VERSION_2_Listing}

cd ${cwd}
