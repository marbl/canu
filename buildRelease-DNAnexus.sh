#!/bin/sh

mkdir -p src/dx-canu/resources/bin/

#  Generate a script to compile Canu using the Holy Build Box.

echo  > build-linux.sh  \#\!/bin/bash
echo >> build-linux.sh  yum install -y git
echo >> build-linux.sh  cd /build/src
echo >> build-linux.sh  gmake -j 12 \> ../Linux-amd64.out 2\>\&1
echo >> build-linux.sh  cd ..
#echo >> build-linux.sh  rm -rf Linux-amd64/obj
#echo >> build-linux.sh  tar -cf canu-linux.Linux-amd64.tar  canu-linux/README* canu-linux/Linux-amd64

chmod 755 build-linux.sh

echo Build Linux and make tarballs.

docker run \
  -v `pwd`:/build \
  -t \
  -i \
  --rm phusion/holy-build-box-64:latest /hbb_exe/activate-exec bash /build/build-linux.sh

rm -f build-linux.sh

tail -n 10 Linux-amd64.out

#  Fetch the Upload Agent and install in our bin/.

if [ ! -e src/dx-canu/resources/bin/ua ] ; then
  curl -L -R -O https://wiki.dnanexus.com/images/files/dnanexus-upload-agent-1.5.31-linux.tar.gz

  tar zxf dnanexus-upload-agent-1.5.31-linux.tar.gz

  mv dnanexus-upload-agent-1.5.31-linux/ua src/dx-canu/resources/bin/

  rm -rf dnanexus-upload-agent-1.5.31-linux.tar.gz
  rm -rf dnanexus-upload-agent-1.5.31-linux
fi

#  Package that up into the DNAnexus app.

rm -rf src/dx-canu/resources/usr/bin/
rm -rf src/dx-canu/resources/usr/lib/
rm -rf src/dx-canu/resources/usr/share/

mkdir -p src/dx-canu/resources/usr/bin/
mkdir -p src/dx-canu/resources/usr/lib/
mkdir -p src/dx-canu/resources/usr/share/

rsync -a Linux-amd64/bin/   src/dx-canu/resources/usr/bin/
rsync -a Linux-amd64/lib/   src/dx-canu/resources/usr/lib/
rsync -a Linux-amd64/share/ src/dx-canu/resources/usr/share/

#rm -fr Linux-amd64/obj Linux-amd64.tar Linux-amd64.tar.gz
#tar -cf Linux-amd64.tar Linux-amd64
#gzip -1v Linux-amd64.tar
#dx rm Linux-amd64.tar.gz
#dx upload Linux-amd64.tar.gz

cd src
dx build -f dx-canu

exit 0
