for file in `cat list |awk '{print $1}'`; do
   path=`echo $file |sed s/utgcns//g |sed s/libboost//g |sed s/boost//g`
   echo $path
   cp --parents /usr/local/BOOST/1_56_0/include/boost/$path ./
done
