for file in `ls *.rst`; do
  cmd=`echo $file |sed s/.rst//g`
  export PATH=$PATH:../../../Darwin-amd64/bin
  OUTPUT=`$cmd 2>&1 |awk '{print "  "$0}'`
  echo "$cmd"     > $file
  echo "~~~~~~"  >> $file
  echo ""        >> $file
  echo "::"      >> $file
  echo ""        >> $file
  echo "$OUTPUT" >> $file
done
