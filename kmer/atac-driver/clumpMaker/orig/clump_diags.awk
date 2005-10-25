$NF==-1{next}
$NF!=prev{
  if(NR>1){
    if(y1>y2){
      tmp=y1;
      y1=y2;
      y2=tmp;
    }
    print prev,xch,ych,x1,y1,x2,y2,ori;
  };
  prev=$NF;
  ori=$12;
  x1=$6;
  if(ori==1){
    y1=$10-x1;
  }else{
    y1=$10+x1;
  }
}
{
  x2=$6+$7;
  if(ori==1){
    y2=$10+$11-x2;
  } else {
    y2=$10+$11+x2;
  }
  xch=$5;
  ych=$9
}
END{
  if(y1>y2){
    tmp=y1;
    y1=y2;
    y2=tmp;
  };
  print prev,xch,ych,x1,y1,x2,y2,ori
}
