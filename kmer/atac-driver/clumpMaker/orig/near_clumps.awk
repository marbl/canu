NF<12{
  xch[$1]=$2;
  ych[$1]=$3;
  x1[$1]=$4;
  y1[$1]=$5;
  x2[$1]=$6;
  y2[$1]=$7;
  ori[$1]=$8;
}
NF>=12{

  if($NF==-1){
  # if not in a clump, do lots of stuff
    o="";

    # test whether near the most-recently encountered clump
        # ie whether match is within parallelogram defined
	# by x-begin and -end of clump and high and low diagonals
	# of clump, with all values expanded by the fudge factor
     # four cases depending on orientation

    if(prev!=""){
      if($12==1){
	if(ori[prev]==1){
          # match and clump both forward
	  diag1=$10-$6;
	  diag2=$10+$11-($6+$7);
	  if(diag1>diag2){
	    tmp=diag1;
	    diag1=diag2;
	    diag2=diag1;
	  }
	} else {
# match forward, clump reversed
	  diag1=$10+$6;
	  diag2=$10+$11+($6+$7);
	  o="reversed";
	}
      } else {
	if(ori[prev]==1){
# match reversed, clump forward
	  diag2=($10+$11)-$6;
	  diag1=$10-($6+$7);
	  o="reversed";
	} else {
# both match and clump reversed
	  diag1=$10+$6+$7;
	  diag2=$10+$11+$6;
	  if(diag1>diag2){
	    tmp=diag1;
	    diag1=diag2;
	    diag2=diag1;
	  }
	}
      }

      if(xch[prev]==$5&&ych[prev]==$9&&
	 x2[prev]+fudge>$6&&x1[prev]-fudge<$6+$7&&
	 diag2>y1[prev]-fudge&&diag1<y2[prev]+fudge){
	if(o!=""){
	  print $0,"near", prev,o;
	} else {
	  print $0,"near", prev;
	}
	next;
      }
    }

    # test whether near the next clump (not yet encountered)
        # ie whether match is within parallelogram defined
	# by x-begin and -end of clump and high and low diagonals
	# of clump, with all values expanded by the fudge factor
     # four cases depending on orientation

    if(prev!=0){
      if($12==1){
	if(ori[prev-1]==1){
          # match and clump both forward
	  diag1=$10-$6;
	  diag2=$10+$11-($6+$7);
	  if(diag1>diag2){
	    tmp=diag1;
	    diag1=diag2;
	    diag2=diag1;
	  }
	} else {
# match forward, clump reversed
	  diag1=$10+$6;
	  diag2=$10+$11+($6+$7);
	  o="reversed";
	}
      } else {
	if(ori[prev-1]==1){
# match reversed, clump forward
	  diag2=($10+$11)-$6;
	  diag1=$10-($6+$7);
	  o="reversed";
	} else {
# both match and clump reversed
	  diag1=$10+$6+$7;
	  diag2=$10+$11+$6;
	  if(diag1>diag2){
	    tmp=diag1;
	    diag1=diag2;
	    diag2=diag1;
	  }
	}
      }

      if(xch[prev-1]==$5&&ych[prev-1]==$9&&
	 x2[prev-1]+fudge>$6&&x1[prev-1]-fudge<$6+$7&&
	 diag2>y1[prev-1]-fudge&&diag1<y2[prev-1]+fudge){
	if(o!=""){
	  print $0,"nearby", prev-1,o;
	} else {
	  print $0,"nearby", prev-1;
	}
	next;
      }
    }

    print;
  } else {
    # in a clump, so just output as is
    prev=$NF;
    print;
  }
}
