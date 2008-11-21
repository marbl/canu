#!/usr/local/bin/perl

my @imps;

while(<>){
 @w=split;
 if($w[0] eq "Seed"){
   if($N>0){
#     $xx=@imps;
     output_IUM();
   }
   $N++;
   undef @imps;
#   $yy=@imps;
   $iid=$w[1];
   $min=0;
   $max=0;
   $n=0;
   next;
 }
 if($w[0] eq "SEEN:"){
   next;
   $f=1;
 } else {
   $f=0;
 }
 $n++;
 if($w[1+$f]<$min){$min=$w[1+$f];}
 if($w[2+$f]>$max){$max=$w[2+$f];}
 if($w[3+$f] eq "<--"){
   $beg=$w[2+$f];
   $end=$w[1+$f];
 } else {
   $beg=$w[1+$f];
   $end=$w[2+$f];
 }
 $mid=$w[0+$f];
 $curr = new implayout($mid,$beg,$end);
 push(@imps,$curr);
}

if($N>0){
  output_IUM();
}

exit(0);


sub output_IUM(){
  if($n<2){return;}

#  $zz=@imps;
#  print "counts: $yy $xx $zz $q\n";
  $len=$max-$min;

  @sortimps = sort {  ($a->{beg} <$a->{end}? $a->{beg} : $a->{end}) <=>  ($b->{beg} <$b->{end}? $b->{beg} : $b->{end}) } @imps;

# also output the IFG records for each read used in this IUM
for($i=0;$i<$n;$i++){
print
"{IAF
acc:$sortimps[$i]->{mid}
typ:R
chi:0
cha:0
clr:-1,-1
mst:U
}\n";
}

print
"{IUM
acc:$iid
cov:0.000
mhp:0.000
sta:X
fur:X
len:$len
cns:
.
qlt:
.
for:0
nfr:$n\n";

for($i=0;$i<$n;$i++){
$sortimps[$i]->{beg}+=-$min;
$sortimps[$i]->{end}+=-$min;

print
"{IMP
typ:R
mid:$sortimps[$i]->{mid}
con:0
pid:0
pos:$sortimps[$i]->{beg},$sortimps[$i]->{end}
ahg:0
bhg:0
dln:0
del:
}\n";

}

print "}\n";

}

package implayout;
sub new {
    my $class = shift;
    my $self = {};

    $self->{mid} = undef;
    $self->{beg} = undef;
    $self->{end} = undef;

    my ($mid,$beg,$end) = @_;
    $self->{beg} = $beg;
    $self->{end} = $end;
    $self->{mid} = $mid;
    return bless $self, ref($class) || $class;
}
