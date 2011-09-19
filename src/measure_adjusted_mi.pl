#!/usr/bin/perl 
# $Id: measure_adjusted_mi.pl,v 1.2 2002/08/14 13:46:18 clark Exp $
# This file takes three arguments 
# measure_mi.pl  gold clusters maxcount 
# gold is a file in format word tag nl
# clusters is the output of the other
# we measure the mi of one with the other 
# BUT we adjust it by dividing by the count
# so it is not biased by the frequent words

$corpus = $ARGV[0];
$clusters = $ARGV[1];
if (scalar(@ARGV) > 2){
    $max = $ARGV[2];
}
else {
    $max = 100000000;
}
if (scalar(@ARGV) > 3){
    $min = $ARGV[3];
}
else {
    $min = -1;
}

print("max is $max, min is $min\n");
$numberLines = 0;
$numberGoodLines = 0;
open(CLUSTERS,"<$clusters") or die "Couldnt open file";

open(COUNT,"<$corpus") or die $!;
while($lineone = <COUNT>){
    $numberLines++;
    if ($lineone =~ m/^\s*([^\s]+)\s+([^\s]+)\s*$/){
	$numberGoodLines++;
	$count{$1}++;
    }
    else {
	#print("Skipping $lineone");
    }
}
print("Read $numberGoodLines good lines  out of $numberLines total\n");
close(COUNT);
open(CORPUS,"<$corpus") or die $!;

while(<CLUSTERS>){
    if (m/^(.+) (.+) (.+)$/){
        #print("Word $1 has tag $2\n");
        $dict{$1} = $2;
    }
}

# now read the corpus
while ($lineone = <CORPUS>){
    if ($lineone =~ m/^\s*([^\s]+)\s+([^\s]+)\s*$/){
	$wordone = $1;
	$tagone = $2;
	$tagtwo = $dict{$wordone};
	$n = $count{$wordone};
	if ($n <= $max && $n >= $min){
	    $inc = 1.0/$n;
	    $lineno + $inc;
	    $bigrams{$tagone . " " . $tagtwo} += $inc;
	    $left{$tagone} += $inc;
	    $right{$tagtwo} += $inc;
	    $count += $inc;
	}
    }
    elsif ($lineone =~ m/^[^\s]+$/){
	print("Unrecognised line $lineone");
    } 
}

print(" $lineno valid lines\n");
$mi = 0;
$hx = 0;
$hy = 0;
print("Starting to calc \n");
$hxy = 0;
foreach $tagone (keys %left) {
 #   print("Tag one $tagone \n");
    $p = $left{$tagone}/$count;
    $hx -= $p * log($p);
    foreach $tagtwo (keys %right){
        $q = $right{$tagtwo}/$count;
        if (exists($bigrams{$tagone . " " . $tagtwo})){
            $pq = $bigrams{$tagone . " " . $tagtwo}/$count;
  #         print("$tagone $tagtwo $pq \n");
            $mi += $pq * log($pq / ($p * $q));
            $hxy -= $pq * log($pq);
            
        }
    }
}
foreach $tagtwo (keys %right){
    $q = $right{$tagtwo}/$count;
    $hy -= $q * log($q);
}
print("X is gold standard, Y is clusters\n");
print("I(X;Y) =  $mi\n");
print("H(X,Y) = $hxy\n");
print("H(X) = $hx\n");
print("H(Y) = $hy\n");
$hxgiveny = $hxy - $hy;
print("H(X|Y) = $hxgiveny\n");
$hygivenx = $hxy - $hx;
print("H(Y|X) = $hygivenx\n");
print("$max $hxgiveny\n");
