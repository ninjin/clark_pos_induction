#!/usr/bin/perl 
# $Id: tag_text_with_clusters.pl,v 1.2 2002/08/14 13:46:18 clark Exp $
# This file takes two arguments 
# tag_text_with_clusters text clusters 


$corpus = $ARGV[0];
$clusters = $ARGV[1];

$numberLines = 0;
$numberGoodLines = 0;

open(CLUSTERS,"<$clusters") or die "Couldnt open file";

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
	if ($count{$wordone} <= $max && $count{$wordone} >= $min){
	    $lineno++;
	    $bigrams{$tagone . " " . $tagtwo}++;
	    $left{$tagone}++;
	    $right{$tagtwo}++;
	    $count++;
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
