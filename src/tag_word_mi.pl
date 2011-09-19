#!/usr/bin/perl -w
# $Id: tag_word_mi.pl,v 1.2 2002/08/14 13:46:18 clark Exp $
# calculates the mutual information and other figures
# between the first symbol and the second symbol



($one) = @ARGV;

open(ONE,"<$one") or die $!;

$count = 0;
$lineno = 0;
while ($lineone = <ONE>){
    $lineno++;
    if ($lineone =~ m/^(.*) (.*)$/){
	$tagone = $1;
	$tagtwo = $2;
	$bigrams{$tagone . " " . $tagtwo}++;
	$left{$tagone}++;
	$right{$tagtwo}++;
	$count++;
    }
}
print("Read $count lines out of $lineno total lines\n");
$mi = 0;
$hx = 0;
$hy = 0;
$hxy = 0;
foreach $tagone (keys %left) {
    $x++;
    $p = $left{$tagone}/$count;
    $hx -= $p * log($p);
    foreach $tagtwo (keys %right){
	$q = $right{$tagtwo}/$count;
	if (exists($bigrams{$tagone . " " . $tagtwo})){
	    $pq = $bigrams{$tagone . " " . $tagtwo}/$count;
#	    print("$tagone $tagtwo $pq \n");
	    $mi += $pq * log($pq / ($p * $q));
	    $hxy -= $pq * log($pq);
	    
	}
    }
}
foreach $tagtwo (keys %right){
    $y++;
    $q = $right{$tagtwo}/$count;
    $hy -= $q * log($q);
}
print("Distinct types X   $x Y $y\n");
print("I(X;Y) =  $mi\n");
print("H(X,Y) = $hxy\n");
print("H(X) = $hx\n");
print("H(Y) = $hy\n");
$hxgiveny = $hxy - $hy;
print("H(X|Y) = $hxgiveny\n");
$hygivenx = $hxy - $hx;
print("H(Y|X) = $hygivenx\n");

