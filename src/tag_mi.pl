#!/usr/bin/perl -w
# $Id: tag_mi.pl,v 1.1 2002/05/30 09:09:15 clark Exp $
# calculates the mutual information between two tag sequences



($one,$two) = @ARGV;
print("$one and $two \n");

open(ONE,"<$one") or die $!;
open(TWO,"<$two") or die $!;
$count = 0;
$lineno = 0;
while ($lineone = <ONE>){
    $lineno++;
    $linetwo = <TWO>;
    if ($lineone =~ m/^\s*$/){
	if (!($linetwo =~ m/^\s*$/)){
	    die("Line mismatch $linetwo at line number $lineno\n");
	}
    }
    else {
	if ($lineone =~ m/^(.*) (.*)$/){
	    $wordone = $1;
	    $tagone = $2;
	    if ($linetwo =~ m/^(.*) (.*)$/){
		$wordtwo = $1;
		$tagtwo = $2;
		if ($wordone eq $wordtwo){
		    $bigrams{$tagone . " " . $tagtwo}++;
		    $left{$tagone}++;
		    $right{$tagtwo}++;
		    $count++;
		}
		else {
		    die("Word mismatch $wordone eq $wordtwo at line number $lineno\n");
		}
	    }
	}
    }
}
print("Read $count lines \n");
$mi = 0;
$hx = 0;
$hy = 0;
$hxy = 0;
foreach $tagone (keys %left) {
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
    $q = $right{$tagtwo}/$count;
    $hy -= $q * log($q);
}
print("I(X;Y) =  $mi\n");
print("H(X,Y) = $hxy\n");
print("H(X) = $hx\n");
print("H(Y) = $hy\n");
$hxgiveny = $hxy - $hy;
print("H(X|Y) = $hxgiveny\n");
$hygivenx = $hxy - $hx;
print("H(Y|X) = $hygivenx\n");

