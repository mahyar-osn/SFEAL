#!/usr/bin/perl

use strict;

my $nfile = $ARGV[0];
my $ofile = $ARGV[1];

open IPNODE, "<$nfile" or die "\033[31mError: Can't open $nfile\033[0m ";
my ($node, $i, @xyz, $nv, @deriv);

my $line = <IPNODE>;
while ($line) {
    if ($line =~ /Node number \[.*\]:\s*(.*)/) {
        $node = $1;
        $nv = 1;
    }
    if ($line =~ /For version number (\d):/) {$nv = $1;}
    if ($line =~ /The Xj\((\d)\) coordinate is \[.*\]:\s*(.*)/){
        $i = $1;
        $xyz[$node][$i-1] = $2;
    }
    if ($line =~ /The derivative wrt direction 1 is \[.*\]:\s*(.*)/) {$deriv[$node][$nv][$i][0] = $1;}
    if ($line =~ /The derivative wrt direction 2 is \[.*\]:\s*(.*)/) {$deriv[$node][$nv][$i][1] = $1;}
    if ($line =~ /The derivative wrt directions 1 \& 2 is \[.*\]:\s*(.*)/) {$deriv[$node][$nv][$i][2] = $1;}
    $line = <IPNODE>;
}
close IPNODE;

my (@versions,@Nodes,@base);
@Nodes = (1..50);

# Nodes with 1 version
for my $i (1..50){
    $versions[$i] = 1;
}

#Nodes with 2 versions
for my $i (2,4..6,8,10..13,15..18,23..24,30,32,34,36..37,39..40,44,49){
    $versions[$i] = 2;
}

# Nodes with 3 versions
for my $i (3,14){
	$versions[$i] = 3;
}

# Nodes with 4 versions
for my $i (9,35){
	$versions[$i] = 4 
}

# Nodes with 5 versions
for my $i (1,7,20,26,29,33,48,50){
	$versions[$i] = 5; 
}

# Nodes with 6 versions
for my $i (31,38){
	$versions[$i] = 6; 
}

open OPNODE, ">$ofile" or die "\033[31mError: Can't open $ofile\033[0m ";

&PrintModifiedNodes;

# --------------------------------------------------------------------------------
sub PrintModifiedNodes {
    foreach my $i (@Nodes){
        if($versions[$i] == 1){
            $nv = 1;
            print OPNODE "$i $xyz[$i][0] $deriv[$i][$nv][1][0] $deriv[$i][$nv][1][1] 0\n$i $xyz[$i][1] $deriv[$i][$nv][2][0] $deriv[$i][$nv][2][1] 0 \n$i $xyz[$i][2] $deriv[$i][$nv][3][0] $deriv[$i][$nv][3][1] 0 \n";
        }else{
            for my $nv (1..$versions[$i]) {
                print OPNODE "$i.$nv $xyz[$i][0] $deriv[$i][$nv][1][0] $deriv[$i][$nv][1][1] 0\n$i.$nv $xyz[$i][1] $deriv[$i][$nv][2][0] $deriv[$i][$nv][2][1] 0 \n$i.$nv $xyz[$i][2] $deriv[$i][$nv][3][0] $deriv[$i][$nv][3][1] 0 \n";
            }
        }
    }
}
my $arraysize = @Nodes;
#print "Total Nodes = $arraysize \n";
