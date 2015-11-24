#!/usr/bin/env perl

use strict;
use warnings;

die "USAGE: $0 pAss directory outputFile\n" if $#ARGV != 1;

open my $out, ">", $ARGV[1]
    || die "cannot open $ARGV[1] $!";

print $out join("\t", qw/ko contig/), "\n";

for (grep { $_ =~ /fna$/ } <$ARGV[0]/*>)
{
    my $i = 0;
    my ($theKO) = $_ =~ m/(K\d{5})/;

    open my $in, "<", $_
        || die "cannot open $_ $!";

    while(<$in>)
    {
        print $out "$theKO\t$1\n" if m/^\>(contig\d+)/
    }
}

__DATA__
>contig00295 ## spanning:89 msaStart:240 msaEND:605 max10BPwindow:143.3
