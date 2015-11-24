#!/usr/bin/env perl

use strict;
use warnings;

die "$0 <contigDir> <targetFile>\n" unless $#ARGV == 1;

my @kos = <"$ARGV[0]/*/454AllContigs.fna">;
print "## Registered all Files!\n";
print "##\tTime to read\n";

open my $out, ">", $ARGV[1] || die "$! cannot open output file: $ARGV[1]";

print $out join("\t", qw/ko contig length numreads/), "\n";
for my $contigFile (@kos)
{
    my ($ko) = $contigFile =~ m/(K\d{5})/;
    open(my $in, "<", $contigFile) || die "$! cannot open contig fasta";
    while(<$in>)
    {
        if(/^>/)
        {
            m/>(\S+)\s+length=(\d+)\s+numreads=(\d+)$/;
            print $out join("\t", $ko, $1, $2, $3), "\n";
        }
    }
}
close $out;

__END__
>contig00002  length=876   numreads=530
