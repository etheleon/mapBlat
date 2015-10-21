#!/usr/bin/env perl

#Mapping to neo4j based on name is going to be very SLOW and inaccurate
#This will parse Daniel's BLAST2LCA mapper output to reflect the NCBI taxonIDs


die "USAGE $0 <ncbi.map> <blast2lca>\n" unless $#ARGV == 2;

my %map;

# kingdom phylum class order family genus species

open my $map, "<", $ARGV[0] || die "$! ncbi mapping file not found";
while(<$map>)
{
    my ($taxid, $name) = (split(/\t/))[0,1];
    $map{$name} = $taxid;
    #say join "\t", $taxid, $name
}


open my $assignment , "<" , $ARGV[1] || die "$! blast2lca output not found";
open my $output     , ">" , $ARGV[2] || die "$! output file not specified";

print $output join(",", qw/ko contig rank taxid score/), "\n";
while(<$assignment>)
{
    chomp;
    my @entry = (split/;/);
    my $id = shift @entry;
    shift @entry;
    my ($contig, $ko) = (split /\|/, $id);
    while(@entry)
    {
        my $id                  = shift @entry;
        my $score               = shift @entry;
        $score =~ s/\s//g;

        my ($classifier, $name) = (split(/__/, $id));

        print $output join(',', $ko, $contig,  $classifier, $map{$name}, $score), "\n";
    }
}

__DATA__
contig00283|K00001; ;k__Bacteria; 100;p__Proteobacteria; 100;c__Betaproteobacteria; 100;o__Burkholderiales; 96;f__Burkholderiaceae; 76;g__Cupriavidus; 60;s__Cupriavidus taiwanensis; 12;
