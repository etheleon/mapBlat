#!/usr/bin/env perl

use Modern::Perl '2015';
use Getopt::Lucid qw/:all/;
use experimental qw/signatures postderef/;
use Pod::Usage;
use Bio::SeqIO;
use autodie;
use Parallel::ForkManager;


my $opt = Getopt::Lucid->getopt([
    Param("contigs|c"),
    Param("output|o"),
    Param("fastq|q"),
    Param("mapper|m")->default("blat"),
    Param("threads|t")->default(1),
    Param("ko|k"),
    Param("help|h")
]);
$opt->validate({'requires' => ['contigs', 'fastqFile', 'outputDIR']});
pod2usage(-verbose=>2) if $opt->get_help;

my $threads = $opt->get_threads;
my $output  = $opt->get_outputDIR;
my $fq      = $opt->get_fastqFile;
my $contigs = $opt->get_contigs;
my $mapper  = $opt->get_mapper;
my $kos     = $opt->get_ko;

open my $kofile, "$kos";
my @ko;
push(@ko, $_) while <$kofile>;

my $pm = Parallel::ForkManager->new($threads);

if($mapper eq "blat")
{
    foreach my $ko (@ko)
    {
        my $pid = $pm->start and next;
        die "$ko does not have contigs\n" unless -e "$contigs/$ko/454AllContigs.fna";
        convert("$fq/$ko.1.fq", "$fq/$ko.2.fq", "$output/$ko.fasta");
        runBlat("$contigs/$ko/454AllContigs.fna", "$output/$ko.fasta", $output, $ko);
        system "rm $output/$ko.fasta";
        $pm->finish;
    }
    $pm->wait_all_children;
}


sub runBWAmem{
say "# Indexing";
    `bwa index $contigs`;
say "# Mapping";
    #Separately because the reads dun come together
    #Trims quality by 20
    my @a = `bwa mem $contigs $fq.1.fq | samtools view -S -q 20 `;
    my @b = `bwa mem $contigs $fq.2.fq`;
    parseBWA(\@a, \@b);
}

sub parseBWA($alnOne, $alnTwo){
my @a = $alnOne->@*;
    foreach(@a){
        chomp;
        unless(m/^\@/)
        {
            #bwa assumes phred+33
            my($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIGAR,$RNEXT,$PNEXT,$TLEN,$SEQ,$QUAL)=(split(/\t/));
            say $QNAME;
        }
    }
}

sub convert($fastq1, $fastq2, $fasta)
{
    say "# mapper::BLAT\n# Changing fq to fa format";

    `rm -rf $fasta` if -e $fasta;
    foreach my $fastQ (($fastq1, $fastq2))
    {
        if(-e $fastQ)
        {
            my $in  = Bio::SeqIO->new(-file => $fastQ,    -format=>'fastq');
            my $out = Bio::SeqIO->new(-file => ">>$fasta",-format=>'fasta');
            print $out $_ while <$in>;
        }
    }
}

sub runBlat($contigFile, $fasta, $output, $ko)
{
    say "# BLAT: Running BLAT at 100% identity, 100% alignment";
    `blat -fastMap -minIdentity=100 -out=blast8 $contigFile $fasta $output/blatouput_$ko`;
#    open my $outputFile , ">" , "$output/trimmed_blatoutput_$ko";
#    say $outputFile join "\t", qw/query subject qStart qEnd sStart send/;
    #open my $input  , "<" , "$output/blatouput_$ko";

    #while(<$input>)
    #{
        #my ($q, $s, undef, $aln, undef, undef, $qstart, $qend, $sstart, $send, undef, undef) = (split /\t/);
        #say $outputFile join("\t", $q, $s, $qstart, $qend, $sstart, $send) if $aln == 101;
    #}
}

=pod

=head1 NAME

 metaMAPPER - Wrapper for mapping cDNA and gDNA metagenomic reads onto gene centric assemblies from the pAss pipeline. Outputs to tabular blast format

=head1 SYPNOPSIS

 metaMAPPER -contigs path/to/contig.fa -q fastq1 -q fastq2

=head1 BLAT
 BLAT AND Calculate RPKM

   -fastMap    Run for fast DNA/DNA remapping - not allowing introns,
                  requiring high %ID
   -minMatch=N sets the number of tile matches.  Usually set from 2 to 4
               Default is 2 for nucleotide, 1 for protein.
   -minScore=N sets minimum score.  This is the matches minus the
               mismatches minus some sort of gap penalty.  Default is 30
   -minIdentity=N Sets minimum sequence identity (in percent).  Default is
               90 for nucleotide searches, 25 for protein or translated
               protein searches.

 BWA parameters

=over 4

=item --contigs, -c

    Path to contig fasta files. Path/KXXXXX/454AllContigs.fna

=item --output, -o

    Path to the folder in which to output files

=item --fastq, -q

    Path to KO fastq files. FASTQ files must follow the following structure /path/to/fastq/KXXXXX.1.fq and /path/to/fastq/KXXXXX.1.fq

=item --mapper, -m

    Mapper to use

=back

=cut

