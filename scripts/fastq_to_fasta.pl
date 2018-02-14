#!/usr/bin/perl
# Kirk Amundson
# QtoA.pl
# 3 October 2016
# Accept fastq formatted file at command line, write to converted fasta file

# USAGE: program.pl input.fastq

# interpret file provided at command line
my $fastq = shift @ARGV; "Usage: program.pl input.fastq\n";
my $fasta = (split /\./, $fastq)[0] . ".fasta";

# open files and process, reading from fastq and writing to fasta
open FASTQ, "<", "$fastq" or die "Usage: program.pl input.fastq: $!";
open FASTA, ">", "$fasta" or die "problem opening $fasta\n";
while (my $line = <FASTQ>) {
    $linecount++;
    if ($linecount%4 == 1) {
        $line =~ s/^@/>/g;
        print FASTA $line;
    }
    elsif ($linecount%4 == 2) {
        print FASTA $line;
    }
}