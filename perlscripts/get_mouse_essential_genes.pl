#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use String::Util qw(trim);

use Cwd 'abs_path';
use File::Basename;
use Bio::OntologyIO;


my $ver = '1.0';
my $last_modified = 'June 17, 2016';

## --------------------------------------------
## -- version history --
## --------------------------------------------

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s", "o:s", "g:s", "obo" );

if (!$opts{o} ){
    print "----------------------------------------------------------------------
    \t\tversion : $ver by Weihua Chen; last modified : $last_modified
----------------------------------------------------------------------
    USAGE: perl $0
        -o output gene to essentiality file
            1. Ensembl geneID
            2. essentiality (E/NE)
      [optinal input files]
        -i input mouse gene to phenotype file, such as ftp://ftp.informatics.jax.org/pub/reports/MGI_GenePheno.rpt
        -g input mgi to ensembl ID file, ftp://ftp.informatics.jax.org/pub/reports/MGI_Gene_Model_Coord.rpt
        -obo obo file
----------------------------------------------------------------------\n";
    exit;
}

## -- get obo terms --
## -- 
my $script_full_path    = dirname(__FILE__);
my $mgi2phenoFile       = defined $opts{i} ? $opts{i} : $script_full_path . "/mouse_mgi_files/MGI_GenePheno.rpt.gz";
my $mgi2EssemblFile     = defined $opts{g} ? $opts{g} : $script_full_path . "/mouse_mgi_files/MGI_Gene_Model_Coord.rpt.gz";
my $mgiOboFile          = defined $opts{obo} ? $opts{obo} : $script_full_path . "/mouse_mgi_files/MPheno_OBO.ontology";
my $hrEssMPs = &oboParser4EssTerms();

##print Dumper $hrEssMPs;
#exit;


## -- load input files --
print STDERR "\tload mgi to pheotype file ... \n";
my %hEss = (); ## essential genes
my %hOthers = (); ## -- 
if( $mgi2phenoFile =~ /\.gz$/ ){
    open IN, "gunzip -c $mgi2phenoFile |" or die "cannot open pipe to $mgi2phenoFile\n";
	     ##                        ^  -- this | is very important!!
} else {
    open IN, $mgi2phenoFile or die "cannot open pipe to $mgi2phenoFile\n";
}
while(<IN>){
    next if(/^#/ or !/\S/);
    chomp;
    my @arr = split(/\t/, $_);
    my $mp = $arr[4];
    my $mgi = $arr[6];
    if( defined $mp and defined $mgi and length($mp) > 0 and length($mgi) > 0 ){
        if( exists $$hrEssMPs{ $mp } ){
            $hEss{ $mgi } ++;
        } else {
            $hOthers{ $mgi } ++;
        }
    }
}
close IN;
print STDERR "\t\t... done, ", scalar keys %hEss ," ess MGIs\n\n";

## -- --
my %hMgi2Locus = ();
print STDERR "\tload mgi to ensembl gene locus file ... \n";
if( $mgi2EssemblFile =~ /\.gz$/ ){
    open GENE, "gunzip -c $mgi2EssemblFile |" or die "cannot open pipe to $mgi2EssemblFile\n";
	     ##                            ^  -- this | is very important!!
} else {
    open GENE, $mgi2EssemblFile or die "cannot open pipe to $mgi2EssemblFile\n";
}
while(<GENE>){
    next if(/^#/ or !/\S/);
    chomp;
    my @arr = split(/\t/, $_);
    my $ensembl = $arr[10];
    my $type = $arr[ 1 ];
    my $mgi = $arr[0];
    if( defined $ensembl and defined $mgi and length($ensembl) > 0 and length($mgi) > 0 ){
        $hMgi2Locus{ $mgi }{$ensembl} ++;
    }
}
close GENE;
print STDERR "\t\t... done\n\n";

## -- post process 
print STDERR "\tpost process and print Ess genes ... \n";
my %hEnsGenes = ();
my %hOtherGenes = ();
foreach my $mgi ( keys %hEss ){
    if( exists $hMgi2Locus{ $mgi } ){
        foreach my $locus ( keys %{ $hMgi2Locus{ $mgi } } ){
            $hEnsGenes{ $locus } ++;
        }
    }
}

foreach my $mgi ( keys %hOthers ){
    if( exists $hMgi2Locus{ $mgi } ){
        foreach my $locus ( keys %{ $hMgi2Locus{ $mgi } } ){
            $hOtherGenes{ $locus } ++ if( !exists $hEnsGenes{ $locus } );
        }
    }
}

## -- print --
open OUT, ">$opts{o}" or die;
foreach my $locus (keys %hEnsGenes){
    print OUT join("\t", $locus, "E"), "\n";
}
foreach my $locus (keys %hOtherGenes){
    print OUT join("\t", $locus, "NE"), "\n";
}
close OUT;
print STDERR "\t\tdone ... \n\n";
print STDERR "\tall jobs done ... \n\n";



##########################
## -- sub
sub oboParser4EssTerms{
    my %hash = ();
    my $parser = Bio::OntologyIO->new(-format => "obo",-file => $mgiOboFile);
    my $ont = $parser->next_ontology;
    
    foreach my $term ($ont->get_all_terms){
        my $mp = $term->identifier();
        my $name = $term->name;
        
        if( $name =~ /lethality/ or
           $name =~ /infertility/ or
           $name =~ /premature death/){
            $hash{ $mp } = $name;
        }
        
        #print $mp, "\t", $name, "\n";
        
    }

    return \%hash;
}