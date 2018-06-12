#!/usr/bin/perl

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding

use Bio::Perl;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

my $taxDB = Bio::LITE::Taxonomy::NCBI-> new (
                                             db=>"NCBI",
                                             names=>
                                             "/SAN/db/taxonomy/names.dmp",
                                             nodes=>
                                             "/SAN/db/taxonomy/nodes.dmp",
                                            );

my $cmd = "blast_formatter -archive ".$ARGV[0].
  " -outfmt \'6 qseqid sacc bitscore staxids pident length sstart send\' -max_target_seqs 5"; 

print STDERR "running $cmd \n";

my @blast = `$cmd`; 

print STDERR "obtaining species, genus, family, phylum, kingdom and superkingdom\n\n";

print  "query,subject,bitscore,taxid,pident,length,substart,subend,species,genus,family,order,class,phylum,kingdom,superkingdom\n";

foreach (@blast) {
  next if /#/;
  chomp($_);
  my @blt = split (" ", $_);
  my $species = $taxDB->get_term_at_level($blt[3],"species");
  my $genus = $taxDB->get_term_at_level($blt[3],"genus");
  my $family = $taxDB->get_term_at_level($blt[3],"family");
  my $order = $taxDB->get_term_at_level($blt[3],"order");
  my $class = $taxDB->get_term_at_level($blt[3],"class");
  my $phylum = $taxDB->get_term_at_level($blt[3],"phylum");
  my $kingdom = $taxDB->get_term_at_level($blt[3],"kingdom");
  my $superkingdom = $taxDB->get_term_at_level($blt[3],"superkingdom");

  print "$blt[0],$blt[1],$blt[2],$blt[3],$blt[4],$blt[5],$blt[6],$blt[7],$species,$genus,$family,$order,$class,$phylum,$kingdom,$superkingdom\n";
}

