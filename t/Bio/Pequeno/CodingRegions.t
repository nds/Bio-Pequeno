#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::Files;

BEGIN { unshift( @INC, './lib' ) }
$ENV{PATH} .= ":./bin";

BEGIN {
    use Test::Most;
    use_ok('Bio::Pequeno::CodingRegions');
}

my $obj;
ok($obj = Bio::Pequeno::CodingRegions->new(gff_file => 't/data/query_1.gff'), 'initialise obj');
is($obj->fasta_file,'query_1.fa', 'fasta file created');
ok(-e $obj->fasta_file, 'fasta file created');

compare_ok( 'query_1.fa', 't/data/expected_query_1.fa', 'extracted nuc sequences as expected' );

done_testing();

unlink('query_1.fa');