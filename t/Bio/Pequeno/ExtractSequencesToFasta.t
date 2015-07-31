#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::Files;

BEGIN { unshift( @INC, './lib' ) }
$ENV{PATH} .= ":./bin";

BEGIN {
    use Test::Most;
    use_ok('Bio::Pequeno::ExtractSequencesToFasta');
}

my $obj;
ok(
    $obj = Bio::Pequeno::ExtractSequencesToFasta->new(
        fasta_file   => 't/data/query_1.fa',
        sequence_ids => { 'abc_00004' => 0, 'abc_00006' => 0.01, }
    ),
    'initialise obj'
);

ok( $obj->extract_sequence_ids_to_file, 'extract nuc sequences from blocks' );
is( $obj->output_filename, 'query_1.nohits.fa', 'nohits fasta filename' );
ok( -e $obj->output_filename, 'nohits fasta file created' );

compare_ok( 'query_1.nohits.fa', 't/data/expected_query_1.nohits.fa', 'extracted nuc sequences as expected' );
unlink('query_1.nohits.fa');

done_testing();

