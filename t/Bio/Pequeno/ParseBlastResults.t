#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::Files;

BEGIN { unshift( @INC, './lib' ) }
$ENV{PATH} .= ":./bin";

BEGIN {
    use Test::Most;
    use_ok('Bio::Pequeno::ParseBlastResults');
    use_ok('Bio::Pequeno::Blastn');
}

calculate_coverage(
    't/data/query_1.fa',
    't/data/blast_results_query_1',
    {
        'abc_00012' => '0',
        'abc_00004' => '0',
        'abc_00011' => '0',
        'abc_00006' => '0',
        '1_6'       => '0',
        'abc_00014' => '0',
        'abc_00002' => '0',
        'abc_00003' => '0',
        'abc_00010' => '0'
    }
);

calculate_coverage(
    't/data/query_1.fa',
    't/data/blast_results_query_1_small_hits',
    {
        'abc_00012' => '0',
        'abc_00004' => '0',
        'abc_00011' => '0',
        'abc_00006' => '0',
        '1_6'       => '0',
        'abc_01705' => '0',
        'abc_00014' => '0',
        'abc_00002' => '0',
        'abc_00003' => '0',
        'abc_00010' => '0'
    }
);

done_testing();

sub calculate_coverage {
    my ( $fasta_file, $blast_results, $expected_coverage ) = @_;

    ok( my $blast_obj = Bio::Pequeno::Blastn->new( fasta_file => $fasta_file, _blast_command => "cat $blast_results" ) );
    ok( my $obj = Bio::Pequeno::ParseBlastResults->new( fasta_file => $fasta_file, blast_results => $blast_obj->blast_results ),
        'initialise obj' );
    is_deeply( $obj->sequence_calculate_coverage, $expected_coverage, 'sequence coverage' );
}

