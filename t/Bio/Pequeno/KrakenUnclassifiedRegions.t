#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::Files;

BEGIN { unshift( @INC, './lib' ) }
$ENV{PATH} .= ":./bin";

BEGIN {
    use Test::Most;
    use_ok('Bio::Pequeno::KrakenUnclassifiedRegions');
}

my $obj;
ok(
    $obj = Bio::Pequeno::KrakenUnclassifiedRegions->new( fasta_file => 't/data/query_1.fa', kraken_db => 't/data/kraken_test/',minimum_gene_id_threshold => 1),
    'initialise obj'
);

is(
    $obj->_kraken_cmd,
    'kraken --fasta-input --threads 1 --db t/data/kraken_test/ --quick --preload t/data/query_1.fa > t/data/query_1.fa.report',
    'kraken command'
);
is_deeply( $obj->unclassified_gene_ids, [ '1_6', 'abc_00012' ], 'unclassified gene ids' );

done_testing();
