#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::Files;

BEGIN { unshift( @INC, './lib' ) }
$ENV{PATH} .= ":./bin";

BEGIN {
    use Test::Most;
    use_ok('Bio::Pequeno::GenesToContigs');
}

my $obj;
ok( $obj = Bio::Pequeno::GenesToContigs->new( gff_file => 't/data/query_1.gff',gene_ids => 
['abc_00002', 'abc_00004'], min_sequence_length => 200, min_genes_on_contig => 2
 ), 'initialise obj' );
ok( $obj->extract_nuc_sequences_from_blocks, 'extract nuc sequences from blocks' );
is($obj->output_filename,'query_1.novel.fa', 'novel fasta filename' );
ok( -e $obj->output_filename, 'novel fasta file created' );

compare_ok( 'query_1.novel.fa', 't/data/expected_query_1.novel.fa', 'extracted nuc sequences as expected' );

done_testing();

unlink('query_1.novel.fa');
