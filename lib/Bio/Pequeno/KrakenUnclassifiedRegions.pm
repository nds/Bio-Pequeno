package Bio::Pequeno::KrakenUnclassifiedRegions;

# ABSTRACT:  Find unclassified sequences with kraken

=head1 SYNOPSIS

Given a FASTA of nucleotide sequences of coding regions, find unclassified sequences with kraken

=cut

use Moose;

has 'kraken_db'                 => ( is => 'ro', isa => 'Str', default  => "/lustre/scratch108/pathogen/pathpipe/kraken/pi_qc_2015521/" );
has 'fasta_file'                => ( is => 'ro', isa => 'Str', required => 1 );
has 'cpus'                      => ( is => 'ro', isa => 'Int', default  => 1 );
has 'minimum_gene_id_threshold' => ( is => 'ro', isa => 'Int', default  => 3 );
has 'kraken_exec'               => ( is => 'ro', isa => 'Str', default  => 'kraken' );

has '_report_filename'      => ( is => 'ro', isa => 'Str',             lazy => 1, builder => '_build__report_filename' );
has 'unclassified_gene_ids' => ( is => 'ro', isa => 'Maybe[ArrayRef]', lazy => 1, builder => '_build_unclassified_gene_ids' );

sub _build__report_filename {
    my ($self) = @_;
    return $self->fasta_file . ".report";
}

sub _kraken_cmd {
    my ($self) = @_;
    return
        $self->kraken_exec
      . " --fasta-input --threads "
      . $self->cpus
      . " --db "
      . $self->kraken_db
      . " --quick --preload "
      . $self->fasta_file . " > "
      . $self->_report_filename;
}

sub _get_unclassified_coding_regions {
    my ( $self, $filename ) = @_;
    open( my $fh, $filename );
    my @unclassified;
    while (<$fh>) {
        chomp;
        my $line = $_;

        # get all unclassified
        next unless ( $line =~ /^U/ );
        my @classification_details = split( /\t/, $line );
        push( @unclassified, $classification_details[1] );

    }
    my @sorted_unclassified = sort @unclassified;
    return \@sorted_unclassified;
}

sub _build_unclassified_gene_ids {
    my ($self) = @_;
    system( $self->_kraken_cmd );
    my $unclassified = $self->_get_unclassified_coding_regions( $self->_report_filename );
    unlink( $self->_report_filename );
    return undef unless ( @{$unclassified} >= $self->minimum_gene_id_threshold );
    return $unclassified;
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
