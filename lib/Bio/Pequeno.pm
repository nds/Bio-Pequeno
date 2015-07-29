undef $VERSION;

package Bio::Pequeno;

# ABSTRACT: Find novel sequences

=head1 SYNOPSIS

Find novel sequences

=cut

use Moose;
use Getopt::Long qw(GetOptionsFromArray);
use Bio::Tools::GFF;
use Bio::Pequeno::CodingRegions;
use Bio::Pequeno::KrakenUnclassifiedRegions;
use Bio::Pequeno::GenesToContigs;

has 'args'        => ( is => 'ro', isa => 'ArrayRef', required => 1 );
has 'script_name' => ( is => 'ro', isa => 'Str',      required => 1 );
has 'help'        => ( is => 'rw', isa => 'Bool',     default  => 0 );

has 'input_files'               => ( is => 'rw', isa => 'ArrayRef' );
has 'kraken_db'                 => ( is => 'rw', isa => 'Str', default => "/lustre/scratch108/pathogen/pathpipe/kraken/pi_qc_2015521/" );
has 'cpus'                      => ( is => 'rw', isa => 'Int', default => 1 );
has 'minimum_gene_id_threshold' => ( is => 'rw', isa => 'Int', default => 3 );

sub BUILD {
    my ($self) = @_;

    my ( $input_files, $kraken_db, $help, $cpus, $minimum_gene_id_threshold, );

    GetOptionsFromArray(
        $self->args,
        'd|kraken_db=s'  => \$kraken_db,
        'p|processors=i' => \$cpus,
        'u|min_genes=i'  => \$minimum_gene_id_threshold,
        'h|help'         => \$help,
    );

    $self->kraken_db($kraken_db)                                 if ( defined($kraken_db) );
    $self->cpus($cpus)                                           if ( defined($cpus) );
    $self->minimum_gene_id_threshold($minimum_gene_id_threshold) if ( defined($minimum_gene_id_threshold) );
    $self->input_files( $self->args );

}

sub run {
    my ($self) = @_;

    ( !$self->help ) or die $self->usage_text;
    for my $file ( @{ $self->input_files } ) {
        my $coding = Bio::Pequeno::CodingRegions->new( gff_file => $file );
        my $unclassified_genes_obj = Bio::Pequeno::KrakenUnclassifiedRegions->new(
            fasta_file                => $coding->fasta_file,
            cpus                      => $self->cpus,
            kraken_db                 => $self->kraken_db,
            minimum_gene_id_threshold => $self->minimum_gene_id_threshold
        );

        my $genes_to_contigs_obj =
          Bio::Pequeno::GenesToContigs->new( gff_file => $file, gene_ids => $unclassified_genes_obj->unclassified_gene_ids );
		$genes_to_contigs_obj->print_gene_products_on_contigs;
		$genes_to_contigs_obj->extract_nuc_sequences_from_blocks;

        #extract annotation thats not hypo
        #blast coding regions against nt
    }

}

sub usage_text {
    my ($self) = @_;

    return <<USAGE;
    Usage: pequeno [options]
    Find novel sequences
	
 	 -c : number of cpus
	 -d : path to kraken database
	 -u : minimum number of unclassified genes in a genome to consider

    # This help message
    pequeno -h

USAGE
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
