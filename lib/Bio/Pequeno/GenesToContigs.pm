package Bio::Pequeno::GenesToContigs;

# ABSTRACT:  Given a GFF file and a list of IDs extract contigs names

=head1 SYNOPSIS

Given a GFF file and a list of IDs extract contigs names

=cut

use Moose;
use Bio::Tools::GFF;

has 'gff_file'  => ( is => 'ro', isa => 'Str',  required => 1 );
has 'gene_ids'  => ( is => 'ro', isa => 'ArrayRef',  required => 1 );

has 'sequence_ids_to_genes' => ( is => 'ro', isa => 'HashRef', lazy => 1, builder => '_build_sequence_ids_to_genes' );
has 'contigs' => ( is => 'ro', isa => 'HashRef', lazy => 1, builder => '' );

#
#loop over file to link gene_ids to sequenc_ids
#loop again - store all annonotation for a contig, flag which genes are unclassified
#add the sequence to a contig



sub _build_sequence_ids_to_genes
{
	my ($self) = @_;
	my %seq_ids_to_genes;
    my $gffio = Bio::Tools::GFF->new( -file => $self->gff_file, -gff_version => 3 );
    while ( my $feature = $gffio->next_feature() ) {
        my $gene_id = $self->_get_feature_id($feature);
        next unless ($gene_id);
		
		push(@{$seq_ids_to_genes{$feature->seq_id}},$gene_id);
	}
	return \%seq_ids_to_genes;
}


sub _get_feature_id {
    my ( $self, $feature ) = @_;
    my ( $gene_id, @junk );
    if ( $feature->has_tag('ID') ) {
        ( $gene_id, @junk ) = $feature->get_tag_values('ID');
    }
    elsif ( $feature->has_tag('locus_tag') ) {
        ( $gene_id, @junk ) = $feature->get_tag_values('locus_tag');
    }
    else {
        return undef;
    }
    $gene_id =~ s!["']!!g;
    return undef if ( $gene_id eq "" );
    return $gene_id;
}


__PACKAGE__->meta->make_immutable;
no Moose;
1;
