package Bio::Pequeno::ExtractSequencesToFasta;

# ABSTRACT:  given a list of sequence ids, pull them out of a fasta file and output them to a new fasta file

=head1 SYNOPSIS

given a list of sequence ids, pull them out of a fasta file and output them to a new fasta file

=cut

use Moose;
use Bio::Tools::GFF;
use Bio::SeqIO;
use File::Basename;

has 'fasta_file'      => ( is => 'ro', isa => 'Str',            required => 1 );
has 'sequence_ids'    => ( is => 'ro', isa => 'Maybe[HashRef]', required => 1 );
has 'output_filename' => ( is => 'ro', isa => 'Str',            lazy     => 1, builder => '_build_output_filename' );

sub _build_output_filename {
    my ($self) = @_;
    my ( $filename, $directories, $suffix ) = fileparse( $self->fasta_file, qr/\.[^.]*/ );
    return $filename . '.nohits.fa';
}

sub extract_sequence_ids_to_file {
    my ($self) = @_;
    return undef unless ( defined( $self->sequence_ids ) );

    my $out_seq_io = Bio::SeqIO->new( -file => ">" . $self->output_filename, -format => 'Fasta' );
    my $seq_io     = Bio::SeqIO->new( -file => $self->fasta_file,            -format => 'Fasta' );
    while ( my $input_seq = $seq_io->next_seq() ) {
        if ( defined( $self->sequence_ids->{ $input_seq->display_id() } ) ) {
            $out_seq_io->write_seq($input_seq);
        }
    }
    return 1;
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
