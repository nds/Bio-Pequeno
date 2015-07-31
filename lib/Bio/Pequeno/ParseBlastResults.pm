package Bio::Pequeno::ParseBlastResults;

# ABSTRACT:  filter blast results

=head1 SYNOPSIS

filter blast results

=cut

use Moose;
use Bio::SeqIO;

has 'blast_results'            => ( is => 'ro', isa => 'ArrayRef',       required => 1 );
has 'fasta_file'               => ( is => 'ro', isa => 'Str',            required => 1 );
has '_filtered_results'        => ( is => 'ro', isa => 'Maybe[HashRef]', lazy     => 1, builder => '_build__filtered_results' );
has 'minimum_alignment_length' => ( is => 'ro', isa => 'Int',            default  => 400 );
has 'max_percentage_coverage'  => ( is => 'ro', isa => 'Num',            default  => 0.5 );

sub sequence_calculate_coverage {
    my ($self) = @_;

    my %sequence_coverage;
    my $seq_io = Bio::SeqIO->new( -file => $self->fasta_file, -format => 'Fasta' );
    while ( my $input_seq = $seq_io->next_seq() ) {
        my $sequence_name     = $input_seq->display_id();
        my $sequence_length   = $input_seq->length();
        my @sequence_coverage = 0 x $sequence_length;

        for my $blast_result ( @{ $self->_filtered_results->{$sequence_name} } ) {

            my $start_coord = $blast_result->[6];
            my $end_coord   = $blast_result->[7];
            if ( $start_coord > $end_coord ) {
                my $tmp = $start_coord;
                $start_coord = $end_coord;
                $end_coord   = $tmp;
            }
            for ( my $i = $start_coord ; $i < $end_coord ; $i++ ) {
                $sequence_coverage[ $i - 1 ] = 1;
            }
        }

        # count the number of 1's.
        my $sum = 0;
        for ( my $i = 0 ; $i < @sequence_coverage ; $i++ ) {
            $sum += $sequence_coverage[$i] if ( defined( $sequence_coverage[$i] ) );
        }

        my $coverage = $sum / $sequence_length;
        next if ( $coverage > $self->max_percentage_coverage );
        $sequence_coverage{$sequence_name} = $coverage;
    }
    return \%sequence_coverage;
}

sub _build__filtered_results {
    my ($self) = @_;
    my %filtered_results;

    for my $blast_result ( @{ $self->blast_results } ) {
        next if ( !defined($blast_result) || $blast_result eq "" );
        my @blast_result_details = split( /\t/, $blast_result );

        next if ( $blast_result_details[3] < $self->minimum_alignment_length );
        push( @{ $filtered_results{ $blast_result_details[0] } }, \@blast_result_details );
    }

    return \%filtered_results;
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
