package Bio::Pequeno::Blastn;

# ABSTRACT:  take in a fasta file of nucleotides and run blastn

=head1 SYNOPSIS

take in a fasta file of nucleotides and run blastn

=cut

use Moose;
use Bio::SeqIO;
use File::Path qw(make_path);
use File::Basename;
use Bio::Tools::GFF;

has 'fasta_file'      => ( is => 'ro', isa => 'Str',      required => 1 );
has 'blast_database'  => ( is => 'ro', isa => 'Str',      default  => '/data/blastdb/Supported/NT/nt' );
has 'blastn_exec'     => ( is => 'ro', isa => 'Str',      default  => 'blastn' );
has 'evalue'          => ( is => 'ro', isa => 'Num',      default  => 0.00001 );
has 'cpus'            => ( is => 'ro', isa => 'Int',      default  => 1 );
has 'perc_identity'   => ( is => 'ro', isa => 'Num',      default  => 70 );
has 'max_target_seqs' => ( is => 'ro', isa => 'Int',      default  => 2000 );
has 'outfmt'          => ( is => 'ro', isa => 'Int',      default  => 6 );
has 'blast_results'   => ( is => 'ro', isa => 'ArrayRef', lazy     => 1, builder => '_build_blast_results' );

sub _blast_command {
    my ($self) = @_;
    return join(
        ' ',
        (
            $self->blastn_exec,   '-db',              $self->blast_database, '-evalue',
            $self->evalue,        '-num_threads',     $self->cpus,           '-outfmt',
            $self->outfmt,        '-query',           $self->fasta_file,     '-perc_identity',
            $self->perc_identity, '-max_target_seqs', $self->max_target_seqs
        )
    );
}

sub _build_blast_results {
    my ($self) = @_;
    my @results;
    open( my $results_fh, '-|', $self->_blast_command );
    while (<$results_fh>) {
		chomp();
        push( @results, $_ );
    }
    return \@results;
}

__PACKAGE__->meta->make_immutable;
no Moose;
1;
