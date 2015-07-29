package Bio::Pequeno::CodingRegions;

# ABSTRACT:  take in a gff and output a fasta with nucleotide sequences of coding regions

=head1 SYNOPSIS

take in a gff and output a fasta with nucleotide sequences of coding regions

=cut

use Moose;
use Bio::SeqIO;
use File::Path qw(make_path);
use File::Basename;
use Bio::Tools::GFF;

has 'gff_file'         => ( is => 'ro', isa => 'Str',                           required => 1 );

has 'fasta_file'   => ( is => 'ro', isa => 'Str',        lazy => 1, builder => '_build_fasta_file' );
has '_input_seqio' => ( is => 'ro', isa => 'Bio::SeqIO', lazy => 1, builder => '_build__input_seqio' );

has '_tags_to_filter'   => ( is => 'ro', isa => 'Str', default => '(CDS|ncRNA|tRNA|tmRNA|rRNA)' );
has 'min_gene_size_in_nucleotides'   => ( is => 'ro', isa => 'Int',  default  => 120 );
has 'output_filename' => ( is => 'ro', isa => 'Str', lazy => 1, builder => '_build_output_filename' );


sub _build_output_filename
{
  my ($self) = @_;
  my ( $filename, $directories, $suffix ) = fileparse($self->gff_file, qr/\.[^.]*/);
  return  $filename.'.fa' ;
}

sub _bed_output_filename {
    my ($self) = @_;
    return join( '.', ( $self->output_filename, 'intermediate.bed' ) );
}

sub _create_bed_file_from_gff {
    my ($self) = @_;

    open( my $bed_fh, '>', $self->_bed_output_filename );
    my $gffio = Bio::Tools::GFF->new( -file => $self->gff_file, -gff_version => 3 );
    while ( my $feature = $gffio->next_feature() ) {

        next unless defined($feature);

        # Only interested in a few tags
        my $tags_regex = $self->_tags_to_filter;
        next if !( $feature->primary_tag =~ /$tags_regex/ );

        # Must have an ID tag
        my $gene_id = $self->_get_feature_id($feature);
        next unless($gene_id);

        #filter out small genes
        next if ( ( $feature->end - $feature->start ) < $self->min_gene_size_in_nucleotides );

        my $strand = ($feature->strand > 0)? '+':'-' ;
        print {$bed_fh} join( "\t", ( $feature->seq_id, $feature->start -1, $feature->end, $gene_id, 1, $strand ) ) . "\n";
    }
    $gffio->close();
}

sub _get_feature_id
{
    my ($self, $feature) = @_;
    my ( $gene_id, @junk ) ;
    if ( $feature->has_tag('ID') )
    {
         ( $gene_id, @junk ) = $feature->get_tag_values('ID');
    }
    elsif($feature->has_tag('locus_tag'))
    {
        ( $gene_id, @junk ) = $feature->get_tag_values('locus_tag');
    }
    else
    {
        return undef;
    }
    $gene_id =~ s!["']!!g;
    return undef if ( $gene_id eq "" );
    return $gene_id ;
}


sub _create_nucleotide_fasta_file_from_gff {
    my ($self) = @_;
    my $cmd =
        'sed -n \'/##FASTA/,//p\' '
      . $self->gff_file
      . ' | grep -v \'##FASTA\' > '
      . $self->_nucleotide_fasta_file_from_gff_filename;
    system($cmd);
}

sub _nucleotide_fasta_file_from_gff_filename {
    my ($self) = @_;
    return join( '.', ( $self->output_filename, 'intermediate.fa' ) );
}

sub _extract_nucleotide_regions {
    my ($self) = @_;

    $self->_create_nucleotide_fasta_file_from_gff;
    $self->_create_bed_file_from_gff;

    my $cmd =
        'bedtools getfasta -s -fi '
      . $self->_nucleotide_fasta_file_from_gff_filename
      . ' -bed '
      . $self->_bed_output_filename . ' -fo '
      . $self->_extracted_nucleotide_fasta_file_from_bed_filename
      . ' -name > /dev/null 2>&1';
    system($cmd);
    unlink( $self->_nucleotide_fasta_file_from_gff_filename );
    unlink( $self->_bed_output_filename );
    unlink( $self->_nucleotide_fasta_file_from_gff_filename . '.fai' );
	system("mv ".$self->_extracted_nucleotide_fasta_file_from_bed_filename."  ".$self->output_filename);
    return $self->output_filename;
}
sub _build_fasta_file {
    my ($self) = @_;
    return $self->_extract_nucleotide_regions;
}


sub _extracted_nucleotide_fasta_file_from_bed_filename {
    my ($self) = @_;
    return join( '.', ( $self->output_filename, 'extracted.fa' ) );
}


__PACKAGE__->meta->make_immutable;
no Moose;
1;
