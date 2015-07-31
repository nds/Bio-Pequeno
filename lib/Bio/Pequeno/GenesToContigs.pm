package Bio::Pequeno::GenesToContigs;

# ABSTRACT:  Given a GFF file and a list of IDs extract contigs names

=head1 SYNOPSIS

Given a GFF file and a list of IDs extract contigs names

=cut

use Moose;
use Bio::Tools::GFF;
use Bio::SeqIO;
use File::Basename;

has 'gff_file'         => ( is => 'ro', isa => 'Str',      required => 1 );
has 'gene_ids'         => ( is => 'ro', isa => 'ArrayRef', required => 1 );
has 'genes_to_feature' => ( is => 'ro', isa => 'HashRef',  default  => sub { {} } );
has 'genes_to_product' => ( is => 'ro', isa => 'HashRef',  default  => sub { {} } );
has '_sequences' => ( is => 'rw', isa => 'ArrayRef' );
has 'output_filename' => ( is => 'ro', isa => 'Str', lazy => 1, builder => '_build_output_filename' );
has 'max_gap_between_genes'            => ( is => 'ro', isa => 'Int',     default => 20 );
has 'min_genes_on_contig'              => ( is => 'ro', isa => 'Int',     default => 5 );
has '_percentage_to_take_whole_contig' => ( is => 'ro', isa => 'Num',     default => 0.7 );
has 'min_sequence_length'              => ( is => 'ro', isa => 'Int',     default => 3000 ); #Phix 5386
has 'max_sequence_length'              => ( is => 'ro', isa => 'Int',     default => 1800000 );
has 'sequence_ids_to_genes'            => ( is => 'ro', isa => 'HashRef', lazy    => 1, builder => '_build_sequence_ids_to_genes' );
has '_blocks_to_sequences'             => ( is => 'ro', isa => 'HashRef', lazy    => 1, builder => '_build__blocks_to_sequences' );

sub _build_output_filename {
    my ($self) = @_;
    my ( $filename, $directories, $suffix ) = fileparse( $self->gff_file, qr/\.[^.]*/ );
    return $filename . '.novel.fa';
}

sub print_gene_products_on_contigs {
    my ($self) = @_;

    for my $seq_id ( keys %{ $self->sequence_ids_to_genes } ) {
        for my $gene ( @{ $self->sequence_ids_to_genes->{$seq_id} } ) {
            print $self->gff_file . "\t" . $seq_id . "\t" . $gene . "\t" . $self->genes_to_product->{$gene} . "\n";
        }
        print "\n";
    }
}

sub extract_nuc_sequences_from_blocks {
    my ($self) = @_;
    $self->_blocks_to_sequences;

    my $out_seq_io = Bio::SeqIO->new( -file => ">" . $self->output_filename, -format => 'Fasta' );
    my ( $sample_name, $directories, $suffix ) = fileparse( $self->gff_file, qr/\.[^.]*/ );

    for my $seq_obj ( @{ $self->_sequences } ) {
        next unless ( defined( $self->_blocks_to_sequences->{ $seq_obj->display_id } ) );

        for my $blocks ( @{ $self->_blocks_to_sequences->{ $seq_obj->display_id } } ) {
            next unless ( defined( $blocks->[0] ) && defined( $blocks->[1] ) );
            my $start_coord = $self->genes_to_feature->{ $blocks->[0] }->start;
            my $end_coord   = $self->genes_to_feature->{ $blocks->[0] }->end;

            $start_coord = $self->genes_to_feature->{ $blocks->[1] }->start
              if ( $self->genes_to_feature->{ $blocks->[1] }->start < $start_coord );
            $end_coord = $self->genes_to_feature->{ $blocks->[1] }->end if ( $self->genes_to_feature->{ $blocks->[1] }->end > $end_coord );

			my $whole_sequence_used = 0 ;
            # If the block makes up a big percentage of the contig, then just take the whole thing
            if ( ( $end_coord - $start_coord ) / $seq_obj->length > $self->_percentage_to_take_whole_contig ) {
                $start_coord = 1;
                $end_coord   = $seq_obj->length;
				$whole_sequence_used =1 ;
            }
			
			next if( ($end_coord - $start_coord)  < $self->min_sequence_length);
			next if( ($end_coord - $start_coord)  > $self->max_sequence_length);

            my $sample_name = join( '___', ( $sample_name, $seq_obj->display_id, $start_coord, $end_coord ) );
            $sample_name =~ s!\W!_!gi;
            $out_seq_io->write_seq( Bio::Seq->new( -display_id => $sample_name, -seq => $seq_obj->subseq( $start_coord, $end_coord ) ) );
			
			last if($whole_sequence_used == 1);
        }

    }
    return 1;
}

sub _build__blocks_to_sequences {
    my ($self) = @_;
    my %blocks;
    for my $seq_id ( keys %{ $self->sequence_ids_to_genes } ) {
        my $smallest_gene;
        my $previous_gene_number;
        my $largest_gene;

        # big assumption here is that the genes are sequentially numbered and that there is an underscore plus digits
        my @gene_ids = sort @{ $self->sequence_ids_to_genes->{$seq_id} };
        next if ( @gene_ids < $self->min_genes_on_contig );

        for ( my $i = 0 ; $i < @gene_ids ; $i++ ) {
            if ( $gene_ids[$i] =~ /\_([\d]+)$/ ) {
                my $gene_number = $1;

                if ( !defined($smallest_gene) ) {
                    $previous_gene_number = $gene_number;
                    $smallest_gene        = $gene_ids[$i];
                    next;
                }

                if ( $gene_ids[$i] ne $smallest_gene && ( $previous_gene_number + $self->max_gap_between_genes ) > $gene_number ) {
                    $largest_gene = $gene_ids[$i];
                }
                else {
                    my @current_block = ( $smallest_gene, $largest_gene );
                    push( @{ $blocks{$seq_id} }, \@current_block );

                    # new block
                    $smallest_gene        = undef;
                    $previous_gene_number = undef;
                    $largest_gene         = undef;
                }

                $previous_gene_number = $gene_number;
            }
        }
        if ( defined($largest_gene) && defined($smallest_gene) ) {
            my @current_block = ( $smallest_gene, $largest_gene );
            push( @{ $blocks{$seq_id} }, \@current_block );
        }

    }
    return \%blocks;

}

sub _build_sequence_ids_to_genes {
    my ($self) = @_;
    my %seq_ids_to_genes;
    my $gffio = Bio::Tools::GFF->new( -file => $self->gff_file, -gff_version => 3 );
    while ( my $feature = $gffio->next_feature() ) {
        my $gene_id = $self->_get_feature_id($feature);
        next unless ($gene_id);

        for my $unclassified_gene_id ( @{ $self->gene_ids } ) {
            if ( $gene_id eq $unclassified_gene_id ) {
                push( @{ $seq_ids_to_genes{ $feature->seq_id } }, $gene_id );
                $self->genes_to_feature->{$gene_id} = $feature;

                if ( $feature->has_tag('product') ) {
                    my ( $product, @junk ) = $feature->get_tag_values('product');
                    next unless ( defined($product) );
                    $self->genes_to_product->{$gene_id} = $product;
                }

                last;
            }
        }
    }
    my @seq_objects = $gffio->get_seqs();
    $self->_sequences( \@seq_objects );

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
