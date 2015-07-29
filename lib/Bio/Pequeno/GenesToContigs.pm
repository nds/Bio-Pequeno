package Bio::Pequeno::GenesToContigs;

# ABSTRACT:  Given a GFF file and a list of IDs extract contigs names

=head1 SYNOPSIS

Given a GFF file and a list of IDs extract contigs names

=cut

use Moose;
use Bio::Tools::GFF;

has 'gff_file'  => ( is => 'ro', isa => 'Str',  required => 1 );
has 'gene_ids'  => ( is => 'ro', isa => 'ArrayRef',  required => 1 );

has 'genes_to_product' => ( is => 'rw', isa => 'HashRef');
has 'max_gap_between_genes' => ( is => 'ro', isa => 'Int', required => 20);
has 'min_genes_on_contig' => ( is => 'ro', isa => 'Int', required => 2);

has 'sequence_ids_to_genes' => ( is => 'ro', isa => 'HashRef', lazy => 1, builder => '_build_sequence_ids_to_genes' );
has 'contigs' => ( is => 'ro', isa => 'HashRef', lazy => 1, builder => '' );

#loop over file to link gene_ids to sequenc_ids
#loop again - store all annonotation for a contig, flag which genes are unclassified
#add the sequence to a contig

sub calculate_blocks
{
	my ($self) = @_;
	my @blocks;
	for my $seq_id (keys %{$self->sequence_ids_to_genes})
	{
		my $smallest_gene;
		my $previous_gene_number;
		my $largest_gene;
		# big assumption here is that the genes are sequentially numbered and that there is an underscore plus digits
		my @gene_ids = sort @{$self->sequence_ids_to_genes->{$seq_id}};
		next if( @gene_ids < $self->min_genes_on_contig);
		
		for(my $i = 0; $i< @gene_ids; $i++ )
		{
			if( $gene_ids[$i] =~ /\_([\d]+)$/)
			{
			  my $gene_number = $1;	
			
			  if(!defined($smallest_gene))
			  {
			  	$previous_gene_number = $gene_number;
				$smallest_gene = $gene_ids[$i];
				next;
			  }

			  if( $gene_ids[$i] ne $smallest_gene && ($previous_gene_number + $self->max_gap_between_genes) > $gene_number)
			  {
			  	$largest_gene =  $gene_ids[$i];
			  }
			  else
			  {
				my @current_block = ($smallest_gene,$largest_gene);
				push(@blocks, \@current_block);
				# new block 
		  		$smallest_gene        = undef;
		  		$previous_gene_number = undef;
		  		$largest_gene         = undef;
			  }
				  
			  $previous_gene_number = $gene_number;
			}			
		}
		
	}
		
}

sub _build_sequence_ids_to_genes
{
	my ($self) = @_;
	my %seq_ids_to_genes;
    my $gffio = Bio::Tools::GFF->new( -file => $self->gff_file, -gff_version => 3 );
    while ( my $feature = $gffio->next_feature() ) {
        my $gene_id = $self->_get_feature_id($feature);
        next unless ($gene_id);
		
		for my $unclassified_gene_id( @{$self->gene_ids} )
		{
			if($gene_id eq $unclassified_gene_id)
			{
				push(@{$seq_ids_to_genes{$feature->seq_id}},$gene_id);
				
		        if ( $feature->has_tag('product') ) {
		            my ( $product, @junk ) = $feature->get_tag_values('product');
		            $self->genes_to_product->{$gene_id} = $product;
		        }
				
				last;
			}
		}
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
