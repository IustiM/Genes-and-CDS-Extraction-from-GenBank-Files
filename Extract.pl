#M.I.
use strict;
use warnings;
use Bio::SeqIO;

sub extract_genes_and_cds {
    my ($input_file, $output_file, $start, $end) = @_;

    # Convertim $start si $end la numere
    $start = int($start);
    $end = int($end);

    # Deschidem fisierul de intrare pentru citire header
    open(my $gb_file, "<", $input_file) or die "Nu se poate deschide $input_file: $!";
    my @header;
    while (my $line = <$gb_file>) {
        push @header, $line;
        last if $line =~ /^FEATURES/;
    }
    close($gb_file);

    # Citim fisierul GenBank
    my $seqio = Bio::SeqIO->new(-file => $input_file, -format => 'genbank');
    my $record = $seqio->next_seq;

    # Deschidem fisierul de iesire pentru scriere#M.I.
    open(my $out_file, ">", $output_file) or die "Nu se poate deschide $output_file: $!";

    # Scriem header-ul in fisierul de iesire
    print $out_file @header;

    # Initializam dictionare pentru a stoca informatiile gene si CDS
    my %gene_info_dict;
    my %cds_info_dict;

    # Iteram prin toate caracteristicile (features) din fisierul GenBank
    for my $feature ($record->get_SeqFeatures) {
        my $feature_start = $feature->location->start;
        my $feature_end = $feature->location->end;

        # Filtram doar genele si CDS-urile care se afla in intervalul specificat
        if (($feature_start >= $start && $feature_start <= $end) ||
            ($feature_end >= $start && $feature_end <= $end)) {
            if ($feature->primary_tag eq "gene") {
                my $gene_info = extract_feature_info($feature, $record);
                my $gene_id = $gene_info->{gene_id};
                $gene_info_dict{$gene_id} = $gene_info;
            } elsif ($feature->primary_tag eq "CDS") {
                my $cds_info = extract_feature_info($feature, $record);
                my $gene_id = $cds_info->{gene_id};
                $cds_info_dict{$gene_id} = $cds_info;
            }
        }
    }

    # Scriem informatiile despre gene si CDS-uri in fisierul de iesire
    for my $gene_id (keys %gene_info_dict) {
        my $gene_info = $gene_info_dict{$gene_id};
		
        # Scriem informatiile despre gene#M.I.
        print $out_file "Gene: $gene_info->{location}\n";
        print $out_file " /gene=$gene_id\n";
        print $out_file " /locus_tag=\"$gene_info->{locus_tag}\"\n" if $gene_info->{locus_tag};
        print $out_file " /note=\"$gene_info->{note}\"\n" if $gene_info->{note};
        print $out_file " /sequence:\n$gene_info->{sequence}\n\n";

        # Scriem informatiile despre CDS daca exista
        if (exists $cds_info_dict{$gene_id}) {
            my $cds_info = $cds_info_dict{$gene_id};
            print $out_file "CDS: $cds_info->{location}\n";
            print $out_file " /gene=$gene_id\n";
            print $out_file " /protein_id=\"$cds_info->{protein_id}\"\n" if $cds_info->{protein_id};
            print $out_file " /note=\"$cds_info->{note}\"\n" if $cds_info->{note};
            print $out_file " /sequence:\n$cds_info->{sequence}\n\n";
        }
    }

    close($out_file);
    print "Procesul de extragere si scriere a fost finalizat cu succes de catre Mutescu Iustinian_GM1!\n";
}

sub extract_feature_info {
    my ($feature, $record) = @_;
    my %info;

    # Extragem informatiile relevante din feature
    $info{location} = $feature->location->to_FTstring;
    for my $tag ($feature->get_all_tags) {
        if ($tag eq "gene") {
            ($info{gene_id}) = $feature->get_tag_values($tag);
        } elsif ($tag eq "locus_tag") {
            ($info{locus_tag}) = $feature->get_tag_values($tag);
        } elsif ($tag eq "note") {
            ($info{note}) = $feature->get_tag_values($tag);
        } elsif ($tag eq "protein_id") {
            ($info{protein_id}) = $feature->get_tag_values($tag);
        }
    }

    # Extragem secventa ADN pentru gene si CDS
    if ($feature->primary_tag eq "gene" || $feature->primary_tag eq "CDS") {
        my $seq = $feature->spliced_seq->seq;
        $info{sequence} = $seq;
    }

    return \%info;
}

# Specificam fisierul de iesire
my $output_file = "C:\\Users\\IUSTI\\Desktop\\output.txt";

# Apelam functia pentru extragerea genelor si CDS-urilor, utilizand argumentele din linia de comanda#M.I.
extract_genes_and_cds($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);