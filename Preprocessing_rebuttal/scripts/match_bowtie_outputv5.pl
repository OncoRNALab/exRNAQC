use strict;
use warnings;

use Cwd "abs_path";
use File::Basename;
use Time::HiRes qw(time);
use List::Util qw(min);
use Data::Dumper;

# Get current working directory -----------------
#my $scriptdir = "/user/data/gent/gvo000/gvo00027/RNA_seq_pipeline/resources_Annelien/";

# Get parameters --------------------------------
$0 =~ m/^.*\/(.*[.].*)$/;
my $starttime = time;
print "=================\n";
print "Script '".$1."' started on ".localtime($starttime)."\n\n";

# Initialize variables ===================================================
print "Initializing variables... ";

#my $project_name		= $ARGV[1];
#my $analysis_name		= $ARGV[2];
#my $ensembl_source		= $ARGV[7];

my $wdir				= $ARGV[0];
my $sample_name			= $ARGV[1];
my $organism			= "human";
my $mirbase_version		= "v22";
my $ensembl_version		= "v91";
my $ucsc_version		= "hg38";
my $error_flag			= 0;
my $annotation_path 	= "../../../resources/";

my $curr_time = time;
my $time_diff = $curr_time - $starttime;
print "[".sprintf("%.6f", $time_diff)."s]\n";
my $prev_time = $curr_time;

# Create UCSC-Ensembl chromosome conversion hash ==================================================
print "Creating UCSC-Ensembl chromosome conversion hash... ";
my $conversion_file = $annotation_path."Ensembl_" . $ensembl_version . "_chromosome_names.tsv";
open(FH, "<", $conversion_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$conversion_file."': $!";
	print STDERR $err_message;
	exit(1);
}

my $conv_header = <FH>;
my %chrom_conversion;
while (my $line = <FH>) {
	chomp $line;
	my (undef, $ucsc_id, $ensembl_id) = split("\t", $line);
	
	### If no UCSC chromosome
	if($ucsc_id eq "\\N"){
		next;
	}
	
	### If no Ensembl chromosome
	if($ensembl_id eq "\\N"){
		$chrom_conversion{$ucsc_id} = undef;
	} else {
		$chrom_conversion{$ucsc_id} = $ensembl_id;
	}
}

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Read miRbase file ===============================================================================
# Inlezen van miRbase-file met genomische coordinaten van mature en precursor miRNAs.
# Er worden 2 hashes aangemaakt, één voor mature miRNAs en één voor precursors.
# Keys van deze hashes zijn MIMAT- of MI-ids die op hun beurt gelinkt zijn aan chromosoom, start en stop via nieuwe hash (hash of hashes).
print "Reading miRbase file... ";

my %mature;
my %precursor;
my @ids;

my $mirbase_file = $annotation_path."annotations.gff3";
open(FH, "<", $mirbase_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$mirbase_file."': $!";
	print STDERR $err_message;
	exit(1);
}

while (my $line = <FH>) {
	chomp($line);
	$line =~ s/^\s+//;
	if($line =~ /^#/){
		next;
	}
	my ($chr, undef, $type, $start, $end, undef, $strand, undef, $attribute) = split("\t", $line);
	
	if(defined($chrom_conversion{$chr})){
		
		$chr = $chrom_conversion{$chr};
	
		if($type eq "miRNA_primary_transcript"){
			@ids	= split(";", $attribute);
			$ids[0]	=~ /MI[0-9]+/;
			
			$precursor{$chr}{$strand}{$&} = {'start'=> $start,
											 'end'	=> $end};
			
		} else {
			@ids	= split(";", $attribute);
			$ids[0]	=~ /MIMAT[0-9]+/;
			my $id_pre = $&;
			
			$ids[3]	=~ /MI[0-9]+/;
			my $id_suff	= $&;
			
			my $id_final = $id_pre."_".$id_suff;
			
			$mature{$chr}{$strand}{$id_final} = {'start'=> $start,
												 'end'	=> $end};
		}
	}
}
close(FH);


$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Read Ensembl contaminants file ==================================================================
# Inlezen van contaminants file met genomische coordinaten van snoRNA, snRNA, MT_tRNA, MT_rRNA, rRNA en miscRNA uit Ensembl.
# Er wordt 1 hash aangemaakt, keys van deze hash zijn contaminant-ids die op hun beurt gelinkt zijn aan chromosoom, start, stop, strand en type via nieuwe hash (hash of hashes).
print "Reading Ensembl contaminants file... ";

my %contam_ens;
my $suffix = 0;

my $ensembl_contaminants_file = $annotation_path."biomart_contaminants.txt";
open(FH, "<", $ensembl_contaminants_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$ensembl_contaminants_file."': $!";
	print STDERR $err_message;
	exit(1);
}

my $ens_header = <FH>;
while (my $line = <FH>) {
	$suffix ++;
	chomp($line);
	
	my ($gene_id, $start, $end, $chr, $strand, $gene_type) = split("\t", $line);
	
	if($strand eq "1"){
		$strand = "+";
		
	} elsif($strand eq "-1"){
		$strand = "-";
	}
	
	if(exists $contam_ens{$gene_id}){
		$gene_id = $gene_id."_".$suffix;
		$contam_ens{$chr}{$strand}{$gene_id} = {'start'	=> $start,
												'end'	=> $end,
												'type'	=> $gene_type};
	} else {
		$contam_ens{$chr}{$strand}{$gene_id} = {'start'	=> $start,
												'end'	=> $end,
												'type'	=> $gene_type};
	}
}
close(FH);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;


# Read UCSC tRNA contaminants file ================================================================
# Inlezen van contaminants file met genomische coordinaten van tRNA uit UCSC.
# Er wordt 1 hash aangemaakt, keys van deze hash zijn contaminant-ids die op hun beurt gelinkt zijn aan chromosoom, start, stop, strand en type via nieuwe hash (hash of hashes).
print "Reading UCSC tRNA contaminants file... ";

my %contam_ucsc;
$suffix = 0;

my $ucsc_contaminants_file = $annotation_path."contaminants_tRNA.gtf";
open(FH, "<", $ucsc_contaminants_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$ucsc_contaminants_file."': $!";
	print STDERR $err_message;
	exit(1);
}

while (my $line = <FH>) {
	$suffix ++;
	chomp($line);
	
	my ($chr, undef, undef, $start, $end, undef, $strand, undef, $attribute) = split("\t", $line);
	
	my @ids = split(";", $attribute);
	$ids[1] =~ /tRNA\w+-\w+/;
	my $key = $&;
	
	if(defined($chrom_conversion{$chr})){
		
		$chr = $chrom_conversion{$chr};
		
		if(exists $contam_ucsc{$key}){
			$key = $key."_".$suffix;
			
			$contam_ucsc{$chr}{$strand}{$key} = {'start'=> $start,
												 'end'	=> $end};
												
		} else {
			$contam_ucsc{$chr}{$strand}{$key} = {'start'=> $start,
												 'end'	=> $end};
												
		}
	}
}
close(FH);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Read UCSC piRNA contaminants file ===============================================================
# Inlezen van contaminants file met genomische coordinaten van piRNAs uit UCSC.
# Er wordt 1 hash aangemaakt, keys van deze hash zijn contaminant ids die op hun beurt gelinkt zijn aan chromosoom, start, stop, strand en type via nieuwe hash (hash of hashes).
print "Reading UCSC piRNA contaminants file... ";

my %piRNA = ();
$suffix = 0;

my $ucsc_pirna_file = $annotation_path."contaminants_piRNA.gtf";
open(FH, "<", $ucsc_pirna_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$ucsc_pirna_file."': $!";
	print STDERR $err_message;
	exit(1);
}

while (my $line = <FH>) {
	$suffix ++;
	chomp($line);
	
	my ($chr, undef, undef, $start, $end, undef, $strand, undef, $attribute) = split("\t", $line);
	
	@ids = split(";", $attribute);
	$ids[1]	=~ /uc.+[^"]/;
	my $key	= $&;
	
	if(defined($chrom_conversion{$chr})){
		
		$chr = $chrom_conversion{$chr};
		
		if(exists $piRNA{$key}){
			$key = $key."_".$suffix;
			
			$piRNA{$chr}{$strand}{$key} = {'start'=> $start,
										   'end'  => $end};
		
		} else {
			$piRNA{$chr}{$strand}{$key} = {'start'=> $start,
										   'end'  => $end};
		}
	}
}
close(FH);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Read Bowtie results =============================================================================
print "Reading Bowtie results... ";

my @mature_mir_reads	= ();
my %mature_mir_index	= ();
my @contam_reads		= ();
my %contam_index		= ();
my @not_annot_reads		= ();
my %not_annot_index		= ();
my @tech_artifacts		= ();
my @read_id_stop_oligo	= ();
my $i					= 0;

my $sam_file = $wdir."/".$sample_name."_mapped.sam";
open(FH, "<", $sam_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$sam_file."': $!";
	print STDERR $err_message;
	exit(1);
}

while (my $line = <FH>) {
	$i ++;
	if($line =~ /^@/){
		next;
	}
	chomp $line;
	my @bowtie_l	= split("\t", $line);
	my @tmp			= split("-", $bowtie_l[0]);
	my $read_id		= $tmp[0];
	my $read_count	= $tmp[1];
	my $read_start	= $bowtie_l[3] + 1;
	my $read_end	= $bowtie_l[3] + length($bowtie_l[4]);
	# my $read_chr	= substr($bowtie_l[2],3); #remove chr notation
	my $read_chr	= $bowtie_l[2]; #dont remove chr notation
	#my $read_chr	= $bowtie_l[2];
	my $read_strand = $bowtie_l[1];
	my $read_seq	= $bowtie_l[4];
	
	#if($read_strand == 0){
		#$read_strand = "+";
	#} elsif ($read_strand == 16){
	#	$read_strand = "-";
	if ($read_strand eq "-"){		
		# If strand is reverse, calculate reverse complement of sequence
		my $rev_comp_seq = reverse($read_seq);
		$rev_comp_seq =~ tr/ACGTacgt/TGCAtgca/;
		$read_seq = $rev_comp_seq;
	}	
	#} elsif ($read_strand == 4){
	#	next;
	#}

	
	my @keys_mature;
	my @keys_contam_ens;
	my @keys_contam_ucsc;
	my @keys_piRNA;
	if(exists($mature{$read_chr}{$read_strand})){
		@keys_mature			= grep {$mature		{$read_chr}{$read_strand}{$_}{'start'}	<= $read_end &&
										#$mature		{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys $mature		{$read_chr}{$read_strand};
										$mature		{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys %{$mature{$read_chr}{$read_strand}};
	}
	
	if(exists($contam_ens{$read_chr}{$read_strand})){
		@keys_contam_ens		= grep {$contam_ens	{$read_chr}{$read_strand}{$_}{'start'}	<= $read_end &&
										#$contam_ens	{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys $contam_ens	{$read_chr}{$read_strand};
										$contam_ens	{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys %{$contam_ens{$read_chr}{$read_strand}};
	}
	
	if(exists($contam_ucsc{$read_chr}{$read_strand})){
		@keys_contam_ucsc	= grep {$contam_ucsc{$read_chr}{$read_strand}{$_}{'start'}	<= $read_end &&
										#$contam_ucsc{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys $contam_ucsc	{$read_chr}{$read_strand};
										$contam_ucsc{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys %{$contam_ucsc{$read_chr}{$read_strand}};
	}
	
	if(exists($piRNA{$read_chr}{$read_strand})){
		@keys_piRNA			= grep {$piRNA		{$read_chr}{$read_strand}{$_}{'start'}	<= $read_end &&
										#$piRNA		{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys $piRNA			{$read_chr}{$read_strand};
										$piRNA		{$read_chr}{$read_strand}{$_}{'end'}	>= $read_start} keys %{$piRNA{$read_chr}{$read_strand}};
	}
	
	# Collapse isomiRs, only report isomiRs: per read check offset and check if only miR is mapped and no contaminants
	my $size_keys_mature	  = @keys_mature;
	my $size_keys_contam_ens  = @keys_contam_ens;
	my $size_keys_contam_ucsc = @keys_contam_ucsc;
	my $size_keys_piRNA		  = @keys_piRNA;

	### Mature miRNA (miRBase)
	if($size_keys_contam_ens == 0 && $size_keys_contam_ucsc == 0 && $size_keys_piRNA == 0 && $size_keys_mature > 0){
		my $offset5 = $mature{$read_chr}{$read_strand}{$keys_mature[0]}{'start'} - $read_start;
		my $offset3 = $read_end - $mature{$read_chr}{$read_strand}{$keys_mature[0]}{'end'};
		my %tmp = (
			'read_id'		=> $read_id,
			'read_count'	=> $read_count,
			'transcript_id'	=> $keys_mature[0],
			'read_strand'	=> $read_strand,
			'offset5'		=> $offset5,
			'offset3'		=> $offset3,
			'read_seq'		=> $read_seq,
		);
		my $tmp_index = scalar(@mature_mir_reads);
		push(@mature_mir_reads, \%tmp );
		push(@{$mature_mir_index{$read_id}}, $tmp_index);
		
	### Contaminants (Ensembl)
	} elsif($size_keys_contam_ens > 0){
		my $offset5 = $contam_ens{$read_chr}{$read_strand}{$keys_contam_ens[0]}{'start'} - $read_start;
		my $offset3 = $read_end - $contam_ens{$read_chr}{$read_strand}{$keys_contam_ens[0]}{'end'};
		my %tmp = (
			'read_id'		=> $read_id,
			'read_count'	=> $read_count,
			'transcript_id'	=> $keys_contam_ens[0],
			'read_strand'	=> $read_strand,
			'offset5'		=> $offset5,
			'offset3'		=> $offset3,
			'annot_type'	=> $contam_ens{$read_chr}{$read_strand}{$keys_contam_ens[0]}{'type'},
		);
		my $tmp_index = scalar(@contam_reads);
		push(@contam_reads, \%tmp );
		push(@{$contam_index{$read_id}}, $tmp_index);
	### tRNA contaminants (UCSC)
	} elsif($size_keys_contam_ucsc > 0){
		my $offset5 = $contam_ucsc{$read_chr}{$read_strand}{$keys_contam_ucsc[0]}{'start'} - $read_start;
		my $offset3 = $read_end - $contam_ucsc{$read_chr}{$read_strand}{$keys_contam_ucsc[0]}{'end'};
		my %tmp = (
			'read_id'		=> $read_id,
			'read_count'	=> $read_count,
			'transcript_id'	=> $keys_contam_ucsc[0],
			'read_strand'	=> $read_strand,
			'offset5'		=> $offset5,
			'offset3'		=> $offset3,
			'annot_type'	=> "tRNA",
		);
		my $tmp_index = scalar(@contam_reads);
		push(@contam_reads, \%tmp );
		push(@{$contam_index{$read_id}}, $tmp_index);
		
	### piRNA contaminants (UCSC)
	} elsif($size_keys_piRNA > 0){
		my $offset5 = $piRNA{$read_chr}{$read_strand}{$keys_piRNA[0]}{'start'} - $read_start;
		my $offset3 = $read_end - $piRNA{$read_chr}{$read_strand}{$keys_piRNA[0]}{'end'};
		my %tmp = (
			'read_id'		=> $read_id,
			'read_count'	=> $read_count,
			'transcript_id'	=> $keys_piRNA[0],
			'read_strand'	=> $read_strand,
			'offset5'		=> $offset5,
			'offset3'		=> $offset3,
			'annot_type'	=> "piRNA",
		);
		my $tmp_index = scalar(@contam_reads);
		push(@contam_reads, \%tmp );
		push(@{$contam_index{$read_id}}, $tmp_index);
		
	} elsif($size_keys_contam_ens == 0 && $size_keys_contam_ucsc == 0 && $size_keys_piRNA == 0 && $size_keys_mature == 0){
		
		### Filter reads that map to stop oligo
		if($read_chr eq "stop_oligo"){
			my %tmp = (
				'read_id'	  => $read_id,
				'read_count'  => $read_count,
				'read_seq'	  => $read_seq,
				'read_chr'	  => $read_chr,
				'read_start'  => $read_start,
				'read_end'	  => $read_end,
				'read_strand' => $read_strand,
			);
			push(@tech_artifacts, \%tmp );
			push(@read_id_stop_oligo, $read_id);
			
		### Not annotated reads
		} else {
			my %tmp = (
				'read_id'	  => $read_id,
				'read_count'  => $read_count,
				'read_seq'	  => $read_seq,
				'read_chr'	  => $read_chr,
				'read_start'  => $read_start,
				'read_end'	  => $read_end,
				'read_strand' => $read_strand,
			);
			my $tmp_index = scalar(@not_annot_reads);
			push(@not_annot_reads, \%tmp );
			push(@{$not_annot_index{$read_id}}, $tmp_index);
		}
	}

}
close(FH);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Remove reads that map on stop_oligo from all other categories ===================================
print "Removing reads that map on stop_oligo from all other categories... ";

foreach my $read_id (@read_id_stop_oligo) {
	
	### Removing reads from mature miRNA
	if(defined $mature_mir_index{$read_id}){
		delete(@mature_mir_reads[@{$mature_mir_index{$read_id}}]);
	}
	
	### Removing reads from contaminants
	if(defined $contam_index{$read_id}){
		delete(@contam_reads[@{$contam_index{$read_id}}]);
	}
	
	### Removing reads from non-annotated reads
	if(defined $not_annot_index{$read_id}){
		delete(@not_annot_reads[@{$not_annot_index{$read_id}}]);
	}
}

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Write technical artifacts to file ===============================================================
print "Writing technical artifacts to file... ";

my $sample_tech = $sample_name."_technical_artifacts.txt";
open(FH_tech, ">", $sample_tech) and truncate(FH_tech, 0) or $error_flag = 1;

if ($error_flag == 1) {
	my $err_message = "Cannot open or clear file '".$sample_tech."': $!";
	print STDERR $err_message;
	exit(1);
}

foreach my $read (@tech_artifacts){
	print FH_tech $read->{'read_id'}."\t".$read->{'read_count'}."\t".$read->{'read_seq'}."\t".$read->{'read_chr'}."\t".$read->{'read_start'}."\t".$read->{'read_end'}."\t".$read->{'read_strand'}."\n";
}

close(FH_tech);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Write non-annotated reads to file ===============================================================
print "Writing non-annotated reads to file... ";

my $sample_not_annot = $sample_name."_not_annotated.txt";
open(FHno_annot, ">", $sample_not_annot) and truncate(FHno_annot, 0) or $error_flag = 1;

if ($error_flag == 1) {
	my $err_message = "Cannot open or clear file '".$sample_not_annot."': $!";
	print STDERR $err_message;
	exit(1);
}

foreach my $read (@not_annot_reads){
	if( ! defined($read)){
		next;
	}
	print FHno_annot $read->{'read_id'}."\t".$read->{'read_count'}."\t".$read->{'read_seq'}."\t".$read->{'read_chr'}."\t".$read->{'read_start'}."\t".$read->{'read_end'}."\t".$read->{'read_strand'}."\n";
}

close(FHno_annot);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Select unique reads that mapped on contam =======================================================
# For each unique read, write type to file. If multiple types, write 'multiple' to file
print "Selecting unique reads that mapped on contam... ";

my $sample_contam_file = $sample_name."_contam.txt";
open(FHcontam, ">", $sample_contam_file) and truncate(FHcontam, 0) or $error_flag = 1;

if ($error_flag == 1) {
	my $err_message = "Cannot open or clear file '".$sample_contam_file."': $!";
	print STDERR $err_message;
	exit(1);
}

### Create hash of reads grouped by read-ID
my %unique_contam_read = ();
foreach my $read (@contam_reads){
	if( ! defined($read)){
		next;
	}
	
	if(!defined($unique_contam_read{$read->{'read_id'}})){
		$unique_contam_read{$read->{'read_id'}} = [];
	}
	
	push(@{$unique_contam_read{$read->{'read_id'}}}, $read);
}

foreach my $read_id (keys %unique_contam_read){
	my @reads_selected = @{$unique_contam_read{$read_id}};
	
	if (scalar @reads_selected == 1){
		print FHcontam $sample_name."\t".$reads_selected[0]->{'read_id'}."\t".$reads_selected[0]->{'read_count'}."\t".$reads_selected[0]->{'annot_type'}."\t".$reads_selected[0]->{'transcript_id'}."\n";
		
	} else {
		### Get all annotation types for this read_id. If multiple, join them together
		my %unique_contam = ();
		foreach my $item (@reads_selected){
			$unique_contam{$item->{'annot_type'}} ++;
		}
		
		if(keys %unique_contam == 1){
			print FHcontam $sample_name."\t".$reads_selected[0]->{'read_id'}."\t".$reads_selected[0]->{'read_count'}."\t".$reads_selected[0]->{'annot_type'}."\t".$reads_selected[0]->{'transcript_id'}."\n";
			
		} else {
			my $type_multiple = join("_", keys %unique_contam);
			print FHcontam $sample_name."\t".$reads_selected[0]->{'read_id'}."\t".$reads_selected[0]->{'read_count'}."\t"."multiple_".$type_multiple."\t"."none"."\n";
		}
	}
}
close(FHcontam);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Select uniqe reads that mapped on miRNA =========================================================
# For each unique read, check if mapped miR has same mimat and same offset, if different mimat, calculate sumof absolute offsets per mimat and include miR with lowest offset_sum
print "Selecting uniqe reads that mapped on miRNA... ";

my $sample_isomir_file = $sample_name."_isomiRs.txt";
open(FHisomir, ">", $sample_isomir_file) and truncate(FHisomir, 0) or $error_flag = 1;

if ($error_flag == 1) {
	my $err_message = "Cannot open or clear file '".$sample_isomir_file."': $!";
	print STDERR $err_message;
	exit(1);
}

### Create hash of reads grouped by read-ID
my %unique_mature_read = ();
foreach my $read (@mature_mir_reads){
	if( ! defined($read)){
		next;
	}
	
	if(!defined($unique_mature_read{$read->{'read_id'}})){
		$unique_mature_read{$read->{'read_id'}} = [];
	}
	
	push(@{$unique_mature_read{$read->{'read_id'}}}, $read);
	
}

my @isomirs = ();
foreach my $read_id (keys %unique_mature_read){
	my $check_id		= 0;
	my $check_offset	= 1;
	my @reads_selected	= @{$unique_mature_read{$read_id}};
	
	if(scalar @reads_selected == 1){
		my $offset5;
		my $offset3;
		
		### If strand is '-', Switch 3' and 5' offset
		if($reads_selected[0]->{'read_strand'} eq "-"){
			$offset5 = $reads_selected[0]->{'offset3'};
			$offset3 = $reads_selected[0]->{'offset5'};
			
		} else {
			$offset5 = $reads_selected[0]->{'offset5'};
			$offset3 = $reads_selected[0]->{'offset3'};
		}
		
		### Filter mature-ID from the string 'mature-ID_precursor'
		$reads_selected[0]->{'transcript_id'} =~ /MIMAT[0-9]+/;
		my $mature_mir_id = $&;
		
		print FHisomir $sample_name."\t".$reads_selected[0]->{'read_id'}."\t".$reads_selected[0]->{'read_count'}."\t".$mature_mir_id."\t".$offset5."\t".$offset3."\t".$reads_selected[0]->{'read_seq'}."\n";
		
		my %tmp = (
			'read_id'		=> $reads_selected[0]->{'read_id'},
			'read_count'	=> $reads_selected[0]->{'read_count'},
			'transcript_id'	=> $mature_mir_id,
			'offset5'		=> $offset5,
			'offset3'		=> $offset3,
		);
		
		push(@isomirs, \%tmp);
		
	} else {
		my %unique_mature_mir = ();
		foreach my $read (@reads_selected){
			$read->{'transcript_id'} =~ /MIMAT[0-9]+/;
			$unique_mature_mir{$&} ++;
		}
		
		my @read1;
		my $read1_offset5;
		my $read1_offset3;
		
		if(keys %unique_mature_mir == 1){
			$check_id = 1;
		
			### Loop over isomiRs and check if 5' and 3' offsets are all the same
			for($i = 1; $i < scalar @reads_selected; $i ++){
				my $read2_offset5;
				my $read2_offset3;
				
				@read1 = $reads_selected[$i-1];
				my @read2 = $reads_selected[$i];
				
				### Switch 3' and 5' offset for read1 if strand is reverse
				if($read1[0]->{'read_strand'} eq "-"){
					$read1_offset5 = $read1[0]->{'offset3'};
					$read1_offset3 = $read1[0]->{'offset5'};
					
				} else {
					$read1_offset5 = $read1[0]->{'offset5'};
					$read1_offset3 = $read1[0]->{'offset3'};
					
				}
			
				### Switch 3' and 5' offset for read2 if strand is reverse
				if($read2[0]->{'read_strand'} eq "-"){
					$read2_offset5 = $read2[0]->{'offset3'};
					$read2_offset3 = $read2[0]->{'offset5'};
					
				} else {
					$read2_offset5 = $read2[0]->{'offset5'};
					$read2_offset3 = $read2[0]->{'offset3'};
					
				}
				
				if($read1_offset5 != $read2_offset5 || $read1_offset3 != $read2_offset3){
					$check_offset = 0;
					last;
					
				}
			}
			
		} else {
			
			my @sum_abs_offset;
			### Get total offset for each read
			foreach my $read (@reads_selected){
				push(@sum_abs_offset, abs($read->{'offset5'}) + abs($read->{'offset3'}));
			}
			
			### Get minimal offset and frequency for each offset
			my $min_offset = min(@sum_abs_offset);
			my %freq;
			foreach (@sum_abs_offset){
				$freq{$_} ++;
			}
			
			my $array_index = 0;
			if($freq{$min_offset} == 1){
				until($sum_abs_offset[$array_index] == $min_offset){
					$array_index++;
				}
				
				@read1 = $reads_selected[$array_index];
				
				### Switch 3' and 5' offset for read2 if strand is reverse
				if($reads_selected[$array_index]->{'read_strand'} eq "-"){
					$read1_offset5 = $reads_selected[$array_index]->{'offset3'};
					$read1_offset3 = $reads_selected[$array_index]->{'offset5'};
				} else {
					$read1_offset5 = $reads_selected[$array_index]->{'offset5'};
					$read1_offset3 = $reads_selected[$array_index]->{'offset3'};
				}
				
				$check_id = 1;
			}
		}
		
		if($check_id == 1 && $check_offset == 1){
			$read1[0]->{'transcript_id'} =~ /MIMAT[0-9]+/;
			my $mature_mir_id = $&;
			
			print FHisomir $sample_name."\t".$read1[0]->{'read_id'}."\t".$read1[0]->{'read_count'}."\t".$mature_mir_id."\t".$read1_offset5."\t".$read1_offset3."\t".$read1[0]->{'read_seq'}."\n";
			
			my %tmp = (
				'read_id'		=> $read1[0]->{'read_id'},
				'read_count'	=> $read1[0]->{'read_count'},
				'transcript_id'	=> $mature_mir_id,
				'offset5'		=> $read1_offset5,
				'offset3'		=> $read1_offset3,
			);
			
			push(@isomirs, \%tmp);
		}
	}
}
close(FHisomir);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Collapse isomiR file to mature miRNA file =======================================================
print "Collapsing isomiR file to mature miRNA file... ";

my $sample_mir_file = $sample_name."_miRs.txt";
open(FHmir, ">", $sample_mir_file) and truncate(FHmir, 0) or $error_flag = 1;

if ($error_flag == 1) {
	my $err_message = "Cannot open or clear file '".$sample_mir_file."': $!";
	print STDERR $err_message;
	exit(1);
}

### Create hash of reads grouped by read-ID
my %unique_mir = ();
foreach my $read (@isomirs){
	if(!defined($unique_mir{$read->{'transcript_id'}})){
		$unique_mir{$read->{'transcript_id'}} = [];
	}
	
	push(@{$unique_mir{$read->{'transcript_id'}}}, $read);
}

my @mature_mirs = ();
foreach my $mimat_id (keys %unique_mir){
	my @reads_selected = @{$unique_mir{$mimat_id}};
	
	### Get the total count for each mimat-ID
	my $mimat_count = 0;
	foreach my $isomir (@reads_selected){
		$mimat_count = $mimat_count + $isomir->{'read_count'};
	}
	print FHmir $sample_name."\t".$mimat_id."\t".$mimat_count."\n";
	my @tmp = ($mimat_id, $mimat_count);
	push(@mature_mirs, join("\t", @tmp));
}
close(FHmir);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Total running time ==============================================================================
print "\nScript finished on ".localtime($curr_time)."\n";
print "Running time: ";
my $total_time = $curr_time - $starttime;
print sprintf("%.6f", $total_time)."s\n";
