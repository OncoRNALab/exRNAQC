use strict;
use warnings;

use Cwd "abs_path";
use File::Basename;

# Get current working directory =================
my $scriptdir = "/user/data/gent/gvo000/gvo00027/TCGA_data/small_RNA_seq/finalscripts";

# Get parameters --------------------------------
my $starttime = time;
print "=================\n";
print "Script started on ".localtime($starttime)."\n\n";

# Initialize variables ===================================================

my $wdir			= $ARGV[0];
my $sample			= $ARGV[1];

my $error_flag		= 0;

my $read_count = $wdir."/".$sample."_read_count_new.txt";
open(FH_read_count, '>', $read_count) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$read_count."': $!";
	print STDERR $err_message;
	exit(1);
}

my $curr_time = time;
my $time_diff = $curr_time - $starttime;
print "[".sprintf("%.6f", $time_diff)."s]\n";
my $prev_time = $curr_time;

# Count total number of reads in original fastq ===================================================
print "Counting total number of reads in original fastq-file... ";

#my $fastq_file = $wdir."/".$sample."_R1.fastq.gz";
#my $fastq_count = get_fastq_read_count($fastq_file);
#print FH_read_count "fastq\t".$fastq_count."\n";

#$curr_time = time;
#$time_diff = $curr_time - $prev_time;
#print "[".sprintf("%.6f", $time_diff)."s]\n";
#$prev_time = $curr_time;

# Count total number of reads after clipping ======================================================
print "Counting total number of reads after clipping... ";

#my $clipped_fastq_file = $wdir."/".$sample."_clipped.fastq.gz";
#my $clipped_count = get_fastq_read_count($clipped_fastq_file);
#print FH_read_count "clipped\t".$clipped_count."\n";

#$curr_time = time;
#$time_diff = $curr_time - $prev_time;
#print "[".sprintf("%.6f", $time_diff)."s]\n";
#$prev_time = $curr_time;

# Count total number of reads after qc ============================================================
print "Counting total number of reads after QC-check... ";

#my $qc_clipped_fastq_file = $wdir."/".$sample."_qc_clipped.fastq.gz";
#my $qc_clipped_count = get_fastq_read_count($qc_clipped_fastq_file);
#print FH_read_count "qc_clipped\t".$qc_clipped_count."\n";

#$curr_time = time;
#$time_diff = $curr_time - $prev_time;
#print "[".sprintf("%.6f", $time_diff)."s]\n";
#$prev_time = $curr_time;

# Count total number of mapped reads ==============================================================
print "Counting total number of mapped reads... ";

my $sam_file = $wdir."/".$sample."_mapped.sam";
open (FH, '<', $sam_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$sam_file."': $!";
	print STDERR $err_message;
	exit(1);
}

my %reads = ();
my %stop_oligo_read_ids = ();
while (my $line = <FH>) {
	chomp $line;
	
	### Skip header
	if($line =~ /^@.+/){
		next;
	}
	my @bowtie_l = split("\t", $line);
	
	### Skip unmapped reads
#	if($bowtie_l[1] == 4){
#		next;
#	}
	
	### Get read_ids for reads mapping on stop-oligo
	if($bowtie_l[2] eq "stop_oligo"){
		my @tmp = split("-", $bowtie_l[0]);
		$stop_oligo_read_ids{$tmp[0]} ++;
	}
	
	$reads{$bowtie_l[0]} ++;
}
close(FH);

### Get unique read
my @stop_oligo_read_ids_unique = keys %stop_oligo_read_ids;
my @reads_unique = keys %reads;

### Make sum of mapped reads
my $total_read_count = 0;
my $stop_oligo_count = 0;
foreach my $item (@reads_unique){
	my @tmp = split("-", $item);
	
	# Reads mapped to stop-oligo
	if(grep {$_ eq $tmp[0]} @stop_oligo_read_ids_unique){
		$stop_oligo_count = $stop_oligo_count + $tmp[1];
		
	# Actual mapped reads
	} else {
		$total_read_count = $total_read_count + $tmp[1];
	}
}

print FH_read_count "technical\t".$stop_oligo_count."\n";
print FH_read_count "mapped\t".$total_read_count."\n";

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Count total number of mapped reads identifying miRNAs ===========================================
print "Counting total number of mapped reads identifying miRNAs... ";

my $mir_file = $wdir."/".$sample."_miRs.txt";
open (FH, '<', $mir_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$mir_file."': $!";
	print STDERR $err_message;
	exit(1);
}

my $total_miRNA_count = 0;
while (my $line = <FH>) {
	my @tmp = split("\t", $line);
	$total_miRNA_count = $total_miRNA_count + $tmp[2];
}
close(FH);

print FH_read_count "mapped_miR\t".$total_miRNA_count."\n";

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Count total number of mapped reads identifying contaminants (piRNAs, tRNAs, rRNAs) ===========================================
print "Counting total number of mapped reads identifying contaminants (piRNAs, tRNAs, rRNAs)... ";

my $contam_file = $wdir."/".$sample."_contam.txt";
open (FH, '<', $contam_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$contam_file."': $!";
	print STDERR $err_message;
	exit(1);
}

my $total_contam_count = 0;
while (my $line = <FH>) {
	my @tmp = split("\t", $line);
	$total_contam_count = $total_contam_count + $tmp[2];
}
close(FH);

print FH_read_count "mapped_contam\t".$total_contam_count."\n";

close(FH_read_count);

# Total running time ==============================================================================
print "\nScript finished on ".localtime($curr_time)."\n";
print "Running time: ";
my $total_time = $curr_time - $starttime;
print sprintf("%.6f", $total_time)."s\n";



sub get_fastq_read_count{
	my ($file) = @_;
	
	if(! -e $file || ! -r $file ){
		my $err_message = "Cannot open file '".$file."': $!";
		print STDERR $err_message;
		exit(1);
	}
	
	my $fastq_count = `zcat $file | grep "^@"|wc -l` or $error_flag = 1;
	if ($error_flag == 1) {
		my $err_message = "Cannot read counts in fastq-file ('".$file."'): $!";
		print STDERR $err_message;
		exit(1);
	}
	chomp $fastq_count;
	return $fastq_count;
}
