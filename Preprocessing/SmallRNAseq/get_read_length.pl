use strict;
use warnings;

use Cwd "abs_path";
use File::Basename;
use Bio::SeqIO;

# Get current working directory -----------------
my $scriptdir = "/user/data/gent/gvo000/gvo00027/TCGA_data/small_RNA_seq/finalscripts";

# Get parameters --------------------------------
my $starttime = time;
print "=================\n";
print "Script started on ".localtime($starttime)."\n\n";


# Initialize variables --------------------------
print "Initializing variables... ";
my $wdir			= $ARGV[0];
my $sample_name			= $ARGV[1];

my $error_flag			= 0;

my $curr_time = time;
my $time_diff = $curr_time - $starttime;
print "[".sprintf("%.6f", $time_diff)."s]\n";
my $prev_time = $curr_time;

# Create SeqIO-object ---------------------------
print "Creating SeqIO-object... ";

my $sequence_object;
my $sample_fasta = $wdir."/".$sample_name."_collapse.fa"; #make sure no additional suffixes are present in sample_name
#my $sample_fasta = $wdir."/".$sample_name."_qc_collapse.fa";
my $seqio = Bio::SeqIO ->  new ('-format' => 'fasta', 
								'-file'	  => $sample_fasta);
$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Write to file ---------------------------------
print "Writing to file... ";

my $read_length_file = $wdir."/".$sample_name."_read_length_new.txt";
open(FH, ">", $read_length_file) or $error_flag = 1;
if ($error_flag == 1) {
	my $err_message = "Cannot open file '".$read_length_file."': $!";
	print STDERR $err_message;
	exit(1);
}

print FH "read_id\tread_length\tread_count\n";
                                
while ($sequence_object = $seqio -> next_seq) { 
	my $read_length = $sequence_object -> length();
	my $id			 = $sequence_object -> id();
	my @id_items	 = split("-", $id);
	
	print FH $id_items[0]."\t".$read_length."\t".$id_items[1]."\n"; 
}
close(FH);

$curr_time = time;
$time_diff = $curr_time - $prev_time;
print "[".sprintf("%.6f", $time_diff)."s]\n";
$prev_time = $curr_time;

# Total running time ----------------------------
print "\nScript finished on ".localtime($curr_time)."\n";
print "Running time: ";
my $total_time = $curr_time - $starttime;
print sprintf("%.6f", $total_time)."s\n";
