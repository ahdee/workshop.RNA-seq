use strict;
	
my $unstrand = 0; 	# most sequences this is 0 so becareful with this. 
my $post = "/data/post_analysis/";
my $fastin = "/data/RNA-Seq_Raw/"; 
my $gofast = 0; 
my %fast; 
# will store into fast has as basename_1 or _2 you can get this value whenever u want. 
if ($gofast == 1){
	opendir(DIR,$fastin);
	my @f = readdir(DIR); 
	foreach ( @f ) {
		my $file = $_; 
		next if $file !~ /(.+?)[-|_].*_R*(\d)_.+html/; 
		my $base = $1; 
		my $or = $2; 
		my $k = $base."_".$or; 
		my @fq_s = openfile($fastin.$file);
        my $temp = "@fq_s";
        $temp =~ /(Total Sequences).+?>(\d\d+)</; # Total Sequences</td><td>95376061</td>
 		$fast{$k} = $2; 
       		
		}

}


opendir(DIR,$post);

my $header = "sample,Total-reads,Mapped-reads,Unmapped-reads,Uniquely-mapped-reads,Multiple-mapped-reads,Unmapped-reads/Mapped-reads ratio, Uniquely-mapped-reads/Multiple-mapped-reads ratio,Q30-pct,Sense/anti-sense ratio,Orientation,Exonic/Intronic ratio,Library-size (Mil),Uniquely-Mapped-NonDupeReads,dataset\n"; 

my @dataset = readdir(DIR); 

 


#@dataset = ('ASC20', 'ASC21'); # you can overide it by looping through your own dataset
#@dataset = ('ASC15');
my $bigOut = $header; 


my $single_ds = 'ASC48';
#@dataset = ($single_ds);



my $all = 1; 
if ($all == 1){
	foreach ( @dataset ) {
	my $dataset = $_; 
	next if $dataset !~ /\w/; 
	next if $dataset !~ /ASC|^p/; 
	#print "$dataset \n ";
	if ($dataset){
		print "\n ----- $dataset -----\n";
		godoit($dataset); 
	} 
	}
}


print "done \n"; 
writeToo("$post".""."_summaryQC.csv",$bigOut,0); 



sub godoit{
	my $dataset = shift @_; 
	my $ds = "$post$dataset/"; # main directory - make sure to create a folder 'Star_counts' below this. 
	my $qc=$ds."QC/";
	my $counts=$ds."Star_counts/"; # below t
	my $file_out=$qc.$dataset."_QC.csv"; 
	my $umem=$ds."star.align/"; 

    return if ! -d $qc;
	opendir(DIR,$qc);

	# main header
	my $output = $header;  
	   

	QC: while (my $file = readdir(DIR)) {
			
			next if $file !~ /(.+?)\.stats.txt/;
			my $basename  = $1; 
			print "$basename \n";
			# open file and then parse
			# there should be a corresponding count data.
			# warn if count file does not exists.
			#print "$file \n"; 

			my $qc_stats = getQC($qc.$file);
			my $lib_size = ''; 
			my $umem_file = $umem.$basename.".QC.txt"; 
			my $umem_score = 'NA';  
			### check ummem scores
			
			if (-e $umem_file){
				$umem_score = getUmem($umem_file);  
				
			   }
			else {
				print "$umem_file NO UMEM data check! \n"; 
			}
			
			
			if (-e $counts.$basename."ReadsPerGene.out.tab"){
				$lib_size = getCounts($counts.$basename."ReadsPerGene.out.tab");  
			   }
			else {
				print "NO count data check!! ". $counts.$basename."ReadsPerGene.out.tab \n"; 
			}
			
			### these contains fastqc total reads if you want. 
			my $r1 = $fast{$basename."_1"}; 
		    my $r2 = $fast{$basename."_2"}; 
		
		
			$output .= "$basename,$qc_stats,$lib_size,$umem_score,$dataset\n"; 
			$bigOut .= "$basename,$qc_stats,$lib_size,$umem_score,$dataset\n";
			
		}


	print "$output \n";
	writeToo("$file_out",$output,0); 
    
	
	
	

}


	sub getCounts{
		my $file = shift;
		my @counts = openfile($file);
		my $sum = 0; my $sumSense = 0; 
		my $last = @counts * 1 -1; 
		foreach (4..$last ){
		 my @temp = split ("\t",$counts[$_]);
		 $sum += $temp[3];
         $sumSense += $temp[2];		 
		}
        $sum = $sumSense if $sumSense > $sum; 
		$sum = $sumSense  + $sum if $unstrand == 1; 
		return sprintf("%.1f", $sum/1000000);
		
	}

	
	sub getUmem{
		my $file = shift; 
		my @u = openfile($file);
		@u = split ("\t",$u[1]);
		return $u[2]; 
		
		
	}
	

	sub getQC{
		my $file = shift;
		my @qc_s = openfile($file);
		$file =~ /.+\/(.+?).stats.txt/;
		my $base = $1; 
		## loop through and find keywords for the following
		# Sample,Mapped-reads,unmapped/mapped,unique/multiple,unique/total,Q30,anti/sense,Orientation,exon/intron,library size (counts)'."\n"; 
		my ($total_reads, $Mapped_reads,$Unmapped_reads,$Uniquely_mapped_reads,$Multiple_mapped_reads,$Q30_pct,$anti2sense,$Orientation,$exon2intron); 
		foreach (@qc_s){
		    $total_reads = $1 if $_ =~ /Total-reads:\t([\d|\.]+)/;
			$Mapped_reads = $1 if $_ =~ /Mapped-reads:\t([\d|\.]+)/; 
			$Unmapped_reads =$1 if $_ =~ /Unmapped-reads:\t([\d|\.]+)/; 
			$Uniquely_mapped_reads  =$1 if $_ =~ /Uniquely-mapped-reads:\t([\d|\.]+)/; 
			$Multiple_mapped_reads  =$1 if $_ =~ /Multiple-mapped-reads:\t([\d|\.]+)/; 
			$Q30_pct  =$1 if $_ =~ /Q30-pct:\t([\d|\.]+)/; 
			$anti2sense =$1 if $_ =~ /Sense\/anti-sense ratio\s+(.+)/; 
			$Orientation =$1 if $_ =~ /Orientation\t([A-Z]+)/; 
			$exon2intron = $1 if $_ =~ /Exonic\/intronic ratio\s+(.+)/; 
		}
		# calcualte ratios 
		my $unmap2map = sprintf("%.2f", $Unmapped_reads / $Mapped_reads );
		my $unique2multiple = sprintf("%.2f", $Uniquely_mapped_reads / $Multiple_mapped_reads );
		my $unique2total = sprintf("%.2f", $Uniquely_mapped_reads / $Mapped_reads );
		
		#print "$base\n mapped reads: $Mapped_reads \n Unmap_reads: $Unmapped_reads \n unique mapped reads : $Uniquely_mapped_reads \n multiple reads: $Multiple_mapped_reads \n Q30: $Q30_pct \n antise-sense:  $anti2sense \n orientation: $Orientation \n exon2intron: $exon2intron \n ";
		#print "unmape2map: $unmap2map\n unique2multiple: $unique2multiple \n unique2total: $unique2total\n "; 
		
		return "$total_reads,$Mapped_reads,$Unmapped_reads,$Uniquely_mapped_reads,$Multiple_mapped_reads,$unmap2map,$unique2multiple,$Q30_pct,$anti2sense,$Orientation,$exon2intron";

	}



sub writeToo{
# file name, input, append (1) or not (0)
my $file = shift;
my $input = shift;
my $append = shift;

$append = 0 if !$append;

if ($append==0){$append = '>';}
else {$append = '>>';}

open (MYFILE, $append.$file) or die print "$! ERROR: $file \n";
print MYFILE "$input";
close (MYFILE);




}


sub openfile{
my $file = shift;
open (FILE, $file) or die print "cannot open $file \n";
my @content = <FILE>;
close FILE;

return @content; 
}