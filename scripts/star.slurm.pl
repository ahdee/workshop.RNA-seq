
use strict; 

my $in_dir = '/scratch/users/tree101/human/';
my $a = "p2b";
 
my $mus = 0;  
my $star_index = '/srv/gsfs0/projects/sweet-cordero/alex/bin/index/STAR/hg38.p5/';
my $gtf = $star_index.'genecode_nopseudo.gtf';


if ($mus == 1){
	$star_index = '/srv/gsfs0/projects/sweet-cordero/alex/bin/index/STAR/mm10/';
	$gtf = $star_index.'gencode.vM13.annotation.gtf';
}

my %sample = ( );

####################### creating sh directory 
my $sh  = `pwd`;
$sh =~ s/[\n\r\f\t]//g;
$sh .= "/star.count/" ; 
if (! -e $sh  and ! -d $sh ) {
print " sh directory not exist so making with ($sh ) 755 \n";
mkdir $sh , oct('0755') or die print "\n $! #### Cannot make dir ($sh) ##### \n" ; 
}

my $dataset = $a; 
$in_dir = $in_dir.$a.'/'; 

my $go = 0; 

while ($go != 1){

    print "is this the correct input directory? \n===(--|--)===\n ($in_dir) if yes then type 'yes' else type full directory with trailing \/ \n===(--|--)===\n ";
	 print "  $star_index \/ \n===(--|--)===\n $gtf \n";
	my $a1 =   <STDIN>; 
	$a1  =~ s/[\n\r\f\t]//g;
	
	if ($a1 ne 'yes'){
		$in_dir = $a1; 
	}
	else{
	$go = 1; 
		 
		print "working on $a now \n"; sleep 1; 
	}
}



my $fastq = $in_dir."trim/";
my $star_out = $fastq.'star_aligned/'; 
my $reports = $star_out.'reports/'; 
my $error_name = 'st_';
my $qc = $star_out.'QC/'; 



### java programs  

my $picard = "/srv/gsfs0/projects/sweet-cordero/alex/bin/picard/picard-tools-2.2.4/"; 
my $ng = '/srv/gsfs0/projects/sweet-cordero/alex/bin/ngsutilsj';

# config scheduler
my ($slots,$l,$vmem,$email) = (8,'--partition=nih_s10 --time=12:00:00','--mem=71GB','tree101@stanford.edu');
my $jmem = 4 * 4;



my $c = "STAR --genomeDir $star_index --runThreadN $slots --readFilesIn IN --readFilesCommand zcat --outFileNamePrefix $star_out/BASENAME --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile $gtf" ;

my $c4 = "$ng bam-stats --silent --unique --gtf $gtf -o "; 

if (! -e $star_out  and ! -d $star_out ) {
print " output directory not exist so making with ($star_out ) 755 \n";
mkdir $star_out , oct('0755') or die print "\n $! #### Cannot make dir ($star_out) ##### \n" ; 
}

if (! -e $reports and ! -d $reports) {
print " output does not exist so making with ($reports) 755 \n";
mkdir $reports, oct('0755') or die print "\n #### $! Cannot make dir ($reports) ##### \n" ; 
}

if (! -e $qc and ! -d $qc) {
print " output does not exist so making with ($qc) 755 \n";
mkdir $qc, oct('0755') or die print "\n #### $! Cannot make dir ($qc) ##### \n" ; 
}

if (! -d $star_index) {
print " missing index: please check $star_index  \n";
exit; 
} 
 


opendir(DIR,$fastq);



my $log = 'module purge; module load STAR/2.5.1b; module load java/latest; module load samtools/1.5; alias STAR=\'/home/tree101/scratch/STAR/STAR/bin/Linux_x86_64/STAR\'; ' . "\n"; 
   $log = $log. 'PATH=/home/tree101/scratch/STAR/STAR/bin/Linux_x86_64:$PATH;' . "\n";  

   
   

D: while (my $file = readdir(DIR)) {
     # print "$file \n";
    next if $file !~ /_\d.*.f/i;  
	#next if $file !~ /U11/i;
	#print "** $file\n"; next; 
    if ($file =~ /(.+)_1\./){
        my $file2 = $file;
        my $basename = $1;
		 
        $file2 =~ s/_1\./_2\./;
		$log .= "#### for $basename #### \n"; 
		$file = $fastq.$file; 
		$file2 = $fastq.$file2; 
		#print "$file and $file2 \n";  next; 
		#next; 
		
		my $j1 = $c; 
		$j1 =~ s/IN/$file $file2/; 
		$j1 =~ s/BASENAME/$basename/;
		
		# error reports and name of run  
		my $j_name = "st_$basename"; 
		my $e = "-e $reports".'s1_'.$basename.".error -o ".$reports.'s1_'.$basename."_.o";
		## align
		
		
		#my  $submit = "qsub $e -N $j_name -V $vmem -m ea -M $email -pe shm $slots -b y '$j1';";
		# JOBID=$(sbatch my_job_script.sh | awk '{print $4}')
		
		my $submit = "JOBID=\$(sbatch $l --workdir=$star_out --account=ascor --export=ALL $e --job-name=$j_name $vmem -n $slots --wrap '$j1' | awk '{print \$4}')";
		
		
		
		#### log #system($submit) if ! -e "$star_out/$basename"."Aligned.sortedByCoord.out.bam"; 
		$log .= "$submit\n" if ! -e "$star_out/$basename"."Aligned.sortedByCoord.out.bam"; 
		# do QC ngutil
	
		my $tout = $basename."Aligned.sortedByCoord.out.bam"; 
		my $bam = "$star_out/$basename"."Aligned.sortedByCoord.out.bam";
        my $j4 = $c4.$qc."$basename.stats.txt ".$star_out.$tout; 
		
		my $loc_mem = 30; 
		$j4 = "java -Xmx$loc_mem"."g -jar $j4";
		
		my ($slots2,$l2,$vmem2,$email) = (1,'--partition=nih_s10 --time=12:00:00','--mem=30GB','tree101@stanford.edu');
		
		#my  $submit4 = "qsub -hold_jid '$j_name' -wd $star_out $e -N ng_$basename -V -l h_vmem=$loc_mem2 -m ea -M $email -pe shm 2 -b y '$j4';";
        my  $submit4 = "JOBID=\$(sbatch $l2 --dependency=afterany:\$JOBID --workdir=$star_out --account=ascor --export=ALL $e --job-name=ng$j_name $vmem2 -n $slots2 --wrap '$j4' | awk '{print \$4}');";		
		#### log ### system($submit4); 
		$log .= "$submit4\n" if ! -e "$star_out/$basename"."Aligned.sortedByCoord.out.bam"; ; 
		
		print "--- $basename \n";
		 
		 
		
		}
	
	}
	
my $fileO = $sh.$dataset.".align".time.".sh"; 
writeToo($fileO,$log);
#writeToo("$star_out"."star.log.txt",$log);
print "completed.. please see sh under star.count \nsh $fileO\n"; 
sub writeToo{
# file name, input, append (1) or not (0)
my $file = shift;
my $input = shift;
my $append = shift;

$append = 0 if !$append;

if ($append==0){$append = '>';}
else {$append = '>>';}

open (MYFILE, $append.$file);
print MYFILE "$input";
close (MYFILE);




}
