# login to HPC
# clone git clone https://github.com/ahdee/workshop.RNA-seq.git
cd workshop.RNA-seq
pwd
ls
man ls # get info about command.  Sometimes just running the program is enough
# to end a command quickly do c+ ESC key 
history
clear
# use this \ to write commands multiple lines
echo \
"hello world"
ls -lS # list and sort by size based on manual 
# lets try to change some permission
cd ./bin
./tabutils # why is there an error
stat ./tabutils # find more info 
chmod 755 ./tabutils # lets change its permission
stat ./tabutils # how does it look like now
./tabutils # good to go 
# some bit of navigating around
cd ~/ # this is your home directory
cd bin # this is your executables
# when lost go back to ~/ or 
cd - # go back to original 
# lets play around a bit 
cd ~/workshop.RNA-seq
mkdir practice # create a new directory
cd practice
cp ../full_length.tab ./ # copy a file: recall . vs ..
ls –lS
mv ./full_length.tab ./fun.tab # rename or move 
# copying is fine however this could be a bit troublesome for BIG files
# use rysnc instead 
# extra benefit can copy files from outside 
rsync -avP /asclab/data1/backup/projects/spcg/working/alex.rna/DATA/UCSF_747/BC6/raw/* ./
# if external then normally its username@<host>:<dir/target>
cd ..
rm –rf practice # be careful of –rf there are interesting stories out there
# storage viewing 
df -h 
du -hs
# concept of variables 
test="hello"
echo $test 
echo "${test} World"
# special variables 

echo $HOME # home directory
# $PATH shows you were to look for your executable 
echo $PATH | sed 's/:/\n/g' # we can explain a bit later
tree
export PATH="${HOME}/workshop.RNA-seq/bin/:$PATH"
# find out where your program is 
which tree
cd ~/workshop.RNA-seq/
# lets look at files now 
head -10 full_length.tab 
less full_length.tab # ESC to exit
# concept of pipes
cd ~/workshop.RNA-seq/seq
gunzip -c test_R1.fastq.gz | less -S
# viewing things in tab instead 
module load tabutils # this will load the system's program
cd ~/workshop.RNA-seq/
tabutils view full_length.tab | less -S
# searching for stuff
grep -irw './' -e 'kras'
# searching for stuff for a specific column 
awk '$2~/KRAS/' full_length.tab
# wild card
awk '$2~/KR.S/' full_length.tab
awk '$2~/^AS./' full_length.tab
# challenge question to pipe and save file
awk -F "\t" '$2~/^S.*[0..9]$/' full_length.tab | tabutils view > test.txt
# test this out
nano test.txt
# search and replace
echo "KRAS is a bad gene"
echo "KRAS is a bad gene" | sed s/bad/good/
head -10 full_length.tab
# lets replace SDF4 with KRAS
head -25 full_length.tab | sed s/SDF4/KRAS/
# the concept of greedy searches
echo "aaaaE_KRASaaaaaEEEEE" | sed s/aE.*//
echo "aaaaE_KRASaaaaaEEEEE" | sed s/aE.*?//
# limitation of sed
head -25 full_length.tab | sed s/E/_/g
# fix this with awk 
head -25 full_length.tab | awk '{ gsub("E","_",$2); print $2 "\t" $1 }' | less -S
# now lets take a peak at some files we need to know about
# what is a fastq 
cd ~/workshop.RNA-seq/seq
# what is a fasta file 
cd ~/workshop.RNA-seq/index
ls -lS
less -S cdna.fa
# GTF files
tabutils view gencode.v24.annotation.gtf | less -S
# do interactive 
srun -N 1 -n 1 --mem=6GB --time=6:0:0 --pty bash
squeue -u <username>
hostname # check if you are in the interactive mode
# use screen so that you you can go back to your session 
screen --ls
screen -S alex
# check resources
module avail
# load module 
module load STAR/2.5.3a
# submit your first job 
cd ~/workshop.RNA-seq/scripts
JOBID=$(sbatch --time=12:00:00 \
--export=ALL \
-e ./test/test.error \
-o ./test/test.o \
--mem=7GB \
-n 1 --wrap \
'gunzip -c ../seq/test_R1.fastq.gz | head -4000 | grep 0:TTAGGC | wc -l; echo "hello world"' \
| awk '{print $4}')
# submit with scripts
cd ~/workshop.RNA-seq/scripts
sbatch s.sh
squeue -u <user>
# final real world example 
nano run1.sh







