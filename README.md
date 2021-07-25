# denovo-RNA-Seq-from-raw-reads-to-assembly-to-count-matrix-and-annotation
This is a compilation of scripts I have used for denovo RNA-Seq assembly and annotation and getting count metrics for differential gene expression analysis in New Zealand eScience Infrastructure **(Nesi)**.

Workflow follows the softwares listed below:
1. Fastqc on raw reads
2. Trimmomatic
3. Fastqc after QC
4. RCorrector
5. Running Trinity in two phases (Phase I & Phase II)
6. Salmon
7. Corset


## 1. Fastqc on raw reads
Fastqc is a program that gives you the report on the quality of your sequencing data. This is the first step we want to do with most of the sequencing data types after we receive them.

To run Fastqc, Let's first change our working directory to the directory where we have our sequencing data using `cd`.
```
cd /nesi/nobackup/uoo002752/RNA-seq  # Path can vary
```

### Lets create a slurm script
```
nano fastqc.sl
```
Above script will open a nano text editor with name fastqc.sl.
Copy and paste the script below, edit account, email address in the script below as needed.
We can save and close the nano text editor with `ctrl+o` command

```
#!/bin/bash -e
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name fastqc
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

mkdir fastqc_raw
module load FastQC
fastqc *.fastq.gz -o ./fastqc_raw
```
Now we can submit the job to **Nesi** hpc using `sbatch` command
```
sbatch fastqc.sl
```
We will get `html` files for each of our `fastq.gz` files in a sub folder `fastqc_raw` inside our working directory.
we can copy those `html` files to our local computer and open in our web browser and see how the data looks like, depending on the result we will perform `QC` using `trimmomatic` in our next step.

To copy those `html` files from `Nesi` to our Personal computer. lets open another terminal window and type.

```
scp mahuika:/nesi/nobackup/uoo002752/RNA-seq/fastqc_raw/*.html /path/to/a/folder/in/our/personal/computer
```
## 2. Trimmomatic
We checked how our data looks like (quality, adapter contamination, GC contents etc...) using fastqc tool above now its time to clean it up remove adapters, and low quality reads we will be using trimmomatic for this purpose, however  there are many other tools to do the same job.

We will have multiple fastq.gz files (raw data) to clean using trimmomatic, so we will be running it on a loupe. 
There are other ways (may be more efficient) to do it, but for now we will first create a text file called `names` listing the unique part of names of all the fastq.gz files. 
We can create this file using `ls` followed by piping `|` it to `sed` to edit the names to keep the unique part and saving it to a text file called `names`.
Basically, if the name of our fastq.gz files are as below:
```
EW459-W1-M_S20_L001_R1_001.fastq.gz
EW459-W1-M_S20_L001_R2_001.fastq.gz
EW459-W1-M_S20_L002_R1_001.fastq.gz
EW459-W1-M_S20_L002_R2_001.fastq.gz
EW460-W1-M_S21_L001_R1_001.fastq.gz
```
We will prepare a have following lines in `names'
```
EW459-W1-M_S20_L001
EW459-W1-M_S20_L001
EW459-W1-M_S20_L002
EW459-W1-M_S20_L002
EW460-W1-M_S21_L001
```

### Lets make slurm job script for trimmomatic now.
```
nano trimmomatic.sl
```
Copy, paste, and save `ctrl+o` the command below
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name trimmomatic
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

mkdir paired
mkdir unpaired
module load Trimmomatic/0.38-Java-1.8.0_144
for f in $(<names)
do
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -phred33 -threads 10 "${f}_R1_001.fastq.gz" "${f}_R2_001.fastq.gz" \
"paired/${f}_R1_001_trim.fastq.gz" "unpaired/${f}_R1_001_trim_unpaired.fastq.gz" \
"paired/${f}_R2_001_trim.fastq.gz" "unpaired/${f}_R2_001_trim_unparied.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
done
```
We can submit the job to `Nesi` using `sbatch` command as we did for fastqc job.
```
sbatch trimmomatic.sl
```
This will create two folders paired and unpaired you will have your clean data in paired folder for further use.
Check [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) for options and parameters we used and to fine tune them according to your need.
[Note] When I ran trimmomatic, I don't know why but, it had problem finding the illumina adapter file `TruSeq3-PE-2.fa` So I had to copy that file to my working directory. If you come across this problem, you can find adapter file in Trimmomatic installation directory.

## 3. Fastqc after QC
Now we want to see how our cleaned data looks like so for that purpose we will have to run fastqc on paired folder that was created after we ran trimmomatic. For this purpose we can use the same slurm script as we used for fastqc above. may be we just want to modify the later part to make it more meaningful like change
```
mkdir fastqc_raw
module load FastQC
fastqc *.fastq.gz -o ./fastqc_raw
```
to 
```
mkdir fastqc_trimmed
module load FastQC
fastqc *.fastq.gz -o ./fastqc_trimmed
```
So that the reports will be in new subfolder `fastqc_trimmed` inside `paired` folder
Now we can copy all the `.html` files in our personal computer as above and open them in web browser to compare how the data looks like before and after quality control with trimmomatic, make sure all the adapters and low quality reads have been removed.

For easy comparision of fastqc results we can additionally run `MultiQC` program. this will prepare a single report for all the report files in a directory so it will be easy to compare them.
To do that, we need to have all the `.html` files from first fastqc and second fastqc in one folder, we can move them to a new folder `multiqc`.
Lets create a folder `multiqc` in our parent working directory `RNA-seq`
```
cd /nesi/nobackup/uoo002752/RNA-seq # we can also use relative path cd ../..
mkdir multiqc
mv /fastqc_raw/*.html ./
mv /paired/fastq_trim/*.html ./
```
now we can run multiqc in this directory to summarize all the reports. 
```
nano multiqc.sl
```
copy, paste, and save the script below.
```
#!/bin/bash -e
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name multiqc
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load MultiQC
multiqc ./
```
run multiqc with `sbatch`
```
sbatch multiqc.sl
```
This will create a single `html` file which we can copy to our personal computer and view in web browser (as above).

## 4.Rcorrector
Rcorrector uses a k-mer based method to correct random sequencing errors in Illumina RNA-seq reads. It is used on the sequencing data where the read coverage is non-uniform.
Please read about [Rcorrector](https://academic.oup.com/gigascience/article/4/1/s13742-015-0089-y/2707778) before applying it on new data.
To run Rcorrector on our trimmomatic processed data, lets go to paired folder
```
cd /nesi/nobackup/uoo002752/RNA-seq/paired/
```
```
nano rcorrector.sl
```
copy, paste and save the below script to run rcorrector. This runs in a loop and for this we have to create a text file `names` same as the one for trimmomatic (see above)
```
#!/bin/bash -e

#SBATCH --job-name=rcorrector
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 10
#SBATCH --mem=50G
#SBATCH --partition=bigmem,large
#SBATCH --time=48:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Rcorrector
module load Jellyfish
for f in $(<names)
do
run_rcorrector.pl -t 10 -p “{f}_R1_001_trim.fastq.gz” “${f}_R2_001_trim.fastq.gz” -od ./rcorrector -verbose
done
```
## 5.Trinity
Now we will find our filtered (trimmomatic) and corrected (rcorrector) data files in a folder `rcorrector`.
We can now assemble this data using trinity.
We usually perform trinity in two phases in Nesi to make the process more efficient in the cluster in terms of resource usage. 
To run trinity we first need to create a configuration file `SLUM.conf`
Let's change our directory to rcorrector
```
cd rcorrector
```
```
nano SLUM.conf
```
Copy, paste, edit (eg, account) and save `ctrl+o` the following script
```
[GRID]
# grid type:
gridtype=SLURM
# template for a grid submission
# make sure:
# --partition is chosen appropriately for the resource requirement
# (here we choose either large or bigmem, whichever is available first)
# --ntasks and --cpius-per-task should always be 1
# --mem may need to be adjusted
# --time may need to adjusted
# (must be enough time for a batch of commands to finish)
# --account should be your NeSI project code
# add other sbatch options as required
cmd=sbatch --partition=large,bigmem --mem=5G --ntasks=1 --cpus-per-task=1 --time=10:00:00 --
account=uoo02752
#note -e error.file -o out.file are set internally, so dont set them in the above cmd.
####################################################################################
# settings below configure the Trinity job submission system, not tied to the grid itself.
####################################################################################
# number of grid submissions to be maintained at steady state by the Trinity submission system
max_nodes=100
# number of commands that are batched into a single grid submission job.
cmds_per_node=100
```
Lets create a slurm script for trinity phase I
```
nano trinity_phase_1.sl
```
Copy, paste, edit (eg, account, file names, SS_lib_type and other trinity parameters) and save `ctrl+o` the following script and run it using `sbatch trinity_phase_1.sl`
```
#!/bin/bash -e
#SBATCH --nodes 1
#SBATCH --job-name=trinity-phase1
#SBATCH --account=uoo02752
#SBATCH --time=30:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --partition=bigmem
#SBATCH --mem=220G
#SBATCH --hint=nomultithread
##SBATCH --output=%x.%j.out
##SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz

# load a Trinity module
#module load Trinity/2.11.0-gimkl-2020a
module load Trinity/2.8.6-gimkl-2020a #above version of trinity didn’t work for one of my data, so I had to go with older version, I don't know why it didn't work.
# run trinity, stop before phase 2
srun Trinity --no_distributed_trinity_exec \
--CPU ${SLURM_CPUS_PER_TASK} --max_memory 200G \
--seqType fq \
--left (comma separated list of R1 of paired end files here) \
--right (comma separated list of R2 of paired end files here) \
--include_supertranscripts \
--min_contig_length 200 \
--SS_lib_type RF \
--output trinity_out
```
When the run is finished you can run trinity Phase II. using following script, save the following script as `trinity_phase_2.sl` and submit the job using `sbatch`
```
nano trinity_phase_2.sl
```
```
#!/bin/bash -e
#SBATCH --nodes 1
#SBATCH --job-name=trinity-phase2
#SBATCH --account=uoo02752
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bigmem
#SBATCH --mem=20G
#SBATCH --hint=nomultithread
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz

# load a Trinity module
#module load Trinity/2.11.0-gimkl-2020a
module load Trinity/2.8.6-gimkl-2020a
module load HpcGridRunner/20181005
srun Trinity --CPU ${SLURM_CPUS_PER_TASK} --max_memory 20G \
--grid_exec "hpc_cmds_GridRunner.pl --grid_conf ${SLURM_SUBMIT_DIR}/SLURM.conf -c" \
--seqType fq \
--left (comma separated list of R1 of paired end files here) \
--right (comma separated list of R2 of paired end files here) \
--include_supertranscripts \
--min_contig_length 200 \
--SS_lib_type RF \
--output trinity_out
```
Upon completion of this job you will get Trinity.fasta as a final RNA-seq assembly in `trinity_out` folder.

## 6. Salmon to get count metrix
After denovo assembly we can use salmon to get the count metrics for genes and or trascripts and cluster them with corset.
Salmon first indexes the assembly and quantifies each pairdend files individually, we will be using for loop for it to go through all the samples
below is the script:

```
#!/bin/bash -e

##SBATCH --job-name=salmon
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 10
#SBATCH --mem=50G
#SBATCH --partition=bigmem
#SBATCH --time=12:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

#module load Salmon/1.3.0-gimkl-2020a

salmon index -t path/to/Trinity.fasta -i RNA_index

FILES=`ls *_L00A_R1_001_trim.fastq.gz | sed 's/_R1_001_trim.fastq.gz//g'`
for F in $FILES ; do
    R1=${F}_R1_001_trim.fastq.gz
    R2=${F}_R2_001_trim.fastq.gz
    salmon quant --threads 10 --gcBias --index 
RNA_index --libType A --dumpEq --hardFilter --skipQuant -1 $R1 -2 $R2 --output ${F}.out
done
```
## 7. corset
We now will use corset to cluster the genes/transcripts we can do so using script below.
This will give you a count matrix which you can use with DESeq2 or edgeR or other softwares in R to do differential gene expression analysis.
-g and -n flag above should be according your sample composition different numbers in -g represents different groups of samples to compare and -n are their corresponding names above I have four groups to compare group 1 has two samples AM1 and AM2 and like wise.

```
#!/bin/bash -e

##SBATCH --job-name=corset
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 10
#SBATCH --mem=50G
#SBATCH --partition=bigmem
#SBATCH --time=12:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread


module load Corset/1.09-GCC-9.2.0

corset -g 1,1,2,2,2,3,3,4,4,4,4 -n AM1,AM2,BML1,BML2,BML3,BMS1,BMS2,DM1,DM2,DM3,DM4 -i salmon_eq_classes *_L00A.out/aux_info/eq_classes.txt
```
## 8.Annotation

### 8.1. Transdecoder 

To annotate our denovo assembly `Trinity.fasta` we will first run transdecoder with the script below

```
#!/bin/bash -e
#SBATCH --job-name=transdec.nem
#SBATCH --account=uoo02752
#SBATCH --time=60:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --partition=bigmem
#SBATCH --error=%x.%j.err
#SBATCH --output=%x.%j.out
#SBATCH --hint=nomultithread
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --mail-type=ALL

module load TransDecoder/5.5.0-GCC-9.2.0-Perl-5.30.1

TransDecoder.LongOrfs -t path/to/Trinity.fasta
```
### 8.2 series of other programs to run to get the annotation, comments on the script below are quite explanatory. You have to run them one by one in order, not all at the same time. use comments `#` to stop scripts from running.

```
#!/bin/bash -e

#SBATCH --job-name=Trinotate
#SBATCH --account=uoo02752
#SBATCH --nodes 1 
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks 16
#SBATCH --mem=50G
##SBATCH --qos=debug
##SBATCH --time=00:15:00
#SBATCH --partition=bigmem,large
#SBATCH --time=72:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=All
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread


# modules to be used
module load Trinotate/3.2.1-GCC-9.2.0
module load SQLite/3.31.1-GCCcore-9.2.0
module load TransDecoder/5.5.0-GCC-9.2.0-Perl-5.30.1
module load HMMER/3.3-GCC-9.2.0

# run below scripts one by one
# Download databases (Uniprot)
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_sprot.trinotate_v2.0.pep.gz

# rename, uncompress and index
mv uniprot_sprot.trinotate_v2.0.pep.gz uniprot_sprot.trinotate.pep.gz
gunzip uniprot_sprot.trinotate.pep.gz

#Install ncbi blast+ with wget following ncbi instruction and prepare the protein database for blast search by:
/nesi/nobackup/uoo02752/bin/ncbi-blast-2.11.0+/bin/makeblastdb -in uniprot_sprot.pep -dbtype prot

# Download databases (Pfam)
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz

# gunzip and prepare for use with hmmscan
gunzip Pfam-A.hmm.gz

module load HMMER/3.3-GCC-9.2.0
hmmpress Pfam-A.hmm

# Download  and uncompress the Trinotate SQLite database
wget "https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Trinotate.sprot_uniref90.20150131.boilerplate.sqlite.gz" -O Trinotate.sqlite.gz
gunzip Trinotate.sqlite.gz

##Run transdecoder with trinity.fasta that will produce long_orf.pep file containing longest ORF peptide candidates from trinity assembly to be used under
module load TransDecoder/5.5.0-GCC-9.2.0-Perl-5.30.1
TransDecoder.LongOrfs -t path/to/Trinity.fasta
TransDecoder.Predict -t Trinity.fasta

#Run blastx for your trinity.fasta
/nesi/nobackup/uoo02752/bin/ncbi-blast-2.11.0+/bin/blastx -query Trinity.fasta -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6

#if you split your trinity.fasta using split_fasta.pl script you have to merge the above results with the script below. but not our case here
cat blastx.vol.*.outfmt6 > blastx.outfmt6

# now run blastp
/nesi/nobackup/uoo02752/bin/ncbi-blast-2.11.0+/bin/blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6

#same as above, if you had split the transdecoder.pep you have to merge above output with cat as below, not our case
cat blastp.vol.*.outfmt6 > blastp.outfmt6

# now run hmmer against Pfam to identify protein domains
hmmscan --cpu 16 --domtblout TrinotatePFAM.out ./Pfam-A.hmm Trinity.fasta.transdecoder.pep > pfam.log

# There are other optional steps to run blast with Uniref90 database and others I have not done that.
# Now combine our results and create annotation report # gene_transcript_map.txt used below should be created when you run Trinity assembly. so look in that out put folder for that file.
module load Trinotate/3.2.1-GCC-9.2.0
Trinotate Trinotate.sqlite init --gene_trans_map path/to/trinity/output/folder/gene_transcript_map.txt --transcript_fasta path/to/Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep 

# load blast results
module load Trinotate/3.2.1-GCC-9.2.0
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6

# load Pfam results
module load Trinotate/3.2.1-GCC-9.2.0
Trinotate Trinotate.sqlite Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

# Load other results from Uniref90, signalP and rnammer if you have run them. I have not.

# Now its time to create annotation report
module load Trinotate/3.2.1-GCC-9.2.0
Trinotate Trinotate.sqlite report -E 0.00001 > trinotate_annotation_report.xls

# trinotate_annotation_report.xls will have all the ids and annotation report for further use.
```
