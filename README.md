# denovo-RNA-Seq-on-Nesi
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
