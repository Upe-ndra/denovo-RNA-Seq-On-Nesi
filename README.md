# denovo-RNA-Seq-in-Nesi
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



