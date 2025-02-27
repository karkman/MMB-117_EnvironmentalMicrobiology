# Quality control and trimming

During this part we will do quality control on the amplicon sequence data and remove primer sequences from each amplicon.  
But first we need some setting up.  

## Setup

Log into Puhti using the web interface (you can use either "Login node shell" or "Visual Studio Code"). Go to the course project scratch folder and create your own folder with your username.  
You can check the path using command `csc-workspaces` that lists all projects yor affiliated with.  

```bash
csc-workspaces

cd /scratch/project_XXXXXX/
mkdir $USER
```

Check what was created and then go to your own folder and clone this Github repository there.  
We will run all analysis in this folder.  

```bash 
cd $USER
git clone https://github.com/karkman/MMB-117_EnvironmentalMicrobiology.git
cd MMB-117_EnvironmentalMicrobiology
```

Now we have a workig directory for the course that has all the instructions as well. So if you're using VS Code, you can open the instructions in VS Code. You can also check the full path to this folder to make it easier to navigate here in the future.  

Next we will make some folders and copy the raw sequencing data to the `01_RAW` folder.  

```bash
mkdir -p 00_LOGS
mkdir -p 01_RAW
cp /scratch/project_XXX/DATA/*.fastq.gz 01_RAW/
```

Check how many files were copied and does that correspond to what we should have.  

If all is good, make a list of all sample names. The one-liner first lists all R1 files, then splits each name by "-" and prints the 2nd and 3rd fields to file called `sample_names.txt`.  
After creating the file, you can check the content of it.  

```bash
ls 01_RAW/*R1*.fastq.gz | cut -d "-" -f 2,3 > sample_names.txt
```

## Sequence data quality control

For the next steps we need some computing resources (4h, 5G and 4 CPUs), so first applu the resources and then run the QC steps below.  
During the QC steps we first use `FastQC` to analyse each file separately and then combine all reports with `MultiQC`.  

```bash
module load biokit/11.3.0
mkdir 01_RAW/FASTQC

fastqc --outdir 01_RAW/FASTQC 01_RAW/*.fastq.gz --threads $SLURM_CPUS_PER_TASK
module purge

module load multiqc/1.19

multiqc --interactive --outdir 01_RAW/FASTQC 01_RAW/FASTQC/
module purge
```

Before going to the primer removal step, we'll go thru the multiQC reports together.  

## Primer removal

Since the amplification will create artificial sequences at the primer sites (all identical to the primers), we need to remove those before running any further analyses. You will need the primer sequences for this step. Also check the `cutadapt` manual for the other options were using.  

```bash
module load cutadapt/4.9
mkdir -p 02_TRIMMED

for sample in $(cat sample_names.txt); do 
    cutadapt \
        -g XXXX \
        -G YYYY \
        -O 10 \
        --discard-untrimmed \
        --cores $SLURM_CPUS_PER_TASK \
        01_RAW/*${sample}*_R1_001.fastq.gz \
        01_RAW/*${sample}*_R2_001.fastq.gz \
        -o 02_TRIMMED/${sample}_trimmed_1.fastq.gz \
        -p 02_TRIMMED/${sample}_trimmed_2.fastq.gz \
        > 02_TRIMMED/${sample}.log
done

module purge
```

Have a look at the cutadapt log files, did the reads contain the primers? How much data was removed?  

## Quality control for trimmed data

It's good practice to check the data quality once more to make sure everything is OK. Run the same QC steps to the trimmmed data. You need to fill in the correct paths to input and output files, and any relevant options for both steps.  

```bash
module load biokit/11.3.0
mkdir 02_TRIMMED/FASTQC

fastqc 
module purge

module load multiqc/1.19
multiqc
module purge
```

Now we should have read files with primers removed and we can proceed to the ASV part with DADA2 in R.
