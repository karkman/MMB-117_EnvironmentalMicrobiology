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

Next we will make some folders and copy the raw sequencing data to the right folder.  

```bash
mkdir -p 00_LOGS
mkdir -p 01_RAW
cp /scratch/project_XXX/DATA/*.fastq.gz 01_RAW/
```

Make a list of all sample names. After creating the file, you can check the content of it.  

```bash
ls 01_RAW/*R1*.fastq.gz | cut -d "-" -f 2,3 > sample_names.txt
```

## Sequence data quality control

For the next steps we need some computing resources, so apply for resources and run the QC steps.  


```bash
module load biokit/11.3.0
mkdir 01_RAW/FASTQC

fastqc --outdir 01_RAW/FASTQC 01_RAW/*.fastq.gz --threads $SLURM_CPUS_PER_TASK
module purge

module load multiqc/1.19

multiqc --interactive --outdir 01_RAW/FASTQC 01_RAW/FASTQC/
module purge
```

Primer removal 

```bash
module load cutadapt/4.9
mkdir -p 02_TRIMMED

for sample in $(cat sample_names.txt); do 
    cutadapt \
    -g CCTACGGGNGGCWGCAG \
    -G GACTACHVGGGTATCTAATCC \
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

Trimmed QC

```bash
module load biokit/11.3.0
mkdir 02_TRIMMED/FASTQC

fastqc --outdir 02_TRIMMED/FASTQC 02_TRIMMED/*.fastq.gz --threads $SLURM_CPUS_PER_TASK
module purge

module load multiqc/1.19
multiqc --interactive --outdir 02_TRIMMED/FASTQC 02_TRIMMED/FASTQC
module purge
```
