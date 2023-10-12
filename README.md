# MELFUwgrs!
Pipeline to call snps using a reference genome

# Download
        git clone https://github.com/Yuma248/MELFUwgrs-.git
  
# Dependecies 

Perl Parallel:::Loops

GNU Parallel

AdapterRemoval

Bowtie2

BWA

samtools

bcftools

bedtools

vcftools

SNAP

ANGSD (including NGSadmix, PCAngsd, ngsLD)

The easiest way to install the dependencies is using conda 

conda create --name MELFUwgrs -c conda-forge -c bioconda perl-parallel-loops parallel stacks adapterremoval bowtie2 dragmap bwa samtools bedtools bcftools vcftools angsd

For SNAP 

    git clone https://github.com/amplab/snap.git 

    cd snap 

    make  

Then copy snap-aligment to your path or include snap folder in your $PATH 




# Basic usage:


Usage:

MELFUwgrs.pl 
        -stp <You need at least determine what steps you want to run>
        
                indref: <Indexs the reference genome with samtools, picard, bowtie2 and snap> 
                
                trim: <It will use AdapterRemoval to trim and filter reads> 
                
                concat: <It will concatenate fastq files of the same sample but different runs in one file> 
                
                alignment: <It will use bowtie2, bwa or snap to align reads to a reference genome> 
                
                dedup: <This step will sort sam/bam files, convert sam to bam (if necessary) and mask duplicates> 
                
                indelrea: <This step will locally realign indels, although this is not recommended any more> 
                
                bedmarkrep: <This step will mask repeat regions in the genome> 
                
                snpcalling: <This step will use ANGSD to simultaneously call and genotype SNP> 
                
                filtering: <This step will use vcftools to filter SNPs, I recommend to use this automatically to have an idea of your data, but play whit the parameters if you have the time> 
                
# Scripts

## callSFS
This script uses ANGSD to calculate SFS per individual and population, and then calculate He, and Fst based on these SFS. The required arguments are a popmap that assigns individuals to each pop, an input folder with the bam files of all the samples, and an output folder to save all the outputs.

Usage:

callSFS.pl 
        -i <inputfolder with all bam files to used>
        
         -o <output folder to save all the resutls>
         
         -pm <popmap, tab delimited files with sample names in the first column and corresponding pop in the second column>
         
         -rg <path to reference genome used to map the bam files>
         
Optional:
        -s <liste of sites to be considere, recomended to use pruned SNPs to avoid LD effects, this list should be in the site format from ANGSD>
        
        -lnc <number of runs in parallel, it works per sample pop and pairs of pop, default 10>
        
        -snc <number of CPU to use in pairwise Fst calculations, default 4>
        

Example:

calSFS -i /Yuma/mapped/ -o /Yuma/SFSout/ -pm /Yuma/popmap -rg /Yuma/genome/REFERENCE.fasta -s /Yuma/pruned/Somatic_prunedSNPs.list -snc 4 -lnc 15.
