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

vcftools

SNAP

ANGSD

The easiest way to install the dependencies is using conda 

conda create --name MELFUwgrs -c conda-forge -c bioconda perl-parallel-loops parallel stacks adapterremoval bowtie2 bwa samtools bcftools vcftools angsd

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
                

